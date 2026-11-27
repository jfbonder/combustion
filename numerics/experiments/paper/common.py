"""Shared setup, I/O and plotting helpers for the scripts of experiments/paper
(not a script itself): the test problems, save_json / load_json (results/paper,
with the provenance record), the node histories (*.npz), the mesh analysis and
the plot helpers.

The scripts are run from numerics/, e.g.
    OMP_NUM_THREADS=1 python -m experiments.paper.meshes
experiments/paper/README.md lists them with the items of the paper they produce.
"""

import json
import time
from pathlib import Path

import numpy as np

import combustion1d.ale as _ale
from combustion1d.exact import SelfSimilar, tw_profile_limit, tw_speed_limit
from combustion1d.model import fb_slope
from combustion1d.solver import MeshParams, NumericalFailure, Result
from experiments.util import to_jsonable, with_provenance, worker_pool  # noqa: F401

# palette / rcParams of the existing figures (experiments/plots.py)
from experiments.plots import BLUES, GRID, INK, MARKERS, MUTED, SERIES, plt  # noqa: F401

ROOT = Path(__file__).resolve().parent.parent.parent
RES = ROOT / "results" / "paper"
FIG = ROOT / "figures" / "paper"
RES.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
DPI = 200


# ---------------------------------------------------------------------------
# Test problems
# ---------------------------------------------------------------------------

A_TW = 1.0
X0_TW = 2.0


def tw_setup(w0, T, A=A_TW, x0=X0_TW):
    """Planar travelling wave: limit profile as initial datum, domain long
    enough that the fresh end does not feel the front."""
    c0 = float(tw_speed_limit(A, w0))
    u0 = lambda x: tw_profile_limit(x, 0.0, A, c0, x0=x0)
    L = x0 + c0 * T + 12.0 / c0 + 2.0
    return u0, L, c0


def extinction_setup(N, w0=0.5, T=1.0):
    """Self-similar shrinking ball of the limit problem (exact)."""
    V = SelfSimilar(N, w0, T)
    L = 1.5 * V.R * np.sqrt(T) + 1.0
    u0 = lambda r: V(r, 0.0)
    return V, u0, L


def expanding_u0(r0, w0=0.5, A=1.0):
    lam = float(fb_slope(w0))
    return lambda r: np.where(r > r0, A * np.tanh(lam * (r - r0) / A), 0.0)


# Rigid-mesh experiments (the whole mesh translates: ell0 = 1e6,
# remesh_distance = 1e9). By construction the cells in the tapered zones at
# the ends of [0, L] are compressed far below h_f (the mesh moves, the ends do
# not), and the experiments rely on the mesh never being rebuilt. They pass
# min_cell = RIGID_MIN_CELL to simulate_ale, which switches off the remeshing
# of compressed cells (crossing nodes still raise MeshError); the travel
# distance is limited by hand in those scripts.
RIGID_MIN_CELL = 0.0


# ---------------------------------------------------------------------------
# ALE runs that may end with the disappearance of the front
# ---------------------------------------------------------------------------

class _KeepResult(Result):
    """Result that also records the step index of every remeshing
    (res.remesh_steps: indices into res.t), by intercepting the updates of
    the counter n_remesh."""

    def __init__(self, *a, **k):
        self._nrm = 0
        self.remesh_steps = []
        super().__init__(*a, **k)

    @property
    def n_remesh(self):
        return self._nrm

    @n_remesh.setter
    def n_remesh(self, val):
        if val > self._nrm:
            self.remesh_steps.append(len(self.n_nodes) + 1)   # = index into res.t
        self._nrm = val


def simulate_ale_safe(*args, **kw):
    """simulate_ale with remeshing times recorded (res.remesh_steps: indices
    into res.t). Returns (res, error):

    * error is None when the run ended normally: res.status is "ok", or
      "extinct" if the front disappeared (a physical outcome; res.reason
      says when, and res.front ends with nan);
    * a numerical failure of the scheme (combustion1d.solver.NumericalFailure:
      Peclet condition, crossing nodes) is returned as (partial res with
      status "failed", message). The caller must treat it as a failed run,
      never as an outcome of the model.

    Any other exception (programming errors, invalid arguments, failures of
    the linear solver) is not caught.
    """
    old = _ale.Result
    _ale.Result = _KeepResult
    try:
        return _ale.simulate_ale(*args, **kw), None
    except NumericalFailure as ex:
        if ex.result is None:
            raise
        return ex.result, f"{type(ex).__name__}: {ex}"
    finally:
        _ale.Result = old


def cycle_speed(res, t_from):
    """Mean front speed between the first and the last remeshing after t_from
    (whole remeshing cycles: the small periodic oscillation of the front
    position over a cycle cancels). Falls back to least squares."""
    idx = [i for i in res.remesh_steps if res.t[i] >= t_from]
    if len(idx) >= 2:
        a, b = idx[0], idx[-1]
        return (res.front[b] - res.front[a]) / (res.t[b] - res.t[a]), len(idx) - 1
    m = res.t >= t_from
    return np.polyfit(res.t[m], res.front[m], 1)[0], 0


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def save_json(name, obj):
    """Write obj to results/paper/<name> with its provenance (code commit and
    library versions; see experiments.util.with_provenance)."""
    (RES / name).write_text(json.dumps(with_provenance(obj), indent=1, default=_conv))


def load_json(name):
    return json.loads((RES / name).read_text())


_conv = to_jsonable


def save_mesh_history(name, hist, **extra):
    """[(t, x)] -> npz with concatenated nodes and offsets (float64)."""
    ts = np.array([h[0] for h in hist])
    sizes = np.array([h[1].size for h in hist])
    xs = np.concatenate([h[1] for h in hist])
    np.savez_compressed(RES / name, t=ts, sizes=sizes, x=xs, **extra)


def load_mesh_history(name):
    d = np.load(RES / name)
    off = np.concatenate(([0], np.cumsum(d["sizes"])))
    hist = [(d["t"][k], d["x"][off[k]:off[k + 1]]) for k in range(d["t"].size)]
    extra = {k: d[k] for k in d.files if k not in ("t", "sizes", "x")}
    return hist, extra


class Timer:
    def __enter__(self):
        self.t0 = time.time()
        return self

    def __exit__(self, *a):
        self.seconds = time.time() - self.t0


# ---------------------------------------------------------------------------
# Mesh analysis
# ---------------------------------------------------------------------------

def mesh_parts(x, eps, mp):
    """Number of nodes in the band (spacing h_f), transition and coarse zone."""
    h = np.diff(x)
    h_f = eps / mp.pts_per_eps
    band = np.abs(h - h_f) < 1e-6 * h_f
    coarse = h > 0.999 * mp.h_coarse
    nb = int(band.sum()) + 1
    nc = int(coarse.sum())
    return dict(total=int(x.size), band=nb, coarse=nc, transition=int(x.size - nb - nc))


def mesh_parts_predicted(eps, mp, L):
    h_f = eps / mp.pts_per_eps
    band = int(round((mp.K_burnt + mp.K_fresh + 2 * mp.K_margin) * mp.pts_per_eps)) + 1
    m = int(np.ceil(np.log(mp.h_coarse / h_f) / np.log(mp.growth)))
    trans_len = 2 * sum(min(h_f * mp.growth ** k, mp.h_coarse) for k in range(1, m + 1))
    band_len = (mp.K_burnt + mp.K_fresh + 2 * mp.K_margin) * eps
    coarse = (L - band_len - trans_len) / mp.h_coarse
    return dict(band=band, transition=2 * m, coarse=coarse, total=band + 2 * m + coarse)


def remesh_times(hist, vmax_dt):
    """Times of the records where the node set was rebuilt (size change or a
    node displacement larger than the ALE step can produce)."""
    out = []
    for (t0, x0), (t1, x1) in zip(hist[:-1], hist[1:]):
        if x0.size != x1.size or np.max(np.abs(x1 - x0)) > 1.5 * vmax_dt + 1e-12:
            out.append(t1)
    return np.array(out)


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def slope_ref(ax, x, y0, label="slope 1", p=1.0, **kw):
    x = np.array([min(x), max(x)])
    y = y0 * (x / x[0]) ** p
    ax.plot(x, y, ls="--", lw=1, color=MUTED)
    ax.annotate(label, (x[1], y[1]), xytext=(4, 0), textcoords="offset points",
                color=MUTED, va="center", **kw)


def plain_ticks(ax, ticks, axis="x", fmt="{:g}"):
    """Explicit ticks with plain labels on a log axis (no minor tick labels)."""
    getattr(ax, f"set_{axis}ticks")(ticks)
    getattr(ax, f"set_{axis}ticklabels")([fmt.format(t) for t in ticks])
    getattr(ax, f"{axis}axis").set_minor_locator(plt.NullLocator())


def title(ax, s):
    ax.set_title(s, fontsize=9, color=INK, loc="left")


def savefig(fig, name):
    fig.savefig(FIG / name, dpi=DPI)
    plt.close(fig)
    print("wrote", FIG / name, flush=True)
