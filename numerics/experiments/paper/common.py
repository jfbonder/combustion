"""Shared setup, I/O and plotting helpers for the paper experiments.

Run the group scripts from numerics/, e.g.
    OMP_NUM_THREADS=1 ../.venv/bin/python -m experiments.paper.meshes
"""

import json
import time
from pathlib import Path

import numpy as np

import combustion1d.ale as _ale
from combustion1d.exact import SelfSimilar, tw_profile_limit, tw_speed_limit
from combustion1d.model import fb_slope
from combustion1d.solver import MeshParams, Result

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


# ---------------------------------------------------------------------------
# ALE runs that may end with the disappearance of the front
# ---------------------------------------------------------------------------

class _KeepResult(Result):
    """Result that remembers the last instance (combustion1d.ale.simulate_ale
    raises a TypeError when the front disappears; with this we keep the
    snapshots / mesh history recorded up to that step)."""
    last = None

    def __post_init__(self):
        _KeepResult.last = self


def simulate_ale_safe(*args, **kw):
    """simulate_ale; on failure returns (partial Result, error string).

    In the partial result t / front are None: only snapshots, n_nodes and
    mesh_history are available."""
    old = _ale.Result
    _ale.Result = _KeepResult
    try:
        return _ale.simulate_ale(*args, **kw), None
    except Exception as ex:          # front disappeared (extinction), etc.
        return _KeepResult.last, repr(ex)
    finally:
        _ale.Result = old


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def save_json(name, obj):
    (RES / name).write_text(json.dumps(obj, indent=1, default=_conv))


def load_json(name):
    return json.loads((RES / name).read_text())


def _conv(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    raise TypeError(type(o))


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


def title(ax, s):
    ax.set_title(s, fontsize=9, color=INK, loc="left")


def savefig(fig, name):
    fig.savefig(FIG / name, dpi=DPI)
    plt.close(fig)
    print("wrote", FIG / name, flush=True)
