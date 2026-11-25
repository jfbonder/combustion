"""Moving-mesh finite volume solver for (P_eps) in 1D / radial symmetry.

    u_t = Delta u - (u + w) f_eps(u),      w_t = Delta w,

on [0, L] (planar: x in [0, L]; radial: r = |x| in [0, L], Delta = r^{1-N}
(r^{N-1} u_r)_r), with zero flux at both ends. w = Y - u.

Space: vertex-centred finite volumes on a non-uniform mesh (conservative,
second order on smoothly graded meshes, handles r = 0 without special
cases because the flux weight r^{N-1} vanishes there).

Time: scipy's BDF (variable order, adaptive step, analytic sparse Jacobian).

Mesh adaptivity ("chasing the flame"): every remeshing interval we
  1. locate the front x_f (where u = eps; f_eps vanishes for u >= eps),
  2. build a band around it that covers the reaction zone plus a safety
     margin the front cannot cross before the next remesh,
  3. build a new mesh: uniform size h_fine inside the band, cells growing
     geometrically (ratio `growth`) away from it up to h_coarse; this costs
     only O(log(h_coarse/h_fine)/log(growth)) extra nodes, and refinement
     near the new front and coarsening behind the old one happen together,
  4. transfer u and w to the new nodes by monotone cubic interpolation.
The remeshing interval is adapted to the measured front speed so the front
moves at most half the margin between remeshes, and it is at most
max_dt_remesh (default 0.05 * t_end; see simulate).

Assumption: a single front (the case for the 1D / radial tests).
"""

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

from .model import resolve_reaction


# ---------------------------------------------------------------------------
# Failures
# ---------------------------------------------------------------------------

class NumericalFailure(RuntimeError):
    """A run could not be continued for a numerical reason (as opposed to a
    physical outcome such as the extinction of the flame, which is reported
    in Result.status). `result` holds the partial Result (status "failed",
    reason = the message), when the solver could build it."""
    result = None


class PecletError(NumericalFailure):
    """The cell Peclet condition of the ALE scheme is violated: the matrix of
    the step is not an M-matrix (see combustion1d.ale._matrix)."""


class MeshError(NumericalFailure):
    """The moving mesh is not valid (two nodes cross within one step)."""


# ---------------------------------------------------------------------------
# Discrete Laplacian
# ---------------------------------------------------------------------------

def laplacian(x, N):
    """Finite volume Laplacian (zero flux at both ends) on nodes x."""
    n = x.size
    h = np.diff(x)
    xm = 0.5 * (x[1:] + x[:-1])               # interior faces
    faces = np.concatenate(([x[0]], xm, [x[-1]]))
    if N == 1:
        vol = np.diff(faces)
        a = 1.0 / h
    else:
        vol = np.diff(faces**N) / N
        a = xm ** (N - 1) / h
    main = np.zeros(n)
    main[:-1] -= a
    main[1:] -= a
    lower = a / vol[1:]
    upper = a / vol[:-1]
    main /= vol
    return sp.diags([lower, main, upper], [-1, 0, 1], format="csc")


# ---------------------------------------------------------------------------
# Mesh
# ---------------------------------------------------------------------------

@dataclass
class MeshParams:
    pts_per_eps: float = 16.0   # h_fine = eps / pts_per_eps
    h_coarse: float = 0.01
    K_fresh: float = 5.0        # band extent on the fresh side (units of eps)
    K_burnt: float = 20.0       # band extent on the burnt side (units of eps)
    K_margin: float = 15.0      # safety margin on both sides (units of eps)
    growth: float = 1.1         # cell size ratio outside the band


def find_front(x, u, level):
    """Position of the first crossing of u = level (linear interpolation), and
    the side of the burnt region {u < level}: -1 = to the left of the front,
    +1 = to the right. Returns (None, 0) if u - level does not change sign.

    A node where u equals the level exactly counts as a crossing if u - level
    has opposite signs on its two sides (the front is then that node); a node
    that only touches the level (same sign on both sides) or a zero at an end
    of the mesh is not a crossing.
    """
    s = np.sign(u - level)
    if np.any(s == 0):
        # nodes exactly on the level: look for the sign change between the
        # neighbouring nonzero values
        nz = np.nonzero(s)[0]
        k = np.nonzero(s[nz[:-1]] * s[nz[1:]] < 0)[0]
        if k.size == 0:
            return None, 0
        i, j = nz[k[0]], nz[k[0] + 1]
        if j > i + 1:                   # x[i+1] (.. x[j-1]) lie on the level
            return float(x[i + 1]), (-1 if s[i] < 0 else +1)
    else:
        idx = np.nonzero(s[:-1] * s[1:] < 0)[0]
        if idx.size == 0:
            return None, 0
        i = idx[0]
    xf = x[i] + (level - u[i]) * (x[i + 1] - x[i]) / (u[i + 1] - u[i])
    burnt_left = u[i] < level
    return xf, (-1 if burnt_left else +1)


MIN_END_CELL = 0.5      # smallest cell at an end of [0, L], in units of h_fine


def build_mesh(x, u, eps, L, mp, anchor=0.0):
    """Band mesh around the front of (x, u): returns (nodes, front position).

    The band has spacing h_fine = eps / pts_per_eps and lies on the lattice
    anchor + h_fine * Z; outside it the cells grow geometrically up to
    h_coarse. All cells are >= MIN_END_CELL * h_fine (for h_coarse >= h_fine).
    Without a front the mesh is uniform with spacing h_coarse and the second
    value is None.
    """
    xf, burnt_dir = find_front(x, u, eps)
    if xf is None:
        n = int(np.ceil(L / mp.h_coarse)) + 1
        return np.linspace(0.0, L, n), None
    lo_K = mp.K_burnt if burnt_dir < 0 else mp.K_fresh
    hi_K = mp.K_fresh if burnt_dir < 0 else mp.K_burnt
    a = max(0.0, xf - (lo_K + mp.K_margin) * eps)
    b = min(L, xf + (hi_K + mp.K_margin) * eps)
    h_fine = eps / mp.pts_per_eps
    # snap the band to the lattice anchor + h_fine * Z: consecutive bands share
    # their nodes, so the transfer of the solution is exact inside the band
    # (interpolation only acts in the graded transition, where u is smooth)
    ia = int(np.floor((a - anchor) / h_fine))
    ib = int(np.ceil((b - anchor) / h_fine))
    ia = max(ia, int(np.ceil((0.0 - anchor) / h_fine)))
    ib = min(ib, int(np.floor((L - anchor) / h_fine)))
    band = anchor + h_fine * np.arange(ia, ib + 1)
    # A band that reaches an end of [0, L] leaves there a remainder shorter
    # than h_fine (the anchor is arbitrary with respect to 0 and L). A
    # remainder below MIN_END_CELL * h_fine is merged into the adjacent band
    # cell (the last lattice node is dropped), so that every cell of the mesh
    # is >= MIN_END_CELL * h_fine.
    if band.size > 2 and 0.0 < band[0] < MIN_END_CELL * h_fine:
        band = band[1:]
    if band.size > 2 and 0.0 < L - band[-1] < MIN_END_CELL * h_fine:
        band = band[:-1]
    a, b = band[0], band[-1]
    h0 = h_fine
    left = a - _graded(a, h0, mp.h_coarse, mp.growth)[::-1]
    right = b + _graded(L - b, h0, mp.h_coarse, mp.growth)
    xn = np.concatenate((left[:-1], band, right[1:]))
    xn[0], xn[-1] = 0.0, L          # remove round-off at the ends
    return xn, xf


def _graded(length, h0, h_max, growth):
    """Offsets 0 = d_0 < ... < d_m = length with cell sizes growing
    geometrically from h0 (factor `growth`) up to h_max."""
    if length <= 0.0:
        return np.array([0.0])
    d, h = [0.0], h0
    while d[-1] < length:
        h = min(h * growth, h_max)
        d.append(d[-1] + h)
    d = np.array(d)
    if len(d) > 2 and length - d[-2] < 0.5 * (d[-1] - d[-2]):
        d = d[:-1]          # merge a too-small last cell into the previous one
    return d * (length / d[-1])


# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------

@dataclass
class Result:
    t: np.ndarray
    front: np.ndarray
    snapshots: dict = field(default_factory=dict)  # t -> (x, u, w)
    n_remesh: int = 0
    n_steps: int = 0
    n_nodes: list = field(default_factory=list)
    mesh_history: list = field(default_factory=list)   # [(t, nodes)] if requested
    extinct: bool = False      # True iff status == "extinct"
    # "ok": the run reached t_end with a front;
    # "extinct": the front {u = eps} disappeared (ALE stops at that step; the
    #     Eulerian solvers continue on a uniform mesh), reason says when;
    # "failed": numerical failure, reason says which (the ALE solver raises
    #     the corresponding NumericalFailure, with this Result attached)
    status: str = "ok"
    reason: str = None
    n_secant_fail: int = 0     # ALE: steps accepted with |front - reference| >= tolerance
    n_mesh_remesh: int = 0     # ALE: remeshings forced by a cell below min_cell * h_f
    # ALE diagnostic: number of steps that ended with the front outside the
    # fine band (in a transition or coarse cell, where the reaction layer is
    # not resolved), and the time of the first one. Not acted upon.
    n_front_off_band: int = 0
    t_front_off_band: float = None
    n_solves: int = 0          # tridiagonal solves (scheme (S): ALE and simple.py)
    n_trials: int = 0          # ALE: trial velocities (each = 2 tridiagonal solves)
    n_lu: int = 0              # BDF solver: LU factorisations (2n x 2n sparse)
    n_fev: int = 0             # BDF solver: right-hand side evaluations
    v: list = field(default_factory=list)              # ALE: mesh velocity of each step
    dts: list = field(default_factory=list)            # simple.py: distinct step sizes used


def initial_data(x, u0, w0, eps):
    """u0^eps = u0 - eps w0 (1 - u0/eps)^+ ,  w = eps w0.

    Then (Y0 - u0^eps)/eps = w0 exactly, Y0 = u0^eps + eps w0 >= 0, the
    burnt region {u0 = 0} has Y0 = 0 and u0^eps -> u0 uniformly (Lipschitz).
    """
    U = u0(x)
    W0 = w0(x) if callable(w0) else np.full_like(x, w0)
    u = U - eps * W0 * np.maximum(0.0, 1.0 - U / eps)
    return u, eps * W0


def _mark_extinct(res, t, u, eps):
    """Record in res that no front was found at time t (first time only)."""
    if res.status != "ok":
        return
    if not np.all(np.isfinite(u)):
        res.status, res.reason = "failed", f"non-finite values in the solution at t = {t:.6g}"
        return
    side = "u > eps" if u[0] > eps else "u < eps"
    res.status, res.extinct = "extinct", True
    res.reason = f"front disappeared: {side} on the whole mesh at t = {t:.6g}"


def simulate(u0, w0, eps, L, t_end, N=1, mesh=None, speed_guess=1.0,
             snapshot_times=(), rtol=1e-5, atol_rel=1e-4, max_dt_remesh=None,
             verbose=False, record_mesh=False, reaction=None):
    """Reference solver: adaptive BDF on band meshes rebuilt around the front.

    The mesh is rebuilt (and the BDF integrator restarted) after every time
    interval of length
        min(0.5 * K_margin * eps / v_est, max_dt_remesh, t_end - t),
    shortened to land on the next snapshot time; v_est is the front speed
    measured over the previous interval (speed_guess for the first one,
    never below 1e-3). max_dt_remesh defaults to 0.05 * t_end, so the
    sequence of meshes depends on the final time unless max_dt_remesh is
    given. res.t / res.front hold one sample per interval.

    reaction: None (the rate installed with model.set_reaction) or f or
        (f, df) for this call only (see model.resolve_reaction).
    """
    mp = MeshParams() if mesh is None else mesh
    if max_dt_remesh is None:
        max_dt_remesh = 0.05 * t_end
    if not max_dt_remesh > 0:
        raise ValueError("max_dt_remesh must be positive")
    f_r, df_r = resolve_reaction(reaction)

    def f_eps(s, e):
        return f_r(s / e) / e**2

    def df_eps(s, e):
        return df_r(s / e) / e**3

    x = np.linspace(0.0, L, 20001)
    u, w = initial_data(x, u0, w0, eps)
    x_new, xf = build_mesh(x, u, eps, L, mp)
    u, w = initial_data(x_new, u0, w0, eps)   # exact initial data on mesh
    x = x_new
    # front of the initial datum on the actual mesh (build_mesh located it on
    # the auxiliary uniform grid, whose spacing L/20000 can exceed eps)
    xf = find_front(x, u, eps)[0]

    t = 0.0
    ts, fronts = [0.0], [xf if xf is not None else np.nan]
    res = Result(t=None, front=None)
    pending = sorted(snapshot_times)
    v_est = speed_guess
    h_last = None
    if pending and pending[0] == 0.0:
        res.snapshots[0.0] = (x.copy(), u.copy(), w.copy())
        pending.pop(0)

    while t < t_end - 1e-14:
        n = x.size
        Lap = laplacian(x, N)

        def rhs(_, y):
            uu, ww = y[:n], y[n:]
            return np.concatenate((Lap @ uu - (uu + ww) * f_eps(uu, eps), Lap @ ww))

        def jac(_, y):
            uu, ww = y[:n], y[n:]
            fe = f_eps(uu, eps)
            d = fe + (uu + ww) * df_eps(uu, eps)
            return sp.bmat([[Lap - sp.diags(d), -sp.diags(fe)],
                            [None, Lap]], format="csc")

        dt = min(0.5 * mp.K_margin * eps / max(v_est, 1e-12), max_dt_remesh, t_end - t)
        if pending and t + dt >= pending[0]:
            dt = pending[0] - t
        sol = solve_ivp(rhs, (t, t + dt), np.concatenate((u, w)), method="BDF",
                        jac=jac, rtol=rtol, atol=atol_rel * eps,
                        first_step=min(h_last, dt) if h_last is not None else None)
        if not sol.success:
            res.status, res.reason = "failed", f"BDF integrator at t = {t:.6g}: {sol.message}"
            res.t, res.front = np.array(ts), np.array(fronts)
            err = NumericalFailure(res.reason)
            err.result = res
            raise err
        t = sol.t[-1]
        if sol.t.size > 2:
            h_last = sol.t[-2] - sol.t[-3]   # last full step
        res.n_steps += sol.t.size - 1
        res.n_lu += sol.nlu
        res.n_fev += sol.nfev
        u, w = sol.y[:n, -1], sol.y[n:, -1]

        xf_new, _ = find_front(x, u, eps)
        if xf_new is not None and xf is not None:
            v_est = max(abs(xf_new - xf) / dt, 1e-3)
        xf = xf_new
        if xf is None:
            _mark_extinct(res, t, u, eps)
        ts.append(t)
        fronts.append(xf if xf is not None else np.nan)
        if pending and abs(t - pending[0]) < 1e-12:
            res.snapshots[pending.pop(0)] = (x.copy(), u.copy(), w.copy())

        # remesh: refine around the new front position, coarsen elsewhere
        x_new, _ = build_mesh(x, u, eps, L, mp)
        with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
            u = PchipInterpolator(x, u)(x_new)   # flat regions trigger
            w = PchipInterpolator(x, w)(x_new)   # harmless 1/0 warnings
        x = x_new
        res.n_remesh += 1
        res.n_nodes.append(x.size)
        if record_mesh:
            res.mesh_history.append((t, x.copy()))
        if verbose and res.n_remesh % 200 == 0:
            print(f"  t={t:.4f} front={xf} nodes={x.size} v={v_est:.3f}")

    res.t, res.front = np.array(ts), np.array(fronts)
    return res


def front_speed(res, t_from):
    """Least squares slope of the front position for t >= t_from."""
    m = (res.t >= t_from) & np.isfinite(res.front)
    return np.polyfit(res.t[m], res.front[m], 1)[0]
