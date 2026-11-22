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
moves at most half the margin between remeshes.

Assumption: a single front (the case for the 1D / radial tests).
"""

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

from .model import df_eps, f_eps


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
    pts_per_eps: float = 8.0    # h_fine = eps / pts_per_eps
    h_coarse: float = 0.01
    K_fresh: float = 5.0        # band extent on the fresh side (units of eps)
    K_burnt: float = 20.0       # band extent on the burnt side (units of eps)
    K_margin: float = 15.0      # safety margin on both sides (units of eps)
    growth: float = 1.1         # cell size ratio outside the band


def find_front(x, u, level):
    """Position of the (single) crossing u = level, or None."""
    s = np.sign(u - level)
    idx = np.nonzero(s[:-1] * s[1:] < 0)[0]
    if idx.size == 0:
        return None, 0
    i = idx[0]
    xf = x[i] + (level - u[i]) * (x[i + 1] - x[i]) / (u[i + 1] - u[i])
    burnt_left = u[i] < level
    return xf, (-1 if burnt_left else +1)


def build_mesh(x, u, eps, L, mp):
    xf, burnt_dir = find_front(x, u, eps)
    if xf is None:
        n = int(np.ceil(L / mp.h_coarse)) + 1
        return np.linspace(0.0, L, n), None
    lo_K = mp.K_burnt if burnt_dir < 0 else mp.K_fresh
    hi_K = mp.K_fresh if burnt_dir < 0 else mp.K_burnt
    a = max(0.0, xf - (lo_K + mp.K_margin) * eps)
    b = min(L, xf + (hi_K + mp.K_margin) * eps)
    h_fine = eps / mp.pts_per_eps
    # snap the band to the global lattice h_fine * Z: consecutive bands share
    # their nodes, so the transfer of the solution is exact inside the band
    # (interpolation only acts in the graded transition, where u is smooth)
    ia = int(np.floor(a / h_fine))
    ib = int(np.ceil(b / h_fine))
    ib = min(ib, int(np.floor(L / h_fine)))
    band = h_fine * np.arange(ia, ib + 1)
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


def initial_data(x, u0, w0, eps):
    """u0^eps = u0 - eps w0 (1 - u0/eps)^+ ,  w = eps w0.

    Then (Y0 - u0^eps)/eps = w0 exactly, Y0 = u0^eps + eps w0 >= 0, the
    burnt region {u0 = 0} has Y0 = 0 and u0^eps -> u0 uniformly (Lipschitz).
    """
    U = u0(x)
    W0 = w0(x) if callable(w0) else np.full_like(x, w0)
    u = U - eps * W0 * np.maximum(0.0, 1.0 - U / eps)
    return u, eps * W0


def simulate(u0, w0, eps, L, t_end, N=1, mesh=None, speed_guess=1.0,
             snapshot_times=(), rtol=1e-5, atol_rel=1e-4, max_dt_remesh=None,
             verbose=False):
    mp = mesh or MeshParams()
    max_dt_remesh = max_dt_remesh or 0.05 * t_end
    x = np.linspace(0.0, L, 20001)
    u, w = initial_data(x, u0, w0, eps)
    x_new, xf = build_mesh(x, u, eps, L, mp)
    u, w = initial_data(x_new, u0, w0, eps)   # exact initial data on mesh
    x = x_new

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
                        first_step=min(h_last, dt) if h_last else None)
        if not sol.success:
            raise RuntimeError(sol.message)
        t = sol.t[-1]
        if sol.t.size > 2:
            h_last = sol.t[-2] - sol.t[-3]   # last full step
        res.n_steps += sol.t.size - 1
        u, w = sol.y[:n, -1], sol.y[n:, -1]

        xf_new, _ = find_front(x, u, eps)
        if xf_new is not None and xf is not None:
            v_est = max(abs(xf_new - xf) / dt, 1e-3)
        xf = xf_new
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
        if verbose and res.n_remesh % 200 == 0:
            print(f"  t={t:.4f} front={xf} nodes={x.size} v={v_est:.3f}")

    res.t, res.front = np.array(ts), np.array(fronts)
    return res


def front_speed(res, t_from):
    """Least squares slope of the front position for t >= t_from."""
    m = (res.t >= t_from) & np.isfinite(res.front)
    return np.polyfit(res.t[m], res.front[m], 1)[0]
