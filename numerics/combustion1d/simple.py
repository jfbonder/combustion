"""The "simple" (provably monotone) version of the method.

Linearly implicit Euler with the reaction evaluated as Y^{n+1} f_eps(U^n):

    (Y^{n+1} - Pi Y^n)/dt = Delta_h Y^{n+1} - f_eps(U^n) Y^{n+1},
    (U^{n+1} - Pi U^n)/dt = Delta_h U^{n+1} - f_eps(U^n) Y^{n+1},

on the same front-following meshes as the main solver, with Pi = piecewise
LINEAR interpolation at remeshing (a positive operator). W = Y - U then
solves the discrete heat equation exactly.

If f is nonincreasing, the scheme satisfies a discrete comparison principle
for every dt > 0: the matrix (1/dt - Delta_h + diag f_eps(U^n)) is an
M-matrix and the reaction term is order preserving as long as Y >= 0 (see
the paper). Each step costs two tridiagonal solves.
"""

import numpy as np
from scipy.linalg import solve_banded

from .model import f_eps
from .solver import MeshParams, Result, build_mesh, find_front, initial_data, laplacian


def _banded(Lap, diag_extra, dt):
    """Banded storage of (1/dt) I - Lap + diag(diag_extra)."""
    n = Lap.shape[0]
    ab = np.zeros((3, n))
    ab[0, 1:] = -Lap.diagonal(1)
    ab[1, :] = 1.0 / dt - Lap.diagonal(0) + diag_extra
    ab[2, :-1] = -Lap.diagonal(-1)
    return ab


def simulate_simple(u0, w0, eps, L, t_end, N=1, dt=None, dt_factor=1.0, mesh=None,
                    speed_guess=1.0, snapshot_times=(), max_dt_remesh=None):
    """dt defaults to dt_factor * eps^2."""
    mp = mesh or MeshParams()
    dt = dt or dt_factor * eps**2
    max_dt_remesh = max_dt_remesh or 0.05 * t_end
    x0 = np.linspace(0.0, L, 20001)
    u, w = initial_data(x0, u0, w0, eps)
    x, xf = build_mesh(x0, u, eps, L, mp)
    u, w = initial_data(x, u0, w0, eps)
    y = u + w

    res = Result(t=None, front=None)
    t, ts, fronts = 0.0, [0.0], [xf]
    pending = sorted(snapshot_times)
    v_est = speed_guess
    while t < t_end - 1e-14:
        Lap = laplacian(x, N)
        heat = _banded(Lap, np.zeros(x.size), dt)
        interval = min(0.5 * mp.K_margin * eps / max(v_est, 1e-12), max_dt_remesh, t_end - t)
        if pending and t + interval >= pending[0]:
            interval = pending[0] - t
        nsteps = max(1, int(np.ceil(interval / dt - 1e-9)))
        h = interval / nsteps
        heat = _banded(Lap, np.zeros(x.size), h)
        for _ in range(nsteps):
            fe = f_eps(u, eps)
            y = solve_banded((1, 1), _banded(Lap, fe, h), y / h)
            u = solve_banded((1, 1), heat, u / h - fe * y)
            res.n_steps += 1
        t += interval
        xf_new, _ = find_front(x, u, eps)
        if xf_new is not None and xf is not None:
            v_est = max(abs(xf_new - xf) / interval, 1e-3)
        xf = xf_new
        ts.append(t)
        fronts.append(xf if xf is not None else np.nan)
        if pending and abs(t - pending[0]) < 1e-12:
            res.snapshots[pending.pop(0)] = (x.copy(), u.copy(), (y - u).copy())
        x_new, _ = build_mesh(x, u, eps, L, mp)
        u, y = np.interp(x_new, x, u), np.interp(x_new, x, y)   # positive operator
        x = x_new
        res.n_remesh += 1
        res.n_nodes.append(x.size)
    res.t, res.front = np.array(ts), np.array(fronts, dtype=float)
    return res
