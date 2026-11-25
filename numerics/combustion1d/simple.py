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

from .model import resolve_reaction
from .solver import (MeshParams, Result, _mark_extinct, build_mesh, find_front,
                     initial_data, laplacian)


def _banded(Lap, diag_extra, dt):
    """Banded storage of (1/dt) I - Lap + diag(diag_extra)."""
    n = Lap.shape[0]
    ab = np.zeros((3, n))
    ab[0, 1:] = -Lap.diagonal(1)
    ab[1, :] = 1.0 / dt - Lap.diagonal(0) + diag_extra
    ab[2, :-1] = -Lap.diagonal(-1)
    return ab


def simulate_simple(u0, w0, eps, L, t_end, N=1, dt=None, dt_factor=1.0, mesh=None,
                    speed_guess=1.0, snapshot_times=(), max_dt_remesh=None,
                    exact_initial=False, monitor=None, reaction=None):
    """dt defaults to dt_factor * eps^2.

    Every step has exactly the nominal size dt, except possibly the last one
    (to land on t_end) and the steps that land on a snapshot time; the step
    sizes used are recorded in res.dts. The mesh is rebuilt every
    nsteps = floor(interval/dt) steps (at least one step), where
        interval = min(0.5 K_margin eps / v_est, max_dt_remesh)
    and v_est is the front speed measured over the previous interval
    (speed_guess for the first one). max_dt_remesh defaults to 0.05 * t_end,
    so the sequence of meshes depends on the final time unless max_dt_remesh
    is given. res.t / res.front hold one sample per interval (res.front[0]
    is the front of the initial datum on the initial mesh).

    exact_initial: if True, u(x, 0) = u0(x) and w = eps*w0 (no modification).
    monitor: optional callable monitor(t, x, u, v=None), called at t = 0 and
        after every step (before remeshing).
    reaction: None (the rate installed with model.set_reaction) or a
        callable f for this call only (see model.resolve_reaction).
    """
    mp = MeshParams() if mesh is None else mesh
    if dt is None:
        dt = dt_factor * eps**2
    if max_dt_remesh is None:
        max_dt_remesh = 0.05 * t_end
    if not dt > 0 or not max_dt_remesh > 0:
        raise ValueError("dt and max_dt_remesh must be positive")
    f_r = resolve_reaction(reaction)[0]

    def f_eps(s, e):
        return f_r(s / e) / e**2

    def init(xx):
        if exact_initial:
            return u0(xx), np.full_like(xx, eps * w0)
        return initial_data(xx, u0, w0, eps)
    x0 = np.linspace(0.0, L, 20001)
    u, w = init(x0)
    x, xf = build_mesh(x0, u, eps, L, mp)
    u, w = init(x)
    y = u + w
    # front of the initial datum on the actual mesh (build_mesh located it on
    # the auxiliary uniform grid x0, whose spacing L/20000 can exceed eps)
    xf = find_front(x, u, eps)[0]

    res = Result(t=None, front=None)
    t, ts, fronts = 0.0, [0.0], [xf if xf is not None else np.nan]
    pending = sorted(snapshot_times)
    v_est = speed_guess
    dts = set()
    if monitor is not None:
        monitor(0.0, x, u, None)
    while t < t_end - 1e-12:
        Lap = laplacian(x, N)
        interval = min(0.5 * mp.K_margin * eps / max(v_est, 1e-12), max_dt_remesh)
        nsteps = max(1, int(np.floor(interval / dt + 1e-9)))
        t_stop = min(t + nsteps * dt, t_end)
        if pending and t_stop >= pending[0] - 1e-12:
            t_stop = pending[0]
        t_start = t
        heat, h_heat = None, None
        while t < t_stop - 1e-12:
            h = min(dt, t_stop - t)
            if h != h_heat:
                heat, h_heat = _banded(Lap, np.zeros(x.size), h), h
            fe = f_eps(u, eps)
            y = solve_banded((1, 1), _banded(Lap, fe, h), y / h)
            u = solve_banded((1, 1), heat, u / h - fe * y)
            res.n_steps += 1
            res.n_solves += 2
            dts.add(float('%.10g' % h))
            t = t + h if t_stop - (t + h) > 1e-12 else t_stop
            if monitor is not None:
                monitor(t, x, u, None)
        xf_new, _ = find_front(x, u, eps)
        if xf_new is not None and xf is not None:
            v_est = max(abs(xf_new - xf) / (t - t_start), 1e-3)
        xf = xf_new
        if xf is None:
            _mark_extinct(res, t, u, eps)
        ts.append(t)
        fronts.append(xf if xf is not None else np.nan)
        if pending and abs(t - pending[0]) < 1e-12:
            res.snapshots[pending.pop(0)] = (x.copy(), u.copy(), (y - u).copy())
        x_new, _ = build_mesh(x, u, eps, L, mp)
        u, y = np.interp(x_new, x, u), np.interp(x_new, x, y)   # positive operator
        x = x_new
        res.n_remesh += 1
        res.n_nodes.append(x.size)
    res.dts = sorted(dts)
    res.t, res.front = np.array(ts), np.array(fronts, dtype=float)
    return res
