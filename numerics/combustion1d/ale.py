"""Monotone scheme (S) in ALE form: the mesh follows the flame (planar, N = 1).

Nodes move with velocity  xdot_i = v(t) chi(x_i),  chi = 1 on a neighbourhood of the fine band,
decaying linearly to 0 over a distance ell0 (nodes far from the front do not
move). Along node trajectories  D_t u = u_t + xdot u_x, so

    (Y^{n+1} - Y^n)/dt = Delta_h Y^{n+1} + xdot D_x Y^{n+1} - f_eps(U^n) Y^{n+1},
    (U^{n+1} - U^n)/dt = Delta_h U^{n+1} + xdot D_x U^{n+1} - f_eps(U^n) Y^{n+1},

on the mesh x^{n+1} = x^n + dt xdot (no interpolation between steps). D_x is
the second order three-point derivative on the non-uniform mesh; the matrix
(1/dt - Delta_h - xdot D_x + diag f) is an M-matrix as long as the cell Peclet
number satisfies |xdot| h <= 2 (checked), so the discrete comparison principle
of the Eulerian scheme still holds.

The band travels rigidly with the front, so inside the band the reaction layer
is (nearly) steady in the node frame and the time step is not limited by
eps. The mesh is rebuilt (on the lattice of the current band, hence exactly in
the band) only when the transition region has been distorted.
"""

import numpy as np
from scipy.linalg import solve_banded

from .model import f_eps
from .solver import MeshParams, Result, build_mesh, find_front, initial_data, laplacian


def _derivative_coeffs(x):
    """Coefficients (lo, di, up) of the 3-point second-order u_x at interior nodes."""
    hm = np.diff(x)[:-1]
    hp = np.diff(x)[1:]
    lo = -hp / (hm * (hm + hp))
    up = hm / (hp * (hm + hp))
    di = -(lo + up)
    return lo, di, up


def _matrix(Lap, x, xdot, react, dt):
    """Banded storage of (1/dt) I - Lap - diag(xdot) D_x + diag(react)."""
    n = x.size
    ab = np.zeros((3, n))
    ab[0, 1:] = -Lap.diagonal(1)
    ab[1, :] = 1.0 / dt - Lap.diagonal(0) + react
    ab[2, :-1] = -Lap.diagonal(-1)
    lo, di, up = _derivative_coeffs(x)
    i = np.arange(1, n - 1)
    ab[1, i] -= xdot[i] * di
    ab[0, i + 1] -= xdot[i] * up          # entry (i, i+1)
    ab[2, i - 1] -= xdot[i] * lo          # entry (i, i-1)
    if np.any(ab[0, 1:] > 1e-12) or np.any(ab[2, :-1] > 1e-12):
        raise RuntimeError("cell Peclet number too large: M-matrix property lost")
    return ab


def simulate_ale(u0, w0, eps, L, t_end, dt, mesh=None, speed_guess=1.0, ell0=None,
                 control_steps=10, recenter_time=None, remesh_distance=None,
                 snapshot_times=()):
    mp = mesh or MeshParams()
    # the rigidly moving zone must contain several decay lengths 1/c of the
    # profile ahead of the front; otherwise the end of the zone, where u is not
    # steady in the node frame, adds an O(dt) error that grows as eps -> 0
    ell0 = ell0 or max(2.0, 3.0 / max(speed_guess, 1e-6))
    remesh_distance = remesh_distance or 0.25 * ell0
    x0 = np.linspace(0.0, L, 20001)
    u, w = initial_data(x0, u0, w0, eps)
    x, xf = build_mesh(x0, u, eps, L, mp)
    u, w = initial_data(x, u0, w0, eps)
    y = u + w
    h_f = eps / mp.pts_per_eps

    # chi = 1 on a neighbourhood of the band (centred at the reference point xc),
    # decaying linearly to 0 over a distance ell0; nodes far away do not move
    R1 = (max(mp.K_burnt, mp.K_fresh) + mp.K_margin) * eps + 0.5 * ell0

    def chi(z):
        return np.clip(1.0 - (np.abs(z - xc) - R1) / ell0, 0.0, 1.0)

    res = Result(t=None, front=None)
    nsteps = int(round(t_end / dt))
    dt = t_end / nsteps
    recenter_time = recenter_time or 20 * dt
    ts, fronts = [0.0], [xf]
    v = speed_guess
    xc = xf                         # where we want the front (band reference point)
    moved = 0.0
    t_prev, xf_prev = 0.0, xf
    pending = sorted(snapshot_times)
    for n in range(nsteps):
        xdot = v * chi(x)
        xn = x + dt * xdot
        Lap = laplacian(xn, 1)
        F = f_eps(u, eps)
        y = solve_banded((1, 1), _matrix(Lap, xn, xdot, F, dt), y / dt)
        u = solve_banded((1, 1), _matrix(Lap, xn, xdot, np.zeros_like(F), dt), u / dt - F * y)
        x = xn
        xc += dt * v
        moved += dt * v
        t = (n + 1) * dt
        res.n_steps += 1
        if (n + 1) % control_steps == 0 or n == nsteps - 1:
            xf, _ = find_front(x, u, eps)
            ts.append(t)
            fronts.append(xf)
            # speed estimate + relaxation of the front towards the band centre
            v_meas = (xf - xf_prev) / (t - t_prev)
            v = v_meas + (xf - xc) / recenter_time
            t_prev, xf_prev = t, xf
            if pending and t >= pending[0] - 1e-12:
                res.snapshots[pending.pop(0)] = (x.copy(), u.copy(), (y - u).copy())
            if moved > remesh_distance:
                band_nodes = x[np.abs(np.diff(x, prepend=x[0] - 1)) < 1.5 * h_f]
                anchor = band_nodes[0] if band_nodes.size else 0.0
                x_new, _ = build_mesh(x, u, eps, L, mp, anchor=anchor)
                u, y = np.interp(x_new, x, u), np.interp(x_new, x, y)
                x = x_new
                xc = xf
                moved = 0.0
                res.n_remesh += 1
            res.n_nodes.append(x.size)
    res.t, res.front = np.array(ts), np.array(fronts, dtype=float)
    return res
