"""Monotone scheme (S) in ALE form: the mesh follows the flame (planar or radial).

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
eps. By default the mesh velocity is IMPLICIT: at each step v is chosen (secant
method, 2-4 solves) so that the front ends the step at the band reference
point; this keeps the mesh-velocity error small also for accelerating fronts
(e.g. extinction), where an explicit estimate lags by (acceleration)*dt. The mesh is rebuilt (on the lattice of the current band, hence exactly in
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
                 recenter_time=None, remesh_distance=None,
                 snapshot_times=(), fixed_velocity=None, N=1, record_mesh=False,
                 implicit=True):
    mp = mesh or MeshParams()
    # the rigidly moving zone must contain several decay lengths 1/c of the
    # profile ahead of the front; otherwise the end of the zone, where u is not
    # steady in the node frame, adds an O(dt) error that grows as eps -> 0
    ell0 = ell0 or max(2.0, 3.0 / max(abs(speed_guess), 1e-6))
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

    ell_end = min(0.5, 0.1 * L)       # taper near the ends (no node may cross r = 0)

    def chi(z):
        c = np.clip(1.0 - (np.abs(z - xc) - R1) / ell0, 0.0, 1.0)
        return c * np.clip(z / ell_end, 0.0, 1.0) * np.clip((L - z) / ell_end, 0.0, 1.0)

    res = Result(t=None, front=None)
    nsteps = int(round(t_end / dt))
    dt = t_end / nsteps
    recenter_time = recenter_time or 20 * dt
    ts, fronts = [0.0], [xf]
    v = speed_guess if fixed_velocity is None else fixed_velocity
    xc = xf                         # where we want the front (band reference point)
    moved = 0.0
    t_prev, xf_prev = 0.0, xf
    pending = sorted(snapshot_times)
    if record_mesh:
        res.mesh_history.append((0.0, x.copy()))
    def advance(v):
        """One step of (S) on the mesh moving with velocity v*chi."""
        xdot = v * chi(x)
        xdot[0] = xdot[-1] = 0.0          # the ends of [0, L] stay fixed
        xn = x + dt * xdot
        Lap = laplacian(xn, N)
        F = f_eps(u, eps)
        yn = solve_banded((1, 1), _matrix(Lap, xn, xdot, F, dt), y / dt)
        un = solve_banded((1, 1), _matrix(Lap, xn, xdot, np.zeros_like(F), dt), u / dt - F * yn)
        xfn, _ = find_front(xn, un, eps)
        return xn, un, yn, xfn

    for n in range(nsteps):
        if implicit and fixed_velocity is None:
            # implicit mesh velocity: choose v so that, at the end of the step,
            # the front sits at the band reference point xc + dt v (secant method)
            def gap(vv):
                try:
                    out = advance(vv)
                except RuntimeError:          # trial velocity violates the Peclet bound
                    return np.inf, None
                return (np.inf if out[3] is None else out[3] - (xc + dt * vv)), out
            v0, (g0, out0) = v, gap(v)
            v1 = v0 + 0.05 * abs(v0) + 0.01
            g1, out1 = gap(v1)
            for _ in range(8):
                if not np.isfinite(g0) or not np.isfinite(g1) or g1 == g0:
                    break
                v2 = v1 - g1 * (v1 - v0) / (g1 - g0)
                v2 = float(np.clip(v2, v1 - 2 * (abs(v1) + 1), v1 + 2 * (abs(v1) + 1)))
                v0, g0, out0 = v1, g1, out1
                v1 = v2
                g1, out1 = gap(v1)
                if abs(g1) < 1e-3 * h_f:
                    break
            if np.isfinite(g1):
                v, (x, u, y, xf) = v1, out1
            elif out0 is not None:
                v, (x, u, y, xf) = v0, out0
            else:                             # fall back to the previous velocity
                x, u, y, xf = advance(v)
        else:
            x, u, y, xf = advance(v)
        if xf is None:                    # the flame has gone out
            ts.append((n + 1) * dt)
            fronts.append(np.nan)
            res.extinct = True
            break
        xc += dt * v
        moved += dt * abs(v)
        t = (n + 1) * dt
        res.n_steps += 1
        ts.append(t)
        fronts.append(xf)
        if not implicit and fixed_velocity is None:
            v_meas = (xf - xf_prev) / (t - t_prev)
            v = v_meas + (xf - xc) / recenter_time
        t_prev, xf_prev = t, xf
        if pending and t >= pending[0] - 1e-12:
            res.snapshots[pending.pop(0)] = (x.copy(), u.copy(), (y - u).copy())
        # rebuild (exactly on the band lattice) when the mesh has travelled far
        # or when the front has drifted from the band reference point
        if moved > remesh_distance or abs(xf - xc) > 0.25 * mp.K_margin * eps:
            # anchor the new band on the lattice of the current one: a node whose
            # cell to the right has exactly the band spacing h_f (transition cells
            # are >= growth * h_f, so they are excluded)
            exact = np.nonzero(np.abs(np.diff(x) - h_f) < 1e-6 * h_f)[0]
            anchor = x[exact[0]] if exact.size else 0.0
            x_new, _ = build_mesh(x, u, eps, L, mp, anchor=anchor)
            u, y = np.interp(x_new, x, u), np.interp(x_new, x, y)
            x = x_new
            xc = xf
            moved = 0.0
            res.n_remesh += 1
        res.n_nodes.append(x.size)
        if record_mesh:
            res.mesh_history.append((t, x.copy()))
    res.t, res.front = np.array(ts), np.array(fronts, dtype=float)
    return res
