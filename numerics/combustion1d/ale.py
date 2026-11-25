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
method; each trial velocity costs two tridiagonal solves, counted in
res.n_trials / res.n_solves) so that the front ends the step at the band reference
point; this keeps the mesh-velocity error small also for accelerating fronts
(e.g. extinction), where an explicit estimate lags by (acceleration)*dt.
The mesh is rebuilt (on the lattice of the current band, hence exactly in
the band) when it has travelled remesh_distance, when the front has drifted
from the band reference point, or when a cell has been compressed below
min_cell * h_f.

Outcome of a run (Result.status): "ok", "extinct" (the front {u = eps}
disappeared: the run stops at that step) or "failed" (a numerical failure:
PecletError or MeshError is raised, with the partial Result attached to the
exception as .result).
"""

import numpy as np
from scipy.linalg import solve_banded

from .model import resolve_reaction
from .solver import (MeshError, MeshParams, NumericalFailure, PecletError, Result,
                     build_mesh, find_front, initial_data, laplacian)

SECANT_TOL = 1e-3       # tolerance of the implicit mesh velocity, in units of h_f
SECANT_MAXIT = 8        # secant updates per step (after the two starting trials)


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
    # Peclet condition: the off-diagonal entries must be <= 0 (M-matrix)
    if np.any(ab[0, 1:] > 1e-12) or np.any(ab[2, :-1] > 1e-12):
        raise PecletError("cell Peclet number too large: M-matrix property lost")
    return ab


def _band_anchor(x, h_f):
    """Anchor of the lattice of the current band: a node whose cell to the
    right has exactly the band spacing h_f (transition cells are >= growth *
    h_f, so they are excluded). A band built on anchor + h_f * Z shares its
    nodes with the current one, so the transfer by interpolation is exact
    there. 0.0 if no such cell exists (e.g. a band deformed near r = 0)."""
    exact = np.nonzero(np.abs(np.diff(x) - h_f) < 1e-6 * h_f)[0]
    return x[exact[0]] if exact.size else 0.0


def _band_cells(x, h_f):
    """(first, last) index of the cells of the band (the cells of spacing
    h_f) of a mesh that has just been built; a band that reaches an end of
    the domain is extended to that end. Cell j is [x[j], x[j+1]]."""
    cells = np.nonzero(np.abs(np.diff(x) - h_f) < 1e-6 * h_f)[0]
    if cells.size == 0:
        return 0, x.size - 2
    lo, hi = cells[0], cells[-1]
    return (0 if lo <= 1 else lo), (x.size - 2 if hi >= x.size - 3 else hi)


def simulate_ale(u0, w0, eps, L, t_end, dt, mesh=None, speed_guess=1.0, ell0=None,
                 recenter_time=None, remesh_distance=None,
                 snapshot_times=(), fixed_velocity=None, N=1, record_mesh=False,
                 implicit=True, rigid_extra=1.0, exact_initial=False, monitor=None,
                 remesh=True, min_cell=0.5, remesh_off_band=False, reaction=None):
    """ALE version of scheme (S).

    ell0: width of the zone where the mesh velocity decays linearly from v
        to 0. Default max(2, 3/|speed_guess|).
    rigid_extra: the mesh velocity is v (rigid zone) for
        |x - x_c| <= R1 = (max(K_burnt, K_fresh) + K_margin) eps + rigid_extra*ell0,
        measured symmetrically from the band reference point x_c, and decays
        linearly to 0 over the next ell0. With the default mesh the rigid
        zone therefore extends rigid_extra*ell0 beyond the band on the burnt
        side and (K_burnt - K_fresh) eps + rigid_extra*ell0 beyond it on the
        fresh side. Default 1.0. (With 0.5 the non-rigid zone leaves an
        O(dt e^{-c R1}) error of about 0.06 dt for w0 = 1 and the default
        ell0 = 2; with 1.0 about 0.011 dt.)
    remesh_distance: the mesh is rebuilt when it has travelled this distance
        (default ell0/4), or when |front - x_c| > K_margin eps / 4.
    min_cell: the mesh is also rebuilt when a cell has been compressed below
        min_cell * h_f (h_f = eps/pts_per_eps); these remeshings are counted
        in res.n_mesh_remesh. With remesh=False the mesh is not rebuilt. A
        trial step in which two nodes would cross raises MeshError.
        The run also counts the steps that end with the front outside the
        fine band (res.n_front_off_band, first time res.t_front_off_band).
        By default this is only a diagnostic. It happens where the mesh
        velocity is tapered (near r = 0 or x = L), because there the band
        nodes move more slowly than the front; once the front is in the
        coarse cells the reaction layer is not resolved and the front can
        stall.
    remesh_off_band: if True, the mesh is also rebuilt (around the front)
        after a step that ends with the front outside the fine band; these
        remeshings are counted in res.n_mesh_remesh too. Default False.
    fixed_velocity: prescribed mesh velocity, a number or a function of t
        (evaluated at t_n for the step n -> n+1); None = computed (implicit).
    implicit: with fixed_velocity=None, True (default) = the mesh velocity of
        each step is found by a secant iteration so that the front ends the
        step at x_c + dt v (tolerance SECANT_TOL * h_f, two starting trials
        plus at most SECANT_MAXIT updates). If the iteration ends without
        meeting the tolerance, the last trial with a front is accepted and
        the step is counted in res.n_secant_fail. False = explicit estimate
        from the measured front speed (recenter_time, default 20 dt).
    remesh: if False the mesh is never rebuilt (setting of Theorem cv:thm:main).
    exact_initial: if True, u(x, 0) = u0(x) and w = eps*w0 exactly (no
        modification of u0 near the front; use with the eps-travelling wave).
    monitor: optional callable monitor(t, x, u, v) called at t = 0 (v = None)
        and after every step, with the state carried to the next step (i.e.
        after remeshing, if any).
    reaction: None (the rate installed with model.set_reaction) or a
        callable f for this call only (see model.resolve_reaction).

    Returns a Result. res.status is "ok" or "extinct" (no crossing of u = eps
    after a step: the run stops there, res.front ends with nan, res.reason
    says when). A numerical failure (Peclet condition violated for every
    usable trial velocity, crossing nodes) raises PecletError / MeshError;
    the exception carries the partial Result (status "failed") as .result.
    """
    mp = MeshParams() if mesh is None else mesh
    # the rigidly moving zone must contain several decay lengths 1/c of the
    # profile ahead of the front; otherwise the end of the zone, where u is not
    # steady in the node frame, adds an O(dt) error that grows as eps -> 0
    if ell0 is None:
        ell0 = max(2.0, 3.0 / max(abs(speed_guess), 1e-6))
    if remesh_distance is None:
        remesh_distance = 0.25 * ell0
    if not dt > 0 or not ell0 > 0:
        raise ValueError("dt and ell0 must be positive")
    f_r = resolve_reaction(reaction)[0]

    def f_eps(s, e):
        return f_r(s / e) / e**2

    x0 = np.linspace(0.0, L, 20001)
    def init(xx):
        if exact_initial:
            return u0(xx), np.full_like(xx, eps * w0)
        return initial_data(xx, u0, w0, eps)
    u, w = init(x0)
    x, xf = build_mesh(x0, u, eps, L, mp)
    u, w = init(x)
    y = u + w
    # reference point = front on the actual mesh (build_mesh located it on the
    # auxiliary uniform grid x0, whose spacing L/20000 can exceed eps; using
    # that value made the first implicit step jump by (xf error)/dt)
    xf = find_front(x, u, eps)[0]
    if xf is None:
        raise ValueError("the initial datum has no front (no crossing of u = eps)")
    h_f = eps / mp.pts_per_eps

    # chi = 1 for |x - xc| <= R1 (xc = band reference point), then decays
    # linearly to 0 over a distance ell0; nodes far away do not move. R1 is
    # the same on both sides of xc: it covers the longer side of the band
    # (max(K_burnt, K_fresh) + K_margin) eps plus rigid_extra*ell0
    R1 = (max(mp.K_burnt, mp.K_fresh) + mp.K_margin) * eps + rigid_extra * ell0

    ell_end = min(0.5, 0.1 * L)       # taper near the ends (no node may cross r = 0)

    def chi(z):
        c = np.clip(1.0 - (np.abs(z - xc) - R1) / ell0, 0.0, 1.0)
        return c * np.clip(z / ell_end, 0.0, 1.0) * np.clip((L - z) / ell_end, 0.0, 1.0)

    res = Result(t=None, front=None)
    nsteps = int(round(t_end / dt))
    dt = t_end / nsteps
    if recenter_time is None:
        recenter_time = 20 * dt
    if not recenter_time > 0:
        raise ValueError("recenter_time must be positive")
    ts, fronts = [0.0], [xf]
    v_of_t = fixed_velocity if callable(fixed_velocity) else None
    v = speed_guess if fixed_velocity is None else (v_of_t(0.0) if v_of_t else fixed_velocity)
    xc = xf                         # where we want the front (band reference point)
    moved = 0.0
    t_prev, xf_prev = 0.0, xf
    pending = sorted(snapshot_times)
    if record_mesh:
        res.mesh_history.append((0.0, x.copy()))
    if monitor is not None:
        monitor(0.0, x, u, None)

    def advance(v):
        """One step of (S) on the mesh moving with velocity v*chi."""
        xdot = v * chi(x)
        xdot[0] = xdot[-1] = 0.0          # the ends of [0, L] stay fixed
        xn = x + dt * xdot
        if np.diff(xn).min() <= 0.0:
            raise MeshError("mesh nodes cross within one step (dt |v| |chi'| >= 1)")
        Lap = laplacian(xn, N)
        F = f_eps(u, eps)
        res.n_trials += 1
        res.n_solves += 2
        yn = solve_banded((1, 1), _matrix(Lap, xn, xdot, F, dt), y / dt)
        un = solve_banded((1, 1), _matrix(Lap, xn, xdot, np.zeros_like(F), dt), u / dt - F * yn)
        xfn, _ = find_front(xn, un, eps)
        return xn, un, yn, xfn

    tol = SECANT_TOL * h_f
    band = _band_cells(x, h_f)      # cells of the band (the nodes keep their indices)

    def step(n, v):
        """Advance one step; returns (v, x, u, y, xf) of the accepted state."""
        if implicit and fixed_velocity is None:
            # implicit mesh velocity: choose v so that, at the end of the step,
            # the front sits at the band reference point xc + dt v (secant method)
            def gap(vv):
                try:
                    out = advance(vv)
                except (PecletError, MeshError):   # unusable trial velocity
                    return np.inf, None
                return (np.inf if out[3] is None else out[3] - (xc + dt * vv)), out
            v0, (g0, out0) = v, gap(v)
            v1 = v0 + 0.05 * abs(v0) + 0.01
            g1, out1 = gap(v1)
            for _ in range(SECANT_MAXIT):
                if not np.isfinite(g0) or not np.isfinite(g1) or g1 == g0:
                    break
                v2 = v1 - g1 * (v1 - v0) / (g1 - g0)
                v2 = float(np.clip(v2, v1 - 2 * (abs(v1) + 1), v1 + 2 * (abs(v1) + 1)))
                v0, g0, out0 = v1, g1, out1
                v1 = v2
                g1, out1 = gap(v1)
                if abs(g1) < tol:
                    break
            if np.isfinite(g1):
                v, g, (x_, u_, y_, xf_) = v1, g1, out1
            elif out0 is not None:
                v, g, (x_, u_, y_, xf_) = v0, g0, out0
            else:                             # fall back to the previous velocity
                x_, u_, y_, xf_ = advance(v)  # (raises if it is not usable either)
                g = np.inf if xf_ is None else xf_ - (xc + dt * v)
            # accepted without meeting the tolerance (a state without a front
            # is reported as extinction by the caller, not counted here)
            if xf_ is not None and not abs(g) < tol:
                res.n_secant_fail += 1
            return v, x_, u_, y_, xf_
        if v_of_t is not None:
            v = v_of_t(n * dt)
        return (v,) + advance(v)

    def fail(n, msg, exc=None):
        """Close the partial result and raise the failure (never returns)."""
        res.status = "failed"
        res.reason = f"step {n + 1} (t = {(n + 1) * dt:.6g}): {msg}"
        res.t, res.front = np.array(ts), np.array(fronts, dtype=float)
        err = exc if exc is not None else NumericalFailure(res.reason)
        err.args = (res.reason,)
        err.result = res
        raise err

    for n in range(nsteps):
        try:
            v, x, u, y, xf = step(n, v)
        except NumericalFailure as ex:
            fail(n, f"{type(ex).__name__}: {ex}", ex)
        if xf is None:                    # no crossing of u = eps on the mesh
            if not (np.all(np.isfinite(u)) and np.all(np.isfinite(y))):
                fail(n, "non-finite values in the solution")
            ts.append((n + 1) * dt)
            fronts.append(np.nan)
            res.extinct = True            # the flame has gone out
            res.status = "extinct"
            res.reason = ("front disappeared: {} on the whole mesh after step {} (t = {:.6g})"
                          .format("u > eps" if u[0] > eps else "u < eps", n + 1, (n + 1) * dt))
            break
        res.v.append(v)
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
        # diagnostic: is the front still inside the fine band?
        off_band = not band[0] <= np.searchsorted(x, xf) - 1 <= band[1]
        if off_band:
            res.n_front_off_band += 1
            if res.t_front_off_band is None:
                res.t_front_off_band = t
        if pending and t >= pending[0] - 1e-12:
            res.snapshots[pending.pop(0)] = (x.copy(), u.copy(), (y - u).copy())
        # rebuild (exactly on the band lattice) when the mesh has travelled far,
        # when the front has drifted from the band reference point, or when the
        # mesh validation asks for it: a cell compressed below min_cell * h_f
        # (or, optionally, the front outside the band)
        squeezed = remesh and (np.diff(x).min() < min_cell * h_f
                               or (remesh_off_band and off_band))
        if remesh and (moved > remesh_distance or abs(xf - xc) > 0.25 * mp.K_margin * eps
                       or squeezed):
            x_new, _ = build_mesh(x, u, eps, L, mp, anchor=_band_anchor(x, h_f))
            u, y = np.interp(x_new, x, u), np.interp(x_new, x, y)
            x = x_new
            band = _band_cells(x, h_f)
            xc = xf
            moved = 0.0
            res.n_remesh += 1
            if squeezed:
                res.n_mesh_remesh += 1
        res.n_nodes.append(x.size)
        if record_mesh:
            res.mesh_history.append((t, x.copy()))
        if monitor is not None:
            monitor(t, x, u, v)
    res.t, res.front = np.array(ts), np.array(fronts, dtype=float)
    return res
