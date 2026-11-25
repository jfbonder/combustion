import numpy as np
import pytest

from combustion1d.exact import (SelfSimilar, tw_profile_limit, tw_speed_eps,
                                tw_speed_limit)
from combustion1d.model import M, fb_slope
from combustion1d.solver import (MeshParams, build_mesh, front_speed,
                                 laplacian, simulate)


def test_M_closed_form():
    from scipy.integrate import quad
    for w0 in (-0.5, 0.0, 1.0):
        num = quad(lambda s: (s + w0) * max(1 - s, 0), -w0, 1)[0]
        assert M(w0) == pytest.approx(num, rel=1e-12)


@pytest.mark.parametrize("N", [1, 2, 3])
def test_laplacian_exact_on_quadratics(N):
    x = np.linspace(0, 1, 301) ** 1.7          # non-uniform mesh
    v = laplacian(x, N) @ x**2
    assert np.allclose(v[:-1], 2 * N, atol=1e-8)   # last node: zero-flux BC


@pytest.mark.parametrize("N", [1, 2, 3])
def test_self_similar_is_exact_solution(N):
    V = SelfSimilar(N, w0=0.5, T=1.0)
    t, h, k = 0.3, 1e-4, 1e-5
    for r in np.linspace(0.2, 0.9, 5) * V.radius(t):
        lap = ((V(r + h, t) - 2 * V(r, t) + V(r - h, t)) / h**2
               + (N - 1) / r * (V(r + h, t) - V(r - h, t)) / (2 * h))
        vt = (V(r, t + k) - V(r, t - k)) / (2 * k)
        assert abs(lap - vt) < 1e-4
    rf = V.radius(t)
    grad = (V(rf - 2 * h, t) - V(rf - h, t)) / h      # one-sided, inside
    assert grad == pytest.approx(fb_slope(0.5), rel=1e-3)


def test_travelling_wave_speed_converges_linearly():
    for w0 in (-0.5, 0.0, 1.0):
        c0 = tw_speed_limit(1.0, w0)
        errs = [abs(tw_speed_eps(1.0, w0, e) - c0) for e in (0.02, 0.01, 0.005)]
        assert errs[0] > errs[1] > errs[2]
        assert errs[1] / errs[2] == pytest.approx(2.0, rel=0.15)


def test_mesh_refines_band_and_coarsens_elsewhere():
    eps, L = 0.01, 10.0
    x = np.linspace(0, L, 5001)
    u = tw_profile_limit(x, 0, 1.0, 0.5, x0=4.0)
    mp = MeshParams()
    xn, xf = build_mesh(x, u, eps, L, mp)
    h = np.diff(xn)
    assert xn[0] == 0 and xn[-1] == L and np.all(h > 0)
    assert h.min() == pytest.approx(eps / mp.pts_per_eps, rel=0.05)
    assert h.max() <= mp.h_coarse * 1.5
    assert np.max(h[1:] / h[:-1]) < 1.6 and np.min(h[1:] / h[:-1]) > 1 / 1.6
    near = np.abs(xn - xf) < 5 * eps
    assert np.all(np.diff(xn[near]) < 1.01 * eps / mp.pts_per_eps)


def test_pde_speed_matches_travelling_wave():
    A, w0, eps = 1.0, 0.0, 0.05
    c0 = tw_speed_limit(A, w0)
    u0 = lambda x: tw_profile_limit(x, 0.0, A, c0, x0=2.0)
    T = 6.0 / c0
    r = simulate(u0, w0, eps, 2 + 6 + 12 / c0, T, speed_guess=c0)
    assert front_speed(r, 0.5 * T) == pytest.approx(tw_speed_eps(A, w0, eps), abs=2e-3)


# --- limit problem solver ---------------------------------------------------

from combustion1d.limit import solve_planar, solve_radial


@pytest.mark.parametrize("N", [2, 3])
def test_limit_radial_matches_self_similar_second_order(N):
    V = SelfSimilar(N, w0=0.5, T=1.0)
    errs = []
    for n in (100, 200):
        r = solve_radial(lambda x: V(x, 0.0), 0.5, N, V.R, 0.5, n=n,
                         t_eval=np.linspace(0, 0.5, 11))
        errs.append(np.max(np.abs(r.s - V.radius(r.t))))
    assert errs[1] < 1e-4
    assert errs[0] / errs[1] == pytest.approx(4.0, rel=0.3)


@pytest.mark.parametrize("w0", [-0.5, 0.0, 1.0])
def test_limit_planar_matches_travelling_wave(w0):
    c0 = tw_speed_limit(1.0, w0)
    L = 2 + 8 + 12 / c0 + 2
    r = solve_planar(lambda x: tw_profile_limit(x, 0, 1.0, c0, x0=2.0), w0, L, 2.0,
                     8 / c0, n=400, t_eval=np.linspace(0, 8 / c0, 20))
    assert np.max(np.abs(r.s - (2.0 + c0 * r.t))) < 5e-5


@pytest.mark.parametrize("w0", [0.0, 1.0])
def test_tw_profile_eps_solves_the_ode(w0):
    from combustion1d.exact import tw_profile_eps
    from combustion1d.model import f_eps
    A, eps = 1.0, 0.01
    c = tw_speed_eps(A, w0, eps)
    P = lambda z: tw_profile_eps(z, A, w0, eps)
    assert P(np.array([0.0]))[0] == pytest.approx(eps, rel=1e-12)
    assert P(np.array([-60 * eps]))[0] == pytest.approx(-eps * w0, abs=1e-8 * eps)
    assert np.all(np.diff(P(np.linspace(-0.5, 2, 2001))) >= 0)
    # residual of U'' + c U' = (U + eps w0) f_eps(U), by centred differences
    h = 1e-4 * eps
    for z in np.linspace(-30 * eps, 3 * eps, 23):
        if abs(z) < 2 * h:
            continue
        p = P(np.array([z - h, z, z + h]))
        res = (p[2] - 2 * p[1] + p[0]) / h**2 + c * (p[2] - p[0]) / (2 * h) \
            - (p[1] + eps * w0) * f_eps(p[1], eps)
        assert abs(res) < 1e-3 / eps   # terms are O(1/eps)
    # C^1 matching at z = 0
    d = 1e-7
    gl = (P(np.array([0.0]))[0] - P(np.array([-d]))[0]) / d
    gr = (P(np.array([d]))[0] - P(np.array([0.0]))[0]) / d
    assert gl == pytest.approx(gr, rel=1e-4)


def test_set_reaction_custom_f():
    """A user-supplied f: M by quadrature, travelling-wave speed -> sqrt(2M)/A,
    and restoring the default."""
    from combustion1d import model
    from combustion1d.exact import tw_speed_eps, tw_speed_limit

    try:
        # same f as the default, given as a plain function: M by quadrature
        model.set_reaction(lambda s: np.maximum(1.0 - s, 0.0))
        assert abs(model.M(0.5) - 1.5**3 / 6) < 1e-10
        # a different rate: f(s) = ((1 - s)^+)^2, M(0) = int_0^1 s (1-s)^2 = 1/12
        model.set_reaction(lambda s: np.maximum(1.0 - s, 0.0) ** 2,
                           lambda s: np.where(s < 1.0, -2.0 * (1.0 - s), 0.0))
        assert abs(model.M(0.0) - 1.0 / 12.0) < 1e-10
        c0 = tw_speed_limit(1.0, 0.0)
        assert abs(c0 - np.sqrt(1.0 / 6.0)) < 1e-10
        c = tw_speed_eps(1.0, 0.0, 1e-3)
        assert abs(c - c0) < 5e-3
    finally:
        model.set_reaction()
    assert abs(model.M(0.0) - 1.0 / 6.0) < 1e-14


def test_ale_travelling_wave_from_exact_profile():
    """ALE scheme started from Phi_eps: profile and speed errors are small
    and do not change with the time step (Theorem an:thm:ale)."""
    from combustion1d.ale import simulate_ale
    from combustion1d.exact import tw_profile_eps, tw_speed_eps

    A, w0, eps, T, x0 = 1.0, 0.0, 0.01, 2.0, 2.0
    c = tw_speed_eps(A, w0, eps)
    L = x0 + c * T + 12.0 / c
    errs = []
    for dt in (0.05, 0.01):
        res = simulate_ale(lambda x: tw_profile_eps(x - x0, A, w0, eps), w0, eps, L, T, dt,
                           speed_guess=c, exact_initial=True, snapshot_times=(T,))
        x, u, _ = res.snapshots[T]
        errs.append(np.max(np.abs(u - tw_profile_eps(x - x0 - c * T, A, w0, eps))))
        assert abs(res.front[-1] - (x0 + c * T)) < 2e-3
    assert max(errs) < 2e-3
    assert abs(errs[0] - errs[1]) < 0.2 * max(errs)


# --- robustness of the front detector and of the mesh generator -------------

def test_find_front_node_exactly_on_the_level():
    from combustion1d.solver import find_front
    x = np.linspace(0.0, 1.0, 11)
    # the node x = 0.5 has u = level exactly
    xf, side = find_front(x, x.copy(), 0.5)
    assert xf == pytest.approx(0.5, abs=1e-15) and side == -1
    xf, side = find_front(x, 1.0 - x, 0.5)
    assert xf == pytest.approx(0.5, abs=1e-15) and side == +1
    # two consecutive nodes on the level: the first one is the front
    u = np.array([0.0, 0.2, 0.5, 0.5, 0.8, 1.0])
    xf, side = find_front(np.arange(6.0), u, 0.5)
    assert xf == 2.0 and side == -1
    # generic crossing: linear interpolation, as before
    xf, side = find_front(x, x**2, 0.3)
    assert 0.5 < xf < 0.6 and side == -1
    # touching the level without crossing, level at an end node, no crossing
    assert find_front(x, np.abs(x - 0.5) + 0.5, 0.5) == (None, 0)
    assert find_front(x, x.copy(), 0.0) == (None, 0)
    assert find_front(x, x + 1.0, 0.5) == (None, 0)
    assert find_front(x, np.full_like(x, np.nan), 0.5) == (None, 0)


def test_mesh_has_no_tiny_cells_when_the_band_touches_an_end():
    """Random lattice anchors, front close to x = 0 and to x = L, both
    orientations: every cell is >= h_f/2 and neighbouring cells are graded."""
    rng = np.random.default_rng(12345)
    mp = MeshParams()
    L = 3.0
    xx = np.linspace(0.0, L, 4001)
    for eps in (0.01, 0.003):
        h_f = eps / mp.pts_per_eps
        special = (1e-13, 0.3 * h_f, 0.49 * h_f, L - 1e-13, L - 0.3 * h_f)
        for side in (-1.0, 1.0):
            for near_zero in (True, False):
                for k in range(40):
                    d = rng.uniform(0.5, 45.0) * eps
                    xf = d if near_zero else L - d
                    u = eps + side * (xx - xf)          # one crossing of u = eps, at xf
                    anchor = special[k] if k < len(special) else rng.uniform(0.0, L)
                    x, xf_found = build_mesh(xx, u, eps, L, mp, anchor=anchor)
                    h = np.diff(x)
                    assert x[0] == 0.0 and x[-1] == L
                    assert xf_found == pytest.approx(xf, abs=1e-9)
                    assert h.min() >= 0.5 * h_f * (1 - 1e-9), (eps, side, d / eps, anchor)
                    ratio = h[1:] / h[:-1]
                    assert ratio.max() < 2.2 and ratio.min() > 1 / 2.2
                    # the band keeps the lattice anchor + h_f Z around the front
                    band = x[1:-1][np.abs(x[1:-1] - xf) < 3 * eps]
                    off = (band - anchor) / h_f
                    assert np.max(np.abs(off - np.round(off))) < 1e-6


def test_remeshing_transfer_is_exact_on_the_band():
    """A band rebuilt on the lattice of the previous one (as simulate_ale
    does) shares its nodes with it, so linear interpolation reproduces the
    nodal values there exactly; interpolation acts only outside."""
    from combustion1d.ale import _band_anchor
    eps, L = 0.01, 10.0
    mp = MeshParams()
    h_f = eps / mp.pts_per_eps
    prof = lambda x, xf: tw_profile_limit(x, 0.0, 1.0, 0.6, x0=xf) + eps * (x > xf) * 0 \
        - eps * np.exp(-np.abs(x - xf) / eps) * (x <= xf)
    xf0 = 4.0 + 0.123 * h_f
    grid = np.linspace(0.0, L, 20001)
    x1, _ = build_mesh(grid, prof(grid, xf0) + eps, eps, L, mp, anchor=0.37 * h_f)
    for shift in (3.3 * h_f, 7.9 * eps, -2.6 * eps):
        u1 = prof(x1, xf0 + shift) + eps         # the front has moved by `shift`
        x2, xf2 = build_mesh(x1, u1, eps, L, mp, anchor=_band_anchor(x1, h_f))
        u2 = np.interp(x2, x1, u1)
        h2 = np.diff(x2)
        in_band = np.nonzero(np.abs(h2 - h_f) < 1e-6 * h_f)[0]
        band2 = x2[in_band[0]:in_band[-1] + 2]
        # nodes of the new band that lie inside the old band
        old = np.nonzero(np.abs(np.diff(x1) - h_f) < 1e-6 * h_f)[0]
        lo, hi = x1[old[0]], x1[old[-1] + 1]
        common = band2[(band2 >= lo - 1e-12) & (band2 <= hi + 1e-12)]
        assert common.size > 20 * mp.pts_per_eps      # the bands overlap widely
        j = np.searchsorted(x1, common - 1e-12)
        assert np.max(np.abs(x1[j] - common)) < 1e-12 * L
        exact = prof(common, xf0 + shift) + eps
        assert np.max(np.abs(np.interp(common, x2, u2) - exact)) < 1e-13
        # the reaction zone is inside the common part
        assert lo < xf2 - 5 * eps and xf2 + 5 * eps < hi


# --- the Eulerian scheme (S), simple.py --------------------------------------

def test_simple_scheme_error_is_first_order_in_dt_over_eps():
    """Eulerian scheme (S) from the exact wave: nominal step used exactly, two
    tridiagonal solves per step, and a speed error proportional to dt/eps."""
    from combustion1d.exact import tw_profile_eps
    from combustion1d.simple import simulate_simple
    A, w0, eps, x0, T = 1.0, 0.0, 0.05, 2.0, 3.0
    c = tw_speed_eps(A, w0, eps)
    L = x0 + c * T + 12.0 / c
    rel = []
    for k in (1.0, 0.25):
        dt = k * eps**2
        r = simulate_simple(lambda x: tw_profile_eps(x - x0, A, w0, eps), w0, eps, L, T,
                            dt=dt, speed_guess=c, exact_initial=True)
        assert r.status == "ok" and r.dts == [pytest.approx(dt)]
        assert r.n_solves == 2 * r.n_steps and r.n_steps == round(T / dt)
        assert r.front[0] == pytest.approx(x0, abs=1e-12)     # exact front of the datum
        rel.append((front_speed(r, 0.5 * T) - c) / c)
    assert rel[0] < rel[1] < 0                      # the scheme is slower than the wave
    assert rel[0] / rel[1] == pytest.approx(4.0, rel=0.1)
    assert rel[1] / (0.25 * eps) == pytest.approx(-0.5, abs=0.1)


def test_initial_front_is_located_on_the_actual_mesh():
    """res.front[0] of the Eulerian solvers is the front of the initial datum
    on the mesh of the run, not on the auxiliary grid linspace(0, L, 20001)
    (whose spacing exceeds eps here)."""
    from combustion1d.exact import tw_profile_eps
    from combustion1d.simple import simulate_simple
    A, eps, L = 1.0, 1e-3, 27.0
    assert L / 20000 > eps
    x0 = 2.0 + 0.37 * eps / 16
    c = tw_speed_eps(A, 0.0, eps)
    r = simulate_simple(lambda x: tw_profile_eps(x - x0, A, 0.0, eps), 0.0, eps, L, 2e-6,
                        dt=1e-6, speed_guess=c, exact_initial=True)
    assert abs(r.front[0] - x0) < 1e-2 * eps
    # BDF solver, data of the paper with w0 = 1: u0^eps = 2 u0 - eps below eps,
    # so the front is where A (1 - exp(-c (x - 2))) = eps
    r = simulate(lambda x: tw_profile_limit(x, 0.0, A, c, x0=2.0), 1.0, eps, L, 1e-7,
                 speed_guess=c)
    assert abs(r.front[0] - (2.0 - np.log(1.0 - eps / A) / c)) < 1e-2 * eps


# --- Lemma lem:bounds and Lemma lem:comparison of the paper -------------------

def _check_bounds(rec, snapshots, eps, w0):
    """Y >= 0, W = Y - U = eps w0 exactly, max U nonincreasing."""
    lo = np.array([a for a, _ in rec])
    hi = np.array([b for _, b in rec])
    assert lo.min() + eps * w0 >= -1e-11          # Y = U + W >= 0
    assert np.max(np.diff(hi)) <= 1e-11           # max U^{n+1} <= max U^n
    assert snapshots
    for x, u, w in snapshots.values():
        assert np.max(np.abs(w - eps * w0)) < 1e-11
        assert np.min(u + w) >= -1e-11


def test_bounds_of_the_scheme_eulerian():
    from combustion1d.simple import simulate_simple
    eps, w0 = 0.02, 0.5
    c0 = tw_speed_limit(1.0, w0)
    rec = []
    r = simulate_simple(lambda x: tw_profile_limit(x, 0.0, 1.0, c0, x0=2.0), w0, eps, 12.0, 2.0,
                        dt=0.02, speed_guess=c0, snapshot_times=(0.5, 1.0, 2.0),
                        monitor=lambda t, x, u, v: rec.append((u.min(), u.max())))
    assert r.n_remesh > 3 and r.status == "ok"     # the bounds survive the remeshings
    _check_bounds(rec, r.snapshots, eps, w0)


@pytest.mark.parametrize("dt", [0.05, 0.01])
def test_bounds_of_the_scheme_ale_radial(dt):
    """Radial extinction with large steps (dt = 5 eps and eps): the bounds of
    Lemma lem:bounds hold on the moving mesh and across remeshings, and the
    front follows the exact self-similar solution to O(dt)."""
    from combustion1d.ale import simulate_ale
    N, eps, w0, t_end = 2, 0.01, 0.5, 0.75
    V = SelfSimilar(N, w0, 1.0)
    L = 1.5 * V.R + 1.0
    rec = []
    r = simulate_ale(lambda x: V(x, 0.0), w0, eps, L, t_end, dt, speed_guess=-V.R / 2, N=N,
                     ell0=1.0, remesh_distance=0.1, snapshot_times=(0.25, 0.5, t_end),
                     monitor=lambda t, x, u, v: rec.append((u.min(), u.max())))
    assert r.status == "ok" and not r.extinct and r.n_remesh > 0
    assert r.n_mesh_remesh == 0 and r.n_front_off_band == 0 and r.n_solves == 2 * r.n_trials
    _check_bounds(rec, r.snapshots, eps, w0)
    m = r.t >= 0.05
    err = np.abs(r.front[m] + eps / fb_slope(w0) - V.radius(r.t[m]))
    assert err.max() < 2.0 * dt
    if dt == 0.01:
        assert r.n_secant_fail == 0 and r.n_trials <= 3.3 * r.n_steps


def _ordered_data(w0):
    lam = float(fb_slope(w0))
    base = lambda x: np.where(x > 2.0, np.tanh(lam * (x - 2.0)), 0.0)
    # compactly supported bump away from the front: same front, same meshes
    bump = lambda x: 0.3 * np.maximum(0.0, 1.0 - ((x - 4.0) / 0.8) ** 2) ** 2
    return base, (lambda x: base(x) + bump(x))


def _compare(runs):
    (a, b) = runs
    assert len(a) == len(b) > 5
    worst = 0.0
    for (ta, xa, ua), (tb, xb, ub) in zip(a, b):
        assert ta == tb and np.array_equal(xa, xb)       # same sequence of meshes
        worst = min(worst, float(np.min(ub - ua)))
    assert worst > -1e-10
    assert np.max(b[-1][2] - a[-1][2]) > 0.05            # the data really differ


def test_discrete_comparison_principle_eulerian():
    """Two ordered initial data on the same mesh give ordered solutions, for
    a time step 1000 times larger than eps^2 (Lemma lem:comparison)."""
    from combustion1d.simple import simulate_simple
    eps, w0, L, T, dt = 0.01, 0.5, 6.0, 2.0, 0.1
    runs = []
    for u0 in _ordered_data(w0):
        rec = []
        r = simulate_simple(u0, w0, eps, L, T, dt=dt, speed_guess=1e-12, max_dt_remesh=T,
                            monitor=lambda t, x, u, v: rec.append((t, x.copy(), u.copy())))
        assert r.n_steps == 20 and r.n_remesh == 1       # one mesh for the whole run
        runs.append(rec)
    _compare(runs)


@pytest.mark.parametrize("N", [1, 2])
def test_discrete_comparison_principle_moving_mesh(N):
    """The same on a mesh that moves with a prescribed velocity (ALE form of
    the scheme, no remeshing), dt = 5 eps."""
    from combustion1d.ale import simulate_ale
    eps, w0, L = 0.01, 0.5, 6.0
    runs = []
    for u0 in _ordered_data(w0):
        rec = []
        r = simulate_ale(u0, w0, eps, L, 1.0, 0.05, speed_guess=1.0, fixed_velocity=1.0,
                         remesh=False, ell0=1.0, N=N,
                         monitor=lambda t, x, u, v: rec.append((t, x.copy(), u.copy())))
        assert r.status == "ok" and r.n_remesh == 0
        runs.append(rec)
    _compare(runs)
    assert not np.array_equal(runs[0][0][1], runs[0][-1][1])     # the mesh did move


# --- failures of the ALE scheme ----------------------------------------------

def test_peclet_error_is_raised():
    from combustion1d.ale import _matrix, simulate_ale
    from combustion1d.solver import NumericalFailure, PecletError
    assert issubclass(PecletError, NumericalFailure) and issubclass(NumericalFailure, RuntimeError)
    # planar, uniform mesh h = 0.1: the condition is |xdot| h <= 2
    x = np.linspace(0.0, 1.0, 11)
    Lap = laplacian(x, 1)
    for v in (19.0, -19.0, 21.0, -21.0):
        xdot = np.full(x.size, v)
        xdot[0] = xdot[-1] = 0.0
        if abs(v) * 0.1 <= 2.0:
            ab = _matrix(Lap, x, xdot, np.zeros(x.size), 0.1)
            assert np.all(ab[0, 1:] <= 0.0) and np.all(ab[2, :-1] <= 0.0)
        else:
            with pytest.raises(PecletError, match="Peclet"):
                _matrix(Lap, x, xdot, np.zeros(x.size), 0.1)
    # a run with a prescribed mesh velocity that violates it (h_c = 0.05, v = 60)
    c0 = tw_speed_limit(1.0, 0.0)
    with pytest.raises(PecletError) as info:
        simulate_ale(lambda x: tw_profile_limit(x, 0.0, 1.0, c0, x0=2.0), 0.0, 0.02, 12.0, 0.1,
                     0.01, mesh=MeshParams(h_coarse=0.05), speed_guess=c0, fixed_velocity=60.0)
    res = info.value.result
    assert res.status == "failed" and "PecletError" in res.reason and "step 1 " in res.reason
    assert not res.extinct and res.n_steps == 0 and res.t.size == res.front.size == 1
    with pytest.raises(RuntimeError):                   # still a RuntimeError for old callers
        simulate_ale(lambda x: tw_profile_limit(x, 0.0, 1.0, c0, x0=2.0), 0.0, 0.02, 12.0, 0.1,
                     0.01, mesh=MeshParams(h_coarse=0.05), speed_guess=c0, fixed_velocity=60.0)


def test_extinction_is_reported_as_status_not_as_failure():
    """A burnt kernel that is too small collapses: the run stops with status
    'extinct' (a physical outcome), without raising."""
    from combustion1d.ale import simulate_ale
    eps, w0, r0 = 0.02, 0.5, 0.4
    lam = float(fb_slope(w0))
    u0 = lambda r: np.where(r > r0, np.tanh(lam * (r - r0)), 0.0)
    r = simulate_ale(u0, w0, eps, 12.0, 6.0, 0.002, speed_guess=lam, N=2)
    assert r.status == "extinct" and r.extinct and "front disappeared" in r.reason
    assert np.isnan(r.front[-1]) and np.all(np.isfinite(r.front[:-1]))
    assert r.t[-1] < 1.0 and r.front[-2] < 0.25 * r0
    assert r.n_steps == r.t.size - 2


def test_mesh_validation_triggers_a_remesh():
    """Rigid translation of the whole mesh: the cells near the fixed end
    x = L are compressed. With the default min_cell the mesh is rebuilt
    before a cell falls below h_f/2; with min_cell = 0 it is left alone."""
    from combustion1d.ale import simulate_ale
    from combustion1d.exact import tw_profile_eps
    A, eps, x0, T = 1.0, 0.02, 2.0, 3.0
    h_f = eps / MeshParams().pts_per_eps
    c = tw_speed_eps(A, 0.0, eps)
    L = x0 + c * T + 12.0 / c
    kw = dict(speed_guess=c, exact_initial=True, fixed_velocity=c, ell0=1e6, remesh_distance=1e9,
              record_mesh=True)
    r = simulate_ale(lambda x: tw_profile_eps(x - x0, A, 0.0, eps), 0.0, eps, L, T, 0.05, **kw)
    assert r.status == "ok" and r.n_mesh_remesh >= 1 and r.n_remesh == r.n_mesh_remesh
    assert min(np.diff(x).min() for _, x in r.mesh_history) >= 0.5 * h_f
    r = simulate_ale(lambda x: tw_profile_eps(x - x0, A, 0.0, eps), 0.0, eps, L, T, 0.05,
                     min_cell=0.0, **kw)
    assert r.status == "ok" and r.n_remesh == 0 and r.n_mesh_remesh == 0
    assert min(np.diff(x).min() for _, x in r.mesh_history) < 0.5 * h_f


def test_front_leaving_the_band_is_detected():
    """A collapsing kernel with a narrow band (small eps) in the zone near
    r = 0 where the mesh velocity is tapered: the front runs ahead of the band
    nodes. The run records it (n_front_off_band); with remesh_off_band=True
    the mesh is rebuilt around the front instead."""
    from combustion1d.ale import simulate_ale
    eps, w0, r0 = 0.004, 0.5, 0.4
    lam = float(fb_slope(w0))
    u0 = lambda r: np.where(r > r0, np.tanh(lam * (r - r0)), 0.0)
    a = simulate_ale(u0, w0, eps, 12.0, 2.0, 0.004, speed_guess=lam, N=2)
    assert a.status == "extinct" and a.n_front_off_band > 10
    assert 0.0 < a.t_front_off_band < a.t[-1]
    b = simulate_ale(u0, w0, eps, 12.0, 2.0, 0.004, speed_guess=lam, N=2, remesh_off_band=True)
    assert b.status == "extinct" and b.n_front_off_band <= 3
    assert b.n_mesh_remesh > a.n_mesh_remesh
    assert b.t_front_off_band == a.t_front_off_band       # identical until then
    # a front that stays in the band is never flagged
    c0 = tw_speed_limit(1.0, 0.0)
    c = simulate_ale(lambda x: tw_profile_limit(x, 0.0, 1.0, c0, x0=2.0), 0.0, 0.02, 12.0, 1.0,
                     0.02, speed_guess=c0)
    assert c.status == "ok" and c.n_front_off_band == 0 and c.t_front_off_band is None


def test_explicit_zero_arguments_are_not_replaced_by_defaults():
    from combustion1d.ale import simulate_ale
    from combustion1d.simple import simulate_simple
    u0 = lambda x: tw_profile_limit(x, 0.0, 1.0, 0.5, x0=2.0)
    with pytest.raises(ValueError):
        simulate_simple(u0, 0.0, 0.05, 8.0, 0.1, dt=0.0)
    with pytest.raises(ValueError):
        simulate_simple(u0, 0.0, 0.05, 8.0, 0.1, max_dt_remesh=0.0)
    with pytest.raises(ValueError):
        simulate(u0, 0.0, 0.05, 8.0, 0.1, max_dt_remesh=0.0)
    with pytest.raises(ValueError):
        simulate_ale(u0, 0.0, 0.05, 8.0, 0.1, 0.01, ell0=0.0)
    # remesh_distance = 0 now means "rebuild after every step"
    r = simulate_ale(u0, 0.0, 0.05, 8.0, 0.1, 0.01, speed_guess=0.5, remesh_distance=0.0)
    assert r.n_remesh == r.n_steps == 10


# --- the reaction rate: per-call argument and worker processes ---------------

def _f_quadratic(s):
    return np.maximum(1.0 - s, 0.0) ** 2


def _df_quadratic(s):
    return np.where(s < 1.0, -2.0 * (1.0 - s), 0.0)


def _M_in_worker(w0):
    from combustion1d import model
    return float(model.M(w0))


def test_reaction_argument_matches_set_reaction():
    from combustion1d import model
    from combustion1d.ale import simulate_ale
    from combustion1d.simple import simulate_simple
    u0 = lambda x: tw_profile_limit(x, 0.0, 1.0, 0.4, x0=2.0)
    args = (u0, 0.0, 0.05, 8.0, 0.2)
    a = simulate_simple(*args, dt=0.01, reaction=_f_quadratic)
    b = simulate_ale(*args, 0.01, speed_guess=0.4, reaction=_f_quadratic)
    c = simulate(*args, reaction=(_f_quadratic, _df_quadratic))
    default = simulate_simple(*args, dt=0.01)
    assert model.get_reaction() is None                  # the module state is untouched
    try:
        model.set_reaction(_f_quadratic, _df_quadratic)
        assert model.get_reaction()[0] is _f_quadratic
        a2 = simulate_simple(*args, dt=0.01)
        b2 = simulate_ale(*args, 0.01, speed_guess=0.4)
        c2 = simulate(*args)
    finally:
        model.set_reaction()
    assert np.array_equal(a.front, a2.front) and np.array_equal(b.front, b2.front)
    assert np.array_equal(c.front, c2.front)
    assert abs(a.front[-1] - default.front[-1]) > 1e-4   # and it is not the default rate


def test_worker_pool_carries_the_reaction_to_spawned_workers():
    """set_reaction is state of one process: a 'spawn' worker would use the
    default f. experiments.util.worker_pool installs the parent's reaction in
    every worker (and refuses reactions that cannot be pickled)."""
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor

    from combustion1d import model
    from experiments.util import worker_pool
    ctx = mp.get_context("spawn")
    try:
        model.set_reaction(_f_quadratic, _df_quadratic)
        expected = float(model.M(0.0))                   # 1/12 for the quadratic rate
        assert expected == pytest.approx(1.0 / 12.0, abs=1e-10)
        with ProcessPoolExecutor(1, mp_context=ctx) as ex:     # plain pool: default rate
            assert ex.submit(_M_in_worker, 0.0).result() == pytest.approx(1.0 / 6.0, abs=1e-12)
        with worker_pool(1, mp_context=ctx) as ex:
            assert ex.submit(_M_in_worker, 0.0).result() == pytest.approx(expected, abs=1e-10)
        model.set_reaction(lambda s: np.maximum(1.0 - s, 0.0))
        with pytest.raises(ValueError, match="worker processes"):
            worker_pool(1, mp_context=ctx)
    finally:
        model.set_reaction()
    with worker_pool(1, mp_context=ctx) as ex:           # default rate again
        assert ex.submit(_M_in_worker, 0.0).result() == pytest.approx(1.0 / 6.0, abs=1e-12)


# --- stored results -----------------------------------------------------------

def test_regression_against_stored_radial_extinction_result():
    """One cheap case of results/paper/B_convergence.json (group B4: N = 2,
    eps = 0.01, dt = 0.04, 19 steps) recomputed with the script that produced
    it."""
    import json
    from pathlib import Path
    pytest.importorskip("matplotlib")
    stored_file = Path(__file__).resolve().parent.parent / "results" / "paper" / "B_convergence.json"
    if not stored_file.exists():
        pytest.skip("stored results not available")
    from experiments.paper.convergence import run_radial
    job = dict(group="B4", N=2, eps=0.01, dt=0.04)
    stored = [r for r in json.loads(stored_file.read_text())
              if all(r.get(k) == v for k, v in job.items())]
    assert len(stored) == 1
    new = run_radial(job)
    assert new["steps"] == stored[0]["steps"] == 19
    for key in ("front_err", "front_err_shifted", "sup_err", "sup_err_int"):
        assert new[key] == pytest.approx(stored[0][key], rel=1e-6), key


def test_provenance_is_attached_to_results():
    import json

    from experiments.util import provenance, to_jsonable, with_provenance
    p = provenance()
    assert {"git", "git_dirty", "python", "numpy", "scipy", "platform"} <= set(p)
    assert p["numpy"] == np.__version__
    d = with_provenance({"a": np.float64(1.0)}, p)
    assert d["provenance"] is p and "provenance" not in {"a": 1}
    rows = with_provenance([{"a": 1}, {"a": 2, "provenance": "old"}], p)
    assert rows[0]["provenance"] is p and rows[1]["provenance"] == "old"
    json.dumps(d, default=to_jsonable)
