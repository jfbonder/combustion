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
