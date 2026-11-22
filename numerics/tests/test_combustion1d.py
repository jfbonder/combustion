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
