"""Exact / reference solutions used to validate the solver.

1. Planar travelling waves.
   * Limit problem: u = A (1 - exp(-c (x - c t)))^+ moves to the right
     (into the fresh region u > 0) with speed c_0 = sqrt(2 M) / A.
   * Problem (P_eps) with constant w = eps w_0: the speed c_eps is computed
     by shooting on the travelling wave ODE (high accuracy, independent of
     the PDE solver).

2. Radial self-similar solution of the limit problem (Caffarelli--Vazquez):
       V(x, t) = kappa (T - t)^{1/2} phi(|x| (T - t)^{-1/2}),
   where phi'' + ((N-1)/r - r/2) phi' + phi/2 = 0, phi'(0) = 0, i.e.
       phi(r) = M(-1/2, N/2, r^2/4)   (Kummer's function),
   R is the first zero of phi and kappa = sqrt(2M)/|phi'(R)|.  The support
   {V > 0} is the ball of radius R (T - t)^{1/2}, which shrinks to a point
   at t = T.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import hyp1f1

from . import model


# ---------------------------------------------------------------------------
# Travelling waves
# ---------------------------------------------------------------------------

def tw_speed_limit(A, w0):
    return model.fb_slope(w0) / A


def tw_profile_limit(x, t, A, c, x0=0.0):
    """Limit travelling wave with free boundary at x0 + c t."""
    xi = x - x0 - c * t
    return np.where(xi > 0, A * (1.0 - np.exp(-c * np.maximum(xi, 0.0))), 0.0)


def tw_speed_eps(A, w0, eps, tol=1e-13):
    """Speed of the travelling wave of (P_eps) with w = eps*w0 constant.

    In the variables U = eps*phi, xi = eps*z, sigma = eps*c the profile solves
        phi'' + sigma phi' = (phi + w0) f(phi),
    phi -> -w0 as z -> -inf. Ahead of the reaction zone (phi >= 1) the profile
    is explicit, U = A - (A - eps) exp(-c (xi - xi_1)), which gives
    phi(0) = 1, phi'(0) = sigma (A/eps - 1). We integrate backwards in z and
    bisect on sigma: overshooting -w0 means sigma too large, turning back
    means sigma too small (the energy p^2/2 - G(phi) is monotone in z).
    """
    if not A > eps:
        raise ValueError("need A > eps")
    gamma = A / eps - 1.0
    k = np.sqrt(1.0 + w0)  # decay rate of the tail towards -w0
    s_max = 60.0 / k + 60.0

    def classify(sigma):
        # s = -z ; y = (phi, p) with p = dphi/dz
        def rhs(s, y):
            phi, p = y
            return [-p, -((phi + w0) * model.f(phi) - sigma * p)]

        def overshoot(s, y):
            return y[0] + w0
        overshoot.terminal = True
        overshoot.direction = -1

        def turn(s, y):
            return y[1]
        turn.terminal = True
        turn.direction = -1

        sol = solve_ivp(rhs, (0.0, s_max), [1.0, sigma * gamma],
                        events=(overshoot, turn), rtol=1e-12, atol=1e-14,
                        method="DOP853")
        if sol.t_events[0].size:
            return +1  # sigma too large
        if sol.t_events[1].size:
            return -1  # sigma too small
        return +1 if sol.y[0, -1] + w0 < 0 else -1

    hi = model.fb_slope(w0) / gamma  # zero initial energy: always overshoots
    lo = hi * 1e-3
    while classify(lo) > 0:
        lo *= 0.1
    while hi - lo > tol * hi:
        mid = 0.5 * (lo + hi)
        if classify(mid) > 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi) / eps


# ---------------------------------------------------------------------------
# Radial self-similar solution of the limit problem
# ---------------------------------------------------------------------------

class SelfSimilar:
    """Exact radial solution V of the limit problem, extinct at time T."""

    def __init__(self, N, w0, T=1.0):
        self.N, self.w0, self.T = N, w0, T
        a, b = -0.5, N / 2.0
        self._a, self._b = a, b
        g = lambda z: hyp1f1(a, b, z)
        z_hi = 1.0
        while g(z_hi) > 0:
            z_hi *= 2.0
        z0 = brentq(g, 0.0, z_hi, xtol=1e-15)
        self.R = 2.0 * np.sqrt(z0)
        self.kappa = model.fb_slope(w0) / abs(self.dphi(self.R))

    def phi(self, r):
        return hyp1f1(self._a, self._b, np.asarray(r) ** 2 / 4.0)

    def dphi(self, r):
        r = np.asarray(r)
        return (r / 2.0) * (self._a / self._b) * hyp1f1(self._a + 1, self._b + 1, r**2 / 4.0)

    def radius(self, t):
        return self.R * np.sqrt(self.T - t)

    def __call__(self, r, t):
        tau = self.T - t
        eta = np.asarray(r, dtype=float) / np.sqrt(tau)
        v = self.kappa * np.sqrt(tau) * self.phi(np.minimum(eta, self.R))
        return np.where(eta < self.R, v, 0.0)


# ---------------------------------------------------------------------------
# Travelling wave profile of (P_eps)
# ---------------------------------------------------------------------------

_PROFILE_CACHE = {}


def _tw_inner(A, w0, eps):
    """Inner profile phi on [-s_end, 0] (inner variable z = xi/eps) plus the
    exponential tail beyond. Returns (c_eps, evaluator of phi(z) for z <= 0)."""
    key = (float(A), float(w0), float(eps), model.version)
    if key in _PROFILE_CACHE:
        return _PROFILE_CACHE[key]
    c = tw_speed_eps(A, w0, eps)
    sigma = c * eps
    gamma = A / eps - 1.0
    a = float(model.f(-w0))            # f(-w0) > 0 (needs -w0 < 1)
    # decaying rate of psi = phi + w0 as z -> -inf: psi'' + sigma psi' = a psi
    lam = 0.5 * (-sigma + np.sqrt(sigma**2 + 4.0 * a))
    # integrate backwards (s = -z) until psi is small; beyond that point the
    # shooting error (an exponentially growing mode) would take over, so we
    # continue with the exact linear tail psi(z) = psi* exp(lam (z - z*)).
    tol_psi = 1e-9

    def rhs(s, y):
        phi, p = y
        return [-p, -((phi + w0) * model.f(phi) - sigma * p)]

    def small(s, y):
        return y[0] + w0 - tol_psi
    small.terminal = True
    small.direction = -1

    def turn(s, y):
        return y[1]
    turn.terminal = True
    turn.direction = -1

    sol = solve_ivp(rhs, (0.0, 200.0 / np.sqrt(a) + 200.0), [1.0, sigma * gamma],
                    events=(small, turn), rtol=1e-12, atol=1e-15, method="DOP853",
                    dense_output=True)
    s_end = sol.t[-1]
    psi_end = sol.y[0, -1] + w0
    # tabulate on a fine grid for fast evaluation (cubic Hermite)
    from scipy.interpolate import CubicHermiteSpline
    s = np.linspace(0.0, s_end, int(np.ceil(s_end / 2e-3)) + 1)
    Y = sol.sol(s)
    # phi as a function of z = -s, increasing order in z
    z = -s[::-1]
    spl = CubicHermiteSpline(z, Y[0][::-1], Y[1][::-1])
    z_end = -s_end

    def phi(zz):
        zz = np.asarray(zz, dtype=float)
        out = np.empty_like(zz)
        tail = zz < z_end
        out[tail] = -w0 + psi_end * np.exp(lam * (zz[tail] - z_end))
        out[~tail] = spl(np.minimum(zz[~tail], 0.0))
        return out

    _PROFILE_CACHE[key] = (c, phi)
    return c, phi


def tw_profile_eps(z, A, w0, eps):
    """Travelling wave Phi_eps of (P_eps) with w = eps*w0 constant,
    normalised by Phi_eps(0) = eps (the front {u = eps} at z = 0):

        Phi_eps(z) = A - (A - eps) exp(-c_eps z)          for z >= 0,
        Phi_eps(z) = eps phi(z/eps)                        for z < 0,

    where phi'' + sigma phi' = (phi + w0) f(phi), sigma = eps c_eps,
    phi(0) = 1, phi'(0) = sigma (A/eps - 1), phi -> -w0 as z -> -inf.
    The function U(x, t) = Phi_eps(x - c_eps t) solves
    u_t = u_xx - (u + eps w0) f_eps(u).  Returns an array like z.
    """
    c, phi = _tw_inner(A, w0, eps)
    z = np.asarray(z, dtype=float)
    out = np.empty_like(z)
    pos = z >= 0
    out[pos] = A - (A - eps) * np.exp(-c * z[pos])
    out[~pos] = eps * phi(z[~pos] / eps)
    return out
