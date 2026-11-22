"""Direct solver for the limit free boundary problem (P) in 1D / radial symmetry.

    u_t = Delta u  in {u > 0},   u = 0 and |grad u| = lam = sqrt(2M(w0)) on the
    free boundary.

Front fixing (Landau transformation) on a uniform grid in xi in [0, 1]:

* radial (fresh region = ball {r < s(t)}): r = s xi,
      v_t = (v_xixi + (N-1)/xi v_xi)/s^2 + xi (s'/s) v_xi,
      v_xi(0) = 0,  v(1) = 0,  v_xi(1) = -s lam;
* planar (fresh region = {s(t) < x < L}, zero flux at L): x = s + xi (L - s),
      v_t = v_xixi/l^2 + s' (1 - xi)/l v_xi,   l = L - s,
      v(0) = 0,  v_xi(0) = l lam,  v_xi(1) = 0.

Both free boundary conditions are imposed: u = 0 as a Dirichlet condition and
|grad u| = lam through a ghost node. The front velocity follows from keeping
u = 0 at the front, d/dt v(front) = 0, which gives s' = Delta u / lam (radial)
and s' = -u_xx / lam (planar) evaluated at the front with the ghost node.
Second order in space; BDF in time.
"""

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
from scipy.integrate import solve_ivp

from .model import fb_slope


@dataclass
class LimitResult:
    t: np.ndarray
    s: np.ndarray
    xi: np.ndarray           # mapped coordinate in [0, 1] (g(xi) for planar)
    v: np.ndarray            # v[:, k] on the grid xi at time t[k] (front value 0 included)
    geometry: str
    L: float = None

    def x_of(self, k):
        if self.geometry == "radial":
            return self.s[k] * self.xi
        return self.s[k] + self.xi * (self.L - self.s[k])

    def u_on(self, x, k):
        """Solution at time t[k] evaluated at points x (0 in the burnt region)."""
        xs, vs = self.x_of(k), self.v[:, k]     # xs is increasing in both cases
        return np.interp(x, xs, vs, left=0.0, right=0.0 if self.geometry == "radial" else vs[-1])


def _sparsity(m):
    """Tridiagonal in v, plus full coupling with s (last unknown) and with the
    unknown that enters s' (all rows depend on s')."""
    S = sp.lil_matrix((m + 1, m + 1))
    for i in range(m):
        for j in (i - 1, i, i + 1):
            if 0 <= j < m:
                S[i, j] = 1
    S[:, m - 1] = 1    # s' depends on v_{n-1}
    S[:, m] = 1
    S[m, :] = 1
    return S


def solve_radial(u0, w0, N, s0, t_end, n=400, t_eval=None, rtol=1e-8, atol=1e-10):
    lam = float(fb_slope(w0))
    xi = np.linspace(0.0, 1.0, n + 1)
    h = xi[1] - xi[0]
    v0 = u0(s0 * xi[:-1])                 # unknowns v_0 .. v_{n-1}; v_n = 0

    def front_speed(v, s):
        vxx = (2.0 * v[-1] - 2.0 * h * s * lam) / h**2      # ghost from v_xi = -s lam
        lap = (vxx + (N - 1) * (-s * lam)) / s**2
        return lap / lam

    def rhs(_, y):
        v, s = y[:-1], y[-1]
        sp_ = front_speed(v, s)
        vf = np.append(v, 0.0)
        dv = np.empty(n)
        # xi = 0: symmetry, Delta = N v_xixi
        dv[0] = N * 2.0 * (vf[1] - vf[0]) / h**2 / s**2
        i = np.arange(1, n)
        vx = (vf[i + 1] - vf[i - 1]) / (2 * h)
        vxx = (vf[i + 1] - 2 * vf[i] + vf[i - 1]) / h**2
        dv[1:] = (vxx + (N - 1) / xi[i] * vx) / s**2 + xi[i] * (sp_ / s) * vx
        return np.append(dv, sp_)

    def stop(_, y):                       # extinction
        return y[-1] - 0.05 * s0
    stop.terminal = True

    sol = solve_ivp(rhs, (0.0, t_end), np.append(v0, s0), method="BDF",
                    jac_sparsity=_sparsity(n), rtol=rtol, atol=atol,
                    t_eval=t_eval, events=stop)
    v = np.vstack((sol.y[:-1], np.zeros(sol.t.size)))
    return LimitResult(t=sol.t, s=sol.y[-1], xi=xi, v=v, geometry="radial")


def solve_planar(u0, w0, L, s0, t_end, n=800, beta=4.0, t_eval=None, rtol=1e-8, atol=1e-10):
    """Planar case on a grid clustered at the front: x = s + l g(xi),
    g(xi) = (e^{beta xi} - 1)/(e^beta - 1), l = L - s. Then
        u_x = v_xi/(l g'),  u_xx = (v_xixi - (g''/g') v_xi)/(l g')^2,
        v_t = u_xx + s' (1 - g)/(l g') v_xi,
    with v(0) = 0, v_xi(0) = l g'(0) lam (ghost node) and v_xi(1) = 0."""
    lam = float(fb_slope(w0))
    xi = np.linspace(0.0, 1.0, n + 1)
    h = xi[1] - xi[0]
    E = np.expm1(beta)
    g = np.expm1(beta * xi) / E
    g1 = beta * np.exp(beta * xi) / E
    g2 = beta * g1
    v0 = u0(s0 + g[1:] * (L - s0))         # unknowns v_1 .. v_n; v_0 = 0

    def uxx_front(v, ell):
        vx0 = ell * g1[0] * lam
        vxx = (2.0 * v[0] - 2.0 * h * vx0) / h**2
        return (vxx - g2[0] / g1[0] * vx0) / (ell * g1[0]) ** 2

    def rhs(_, y):
        v, s = y[:-1], y[-1]
        ell = L - s
        sp_ = -uxx_front(v, ell) / lam
        vf = np.concatenate(([0.0], v, [v[-2]]))             # ghost: v_xi(1) = 0
        i = np.arange(1, n + 1)
        vx = (vf[i + 1] - vf[i - 1]) / (2 * h)
        vxx = (vf[i + 1] - 2 * vf[i] + vf[i - 1]) / h**2
        uxx = (vxx - g2[i] / g1[i] * vx) / (ell * g1[i]) ** 2
        dv = uxx + sp_ * (1.0 - g[i]) / (ell * g1[i]) * vx
        return np.append(dv, sp_)

    def stop(_, y):                        # front too close to the far boundary
        return (L - y[-1]) - 0.2 * L
    stop.terminal = True

    # every row depends on v_1 (through s') and on s
    S = sp.lil_matrix((n + 1, n + 1))
    for i in range(n):
        for j in (i - 1, i, i + 1):
            if 0 <= j < n:
                S[i, j] = 1
    S[:, 0] = 1
    S[:, n] = 1
    S[n, :] = 1
    sol = solve_ivp(rhs, (0.0, t_end), np.append(v0, s0), method="BDF",
                    jac_sparsity=S, rtol=rtol, atol=atol, t_eval=t_eval, events=stop)
    v = np.vstack((np.zeros(sol.t.size), sol.y[:-1]))
    res = LimitResult(t=sol.t, s=sol.y[-1], xi=g, v=v, geometry="planar", L=L)
    return res
