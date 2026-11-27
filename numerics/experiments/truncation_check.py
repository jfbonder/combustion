"""Equation (an:eq:solv) of the article (Section 5.4; (37) in the draft in
paper/): truncation error of the monotone scheme (S) on the exact travelling
wave, and its average with the adjoint eigenfunction.

For the travelling wave Phi of (P_eps) with speed c (w = eps*w0 constant) we
evaluate, on a grid x_j = (j + s) h,
    tau_j = (V^{n+1}_j - V^n_j)/dt - delta^2 V^{n+1}_j + f_eps(V^n_j)(V^{n+1}_j + eps w0),
    V^n_j = Phi(x_j),  V^{n+1}_j = Phi(x_j - c dt),
and report
  * max_j |tau_j|  (what a pointwise barrier must absorb),
  * the drift predicted by the solvability condition of the linearized
    problem, delta_c = sum tau_j psi*_j / sum Phi'_j psi*_j, psi* = e^{cz} Phi',
    averaged over 16 offsets s in [0, 1)
(the numerical front then moves at c + delta_c). No PDE is solved.

Prints two tables and writes them to results/paper/truncation_check.json:
  "time_step"  w0 = 0, h = eps/64, dt = kappa eps^2, eps in {0.02, 0.01, 0.005},
               kappa in {1, 0.25}; keys eps, kappa, dt, h, max_tau, dc and the
               normalised values max_tau_eps2_over_dt, dc_over_dt_eps. The three
               values of dc_over_dt_eps with kappa = 1 are those of (an:eq:solv).
  "mesh_size"  dt = 1e-3 eps^2, h = eps/m, m in {8, 16}, w0 in {0, 1},
               eps in {0.01, 0.005}; keys w0, eps, m, h, dt, max_tau, dc and
               max_tau_eps2_over_h, dc_over_h_eps_sq (not quoted in the paper).

Run from numerics/:  python -m experiments.truncation_check   (about 10 s)
"""

from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from combustion1d.exact import tw_speed_eps
from combustion1d.model import f, f_eps
from experiments.util import write_json

OUT = Path(__file__).resolve().parent.parent / "results" / "paper"


class Profile:
    """Travelling wave Phi(z), Phi(0) = eps, z = x - c t, with derivative."""

    def __init__(self, A, w0, eps):
        self.A, self.w0, self.eps = A, w0, eps
        self.c = tw_speed_eps(A, w0, eps)
        sig = eps * self.c
        gam = A / eps - 1.0

        def rhs(s, y):          # s = -zeta
            phi, p = y
            return [-p, -((phi + w0) * f(phi) - sig * p)]

        def near_saddle(s, y):
            return y[0] + w0 - 1e-6
        near_saddle.terminal = True

        def turn(s, y):          # numerical departure from the unstable manifold
            return y[1]
        turn.terminal = True
        sol = solve_ivp(rhs, (0, 400), [1.0, sig * gam], rtol=1e-12, atol=1e-14,
                        dense_output=True, events=(near_saddle, turn), method="DOP853")
        if sol.t_events[0].size == 0:
            raise RuntimeError("profile did not reach the saddle neighbourhood")
        self.s_end = sol.t[-1]
        self.inner = sol.sol
        lam = 0.5 * (-sig + np.sqrt(sig**2 + 4 * f(-w0)))   # decay rate (in zeta)
        self.lam = lam
        self.tail_amp = (sol.y[0, -1] + w0) * np.exp(lam * self.s_end)

    def __call__(self, z):
        z = np.asarray(z, dtype=float)
        e, A, c, w0 = self.eps, self.A, self.c, self.w0
        out = np.empty_like(z)
        ahead = z >= 0
        out[ahead] = A - (A - e) * np.exp(-c * z[ahead])
        s = -z[~ahead] / e
        inner = s <= self.s_end
        vals = np.empty(s.size)
        vals[inner] = self.inner(s[inner])[0]
        vals[~inner] = -w0 + self.tail_amp * np.exp(-self.lam * s[~inner])
        out[~ahead] = e * vals
        return out

    def dphi(self, z, d=1e-7):
        return (self(z + d) - self(z - d)) / (2 * d)


def truncation(P, h, dt, s, Zl=60.0, Zr=None):
    e, c, w0 = P.eps, P.c, P.w0
    if Zr is None:
        Zr = 40.0 / c
    j = np.arange(int(np.floor(-Zl * e / h)) - 1, int(np.ceil(Zr / h)) + 2)
    x = (j + s) * h
    Vn = P(x)
    Vn1 = P(x - c * dt)
    lap = (Vn1[2:] - 2 * Vn1[1:-1] + Vn1[:-2]) / h**2
    tau = (Vn1[1:-1] - Vn[1:-1]) / dt - lap + f_eps(Vn[1:-1], e) * (Vn1[1:-1] + e * w0)
    xi = x[1:-1]
    psi = np.exp(c * np.minimum(xi, 50.0 / c)) * P.dphi(xi)
    num = np.sum(tau * psi) * h
    den = np.sum(P.dphi(xi) * psi) * h
    return np.max(np.abs(tau)), num / den


def report(A, w0, eps, h, dt, nshift=16):
    P = Profile(A, w0, eps)
    res = [truncation(P, h, dt, s) for s in (np.arange(nshift) + 0.5) / nshift]
    return max(r[0] for r in res), float(np.mean([r[1] for r in res]))


def main():
    A = 1.0
    rows_dt, rows_h = [], []
    print("time step: w0=0, h = eps/64, dt = kappa eps^2")
    for eps in (0.02, 0.01, 0.005):
        for kappa in (1.0, 0.25):
            mx, dc = report(A, 0.0, eps, eps / 64, kappa * eps**2)
            print(f"  eps={eps:<6} kappa={kappa:<5} max|tau|*eps^2/dt={mx*eps**2/(kappa*eps**2):7.3f}"
                  f"   predicted dc/(dt/eps)={dc/(kappa*eps):+.4f}")
            rows_dt.append(dict(A=A, w0=0.0, eps=eps, kappa=kappa, dt=kappa * eps**2, h=eps / 64,
                                max_tau=mx, dc=dc, max_tau_eps2_over_dt=mx * eps**2 / (kappa * eps**2),
                                dc_over_dt_eps=dc / (kappa * eps)))
    print("mesh size: dt = 1e-3 eps^2, h = eps/m")
    for w0 in (0.0, 1.0):
        for eps in (0.01, 0.005):
            for m in (8, 16):
                h = eps / m
                mx, dc = report(A, w0, eps, h, 1e-3 * eps**2)
                print(f"  w0={w0} eps={eps:<6} h=eps/{m:<3} max|tau|*eps^2/h={mx*eps**2/h:7.3f}"
                      f"   predicted dc={dc:+.2e}  dc/(h/eps)^2={dc/(h/eps)**2:+.4f}")
                rows_h.append(dict(A=A, w0=w0, eps=eps, m=m, h=h, dt=1e-3 * eps**2, max_tau=mx,
                                   dc=dc, max_tau_eps2_over_h=mx * eps**2 / h,
                                   dc_over_h_eps_sq=dc / (h / eps)**2))
    write_json(OUT / "truncation_check.json",
               dict(nshift=16, time_step=rows_dt, mesh_size=rows_h))
    print("wrote", OUT / "truncation_check.json")


if __name__ == "__main__":
    main()
