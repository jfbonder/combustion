"""Example 2: a shrinking radial flame, compared with the limit problem.

Initial temperature u0(r) = a (1 - r^2/R0^2)^+ in N = 2 dimensions, with
w0 = 0.5. Solves (P_eps) with the ALE scheme (eps = 0.005, dt = 2.5e-4) and
the free boundary problem (P) with the front-fixing solver up to t = 0.4, and
prints the difference between the two fronts, after adding the distance
eps/sqrt(2M) between the level set {u = eps} and the free boundary.

Not an item of the paper; it is one case of the comparison of Section 6.2
(experiments/paper/limit_examples.py, group C1), with other eps and dt.

Run from numerics/:  python -m examples.ex2_radial_vs_limit      (about 3 s)
"""

import numpy as np

from combustion1d.ale import simulate_ale
from combustion1d.limit import solve_radial
from combustion1d.model import fb_slope

N, w0, eps, dt = 2, 0.5, 0.005, 2.5e-4
lam = fb_slope(w0)                     # free boundary gradient sqrt(2 M(w0))
R0 = 2.0
a = lam * R0 / 2.0                     # makes |u0'(R0)| = lam (compatible data)
u0 = lambda r: a * np.maximum(1.0 - r**2 / R0**2, 0.0)
T = 0.4

lim = solve_radial(u0, w0, N, R0, T, n=800, t_eval=np.linspace(0.0, T, 81))
res = simulate_ale(u0, w0, eps, R0 + 0.5, T, dt, N=N, speed_guess=1.0,
                   ell0=1.0, remesh_distance=0.1)

s_lim = np.interp(res.t, lim.t, lim.s)
shift = eps / lam                      # level set {u = eps} vs free boundary {u = 0}
err = np.max(np.abs(res.front + shift - s_lim))
print(f"limit radius at t = {T}: {lim.s[-1]:.5f};  eps-front + eps/lam: {res.front[-1] + shift:.5f}")
print(f"max_t |front + eps/lam - s(t)| = {err:.2e}  (eps = {eps}, dt = {dt})")
