"""Example 3: another reaction rate f.

Installs f(s) = ((1 - s)^+)^2 (Lipschitz, nonincreasing, positive for s < 1,
zero for s >= 1); M(w0) is then computed by quadrature. Computes the
travelling wave speed of (P_eps) by shooting and by the ALE scheme (eps = 0.01
and 0.002, dt = 0.02), and the limit speed sqrt(2 M)/A. The default rate is
restored at the end.

Not an item of the paper (all its experiments but one use f(s) = (1 - s)^+).

Run from numerics/:  python -m examples.ex3_custom_reaction      (about 3 s)
"""

import numpy as np

from combustion1d import model
from combustion1d.ale import simulate_ale
from combustion1d.exact import tw_profile_eps, tw_speed_eps, tw_speed_limit

model.set_reaction(lambda s: np.maximum(1.0 - s, 0.0) ** 2,
                   lambda s: np.where(s < 1.0, -2.0 * (1.0 - s), 0.0))
try:
    A, w0, T, dt, x0 = 1.0, 0.0, 4.0, 0.02, 2.0
    print(f"M(0) = {model.M(0.0):.6f} (exact 1/12),  c_0 = {tw_speed_limit(A, w0):.6f}")
    for eps in (0.01, 0.002):
        c_eps = tw_speed_eps(A, w0, eps)
        L = x0 + c_eps * T + 12.0 / c_eps
        res = simulate_ale(lambda x: tw_profile_eps(x - x0, A, w0, eps), w0, eps, L, T, dt,
                           speed_guess=c_eps, exact_initial=True)
        m = res.t >= T / 2
        speed = np.polyfit(res.t[m], res.front[m], 1)[0]
        print(f"eps = {eps:<6} c_eps = {c_eps:.6f}  ALE speed = {speed:.6f}")
finally:
    model.set_reaction()               # back to f(s) = (1 - s)^+
