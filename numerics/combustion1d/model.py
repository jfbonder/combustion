"""Reaction term and free boundary constant of the model.

The system (Fernandez Bonder--Wolanski) is

    u_t - Delta u = -Y f_eps(u),      Y_t - Delta Y = -Y f_eps(u),

with f_eps(s) = f(s/eps)/eps^2.  Writing w = Y - u (a caloric function),
it reduces to the single equation

    u_t = Delta u - (u + w) f_eps(u).

If (Y_0 - u_0)/eps -> w_0, the limit solves the free boundary problem

    u_t = Delta u  in {u > 0},     |grad u| = sqrt(2 M)  on d{u > 0},

    M(w_0) = int_{-w_0}^1 (s + w_0) f(s) ds.

We take f(s) = (1 - s)^+ : Lipschitz, positive for s < 1, zero for s >= 1.
Then M(w_0) = (1 + w_0)^3 / 6 in closed form.
"""

import numpy as np


def f(s):
    return np.maximum(1.0 - s, 0.0)


def df(s):
    return np.where(s < 1.0, -1.0, 0.0)


def f_eps(s, eps):
    return f(s / eps) / eps**2


def df_eps(s, eps):
    return df(s / eps) / eps**3


def M(w0):
    """int_{-w0}^1 (s + w0) f(s) ds for f(s) = (1 - s)^+ (needs w0 > -1)."""
    return (1.0 + np.asarray(w0, dtype=float)) ** 3 / 6.0


def fb_slope(w0):
    """Free boundary gradient sqrt(2 M(w0))."""
    return np.sqrt(2.0 * M(w0))
