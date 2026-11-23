"""Reaction term and free boundary constant of the model.

The system (Fernandez Bonder--Wolanski) is

    u_t - Delta u = -Y f_eps(u),      Y_t - Delta Y = -Y f_eps(u),

with f_eps(s) = f(s/eps)/eps^2.  Writing w = Y - u (a caloric function),
it reduces to the single equation

    u_t = Delta u - (u + w) f_eps(u).

If (Y_0 - u_0)/eps -> w_0, the limit solves the free boundary problem

    u_t = Delta u  in {u > 0},     |grad u| = sqrt(2 M)  on d{u > 0},

    M(w_0) = int_{-w_0}^1 (s + w_0) f(s) ds.

By default f(s) = (1 - s)^+ : Lipschitz, positive for s < 1, zero for s >= 1,
and then M(w_0) = (1 + w_0)^3 / 6 in closed form. Another reaction rate can be
installed with set_reaction(f, df); every solver in the package picks it up.
"""

import numpy as np
from scipy.integrate import quad


def _f_default(s):
    return np.maximum(1.0 - s, 0.0)


def _df_default(s):
    return np.where(s < 1.0, -1.0, 0.0)


def _M_default(w0):
    """int_{-w0}^1 (s + w0) f(s) ds for f(s) = (1 - s)^+ (needs w0 > -1)."""
    return (1.0 + np.asarray(w0, dtype=float)) ** 3 / 6.0


f, df, _M = _f_default, _df_default, _M_default
version = 0          # incremented by set_reaction; used to invalidate caches


def f_eps(s, eps):
    return f(s / eps) / eps**2


def df_eps(s, eps):
    return df(s / eps) / eps**3


def M(w0):
    """M(w0) = int_{-w0}^1 (s + w0) f(s) ds (needs w0 > -1)."""
    return _M(w0)


def set_reaction(f_new=None, df_new=None, M_new=None):
    """Install the reaction rate f (and its derivative df, used by the
    Jacobian of the BDF solver). f must be vectorised, Lipschitz,
    nonincreasing, positive on (-inf, 1) and zero on [1, inf) (these are the
    hypotheses of the theory; the monotone schemes need f nonincreasing).
    If M_new is not given, M(w0) is computed by quadrature.
    Call set_reaction() with no arguments to restore f(s) = (1 - s)^+.
    """
    global f, df, _M, version
    if f_new is None:
        f, df, _M = _f_default, _df_default, _M_default
    else:
        f = f_new
        if df_new is None:
            h = 1e-7
            df = lambda s: (f_new(s + h) - f_new(s - h)) / (2 * h)
        else:
            df = df_new
        if M_new is None:
            def _M_quad(w0):
                w = np.atleast_1d(np.asarray(w0, dtype=float))
                out = np.array([quad(lambda s: (s + wi) * f_new(s), -wi, 1.0,
                                     limit=200)[0] for wi in w])
                return out.reshape(np.shape(w0)) if np.ndim(w0) else out[0]
            _M = _M_quad
        else:
            _M = M_new
    version += 1


def fb_slope(w0):
    """Free boundary gradient sqrt(2 M(w0))."""
    return np.sqrt(2.0 * M(w0))
