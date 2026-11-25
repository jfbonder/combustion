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
and then M(w_0) = (1 + w_0)^3 / 6 in closed form.

Another reaction rate can be used in two ways:

* set_reaction(f, df, M) installs it for the whole package. This is STATE OF
  THIS MODULE IN THE CURRENT PROCESS: every solver and exact solution called
  afterwards in this process uses it, until set_reaction() restores the
  default. It is not thread safe and it is NOT inherited by worker processes
  started with the "spawn" method (the default on macOS and Windows): a child
  process re-imports this module and gets the default f. To use a custom rate
  in workers, install it in each worker, e.g.

      ProcessPoolExecutor(4, initializer=model.install_reaction,
                          initargs=(model.get_reaction(),))

  (experiments/util.py: worker_pool does exactly this). The functions must
  then be picklable, i.e. defined at module level (no lambdas).

* the solvers simulate, simulate_ale and simulate_simple accept
  reaction=f or reaction=(f, df) for one call, without touching the module
  state (see resolve_reaction). The exact solutions and M(w0) always use the
  installed rate.
"""

import numpy as np
from scipy.integrate import quad


def _f_default(s):
    return np.maximum(1.0 - s, 0.0)


def _df_default(s):
    # derivative of (1 - s)^+ ; at the kink s = 1 we take the right value 0
    # (only used by the Jacobian of the BDF solver, where any value in
    # [-1, 0] is an admissible generalised derivative)
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

    This changes module-level state of the current process only: see the
    module docstring (worker processes started with "spawn" do not inherit
    it; use install_reaction as pool initializer, or the per-call argument
    `reaction` of the solvers).
    """
    global f, df, _M, version
    if f_new is None:
        f, df, _M = _f_default, _df_default, _M_default
    else:
        f = f_new
        if df_new is None:
            df = _fd_derivative(f_new)
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


def _fd_derivative(f_new, h=1e-7):
    return lambda s: (f_new(s + h) - f_new(s - h)) / (2 * h)


def get_reaction():
    """The installed reaction as a tuple (f, df, M), or None for the default
    f(s) = (1 - s)^+. The value can be passed to install_reaction (in another
    process, if the functions are picklable)."""
    if f is _f_default:
        return None
    return (f, df, _M)


def install_reaction(reaction=None):
    """Install the value returned by get_reaction() (None = default rate).
    Meant as `initializer` of a process pool, so that the workers use the
    same reaction as the parent whatever the start method. Does nothing if
    the requested reaction is already installed."""
    if reaction is None:
        if f is not _f_default:
            set_reaction()
        return
    f_new, df_new, M_new = (tuple(reaction) + (None, None))[:3]
    if f_new is f and (df_new is None or df_new is df) and (M_new is None or M_new is _M):
        return
    set_reaction(f_new, df_new, M_new)


def resolve_reaction(reaction=None):
    """(f, df) for the per-call argument `reaction` of the solvers.

    reaction=None: the rate installed in this module at the time of the call
    (set_reaction); reaction=f or reaction=(f, df): that rate, for this call
    only (df by centred differences if it is not given). f must satisfy the
    hypotheses listed in set_reaction."""
    if reaction is None:
        return f, df
    if callable(reaction):
        return reaction, _fd_derivative(reaction)
    f_new, df_new = (tuple(reaction) + (None,))[:2]
    return f_new, (df_new if df_new is not None else _fd_derivative(f_new))


def fb_slope(w0):
    """Free boundary gradient sqrt(2 M(w0))."""
    return np.sqrt(2.0 * M(w0))
