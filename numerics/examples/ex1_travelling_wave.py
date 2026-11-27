"""Example 1: a planar travelling wave with the ALE scheme.

Starts from the exact travelling wave Phi_eps of (P_eps) (w0 = 0, eps = 0.002,
h_f = eps/16, dt = 0.02 = 10 eps) and prints the computed speed, the sup error
against the exact solution at the final time t = 4, and the cost (nodes, steps,
remeshings, tridiagonal solves).

Not an item of the paper; it is a small version of the experiments of
Section S1 of the supplement (experiments/paper/tw_error.py).

Run from numerics/:  python -m examples.ex1_travelling_wave      (about 2 s)
"""

import numpy as np

from combustion1d.ale import simulate_ale
from combustion1d.exact import tw_profile_eps, tw_speed_eps, tw_speed_limit
from combustion1d.solver import MeshParams

A, w0, eps = 1.0, 0.0, 0.002          # wave height, w0 = lim (Y - u)/eps, small parameter
c_eps = tw_speed_eps(A, w0, eps)      # exact speed of the eps-wave (ODE shooting)
c0 = tw_speed_limit(A, w0)            # speed of the limit free boundary problem
x0, T, dt = 2.0, 4.0, 0.02            # initial front position, final time, time step
L = x0 + c_eps * T + 12.0 / c_eps     # domain [0, L], long enough ahead of the front

u0 = lambda x: tw_profile_eps(x - x0, A, w0, eps)
res = simulate_ale(u0, w0, eps, L, T, dt,
                   mesh=MeshParams(pts_per_eps=16, h_coarse=0.01),   # h_f = eps/16
                   speed_guess=c_eps, exact_initial=True, snapshot_times=(T,))

x, u, w = res.snapshots[T]
err = np.max(np.abs(u - tw_profile_eps(x - x0 - c_eps * T, A, w0, eps)))
speed = np.polyfit(res.t[res.t >= T / 2], res.front[res.t >= T / 2], 1)[0]
print(f"c_0 = {c0:.6f}   c_eps = {c_eps:.6f}   computed speed = {speed:.6f}")
print(f"max |U - Phi_eps| at t = {T}: {err:.2e}")
# The speed error (here about -2e-4) is dominated by "lattice pinning": the front
# stays at a fixed position inside a band cell (here exactly on a node, since
# x0 is a multiple of h_f), and the speed depends on that position with an
# amplitude O((h_f/eps)^2); see the supplement, Section S1 (figure
# D3_phase_pinning.png).
print(f"{x.size} nodes, {res.n_steps} steps of size {dt}, {res.n_remesh} remeshings, "
      f"{res.n_solves} tridiagonal solves")
