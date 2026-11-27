"""Figure 4 of the article (label fig:instability, figures/tw_instability.png,
Section 6.4): planar travelling waves with the BDF solver.

A. Convergence of the flame speed as eps -> 0, for w0 in {-0.5, 0, 1} and
   eps from 0.1 to 0.001 (not in the paper; figures/tw_speed_convergence.png):
   c_pde  : speed measured from the PDE solution of (P_eps) (BDF solver),
   c_eps  : exact speed of the (P_eps) travelling wave (ODE shooting),
   c_0    : speed of the limit problem, sqrt(2 M(w0)) / A.
   |c_pde - c_eps| measures the discretization error, |c_eps - c_0| the
   O(eps) distance to the free boundary problem.

B. Dependence of the limit on w0 (this is Figure 4): the SAME initial
   temperature u_0 (the limit travelling wave for w0 = 0) and
   Y_0 = u_0 + eps w0 with three different w0 produce three different flames,
   with speeds approaching sqrt(2M(w0))/A. eps = 0.005, 0 <= t <= 20.

Output: results/traveling_wave.json = {"A": [...], "B": [...]}; records of A:
w0, eps, c0, c_eps, c_pde, steps, nodes, seconds; records of B: w0, c0, the
front (t, front) and the solution at t = 20 (x, u, Y). The figures are drawn
by experiments/plots.py.

Run from numerics/:  python -m experiments.traveling_wave   (about 21 min, 4 processes;
                     most of it is part A with eps = 0.002 and 0.001)
then                 python -m experiments.plots

Numbers of figures and sections are those of the draft in paper/.
"""

import time
from pathlib import Path

import numpy as np

from combustion1d.exact import tw_profile_limit, tw_speed_eps, tw_speed_limit
from combustion1d.solver import MeshParams, front_speed, simulate
from experiments.util import worker_pool, write_json

A = 1.0
X0 = 2.0
DIST = 8.0
W0S = (-0.5, 0.0, 1.0)
EPS = (0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001)
EPS_B = 0.005
OUT = Path(__file__).resolve().parent.parent / "results"


def run_A(args):
    w0, eps = args
    c0 = tw_speed_limit(A, w0)
    u0 = lambda x: tw_profile_limit(x, 0.0, A, c0, x0=X0)
    T = DIST / c0
    L = X0 + DIST + 12.0 / c0 + 2.0
    t = time.time()
    # the error ahead of the front is O((c0 h_coarse)^2): scale h_coarse with 1/c0
    mesh = MeshParams(h_coarse=float(np.clip(0.01 / c0, 0.01, 0.05)))
    r = simulate(u0, w0, eps, L, T, N=1, speed_guess=c0, mesh=mesh)
    return dict(w0=w0, eps=eps, c0=c0, c_eps=tw_speed_eps(A, w0, eps),
                c_pde=front_speed(r, 0.5 * T), steps=r.n_steps,
                nodes=float(np.mean(r.n_nodes)), seconds=time.time() - t)


def run_B(w0):
    c_ref = tw_speed_limit(A, 0.0)          # same initial data for every w0
    u0 = lambda x: tw_profile_limit(x, 0.0, A, c_ref, x0=X0)
    T = 20.0
    L = X0 + tw_speed_limit(A, max(W0S)) * T + 15.0
    r = simulate(u0, w0, EPS_B, L, T, N=1, speed_guess=tw_speed_limit(A, w0),
                 snapshot_times=(T,), mesh=MeshParams(h_coarse=0.02))
    x, u, w = r.snapshots[T]
    return dict(w0=w0, t=r.t.tolist(), front=r.front.tolist(),
                c0=tw_speed_limit(A, w0), x=x.tolist(), u=u.tolist(), Y=(u + w).tolist())


def main():
    OUT.mkdir(exist_ok=True)
    jobs = [(w0, eps) for eps in EPS for w0 in W0S]
    with worker_pool(4) as ex:
        resB = list(ex.map(run_B, W0S))
        resA = []
        for res in ex.map(run_A, jobs):
            print("A: w0={w0:+.1f} eps={eps:<6} c_pde={c_pde:.6f} c_eps={c_eps:.6f} "
                  "c0={c0:.6f} nodes={nodes:.0f} steps={steps} {seconds:.0f}s".format(**res),
                  flush=True)
            resA.append(res)
    write_json(OUT / "traveling_wave.json", dict(A=resA, B=resB))


if __name__ == "__main__":
    main()
