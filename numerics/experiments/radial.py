"""Radial self-similar solution (N = 2, 3) with the BDF solver. Not used in the
paper (the corresponding experiment of the paper, with the ALE scheme, is group
B4 of experiments/paper/convergence.py).

The limit problem has the exact solution
    V(x, t) = kappa (T - t)^{1/2} phi(|x| (T - t)^{-1/2}),
whose support is a ball of radius R (T - t)^{1/2} that shrinks to a point at
t = T (cf. Proposition 4.3 of the SIAM paper). V is radially decreasing, so
u_0 = V(., 0) is starshaped, and w0 is constant: the hypotheses of the
uniqueness theorem hold and u^eps must converge to V. We start (P_eps) from
u_0^eps built from V(., 0) (w0 = 0.5, T = 1) and measure, for eps from 0.05
to 0.001,
    * the front radius error  max_t |r_eps(t) - R (T - t)^{1/2}|,
    * the sup error           max_x |(u^eps)^+(x, t) - V(x, t)| at t = 0.75 T
      (in the burnt region u^eps = -eps w0, which is trivially O(eps)), and
      the same error restricted to r < R(t) - 20 eps (outside the reaction
      layer, where (u^eps)^+ ~ eps psi(d/eps) differs from the linear profile
      of V by O(eps) by construction).

Output: results/radial.json, a list of 12 records with the keys N, eps, w0, R,
kappa, front_err, sup_err, sup_err_int, t, front, snapshots, steps, nodes,
seconds. experiments/plots.py draws from it figures/radial_convergence.png
and figures/radial_snapshots.png.

Run from numerics/:  python -m experiments.radial   (about 1 min, 4 processes)
"""

import time
from pathlib import Path

import numpy as np

from combustion1d.exact import SelfSimilar
from combustion1d.solver import simulate
from experiments.util import worker_pool, write_json

T = 1.0
T_END = 0.75
W0 = 0.5
NS = (2, 3)
EPS = (0.05, 0.02, 0.01, 0.005, 0.002, 0.001)
SNAP = (0.0, 0.25, 0.5, 0.75)
OUT = Path(__file__).resolve().parent.parent / "results"


def run(args):
    N, eps = args
    V = SelfSimilar(N, W0, T)
    L = 1.5 * V.R * np.sqrt(T) + 1.0
    t = time.time()
    r = simulate(lambda x: V(x, 0.0), W0, eps, L, T_END, N=N,
                 speed_guess=V.R / (2 * np.sqrt(T)), snapshot_times=SNAP,
                 max_dt_remesh=0.01)
    m = r.t >= 0.05                        # skip the initial layer
    front_err = float(np.max(np.abs(r.front[m] - V.radius(r.t[m]))))
    x, u, w = r.snapshots[T_END]
    err = np.abs(np.maximum(u, 0.0) - V(x, T_END))
    sup_err = float(np.max(err))
    interior = x < V.radius(T_END) - 20 * eps   # away from the reaction layer
    sup_err_int = float(np.max(err[interior]))
    snaps = {str(s): dict(x=r.snapshots[s][0].tolist(), u=r.snapshots[s][1].tolist())
             for s in SNAP}
    return dict(N=N, eps=eps, w0=W0, R=V.R, kappa=V.kappa, front_err=front_err,
                sup_err=sup_err, sup_err_int=sup_err_int, t=r.t.tolist(), front=r.front.tolist(),
                snapshots=snaps if eps in (0.05, 0.01) else None,
                steps=r.n_steps, nodes=float(np.mean(r.n_nodes)),
                seconds=time.time() - t)


def main():
    OUT.mkdir(exist_ok=True)
    jobs = [(N, eps) for eps in EPS for N in NS]
    res = []
    with worker_pool(4) as ex:
        for r in ex.map(run, jobs):
            print("N={N} eps={eps:<6} front_err={front_err:.3e} sup_err={sup_err:.3e} int={sup_err_int:.3e} "
                  "nodes={nodes:.0f} steps={steps} {seconds:.0f}s".format(**r), flush=True)
            res.append(r)
    write_json(OUT / "radial.json", res)


if __name__ == "__main__":
    main()
