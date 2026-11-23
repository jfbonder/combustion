"""An experiment inside the hypotheses of Theorem an:thm:ale (uniform rigid mesh, h << eps^2).

Uniform mesh of size h = eps^2/k translating rigidly with a prescribed velocity
v = c_eps + mu, from U^0 = Phi_eps. We report
    E   = max_{t_n <= T} max_j |U^n_j - Phi_eps(x^n_j - x0 - c_eps t_n)|,
    E_s = the same with the front shift removed,
and compare E with the quantities of the theorem, rho_h T = (h/eps^2) T and mu T.

Run from numerics/:  python -m experiments.paper.in_hypotheses     (~2 min, 2 processes)
Results: results/paper/G_in_hypotheses.json
"""

from concurrent.futures import ProcessPoolExecutor

import numpy as np

from combustion1d.ale import simulate_ale
from combustion1d.solver import MeshParams

from .common import Timer, save_json
from .tw_error import Monitor, setup

T = 3.0


def run(job):
    eps, k, dt, mu = job["eps"], job["k"], job["dt"], job["mu"]
    pts = k / eps                       # h_f = eps/pts = eps^2/k
    h = eps / pts
    mp = MeshParams(pts_per_eps=pts, h_coarse=h, growth=1.0000001)
    c, x0, L, P = setup(0.0, eps, T, pts)
    mon = Monitor(P, x0, c, eps)
    with Timer() as tm:
        r = simulate_ale(lambda x: P(x - x0), 0.0, eps, L, T, dt, mesh=mp, speed_guess=c,
                         exact_initial=True, monitor=mon, fixed_velocity=c + mu, ell0=1e6,
                         remesh_distance=1e9)
    s = mon.summary(T)
    return dict(job, h=h, rho_h=h / eps**2, nodes=int(np.mean(r.n_nodes)), steps=r.n_steps,
                E=s["E"], E_shift=s["E_shift"], s_max=s["s_max"], seconds=tm.seconds)


def main():
    jobs = [dict(eps=e, k=k, dt=dt, mu=mu)
            for e in (0.1, 0.05) for k in (8, 16) for dt in (0.1, 0.01) for mu in (0.0, 0.01)]
    out = []
    with ProcessPoolExecutor(2) as ex:
        for r in ex.map(run, jobs):
            print("eps={eps} h=eps^2/{k} dt={dt} mu={mu}: E={E:.3e} E_s={E_shift:.3e} "
                  "s_max={s_max:.3e}  rho_h*T={rt:.3f} nodes={nodes} {seconds:.0f}s".format(
                      rt=r["rho_h"] * T, **r), flush=True)
            out.append(r)
    save_json("G_in_hypotheses.json", out)


if __name__ == "__main__":
    main()
