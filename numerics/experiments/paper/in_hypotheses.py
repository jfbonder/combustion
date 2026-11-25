"""An experiment inside the hypotheses of Theorem an:thm:ale (uniform rigid mesh, h << eps^2).

Uniform mesh of size h = eps^2/k translating rigidly with a prescribed velocity
v = c_eps + mu, from U^0 = Phi_eps. We report
    E   = max_{t_n <= T} max_j |U^n_j - Phi_eps(x^n_j - x0 - c_eps t_n)|,
    E_s = the same with the front shift removed,
and compare E with the quantities of the theorem, rho_h T = (h/eps^2) T and mu T.

Two more groups are run by this script:

* regular: the same experiment with the rate f(s) = 3 ((1 - s)^+)^2, which is
  C^{1,1}(R); then the theorem holds with rho_h = h^2/eps^3. Mesh sizes
  h = eps^{3/2}/k, k = 2, 4, 8, 16 (rho_h = 1/k^2), dt = 0.1 and 0.01, mu = 0.
  Results: results/paper/G_regular_f.json
* large_eps: the method itself (implicit mesh velocity, default band mesh,
  remeshing) from Phi_eps for eps = 0.1 and 0.05, the values that correspond
  to realistic Zeldovich numbers. Results: results/paper/G_large_eps.json

Run from numerics/:  python -m experiments.paper.in_hypotheses     (~5 min, 2 processes)
Results: results/paper/G_in_hypotheses.json
"""

import numpy as np

from combustion1d import model
from combustion1d.ale import simulate_ale
from combustion1d.solver import MeshParams

from .common import RIGID_MIN_CELL, Timer, save_json, worker_pool
from .tw_error import Monitor, run_ale, setup

T = 3.0


def f_regular(s):
    """f(s) = 3 ((1 - s)^+)^2: nonincreasing, C^{1,1}(R), M(0) = 1/4."""
    return 3.0 * np.maximum(1.0 - s, 0.0) ** 2


def df_regular(s):
    return -6.0 * np.maximum(1.0 - s, 0.0)


def run(job):
    eps, k, dt, mu = job["eps"], job["k"], job["dt"], job["mu"]
    pts = k / eps ** job.get("power", 1.0)     # h_f = eps/pts = eps^(1+power)/k
    h = eps / pts
    mp = MeshParams(pts_per_eps=pts, h_coarse=h, growth=1.0000001)
    c, x0, L, P = setup(0.0, eps, T, pts)
    mon = Monitor(P, x0, c, eps)
    with Timer() as tm:
        r = simulate_ale(lambda x: P(x - x0), 0.0, eps, L, T, dt, mesh=mp, speed_guess=c,
                         exact_initial=True, monitor=mon, fixed_velocity=c + mu, ell0=1e6,
                         remesh_distance=1e9, min_cell=RIGID_MIN_CELL)
    s = mon.summary(T)
    rho = h**2 / eps**3 if job.get("regular") else h / eps**2
    return dict(job, h=h, rho_h=rho, nodes=int(np.mean(r.n_nodes)), steps=r.n_steps,
                E=s["E"], E_shift=s["E_shift"], s_max=s["s_max"], seconds=tm.seconds)


def main():
    jobs = [dict(eps=e, k=k, dt=dt, mu=mu)
            for e in (0.1, 0.05) for k in (8, 16) for dt in (0.1, 0.01) for mu in (0.0, 0.01)]
    out = []
    with worker_pool(2) as ex:
        for r in ex.map(run, jobs):
            print("eps={eps} h=eps^2/{k} dt={dt} mu={mu}: E={E:.3e} E_s={E_shift:.3e} "
                  "s_max={s_max:.3e}  rho_h*T={rt:.3f} nodes={nodes} {seconds:.0f}s".format(
                      rt=r["rho_h"] * T, **r), flush=True)
            out.append(r)
    save_json("G_in_hypotheses.json", out)

    # the method itself at the values of eps of realistic Zeldovich numbers
    jobs = [dict(kind="ale", w0=0.0, eps=e, dt=dt, T=8.0) for e in (0.1, 0.05) for dt in (0.1, 0.01)]
    out = []
    with worker_pool(2) as ex:
        for r in ex.map(run_ale, jobs):
            print("large eps: eps={eps} dt={dt}: E={E:.3e} E_s={E_shift:.3e} s_end={s_end:.3e} "
                  "nodes={nodes:.0f}".format(**r), flush=True)
            r.pop("series", None)
            out.append(r)
    save_json("G_large_eps.json", out)

    # a rate of class C^{1,1}(R): rho_h = h^2/eps^3
    model.set_reaction(f_regular, df_regular)
    try:
        jobs = [dict(eps=e, k=k, dt=dt, mu=0.0, power=0.5, regular=True)
                for e in (0.1, 0.05) for k in (2, 4, 8, 16) for dt in (0.1, 0.01)]
        out = []
        with worker_pool(2) as ex:
            for r in ex.map(run, jobs):
                print("regular f: eps={eps} h=eps^1.5/{k} dt={dt}: E={E:.3e} E_s={E_shift:.3e} "
                      "rho_h*T={rt:.3f} nodes={nodes}".format(rt=r["rho_h"] * T, **r), flush=True)
                out.append(r)
    finally:
        model.set_reaction()
    save_json("G_regular_f.json", out)


if __name__ == "__main__":
    main()
