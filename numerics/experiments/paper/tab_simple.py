"""Table tab:simple: Eulerian scheme (S) (simple.py) with dt proportional to eps^2.

Planar travelling wave, w0 = 0, A = 1, eps in {0.02, 0.01}, dt in {4 eps^2,
eps^2, eps^2/4} (exactly the nominal step; see simple.py), default band mesh
(h_f = eps/16, h_c = 0.01). Speed error c_(S) - c_eps with c_eps from the
shooting code.

Two set-ups:
  'B1'  : datum A(1 - e^{-c_eps(x-2)})^+, T = 12 (as group B1); speed by least
          squares on [T/2, T] and by end points (x_f(T) - x_f(T/2))/(T/2);
  'Phi' : datum Phi_eps (exact eps-wave), T = 6; speed from the front
          position error s(T)/T and least squares of s on [T/2, T].

Run from numerics/:
    OMP_NUM_THREADS=1 ../.venv/bin/python -m experiments.paper.tab_simple
"""

import numpy as np

from combustion1d.exact import tw_profile_eps, tw_profile_limit, tw_speed_eps, tw_speed_limit
from combustion1d.simple import simulate_simple
from combustion1d.solver import MeshParams, front_speed

from .common import Timer, save_json, worker_pool

OLD = {(0.02, 4): -1.85e-2, (0.02, 1): -4.60e-3, (0.02, 0.25): -9.30e-4,
       (0.01, 4): -9.41e-3, (0.01, 1): -2.21e-3, (0.01, 0.25): -3.69e-4}


def run(job):
    eps, k, setup = job
    dt = k * eps**2
    c = tw_speed_eps(1.0, 0.0, eps)
    c0 = float(tw_speed_limit(1.0, 0.0))
    if setup == "B1":
        T = 12.0
        u0 = lambda x: tw_profile_limit(x, 0.0, 1.0, c, x0=2.0)
        L = 2.0 + c0 * T + 12.0 / c0 + 2.0
        with Timer() as tm:
            r = simulate_simple(u0, 0.0, eps, L, T, dt=dt, mesh=MeshParams(), speed_guess=c0)
    else:
        T = 6.0
        h_f = eps / 16
        x0 = (round(2.0 / h_f) + 0.37) * h_f
        u0 = lambda x: tw_profile_eps(x - x0, 1.0, 0.0, eps)
        L = x0 + c * T + 12.0 / c + 2.0
        with Timer() as tm:
            r = simulate_simple(u0, 0.0, eps, L, T, dt=dt, mesh=MeshParams(), speed_guess=c,
                                exact_initial=True)
    i2 = np.searchsorted(r.t, 0.5 * T - 1e-9)
    out = dict(eps=eps, k=k, dt=dt, setup=setup, T=T, c_eps=c,
               err_lsq=front_speed(r, 0.5 * T) - c,
               err_end=(r.front[-1] - r.front[i2]) / (r.t[-1] - r.t[i2]) - c,
               steps=r.n_steps, dts=r.dts, seconds=tm.seconds, old=OLD[(eps, k)])
    if setup == "Phi":
        out["s_end"] = r.front[-1] - x0 - c * r.t[-1]
        out["err_pos"] = out["s_end"] / r.t[-1]
    return out


def main():
    jobs = [(e, k, s) for s in ("B1", "Phi") for e in (0.01, 0.02) for k in (0.25, 1, 4)]
    out = []
    with worker_pool(4) as ex:
        for r in ex.map(run, jobs):
            print("{setup:<3} eps={eps} dt={dt:.2e} (k={k}) lsq={err_lsq:+.3e} end={err_end:+.3e} "
                  "old={old:+.2e} steps={steps} dts={dts} {seconds:.0f}s".format(**r), flush=True)
            out.append(r)
    save_json("tab_simple.json", out)


if __name__ == "__main__":
    main()
