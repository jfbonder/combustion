"""Experiment in the setting of Theorem cv:thm:main (convergence of (S) to the limit problem).

Radial self-similar solution v of the limit problem (N = 2, w0 = 0.5, extinction time 1),
data (eq:initdata), 0 <= t <= 0.5, no remeshing, and two kinds of meshes:

  fixed : uniform mesh of size h_f = eps/16 that does not move (v_n = 0, mu = max |R'|);
  band  : band mesh (h_f = eps/16 in the band, h_c = 0.01 initially) whose band moves
          rigidly with the prescribed velocity v_n = R'(t_n) of the free boundary (mu = 0);
          the cells outside the rigid zone are stretched by the motion (hmax in the output).

We report E = max_{t_n <= 0.5} max_j |U^n_j - v(x^n_j, t_n)|, to be compared with
eps log(1/eps) + h_f/eps + h_c + dt + mu min(1, dt/eps), and the error of the front
{u = eps} at t = 0.5 (which contains the shift -eps/sqrt(2M) between that level set and
the free boundary).

Run from numerics/:  python -m experiments.paper.convergence_thm     (~3 min, 2 processes)
Results: results/paper/H_convergence_thm.json
"""

from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.linalg import solve_banded

from combustion1d.ale import simulate_ale
from combustion1d.exact import SelfSimilar
from combustion1d.model import f_eps
from combustion1d.solver import MeshParams, find_front, initial_data, laplacian

from .common import Timer, save_json

N, W0, T_EXT, T, L = 2, 0.5, 1.0, 0.5, 5.0
PTS = 16


def _banded(Lap, extra, dt):
    A = (-Lap).tocsc()
    n = A.shape[0]
    ab = np.zeros((3, n))
    ab[0, 1:] = A.diagonal(1)
    ab[1] = A.diagonal(0) + 1.0 / dt + extra
    ab[2, :-1] = A.diagonal(-1)
    return ab


def run_fixed(eps, dt):
    V = SelfSimilar(N, W0, T_EXT)
    h = eps / PTS
    x = np.linspace(0.0, L, int(round(L / h)) + 1)
    u, w = initial_data(x, lambda r: V(r, 0.0), W0, eps)
    y = u + w
    Lap = laplacian(x, N)
    ab0 = _banded(Lap, 0.0, dt)
    nsteps = int(round(T / dt))
    E = float(np.max(np.abs(u - V(x, 0.0))))
    for n in range(nsteps):
        F = f_eps(u, eps)
        y = solve_banded((1, 1), _banded(Lap, F, dt), y / dt)
        u = solve_banded((1, 1), ab0, u / dt - F * y)
        E = max(E, float(np.max(np.abs(u - V(x, (n + 1) * dt)))))
    return E, x.size, nsteps, h, float(find_front(x, u, eps)[0] - V.radius(nsteps * dt))


def run_band(eps, dt):
    V = SelfSimilar(N, W0, T_EXT)
    mp = MeshParams(pts_per_eps=PTS, h_coarse=0.01)
    stat = dict(E=0.0, hmax=0.0)

    def monitor(t, x, u, v):
        stat["E"] = max(stat["E"], float(np.max(np.abs(u - V(x, t)))))
        stat["hmax"] = max(stat["hmax"], float(np.max(np.diff(x))))

    vel = lambda t: -0.5 * V.R / np.sqrt(T_EXT - t)
    r = simulate_ale(lambda x: V(x, 0.0), W0, eps, L, T, dt, mesh=mp, N=N, ell0=0.5,
                     fixed_velocity=vel, remesh=False, monitor=monitor)
    return stat["E"], int(np.mean(r.n_nodes)), r.n_steps, stat["hmax"], float(r.front[-1] - V.radius(r.t[-1]))


def run(job):
    with Timer() as tm:
        E, nodes, steps, h, ef = (run_fixed if job["mesh"] == "fixed" else run_band)(
            job["eps"], job["dt_over_eps"] * job["eps"])
    return dict(job, E=E, front_error=ef, nodes=nodes, steps=steps, hmax=h, seconds=tm.seconds)


def main():
    jobs = [dict(mesh=m, eps=e, dt_over_eps=q)
            for m in ("fixed", "band") for e in (0.02, 0.01, 0.005) for q in (1 / 64, 1 / 16, 0.25, 1.0, 4.0)]
    out = []
    with ProcessPoolExecutor(2) as ex:
        for r in ex.map(run, jobs):
            print("{mesh:5s} eps={eps} dt/eps={dt_over_eps}: E={E:.3e}  E/eps={q:.2f} front={front_error:+.3e} nodes={nodes} "
                  "steps={steps} hmax={hmax:.4f} {seconds:.0f}s".format(q=r["E"] / r["eps"], **r), flush=True)
            out.append(r)
    save_json("H_convergence_thm.json", out)


if __name__ == "__main__":
    main()
