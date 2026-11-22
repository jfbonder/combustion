"""Experiment 3: (P_eps) versus the limit problem (P), both computed.

The limit problem is solved by front fixing (combustion1d.limit, second
order, validated against the exact solutions). We compare with (P_eps) for
initial data that are NOT exact solutions:

  radial  (N = 2, 3): u0(r) = a (1 - r^2/R0^2)^+,  a = lam R0 / 2,
  planar  (N = 1):    u0(x) = A tanh(lam (x - x0)/A)^+,
with lam = sqrt(2 M(w0)) (so that the data are compatible with the free
boundary condition). We record, as eps -> 0,
  * front error  max_t |x_eps(t) - s(t)|,
  * the same after removing the trivial shift eps/lam between the level set
    {u = eps} that we track and the free boundary,
  * sup error    max_x |(u^eps)^+ - u| at the final time.

Run from numerics/:  python -m experiments.limit_comparison
"""

import json
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from combustion1d.limit import solve_planar, solve_radial
from combustion1d.model import fb_slope
from combustion1d.solver import MeshParams, simulate

W0 = 0.5
EPS = (0.05, 0.02, 0.01, 0.005, 0.002, 0.001)
OUT = Path(__file__).resolve().parent.parent / "results"

R0 = 2.0
A, X0, L_PLANAR = 1.0, 2.0, 16.0


def case(geom):
    lam = float(fb_slope(W0))
    if geom == "planar":
        u0 = lambda x: np.where(x > X0, A * np.tanh(lam * (x - X0) / A), 0.0)
        return dict(N=1, u0=u0, L=L_PLANAR, t_end=6.0, sign=+1.0, s0=X0)
    N = int(geom[-1])
    a = lam * R0 / 2.0
    u0 = lambda r: a * np.maximum(1.0 - (r / R0) ** 2, 0.0)
    return dict(N=N, u0=u0, L=1.5 * R0 + 1.0, t_end=None, sign=-1.0, s0=R0)


def limit_solution(geom, t_end, n=1600):
    c = case(geom)
    t_eval = np.linspace(0.0, t_end, 601)
    if geom == "planar":
        return solve_planar(c["u0"], W0, c["L"], c["s0"], t_end, n=n, t_eval=t_eval)
    return solve_radial(c["u0"], W0, c["N"], c["s0"], t_end, n=n, t_eval=t_eval)


def t_end_of(geom):
    c = case(geom)
    if c["t_end"] is not None:
        return c["t_end"]
    # radial: stop when the limit ball has radius R0/2
    r = solve_radial(c["u0"], W0, c["N"], c["s0"], 10.0, n=400,
                     t_eval=np.linspace(0, 10, 10001))
    return float(r.t[np.argmax(r.s < 0.5 * R0)])


def run(args):
    geom, eps, t_end = args
    c = case(geom)
    lam = float(fb_slope(W0))
    lim = limit_solution(geom, t_end)
    t0 = time.time()
    r = simulate(c["u0"], W0, eps, c["L"], t_end, N=c["N"], speed_guess=0.5,
                 snapshot_times=(t_end,), max_dt_remesh=0.01 * t_end,
                 mesh=MeshParams(h_coarse=0.01))
    m = (r.t >= 0.05 * t_end) & np.isfinite(r.front)
    s = np.interp(r.t[m], lim.t, lim.s)
    d = r.front[m] - s
    x, u, _ = r.snapshots[t_end]
    k = np.argmin(np.abs(lim.t - t_end))
    err = np.abs(np.maximum(u, 0.0) - lim.u_on(x, k))
    return dict(geom=geom, eps=eps, t_end=t_end,
                front_err=float(np.max(np.abs(d))),
                front_err_shifted=float(np.max(np.abs(d - c["sign"] * eps / lam))),
                sup_err=float(np.max(err)),
                t=r.t[m].tolist(), front=r.front[m].tolist(),
                lim_t=lim.t.tolist(), lim_s=lim.s.tolist(),
                steps=r.n_steps, seconds=time.time() - t0)


def main():
    OUT.mkdir(exist_ok=True)
    geoms = ("planar", "radial2", "radial3")
    t_ends = {g: t_end_of(g) for g in geoms}
    print("t_end:", t_ends, flush=True)
    jobs = [(g, eps, t_ends[g]) for eps in EPS for g in geoms]
    res = []
    with ProcessPoolExecutor(4) as ex:
        for r in ex.map(run, jobs):
            print("{geom:<8} eps={eps:<6} front={front_err:.3e} shifted={front_err_shifted:.3e} "
                  "sup={sup_err:.3e} steps={steps} {seconds:.0f}s".format(**r), flush=True)
            res.append(r)
    (OUT / "limit_comparison.json").write_text(json.dumps(res))


if __name__ == "__main__":
    main()
