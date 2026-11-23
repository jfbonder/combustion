"""Group B: convergence (verification of the error estimates).

B1  travelling wave speed error vs eps: Eulerian monotone scheme (simple.py,
    fixed dt) vs ALE, dt in {0.02, 0.01, 0.005}
B2  ALE: error vs h_f = eps/pts_per_eps and vs h_coarse
B3  rigid mesh at the exact speed c_eps: error vs dt (flat)
B4  radial extinction (N = 2, 3) with ALE: shifted front error vs dt

Run from numerics/:
    OMP_NUM_THREADS=1 ../.venv/bin/python -m experiments.paper.convergence [run|plot]
"""

import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from combustion1d.ale import simulate_ale
from combustion1d.exact import tw_speed_eps
from combustion1d.model import fb_slope
from combustion1d.simple import simulate_simple
from combustion1d.solver import MeshParams, front_speed

from .common import (INK, MARKERS, MUTED, SERIES, BLUES, Timer, cycle_speed, extinction_setup,
                     load_json, plt, save_json, savefig, simulate_ale_safe, slope_ref, title,
                     tw_setup)
from combustion1d.exact import tw_profile_limit

T_TW = {0.0: 12.0, 1.0: 8.0}      # final time (front travels ~7 and ~13)
T_RIGID = {0.0: 10.0, 1.0: 3.6}   # rigid mesh: front travels ~6 (see report: longer
                                  # runs collapse the cells at the tapered end x = L)
# width of the decay zone of the mesh velocity: the code default max(2, 3/c)
# gives an O(dt) error for the fast wave w0 = 1 (c = 1.63); we use 4 there
ELL0 = {0.0: None, 1.0: 4.0}
EPS = (0.01, 0.005, 0.002, 0.001)
DTS = (0.02, 0.01, 0.005)
DTS_SIMPLE = (0.02, 0.01, 0.005, 0.002, 0.001)


def _speeds(r, T):
    return front_speed(r, 0.5 * T), front_speed(r, 0.75 * T)


def run_tw(job):
    """job = dict(kind, w0, eps, dt, mesh=dict(...))
    kind: 'simple' (Eulerian monotone scheme), 'ale', 'ale_default' (code
    default ell0), 'rigid' (ALE with the whole mesh translating at c_eps)."""
    w0, eps, dt = job["w0"], job["eps"], job["dt"]
    T = T_RIGID[w0] if job["kind"] == "rigid" else T_TW[w0]
    _, L, c0 = tw_setup(w0, T)
    c_eps = tw_speed_eps(1.0, w0, eps)
    # initial datum: limit profile with the speed c_eps, i.e. exactly the fresh
    # part of the (P_eps) travelling wave (removes a slowly decaying transient)
    u0 = lambda x: tw_profile_limit(x, 0.0, 1.0, c_eps, x0=2.0)
    mp = MeshParams(**job.get("mesh", {}))
    cycles = None
    with Timer() as tm:
        if job["kind"] == "simple":
            r = simulate_simple(u0, w0, eps, L, T, N=1, dt=dt, mesh=mp, speed_guess=c0)
            c_num = front_speed(r, 0.5 * T)
        elif job["kind"] == "rigid":
            r = simulate_ale(u0, w0, eps, L, T, dt, mesh=mp, speed_guess=c_eps,
                             fixed_velocity=c_eps, ell0=1e6, remesh_distance=1e9)
            c_num = front_speed(r, 0.5 * T)
        else:
            ell0 = None if job["kind"] == "ale_default" else ELL0[w0]
            r, fail = simulate_ale_safe(u0, w0, eps, L, T, dt, mesh=mp, speed_guess=c0,
                                        ell0=ell0)
            if fail:
                raise RuntimeError(fail)
            c_num, cycles = cycle_speed(r, 0.5 * T)
    c_ls = front_speed(r, 0.5 * T)
    return dict(job, T=T, c0=c0, c_eps=c_eps, c_num=c_num, err=c_num - c_eps,
                err_lsq=c_ls - c_eps, cycles=cycles, steps=r.n_steps,
                n_remesh=r.n_remesh, nodes=float(np.mean(r.n_nodes)), seconds=tm.seconds)


def run_radial(job):
    N, eps, dt = job["N"], job["eps"], job["dt"]
    w0, t_end = 0.5, 0.75
    lam = float(fb_slope(w0))
    V, u0, L = extinction_setup(N, w0)
    with Timer() as tm:
        r = simulate_ale(u0, w0, eps, L, t_end, dt, mesh=MeshParams(), speed_guess=-V.R / 2,
                         N=N, ell0=1.0, remesh_distance=0.1, snapshot_times=(t_end,))
    m = r.t >= 0.05
    d = r.front[m] - V.radius(r.t[m])
    x, u, _ = r.snapshots[t_end]
    err = np.abs(np.maximum(u, 0.0) - V(x, t_end))
    interior = x < V.radius(t_end) - 20 * eps
    return dict(job, dt_actual=t_end / r.n_steps, front_err=float(np.max(np.abs(d))),
                front_err_shifted=float(np.max(np.abs(d + eps / lam))),
                front_err_shifted_end=float(abs(d[-1] + eps / lam)),
                sup_err=float(np.max(err)), sup_err_int=float(np.max(err[interior])),
                steps=r.n_steps, n_remesh=r.n_remesh, nodes=float(np.mean(r.n_nodes)),
                seconds=tm.seconds)


def jobs_B():
    B1 = [dict(group="B1", kind="ale", w0=w0, eps=e, dt=dt)
          for w0 in (0.0, 1.0) for e in EPS for dt in DTS]
    B1 += [dict(group="B1", kind="ale_default", w0=1.0, eps=e, dt=dt) for e in EPS for dt in DTS]
    B1 += [dict(group="B1", kind="simple", w0=w0, eps=e, dt=dt)
           for w0 in (0.0, 1.0) for e in EPS for dt in DTS_SIMPLE]
    B2 = [dict(group="B2h", kind="ale", w0=0.0, eps=e, dt=0.01,
               mesh=dict(pts_per_eps=p, h_coarse=0.0025))
          for e in EPS for p in (4, 8, 16, 32, 64)]
    B2 += [dict(group="B2c", kind="ale", w0=0.0, eps=e, dt=0.01,
                mesh=dict(pts_per_eps=32, h_coarse=hc))
           for e in (0.01, 0.001) for hc in (0.04, 0.02, 0.01, 0.005, 0.0025)]
    B3 = [dict(group="B3", kind="rigid", w0=w0, eps=e, dt=dt)
          for w0 in (0.0, 1.0) for e in (0.005, 0.001)
          for dt in (0.1, 0.05, 0.02, 0.01, 0.005)]
    B4 = [dict(group="B4", N=N, eps=e, dt=dt) for N in (2, 3) for e in (0.01, 0.002, 0.001)
          for dt in (0.04, 0.02, 0.01, 0.005, 0.0025)]
    return B1 + B2 + B3, B4


def _dispatch(job):
    try:
        return run_radial(job) if job["group"] == "B4" else run_tw(job)
    except Exception as ex:
        return dict(job, failed=repr(ex))


def run_all():
    tw_jobs, rad_jobs = jobs_B()
    jobs = sorted(tw_jobs + rad_jobs, key=lambda j: -j.get("dt", 0) * 0 + j["eps"])  # small eps first
    out = []
    with Timer() as tm, ProcessPoolExecutor(4) as ex:
        for r in ex.map(_dispatch, jobs):
            if "failed" in r:
                print("FAILED", r, flush=True)
            elif r["group"] == "B4":
                print("B4 N={N} eps={eps:<6} dt={dt:<7} shifted={front_err_shifted:.3e} "
                      "sup={sup_err:.3e} {seconds:.1f}s".format(**r), flush=True)
            else:
                print("{group:<4} {kind:<6} w0={w0} eps={eps:<6} dt={dt:<6} mesh={m} err={err:+.3e} "
                      "lsq={err_lsq:+.3e} cyc={cycles} steps={steps} {seconds:.1f}s".format(
                          m=r.get("mesh", {}), **r), flush=True)
            out.append(r)
    print(f"total {tm.seconds:.0f}s")
    save_json("B_convergence.json", out)


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------

def _sel(d, **kw):
    return [r for r in d if all(r.get(k) == v for k, v in kw.items())]


def fig_B1(d):
    fig, axs = plt.subplots(1, 3, figsize=(11.2, 3.6))
    for ax, w0 in zip(axs[:2], (0.0, 1.0)):
        for k, dt in enumerate(DTS):
            for kind, ls, mk in (("simple", "--", "s"), ("ale", "-", "o")):
                rs = sorted(_sel(d, group="B1", kind=kind, w0=w0, dt=dt), key=lambda r: r["eps"])
                e = np.array([r["eps"] for r in rs])
                ax.loglog(e, [abs(r["err"]) for r in rs], color=BLUES[k + 1], ls=ls, marker=mk,
                          ms=5, mfc="white", mew=1.3,
                          label=f"{'Eulerian' if kind == 'simple' else 'ALE'}, $\\Delta t={dt:g}$")
        ax.set_xlabel("$\\varepsilon$")
        title(ax, f"$w_0={w0:g}$: $|c_h-c_\\varepsilon|$")
        ax.set_ylim(3e-6, 2)
    slope_ref(axs[0], [0.002, 0.01], 0.012, label="$\\varepsilon^{-1}$", p=-1)
    axs[0].set_ylabel("speed error")
    axs[1].legend(loc="center right", fontsize=7, ncol=2, bbox_to_anchor=(1.0, 0.55))
    # collapse of the Eulerian error as a function of dt/eps
    ax = axs[2]
    for k, w0 in enumerate((0.0, 1.0)):
        rs = _sel(d, group="B1", kind="simple", w0=w0)
        x = np.array([r["dt"] / r["eps"] for r in rs])
        y = np.array([abs(r["err"]) / r["c_eps"] for r in rs])
        ax.loglog(x, y, ls="none", marker=MARKERS[k], color=SERIES[k], mfc="white", mew=1.3,
                  ms=5, label=f"Eulerian, $w_0={w0:g}$")
        rs = _sel(d, group="B1", kind="ale", w0=w0)
        x = np.array([r["dt"] / r["eps"] for r in rs])
        y = np.array([abs(r["err"]) / r["c_eps"] for r in rs])
        ax.loglog(x, y, ls="none", marker=MARKERS[k], color=SERIES[k], ms=5,
                  label=f"ALE, $w_0={w0:g}$")
    slope_ref(ax, [0.1, 2], 0.03, label="slope 1")
    ax.set_xlabel("$\\Delta t/\\varepsilon$")
    title(ax, "relative error $|c_h-c_\\varepsilon|/c_\\varepsilon$")
    ax.legend(loc="lower right", fontsize=7)
    fig.tight_layout()
    savefig(fig, "B1_tw_speed_eulerian_vs_ale.png")


def fig_B2(d):
    fig, axs = plt.subplots(1, 2, figsize=(8.4, 3.5))
    ax = axs[0]
    for k, e in enumerate(EPS):
        rs = sorted(_sel(d, group="B2h", eps=e), key=lambda r: r["mesh"]["pts_per_eps"])
        x = np.array([1.0 / r["mesh"]["pts_per_eps"] for r in rs])
        ax.loglog(x, [abs(r["err"]) for r in rs], color=BLUES[k], marker="o", ms=5, mfc="white",
                  mew=1.3, label=f"$\\varepsilon={e:g}$")
    slope_ref(ax, [1 / 32, 1 / 4], 3e-5, label="slope 2", p=2)
    ax.set_xlabel("$h_f/\\varepsilon$")
    ax.set_ylabel("$|c_h-c_\\varepsilon|$")
    title(ax, "band mesh size ($h_c=0.0025$)")
    ax.legend(loc="upper left", fontsize=8)
    ax = axs[1]
    for k, e in enumerate((0.01, 0.001)):
        rs = sorted(_sel(d, group="B2c", eps=e), key=lambda r: r["mesh"]["h_coarse"])
        x = np.array([r["mesh"]["h_coarse"] for r in rs])
        ax.loglog(x, [abs(r["err"]) for r in rs], color=SERIES[k], marker=MARKERS[k], ms=5,
                  mfc="white", mew=1.3, label=f"$\\varepsilon={e:g}$")
    slope_ref(ax, [0.0025, 0.04], 2e-5, label="slope 1")
    ax.set_xlabel("$h_c$")
    title(ax, "coarse mesh size ($h_f=\\varepsilon/32$)")
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    savefig(fig, "B2_ale_mesh_sizes.png")


def fig_B3(d):
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    for k, e in enumerate((0.005, 0.001)):
        for j, w0 in enumerate((0.0, 1.0)):
            rs = sorted(_sel(d, group="B3", eps=e, w0=w0), key=lambda r: r["dt"])
            ax.loglog([r["dt"] for r in rs], [max(abs(r["err"]), 1e-9) for r in rs],
                      color=SERIES[k], ls="-" if j == 0 else "--", marker=MARKERS[k], ms=5,
                      mfc="white" if j == 0 else SERIES[k], mew=1.3,
                      label=f"$\\varepsilon={e:g}$, $w_0={w0:g}$")
    ax.set_xlabel("$\\Delta t$")
    ax.set_ylabel("$|c_h-c_\\varepsilon|$")
    ax.set_ylim(1e-8, 1e-3)
    title(ax, "rigid mesh translating with $c_\\varepsilon$")
    ax.legend(loc="upper left", fontsize=8, ncol=2)
    fig.tight_layout()
    savefig(fig, "B3_rigid_mesh_dt.png")


def fig_B4(d):
    fig, axs = plt.subplots(1, 2, figsize=(8.4, 3.5), sharey=True)
    for ax, N in zip(axs, (2, 3)):
        for k, e in enumerate((0.01, 0.002, 0.001)):
            rs = sorted(_sel(d, group="B4", N=N, eps=e), key=lambda r: r["dt"])
            ax.loglog([r["dt_actual"] for r in rs], [r["front_err_shifted"] for r in rs],
                      color=SERIES[k], marker=MARKERS[k], ms=5, mfc="white", mew=1.3,
                      label=f"$\\varepsilon={e:g}$")
        slope_ref(ax, [0.0025, 0.04], 0.0025, label="slope 1")
        ax.set_xlabel("$\\Delta t$")
        title(ax, f"$N={N}$: $\\max_t|r_h(t)+\\varepsilon/\\sqrt{{2M}}-R\\sqrt{{T-t}}|$")
        ax.legend(loc="lower right", fontsize=8)
    axs[0].set_ylabel("front error")
    fig.tight_layout()
    savefig(fig, "B4_radial_extinction_dt.png")


def plot_all():
    d = load_json("B_convergence.json")
    fig_B1(d)
    fig_B2(d)
    fig_B3(d)
    fig_B4(d)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("run", "all"):
        run_all()
    if what in ("plot", "all"):
        plot_all()
