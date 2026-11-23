"""Group E: slow transients and the accumulated mesh-velocity error.

E1  Eulerian scheme (S) from Phi_eps, long runs (T = 60, w0 = 0): the speed
    error measured on successive time windows, divided by dt/eps, against the
    solvability prediction of experiments/truncation_check.py (-0.235).
    This is Table tab:simple, with the time window made explicit.
E2  ALE from Phi_eps (eps = 0.01, dt = 0.02, T = 24): realised mesh velocity
    v_n - c_eps as a function of time for four lattice phases of the front.
E3  ALE from Phi_eps (T = 8, phase 0.37): for dt in {0.1, 0.02, 0.0025} and
    eps in {0.01, 0.001}, the quantities of the ALE theorem:
        M_T = dt * sum_n |v_n - c_eps|   (accumulated mesh-velocity error),
        max_n |v_n - c_eps|, E, E_s, max |s_n|.
E4  ALE from Phi_eps (T = 8, dt = 0.02): E, E_s and the speed error for 16
    lattice phases (eps = 0.01) and 8 phases (eps = 0.001).

Run from numerics/:
    OMP_NUM_THREADS=1 ../.venv/bin/python -m experiments.paper.transient [run|plot] [groups]
Figure: figures/paper/D6_transients.png.  Results: results/paper/E_transient.json.
Runtime: E1 about 25 min on 4 cores (2.4e6 steps for eps = 0.01, dt = eps^2/4); E2-E4 2 min.
"""

import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from combustion1d.ale import simulate_ale
from combustion1d.simple import simulate_simple
from combustion1d.solver import MeshParams

from .common import BLUES, MARKERS, MUTED, SERIES, Timer, load_json, plt, save_json, savefig, title
from .tw_error import Monitor, setup

WINDOW = 6.0


def run_eulerian(job):
    w0, eps, k, T = job["w0"], job["eps"], job["k"], job["T"]
    dt = k * eps**2
    mp = MeshParams()
    c, x0, L, P = setup(w0, eps, T, mp.pts_per_eps)
    with Timer() as tm:
        r = simulate_simple(lambda x: P(x - x0), w0, eps, L, T, dt=dt, mesh=mp, speed_guess=c,
                            exact_initial=True)
    t, s = r.t, r.front - x0 - c * r.t
    win = []
    a = 0.0
    while a + WINDOW <= T + 1e-9:
        m = (t >= a - 1e-9) & (t <= a + WINDOW + 1e-9)
        slope = np.polyfit(t[m], s[m], 1)[0]
        win.append((a, a + WINDOW, slope, slope / (dt / eps)))
        a += WINDOW
    idx = np.unique(np.linspace(0, len(t) - 1, 400).astype(int))
    return dict(job, dt=dt, c_eps=c, steps=r.n_steps, dts=r.dts, seconds=tm.seconds,
                windows=win, t=t[idx], s=s[idx])


def run_ale(job):
    w0, eps, dt, T = job["w0"], job["eps"], job["dt"], job["T"]
    mp = MeshParams()
    c, x0, L, P = setup(w0, eps, T, mp.pts_per_eps, job["theta"])
    mon = Monitor(P, x0, c, eps)
    with Timer() as tm:
        r = simulate_ale(lambda x: P(x - x0), w0, eps, L, T, dt, mesh=mp, speed_guess=c,
                         exact_initial=True, monitor=mon)
    out = dict(job, c_eps=c, steps=r.n_steps, n_remesh=r.n_remesh, seconds=tm.seconds)
    summ = mon.summary(T)
    summ.pop("series")
    out.update(summ)
    v = np.array(mon.v)                      # (t, v - c)
    out["M_T"] = float(dt * np.sum(np.abs(v[:, 1])))
    out["X_T"] = float(dt * np.sum(v[:, 1]))             # mesh displacement error X_N - c T
    if job["group"] == "E2":
        n = max(1, int(round(0.25 / dt)))                # block means over 0.25 time units
        nb = len(v) // n
        out["vt"] = v[:nb * n, 0].reshape(nb, n).mean(1)
        out["vv"] = v[:nb * n, 1].reshape(nb, n).mean(1)
    return out


def jobs(groups):
    J = []
    if "E1" in groups:
        J += [dict(group="E1", w0=0.0, eps=e, k=k, T=60.0) for e in (0.01, 0.02) for k in (0.25, 1, 4)]
    if "E2" in groups:
        J += [dict(group="E2", w0=0.0, eps=0.01, dt=0.02, T=24.0, theta=th)
              for th in (1 / 32, 0.37, 17 / 32, 25 / 32)]
    if "E3" in groups:
        J += [dict(group="E3", w0=0.0, eps=e, dt=dt, T=8.0, theta=0.37)
              for e in (0.01, 0.001) for dt in (0.1, 0.02, 0.0025)]
    if "E4" in groups:
        J += [dict(group="E4", w0=0.0, eps=0.01, dt=0.02, T=8.0, theta=th)
              for th in np.arange(16) / 16 + 1 / 32]
        J += [dict(group="E4", w0=0.0, eps=0.001, dt=0.02, T=8.0, theta=th)
              for th in np.arange(8) / 8 + 1 / 16]
    return J


def _run(job):
    return run_eulerian(job) if job["group"] == "E1" else run_ale(job)


def run_all(groups):
    try:
        old = [r for r in load_json("E_transient.json") if r["group"] not in groups]
    except FileNotFoundError:
        old = []
    J = jobs(groups)
    J.sort(key=lambda j: -(j["T"] / (j["k"] * j["eps"]**2) if j["group"] == "E1" else 0))
    out = []
    with ProcessPoolExecutor(4) as ex:
        for r in ex.map(_run, J):
            if r["group"] == "E1":
                print("E1 eps={eps} k={k} steps={steps} {seconds:.0f}s ".format(**r)
                      + " ".join("%.4f" % w[3] for w in r["windows"]), flush=True)
            else:
                print("{group} eps={eps} dt={dt} theta={theta:.4f} E={E:.3e} Es={E_shift:.3e} "
                      "speed={speed_lsq:+.3e} mu={mu:.3e} M_T={M_T:.3e} smax={s_max:.3e}".format(**r),
                      flush=True)
            out.append(r)
            save_json("E_transient.json", old + out)
    save_json("E_transient.json", old + out)


def plot_all():
    R = load_json("E_transient.json")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 3.8))
    E1 = sorted([r for r in R if r["group"] == "E1"], key=lambda r: (r["eps"], r["k"]))
    for i, r in enumerate(E1):
        w = np.array(r["windows"])
        ls = "-" if r["eps"] == 0.01 else "--"
        col = BLUES[1 + [0.25, 1, 4].index(r["k"])]
        a1.plot(0.5 * (w[:, 0] + w[:, 1]), w[:, 3], ls, color=col, marker=MARKERS[[0.25, 1, 4].index(r["k"])],
                ms=4, label=r"$\varepsilon=%g$, $\Delta t/\varepsilon=%g$" % (r["eps"], r["k"] * r["eps"]))
    a1.axhline(-0.2351, color=MUTED, lw=1, ls=":")
    a1.text(1.5, -0.2351, "averaged truncation error", color=MUTED, fontsize=7, va="bottom", ha="left")
    a1.set_xlabel("$t$ (midpoint of the window)")
    a1.set_ylabel(r"$(c_h-c_\varepsilon)\,/\,(\Delta t/\varepsilon)$")
    a1.legend(fontsize=7, ncol=2)
    title(a1, "Eulerian (S): speed error on windows of length %g" % WINDOW)
    E2 = sorted([r for r in R if r["group"] == "E2"], key=lambda r: r["theta"])
    for i, r in enumerate(E2):
        a2.plot(r["vt"], np.array(r["vv"]) * 1e4, color=(SERIES + [MUTED])[i],
                label=r"$\vartheta=%s$" % {1 / 32: "1/32", 17 / 32: "17/32", 25 / 32: "25/32"}.get(r["theta"], "%.2f" % r["theta"]))
    a2.axhline(0, color=MUTED, lw=0.8)
    a2.set_xlabel("$t$")
    a2.set_ylabel(r"$(v_n-c_\varepsilon)\cdot10^{4}$")
    a2.set_ylim(-3, 3)
    a2.legend(fontsize=7, ncol=2)
    title(a2, r"ALE: mesh velocity, $\varepsilon=0.01$, $\Delta t=0.02$")
    fig.tight_layout()
    savefig(fig, "D6_transients.png")


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    groups = sys.argv[2:] or ["E1", "E2", "E3", "E4"]
    if what in ("run", "all"):
        run_all(groups)
    if what in ("plot", "all"):
        plot_all()
