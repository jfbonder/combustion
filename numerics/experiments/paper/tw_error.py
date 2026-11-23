"""Group D: the error that the theorems bound, for the planar travelling wave.

Initial datum: the exact eps-travelling wave, U^0_j = Phi_eps(x_j - x0)
(`exact.tw_profile_eps`, normalised Phi_eps(0) = eps), w = eps*w0. We record
at every step (ALE) or every few steps (Eulerian)

    E      = max_j |U^n_j - Phi_eps(x^n_j - x0 - c_eps t_n)|,
    s_n    = x_f^n - x0 - c_eps t_n        (front position error, x_f: {U = eps}),
    E_s    = max_j |U^n_j - Phi_eps(x^n_j - x0 - c_eps t_n - s_n)|  (shift removed),
    v_n - c_eps                            (ALE: realised mesh-velocity error),

and report max over t_n <= T of each, mu = max_n |v_n - c_eps| (all steps and
t_n >= 1), and two speed estimators from s: slope of s on [T/2, T] (least
squares) and s(T)/T.

D1  ALE: vs dt (0.1 ... 0.0025), eps in {0.01, 0.003, 0.001}, w0 in {0, 1}
D2  ALE: vs h_f/eps (1/2 ... 1/32; h_c = 0.0025) and vs h_c (h_f = eps/32)
D3  ALE: speed error vs the lattice phase of the front (the "floor")
D4  rigid mesh with prescribed v = c_eps + dv: drag ratio (c_h - c_eps)/dv
D5  Eulerian scheme (S) (simple.py, exact nominal dt): vs dt/eps and h/eps
F   the old "floor" setup (B1 data, T = 12) with three speed estimators

Run from numerics/:
    OMP_NUM_THREADS=1 ../.venv/bin/python -m experiments.paper.tw_error [run|plot] [groups]
"""

import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from combustion1d.ale import simulate_ale
from combustion1d.exact import tw_profile_eps, tw_profile_limit, tw_speed_eps, tw_speed_limit
from combustion1d.simple import simulate_simple
from combustion1d.solver import MeshParams, find_front, front_speed

from .common import (BLUES, MARKERS, MUTED, SERIES, Timer, cycle_speed, load_json, plain_ticks,
                     plt, save_json, savefig, simulate_ale_safe, slope_ref, title)

X0 = 2.0
THETA = 0.37          # lattice phase of the initial front: x0 = (k + THETA) h_f, k h_f ~ X0
T_ALE = {0.0: 8.0, 1.0: 4.0}      # front travels 4.6 and 6.5
T_EUL = {0.0: 4.0, 1.0: 2.5}
EPS = (0.01, 0.003, 0.001)
DTS = (0.1, 0.05, 0.02, 0.01, 0.005, 0.0025)


def setup(w0, eps, T, pts, theta=THETA):
    c = tw_speed_eps(1.0, w0, eps)
    h_f = eps / pts
    x0 = (round(X0 / h_f) + theta) * h_f     # the band lattice is h_f Z (anchor 0)
    L = x0 + c * T + 12.0 / c + 2.0
    P = lambda z: tw_profile_eps(z, 1.0, w0, eps)
    return c, x0, L, P


class Monitor:
    def __init__(self, P, x0, c, eps, stride=1):
        self.P, self.x0, self.c, self.eps, self.stride = P, x0, c, eps, stride
        self.k = 0
        self.rows = []          # t, E, E_s, s
        self.v = []             # t, v - c

    def __call__(self, t, x, u, v):
        if v is not None:
            self.v.append((t, v - self.c))
        self.k += 1
        if (self.k - 1) % self.stride and t > 0:
            return
        xe = self.x0 + self.c * t
        xf, _ = find_front(x, u, self.eps)
        s = np.nan if xf is None else xf - xe
        E = np.max(np.abs(u - self.P(x - xe)))
        Es = np.max(np.abs(u - self.P(x - xe - s))) if xf is not None else np.nan
        self.rows.append((t, E, Es, s))

    def summary(self, T):
        a = np.array(self.rows)
        t, s = a[:, 0], a[:, 3]
        m = t >= 0.5 * T - 1e-12
        out = dict(E=np.nanmax(a[:, 1]), E_shift=np.nanmax(a[:, 2]), s_max=np.nanmax(np.abs(s)),
                   s_end=s[-1], speed_end=s[-1] / t[-1], speed_lsq=np.polyfit(t[m], s[m], 1)[0],
                   E_end=a[-1, 1], E_shift_end=a[-1, 2])
        if self.v:
            v = np.array(self.v)
            out.update(mu=np.max(np.abs(v[:, 1])), mu1=np.max(np.abs(v[v[:, 0] >= 1.0, 1])),
                       mu_first=abs(v[0, 1]), v_mean2=np.mean(v[v[:, 0] >= 0.5 * T, 1]))
        # coarse time series for plots
        idx = np.unique(np.linspace(0, len(a) - 1, 160).astype(int))
        out["series"] = dict(t=a[idx, 0], E=a[idx, 1], Es=a[idx, 2], s=a[idx, 3])
        return out


def run_ale(job):
    """kind: 'ale' (implicit velocity, remeshing), 'drag' (rigid mesh, v = c + dv)."""
    w0, eps, dt = job["w0"], job["eps"], job["dt"]
    mp = MeshParams(**job.get("mesh", {}))
    T = job.get("T", T_ALE[w0])
    c, x0, L, P = setup(w0, eps, T, mp.pts_per_eps, job.get("theta", THETA))
    mon = Monitor(P, x0, c, eps)
    kw = dict(job.get("kw", {}))
    if job["kind"] == "drag":
        kw.update(fixed_velocity=c + job["dv"], ell0=1e6, remesh_distance=1e9)
    elif job["kind"] == "rigid_impl":
        kw.update(ell0=1e6, remesh_distance=1e9)
    with Timer() as tm:
        r = simulate_ale(lambda x: P(x - x0), w0, eps, L, T, dt, mesh=mp, speed_guess=c,
                         exact_initial=True, monitor=mon, **kw)
    out = dict(job, T=T, c_eps=c, steps=r.n_steps, n_remesh=r.n_remesh,
               trials_per_step=r.n_trials / r.n_steps, solves_per_step=r.n_solves / r.n_steps,
               nodes=float(np.mean(r.n_nodes)), seconds=tm.seconds)
    out.update(mon.summary(T))
    if job["kind"] == "drag":
        out["ratio"] = out["speed_lsq"] / job["dv"]
    return out


def run_simple(job):
    w0, eps, dt = job["w0"], job["eps"], job["dt"]
    mp = MeshParams(**job.get("mesh", {}))
    T = job.get("T", T_EUL[w0])
    c, x0, L, P = setup(w0, eps, T, mp.pts_per_eps)
    mon = Monitor(P, x0, c, eps, stride=max(1, int(round(0.005 / dt))))
    with Timer() as tm:
        r = simulate_simple(lambda x: P(x - x0), w0, eps, L, T, dt=dt, mesh=mp, speed_guess=c,
                            exact_initial=True, monitor=mon)
    out = dict(job, T=T, c_eps=c, steps=r.n_steps, n_remesh=r.n_remesh, dts=r.dts,
               nodes=float(np.mean(r.n_nodes)), seconds=tm.seconds)
    out.update(mon.summary(T))
    return out


def run_floor(job):
    """The setup of the old B1 table (data A(1 - e^{-c_eps(x-2)})^+, T = 12),
    speed measured with three estimators."""
    w0, eps, dt = job["w0"], job["eps"], job["dt"]
    T = 12.0
    c = tw_speed_eps(1.0, w0, eps)
    c0 = float(tw_speed_limit(1.0, w0))
    L = 2.0 + c0 * T + 12.0 / c0 + 2.0
    u0 = lambda x: tw_profile_limit(x, 0.0, 1.0, c, x0=2.0)
    kw = dict(job.get("kw", {}))
    if job["kind"] == "floor_rigid":
        kw.update(ell0=1e6, remesh_distance=1e9)
    with Timer() as tm:
        r, fail = simulate_ale_safe(u0, w0, eps, L, T, dt, mesh=MeshParams(), speed_guess=c0, **kw)
    if fail:
        raise RuntimeError(fail)
    cyc, ncyc = cycle_speed(r, 0.5 * T)
    i2 = np.searchsorted(r.t, 0.5 * T - 1e-9)
    v = np.array(r.v)
    return dict(job, T=T, c_eps=c, err_cycle=cyc - c if ncyc else None, cycles=ncyc,
                err_lsq=front_speed(r, 0.5 * T) - c,
                err_end=(r.front[-1] - r.front[i2]) / (r.t[-1] - r.t[i2]) - c,
                err_mesh=float(np.mean(v[len(v) // 2:])) - c,
                mu1=float(np.max(np.abs(v[int(1.0 / dt):] - c))),
                n_remesh=r.n_remesh, seconds=tm.seconds)


# ---------------------------------------------------------------------------

def jobs(groups):
    J = []
    if "D1" in groups:
        J += [dict(group="D1", kind="ale", w0=w0, eps=e, dt=dt)
              for w0 in (0.0, 1.0) for e in EPS for dt in DTS]
        J += [dict(group="D1", kind="ale", w0=1.0, eps=e, dt=dt, kw=dict(ell0=4.0), variant="ell0=4")
              for e in EPS for dt in DTS]
        J += [dict(group="D1", kind="ale", w0=1.0, eps=e, dt=dt, kw=dict(rigid_extra=0.5),
                   variant="old rigid zone") for e in EPS for dt in DTS]
    if "D2" in groups:
        J += [dict(group="D2h", kind="ale", w0=w0, eps=e, dt=0.01,
                   mesh=dict(pts_per_eps=p, h_coarse=0.0025))
              for w0 in (0.0, 1.0) for e in EPS for p in (2, 4, 8, 16, 32)]
        J += [dict(group="D2c", kind="ale", w0=w0, eps=e, dt=0.01,
                   mesh=dict(pts_per_eps=32, h_coarse=hc))
              for w0 in (0.0, 1.0) for e in EPS for hc in (0.04, 0.02, 0.01, 0.005, 0.0025)]
    if "D3" in groups:
        J += [dict(group="D3", kind="ale", w0=w0, eps=0.01, dt=0.02, theta=th,
                   T={0.0: 4.0, 1.0: 2.0}[w0], mesh=dict(pts_per_eps=p))
              for w0 in (0.0, 1.0) for p in (8, 16, 32) for th in np.arange(16) / 16 + 1 / 32]
        J += [dict(group="D3", kind="ale", w0=w0, eps=0.001, dt=0.02, theta=th,
                   T={0.0: 4.0, 1.0: 2.0}[w0], mesh=dict(pts_per_eps=16))
              for w0 in (0.0, 1.0) for th in np.arange(8) / 8 + 1 / 16]
    if "D4" in groups:
        for w0 in (0.0, 1.0):
            for e in EPS:
                for r in (0.03, 0.1, 0.3, 1, 3, 10, 30, 100):
                    dt = r * e
                    if 1e-4 <= dt <= 0.1 + 1e-12:
                        J.append(dict(group="D4", kind="drag", w0=w0, eps=e, dt=dt, dv=1e-3,
                                      T={0.0: 6.0, 1.0: 3.0}[w0]))
    if "D5" in groups:
        for w0 in (0.0, 1.0):
            for e in EPS:
                for r in (1, 0.3, 0.1, 0.03, 0.01):
                    if r * e >= 1e-5 - 1e-15:
                        J.append(dict(group="D5", kind="simple", w0=w0, eps=e, dt=r * e))
        J += [dict(group="D5h", kind="simple", w0=0.0, eps=e, dt=0.01 * e,
                   mesh=dict(pts_per_eps=p)) for e in (0.01, 0.003) for p in (4, 8, 16, 32)]
    if "F" in groups:
        J += [dict(group="F", kind=k, w0=0.0, eps=e, dt=dt)
              for k in ("floor", "floor_rigid") for e in (0.01, 0.001) for dt in (0.02, 0.005)]
    return J


def _dispatch(job):
    try:
        if job["kind"] == "simple":
            return run_simple(job)
        if job["kind"].startswith("floor"):
            return run_floor(job)
        return run_ale(job)
    except Exception as ex:
        return dict(job, failed=repr(ex))


def _cost(j):
    if j["kind"] == "simple":
        return j.get("T", T_EUL[j["w0"]]) / j["dt"] * 0.3
    return j.get("T", T_ALE[j["w0"]]) / j["dt"] * (2 if j["kind"] == "drag" else 3)


def run_all(groups):
    J = sorted(jobs(groups), key=_cost, reverse=True)
    old = []
    try:
        old = [r for r in load_json("D_tw_error.json") if r["group"].rstrip("hc") not in groups
               and not (r["group"] in ("D2h", "D2c") and "D2" in groups)]
    except FileNotFoundError:
        pass
    out = []
    with Timer() as tm, ProcessPoolExecutor(4) as ex:
        for r in ex.map(_dispatch, J):
            if "failed" in r:
                print("FAILED", r, flush=True)
            elif r["group"] == "F":
                print("F {kind:<12} eps={eps:<6} dt={dt:<6} cycle={err_cycle} lsq={err_lsq:+.2e} "
                      "end={err_end:+.2e} mesh={err_mesh:+.2e} mu1={mu1:.2e}".format(**r), flush=True)
            else:
                try:
                    q = {k: v for k, v in r.items() if k not in ("mu1", "mesh", "kw")}
                    print("{group:<4} {kind:<6} w0={w0} eps={eps:<6} dt={dt:<8.3g} {m} {kwd} "
                          "E={E:.2e} Es={E_shift:.2e} smax={s_max:.2e} lsq={speed_lsq:+.2e} "
                          "mu1={mu1s} {seconds:.0f}s".format(
                              m=r.get("mesh", ""), kwd=r.get("kw", ""),
                              mu1s=f"{r['mu1']:.2e}" if "mu1" in r else "-", **q), flush=True)
                except Exception as e:      # never lose a result because of printing
                    print("print failed", repr(e), flush=True)
            out.append(r)
            save_json("D_tw_error.json", old + out)
    print(f"total {tm.seconds:.0f}s")


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------

def _sel(d, **kw):
    return [r for r in d if all(r.get(k) == v for k, v in kw.items())]


def _w(ax, w0):
    title(ax, f"$w_0={w0:g}$")


def fig_D1(d):
    fig, axs = plt.subplots(2, 3, figsize=(11.2, 6.4), sharex=True)
    for i, w0 in enumerate((0.0, 1.0)):
        for k, e in enumerate(EPS):
            rs = sorted([r for r in _sel(d, group="D1", w0=w0, eps=e) if "variant" not in r],
                        key=lambda r: r["dt"])
            dt = [r["dt"] for r in rs]
            kw = dict(color=BLUES[k + 1], marker=MARKERS[k], ms=5, mfc="white", mew=1.3,
                      label=f"$\\varepsilon={e:g}$")
            axs[i, 0].loglog(dt, [r["E"] for r in rs], **kw)
            axs[i, 1].loglog(dt, [r["E_shift"] for r in rs], **kw)
            axs[i, 2].loglog(dt, [r["mu1"] for r in rs], **kw)
            if w0 == 1.0:
                for var, ls in (("old rigid zone", ":"), ("ell0=4", "--")):
                    rv = sorted(_sel(d, group="D1", w0=w0, eps=e, variant=var), key=lambda r: r["dt"])
                    if e == 0.001:
                        lab = {"ell0=4": "$\\ell_0=4$", "old rigid zone": "old rigid zone"}[var]
                    else:
                        lab = None
                    axs[i, 0].loglog([r["dt"] for r in rv], [r["E"] for r in rv], ls=ls,
                                     color=BLUES[k + 1], lw=1, label=lab)
                    axs[i, 2].loglog([r["dt"] for r in rv], [r["mu1"] for r in rv], ls=ls,
                                     color=BLUES[k + 1], lw=1)
        axs[i, 0].set_ylabel(f"$w_0={w0:g}$")
    for ax in axs[1]:
        ax.set_xlabel("$\\Delta t$")
        plain_ticks(ax, [0.0025, 0.005, 0.01, 0.02, 0.05, 0.1])
    title(axs[0, 0], "$E=\\max_{n,j}|U^n_j-\\Phi_\\varepsilon(x^n_j-c_\\varepsilon t_n)|$")
    title(axs[0, 1], "$E$ with the front shift removed")
    title(axs[0, 2], "$\\mu_1=\\max_{t_n\\geq1}|v_n-c_\\varepsilon|$")
    axs[0, 0].legend(fontsize=7, loc="upper left")
    axs[1, 0].legend(fontsize=7, loc="upper left", ncol=2)
    for ax in axs.flat:
        ax.set_ylim(1e-5, 1e-1)
    fig.tight_layout()
    savefig(fig, "D1_ale_error_vs_dt.png")


def fig_D2(d):
    fig, axs = plt.subplots(1, 3, figsize=(11.2, 3.6))
    for j, w0 in enumerate((0.0, 1.0)):
        for k, e in enumerate(EPS):
            rs = sorted(_sel(d, group="D2h", w0=w0, eps=e), key=lambda r: r["mesh"]["pts_per_eps"])
            x = [1 / r["mesh"]["pts_per_eps"] for r in rs]
            kw = dict(color=SERIES[j], marker=MARKERS[k], ms=5, mfc="white", mew=1.3,
                      ls="-" if k == 0 else ("--" if k == 1 else ":"),
                      label=f"$w_0={w0:g}$, $\\varepsilon={e:g}$")
            axs[0].loglog(x, [r["E_shift"] for r in rs], **kw)
            axs[1].loglog(x, [abs(r["speed_lsq"]) for r in rs], **kw)
            rs = sorted(_sel(d, group="D2c", w0=w0, eps=e), key=lambda r: r["mesh"]["h_coarse"])
            axs[2].loglog([r["mesh"]["h_coarse"] for r in rs], [r["E_shift"] for r in rs], **kw)
    for ax in axs[:2]:
        ax.set_xticks([1 / 32, 1 / 16, 1 / 8, 1 / 4, 1 / 2])
        ax.set_xticklabels(["1/32", "1/16", "1/8", "1/4", "1/2"])
        ax.minorticks_off()
        ax.set_xlabel("$h_f/\\varepsilon$")
    slope_ref(axs[0], [1 / 16, 1 / 2], 2e-5, label="slope 2", p=2)
    slope_ref(axs[1], [1 / 16, 1 / 2], 4e-5, label="slope 2", p=2)
    plain_ticks(axs[2], [0.0025, 0.005, 0.01, 0.02, 0.04], fmt="{:g}")
    axs[2].set_xlabel("$h_c$")
    title(axs[0], "$E$, shift removed ($h_c=0.0025$)")
    title(axs[1], "front speed error $|ds/dt|$ ($h_c=0.0025$)")
    title(axs[2], "$E$, shift removed ($h_f=\\varepsilon/32$)")
    axs[0].legend(fontsize=6.5, loc="upper left", ncol=1)
    fig.tight_layout()
    savefig(fig, "D2_ale_error_vs_mesh.png")


def fig_D3(d):
    fig, axs = plt.subplots(1, 2, figsize=(8.4, 3.4), sharey=False)
    for ax, w0 in zip(axs, (0.0, 1.0)):
        for k, p in enumerate((8, 16, 32)):
            rs = sorted([r for r in _sel(d, group="D3", w0=w0, eps=0.01)
                         if r["mesh"]["pts_per_eps"] == p], key=lambda r: r["theta"])
            ax.plot([r["theta"] for r in rs], [r["v_mean2"] for r in rs], color=BLUES[k + 1],
                    marker="o", ms=4, mfc="white", mew=1.2,
                    label=f"$h_f=\\varepsilon/{p}$, $\\varepsilon=0.01$")
        rs = sorted([r for r in _sel(d, group="D3", w0=w0, eps=0.001)], key=lambda r: r["theta"])
        ax.plot([r["theta"] for r in rs], [r["v_mean2"] for r in rs], ls="none", color=SERIES[1],
                marker="x", ms=6, label="$h_f=\\varepsilon/16$, $\\varepsilon=0.001$")
        ax.axhline(0, color=MUTED, lw=0.8)
        ax.set_xlabel("phase $\\theta$ of the front in the band lattice ($x_f=x_k+\\theta h_f$)")
        title(ax, f"$w_0={w0:g}$: mean $v_n-c_\\varepsilon$, $t\\in[T/2,T]$")
        ax.set_xlim(0, 1)
    axs[0].legend(fontsize=7)
    axs[0].set_ylabel("speed error")
    fig.tight_layout()
    savefig(fig, "D3_phase_pinning.png")


def fig_D4(d):
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    for j, w0 in enumerate((0.0, 1.0)):
        for k, e in enumerate(EPS):
            rs = sorted(_sel(d, group="D4", w0=w0, eps=e), key=lambda r: r["dt"])
            ax.semilogx([r["dt"] / e for r in rs], [r["ratio"] for r in rs], color=SERIES[j],
                        marker=MARKERS[k], ms=5, mfc="white", mew=1.3,
                        ls="-" if k == 0 else ("--" if k == 1 else ":"),
                        label=f"$w_0={w0:g}$, $\\varepsilon={e:g}$")
    ax.set_xlabel("$\\Delta t/\\varepsilon$")
    ax.set_ylabel("$(c_h-c_\\varepsilon)/\\delta v$")
    ax.set_ylim(-0.05, 1.05)
    title(ax, "rigid mesh, $v=c_\\varepsilon+\\delta v$, $\\delta v=10^{-3}$")
    ax.legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    savefig(fig, "D4_drag.png")


def fig_D5(d):
    fig, axs = plt.subplots(1, 3, figsize=(11.2, 3.6))
    for j, w0 in enumerate((0.0, 1.0)):
        for k, e in enumerate(EPS):
            rs = sorted(_sel(d, group="D5", w0=w0, eps=e), key=lambda r: r["dt"])
            x = [r["dt"] / e for r in rs]
            kw = dict(color=SERIES[j], marker=MARKERS[k], ms=5, mfc="white", mew=1.3,
                      ls="-" if k == 0 else ("--" if k == 1 else ":"),
                      label=f"$w_0={w0:g}$, $\\varepsilon={e:g}$")
            axs[0].loglog(x, [abs(r["speed_lsq"]) / r["c_eps"] for r in rs], **kw)
            axs[1].loglog(x, [r["E_shift"] for r in rs], **kw)
    for ax in axs[:2]:
        ax.set_xlabel("$\\Delta t/\\varepsilon$")
    slope_ref(axs[0], [0.01, 0.1], 2e-3, label="slope 1")
    title(axs[0], "relative speed error $|ds/dt|/c_\\varepsilon$")
    title(axs[1], "$E$, shift removed")
    axs[0].legend(fontsize=6.5, loc="upper left")
    for k, e in enumerate((0.01, 0.003)):
        rs = sorted(_sel(d, group="D5h", eps=e), key=lambda r: r["mesh"]["pts_per_eps"])
        x = [1 / r["mesh"]["pts_per_eps"] for r in rs]
        axs[2].loglog(x, [abs(r["speed_lsq"]) for r in rs], color=SERIES[0], marker=MARKERS[k],
                      ms=5, mfc="white", mew=1.3, label=f"speed, $\\varepsilon={e:g}$")
        axs[2].loglog(x, [r["E_shift"] for r in rs], color=SERIES[2], marker=MARKERS[k], ms=5,
                      mfc="white", mew=1.3, ls="--", label=f"$E$ shift removed, $\\varepsilon={e:g}$")
    axs[2].set_xticks([1 / 32, 1 / 16, 1 / 8, 1 / 4])
    axs[2].set_xticklabels(["1/32", "1/16", "1/8", "1/4"])
    axs[2].minorticks_off()
    axs[2].set_xlabel("$h_f/\\varepsilon$")
    axs[2].set_ylim(1e-3, 4e-3)
    plain_ticks(axs[2], [1e-3, 1.5e-3, 2e-3, 3e-3, 4e-3], axis="y", fmt="{:.1e}")
    title(axs[2], "$w_0=0$, $\\Delta t=0.01\\varepsilon$")
    axs[2].legend(fontsize=6.5)
    fig.tight_layout()
    savefig(fig, "D5_eulerian_error.png")


def plot_all():
    d = [r for r in load_json("D_tw_error.json") if "failed" not in r]
    for f in (fig_D1, fig_D2, fig_D3, fig_D4, fig_D5):
        try:
            f(d)
        except Exception as ex:          # group not run yet
            print("skip", f.__name__, repr(ex))


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    groups = sys.argv[2:] or ["D1", "D2", "D3", "D4", "D5", "F"]
    if what in ("run", "all"):
        run_all(groups)
    if what in ("plot", "all"):
        plot_all()
