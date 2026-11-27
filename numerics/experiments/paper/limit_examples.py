"""Group C: (P_eps) with ALE against the computed limit problem, and two examples
from the literature.

In the article: C1 is Figure 2 (label fig:C1, figures/paper/C1_ale_vs_limit.png,
Section 6.2) and C2 is Figure 6 (label fig:C2,
figures/paper/C2_self_similar_convergence.png, Section 6.7), with the numbers
quoted in those sections. C3 is not in the paper.

C1  ALE (P_eps) vs the computed limit problem, radial data
    u0 = a (1 - r^2/R0^2)^+, a = sqrt(2M) R0/2, R0 = 2, w0 = 0.5 (N = 2, 3)
C2  convergence to the self-similar extinction profile (Galaktionov-Hulshof-
    Vazquez): s(t)/sqrt(T-t) -> R, rescaled profiles -> kappa*phi
C3  quenching of small burnt kernels (expanding flame, tanh data): critical
    initial radius vs eps, front radius r_f(t) for several r0

Output: results/paper/C1_limit_comparison.json, C2_self_similar.json,
C3_quenching.json and the figures C1_ale_vs_limit.png,
C2_self_similar_convergence.png, C3_quenching.png in figures/paper. The grey
points of C1 (BDF) are read from results/limit_comparison.json: run
experiments.limit_comparison first.

Run from numerics/:
    OMP_NUM_THREADS=1 python -m experiments.paper.limit_examples [run|plot] [C1 C2 C3]
(about 13 min with 4 processes: C1 1 min, C2 2 min, C3 10 min. `plot` draws
the figures whose data exist.)

Failures. A numerical failure of a run (Peclet condition, invalid mesh, BDF
integrator) is never used as a result: C2 and C3 stop with an error. In C3 a
run counts as "quenching" only if the front {u = eps} disappeared (status
"extinct" of the solver); the criterion that decided every bisection step is
stored in the history of C3_quenching.json. A bisection that cannot decide
a step is aborted: its steps are written to C3_quenching_aborted.json and
the script stops without touching C3_quenching.json.

Numbers of figures and sections are those of the draft in paper/.
"""

import json
import sys

import numpy as np

from combustion1d.ale import simulate_ale
from combustion1d.exact import SelfSimilar
from combustion1d.limit import solve_radial
from combustion1d.model import fb_slope
from combustion1d.solver import MeshParams, NumericalFailure, find_front, simulate
from experiments.limit_comparison import case, limit_solution, t_end_of

from .common import (BLUES, INK, MARKERS, MUTED, RES, ROOT, SERIES, Timer, expanding_u0,
                     load_json, plt, save_json, savefig, simulate_ale_safe, slope_ref, title,
                     worker_pool)

W0 = 0.5
LAM = float(fb_slope(W0))
C0 = LAM                    # planar speed sqrt(2M)/A, A = 1
R0 = 2.0
EPS_C1 = (0.05, 0.02, 0.01, 0.005, 0.002, 0.001)
DT_C1 = (0.002, 0.0005, 0.000125)


# ---------------------------------------------------------------------------
# C1
# ---------------------------------------------------------------------------

def run_C1(job):
    geom, eps, dt, t_end = job
    c = case(geom)
    N = c["N"]
    lim = limit_solution(geom, t_end)
    with Timer() as tm:
        r = simulate_ale(c["u0"], W0, eps, c["L"], t_end, dt, mesh=MeshParams(),
                         speed_guess=-N / R0, N=N, ell0=1.0, remesh_distance=0.1,
                         snapshot_times=(t_end,))
    m = r.t >= 0.05 * t_end
    s = np.interp(r.t[m], lim.t, lim.s)
    d = r.front[m] - s
    x, u, _ = r.snapshots[t_end]
    k = np.argmin(np.abs(lim.t - t_end))
    err = np.abs(np.maximum(u, 0.0) - lim.u_on(x, k))
    return dict(geom=geom, N=N, eps=eps, dt=dt, t_end=t_end,
                front_err=float(np.max(np.abs(d))),
                front_err_shifted=float(np.max(np.abs(d + eps / LAM))),
                sup_err=float(np.max(err)),
                t=r.t[m][::5], front=r.front[m][::5], lim_t=lim.t, lim_s=lim.s,
                steps=r.n_steps, n_remesh=r.n_remesh, nodes=float(np.mean(r.n_nodes)),
                seconds=tm.seconds)


# ---------------------------------------------------------------------------
# C2
# ---------------------------------------------------------------------------

def extinction_time(t, s, k=6):
    """T from a linear fit of s^2 on the last k points (s^2 ~ R^2 (T - t))."""
    p = np.polyfit(t[-k:], s[-k:] ** 2, 1)
    return -p[1] / p[0]


def run_C2(job):
    N, eps_dt = job
    a = LAM * R0 / 2.0
    u0 = lambda r: a * np.maximum(1.0 - (r / R0) ** 2, 0.0)
    V = SelfSimilar(N, W0)
    lim = solve_radial(u0, W0, N, R0, 10.0, n=1600)
    m = lim.s < 0.3
    T = extinction_time(lim.t[m], lim.s[m], k=int(m.sum()))
    # profile times: radius 1.5, 1.0, 0.5, 0.25
    prof_t = [float(lim.t[np.argmin(np.abs(lim.s - q))]) for q in (1.5, 1.0, 0.5, 0.25)]
    profiles = []
    for tq in prof_t:
        k = int(np.argmin(np.abs(lim.t - tq)))
        tau = T - lim.t[k]
        rr = lim.x_of(k)
        profiles.append(dict(t=lim.t[k], s=lim.s[k], eta=rr / np.sqrt(tau),
                             v=lim.v[:, k] / np.sqrt(tau)))
    out = dict(N=N, R=V.R, kappa=V.kappa, T=T, lim_t=lim.t, lim_s=lim.s, profiles=profiles,
               eta=np.linspace(0, V.R, 200), phi=V.kappa * V.phi(np.linspace(0, V.R, 200)),
               eps_runs=[])
    # short domain: the burnt region r > R0 is irrelevant (u = -eps w0 there), and
    # the ALE mesh velocity is tapered to 0 on r < min(0.5, 0.1 L) (code
    # constant), which adds an O(dt v/eps) lag of the front close to extinction;
    # L = 2.5 makes that zone 0.25 instead of 0.4 (see README)
    L = R0 + 0.5
    for eps, dt in eps_dt:
        # stop when the limit radius is ~ max(0.15, 30 eps)
        s_end = max(0.15, 30 * eps)
        t_end = float(lim.t[np.argmax(lim.s < s_end)])
        snaps = tuple(np.round(np.array(prof_t) / dt) * dt)
        with Timer() as tm:
            r, err = simulate_ale_safe(u0, W0, eps, L, t_end, dt, mesh=MeshParams(),
                                       speed_guess=-N / R0, N=N, ell0=1.0, remesh_distance=0.1,
                                       snapshot_times=snaps)
        if err is not None:
            raise RuntimeError(f"C2 N={N} eps={eps} dt={dt}: ALE run failed: {err}")
        s_eps = r.front + eps / LAM
        # T_eps from s^2 ~ R^2 (T_eps - t) on s in [s_end, 2 s_end]
        m = (s_eps < 2 * s_end) & (s_eps > 0)
        T_eps = extinction_time(r.t[m], s_eps[m], k=int(m.sum()))
        eprof = []
        for tq, sq in zip(snaps, prof_t):
            if tq in r.snapshots:
                x, u, _ = r.snapshots[tq]
                tau = T_eps - tq
                m = x < 1.2 * np.interp(tq, r.t, s_eps)
                eprof.append(dict(t=tq, eta=x[m] / np.sqrt(tau),
                                  v=np.maximum(u[m], 0) / np.sqrt(tau)))
        out["eps_runs"].append(dict(eps=eps, dt=dt, T_eps=T_eps, t=r.t[::4], s=s_eps[::4],
                                    profiles=eprof, steps=r.n_steps, seconds=tm.seconds,
                                    error=err, status=r.status, reason=r.reason,
                                    n_secant_fail=r.n_secant_fail,
                                    n_mesh_remesh=r.n_mesh_remesh))
    # reference for the smallest eps: Eulerian BDF solver (adaptive in time)
    eps = eps_dt[-1][0]
    s_end = max(0.15, 30 * eps)
    t_end = float(lim.t[np.argmax(lim.s < s_end)])
    with Timer() as tm:
        r = simulate(u0, W0, eps, L, t_end, N=N, speed_guess=N / R0, max_dt_remesh=0.002,
                     snapshot_times=tuple(prof_t))
    s_eps = r.front + eps / LAM
    m = (s_eps < 2 * s_end) & (s_eps > 0)
    T_eps = extinction_time(r.t[m], s_eps[m], k=int(m.sum()))
    eprof = []
    for tq in prof_t:
        x, u, _ = r.snapshots[tq]
        tau = T_eps - tq
        mm = x < 1.2 * np.interp(tq, r.t, s_eps)
        eprof.append(dict(t=tq, eta=x[mm] / np.sqrt(tau), v=np.maximum(u[mm], 0) / np.sqrt(tau)))
    out["bdf_run"] = dict(eps=eps, T_eps=T_eps, t=r.t, s=s_eps, profiles=eprof,
                          steps=r.n_steps, seconds=tm.seconds)
    return out


# ---------------------------------------------------------------------------
# C3
# ---------------------------------------------------------------------------

L_EXP = 12.0
T_EXP = 6.0
# A kernel run is "quenching" when the front disappears. With QUENCH_BELOW = q
# (a number in (0, 1)) a run that ends with a front radius below q*r0 is also
# counted as quenching, although the front is still present; the reason stored
# with the step says which criterion decided. None = extinction only.
QUENCH_BELOW = None
# Collapsing kernels enter the zone near r = 0 where the mesh velocity is
# tapered; there the band nodes move more slowly than the front, and for small
# eps (narrow band) the front leaves the fine band and stalls in cells that do
# not resolve the reaction layer (for eps = 0.00125 it stays at r = 0.26
# indefinitely). With REMESH_OFF_BAND the ALE mesh is rebuilt around the front
# when that happens (simulate_ale(remesh_off_band=True)), so that the collapse
# is followed until the front disappears. It has no effect on runs in which the
# front stays in the band (all the igniting kernels). With False, such runs
# are reported as 'unresolved' and the bisection is aborted.
REMESH_OFF_BAND = True


class KernelRunFailed(RuntimeError):
    """A kernel run ended with a numerical failure (not a physical outcome)."""


def run_kernel(job):
    """Expanding-flame run from a kernel of radius r0; returns the front
    history and the outcome:

    'quenching'  the front disappeared (solver status "extinct");
    'ignition'   the run ended with a front radius above r0 + 0.5;
    'undecided'  the run ended with a front and a radius below r0 + 0.5. This
                 includes a kernel that has shrunk but whose front is still
                 present at the end of the run (unless QUENCH_BELOW is set);
    'unresolved' as 'undecided', and the front had left the fine band of the
                 ALE mesh (res.n_front_off_band > 0): the run cannot decide.

    The entry 'reason' says which criterion decided, 'status' is the status
    of the solver. A numerical failure raises KernelRunFailed: it is never
    turned into an outcome.
    """
    r0, eps, dt, method = job["r0"], job["eps"], job["dt"], job.get("method", "ale")
    t_end = job.get("t_end", T_EXP)
    u0 = expanding_u0(r0, W0)
    L = max(L_EXP, r0 + C0 * t_end + 4.0)
    snaps = tuple(np.round(np.arange(0.02, t_end + 1e-9, 0.02), 10))
    counters = {}
    with Timer() as tm:
        if method == "ale":
            r, err = simulate_ale_safe(u0, W0, eps, L, t_end, dt, mesh=MeshParams(),
                                       speed_guess=C0, N=2, snapshot_times=snaps,
                                       remesh_off_band=REMESH_OFF_BAND)
            if err is not None:
                raise KernelRunFailed(f"r0={r0} eps={eps} dt={dt} (ALE): {err}")
            ts, fs = [0.0], [r0 + eps / LAM]
            for s in sorted(r.snapshots):
                x, u, _ = r.snapshots[s]
                xf, _ = find_front(x, u, eps)
                ts.append(s)
                fs.append(np.nan if xf is None else xf)
            ts, fs = np.array(ts), np.array(fs)
            steps = len(r.n_nodes)
            counters = dict(n_secant_fail=r.n_secant_fail, n_mesh_remesh=r.n_mesh_remesh,
                            n_front_off_band=r.n_front_off_band,
                            t_front_off_band=r.t_front_off_band)
        else:
            # (a failure of the BDF integrator raises NumericalFailure)
            r = simulate(u0, W0, eps, L, t_end, N=2, mesh=MeshParams(), speed_guess=C0,
                         max_dt_remesh=0.02)
            ts, fs, steps = r.t, r.front, r.n_steps
    if r.status not in ("ok", "extinct"):
        raise KernelRunFailed(f"r0={r0} eps={eps} dt={dt} ({method}): {r.status}: {r.reason}")
    ok = np.isfinite(fs)
    last = fs[ok][-1]
    if r.status == "extinct":
        outcome, reason = "quenching", r.reason
    elif not np.all(ok):
        # a recorded sample without a front although the solver reports one
        raise KernelRunFailed(f"r0={r0} eps={eps} dt={dt} ({method}): front samples missing "
                              "but the solver did not report extinction")
    elif last > r0 + 0.5:
        outcome, reason = "ignition", f"front radius {last:.6g} > r0 + 0.5 at t = {ts[-1]:.6g}"
    elif QUENCH_BELOW is not None and last < QUENCH_BELOW * r0:
        outcome = "quenching"
        reason = (f"front radius {last:.6g} < {QUENCH_BELOW:g} r0 at t = {ts[-1]:.6g} "
                  "(front still present: no extinction)")
    else:
        outcome = "undecided"
        reason = (f"front radius {last:.6g} <= r0 + 0.5 at t = {ts[-1]:.6g}, front still present "
                  f"(minimum radius {np.nanmin(fs):.6g})")
        if counters.get("n_front_off_band"):
            # the front has left the fine band (tapered zone near r = 0) and sits
            # in cells that do not resolve the reaction layer: it can stall there
            outcome = "unresolved"
            reason += (f"; the front left the fine band at t = {counters['t_front_off_band']:.6g} "
                       f"({counters['n_front_off_band']} steps outside it)")
    return dict(job, t=ts, front=fs, outcome=outcome, reason=reason, status=r.status,
                r_min=float(np.nanmin(fs)), t_min=float(ts[np.nanargmin(fs)]), steps=steps,
                seconds=tm.seconds, error=None, **counters)


def bisect_critical(args):
    """Bisection of the critical radius. Every step is recorded in `history`
    as a dict (r0, t_end, outcome, reason, status, seconds, ...).

    A step is decided only by 'ignition' or 'quenching'. An undecided run is
    repeated up to t = 2 T_EXP; if it is still undecided, if a run is
    'unresolved' or if it fails numerically, the bisection is aborted: the
    returned dict has r_crit = None and the key 'aborted' with the reason
    (the step is not counted as quenching)."""
    eps, dt, method, lo, hi, tol = args
    hist = []

    def result(aborted=None):
        out = dict(eps=eps, dt=dt, method=method, lo=lo, hi=hi,
                   r_crit=None if aborted else 0.5 * (lo + hi), history=hist)
        if aborted:
            out["aborted"] = aborted
        return out

    def one(mid, t_end):
        r = run_kernel(dict(r0=mid, eps=eps, dt=dt, method=method, t_end=t_end))
        hist.append({k: r[k] for k in ("r0", "t_end", "outcome", "reason", "status", "r_min",
                                       "t_min", "steps", "seconds", "n_secant_fail",
                                       "n_mesh_remesh", "n_front_off_band",
                                       "t_front_off_band") if k in r})
        return r["outcome"]

    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        try:
            outcome = one(mid, T_EXP)
            if outcome == "undecided":   # slow: extend the run
                outcome = one(mid, 2 * T_EXP)
        except (KernelRunFailed, NumericalFailure) as ex:
            hist.append(dict(r0=mid, outcome="failed", reason=f"{type(ex).__name__}: {ex}"))
            return result(aborted=f"run failed at r0 = {mid}: {type(ex).__name__}: {ex}")
        if outcome == "ignition":
            hi = mid
        elif outcome == "quenching":
            lo = mid
        else:
            return result(aborted=f"r0 = {mid} {outcome}: {hist[-1]['reason']}")
    return result()


EPS_C3 = (0.01, 0.005, 0.0025, 0.00125)
DT_C3 = 0.001


def run_all(what=("C1", "C2", "C3")):
    out = {}
    if "C1" in what:
        geoms = ("radial2", "radial3")
        t_ends = {g: t_end_of(g) for g in geoms}
        print("C1 t_end:", t_ends, flush=True)
        jobs = [(g, e, dt, t_ends[g]) for e in EPS_C1 for g in geoms for dt in DT_C1]
        res = []
        with worker_pool(4) as ex:
            for r in ex.map(run_C1, jobs):
                print("C1 {geom} eps={eps:<6} dt={dt:<7} front={front_err:.3e} "
                      "shifted={front_err_shifted:.3e} sup={sup_err:.3e} steps={steps} "
                      "{seconds:.1f}s".format(**r), flush=True)
                res.append(r)
        save_json("C1_limit_comparison.json", res)
    if "C2" in what:
        with worker_pool(2) as ex:
            # dt = 2.5e-5 for eps = 0.002: near r = 0 the mesh velocity is tapered
            # and an O(dt v/eps) lag appears (see README)
            eps_dt = ((0.01, 1e-4), (0.002, 2.5e-5))
            res = list(ex.map(run_C2, [(2, eps_dt), (3, eps_dt)]))
        for r in res:
            print("C2 N={N} T={T:.6f} R={R:.6f}".format(**r), "BDF T_eps", r["bdf_run"]["T_eps"],
                  [(e["eps"], round(e["T_eps"], 5), e["steps"], round(e["seconds"], 1), e["status"],
                    e["n_secant_fail"], e["n_mesh_remesh"]) for e in r["eps_runs"]], flush=True)
        save_json("C2_self_similar.json", res)
    if "C3" in what:
        with worker_pool(4) as ex:
            crit = list(ex.map(bisect_critical,
                               [(e, DT_C3, "ale", 0.6, 0.8, 0.001) for e in EPS_C3]
                               + [(0.00125, 0.00025, "ale", 0.6, 0.8, 0.001),
                                  (0.01, None, "euler", 0.6, 0.8, 0.001)]))
            for c in crit:
                print("C3 crit eps={eps} dt={dt} {method}: r0* in [{lo:.4f}, {hi:.4f}]".format(**c)
                      + ("  ABORTED: " + c["aborted"] if "aborted" in c else ""), flush=True)
                for h in c["history"]:
                    print("     r0={r0:.6f} {outcome}: {reason}".format(**h), flush=True)
            if any("aborted" in c for c in crit):
                # keep the audit trail, do not overwrite the stored results
                save_json("C3_quenching_aborted.json", dict(critical=crit, c0=C0, inv_c0=1 / C0))
                raise RuntimeError(
                    "C3: bisection aborted for " + ", ".join(
                        f"(eps={c['eps']}, dt={c['dt']}, {c['method']})" for c in crit if "aborted" in c)
                    + "; steps written to results/paper/C3_quenching_aborted.json")
            rc = {c["eps"]: c["r_crit"] for c in crit if c["method"] == "ale" and c["dt"] == DT_C3}
            r0s = [0.6, 0.8, rc[0.01] - 0.01, rc[0.01] + 0.01, 1.0, 1.2, 1.5]
            runs = list(ex.map(run_kernel, [dict(r0=r0, eps=0.01, dt=DT_C3, t_end=4.0)
                                            for r0 in r0s]))
            runs += list(ex.map(run_kernel, [dict(r0=r0, eps=0.0025, dt=DT_C3, t_end=4.0)
                                             for r0 in (rc[0.0025] - 0.01, rc[0.0025] + 0.01)]))
        for r in runs:
            print("C3 r0={r0:.4f} eps={eps} {outcome} r_min={r_min:.4f} at t={t_min:.3f} "
                  "{seconds:.1f}s [{reason}]".format(**r), flush=True)
        save_json("C3_quenching.json", dict(critical=crit, runs=runs, c0=C0, inv_c0=1 / C0))


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------

def fig_C1():
    d = load_json("C1_limit_comparison.json")
    old = json.loads((ROOT / "results" / "limit_comparison.json").read_text())
    fig, axs = plt.subplots(1, 3, figsize=(11.2, 4.1))
    ax = axs[0]
    sel = sorted((r for r in d if r["geom"] == "radial2" and r["dt"] == min(DT_C1)
                  and r["eps"] in (0.05, 0.01, 0.002)), key=lambda r: -r["eps"])
    for k, r in enumerate(sel):
        ax.plot(r["t"], r["front"], color=BLUES[k + 1], label=f"$\\varepsilon={r['eps']:g}$")
    ax.plot(sel[0]["lim_t"], sel[0]["lim_s"], ls="--", lw=1.2, color=INK, label="limit $(P)$")
    ax.set_xlabel("$t$")
    ax.set_ylabel("front radius")
    ax.legend(loc="lower left", fontsize=8)
    title(ax, f"$N=2$: fronts, ALE $\\Delta t={min(DT_C1):g}$")
    for ax, key, lab in ((axs[1], "front_err_shifted",
                          "$\\max_t|r_h+\\varepsilon/\\sqrt{2M}-s|$"),
                         (axs[2], "sup_err", "$\\max_r|u_h^+-u|$ at the final time")):
        for k, g in enumerate(("radial2", "radial3")):
            for dt, ls, mfc in zip(DT_C1, (":", "--", "-"), (SERIES[k], SERIES[k], "white")):
                rs = sorted((r for r in d if r["geom"] == g and r["dt"] == dt),
                            key=lambda r: r["eps"])
                ax.loglog([r["eps"] for r in rs], [r[key] for r in rs], color=SERIES[k], ls=ls,
                          marker=MARKERS[k], ms=5, mfc=mfc, mew=1.3,
                          label=f"ALE $N={g[-1]}$, $\\Delta t={dt:g}$")
            rs = sorted((r for r in old if r["geom"] == g), key=lambda r: r["eps"])
            ax.loglog([r["eps"] for r in rs], [r[key] for r in rs], color=MUTED, ls="none",
                      marker=MARKERS[k], ms=7, mfc="none", mew=0.8,
                      label=f"Eulerian BDF $N={g[-1]}$")
        slope_ref(ax, [1e-3, 5e-2], 4e-4, label="$\\varepsilon$")
        ax.set_xlabel("$\\varepsilon$")
        title(ax, lab)
    h, l = axs[1].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, fontsize=7.5, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.14, 1, 1))
    savefig(fig, "C1_ale_vs_limit.png")


def fig_C2():
    d = load_json("C2_self_similar.json")
    fig, axs = plt.subplots(2, 2, figsize=(8.4, 6.6))
    for j, r in enumerate(d):
        ax = axs[0, j]
        t, s, T = np.array(r["lim_t"]), np.array(r["lim_s"]), r["T"]
        m = T - t > 0
        ax.semilogx(T - t[m], s[m] / np.sqrt(T - t[m]), color=INK, lw=1.6, label="limit $(P)$")
        for k, e in enumerate(r["eps_runs"]):
            te, se = np.array(e["t"]), np.array(e["s"])
            m = (e["T_eps"] - te > 0) & (te > 0.01)
            ax.semilogx(e["T_eps"] - te[m], se[m] / np.sqrt(e["T_eps"] - te[m]),
                        color=SERIES[k], lw=1.2, ls="--",
                        label=f"$(P_\\varepsilon)$, $\\varepsilon={e['eps']:g}$")
        if "bdf_run" in r:
            b = r["bdf_run"]
            te, se = np.array(b["t"]), np.array(b["s"])
            m = (b["T_eps"] - te > 0) & (te > 0.01)
            ax.semilogx(b["T_eps"] - te[m], se[m] / np.sqrt(b["T_eps"] - te[m]), color=SERIES[2],
                        lw=1.2, ls=":", label=f"$(P_\\varepsilon)$, $\\varepsilon={b['eps']:g}$, BDF")
        ax.axhline(r["R"], color=MUTED, lw=0.8, ls=":")
        ax.annotate(f"$R={r['R']:.4f}$", (1.0, r["R"]), xycoords=("axes fraction", "data"),
                    xytext=(-2, 4), textcoords="offset points", ha="right", color=MUTED,
                    fontsize=8)
        ax.invert_xaxis()
        ax.set_xlabel("$T-t$")
        ax.set_ylim(r["R"] - 0.2, r["R"] + 0.06)
        title(ax, f"$N={r['N']}$: $s(t)/\\sqrt{{T-t}}$")
        ax.legend(loc="lower right", fontsize=8)
        ax = axs[1, j]
        ax.plot(r["eta"], r["phi"], color=INK, lw=2.4, alpha=0.35, label="$\\kappa\\phi(\\eta)$")
        for k, p in enumerate(r["profiles"]):
            ax.plot(p["eta"], p["v"], color=BLUES[k], lw=1.2, label=f"$s={p['s']:.2f}$")
        e = r["eps_runs"][-1]
        for k, p in enumerate(e["profiles"]):
            ax.plot(p["eta"], p["v"], color=BLUES[k], lw=0.8, ls="--")
        ax.set_xlim(0, r["R"] * 1.05)
        ax.set_ylim(bottom=0)
        ax.set_xlabel("$\\eta=r/\\sqrt{T-t}$")
        title(ax, f"$N={r['N']}$: $u/\\sqrt{{T-t}}$ (dashed: $\\varepsilon={e['eps']:g}$)")
        ax.legend(loc="upper right", fontsize=7.5)
    axs[0, 0].set_ylabel("$s/\\sqrt{T-t}$")
    axs[1, 0].set_ylabel("$u/\\sqrt{T-t}$")
    fig.tight_layout()
    savefig(fig, "C2_self_similar_convergence.png")


def fig_C3():
    d = load_json("C3_quenching.json")
    fig, axs = plt.subplots(1, 2, figsize=(8.4, 3.6), gridspec_kw=dict(width_ratios=[1.4, 1]))
    ax = axs[0]
    runs = sorted((r for r in d["runs"] if r["eps"] == 0.01), key=lambda r: r["r0"])
    y_prev = -1.0
    for k, r in enumerate(runs):
        col = {"quenching": SERIES[0], "ignition": SERIES[1]}.get(r["outcome"], MUTED)
        t, f = np.array(r["t"]), np.array(r["front"])
        ax.plot(t, f, color=col, lw=1.4)
        y_lab = max(f[0], y_prev + 0.18)
        y_prev = y_lab
        ax.annotate(f"$r_0={r['r0']:.3f}$", (t[0], f[0]), xytext=(-0.08, y_lab),
                    textcoords="data", ha="right", va="center", fontsize=7, color=INK,
                    arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.5))
    ax.plot([], [], color=SERIES[1], label="ignition")
    ax.plot([], [], color=SERIES[0], label="quenching")
    ax.axhline(d["inv_c0"], color=MUTED, ls=":", lw=0.8)
    ax.annotate("$1/c_0$", (4.0, d["inv_c0"]), xytext=(0, 3), textcoords="offset points",
                ha="right", color=MUTED, fontsize=8)
    ax.set_xlim(-0.75, 4.0)
    ax.set_ylim(0, 3.5)
    ax.set_xlabel("$t$")
    ax.set_ylabel("front radius $r_f(t)$")
    ax.legend(loc="lower right", fontsize=8)
    title(ax, f"$\\varepsilon=0.01$, ALE $\\Delta t={DT_C3:g}$")
    ax = axs[1]
    crit = sorted((c for c in d["critical"] if c["method"] == "ale" and c["dt"] == DT_C3),
                  key=lambda c: c["eps"])
    e = np.array([c["eps"] for c in crit])
    rc = np.array([c["r_crit"] for c in crit])
    err = np.array([0.5 * (c["hi"] - c["lo"]) for c in crit])
    ax.errorbar(e, rc, yerr=err, color=SERIES[0], marker="o", mfc="white", mew=1.3, ms=6,
                capsize=3, label=f"ALE, $\\Delta t={DT_C3:g}$")
    for c in d["critical"]:
        if c["method"] == "ale" and c["dt"] != DT_C3:
            ax.errorbar([c["eps"]], [c["r_crit"]], yerr=[0.5 * (c["hi"] - c["lo"])],
                        color=SERIES[2], marker="^", ms=6, mfc="none", capsize=3, ls="none",
                        label=f"ALE, $\\Delta t={c['dt']:g}$")
        if c["method"] == "euler":
            ax.errorbar([c["eps"]], [c["r_crit"]], yerr=[0.5 * (c["hi"] - c["lo"])],
                        color=SERIES[1], marker="s", ms=6, mfc="none", capsize=3, ls="none",
                        label="Eulerian BDF")
    # linear extrapolation to eps = 0
    p = np.polyfit(e, rc, 1)
    ee = np.linspace(0, e.max(), 50)
    ax.plot(ee, np.polyval(p, ee), color=MUTED, ls="--", lw=1,
            label=f"linear fit, $r_*(0)\\approx{p[1]:.3f}$")
    ax.axhline(d["inv_c0"], color=INK, ls=":", lw=1, label=f"$1/c_0={d['inv_c0']:.3f}$")
    ax.set_xlim(0, 0.0105)
    ax.set_ylim(0.68, 0.97)
    ax.set_xlabel("$\\varepsilon$")
    ax.set_ylabel("critical radius $r_*$")
    ax.legend(loc="best", fontsize=7.5)
    title(ax, "critical kernel radius")
    fig.tight_layout()
    savefig(fig, "C3_quenching.png")


def plot_all():
    if (RES / "C1_limit_comparison.json").exists():
        fig_C1()
    if (RES / "C2_self_similar.json").exists():
        fig_C2()
    if (RES / "C3_quenching.json").exists():
        fig_C3()


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("run", "all"):
        run_all(tuple(sys.argv[2:]) or ("C1", "C2", "C3"))
    if what in ("plot", "all"):
        plot_all()
