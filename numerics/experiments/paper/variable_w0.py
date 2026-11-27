"""Group I: a non-constant perturbation w0. Figure 3 of the article (label fig:I1,
figures/paper/I1_variable_w0.png) and the numbers of Section 6.3 (label
sec:exp:varw); numbers of the draft in paper/.

All the other experiments use w^eps = eps*w0 with w0 constant, so that M(w0)
is a constant. Here w^eps(x, 0) = eps*W(x, 0) with a non-constant W; then
W(x, t) solves the heat equation and the limit problem is

    u_t = Delta u in {u > 0},   |grad u| = sqrt(2 M(W(x, t))) on the free boundary.

We compare (P_eps), computed with the reference solver (BDF) and with ALE,
against the limit problem solved by front fixing (combustion1d.limit with a
callable w0), for eps in {0.05, 0.02, 0.01, 0.005, 0.002, 0.001}.

  planar   L = 24, W(x, t) = 0.5 + 0.4 exp(-k^2 t) cos(k x), k = pi/6 (zero flux
           at both ends), u0(x) = A tanh(lam0 (x - 2)/A)^+, lam0 = sqrt(2M(W(2, 0)));
           a front that propagates into a mixture stratified in enthalpy.
           t <= 6; ALE with dt = 2.5e-4.
  radial3  N = 3, L = 6, W(r, t) = 0.5 - 0.3 exp(-k^2 t) sin(k r)/(k r) with
           tan(kL) = kL (zero flux at r = L), u0 = a (1 - r^2/R0^2)^+, R0 = 2,
           a = sqrt(2M(W(R0, 0))) R0/2; a pocket of fresh gas, until the radius
           of the limit solution is R0/2; ALE with dt = 1.25e-4. Here
           lambda -> W(lambda x, lambda^2 t) is nondecreasing, which is the
           monotonicity condition of the rate theorem.

In both cases the data satisfy the free boundary condition at t = 0. We record
  * the error in the position of the front {u^eps = eps}, max over t, before
    and after removing the shift eps/sqrt(2M) (with the local value of M),
  * the sup error of (u^eps)^+ at the final time,
  * the error of the discrete W_h = (Y - U)/eps against the exact caloric W
    at the final time (a check of the scheme for W).

Output: results/paper/I_variable_w0.json = {"eps", "t_end", "rows", "limits"};
rows: geom, solver ("bdf" or "ale"), eps, dt, status, front_err,
front_err_shifted, sup_err, W_err, t, front, steps, seconds; limits[geom]: the
fronts of the limit problem with W (t, s) and with w0 = 0.5 (t_const,
s_const), and W at the free boundary (W_front).

Run from numerics/:  python -m experiments.paper.variable_w0 [run|plot]   (about 9 min, 4 processes)
"""

import sys
import time

import numpy as np

from combustion1d.ale import simulate_ale
from combustion1d.limit import solve_planar, solve_radial
from combustion1d.model import fb_slope
from combustion1d.solver import MeshParams, simulate

from .common import FIG, load_json, save_json, worker_pool

EPS = (0.05, 0.02, 0.01, 0.005, 0.002, 0.001)
A, X0, L_PLANAR, K_PLANAR = 1.0, 2.0, 24.0, np.pi / 6.0
R0, L_RADIAL = 2.0, 6.0
K_RADIAL = 4.493409457909064 / L_RADIAL          # first positive root of tan(x) = x
W_MEAN = 0.5
NAME = "I_variable_w0.json"


def W_planar(x, t):
    return W_MEAN + 0.4 * np.exp(-K_PLANAR**2 * t) * np.cos(K_PLANAR * x)


def W_radial(r, t):
    return W_MEAN - 0.3 * np.exp(-K_RADIAL**2 * t) * np.sinc(K_RADIAL * np.asarray(r) / np.pi)


def case(geom):
    if geom == "planar":
        lam0 = float(fb_slope(W_planar(X0, 0.0)))
        u0 = lambda x: np.where(x > X0, A * np.tanh(lam0 * (x - X0) / A), 0.0)
        return dict(N=1, u0=u0, W=W_planar, L=L_PLANAR, sign=+1.0, s0=X0, dt_ale=2.5e-4)
    lam0 = float(fb_slope(W_radial(R0, 0.0)))
    a = lam0 * R0 / 2.0
    u0 = lambda r: a * np.maximum(1.0 - (r / R0) ** 2, 0.0)
    return dict(N=3, u0=u0, W=W_radial, L=L_RADIAL, sign=-1.0, s0=R0, dt_ale=1.25e-4)


def limit_solution(geom, t_end, w=None, n=1600):
    c = case(geom)
    w = c["W"] if w is None else w
    t_eval = np.linspace(0.0, t_end, 601)
    if geom == "planar":
        return solve_planar(c["u0"], w, c["L"], c["s0"], t_end, n=n, t_eval=t_eval)
    return solve_radial(c["u0"], w, c["N"], c["s0"], t_end, n=n, t_eval=t_eval)


def t_end_of(geom):
    if geom == "planar":
        return 6.0
    c = case(geom)                         # radial: until the limit ball has radius R0/2
    r = solve_radial(c["u0"], c["W"], c["N"], c["s0"], 10.0, n=400,
                     t_eval=np.linspace(0, 10, 10001))
    return float(r.t[np.argmax(r.s < 0.5 * R0)])


def run(job):
    geom, solver, eps, t_end = job
    c = case(geom)
    lim = limit_solution(geom, t_end)
    w_init = lambda x: c["W"](x, 0.0)
    t0 = time.perf_counter()
    if solver == "bdf":
        r = simulate(c["u0"], w_init, eps, c["L"], t_end, N=c["N"], speed_guess=0.5,
                     snapshot_times=(t_end,), max_dt_remesh=0.01 * t_end,
                     mesh=MeshParams(h_coarse=0.01))
    else:
        r = simulate_ale(c["u0"], w_init, eps, c["L"], t_end, c["dt_ale"], N=c["N"],
                         mesh=MeshParams(h_coarse=0.01), speed_guess=0.5,
                         snapshot_times=(t_end,))
    seconds = time.perf_counter() - t0
    t, front = np.asarray(r.t), np.asarray(r.front, dtype=float)
    m = (t >= 0.05 * t_end) & np.isfinite(front)
    s = np.interp(t[m], lim.t, lim.s)
    lam = np.array([float(fb_slope(c["W"](si, ti))) for si, ti in zip(s, t[m])])
    d = front[m] - s
    tk = max(r.snapshots)
    x, u, w = r.snapshots[tk]
    k = int(np.argmin(np.abs(lim.t - tk)))
    sup = np.abs(np.maximum(u, 0.0) - lim.u_on(x, k))
    return dict(geom=geom, solver=solver, eps=eps, t_end=t_end, status=r.status,
                dt=c["dt_ale"] if solver == "ale" else None,
                front_err=float(np.max(np.abs(d))),
                front_err_shifted=float(np.max(np.abs(d - c["sign"] * eps / lam))),
                sup_err=float(np.max(sup)),
                W_err=float(np.max(np.abs(w / eps - c["W"](x, tk)))),
                t=t[m][::max(1, m.sum() // 400)].tolist(),
                front=front[m][::max(1, m.sum() // 400)].tolist(),
                steps=int(r.n_steps), seconds=seconds)


def run_all():
    geoms = ("planar", "radial3")
    t_ends = {g: t_end_of(g) for g in geoms}
    print("t_end:", t_ends, flush=True)
    jobs = [(g, sv, eps, t_ends[g]) for sv in ("bdf", "ale") for eps in EPS for g in geoms]
    rows = []
    with worker_pool(4) as ex:
        for r in ex.map(run, jobs):
            print("{geom:<8} {solver:<3} eps={eps:<6} front={front_err:.3e} "
                  "shifted={front_err_shifted:.3e} sup={sup_err:.3e} W={W_err:.1e} "
                  "steps={steps} {seconds:.0f}s [{status}]".format(**r), flush=True)
            rows.append(r)
    limits = {}
    for g in geoms:                        # the limit fronts: variable W and constant W = 0.5
        lv, lc = limit_solution(g, t_ends[g]), limit_solution(g, t_ends[g], w=W_MEAN)
        limits[g] = dict(t=lv.t.tolist(), s=lv.s.tolist(), t_const=lc.t.tolist(), s_const=lc.s.tolist(),
                         W_front=[float(case(g)["W"](si, ti)) for si, ti in zip(lv.s, lv.t)])
    save_json(NAME, dict(eps=list(EPS), t_end=t_ends, rows=rows, limits=limits))


def plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d = load_json(NAME)
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 3.9))
    col = {"planar": "C0", "radial3": "C3"}
    lab = {"planar": "planar", "radial3": "radial, $N=3$"}
    for g in ("planar", "radial3"):
        lim = d["limits"][g]
        ax[0].plot(lim["t"], lim["s"], color=col[g], lw=1.6, label=f"(P), {lab[g]}")
        ax[0].plot(lim["t_const"], lim["s_const"], color=col[g], lw=1.0, ls=":",
                   label="(P) with $w_0\\equiv0.5$")
        for r in d["rows"]:
            if r["geom"] == g and r["solver"] == "ale" and r["eps"] == 0.01:
                ax[0].plot(r["t"][::25], r["front"][::25], "o", ms=3.5, mfc="none", color="k",
                           label="$(P_\\varepsilon)$, $\\varepsilon=0.01$" if g == "planar" else None)
        ax[1].plot(lim["t"], lim["W_front"], color=col[g], label=lab[g])
        for sv, mk in (("bdf", "s-"), ("ale", "o--")):
            rows = sorted((r for r in d["rows"] if r["geom"] == g and r["solver"] == sv),
                          key=lambda r: r["eps"])
            e = np.array([r["eps"] for r in rows])
            ax[2].loglog(e, [r["front_err_shifted"] for r in rows], mk, color=col[g], ms=4,
                         mfc="none", label=f"front$^*$, {lab[g]}, {sv.upper()}")
            ax[2].loglog(e, [r["sup_err"] for r in rows], mk, color=col[g], ms=4, alpha=0.45,
                         label=f"sup, {lab[g]}, {sv.upper()}")
    e = np.array(d["eps"])
    ax[2].loglog(e, 0.5 * e, "k:", lw=1, label="$0.5\\,\\varepsilon$")
    ax[0].set_xlabel("$t$"); ax[0].set_ylabel("position of the front")
    ax[1].set_xlabel("$t$"); ax[1].set_ylabel("$w_0$ at the free boundary")
    ax[2].set_xlabel("$\\varepsilon$"); ax[2].set_ylabel("error")
    ax[0].legend(fontsize=7); ax[1].legend(fontsize=8); ax[2].legend(fontsize=6, ncol=2)
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG / "I1_variable_w0.png", dpi=160)
    print("wrote", FIG / "I1_variable_w0.png")


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("run", "all"):
        run_all()
    if what in ("plot", "all"):
        plot()
