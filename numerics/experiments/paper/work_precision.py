"""Group J: work-precision comparison for the planar travelling wave. Figure S5
(label fig:J1, figures/paper/J1_work_precision.png) and Table S1 (label
tab:work) of the supplement, and the paragraph on the cost in Section 6.6 of
the article; numbers of the drafts in paper/.

Problem: (P_eps) with w0 = 0, A = 1, eps in {0.01, 0.001}, started from the
exact travelling wave Phi_eps, on [0, L], up to T = 4. Error at the final time:
the error in the position of the front {u = eps} and the sup error against
Phi_eps(x - x0 - c_eps T).
Work: wall-clock seconds of the solver call (single thread, pure Python/SciPy
implementations, so only indicative) and implementation-independent counts
(time steps; tridiagonal solves for the scheme (S); LU factorizations and
right-hand side evaluations for the method-of-lines solvers).

Methods
  ale       scheme (S) on the moving band mesh, implicit mesh velocity
            (combustion1d.ale), h_f = eps/8, ..., eps/64, dt from 0.2 to 0.0025;
  eul       scheme (S) on band meshes that do not move (combustion1d.simple),
            h_f = eps/16, dt = eps/4, ..., eps/1024 (up to eps/64 for eps = 0.001);
  uni       scheme (S) on a uniform mesh of size eps/16 (eps = 0.01 only);
  bdf       method of lines on band meshes, SciPy BDF with the sparse analytic
            Jacobian and SuperLU (combustion1d.solver.simulate), rtol from
            1e-3 to 1e-6;
  lsoda     the same method of lines with the unknowns interleaved
            (u_0, w_0, u_1, w_1, ...), which makes the Jacobian banded
            (bandwidth 2), integrated with LSODA and its banded linear algebra.

Output: results/paper/J_work_precision.json = {"T", "w0", "workers", "rows"};
rows: kind, eps, dt or rtol, pts (ale), front_err, sup_err, seconds, steps,
nodes and the counters solves (ale, eul, uni), n_lu and n_fev (bdf), n_fev,
n_jev and n_remesh (lsoda). The column "solves" of Table S1 is solves, n_lu
and n_jev, respectively.

Run from numerics/:  python -m experiments.paper.work_precision [run|plot]
(about 26 min, 2 processes. The measured seconds depend on the machine and on
its load.)
"""

import sys
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

from combustion1d.ale import simulate_ale
from combustion1d.exact import tw_profile_eps, tw_speed_eps
from combustion1d.model import resolve_reaction
from combustion1d.simple import simulate_simple
from combustion1d.solver import MeshParams, build_mesh, find_front, laplacian, simulate

from .common import FIG, load_json, save_json, worker_pool

W0, T, X0, THETA = 0.0, 4.0, 2.0, 0.37
EPS = (0.01, 0.001)
NAME = "J_work_precision.json"

ALE_DTS = (0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.0025)
ALE_PTS = (8, 16, 32, 64)
EUL_RATIOS = {0.01: (4, 16, 64, 256, 1024), 0.001: (4, 16, 64)}      # dt = eps/ratio
UNI_RATIOS = {0.01: (4, 16, 64)}
RTOLS = (1e-3, 1e-4, 1e-5, 1e-6)


def setup(eps, pts=16):
    c = tw_speed_eps(1.0, W0, eps)
    h_f = eps / pts
    x0 = (round(X0 / h_f) + THETA) * h_f
    L = x0 + c * T + 12.0 / c + 2.0
    P = lambda z: tw_profile_eps(z, 1.0, W0, eps)
    return c, x0, L, P


def errors(x, u, eps, c, x0, P):
    xe = x0 + c * T
    xf = find_front(x, u, eps)[0]
    return dict(front_err=float(abs(xf - xe)), sup_err=float(np.max(np.abs(u - P(x - xe)))))


def simulate_lsoda(u0, eps, L, t_end, mesh, speed_guess, rtol, atol_rel=10.0, max_dt_remesh=None):
    """Method of lines on the band meshes of combustion1d.solver.simulate (same
    remeshing rule and PCHIP transfer), with interleaved unknowns and LSODA
    with a banded (finite difference) Jacobian. Planar case, w(x, 0) = eps*W0."""
    mp = mesh
    f_r = resolve_reaction(None)[0]
    if max_dt_remesh is None:
        max_dt_remesh = 0.05 * t_end
    x = np.linspace(0.0, L, 20001)
    x, _ = build_mesh(x, u0(x), eps, L, mp)
    u, w = u0(x), np.full_like(x, eps * W0)
    xf = find_front(x, u, eps)[0]
    t, v_est = 0.0, speed_guess
    out = dict(steps=0, n_fev=0, n_jev=0, n_remesh=0)
    while t < t_end - 1e-14:
        n = x.size
        Lap = laplacian(x, 1)

        def rhs(_, z):
            uu, ww = z[0::2], z[1::2]
            dz = np.empty_like(z)
            dz[0::2] = Lap @ uu - (uu + ww) * f_r(uu / eps) / eps**2
            dz[1::2] = Lap @ ww
            return dz

        dt = min(0.5 * mp.K_margin * eps / max(v_est, 1e-12), max_dt_remesh, t_end - t)
        z0 = np.empty(2 * n)
        z0[0::2], z0[1::2] = u, w
        sol = solve_ivp(rhs, (t, t + dt), z0, method="LSODA", lband=2, uband=2,
                        rtol=rtol, atol=atol_rel * rtol * eps)
        if not sol.success:
            raise RuntimeError(sol.message)
        t = sol.t[-1]
        out["steps"] += sol.t.size - 1
        out["n_fev"] += sol.nfev
        out["n_jev"] += sol.njev
        u, w = sol.y[0::2, -1], sol.y[1::2, -1]
        xf_new = find_front(x, u, eps)[0]
        v_est = max(abs(xf_new - xf) / dt, 1e-3)
        xf = xf_new
        if t < t_end - 1e-14:
            x_new, _ = build_mesh(x, u, eps, L, mp)
            with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
                u = PchipInterpolator(x, u)(x_new)
                w = PchipInterpolator(x, w)(x_new)
            x = x_new
            out["n_remesh"] += 1
    return x, u, out


def run(job):
    kind, eps = job["kind"], job["eps"]
    pts = job.get("pts", 16)
    c, x0, L, P = setup(eps, pts)
    u0 = lambda x: P(x - x0)
    mp = MeshParams(pts_per_eps=pts)
    if kind == "uni":
        mp = MeshParams(pts_per_eps=pts, h_coarse=eps / pts)
    t0 = time.perf_counter()
    if kind == "ale":
        r = simulate_ale(u0, W0, eps, L, T, job["dt"], mesh=mp, speed_guess=c,
                         exact_initial=True, snapshot_times=(T,))
        x, u, _ = r.snapshots[T]
        work = dict(steps=r.n_steps, solves=r.n_solves, nodes=float(np.mean(r.n_nodes)),
                    n_secant_fail=r.n_secant_fail)
    elif kind in ("eul", "uni"):
        r = simulate_simple(u0, W0, eps, L, T, dt=job["dt"], mesh=mp, speed_guess=c,
                            exact_initial=True, snapshot_times=(T,))
        x, u, _ = r.snapshots[T]
        work = dict(steps=r.n_steps, solves=r.n_solves, nodes=float(np.mean(r.n_nodes)))
    elif kind == "bdf":
        r = simulate(u0, W0, eps, L, T, mesh=mp, speed_guess=c, snapshot_times=(T,),
                     rtol=job["rtol"], atol_rel=10.0 * job["rtol"])
        x, u, _ = r.snapshots[T]
        work = dict(steps=r.n_steps, n_lu=r.n_lu, n_fev=r.n_fev, nodes=float(np.mean(r.n_nodes)))
    elif kind == "lsoda":
        x, u, work = simulate_lsoda(u0, eps, L, T, mp, c, job["rtol"])
        work["nodes"] = float(x.size)
    else:
        raise ValueError(kind)
    seconds = time.perf_counter() - t0
    out = dict(job, seconds=seconds, c_eps=c, L=L)
    out.update(work)
    out.update(errors(x, u, eps, c, x0, P))
    return out


def jobs():
    js = []
    for eps in EPS:
        for pts in ALE_PTS:
            js += [dict(kind="ale", eps=eps, pts=pts, dt=dt) for dt in ALE_DTS]
        js += [dict(kind="bdf", eps=eps, rtol=rt) for rt in RTOLS]
        js += [dict(kind="lsoda", eps=eps, rtol=rt) for rt in RTOLS]
        js += [dict(kind="eul", eps=eps, dt=eps / k) for k in EUL_RATIOS[eps]]
        js += [dict(kind="uni", eps=eps, dt=eps / k) for k in UNI_RATIOS.get(eps, ())]
    return js


def run_all(workers=2):
    rows = []
    with worker_pool(workers) as ex:
        for r in ex.map(run, jobs()):
            print("{kind:<5} eps={eps:<6} ".format(**r)
                  + (f"dt={r['dt']:<9.3g} pts={r.get('pts', 16):<3}" if "dt" in r else f"rtol={r['rtol']:<8.0e}")
                  + " front={front_err:.2e} sup={sup_err:.2e} steps={steps} {seconds:.1f}s".format(**r),
                  flush=True)
            rows.append(r)
    save_json(NAME, dict(T=T, w0=W0, workers=workers, rows=rows))


def plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d = load_json(NAME)
    style = {"ale": ("C0", "o-", "ALE, $h_f=\\varepsilon/16$"), "eul": ("C1", "s-", "Eulerian (S), band mesh"),
             "uni": ("C2", "^-", "Eulerian (S), uniform mesh"), "bdf": ("C3", "d-", "BDF (sparse LU)"),
             "lsoda": ("C4", "v-", "LSODA (banded)")}
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for a, eps in zip(ax, EPS):
        for kind, (col, mk, lab) in style.items():
            rows = [r for r in d["rows"] if r["kind"] == kind and r["eps"] == eps
                    and r.get("pts", 16) == 16]
            if not rows:
                continue
            rows.sort(key=lambda r: -r.get("dt", r.get("rtol")))   # coarse to fine
            a.loglog([r["seconds"] for r in rows], [r["front_err"] for r in rows], mk, color=col,
                     ms=4, mfc="none", label=lab)
        for pts, ls in ((8, ":"), (32, "--"), (64, "-.")):
            rows = sorted((r for r in d["rows"] if r["kind"] == "ale" and r["eps"] == eps
                           and r.get("pts") == pts), key=lambda r: -r["dt"])
            a.loglog([r["seconds"] for r in rows], [r["front_err"] for r in rows], "o" + ls,
                     color="C0", ms=3, alpha=0.6, label=f"ALE, $h_f=\\varepsilon/{pts}$")
        a.set_title(f"$\\varepsilon={eps:g}$", fontsize=10)
        a.set_xlabel("seconds")
        a.grid(alpha=0.3, which="both")
    ax[0].set_ylabel("error in the position of the front at $T=4$")
    ax[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "J1_work_precision.png", dpi=160)
    print("wrote", FIG / "J1_work_precision.png")


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("run", "all"):
        run_all()
    if what in ("plot", "all"):
        plot()
