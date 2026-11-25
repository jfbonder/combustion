"""Group A: the meshes (ALE vs Eulerian), how they move, refine and coarsen.

A1  space-time diagram of the nodes, planar travelling wave (ALE | Eulerian)
A2  local mesh size h(x) at several times (travelling wave, radial extinction)
A3  number of nodes vs time and vs eps (band / transition / coarse), cost
A4  node trajectories r(t) for the radial extinction and an expanding flame

Run from numerics/:
    OMP_NUM_THREADS=1 ../.venv/bin/python -m experiments.paper.meshes [run|plot]
"""

import sys

import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm

from combustion1d.ale import simulate_ale
from combustion1d.exact import tw_speed_eps
from combustion1d.model import fb_slope
from combustion1d.solver import MeshParams, build_mesh, front_speed, initial_data, simulate

from .common import (BLUES, INK, MARKERS, MUTED, SERIES, Timer, expanding_u0,
                     extinction_setup, load_json, load_mesh_history, mesh_parts,
                     mesh_parts_predicted, plt, save_json,
                     save_mesh_history, savefig, title, tw_setup, worker_pool)

W0_TW = 0.0
EPS_MESH = 0.01
T_TW = 4.0
DT_TW = 0.01
EPS_NODES = (0.01, 0.005, 0.002, 0.001)
W0_RAD = 0.5
MP = MeshParams()
MP_ILL = MeshParams(h_coarse=0.04)    # coarser far field: node lines visible in pictures
MPS = dict(def_=MP, ill=MP_ILL)


def initial_mesh(u0, w0, eps, L, mp=MP):
    x = np.linspace(0.0, L, 20001)
    u, _ = initial_data(x, u0, w0, eps)
    return build_mesh(x, u, eps, L, mp)[0]


# ---------------------------------------------------------------------------
# runs
# ---------------------------------------------------------------------------

def run_tw(args):
    """ALE and Eulerian (BDF) travelling wave, w0 = 0, A = 1."""
    method, eps, mpn, record = args
    mp = MPS[mpn]
    u0, L, c0 = tw_setup(W0_TW, T_TW)
    c_eps = tw_speed_eps(1.0, W0_TW, eps)
    with Timer() as tm:
        if method == "ale":
            r = simulate_ale(u0, W0_TW, eps, L, T_TW, DT_TW, mesh=mp, speed_guess=c0,
                             record_mesh=record)
            hist = r.mesh_history
            node_t = r.t[1:]
        else:
            r = simulate(u0, W0_TW, eps, L, T_TW, N=1, mesh=mp, speed_guess=c0,
                         record_mesh=record)
            x_init = initial_mesh(u0, W0_TW, eps, L, mp)
            hist = [(0.0, x_init)] + r.mesh_history if record else []
            node_t = r.t[1:]
    x_init = initial_mesh(u0, W0_TW, eps, L, mp)
    out = dict(method=method, eps=eps, mesh=mpn, h_coarse=mp.h_coarse, L=L, c0=c0, c_eps=c_eps,
               c_num=front_speed(r, 0.5 * T_TW), steps=r.n_steps, n_remesh=r.n_remesh,
               seconds=tm.seconds, t=r.t, front=r.front,
               node_t=node_t, n_nodes=np.array(r.n_nodes),
               parts=mesh_parts(x_init, eps, mp),
               parts_pred=mesh_parts_predicted(eps, mp, L),
               trials=r.n_trials, solves=r.n_solves, lu=r.n_lu, fev=r.n_fev)
    if record:
        save_mesh_history(f"A_mesh_tw_{method}_{mpn}.npz", hist)
    return out


def run_radial_mesh(mpn):
    V, u0, L = extinction_setup(2, W0_RAD)
    with Timer() as tm:
        r = simulate_ale(u0, W0_RAD, EPS_MESH, L, 0.9, 0.005, mesh=MPS[mpn],
                         speed_guess=-V.R / 2, N=2, ell0=1.0, remesh_distance=0.1,
                         record_mesh=True)
    save_mesh_history(f"A_mesh_radial_ale_{mpn}.npz", r.mesh_history)
    return dict(case="radial", mesh=mpn, eps=EPS_MESH, L=L, R=V.R, t=r.t, front=r.front,
                steps=r.n_steps, n_remesh=r.n_remesh, seconds=tm.seconds,
                n_nodes=np.array(r.n_nodes))


def run_expanding_mesh(mpn):
    r0, L, t_end = 1.5, 12.0, 6.0
    c0 = float(fb_slope(W0_RAD))
    with Timer() as tm:
        r = simulate_ale(expanding_u0(r0, W0_RAD), W0_RAD, EPS_MESH, L, t_end, 0.01,
                         mesh=MPS[mpn], speed_guess=c0, N=2, record_mesh=True)
    save_mesh_history(f"A_mesh_expanding_ale_{mpn}.npz", r.mesh_history)
    return dict(case="expanding", mesh=mpn, eps=EPS_MESH, L=L, r0=r0, c0=c0, t=r.t, front=r.front,
                steps=r.n_steps, n_remesh=r.n_remesh, seconds=tm.seconds,
                n_nodes=np.array(r.n_nodes))


def run_all():
    jobs = [("ale", EPS_MESH, "ill", True), ("euler", EPS_MESH, "ill", True),
            ("ale", EPS_MESH, "def_", True), ("euler", EPS_MESH, "def_", False)]
    jobs += [(m, e, "def_", False) for e in EPS_NODES if e != EPS_MESH for m in ("ale", "euler")]
    with worker_pool(4) as ex:
        fr = ex.submit(run_radial_mesh, "ill")
        fr2 = ex.submit(run_radial_mesh, "def_")
        fe = ex.submit(run_expanding_mesh, "ill")
        fe2 = ex.submit(run_expanding_mesh, "def_")
        tw = []
        for r in ex.map(run_tw, jobs):
            print("{method:<5} eps={eps:<6} c_num-c_eps={d:+.2e} steps={steps} remesh={n_remesh} "
                  "solves={solves} nodes={n:.0f} {seconds:.1f}s".format(d=r["c_num"] - r["c_eps"],
                                                         n=np.mean(r["n_nodes"]), **r), flush=True)
            tw.append(r)
        rad, exp_ = fr.result(), fe.result()
        rad_d, exp_d = fr2.result(), fe2.result()
    print("radial: steps", rad["steps"], "remesh", rad["n_remesh"], f"{rad['seconds']:.1f}s")
    print("expanding: steps", exp_["steps"], "remesh", exp_["n_remesh"], f"{exp_['seconds']:.1f}s")
    save_json("A_meshes.json", dict(tw=tw, radial=rad, expanding=exp_, radial_def=rad_d,
                                    expanding_def=exp_d))


# ---------------------------------------------------------------------------
# plotting helpers
# ---------------------------------------------------------------------------

def node_segments(hist, xlim, tlim, euler=False):
    """Segments of the node trajectories inside the window, coloured by the
    local mesh size. At a remeshing step the old nodes are drawn up to the
    new record (the jump to the new nodes is not drawn)."""
    segs, hs = [], []
    for (t0, x0), (t1, x1) in zip(hist[:-1], hist[1:]):
        if t1 < tlim[0] or t0 > tlim[1]:
            continue
        h0 = np.gradient(x0)
        xa = x0
        if euler or x0.size != x1.size or np.max(np.abs(x1 - x0)) > 0.1:
            xb = x0
        else:
            xb = x1
        m = ((xa > xlim[0]) | (xb > xlim[0])) & ((xa < xlim[1]) | (xb < xlim[1]))
        segs.append(np.stack([np.stack([xa[m], np.full(m.sum(), t0)], -1),
                              np.stack([xb[m], np.full(m.sum(), t1)], -1)], 1))
        hs.append(h0[m])
    return np.concatenate(segs), np.concatenate(hs)


import matplotlib as _mpl
# viridis without its light-yellow end, so that thin coarse-node lines stay visible
NODE_CMAP = _mpl.colors.LinearSegmentedColormap.from_list(
    "viridis_dark", _mpl.colormaps["viridis"](np.linspace(0.0, 0.78, 256)))


def draw_nodes(ax, hist, xlim, tlim, norm, euler=False, lw=0.35, cmap=NODE_CMAP):
    segs, hs = node_segments(hist, xlim, tlim, euler)
    lc = LineCollection(segs, cmap=cmap, norm=norm, linewidths=lw, rasterized=True)
    lc.set_array(hs)
    ax.add_collection(lc)
    ax.set_xlim(*xlim)
    ax.set_ylim(*tlim)
    ax.grid(False)
    return lc


def ale_remesh_times(hist, t_front, front):
    """Records where the ALE mesh was rebuilt: node count changed, or some
    node moved much more than the front did in that step (in a normal ALE
    step every node moves by dt*v*chi, |chi| <= 1)."""
    out = []
    for (t0, x0), (t1, x1) in zip(hist[:-1], hist[1:]):
        if x0.size != x1.size:
            out.append(t1)
            continue
        df = abs(np.interp(t1, t_front, front) - np.interp(t0, t_front, front))
        if np.max(np.abs(x1 - x0)) > 2.0 * df + 1e-3 * np.min(np.diff(x0)):
            out.append(t1)
    return np.array(out)


def mark_remesh(ax, times, right=True):
    """Small ticks on the right border at the remeshing times."""
    for s in times:
        ax.annotate("", xy=(1.0, s), xycoords=("axes fraction", "data"),
                    xytext=(1.035, s), textcoords=("axes fraction", "data"),
                    arrowprops=dict(arrowstyle="-", color=SERIES[1], lw=1.0))


def _tw(d, method, mesh, eps=EPS_MESH):
    return next(r for r in d["tw"] if r["eps"] == eps and r["method"] == method
                and r["mesh"] == mesh)


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------

def fig_A1(d):
    eps = EPS_MESH
    mp = MP_ILL
    tw = {m: _tw(d, m, "ill") for m in ("ale", "euler")}
    H = {m: load_mesh_history(f"A_mesh_tw_{m}_ill.npz")[0] for m in ("ale", "euler")}
    norm = LogNorm(eps / mp.pts_per_eps, mp.h_coarse)
    t_ale, f_ale = np.array(tw["ale"]["t"]), np.array(tw["ale"]["front"])

    fig, axs = plt.subplots(2, 2, figsize=(7.6, 7.4), sharey="row")
    xlim1 = (f_ale[0] - 1.0, f_ale[-1] + 1.0)
    tlim1 = (0.0, T_TW)
    tm = 2.0
    xf_m = np.interp(tm, t_ale, f_ale)
    edge = xf_m + (mp.K_fresh + mp.K_margin) * eps
    xlim2 = (edge - 0.15, edge + 0.45)
    tlim2 = (tm - 0.3, tm + 0.3)
    for j, m in enumerate(("ale", "euler")):
        hist = H[m]
        r = tw[m]
        t, f = np.array(r["t"]), np.array(r["front"])
        rt = ale_remesh_times(hist, t, f) if m == "ale" else np.array([h[0] for h in hist[1:]])
        for i, (xl, tl) in enumerate(((xlim1, tlim1), (xlim2, tlim2))):
            ax = axs[i, j]
            lc = draw_nodes(ax, hist, xl, tl, norm, euler=(m == "euler"),
                            lw=0.3 if i == 0 else 0.45)
            ax.plot(f, t, color=SERIES[1], lw=1.3, label="front $u=\\varepsilon$")
            mark_remesh(ax, rt[(rt >= tl[0]) & (rt <= tl[1])])
            ax.set_xlabel("$x$")
        axs[0, j].text(0.03, 0.97, f"{r['n_remesh']} remeshings", transform=axs[0, j].transAxes,
                       va="top", ha="left", fontsize=8, color=INK,
                       bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))
        axs[0, j].add_patch(plt.Rectangle((xlim2[0], tlim2[0]), xlim2[1] - xlim2[0],
                                          tlim2[1] - tlim2[0], fill=False, ec=SERIES[0], lw=1.2))
    title(axs[0, 0], "ALE: nodes move with the flame")
    title(axs[0, 1], "Eulerian: fixed nodes, rebuilt meshes")
    title(axs[1, 0], "ALE, zoom (blue box)")
    title(axs[1, 1], "Eulerian, zoom (blue box)")
    axs[0, 0].set_ylabel("$t$")
    axs[1, 0].set_ylabel("$t$")
    axs[0, 0].legend(loc="lower right", fontsize=8, framealpha=0.9, frameon=True,
                     edgecolor="none")
    fig.tight_layout(rect=(0, 0, 0.87, 1), w_pad=2.5)
    cax = fig.add_axes([0.895, 0.25, 0.018, 0.5])
    cb = fig.colorbar(lc, cax=cax)
    cb.set_label("local mesh size $h$")
    savefig(fig, "A1_spacetime_nodes_tw.png")


def _h_profile(x):
    return 0.5 * (x[1:] + x[:-1]), np.diff(x)


def fig_A2(d):
    eps = EPS_MESH
    hist_tw, _ = load_mesh_history("A_mesh_tw_ale_def_.npz")
    hist_rad, _ = load_mesh_history("A_mesh_radial_ale_def_.npz")
    tw = _tw(d, "ale", "def_")
    rad = d["radial_def"]
    fig, axs = plt.subplots(1, 2, figsize=(8.4, 3.6))
    for ax, hist, times, res, lab in (
            (axs[0], hist_tw, (0.0, 1.3, 2.6, 4.0), tw, "planar travelling wave"),
            (axs[1], hist_rad, (0.0, 0.3, 0.6, 0.9), rad, "radial extinction, $N=2$")):
        tt = np.array([h[0] for h in hist])
        t_res, f_res = np.array(res["t"]), np.array(res["front"])
        for k, s in enumerate(times):
            x = hist[int(np.argmin(np.abs(tt - s)))][1]
            xm, h = _h_profile(x)
            ax.semilogy(xm, h, color=BLUES[k], lw=1.4, label=f"$t={s:g}$")
            xf = np.interp(s, t_res, f_res)
            ax.plot([xf], [eps / MP.pts_per_eps * 0.65], marker="^", color=BLUES[k], ms=5,
                    ls="none", clip_on=False)
        ax.axhline(eps / MP.pts_per_eps, color=MUTED, lw=0.6, ls=":")
        ax.axhline(MP.h_coarse, color=MUTED, lw=0.6, ls=":")
        ax.set_ylim(eps / MP.pts_per_eps * 0.5, 0.2)
        title(ax, lab)
        ax.set_ylabel("$h$")
    axs[0].set_xlabel("$x$")
    axs[1].set_xlabel("$r$")
    axs[0].set_xlim(0.0, 7.5)
    axs[1].set_xlim(0.0, rad["L"])
    axs[0].annotate("$h_f=\\varepsilon/16$", (7.4, eps / MP.pts_per_eps), xytext=(0, 3),
                    textcoords="offset points", ha="right", va="bottom", color=MUTED, fontsize=8)
    axs[0].annotate("$h_c$", (7.4, MP.h_coarse), xytext=(0, 3), textcoords="offset points",
                    ha="right", va="bottom", color=MUTED, fontsize=8)
    axs[0].legend(loc="upper right", fontsize=8, ncol=2)
    axs[1].legend(loc="upper right", fontsize=8, ncol=2)
    fig.tight_layout()
    savefig(fig, "A2_mesh_size.png")


def fig_A3(d):
    tw = [r for r in d["tw"] if r["mesh"] == "def_"]
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.4))
    ax = axs[0]
    for k, eps in enumerate((0.01, 0.001)):
        for m, ls in (("ale", "-"), ("euler", "--")):
            r = next(r for r in tw if r["eps"] == eps and r["method"] == m)
            tt, nn = np.array(r["node_t"]), np.array(r["n_nodes"])
            lab = ("ALE" if m == "ale" else "Eulerian") + f", $\\varepsilon={eps:g}$"
            if m == "euler":
                ax.step(tt, nn, where="post", color=SERIES[k], ls=ls, lw=1.2, label=lab)
            else:
                ax.plot(tt, nn, color=SERIES[k], ls=ls, lw=2.2, alpha=0.6, label=lab)
    ax.set_ylim(0, 4500)
    ax.set_xlabel("$t$")
    ax.set_ylabel("number of nodes")
    ax.legend(loc="lower right", fontsize=7.5)
    title(ax, "nodes vs time (travelling wave)")
    ax = axs[1]
    rs = sorted((r for r in tw if r["method"] == "ale"), key=lambda r: r["eps"])
    e = np.array([r["eps"] for r in rs])
    for k, key in enumerate(("band", "transition", "coarse")):
        ax.semilogx(e, [r["parts"][key] for r in rs], color=SERIES[k], marker=MARKERS[k],
                    mfc="white", mew=1.5, ms=6, label=key)
        ax.semilogx(e, [r["parts_pred"][key] for r in rs], color=SERIES[k], ls=":", lw=1)
    ax.set_xlabel("$\\varepsilon$")
    ax.set_ylabel("nodes")
    ax.set_yscale("log")
    ax.legend(loc="center left", fontsize=8)
    title(ax, "nodes per zone (dotted: formula)")
    ax = axs[2]
    for k, m in enumerate(("ale", "euler")):
        rs = sorted((r for r in tw if r["method"] == m), key=lambda r: r["eps"])
        e = np.array([r["eps"] for r in rs])
        ax.loglog(e, [r["steps"] for r in rs], color=SERIES[k], marker=MARKERS[k], mfc="white",
                  mew=1.5, ms=6, label="ALE ($\\Delta t=0.01$)" if m == "ale" else
                  "Eulerian (adaptive BDF)")
    ax.set_xlabel("$\\varepsilon$")
    ax.set_ylabel("time steps to $t=4$")
    ax.legend(loc="upper right", fontsize=8)
    title(ax, "time steps vs $\\varepsilon$")
    fig.tight_layout()
    savefig(fig, "A3_node_count.png")


def fig_A4(d):
    eps = EPS_MESH
    mp = MP_ILL
    norm = LogNorm(eps / mp.pts_per_eps, mp.h_coarse)
    hist_rad, _ = load_mesh_history("A_mesh_radial_ale_ill.npz")
    hist_exp, _ = load_mesh_history("A_mesh_expanding_ale_ill.npz")
    rad, ex = d["radial"], d["expanding"]
    fig, axs = plt.subplots(1, 2, figsize=(8.6, 4.3))
    for ax, hist, res, xl, lab in (
            (axs[0], hist_rad, rad, (0.0, 3.4), "radial extinction ($N=2$)"),
            (axs[1], hist_exp, ex, (0.0, 8.5), "expanding flame ($N=2$, $r_0=1.5$)")):
        t, f = np.array(res["t"]), np.array(res["front"])
        lc = draw_nodes(ax, hist, xl, (0.0, t[-1]), norm, lw=0.3)
        ax.plot(f, t, color=SERIES[1], lw=1.3, label="front $u=\\varepsilon$")
        rt = ale_remesh_times(hist, t, f)
        mark_remesh(ax, rt)
        title(ax, lab)
        ax.set_xlabel("$r$")
        ax.text(0.97, 0.03, f"{res['n_remesh']} remeshings, {res['steps']} steps",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color=INK,
                bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))
    axs[0].set_ylabel("$t$")
    axs[0].legend(loc="upper right", fontsize=8, frameon=True, framealpha=0.9, edgecolor="none")
    fig.tight_layout(rect=(0, 0, 0.88, 1), w_pad=2.5)
    cax = fig.add_axes([0.905, 0.25, 0.018, 0.5])
    cb = fig.colorbar(lc, cax=cax)
    cb.set_label("local mesh size $h$")
    savefig(fig, "A4_radial_meshes.png")


def plot_all():
    d = load_json("A_meshes.json")
    fig_A1(d)
    fig_A2(d)
    fig_A3(d)
    fig_A4(d)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("run", "all"):
        run_all()
    if what in ("plot", "all"):
        plot_all()
