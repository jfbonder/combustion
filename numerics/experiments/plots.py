"""Figures from results/*.json  ->  figures/*.png

Run from numerics/:  python -m experiments.plots
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RES, FIG = ROOT / "results", ROOT / "figures"

# categorical slots in fixed order (validated palette), markers as 2nd channel
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]
MARKERS = ["o", "s", "^"]
INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3de"
BLUES = ["#9ec5f4", "#5598e7", "#256abf", "#0d366b"]   # sequential: time

plt.rcParams.update({
    "figure.dpi": 150, "font.size": 9, "axes.edgecolor": MUTED,
    "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "lines.linewidth": 2,
})


def slope_ref(ax, x, y0, label="pendiente 1"):
    x = np.array([min(x), max(x)])
    ax.plot(x, y0 * x / x[0], ls="--", lw=1, color=MUTED)
    ax.annotate(label, (x[1], y0 * x[1] / x[0]), xytext=(4, 0),
                textcoords="offset points", color=MUTED, va="center")


def tw_figures():
    d = json.loads((RES / "traveling_wave.json").read_text())
    A = d["A"]
    w0s = sorted({r["w0"] for r in A})

    # --- speed convergence -------------------------------------------------
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    for k, w0 in enumerate(w0s):
        rs = sorted((r for r in A if r["w0"] == w0), key=lambda r: r["eps"])
        e = np.array([r["eps"] for r in rs])
        ce = np.array([abs(r["c_eps"] - r["c0"]) for r in rs])
        cp = np.array([abs(r["c_pde"] - r["c0"]) for r in rs])
        ax.loglog(e, ce, color=SERIES[k], label=f"$w_0={w0:+g}$")
        ax.loglog(e, cp, ls="none", marker=MARKERS[k], ms=7, mfc="white",
                  mec=SERIES[k], mew=1.5)
        ax.annotate(f"$w_0={w0:+g}$", (e[-1], ce[-1]), xytext=(6, 0),
                    textcoords="offset points", color=INK, va="center")
    slope_ref(ax, [1e-3, 1e-1], 4e-4)
    ax.set_xlabel(r"$\varepsilon$")
    ax.set_ylabel(r"$|c - c_0|$")
    ax.set_title("Velocidad de la llama: línea = onda viajera exacta de $(P_\\varepsilon)$,\n"
                 "marcadores = EDP con malla móvil", fontsize=9, color=INK, loc="left")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "tw_speed_convergence.png")
    plt.close(fig)

    # --- instability: same u0, different w0 ----------------------------------
    B = sorted(d["B"], key=lambda r: r["w0"])
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(8.4, 3.4),
                                  gridspec_kw=dict(width_ratios=[1.1, 1]))
    for k, r in enumerate(B):
        t, xf = np.array(r["t"]), np.array(r["front"])
        ax.plot(t, xf, color=SERIES[k], label=f"$w_0={r['w0']:+g}$")
        lim = _limit_front_B(r["w0"], t[-1])
        ax.plot(lim.t, lim.s, ls="--", lw=1, color=INK)
        ax.annotate(f"$w_0={r['w0']:+g}$", (t[-1], xf[-1]), xytext=(4, 0),
                    textcoords="offset points", color=INK, va="center")
        x, u, Y = map(np.array, (r["x"], r["u"], r["Y"]))
        m = np.abs(x - xf[-1]) < 0.12
        ax2.plot(x[m] - xf[-1], Y[m], color=SERIES[k], label=f"$Y$, $w_0={r['w0']:+g}$")
        ax2.plot(x[m] - xf[-1], u[m], color=SERIES[k], lw=1, ls="--")
    ax.set_xlabel("$t$")
    ax.set_ylabel("posición del frente")
    ax.set_title("Mismo $u_0$, $Y_0=u_0+\\varepsilon w_0$: $(P_\\varepsilon)$ con $\\varepsilon=0.005$\n"
                 "(color) y problema límite $(P)$ calculado (rayado)", fontsize=9,
                 color=INK, loc="left")
    ax.legend(loc="upper left")
    ax2.set_xlabel("$x - x_f$")
    ax2.set_title("Zona de reacción al final: $Y$ (llena), $u$ (rayada)",
                  fontsize=9, color=INK, loc="left")
    ax2.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "tw_instability.png")
    plt.close(fig)


def radial_figures():
    import sys
    sys.path.insert(0, str(ROOT))
    from combustion1d.exact import SelfSimilar

    rs = json.loads((RES / "radial.json").read_text())
    Ns = sorted({r["N"] for r in rs})

    fig, axs = plt.subplots(1, 2, figsize=(8.4, 3.4), sharey=True)
    for k, N in enumerate(Ns):
        g = sorted((r for r in rs if r["N"] == N), key=lambda r: r["eps"])
        e = np.array([r["eps"] for r in g])
        for ax, key in zip(axs, ("front_err", "sup_err")):
            ax.loglog(e, [r[key] for r in g], color=SERIES[k], marker=MARKERS[k],
                      ms=6, mfc="white", mew=1.5, label=f"$N={N}$")
    slope_ref(axs[0], [1e-3, 5e-2], 1.0e-3, label=r"$\varepsilon$")
    slope_ref(axs[1], [1e-3, 5e-2], 2.5e-4, label=r"$\varepsilon/4$")
    axs[0].set_title(r"$\max_t\,|r_\varepsilon(t) - R\sqrt{T-t}|$", fontsize=9,
                     color=INK, loc="left")
    axs[1].set_title(r"$\max_x\,|(u^\varepsilon)^+ - V|$ en $t=0.75$", fontsize=9,
                     color=INK, loc="left")

    for ax in axs:
        ax.set_xlabel(r"$\varepsilon$")
        ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "radial_convergence.png")
    plt.close(fig)

    r = next(r for r in rs if r["N"] == 2 and r["eps"] == 0.01)
    V = SelfSimilar(2, r["w0"], 1.0)
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    for k, (t, s) in enumerate(sorted(r["snapshots"].items(), key=lambda kv: float(kv[0]))):
        t = float(t)
        x, u = np.array(s["x"]), np.array(s["u"])
        ax.plot(x, u, color=BLUES[k], label=f"$t={t:g}$")
        xx = np.linspace(0, x[-1], 800)
        ax.plot(xx, V(xx, t), color=INK, lw=1, ls="--")
    ax.set_xlim(0, 1.2 * V.R)
    ax.set_xlabel("$r$")
    ax.set_ylabel("$u$")
    ax.set_title("$N=2$, $\\varepsilon=0.01$: $u^\\varepsilon$ (color) y solución exacta $V$ (rayada)",
                 fontsize=9, color=INK, loc="left")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(FIG / "radial_snapshots.png")
    plt.close(fig)


def _limit_front_B(w0, t_end):
    import sys
    sys.path.insert(0, str(ROOT))
    from combustion1d.exact import tw_profile_limit, tw_speed_limit
    from combustion1d.limit import solve_planar
    from experiments.traveling_wave import A, W0S, X0
    c_ref = tw_speed_limit(A, 0.0)
    L = X0 + tw_speed_limit(A, max(W0S)) * t_end + 15.0
    return solve_planar(lambda x: tw_profile_limit(x, 0.0, A, c_ref, x0=X0), w0, L, X0,
                        t_end, n=800, t_eval=np.linspace(0, t_end, 400))


def limit_figures():
    rs = json.loads((RES / "limit_comparison.json").read_text())
    geoms = [("planar", "plano"), ("radial2", "radial $N=2$"), ("radial3", "radial $N=3$")]

    # fronts vs time
    fig, axs = plt.subplots(1, 2, figsize=(8.4, 3.4))
    for ax, g in zip(axs, ("planar", "radial2")):
        sel = sorted((r for r in rs if r["geom"] == g and r["eps"] in (0.05, 0.01, 0.002)),
                     key=lambda r: -r["eps"])
        for k, r in enumerate(sel):
            ax.plot(r["t"], r["front"], color=BLUES[k + 1], label=f"$\\varepsilon={r['eps']:g}$")
        ax.plot(sel[0]["lim_t"], sel[0]["lim_s"], ls="--", lw=1.2, color=INK, label="límite $(P)$")
        ax.set_xlabel("$t$")
        ax.legend(loc="best")
    axs[0].set_ylabel("posición del frente")
    axs[0].set_title("Plano, $u_0=A\\tanh(\\sqrt{2M}(x-x_0)/A)^+$", fontsize=9, color=INK, loc="left")
    axs[1].set_title("Radial $N=2$, $u_0=a(1-r^2/R_0^2)^+$", fontsize=9, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(FIG / "limit_fronts.png")
    plt.close(fig)

    # errors vs eps
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.4), sharey=True)
    keys = [("front_err", r"$\\max_t|x_\\varepsilon(t)-s(t)|$"),
            ("front_err_shifted", r"idem, restando $\\varepsilon/\\sqrt{2M}$"),
            ("sup_err", r"$\\max_x|(u^\\varepsilon)^+-u|$ al final")]
    for k, (g, name) in enumerate(geoms):
        sel = sorted((r for r in rs if r["geom"] == g), key=lambda r: r["eps"])
        e = np.array([r["eps"] for r in sel])
        for ax, (key, _) in zip(axs, keys):
            ax.loglog(e, [r[key] for r in sel], color=SERIES[k], marker=MARKERS[k], ms=6,
                      mfc="white", mew=1.5, label=name)
    for ax, (_, title) in zip(axs, keys):
        slope_ref(ax, [1e-3, 5e-2], 1.0e-3, label=r"$\\varepsilon$")
        ax.set_title(title, fontsize=9, color=INK, loc="left")
        ax.set_xlabel(r"$\\varepsilon$")
        ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "limit_convergence.png")
    plt.close(fig)


if __name__ == "__main__":
    FIG.mkdir(exist_ok=True)
    if (RES / "traveling_wave.json").exists():
        tw_figures()
    if (RES / "radial.json").exists():
        radial_figures()
    if (RES / "limit_comparison.json").exists():
        limit_figures()
