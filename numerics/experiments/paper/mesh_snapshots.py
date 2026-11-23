"""Figure A5: the solution and the mesh at successive times.

For each time the panel shows u and Y, and below them one tick per mesh node
(the fine band appears as a solid block, the geometric transitions as
gradients, the coarse zone as sparse ticks). ALE scheme, h_f = eps/16,
h_c = 0.05 so that individual coarse nodes are visible.

Left:  planar travelling wave (w0 = 0, eps = 0.02), the band travels with the flame.
Right: radial extinction (N = 2, w0 = 0.5, eps = 0.02), self-similar data,
       the band follows the shrinking ball until it disappears.

Run from numerics/:  python -m experiments.paper.mesh_snapshots   (~10 s)
Figure: figures/paper/A5_mesh_snapshots.png
"""

import numpy as np

from combustion1d.ale import simulate_ale
from combustion1d.solver import MeshParams

from .common import BLUES, INK, MUTED, SERIES, extinction_setup, plt, savefig, title, tw_setup

EPS = 0.02
MP = MeshParams(h_coarse=0.05)


def run_tw():
    times = (0.0, 1.5, 3.0, 4.5)
    u0, L, c0 = tw_setup(0.0, times[-1])
    r = simulate_ale(u0, 0.0, EPS, L, times[-1], 0.01, mesh=MP, speed_guess=c0,
                     snapshot_times=times)
    return times, r, (0.0, 7.0), 0.0


def run_radial():
    times = (0.0, 0.3, 0.6, 0.85)
    V, u0, L = extinction_setup(2, 0.5)
    r = simulate_ale(u0, 0.5, EPS, L, times[-1], 0.0025, mesh=MP, speed_guess=-V.R / 2, N=2,
                     ell0=1.0, remesh_distance=0.1, snapshot_times=times)
    return times, r, (0.0, 4.2), 0.5


def panel(ax, x, u, w, xlim, label, label_left):
    m = (x >= xlim[0]) & (x <= xlim[1])
    ax.plot(x[m], u[m], color=SERIES[0], lw=1.8, label="$u$")
    ax.plot(x[m], (u + w)[m], color=SERIES[1], lw=1.2, ls="--", label="$Y$")
    lo, hi = -0.30, -0.06
    ax.vlines(x[m], lo, hi, color=INK, lw=0.35)
    ax.set_xlim(*xlim)
    ax.set_ylim(lo - 0.02, 1.15)
    ax.set_yticks([0, 1])
    if label_left:
        ax.text(0.012, 0.92, label, transform=ax.transAxes, ha="left", va="top", fontsize=8, color=INK)
    else:
        ax.text(0.988, 0.92, label, transform=ax.transAxes, ha="right", va="top", fontsize=8, color=INK)
    ax.text(0.012, 0.035, "nodes", transform=ax.transAxes, fontsize=6.5, color=MUTED,
            bbox=dict(facecolor="white", edgecolor="none", pad=0.6))
    ax.grid(False)


def zoom(ax, x, u, eps, bounds, burnt_left):
    """Inset: the reaction layer with one dot per node."""
    i = np.nonzero((u[:-1] - eps) * (u[1:] - eps) < 0)[0][0]
    xf = x[i]
    a, b = (xf - 9 * eps, xf + 4 * eps) if burnt_left else (xf - 4 * eps, xf + 9 * eps)
    m = (x >= a) & (x <= b)
    ins = ax.inset_axes(bounds)
    ins.plot((x[m] - xf) / eps, u[m] / eps, "-", color=SERIES[0], lw=0.6)
    ins.plot((x[m] - xf) / eps, u[m] / eps, "|", color=INK, ms=3.0, mew=0.5)
    ins.axhline(1.0, color=MUTED, lw=0.6, ls=":")
    ins.set_xlim((a - xf) / eps, (b - xf) / eps)
    ins.set_ylim(-0.8, 3.2)
    ins.set_yticks([0, 1, 2, 3])
    ins.tick_params(labelsize=6, length=2)
    ins.set_xlabel(r"$(x-x_f)/\varepsilon$" if burnt_left else r"$(r-r_f)/\varepsilon$", fontsize=6.5, labelpad=1)
    ins.set_ylabel(r"$u/\varepsilon$", fontsize=6.5, labelpad=1)
    ins.set_title("reaction layer: one tick per node", fontsize=6.5, color=INK, pad=2)
    ins.grid(False)


def main():
    cases = [("planar travelling wave, $w_0=0$", "x") + run_tw(),
             ("radial extinction, $N=2$, $w_0=0.5$", "r") + run_radial()]
    fig, axes = plt.subplots(4, 2, figsize=(10.5, 6.6), sharex="col")
    for col, (name, var, times, r, xlim, w0) in enumerate(cases):
        for row, t in enumerate(times):
            ax = axes[row, col]
            x, u, w = r.snapshots[t]
            panel(ax, x, u, w, xlim, f"$t={t:g}$,  {x.size} nodes", label_left=(col == 0))
            if row == 0:
                title(ax, f"{name}  ($\\varepsilon={EPS:g}$, $h_f=\\varepsilon/16$, $h_c={MP.h_coarse:g}$)")
                ax.legend(loc="center left" if col == 0 else "center right", fontsize=7, ncol=2,
                          bbox_to_anchor=(0.0, 0.55) if col == 0 else (1.0, 0.55))
            if row == len(times) - 1:
                zoom(ax, x, u, EPS, [0.33, 0.50, 0.28, 0.44] if col == 0 else [0.36, 0.50, 0.28, 0.44],
                     burnt_left=(col == 0))
        axes[-1, col].set_xlabel(f"${var}$")
    fig.tight_layout(h_pad=0.4)
    savefig(fig, "A5_mesh_snapshots.png")


if __name__ == "__main__":
    main()
