"""Movie: a circular flame closing in on the fresh region (N = 2).

(P_eps) with w = eps*w0 constant is started from the Caffarelli--Vazquez
self-similar solution V(., 0) of the limit problem (same setup as
experiments/radial.py). The fresh region {u > 0} is a disk of radius
~ R (T - t)^{1/2} that shrinks to a point at t = T: the fresh mixture is
completely burnt in finite time.

Left: disk view of the temperature-like quantity -u (burnt = hot, fresh =
cold), with the numerical front {u = eps} (solid circle) and the exact free
boundary |x| = R (T - t)^{1/2} of the limit problem (dashed circle).
Right: radial profiles u^eps(r, t) (solid) and V(r, t) (dashed).

Run from numerics/:  python -m experiments.movie_radial
Writes figures/flame_extinction.{mp4,gif} and figures/flame_extinction_still.png
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter

from combustion1d.exact import SelfSimilar
from combustion1d.solver import simulate

N, W0, T, EPS = 2, 0.5, 1.0, 0.01
T_END = 0.99          # the solver runs robustly up to here (front radius ~ 25 eps)
N_FRAMES = 151
FPS = 15
FIG = Path(__file__).resolve().parent.parent / "figures"


def run():
    V = SelfSimilar(N, W0, T)
    L = 1.5 * V.R * np.sqrt(T) + 1.0
    times = np.linspace(0.0, T_END, N_FRAMES)
    res = simulate(lambda x: V(x, 0.0), W0, EPS, L, T_END, N=N,
                   speed_guess=V.R / (2 * np.sqrt(T)), snapshot_times=times,
                   max_dt_remesh=0.01)
    return V, times, res


def front_radius(r, u, level):
    """Radius where u crosses `level` (u decreasing through the front)."""
    i = np.nonzero((u[:-1] >= level) & (u[1:] < level))[0]
    if i.size == 0:
        return np.nan
    i = i[-1]
    return r[i] + (level - u[i]) * (r[i + 1] - r[i]) / (u[i + 1] - u[i])


def make_animation(V, times, res):
    R0 = V.radius(0.0)
    Lplot = 1.12 * R0
    n_pix = 501
    xs = np.linspace(-Lplot, Lplot, n_pix)
    X, Y = np.meshgrid(xs, xs)
    RHO = np.hypot(X, Y)
    umax = float(V(0.0, 0.0))
    vmin, vmax = -umax, EPS * W0          # fixed colour limits: -u in [-max u0, eps w0]
    rr = np.linspace(0.0, Lplot, 800)
    th = np.linspace(0.0, 2 * np.pi, 400)

    fig = plt.figure(figsize=(10, 5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], left=0.06, right=0.97,
                          bottom=0.11, top=0.83, wspace=0.28)
    ax = fig.add_subplot(gs[0])
    bx = fig.add_subplot(gs[1])

    def frame(k):
        t = times[k]
        r, u, _ = res.snapshots[t]
        theta = -np.interp(RHO, r, u)
        img.set_data(theta)
        rf = front_radius(r, u, EPS)
        rex = V.radius(t)
        c_num.set_data(rf * np.cos(th), rf * np.sin(th))
        c_ex.set_data(rex * np.cos(th), rex * np.sin(th))
        l_num.set_data(r, u)
        l_ex.set_data(rr, V(rr, t))
        title.set_text(f"t = {t:.3f}")
        info.set_text(f"front radius: (P$_\\varepsilon$) {rf:.3f},  exact {rex:.3f}")
        return img, c_num, c_ex, l_num, l_ex, title, info

    img = ax.imshow(np.zeros_like(RHO), origin="lower", cmap="inferno",
                    extent=(-Lplot, Lplot, -Lplot, Lplot), vmin=vmin, vmax=vmax,
                    interpolation="bilinear")
    c_num, = ax.plot([], [], color="#1f9bff", lw=1.8, label=r"front $\{u^\varepsilon=\varepsilon\}$")
    c_ex, = ax.plot([], [], color="black", lw=1.4, ls=(0, (4, 3)),
                    label=r"exact $|x| = R\sqrt{T-t}$")
    ax.set_aspect("equal")
    ax.set_xlim(-Lplot, Lplot)
    ax.set_ylim(-Lplot, Lplot)
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")
    ax.legend(loc="lower left", fontsize=8, framealpha=0.85)
    cb = fig.colorbar(img, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label(r"$-u$  (temperature-like: burnt hot, fresh cold)", fontsize=8)
    cb.ax.tick_params(labelsize=8)

    l_num, = bx.plot([], [], color="C0", lw=1.8, label=r"$u^\varepsilon(r,t)$, $\varepsilon=0.01$")
    l_ex, = bx.plot([], [], color="k", lw=1.3, ls="--", label=r"exact $V(r,t)$")
    bx.axhline(0.0, color="0.6", lw=0.6)
    bx.set_xlim(0.0, Lplot)
    bx.set_ylim(-0.09, 1.08 * umax)
    bx.set_xlabel("$r = |x|$")
    bx.set_ylabel("$u$")
    bx.legend(loc="upper right", fontsize=8)
    bx.grid(alpha=0.25)
    bx.text(0.97 * Lplot, -0.055, r"burnt region: $u = -\varepsilon w_0$", ha="right",
            va="center", fontsize=8, color="0.35")

    fig.suptitle(r"Circular flame: (P$_\varepsilon$) with $\varepsilon = 0.01$ (solid) vs. "
                 r"exact limit solution (dashed); extinction at $T = 1$",
                 fontsize=11, y=0.975)
    title = fig.text(0.5, 0.895, "", ha="center", fontsize=13, weight="bold")
    info = fig.text(0.5, 0.855, "", ha="center", fontsize=9, color="0.3")
    return fig, frame


def main():
    FIG.mkdir(exist_ok=True)
    V, times, res = run()
    m = res.t >= 0.05
    err = np.abs(res.front[m] - V.radius(res.t[m]))
    print(f"R = {V.R:.6f}, kappa = {V.kappa:.6f}, t_end = {T_END}, frames = {len(times)}")
    print(f"max |front - R sqrt(T-t)| for t in [0.05, {T_END}]: {err.max():.4e} "
          f"(all t: {np.nanmax(np.abs(res.front - V.radius(res.t))):.4e})")

    fig, frame = make_animation(V, times, res)
    frame(len(times) // 2)
    fig.savefig(FIG / "flame_extinction_still.png", dpi=120)
    anim = FuncAnimation(fig, frame, frames=len(times), blit=False)
    anim.save(FIG / "flame_extinction.mp4", dpi=120,
              writer=FFMpegWriter(fps=FPS, bitrate=-1,
                                  extra_args=["-pix_fmt", "yuv420p", "-crf", "20"]))
    print("wrote", FIG / "flame_extinction.mp4")
    anim.save(FIG / "flame_extinction.gif", dpi=60, writer=PillowWriter(fps=FPS))
    print("wrote", FIG / "flame_extinction.gif")


if __name__ == "__main__":
    main()
