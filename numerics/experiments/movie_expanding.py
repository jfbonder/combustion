"""Movie: an expanding circular flame (N = 2).

(P_eps) with w = eps*w0 constant is started from a burnt kernel of radius r0
at the centre surrounded by fresh mixture,

    u0(r) = A tanh(sqrt(2M) (r - r0) / A)  for r > r0,   u0 = 0 for r <= r0,

which satisfies the free boundary condition |grad u0| = sqrt(2M) at r = r0.
The flame propagates outwards; its speed tends (from below, because of the
curvature 1/r) to the planar speed c0 = sqrt(2M)/A.

Note: the kernel must be larger than a critical radius (roughly 1/c0) or it
is quenched; with r0 = 0.5 the front moves inwards and the burnt region
disappears, so r0 = 1.5 is used here.

Left: disk view of -u (burnt = hot, fresh = cold) with the numerical front
{u = eps}. Right: front radius r_f(t) and the planar reference r0 + c0 t;
inset: radial profile u(r, t).

Not used in the paper. Uses the BDF solver (eps = 0.01, t <= 8.5).

Run from numerics/:  python -m experiments.movie_expanding
Writes figures/flame_expanding.mp4 and figures/flame_expanding_still.png.
Needs ffmpeg. Not run again for the current documentation: the run time was
not measured.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter, FuncAnimation

from combustion1d.exact import tw_speed_limit
from combustion1d.model import fb_slope
from combustion1d.solver import simulate

N, W0, EPS, A = 2, 0.5, 0.01, 1.0
R0 = 1.5              # initial kernel radius (0.5 is below the ignition radius)
L = 12.0
T_END = 8.5           # front reaches r ~ 8, i.e. stays > 3 away from r = L
N_FRAMES = 151
FPS = 15
LPLOT = 9.0
FIG = Path(__file__).resolve().parent.parent / "figures"

SLOPE = fb_slope(W0)
C0 = tw_speed_limit(A, W0)


def u0(r):
    return np.where(r > R0, A * np.tanh(SLOPE * (r - R0) / A), 0.0)


def run():
    times = np.linspace(0.0, T_END, N_FRAMES)
    res = simulate(u0, W0, EPS, L, T_END, N=N, speed_guess=C0,
                   snapshot_times=times, max_dt_remesh=0.02)
    return times, res


def make_animation(times, res):
    n_pix = 501
    xs = np.linspace(-LPLOT, LPLOT, n_pix)
    X, Y = np.meshgrid(xs, xs)
    RHO = np.hypot(X, Y)
    vmin, vmax = -A, EPS * W0             # fixed colour limits: -u in [-A, eps w0]
    th = np.linspace(0.0, 2 * np.pi, 400)
    m_ok = np.isfinite(res.front)
    tf, rf_all = res.t[m_ok], res.front[m_ok]

    fig = plt.figure(figsize=(10, 5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], left=0.06, right=0.97,
                          bottom=0.11, top=0.83, wspace=0.28)
    ax = fig.add_subplot(gs[0])
    bx = fig.add_subplot(gs[1])

    def frame(k):
        t = times[k]
        r, u, _ = res.snapshots[t]
        img.set_data(-np.interp(RHO, r, u))
        m = tf <= t + 1e-12
        rf = rf_all[m][-1]
        c_num.set_data(rf * np.cos(th), rf * np.sin(th))
        l_rf.set_data(tf[m], rf_all[m])
        p_rf.set_data([t], [rf])
        l_prof.set_data(r, u)
        p_prof.set_data([rf], [EPS])
        # front speed: least squares slope over the last 0.5 time units
        w = m & (tf >= t - 0.5)
        v = np.polyfit(tf[w], rf_all[w], 1)[0] if t >= 0.3 and w.sum() > 3 else np.nan
        title.set_text(f"t = {t:.2f}")
        info.set_text(f"front radius {rf:.3f},  front speed {v:.3f}  "
                      f"(planar speed $c_0$ = {C0:.3f})" if np.isfinite(v)
                      else f"front radius {rf:.3f}")
        return img, c_num, l_rf, p_rf, l_prof, p_prof, title, info

    img = ax.imshow(np.zeros_like(RHO), origin="lower", cmap="inferno",
                    extent=(-LPLOT, LPLOT, -LPLOT, LPLOT), vmin=vmin, vmax=vmax,
                    interpolation="bilinear")
    c_num, = ax.plot([], [], color="#1f9bff", lw=1.8,
                     label=r"front $\{u^\varepsilon=\varepsilon\}$")
    ax.set_aspect("equal")
    ax.set_xlim(-LPLOT, LPLOT)
    ax.set_ylim(-LPLOT, LPLOT)
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")
    ax.legend(loc="lower left", fontsize=8, framealpha=0.85)
    cb = fig.colorbar(img, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label(r"$-u$  (temperature-like: burnt hot, fresh cold)", fontsize=8)
    cb.ax.tick_params(labelsize=8)

    tt = np.linspace(0.0, T_END, 200)
    bx.plot(tt, R0 + C0 * tt, color="k", lw=1.2, ls="--",
            label=r"$r_0 + c_0 t$,  $c_0=\sqrt{2M}/A$")
    l_rf, = bx.plot([], [], color="C0", lw=1.9, label=r"front radius $r_f(t)$")
    p_rf, = bx.plot([], [], "o", color="C0", ms=5)
    bx.set_xlim(0.0, T_END)
    bx.set_ylim(0.0, R0 + C0 * T_END + 0.3)
    bx.set_xlabel("$t$")
    bx.set_ylabel("$r_f$")
    bx.legend(loc="lower right", fontsize=8)
    bx.grid(alpha=0.25)

    ins = bx.inset_axes([0.13, 0.55, 0.33, 0.31])
    l_prof, = ins.plot([], [], color="C3", lw=1.4)
    p_prof, = ins.plot([], [], "o", color="#1f9bff", ms=3.5)
    ins.axhline(0.0, color="0.6", lw=0.6)
    ins.set_xlim(0.0, L)
    ins.set_ylim(-0.12, 1.08 * A)
    ins.tick_params(labelsize=7)
    ins.text(0.97, 0.06, "$r$", transform=ins.transAxes, ha="right", fontsize=8)
    ins.set_title(r"profile $u^\varepsilon(r,t)$", fontsize=7, pad=2)
    ins.patch.set_alpha(0.9)

    fig.suptitle(r"Expanding circular flame, (P$_\varepsilon$) with $\varepsilon = 0.01$; "
                 r"dashed: planar speed $\sqrt{2M}/A$", fontsize=11, y=0.975)
    title = fig.text(0.5, 0.895, "", ha="center", fontsize=13, weight="bold")
    info = fig.text(0.5, 0.855, "", ha="center", fontsize=9, color="0.3")
    return fig, frame


def main():
    FIG.mkdir(exist_ok=True)
    times, res = run()
    m = np.isfinite(res.front)
    t, f = res.t[m], res.front[m]
    w = t >= T_END - 0.5
    v_end = np.polyfit(t[w], f[w], 1)[0]
    print(f"c0 = {C0:.6f}, r0 = {R0}, t_end = {T_END}, frames = {len(times)}")
    print(f"final front radius {f[-1]:.4f}; speed over [t_end-0.5, t_end] = {v_end:.4f} "
          f"(c0 - 1/r_f = {C0 - 1 / f[-1]:.4f})")

    fig, frame = make_animation(times, res)
    frame(len(times) // 2)
    fig.savefig(FIG / "flame_expanding_still.png", dpi=120)
    anim = FuncAnimation(fig, frame, frames=len(times), blit=False)
    anim.save(FIG / "flame_expanding.mp4", dpi=120,
              writer=FFMpegWriter(fps=FPS, bitrate=-1,
                                  extra_args=["-pix_fmt", "yuv420p", "-crf", "20"]))
    print("wrote", FIG / "flame_expanding.mp4")


if __name__ == "__main__":
    main()
