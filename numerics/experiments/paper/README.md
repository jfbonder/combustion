# Numerical experiments for the paper

Scripts, data and figures for the numerical section. Run everything from `numerics/`:

```bash
OMP_NUM_THREADS=1 ../.venv/bin/python -m experiments.paper.meshes          # A: run + plot (~3 min)
OMP_NUM_THREADS=1 ../.venv/bin/python -m experiments.paper.convergence     # B: run + plot (~2 min)
OMP_NUM_THREADS=1 ../.venv/bin/python -m experiments.paper.limit_examples  # C: run + plot (~12 min)
```

Each script accepts `run` or `plot`; with no argument it does both. `limit_examples run C1 C3` runs only part of group C. Runs use 4 processes.

- `common.py` holds the test problems, I/O, mesh analysis and the palette of `experiments/plots.py`.
- Results are in `results/paper/`: `*.json`, plus `A_mesh_*.npz` with the full node histories.
- Figures are in `figures/paper/` (PNG, 200 dpi).

Model throughout: $f(s)=(1-s)^+$, $M(w_0)=(1+w_0)^3/6$, and $c_0=\sqrt{2M}/A$ is the planar speed of the limit problem.

Default mesh, `MeshParams()`:

- band of width $(K_b+K_f+2K_m)\varepsilon=55\varepsilon$ with $h_f=\varepsilon/16$;
- geometric transition with growth factor 1.1;
- coarse spacing $h_c=0.01$.

"ALE" is `combustion1d.ale.simulate_ale`: implicit mesh velocity, fixed $\Delta t$. "Eulerian" is either `solver.simulate` (adaptive BDF, mesh rebuilt around the front) or `simple.simulate_simple` (monotone scheme, fixed $\Delta t$), as stated in each case.

---

## A. Meshes

### A1 `A1_spacetime_nodes_tw.png`: node trajectories, travelling wave

**Caption.** Space–time diagram of the mesh nodes for the planar travelling wave ($w_0=0$, $A=1$, $\varepsilon=0.01$, $0\le t\le4$; $h_f=\varepsilon/16$, $h_c=0.04$, chosen so that individual coarse nodes are visible). Each node is drawn as its trajectory $x_i(t)$, coloured by the local mesh size; the orange line is the front $\{u=\varepsilon\}$ and the ticks on the right border mark remeshings.

- **Left, ALE ($\Delta t=0.01$).** The fine band, the geometric transition and the coarse nodes near the front all translate with the flame. One remeshing occurs (at $t=2.24$, after the mesh has moved $\tfrac14\ell_0$). The band nodes stay on their lattice; only the transition is rebuilt.
- **Right, Eulerian BDF solver.** The nodes are fixed, and the mesh is rebuilt 31 times. At each rebuild the band jumps ahead: refinement in front of the flame and coarsening behind it happen in the same rebuild.
- **Bottom row.** Zoom (blue box) on the fresh edge of the band and the transition.

Cost (whole run): ALE 400 steps, 0.8 s; Eulerian BDF 1044 steps, 3.4 s.

### A2 `A2_mesh_size.png`: local mesh size

**Caption.** Local mesh size $h_{i+1/2}=x_{i+1}-x_i$ (log scale) for the ALE solver with the default mesh ($h_f=\varepsilon/16$, $h_c=0.01$, $\varepsilon=0.01$) at four times. Triangles mark the front.

- **Left, planar travelling wave ($\Delta t=0.01$).** The plateau $h=h_f$ (band) and the geometric ramps (transition) move with the front.
  - Behind the flame ($x\lesssim1.5$) the cells are stretched up to $6h_c$. This is the coarsening produced by the decay zone of the mesh velocity ($\chi$ goes linearly from 1 to 0 over $\ell_0$), next to the fixed end $x=0$.
  - Ahead of the flame, cells are slightly compressed ($h<h_c$) until the next remeshing.
- **Right, radial extinction ($N=2$, self-similar data, $w_0=0.5$, $\ell_0=1$, remeshing every 0.1, $\Delta t=0.005$).** The band follows the shrinking ball. The burnt side of the band (35ε) faces outwards.

### A3 `A3_node_count.png`: number of nodes, cost

**Caption.** Planar travelling wave, $w_0=0$, domain $[0,27.09]$, default mesh, $t\le4$.

- **(a) Nodes vs time**, ALE (solid) and Eulerian BDF (dashed), $\varepsilon=0.01$ and $0.001$. The count is constant to within one node.
- **(b) Nodes per zone vs ε** (markers), with the formula (dotted):
  - band $=(K_b+K_f+2K_m)\,\mathrm{pts\_per\_eps}+1=881$, independent of ε;
  - transition $=2\lceil\log(h_c/h_f)/\log 1.1\rceil$, growing like $\log(1/\varepsilon)$;
  - coarse $\approx(L-\text{band}-\text{transition})/h_c$.
- **(c) Time steps to reach $t=4$.** ALE ($\Delta t=0.01$) needs 400 steps for every ε. The adaptive BDF Eulerian solver needs a number of steps that grows like $\varepsilon^{-1.4}$.

| ε | total | band | transition | coarse | ALE steps / s | Eulerian BDF steps / remeshings / s |
|---|---|---|---|---|---|---|
| 0.01 | 3573 | 882 | 58 | 2633 | 400 / 1.6 | 932 / 31 / 6.1 |
| 0.005 | 3615 | 882 | 72 | 2661 | 400 / 1.3 | 2018 / 62 / 12.8 |
| 0.002 | 3650 | 882 | 90 | 2678 | 400 / 1.6 | 6533 / 155 / 44.0 |
| 0.001 | 3671 | 882 | 106 | 2683 | 400 / 1.0 | 22687 / 309 / 129.8 |

The formula gives 3573, 3614, 3650 and 3670 nodes. Going from ε = 0.01 to 0.001 adds 98 nodes, all in the transition.

### A4 `A4_radial_meshes.png`: radial node trajectories

**Caption.** Node trajectories $r_i(t)$ for the ALE solver in radial symmetry ($N=2$, $w_0=0.5$, $\varepsilon=0.01$, $h_f=\varepsilon/16$, $h_c=0.04$), coloured by the local mesh size. Orange: front. Ticks on the right: remeshings.

- **Left, extinction of the self-similar ball** ($T=1$, $t\le0.9$, $\Delta t=0.005$, $\ell_0=1$, remeshing distance 0.1). 180 steps, 16 remeshings. The band accelerates with the front. The nodes of the decay zone are dragged and reset at each remeshing, which gives the saw-tooth pattern.
- **Right, expanding flame** from a burnt kernel of radius $r_0=1.5$ ($u_0=\tanh(\sqrt{2M}(r-1.5))^+$, $A=1$, $L=12$, $\Delta t=0.01$, $\ell_0$ default). 600 steps, 6 remeshings. The front first slows down (curvature), then approaches the planar speed $c_0=1.06$.

With $h_c=0.01$ the same runs have about 1340 and 2060 nodes, constant in time.

---

## B. Convergence

How speeds were measured:

- Speeds are compared with the exact speed $c_\varepsilon$ of the travelling wave of $(P_\varepsilon)$, computed by ODE shooting. For example $c_\varepsilon=0.579445,\,0.577559$ ($w_0=0$) and $1.628511,\,1.632539$ ($w_0=1$) for $\varepsilon=0.01,\,0.001$.
- Initial datum: $A(1-e^{-c_\varepsilon(x-2)})^+$, which is exactly the fresh part of the ε-wave. With $c_0$ instead of $c_\varepsilon$ there is a slowly decaying transient (rate about $c^2/4$).
- Final time $T=12$ ($w_0=0$) or $T=8$ ($w_0=1$).
- ALE speed: mean speed between the first and the last remeshing in $[T/2,T]$, i.e. over whole remeshing cycles. The front position oscillates by about $\pm3\cdot10^{-4}$ over a cycle, and a least-squares slope on a window that is not aligned with the cycles is biased by about $10^{-4}$. The least-squares values are also stored (`err_lsq`).
- Eulerian speed: least-squares slope on $[T/2,T]$.

### B1 `B1_tw_speed_eulerian_vs_ale.png`: Eulerian vs ALE

**Caption.** Speed error $|c_h-c_\varepsilon|$ vs ε for the planar travelling wave, $w_0=0$ (left) and $w_0=1$ (middle), $\Delta t\in\{0.02,0.01,0.005\}$. Dashed squares: Eulerian monotone scheme (`simple.py`, same band meshes). Solid circles: ALE.

- **Right:** relative error vs $\Delta t/\varepsilon$, including $\Delta t=0.002$ and $0.001$ for the Eulerian scheme. The Eulerian error is a function of $\Delta t/\varepsilon$ alone (the points collapse onto one curve). It is linear for small $\Delta t/\varepsilon$ and saturates near the stalled front ($|c_h-c_\varepsilon|\to c_\varepsilon$). The ALE error does not depend on ε or $\Delta t$.
- ALE uses $\ell_0=$ default ($3/c_0=5.2$) for $w_0=0$ and $\ell_0=4$ for $w_0=1$; see the remark below the table.

Signed errors $c_h-c_\varepsilon$:

| scheme | w0 | Δt | ε=0.01 | ε=0.005 | ε=0.002 | ε=0.001 |
|---|---|---|---|---|---|---|
| Eulerian monotone | 0 | 0.02 | -2.11e-01 | -2.91e-01 | -3.87e-01 | -4.50e-01 |
| Eulerian monotone | 0 | 0.01 | -1.44e-01 | -2.11e-01 | -3.18e-01 | -3.88e-01 |
| Eulerian monotone | 0 | 0.005 | -9.12e-02 | -1.44e-01 | -2.40e-01 | -3.18e-01 |
| Eulerian monotone | 1 | 0.02 | -9.22e-01 | -1.10e+00 | -1.28e+00 | -1.38e+00 |
| Eulerian monotone | 1 | 0.01 | -7.53e-01 | -9.26e-01 | -1.15e+00 | -1.28e+00 |
| Eulerian monotone | 1 | 0.005 | -5.65e-01 | -7.56e-01 | -9.99e-01 | -1.15e+00 |
| ALE | 0 | 0.02 | +1.66e-04 | +1.78e-04 | +1.82e-04 | +1.83e-04 |
| ALE | 0 | 0.01 | +1.49e-04 | +1.61e-04 | +1.64e-04 | +1.67e-04 |
| ALE | 0 | 0.005 | +1.42e-04 | +1.52e-04 | +1.51e-04 | +1.55e-04 |
| ALE (ℓ0=4) | 1 | 0.02 | +1.16e-04 | -2.31e-04 | +2.77e-05 | +2.67e-05 |
| ALE (ℓ0=4) | 1 | 0.01 | +1.29e-04 | -2.30e-04 | +5.46e-06 | -2.55e-06 |
| ALE (ℓ0=4) | 1 | 0.005 | -1.90e-05 | -5.08e-05 | -1.45e-04 | +4.98e-05 |
| ALE, default ℓ0=2 | 1 | 0.02 | +7.48e-04 | +9.64e-04 | +1.11e-03 | +1.23e-03 |
| ALE, default ℓ0=2 | 1 | 0.01 | +3.89e-04 | +4.68e-04 | +5.16e-04 | +6.00e-04 |
| ALE, default ℓ0=2 | 1 | 0.005 | +1.94e-04 | +2.19e-04 | +2.28e-04 | +2.88e-04 |

The Eulerian errors are identical along diagonals with equal $\Delta t/\varepsilon$. For example, $-0.144$ is obtained for $(\Delta t,\varepsilon)=(0.01,0.01)$, $(0.005,0.005)$, $(0.002,0.002)$ and $(0.001,0.001)$.

**Remark on $\ell_0$ for $w_0=1$.** The code default $\ell_0=\max(2,3/c)$ gives $\ell_0=2$ for $c=1.63$. This leaves an $O(\Delta t)$ speed error of about $0.06\,\Delta t$, which grows mildly as ε decreases. With $\ell_0=4$ the error is flat in $\Delta t$.

### B2 `B2_ale_mesh_sizes.png`: ALE error vs $h_f$ and $h_c$

**Caption.** ALE speed error for the travelling wave ($w_0=0$, $\Delta t=0.01$).

- **Left:** vs $h_f/\varepsilon$ at $h_c=0.0025$, for four values of ε. The curves for different ε coincide: ε enters only through the band mesh size $h_f/\varepsilon$. The error decreases about linearly (ratio 2–2.5 per halving) down to the noise level of the speed measurement, about $2\cdot10^{-5}$ (shaded).
- **Right:** vs $h_c$ at $h_f=\varepsilon/32$. The coarse-mesh error is about $O(h_c^2)$: $5\cdot10^{-5}$ at $h_c=0.04$, $10^{-5}$ at $0.02$. With $\varepsilon=0.001$ it saturates at the $h_f$ level, $1.6\cdot10^{-5}$.

| pts_per_eps ($h_f/\varepsilon$) | ε=0.01 | ε=0.005 | ε=0.002 | ε=0.001 |
|---|---|---|---|---|
| 4 (1/4) | +1.54e-03 | +8.65e-04 | +8.55e-04 | +8.72e-04 |
| 8 (1/8) | +3.70e-04 | +4.01e-04 | +4.12e-04 | +4.20e-04 |
| 16 (1/16) | +1.55e-04 | +1.64e-04 | +1.68e-04 | +1.69e-04 |
| 32 (1/32) | +1.08e-06 | -7.79e-06 | -1.36e-05 | -1.56e-05 |
| 64 (1/64) | +1.52e-05 | +2.16e-05 | +2.51e-05 | +2.61e-05 |

| $h_c$ | ε=0.01 | ε=0.001 |
|---|---|---|
| 0.04 | -5.10e-05 | -6.65e-05 |
| 0.02 | -1.02e-05 | -2.63e-05 |
| 0.01 | -1.18e-07 | -1.57e-05 |
| 0.005 | +7.60e-07 | -1.65e-05 |
| 0.0025 | +1.08e-06 | -1.56e-05 |

The $\approx1.6\cdot10^{-4}$ floor of the default ALE runs (B1, $w_0=0$) therefore comes from $h_f=\varepsilon/16$, not from $h_c$. The background notes attributed it to $h_c$.

### B3 `B3_rigid_mesh_dt.png`: rigid mesh at the exact speed

**Caption.** Whole mesh translating rigidly with $v=c_\varepsilon$ (`fixed_velocity=c_eps`, `ell0=1e6`, `remesh_distance=1e9`), $h_f=\varepsilon/16$, $h_c=0.01$. Front travel is about 6 ($T=10$ for $w_0=0$, $T=3.6$ for $w_0=1$). The speed error shows no trend in $\Delta t$ from 0.1 to 0.005, although $\Delta t=0.1$ is $10^5\varepsilon^2$ for $\varepsilon=10^{-3}$ (Theorem "ALE, travelling wave").

| w0 | ε | Δt=0.1 | 0.05 | 0.02 | 0.01 | 0.005 |
|---|---|---|---|---|---|---|
| 0 | 0.005 | -1.8e-06 | +3.0e-07 | +1.1e-06 | +1.0e-06 | +9.9e-07 |
| 0 | 0.001 | -5.1e-07 | -1.3e-07 | +3.1e-08 | +3.1e-08 | +2.9e-08 |
| 1 | 0.005 | -2.7e-05 | -7.7e-06 | +5.1e-06 | +7.4e-06 | +2.2e-06 |
| 1 | 0.001 | -3.5e-06 | -4.9e-06 | -4.0e-07 | +1.1e-07 | +1.0e-07 |

### B4 `B4_radial_extinction_dt.png`: radial extinction, error vs Δt

**Caption.** ALE, exact self-similar extinction ($w_0=0.5$, $T=1$, $t\le0.75$, $L=1.5R+1$, $\ell_0=1$, remeshing distance 0.1, default mesh). The plot shows the front error after removing the level-set shift $\varepsilon/\sqrt{2M}$,
$\max_{t\ge0.05}|r_h(t)+\varepsilon/\sqrt{2M}-R\sqrt{T-t}|$, vs $\Delta t$. The error is first order in $\Delta t$ with a constant independent of ε: $\approx1.35\,\Delta t$ for $N=2$ and $\approx1.55\,\Delta t$ for $N=3$. Each run takes under 1 s.

| N | ε | Δt=0.04* | 0.02* | 0.01 | 0.005 | 0.0025 | error/Δt at 0.0025 |
|---|---|---|---|---|---|---|---|
| 2 | 0.01 | 6.23e-02 | 3.01e-02 | 1.54e-02 | 8.16e-03 | 4.59e-03 | 1.84 |
| 2 | 0.002 | 5.82e-02 | 2.78e-02 | 1.39e-02 | 6.94e-03 | 3.53e-03 | 1.41 |
| 2 | 0.001 | 5.68e-02 | 2.71e-02 | 1.35e-02 | 6.71e-03 | 3.37e-03 | 1.35 |
| 3 | 0.01 | 6.90e-02 | 3.34e-02 | 1.71e-02 | 8.97e-03 | 4.97e-03 | 1.99 |
| 3 | 0.002 | 6.59e-02 | 3.15e-02 | 1.57e-02 | 7.91e-03 | 4.05e-03 | 1.62 |
| 3 | 0.001 | 6.46e-02 | 3.09e-02 | 1.53e-02 | 7.64e-03 | 3.85e-03 | 1.54 |

\* The actual steps are $0.75/19=0.0395$ and $0.75/38=0.0197$. The small extra error at ε = 0.01 is the $O(\varepsilon)$ part (sup error in the interior: $3.4\cdot10^{-3}$ at Δt = 0.0025). The JSON also has the sup errors `sup_err` and `sup_err_int`.

---

## C. Limit problem and examples from the literature

### C1 `C1_ale_vs_limit.png`: $(P_\varepsilon)$ vs the computed limit problem

**Caption.** Non-self-similar radial data $u_0=a(1-r^2/R_0^2)^+$, $a=\sqrt{2M}R_0/2$, $R_0=2$, $w_0=0.5$, computed up to the time when the limit radius is $R_0/2$ ($t=0.565$ for $N=2$, $0.393$ for $N=3$). The limit problem is solved by front fixing (`limit.solve_radial`, $n=1600$).

- **Left:** fronts for $N=2$.
- **Middle:** front error after removing $\varepsilon/\sqrt{2M}$.
- **Right:** $\max_r|u_h^+-u|$ at the final time.
- Colour: ALE with $\Delta t=0.002$ (dotted), $0.0005$ (dashed) and $0.000125$ (solid). Grey: Eulerian BDF solver (`results/limit_comparison.json`).

The ALE error behaves like $a\varepsilon+b\Delta t$ with $b\approx2.3$–$2.6$. For $\varepsilon\ge0.01$ the ε part dominates and ALE agrees with the BDF solver. For small ε the $\Delta t$ part dominates: with $\Delta t=1.25\cdot10^{-4}$ the floor is about $3.5\cdot10^{-4}$.

| N | ε | Δt | front err | shifted | sup | steps |
|---|---|---|---|---|---|---|
| 2 | 0.05 | 1.25e-4 | 5.66e-02 | 9.44e-03 | 6.80e-03 | 4520 |
| 2 | 0.01 | 1.25e-4 | 1.07e-02 | 1.26e-03 | 1.21e-03 | 4520 |
| 2 | 0.005 | 1.25e-4 | 5.41e-03 | 6.91e-04 | 7.02e-04 | 4520 |
| 2 | 0.002 | 1.25e-4 | 2.32e-03 | 4.39e-04 | 4.49e-04 | 4520 |
| 2 | 0.001 | 1.25e-4 | 1.32e-03 | 3.75e-04 | 3.85e-04 | 4520 |
| 2 | 0.001 | 5e-4 | 2.00e-03 | 1.06e-03 | 1.12e-03 | 1130 |
| 2 | 0.001 | 2e-3 | 4.76e-03 | 3.82e-03 | 4.03e-03 | 282 |
| 3 | 0.01 | 1.25e-4 | 1.07e-02 | 1.29e-03 | 1.23e-03 | 3144 |
| 3 | 0.001 | 1.25e-4 | 1.37e-03 | 4.23e-04 | 4.40e-04 | 3144 |
| BDF, N=2 | 0.01 | adaptive | 1.04e-02 | 9.59e-04 | 9.02e-04 | 1024 |
| BDF, N=2 | 0.001 | adaptive | 1.01e-03 | 6.74e-05 | 9.27e-05 | 14408 |

The unshifted front error is $\approx\varepsilon/\sqrt{2M}$ in all cases, i.e. the distance between the tracked level set $\{u=\varepsilon\}$ and the free boundary.

### C2 `C2_self_similar_convergence.png`: approach to the self-similar extinction profile

**Caption.** Same data as C1, followed until extinction.

- **Top row:** $s(t)/\sqrt{T-t}$ vs $T-t$ (log scale, time running to the right).
  - Black: limit problem, with $T$ fitted from $s^2\sim R^2(T-t)$ on $s<0.3$: $T=0.725090$ ($N=2$), $0.505756$ ($N=3$).
  - Dashed: $(P_\varepsilon)$ with ALE, $\varepsilon=0.01$ ($\Delta t=10^{-4}$) and $0.002$ ($\Delta t=2.5\cdot10^{-5}$). For these, $r_h+\varepsilon/\sqrt{2M}$ is used, with $T_\varepsilon$ fitted in the same way.
  - Dotted: $\varepsilon=0.002$ with the Eulerian BDF solver.
  - Domain $L=R_0+0.5$; see issue 2 below.
  - The ratio tends to the Galaktionov–Hulshof–Vázquez constant $R$ (2.513927 for $N=2$, 3.003951 for $N=3$).
- **Bottom row:** rescaled profiles $u/\sqrt{T-t}$ vs $\eta=r/\sqrt{T-t}$ when the radius is 1.5, 1.0, 0.5 and 0.25 (limit problem solid, ε = 0.002 dashed). They converge to $\kappa\phi(\eta)$ (thick grey).

Limit problem, $s/\sqrt{T-t}$:

| s | N=2 | N=3 |
|---|---|---|
| 1.5 | 2.4474 | 2.9202 |
| 1.0 | 2.4922 | 2.9716 |
| 0.5 | 2.5102 | 2.9971 |
| 0.25 | 2.5130 | 3.0020 |
| R (exact) | 2.513927 | 3.003951 |

Extinction times:

| | $T$ (limit) | BDF, ε=0.002 | ALE, ε=0.01 | ALE, ε=0.002 |
|---|---|---|---|---|
| N=2 | 0.725090 | 0.725024 | 0.72459 | 0.72525 |
| N=3 | 0.505756 | 0.505712 | 0.50540 | 0.50589 |

The small remaining dip of the ALE curve for ε = 0.002 (about 1% at $T-t\approx10^{-2}$) comes from the mesh-velocity taper near $r=0$ (issue 2). The BDF curve lies on the limit curve.

### C3 `C3_quenching.png`: quenching of small burnt kernels

**Caption.** Expanding flame from a burnt kernel of radius $r_0$: $u_0=\tanh(\sqrt{2M}(r-r_0))^+$, $A=1$, $w_0=0.5$, $N=2$, $c_0=1.0607$.

- **Left:** front radius $r_f(t)$ for several $r_0$ ($\varepsilon=0.01$, ALE $\Delta t=0.001$). Kernels that ignite (orange) first shrink, down to $r\approx0.52$ for $r_0$ just above critical, then expand with speed tending to $c_0$. Smaller kernels collapse (blue).
- **Right:** critical radius $r_*$ found by bisection, run to $t=6$ (to $t=12$ when undecided; tolerance $10^{-3}$).
  - $r_*=0.715\pm0.001$ for every ε from 0.01 to 0.00125.
  - The Eulerian BDF solver gives the same value (ε = 0.01), and so does ALE with $\Delta t=2.5\cdot10^{-4}$.
  - The heuristic $1/c_0=0.943$ (stationary radius of the curvature-corrected speed $c_0-1/r$) overestimates $r_*$ by 30%.
  - For these data the threshold is set by the transient (a kernel that barely ignites first shrinks to $r\approx0.52$), so $r_*$ presumably depends on the whole initial profile, not only on $r_0$ (not tested).

| ε | Δt | $r_*\in$ |
|---|---|---|
| 0.01 | 0.001 | [0.7148, 0.7156] |
| 0.005 | 0.001 | [0.7148, 0.7156] |
| 0.0025 | 0.001 | [0.7148, 0.7156] |
| 0.00125 | 0.001 | [0.7141, 0.7148] |
| 0.00125 | 0.00025 | [0.7156, 0.7164] |
| 0.01 | BDF (Eulerian) | [0.7148, 0.7156] |

With $\Delta t=0.005$ (file `C3_quenching_dt0.005.json`) the bisection gives 0.7129, 0.7129, 0.7105 and 0.7059 for ε = 0.01, 0.005, 0.0025 and 0.00125. This apparent drift with ε is an artifact of issue 2: near-critical fronts dip to $r\approx0.5$, inside the taper zone, where the lag is $O(\Delta t\,v/\varepsilon)$.

---

## Issues found in `combustion1d/` (not modified)

1. **`simulate_ale` crashes when the front disappears.**
   - What happens: when the front vanishes (kernel quenching, or extinction), `find_front` returns `None` and the step fails with `TypeError: unsupported operand type(s) for -: 'float' and 'NoneType'` (`abs(xf - xc)` in the remeshing test, `ale.py` line 149). Everything computed so far is lost.
   - Workaround: `common.simulate_ale_safe` swaps in a `Result` subclass that keeps the partial result, including snapshots and mesh history. It also records the remeshing steps.

2. **The mesh velocity is tapered to zero near $r=0$.**
   - Cause: `chi(z)` is multiplied by `clip(z/ell_end)`, with `ell_end = min(0.5, 0.1 L)` (`ale.py` line 80).
   - Effect: for $r_f<$ `ell_end` the scheme is partly Eulerian, and the $O(\Delta t\,v/\varepsilon)$ front lag of the Eulerian scheme comes back. The front speed $v\sim R^2/(2r)$ blows up near extinction.
   - Example: SelfSimilar, $N=2$, $\varepsilon=0.001$, $\Delta t=0.001$, $L=4.77$. At $t=0.99$ the front is at $r_h=0.358$ against the exact $0.251$; the error was $<0.006$ up to $t=0.965$. With $\varepsilon=0.01$ (same $\Delta t$) the error at $t=0.99$ is 0.022.
   - Example: in C2 ($\varepsilon=0.002$, $\Delta t=2.5\cdot10^{-5}$) the extrapolated extinction time is off by $8.5\cdot10^{-4}$ with $L=4$, and by $1.6\cdot10^{-4}$ with $L=2.5$ (taper 0.25).
   - Workaround: shorter domain, smaller $\Delta t$.
   - Suggested fix: taper over a few $h_c$ near $r=0$, or let the band approach the origin.

3. **The default $\ell_0=\max(2,3/c)$ is too short for fast fronts.**
   - $w_0=1$, $c=1.63$ gives $\ell_0=2$ and an $O(\Delta t)$ speed error of about $0.06\,\Delta t$ (B1 table; up to $1.2\cdot10^{-3}$ at $\Delta t=0.02$).
   - $\ell_0=4$ removes it. $\ell_0\approx6.5/c$ seems needed.

4. **The rigid mode is only usable for moderate travel distances.**
   - Rigid mode means `ell0=1e6, remesh_distance=1e9`.
   - The nodes in the tapered end zones are compressed exponentially at $x=L$ and stretched at $x=0$.
   - For $w_0=1$, $T=8$ (travel 13) the cells at $x=L$ collapse and `solve_banded` fails with "array must not contain infs or NaNs".
   - B3 therefore limits the travel to about 6.

5. **Speed measurement.**
   - The ALE front position oscillates by about $\pm3\cdot10^{-4}$ ($w_0=0$, ε = 0.01) over each remeshing cycle, because the decay zone of $\chi$ distorts the mesh and is then reset.
   - `front_speed` (least squares) on a window that is not a whole number of cycles is then biased by about $10^{-4}$, which is the size of the "floor" quoted before.
   - Use `common.cycle_speed`.

## Runtimes (4 processes, `OMP_NUM_THREADS=1`)

**A (~3 min).** ALE runs take 0.3–1.8 s each. The Eulerian BDF runs take 6, 13, 44 and 130 s for ε = 0.01, 0.005, 0.002 and 0.001.

**B (125 s in total, 160 runs).** The largest ALE run (pts_per_eps = 64, $h_c=0.0025$, about 12,000 nodes) takes 13 s.

**C:**

| part | runtime |
|---|---|
| C1 (36 runs) | ~1 min, at most 8.5 s per run |
| C2 | ~1.5 min (ALE ε = 0.002 with Δt = 2.5e-5 takes 48 s) |
| C3 | ~7 min (bisection: 8 runs of about 10 s each per ε; 47 s each at Δt = 2.5e-4) |
