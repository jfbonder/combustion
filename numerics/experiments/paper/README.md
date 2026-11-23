# Numerical experiments for the paper

Scripts, data and figures for the numerical section. Run everything from `numerics/`:

```bash
OMP_NUM_THREADS=1 ../.venv/bin/python -m experiments.paper.meshes          # A: run + plot (~4 min)
OMP_NUM_THREADS=1 ../.venv/bin/python -m experiments.paper.convergence     # B: run + plot (~2 min)
OMP_NUM_THREADS=1 ../.venv/bin/python -m experiments.paper.limit_examples  # C: run + plot (~12 min)
OMP_NUM_THREADS=1 ../.venv/bin/python -m experiments.paper.tw_error        # D, F: run + plot (~6 min)
OMP_NUM_THREADS=1 ../.venv/bin/python -m experiments.paper.tab_simple      # table tab:simple (~3 min)
```

Each group script accepts `run` or `plot`; with no argument it does both. `limit_examples run C1 C3` runs only part of group C, `tw_error run D3 D4` only part of group D. Runs use 4 processes.

- `common.py` holds the test problems, I/O, mesh analysis and the palette of `experiments/plots.py`.
- Results are in `results/paper/`: `*.json`, plus `A_mesh_*.npz` with the full node histories.
- Figures are in `figures/paper/` (PNG, 200 dpi).

Model throughout: $f(s)=(1-s)^+$, $M(w_0)=(1+w_0)^3/6$, and $c_0=\sqrt{2M}/A$ is the planar speed of the limit problem.

Default mesh, `MeshParams()`:

- band of width $(K_b+K_f+2K_m)\varepsilon=55\varepsilon$ with $h_f=\varepsilon/16$;
- geometric transition with growth factor 1.1;
- coarse spacing $h_c=0.01$.

## Which solver produced what

| figure / table | solver |
|---|---|
| A1 left, A2, A4, A3 "ALE" | scheme (S) on the moving mesh, `ale.simulate_ale` |
| A1 right, A3 "Eulerian BDF" | adaptive BDF, PCHIP transfer, `solver.simulate` (not covered by the theorems) |
| B1 "Eulerian", D5, `tab:simple` | scheme (S) on rebuilt band meshes, `simple.simulate_simple`, fixed $\Delta t$ |
| B1 "ALE", B2, B3, B4, D1–D4, F | `ale.simulate_ale` |
| C1, C2, C3 coloured curves / ALE rows | `ale.simulate_ale` |
| C1 grey, C2 dotted "BDF", C3 row "BDF (Eulerian)" | `solver.simulate` |
| `figures/*.png` outside `paper/`, `results_*.log`, `results/*.json` outside `paper/` (Table `tab:tw`, Figs `fig:tw`, `fig:instability`, `fig:limit`, Table `tab:limit`) | `solver.simulate` (not rerun) |

---

## What changed (November 2026, after the referee report)

All ALE-based numbers in this file were recomputed with the code below. Old values are quoted as "(before: …)".

**Code changes in `combustion1d/`.**

1. **Band anchor at remeshing** (`ale.py`; committed just before this rerun). The new band is anchored on a node whose right cell has exactly the spacing $h_f$. Before, it was anchored on a transition node, 0.35 $h_f$ off the lattice, so linear interpolation acted inside the reaction layer at every remeshing (mesh-velocity spikes of 3–4e-3).
2. **Initial reference point** (`ale.py`, new). The reference point $x_c$ of the implicit velocity is now the front of the initial datum *on the actual mesh*. Before, it was the front located on the auxiliary uniform grid `linspace(0, L, 20001)`, whose spacing $L/20000\approx1.4\cdot10^{-3}$ is larger than ε for small ε. The first implicit step then moved the mesh by that location error, i.e. $v_1-c_\varepsilon=(\text{error})/\Delta t$: this was the "first step 0.29" of the referee (M7a). After the fix, starting from $\Phi_\varepsilon$: $|v_1-c_\varepsilon|=5\cdot10^{-4}$ at $\Delta t=0.02$ (any ε).
   This fix also fixes the lattice phase of the front for the whole run (see D3), and it changed all ALE travelling-wave speed errors.
3. **Rigid zone** (`ale.py`, new default `rigid_extra=1.0`). The mesh velocity is $v$ on the band plus $\ell_0$ on each side, then decays linearly to 0 over $\ell_0$: $R_1=\text{band}+\ell_0$ (before: $\text{band}+\tfrac12\ell_0$), $\ell_0=\max(2,3/c)$ unchanged. See D1 and issue 3.
4. **`simple.py` uses the nominal Δt exactly.** The mesh is rebuilt every $\lfloor\text{interval}/\Delta t\rfloor$ steps; only the last step (to land on $T$ or on a snapshot time) can be shorter. The step sizes are recorded in `res.dts`. Before, each interval $7.5\varepsilon/v_{est}$ was split into equal steps $\le\Delta t$ (mean step 0.0091 for nominal 0.01 at ε = 0.001).
5. New: `exact.tw_profile_eps(z, A, w0, eps)`, the travelling wave $\Phi_\varepsilon$ of $(P_\varepsilon)$ with $\Phi_\varepsilon(0)=\varepsilon$ (test `test_tw_profile_eps_solves_the_ode`); options `exact_initial` and `monitor` of `simulate_ale`/`simulate_simple`; counters `res.n_trials`, `res.n_solves` (scheme (S)), `res.n_lu`, `res.n_fev` (BDF); `res.v` (mesh velocity of each ALE step).

**Summary of the effect.**

- The ALE speed errors in B1/B2 changed (e.g. 1.66e-4 → 1.66e-5 for $w_0=0$, ε = 0.01, Δt = 0.02), and the three speed estimators now agree to 2–3 digits (table F).
- The remaining ALE speed error is a **lattice-pinning error**: with implicit velocity the front stays at a fixed position $\theta h_f$ inside a band cell, and the realised speed is $c_\varepsilon+g(\theta)$ with $g$ of zero mean, amplitude $\propto(h_f/\varepsilon)^2$, independent of ε, Δt, $h_c$, remeshing (D3). It is not the band truncation error of the semi-discrete scheme in the usual sense, and it is not $O(h_f/\varepsilon)$.
- Eulerian scheme (S): numbers changed by 0–4 % (exact Δt); the collapse on Δt/ε is now exact to 3 digits.
- C1 improved at small ε (floor 3.75e-4 → 7.7e-5); C2, C3 essentially unchanged.

---

## A. Meshes

### A1 `A1_spacetime_nodes_tw.png`: node trajectories, travelling wave

**Caption.** Space–time diagram of the mesh nodes for the planar travelling wave ($w_0=0$, $A=1$, $\varepsilon=0.01$, $0\le t\le4$; $h_f=\varepsilon/16$, $h_c=0.04$, chosen so that individual coarse nodes are visible). Each node is drawn as its trajectory $x_i(t)$, coloured by the local mesh size; the orange line is the front $\{u=\varepsilon\}$ and the ticks on the right border mark remeshings.

- **Left, ALE ($\Delta t=0.01$).** The fine band, the geometric transition and the coarse nodes near the front all translate with the flame. One remeshing occurs (after the mesh has moved $\tfrac14\ell_0$). The band nodes stay on their lattice; only the transition is rebuilt.
- **Right, Eulerian BDF solver.** The nodes are fixed, and the mesh is rebuilt 31 times. At each rebuild the band jumps ahead: refinement in front of the flame and coarsening behind it happen in the same rebuild.
- **Bottom row.** Zoom (blue box) on the fresh edge of the band and the transition.

Cost (whole run): ALE 400 steps, 2404 tridiagonal solves, 0.7 s; Eulerian BDF 1044 steps, 1305 sparse LU factorisations, 4849 RHS evaluations, 3.4 s. (Unchanged except for the counts, which are new.)

### A2 `A2_mesh_size.png`: local mesh size

**Caption.** Local mesh size $h_{i+1/2}=x_{i+1}-x_i$ (log scale) for the ALE solver with the default mesh ($h_f=\varepsilon/16$, $h_c=0.01$, $\varepsilon=0.01$) at four times. Triangles mark the front.

- **Left, planar travelling wave ($\Delta t=0.01$).** The plateau $h=h_f$ (band) and the geometric ramps (transition) move with the front.
  - Behind the flame the cells are stretched (up to several $h_c$). This is the coarsening produced by the decay zone of the mesh velocity, next to the fixed end $x=0$. With the larger rigid zone the decay zone now starts $\ell_0$ (before: $\ell_0/2$) behind the band.
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

| ε | total | band | transition | coarse | ALE steps / tridiag. solves / s | Eulerian BDF steps / LU / RHS evals / remeshings / s |
|---|---|---|---|---|---|---|
| 0.01 | 3573 | 882 | 58 | 2633 | 400 / 2404 / 1.1 | 932 / 1178 / 4512 / 31 / 5.5 |
| 0.005 | 3615 | 882 | 72 | 2661 | 400 / 2404 / 1.2 | 2018 / 2734 / 9913 / 62 / 11.7 |
| 0.002 | 3650 | 882 | 90 | 2678 | 400 / 2404 / 1.3 | 6533 / 8881 / 31773 / 155 / 39.4 |
| 0.001 | 3671 | 882 | 106 | 2683 | 400 / 2404 / 1.6 | 22687 / 30504 / 107801 / 309 / 122.1 |

The formula gives 3573, 3614, 3650 and 3670 nodes. Going from ε = 0.01 to 0.001 adds 98 nodes, all in the transition.

**Cost of the implicit velocity.** Each trial velocity costs one mesh/matrix assembly and two tridiagonal solves (n ≈ 3600). In every ALE run of this directory the secant iteration used exactly 3 trials per step (3.00–3.01 on average: the initial guess, one perturbed guess, one secant step, which meets the tolerance $10^{-3}h_f$), i.e. **6 tridiagonal solves per step** (2404 = 6·400 + 4). An Eulerian step of scheme (S) costs 2 tridiagonal solves. A BDF step costs on average 1.3 sparse LU factorisations of the $2n\times2n$ Jacobian and 4.8 RHS evaluations (Newton iterations, each with a sparse back-substitution). So "400 ALE steps" corresponds to 2400 tridiagonal solves.

### A4 `A4_radial_meshes.png`: radial node trajectories

**Caption.** Node trajectories $r_i(t)$ for the ALE solver in radial symmetry ($N=2$, $w_0=0.5$, $\varepsilon=0.01$, $h_f=\varepsilon/16$, $h_c=0.04$), coloured by the local mesh size. Orange: front. Ticks on the right: remeshings.

- **Left, extinction of the self-similar ball** ($T=1$, $t\le0.9$, $\Delta t=0.005$, $\ell_0=1$, remeshing distance 0.1). 180 steps, 16 remeshings. The band accelerates with the front. The nodes of the decay zone are dragged and reset at each remeshing, which gives the saw-tooth pattern.
- **Right, expanding flame** from a burnt kernel of radius $r_0=1.5$ ($u_0=\tanh(\sqrt{2M}(r-1.5))^+$, $A=1$, $L=12$, $\Delta t=0.01$, $\ell_0$ default). 600 steps, 6 remeshings. The front first slows down (curvature), then approaches the planar speed $c_0=1.06$.

With $h_c=0.01$ the same runs have about 1340 and 2060 nodes, constant in time.

---

## B. Convergence (speed of the planar wave; radial extinction)

How speeds were measured:

- Speeds are compared with the exact speed $c_\varepsilon$ of the travelling wave of $(P_\varepsilon)$, computed by ODE shooting. For example $c_\varepsilon=0.579445,\,0.577559$ ($w_0=0$) and $1.628511,\,1.632539$ ($w_0=1$) for $\varepsilon=0.01,\,0.001$.
- Initial datum: $A(1-e^{-c_\varepsilon(x-2)})^+$, the fresh part of the ε-wave (not $\Phi_\varepsilon$; group D uses $\Phi_\varepsilon$).
- Final time $T=12$ ($w_0=0$) or $T=8$ ($w_0=1$).
- ALE speed (`err`): mean speed between the first and the last remeshing in $[T/2,T]$. Also stored: least squares on $[T/2,T]$ (`err_lsq`) and end points $(x_f(T)-x_f(T/2))/(T/2)$ (`err_end`). **After the fixes the three agree to 2–3 digits in every B1/B2 run** (before: they differed by up to a factor 4 and in sign, referee M8(c)).
- Eulerian speed: least-squares slope on $[T/2,T]$.

### B1 `B1_tw_speed_eulerian_vs_ale.png`: Eulerian vs ALE

**Caption.** Speed error $|c_h-c_\varepsilon|$ vs ε for the planar travelling wave, $w_0=0$ (left) and $w_0=1$ (middle), $\Delta t\in\{0.02,0.01,0.005\}$. Dashed squares: Eulerian scheme (S) (`simple.py`, same band meshes, exact Δt). Solid circles: ALE. **Right:** relative error vs $\Delta t/\varepsilon$, including $\Delta t=0.002$ and $0.001$ for the Eulerian scheme. The Eulerian error is a function of $\Delta t/\varepsilon$ alone. The ALE error does not depend on Δt (but see below for the default $\ell_0$ with $w_0=1$); its size is set by the lattice phase of the front (D3).

Signed errors $c_h-c_\varepsilon$:

| scheme | w0 | Δt | ε=0.01 | ε=0.005 | ε=0.002 | ε=0.001 |
|---|---|---|---|---|---|---|
| Eulerian (S) | 0 | 0.02 | -2.22e-01 | -3.00e-01 | -4.00e-01 | -4.62e-01 |
| Eulerian (S) | 0 | 0.01 | -1.49e-01 | -2.21e-01 | -3.25e-01 | -4.00e-01 |
| Eulerian (S) | 0 | 0.005 | -9.17e-02 | -1.49e-01 | -2.46e-01 | -3.25e-01 |
| Eulerian (S) | 1 | 0.02 | -9.57e-01 | -1.13e+00 | -1.32e+00 | -1.43e+00 |
| Eulerian (S) | 1 | 0.01 | -7.68e-01 | -9.60e-01 | -1.18e+00 | -1.32e+00 |
| Eulerian (S) | 1 | 0.005 | -5.72e-01 | -7.71e-01 | -1.02e+00 | -1.18e+00 |
| ALE | 0 | 0.02 | +1.66e-05 | +2.85e-05 | +3.52e-05 | +3.74e-05 |
| ALE | 0 | 0.01 | +1.76e-05 | +2.89e-05 | +3.50e-05 | +3.73e-05 |
| ALE | 0 | 0.005 | +1.77e-05 | +2.89e-05 | +3.50e-05 | +3.73e-05 |
| ALE (ℓ0=4) | 1 | 0.02 | -7.30e-04 | -6.31e-04 | -5.71e-04 | -5.51e-04 |
| ALE (ℓ0=4) | 1 | 0.01 | -7.28e-04 | -6.30e-04 | -5.71e-04 | -5.51e-04 |
| ALE (ℓ0=4) | 1 | 0.005 | -7.28e-04 | -6.28e-04 | -5.69e-04 | -5.50e-04 |
| ALE, default ℓ0=2 | 1 | 0.02 | -5.83e-04 | -4.33e-04 | -3.37e-04 | -3.03e-04 |
| ALE, default ℓ0=2 | 1 | 0.01 | -6.53e-04 | -5.29e-04 | -4.52e-04 | -4.25e-04 |
| ALE, default ℓ0=2 | 1 | 0.005 | -6.88e-04 | -5.76e-04 | -5.08e-04 | -4.86e-04 |

Before (old code): ALE $w_0=0$ +1.42e-4 … +1.83e-4; ALE ℓ0=4 $w_0=1$ −2.3e-4 … +1.3e-4 (scattered); default ℓ0 $w_0=1$ +1.9e-4 … +1.23e-3; Eulerian e.g. −0.144 on the diagonal Δt = ε (now −0.149).

- The Eulerian errors are identical along diagonals with equal $\Delta t/\varepsilon$: −0.149 for $(\Delta t,\varepsilon)=(0.01,0.01),(0.005,0.005),(0.002,0.002),(0.001,0.001)$; actual step = nominal step (`dts`).
- ALE, $w_0=0$: $2$–$4\cdot10^{-5}$, flat in Δt. The variation with ε comes from the lattice phase in which the initial front happens to fall (D3), not from ε itself.
- ALE, $w_0=1$, ℓ0 = 4: $-5.5$ to $-7.3\cdot10^{-4}$, flat in Δt to 3 digits; again a phase effect, larger for $w_0=1$ (D3: amplitude ±9e-4 at $h_f=\varepsilon/16$).
- ALE, $w_0=1$, default ℓ0 = 2: the same pinning error plus a remaining $O(\Delta t)$ term from the non-rigid zone of about $+0.007\,\Delta t$ (+1.05e-4 between Δt = 0.005 and 0.02 at ε = 0.01; +1.8e-4 at ε = 0.001), i.e. it still grows mildly as ε decreases. Before the rigid-zone change this term was about $0.06\,\Delta t$. See D1.

### B2 `B2_ale_mesh_sizes.png`: ALE error vs $h_f$ and $h_c$

**Caption.** ALE speed error for the travelling wave ($w_0=0$, $\Delta t=0.01$, B1 initial datum). Left: vs $h_f/\varepsilon$ at $h_c=0.0025$. Right: vs $h_c$ at $h_f=\varepsilon/32$.

| pts_per_eps ($h_f/\varepsilon$) | ε=0.01 | ε=0.005 | ε=0.002 | ε=0.001 |
|---|---|---|---|---|
| 4 (1/4) | -2.10e-03 | -2.03e-03 | -1.99e-03 | -1.97e-03 |
| 8 (1/8) | -2.83e-04 | -2.49e-04 | -2.30e-04 | -2.23e-04 |
| 16 (1/16) | +2.12e-05 | +3.28e-05 | +3.92e-05 | +4.13e-05 |
| 32 (1/32) | +2.34e-05 | +2.35e-05 | +2.32e-05 | +2.30e-05 |
| 64 (1/64) | -1.29e-05 | -8.08e-06 | -5.06e-06 | -4.11e-06 |

| $h_c$ | ε=0.01 | ε=0.001 |
|---|---|---|
| 0.04 | -3.09e-05 | -3.69e-05 |
| 0.02 | +8.06e-06 | +5.68e-06 |
| 0.01 | +1.97e-05 | +1.87e-05 |
| 0.005 | +2.26e-05 | +2.21e-05 |
| 0.0025 | +2.34e-05 | +2.30e-05 |

- The curves for different ε nearly coincide (ε enters through $h_f/\varepsilon$ and the lattice phase). From 1/4 to 1/8 the error drops by a factor 7–9; below $10^{-4}$ the values are dominated by the phase-dependent pinning error and are **not** a clean rate (sign changes). A clean measurement at fixed phase is D2: second order in $h_f/\varepsilon$.
- The $h_c$ dependence is at the level of $5\cdot10^{-5}$ for $h_c=0.04$ and below $2\cdot10^{-5}$ for $h_c\le0.02$; it does not support a rate.
- Before: "≈ first order in $h_f/\varepsilon$, floor 1.6e-4 from $h_f=\varepsilon/16$". That floor came from the two code defects (see F), not from the band.

### B3 `B3_rigid_mesh_dt.png`: rigid mesh at the exact speed

**Caption.** Whole mesh translating rigidly with $v=c_\varepsilon$ (`fixed_velocity=c_eps`, `ell0=1e6`, `remesh_distance=1e9`), $h_f=\varepsilon/16$, $h_c=0.01$. Front travel about 6. **Unchanged** (the fixes do not touch this mode; values identical to before).

| w0 | ε | Δt=0.1 | 0.05 | 0.02 | 0.01 | 0.005 |
|---|---|---|---|---|---|---|
| 0 | 0.005 | -1.8e-06 | +3.0e-07 | +1.1e-06 | +1.0e-06 | +9.9e-07 |
| 0 | 0.001 | -5.1e-07 | -1.3e-07 | +3.1e-08 | +3.1e-08 | +2.9e-08 |
| 1 | 0.005 | -2.7e-05 | -7.7e-06 | +5.1e-06 | +7.4e-06 | +2.2e-06 |
| 1 | 0.001 | -3.5e-06 | -4.9e-06 | -4.0e-07 | +1.1e-07 | +1.0e-07 |

Caveat (referee M8(a), confirmed by D4): with the mesh velocity prescribed, the front is dragged by the mesh. For $\Delta t\gg\varepsilon$ the front speed is mostly the mesh speed, so this test shows locking to the mesh, not accuracy. D4 quantifies the drag; D1 is the meaningful test.

### B4 `B4_radial_extinction_dt.png`: radial extinction, error vs Δt

**Caption.** ALE, exact self-similar extinction ($w_0=0.5$, $T=1$, $t\le0.75$, $L=1.5R+1$, $\ell_0=1$, remeshing distance 0.1, default mesh). Front error after removing the level-set shift $\varepsilon/\sqrt{2M}$, $\max_{t\ge0.05}|r_h(t)+\varepsilon/\sqrt{2M}-R\sqrt{T-t}|$, vs Δt. First order in Δt with a constant 1.4–1.9 that does not depend systematically on ε.

| N | ε | Δt=0.04* | 0.02* | 0.01 | 0.005 | 0.0025 | error/Δt at 0.0025 |
|---|---|---|---|---|---|---|---|
| 2 | 0.01 | 6.41e-02 | 3.07e-02 | 1.55e-02 | 7.98e-03 | 4.28e-03 | 1.71 |
| 2 | 0.002 | 6.19e-02 | 2.95e-02 | 1.47e-02 | 7.42e-03 | 3.81e-03 | 1.52 |
| 2 | 0.001 | 6.19e-02 | 2.96e-02 | 1.49e-02 | 7.59e-03 | 3.98e-03 | 1.59 |
| 3 | 0.01 | 7.08e-02 | 3.40e-02 | 1.72e-02 | 8.88e-03 | 4.78e-03 | 1.91 |
| 3 | 0.002 | 6.90e-02 | 3.31e-02 | 1.67e-02 | 8.56e-03 | 4.55e-03 | 1.82 |
| 3 | 0.001 | 6.80e-02 | 3.21e-02 | 1.57e-02 | 7.57e-03 | 3.54e-03 | 1.42 |

\* The actual steps are $0.75/19=0.0395$ and $0.75/38=0.0197$. Before: 1.35 ($N=2$) and 1.55 ($N=3$) at ε = 0.001, 1.84/1.99 at ε = 0.01. The larger rigid zone accounts for +2–4 % of the change (checked with `rigid_extra=0.5`: 3.84e-3 vs 3.98e-3 for N = 2, ε = 0.001, Δt = 0.0025); the rest comes from fixes 1–2. The statement "constant independent of ε" now holds only roughly (1.4–1.9).

---

## C. Limit problem and examples from the literature

### C1 `C1_ale_vs_limit.png`: $(P_\varepsilon)$ vs the computed limit problem

**Caption.** Non-self-similar radial data $u_0=a(1-r^2/R_0^2)^+$, $a=\sqrt{2M}R_0/2$, $R_0=2$, $w_0=0.5$, computed up to the time when the limit radius is $R_0/2$ ($t=0.565$ for $N=2$, $0.393$ for $N=3$). The limit problem is solved by front fixing (`limit.solve_radial`, $n=1600$). Left: fronts for $N=2$. Middle: front error after removing $\varepsilon/\sqrt{2M}$. Right: $\max_r|u_h^+-u|$ at the final time. Colour: ALE with $\Delta t=0.002$ (dotted), $0.0005$ (dashed) and $0.000125$ (solid). Grey: Eulerian BDF solver (`results/limit_comparison.json`, not rerun).

| N | ε | Δt | front err | shifted | sup | steps |
|---|---|---|---|---|---|---|
| 2 | 0.05 | 1.25e-4 | 5.67e-02 | 9.59e-03 | 6.94e-03 | 4520 |
| 2 | 0.01 | 1.25e-4 | 1.05e-02 | 1.12e-03 | 1.07e-03 | 4520 |
| 2 | 0.005 | 1.25e-4 | 5.19e-03 | 4.80e-04 | 4.80e-04 | 4520 |
| 2 | 0.002 | 1.25e-4 | 2.03e-03 | 1.41e-04 | 2.00e-04 | 4520 |
| 2 | 0.001 | 1.25e-4 | 1.02e-03 | 7.73e-05 | 1.13e-04 | 4520 |
| 2 | 0.001 | 5e-4 | 1.72e-03 | 7.78e-04 | 8.22e-04 | 1130 |
| 2 | 0.001 | 2e-3 | 4.74e-03 | 3.80e-03 | 4.01e-03 | 282 |
| 3 | 0.01 | 1.25e-4 | 1.06e-02 | 1.17e-03 | 1.11e-03 | 3144 |
| 3 | 0.001 | 1.25e-4 | 1.09e-03 | 1.44e-04 | 1.51e-04 | 3144 |
| BDF, N=2 | 0.01 | adaptive | 1.04e-02 | 9.59e-04 | 9.02e-04 | 1024 |
| BDF, N=2 | 0.001 | adaptive | 1.01e-03 | 6.74e-05 | 9.27e-05 | 14408 |

- For $\varepsilon\ge0.01$ the ε part dominates and ALE agrees with the BDF solver. For small ε the Δt part dominates at the larger steps.
- At ε = 0.001, $N=2$ the shifted error is 3.80e-3, 7.78e-4, 7.7e-5 for Δt = 2e-3, 5e-4, 1.25e-4: faster than first order between the last two steps, so "aε + bΔt with b ≈ 2.3–2.6" no longer describes the data well (b = 1.9, 1.6, 0.6 at the three steps).
- Before: floor 3.75e-4 ($N=2$) and 4.23e-4 ($N=3$) at ε = 0.001, Δt = 1.25e-4; now 7.7e-5 and 1.4e-4, comparable to BDF (6.7e-5).
- The unshifted front error is $\approx\varepsilon/\sqrt{2M}$ in all cases, i.e. the distance between the tracked level set $\{u=\varepsilon\}$ and the free boundary.

### C2 `C2_self_similar_convergence.png`: approach to the self-similar extinction profile

**Caption.** Same data as C1, followed until extinction. Top row: $s(t)/\sqrt{T-t}$ vs $T-t$; black: limit problem with $T$ fitted from $s^2\sim R^2(T-t)$ on $s<0.3$: $T=0.725090$ ($N=2$), $0.505756$ ($N=3$). Dashed: ALE, ε = 0.01 (Δt = 1e-4) and 0.002 (Δt = 2.5e-5), using $r_h+\varepsilon/\sqrt{2M}$ with $T_\varepsilon$ fitted in the same way. Dotted: ε = 0.002, BDF. Domain $L=R_0+0.5$. The ratio tends to the Galaktionov–Hulshof–Vázquez constant $R$ (2.513927 for $N=2$, 3.003951 for $N=3$). Bottom row: rescaled profiles; they converge to $\kappa\phi(\eta)$.

Limit problem, $s/\sqrt{T-t}$ (unchanged):

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
| N=2 | 0.725090 | 0.725024 | 0.72464 (before 0.72459) | 0.72533 (before 0.72525) |
| N=3 | 0.505756 | 0.505712 | 0.50543 (before 0.50540) | 0.50594 (before 0.50589) |

The ALE curve for ε = 0.002 dips by 1.4 % ($N=2$) and 1.6 % ($N=3$) below $R$ at $T-t\approx10^{-2}$ (minimum of $s/(R\sqrt{T-t})$ over $0.003<T-t<0.05$; ε = 0.01: 0.3 %, 0.5 %). This is the mesh-velocity taper near $r=0$ (issue 2). The BDF curve lies on the limit curve.

### C3 `C3_quenching.png`: quenching of small burnt kernels

**Caption.** Expanding flame from a burnt kernel of radius $r_0$: $u_0=\tanh(\sqrt{2M}(r-r_0))^+$, $A=1$, $w_0=0.5$, $N=2$, $c_0=1.0607$. Left: front radius $r_f(t)$ for several $r_0$ (ε = 0.01, ALE Δt = 0.001). Kernels that ignite (orange) first shrink, down to $r\approx0.53$ for $r_0$ just above critical (0.52 at ε = 0.0025), then expand with speed tending to $c_0$. Smaller kernels collapse (blue). Right: critical radius $r_*$ by bisection, run to $t=6$ (to 12 when undecided; tolerance $10^{-3}$). **Unchanged by the fixes** (identical brackets).

| ε | Δt | $r_*\in$ |
|---|---|---|
| 0.01 | 0.001 | [0.7148, 0.7156] |
| 0.005 | 0.001 | [0.7148, 0.7156] |
| 0.0025 | 0.001 | [0.7148, 0.7156] |
| 0.00125 | 0.001 | [0.7141, 0.7148] |
| 0.00125 | 0.00025 | [0.7156, 0.7164] |
| 0.01 | BDF (Eulerian) | [0.7148, 0.7156] |

$r_*=0.715\pm0.001$ for this profile and these Δt. The heuristic $1/c_0=0.943$ overestimates it by 30 %. The file `C3_quenching_dt0.005.json` (Δt = 0.005: 0.7129, 0.7129, 0.7105, 0.7059) **predates the fixes and was not rerun**.

---

## D. The error bounded by the theorems: travelling wave from $U^0=\Phi_\varepsilon$ (`tw_error.py`)

Set-up: $U^0_j=\Phi_\varepsilon(x_j-x_0)$, $w=\varepsilon w_0$, $x_0=(k+\theta)h_f$ with $kh_f\approx2$ and lattice phase θ = 0.37 (except D3). $L=x_0+c_\varepsilon T+12/c_\varepsilon+2$, $T=8$ ($w_0=0$) / 4 ($w_0=1$). Recorded at every step:

- $E=\max_{t_n\le T}\max_j|U^n_j-\Phi_\varepsilon(x^n_j-x_0-c_\varepsilon t_n)|$;
- the front position error $s_n=x_f^n-x_0-c_\varepsilon t_n$ (exact front at $\Phi_\varepsilon=\varepsilon$); speed error = least-squares slope of $s$ on $[T/2,T]$ (agrees with $s(T)/T$ and with the mean of $v_n-c_\varepsilon$);
- $E_s$: $E$ with the shift $s_n$ removed, i.e. $\max|U^n_j-\Phi_\varepsilon(x^n_j-x_0-c_\varepsilon t_n-s_n)|$;
- $\mu=\max_n|v_n-c_\varepsilon|$, $\mu_1$ = the same for $t_n\ge1$.

$E\approx\max|s|\cdot\max\Phi'_\varepsilon\approx c_\varepsilon\max|s|$ grows linearly in time with the speed error; $E_s$ is bounded.

### D1 `D1_ale_error_vs_dt.png`: ALE, error vs Δt

Rows $w_0=0,1$; columns $E$, $E_s$, $\mu_1$; curves ε = 0.01, 0.003, 0.001 (default $\ell_0$); for $w_0=1$ dashed: ℓ0 = 4, dotted: old rigid zone ($R_1=\text{band}+\ell_0/2$).

| w0 | ε | quantity | Δt=0.1 | 0.05 | 0.02 | 0.01 | 0.005 | 0.0025 |
|---|---|---|---|---|---|---|---|---|
| 0 | 0.01 | E | 4.7e-4 | 4.5e-4 | 4.5e-4 | 4.4e-4 | 4.4e-4 | 4.4e-4 |
| 0 | 0.001 | E | 4.8e-4 | 4.6e-4 | 4.5e-4 | 4.5e-4 | 4.4e-4 | 4.4e-4 |
| 0 | all | $E_s$ | 6.1–7.3e-5 | 4.1e-5 | 4.3e-5 | 4.4e-5 | 4.4e-5 | 4.4e-5 |
| 0 | all | $\mu_1$ | 1.17e-4 | 1.15e-4 | 1.13e-4 | 1.13e-4 | 1.12e-4 | 1.12e-4 |
| 0 | all | speed err | +9.0–9.4e-5 | +8.5e-5 | +8.1e-5 | +7.9e-5 | +7.9e-5 | +7.8e-5 |
| 0 | all | $\mu$ (all steps) | 2.7e-4 | 3.7e-4 | 5.7e-4 | 7.8–8.1e-4 | 1.1e-3 | 1.5–1.6e-3 |
| 1 | 0.01 | E | 5.3e-3 | 3.8e-3 | 3.1e-3 | 2.8e-3 | 2.7e-3 | 2.6e-3 |
| 1 | 0.001 | E | 7.4e-3 | 4.8e-3 | 3.4e-3 | 2.9e-3 | 2.7e-3 | 2.6e-3 |
| 1 | 0.01 | $E_s$ | 4.3e-4 | 2.0e-4 | 8.6e-5 | 8.7e-5 | 8.9e-5 | 9.0e-5 |
| 1 | 0.001 | $E_s$ | 7.1e-4 | 3.3e-4 | 1.2e-4 | 7.9e-5 | 8.0e-5 | 8.2e-5 |
| 1 | 0.01 | $\mu_1$ | 1.1e-3 | 7.1e-4 | 5.1e-4 | 4.3e-4 | 4.0e-4 | 3.8e-4 |
| 1 | 0.001 | $\mu_1$ | 1.6e-3 | 9.5e-4 | 6.0e-4 | 4.9e-4 | 4.4e-4 | 4.2e-4 |
| 1, ℓ0=4 | all | E / $E_s$ / $\mu_1$ | 2.4–2.5e-3 / 8.4–9.2e-5 / 3.6–3.9e-4 | | flat | | | |
| 1, old $R_1$ | 0.001 | E / $\mu_1$ | 3.4e-2 / 7.1e-3 | 1.8e-2 / 3.6e-3 | 8.4e-3 / 1.6e-3 | 5.4e-3 / 1.0e-3 | 4.0e-3 / 7.0e-4 | 3.2e-3 / 5.5e-4 |

- $w_0=0$: **no Δt dependence and no ε dependence** of $E$, $E_s$, $\mu_1$ from Δt = 0.0025 to 0.1 (Δt/ε up to 100, Δt/ε² up to $10^5$). $E$ is dominated by the drift of the front at the pinned speed $c_\varepsilon+7.9\cdot10^{-5}$ (D3; $\theta=0.37$, $h_f=\varepsilon/16$).
- $\mu$ over all steps is set by the first steps (the discrete profile relaxing from $\Phi_\varepsilon$ to the discrete wave) and *grows* as Δt decreases ($\approx\Delta t^{-0.45}$; the displacement in the first steps is $O(10^{-5})$). So the realised μ of the theorem is $\approx10^{-3}$ including the start, $\approx10^{-4}$ ($w_0=0$) after $t=1$.
- $w_0=1$, default $\ell_0=2$: a remaining Δt term from the non-rigid zone, $\mu_1-\mu_1(\ell_0{=}4)\approx0.007\,\Delta t$ (ε = 0.01) to $0.012\,\Delta t$ (ε = 0.001). With the old rigid zone it was $0.04$–$0.07\,\Delta t$ (factor ≈ 5.5 larger, consistent with $e^{-c\ell_0/2}/e^{-c\ell_0}=e^{1.63}=5.1$). With $\ell_0=4$ (or the rigid mesh) the Δt term is gone (< 1 % change from Δt = 0.0025 to 0.1). **So the larger rigid zone reduces but does not remove the default-ℓ0 problem for $w_0=1$**; it still grows mildly as ε decreases. A default $\ell_0\approx6.5/c$ or a rigid zone of a few $1/c$ would remove it.
- Adopted: `rigid_extra=1.0` (radial B4 changes by +2–4 %, C3 unchanged).

### D2 `D2_ale_error_vs_mesh.png`: ALE, error vs $h_f/\varepsilon$ and $h_c$

Δt = 0.01, θ = 0.37 fixed. Left/middle: vs $h_f/\varepsilon\in\{1/2,\dots,1/32\}$ at $h_c=0.0025$; right: vs $h_c$ at $h_f=\varepsilon/32$.

| $h_f/\varepsilon$ | 1/2 | 1/4 | 1/8 | 1/16 | 1/32 |
|---|---|---|---|---|---|
| $w_0=0$: $E_s$ | 3.3e-3 | 7.7e-4 | 1.9e-4 | 4.5e-5 | 1.1e-5 |
| $w_0=0$: speed err | +5.9e-3 | +1.37e-3 | +3.3e-4 | +8.3e-5 | +2.2e-5 |
| $w_0=1$: $E_s$ | 7.3e-3 | 1.6e-3 | 3.8e-4 | 8.9e-5 | 3.8–6.3e-5 |
| $w_0=1$: speed err | +3.3e-2 | +7.4e-3 | +1.8e-3 | +4.8–5.3e-4 | +1.6–2.1e-4 |

(values for ε = 0.01, 0.003, 0.001 agree to 1–2 digits; ranges given where they differ.)

- **Second order in $h_f/\varepsilon$** (ratio 4.0–4.4 per halving) for $E_s$ and the speed error, independent of ε. No $h_f/\varepsilon^2$ dependence is visible (that would differ by 10× between ε = 0.01 and 0.001).
- $h_c$ ($w_0=0$): $E_s=0.7$–$1.4\cdot10^{-5}$ for all $h_c\in[0.0025,0.04]$, speed error ≤ 4e-5. $w_0=1$: the speed error moves from +1.6e-4 ($h_c=0.0025$) to −5.5e-4 … −7.6e-4 ($h_c=0.04$), differences ≈ 7e-4, 2.3e-4, 6e-5, 1.3e-5 for $h_c$ = 0.04, 0.02, 0.01, 0.005: roughly $O(h_c^2)$; $E_s$ 3.8e-5 → 8.7e-5 (ε = 0.01).

### D3 `D3_phase_pinning.png`: the "floor" is lattice pinning

ALE, ε = 0.01 (lines) and 0.001 (crosses), Δt = 0.02; initial front at phase θ ∈ [0, 1) of a band cell; mean of $v_n-c_\varepsilon$ over $[T/2,T]$.

| $h_f/\varepsilon$ | $w_0=0$: min / max / mean over θ | $w_0=1$: min / max / mean |
|---|---|---|
| 1/8 | −7.5e-4 / +4.5e-4 / −1e-6 | −3.6e-3 / +2.2e-3 / +1e-5 |
| 1/16 | −1.9e-4 / +1.1e-4 / −3e-6 | −8.9e-4 / +5.3e-4 / −3e-6 |
| 1/32 | −5.0e-5 / +2.4e-5 / −4e-6 | −2.3e-4 / +1.3e-4 / −6e-6 |

- With implicit velocity the front stays at its initial phase θ (the reference point moves with the band and the band is now transferred exactly). The realised speed $c_\varepsilon+g(\theta)$: $g$ is smooth in θ, vanishes near θ ≈ 0.2 and 0.8, has **zero mean**, and its amplitude is $\propto(h_f/\varepsilon)^2$ (factor 4 per halving) and independent of ε (crosses on the lines). The same values are obtained with a fully rigid mesh with implicit velocity, with Δt = 0.002 and 2e-4, with $h_c=0.0025$, with a longer band ($K_b=40,K_f=20$) and with growth 1.02 (scratch tests, ε = 0.01, θ = 0.37: all +8.6e-5 … +1.0e-4).
- A rigid mesh at prescribed $v=c_\varepsilon$ (B3) lets the front slide through the lattice and averages $g$: hence the tiny B3 errors.
- **This is the "floor" of the old tables.** Old setup, $w_0=0$, ε = 0.01, Δt = 0.02 (table F): pre-anchor-fix code +1.66e-4 (cycle) / +1.06e-4 (lsq) / +4.4e-5 (end points); after the anchor fix +5.0e-5 with all three; current code +1.6e-5 with all three. The value depends on where the reference point lands in the lattice.

### D4 `D4_drag.png`: drag of the front by the mesh

Rigid mesh with prescribed $v=c_\varepsilon+\delta v$, $\delta v=10^{-3}$, from $\Phi_\varepsilon$; $c_h$ from the slope of $s$ on $[T/2,T]$ ($T=6$, 3). Ratio $(c_h-c_\varepsilon)/\delta v$:

| Δt/ε | 0.03 | 0.1 | 0.3 | 1 | 3 | 10 | 30 | 100 |
|---|---|---|---|---|---|---|---|---|
| $w_0=0$ (ε = 0.01 / 0.003 / 0.001) | 0.04 | 0.07 / 0.08 / 0.09 | 0.14 / 0.15 / 0.15 | 0.33 | 0.60 | 0.83–0.85 | 0.94 | 0.98 |
| $w_0=1$ | 0.21 | 0.34 / 0.34 / 0.37 | 0.46–0.48 | 0.67–0.68 | 0.82–0.85 | 0.95 | 0.98 | 0.99 |

- The ratio is a function of Δt/ε alone (ε-independent to ±0.02 for $w_0=0$): this is the μ-term of the ALE theorem; the front-speed error is $\alpha(\Delta t/\varepsilon)\,\mu$ with α → 1 for Δt ≫ ε (referee: 0.03…0.91; reproduced).
- For $w_0=1$ at small Δt/ε the ratio is presumably biased upwards by lattice pinning (pinning amplitude ±9e-4 at $h_f=\varepsilon/16$ is comparable to δv, so the front does not slide freely); these values are not conclusive.
- Consequence for the implicit velocity: $\partial(\text{gap})/\partial v=-(1-\alpha)\Delta t$, nearly degenerate for Δt ≫ ε; a phase error $g(\theta)$ is then not corrected (D3).

### D5 `D5_eulerian_error.png`: Eulerian scheme (S), error vs Δt/ε and $h/\varepsilon$

`simple.py` with exact nominal Δt, from $\Phi_\varepsilon$, $T=4$ ($w_0=0$) / 2.5 ($w_0=1$), monitor every max(1, 0.005/Δt) steps.

| Δt/ε | 0.01 | 0.03 | 0.1 | 0.3 | 1 |
|---|---|---|---|---|---|
| $w_0=0$: $\lvert c_h-c_\varepsilon\rvert/c_\varepsilon$ | 4.2–4.8e-3 | 1.4e-2 | 4.5e-2 | 0.122 | 0.30 |
| $w_0=0$: $E_s$ | 1.4e-3 | 4.0e-3 | 1.3e-2 | 3.4e-2 | 8.3e-2 |
| $w_0=1$: relative speed error | 1.4–1.6e-2 | 4.4–4.5e-2 | 0.128 | 0.277 | 0.50 |
| $w_0=1$: $E_s$ | 5.8e-3 | 1.7e-2 | 4.9e-2 | 0.11 | 0.22 |

(ranges over ε = 0.01, 0.003, 0.001.) The scheme is slower than $c_\varepsilon$.

- The error is a function of **Δt/ε alone** and **linear** for Δt/ε ≲ 0.1: relative speed error ≈ 0.46 Δt/ε ($w_0=0$), 1.5 Δt/ε ($w_0=1$); $E_s\approx0.13\,\Delta t/\varepsilon$, $0.58\,\Delta t/\varepsilon$. This is the conjectured sharp rate Δt/ε (not Δt/ε²).
- $h_f/\varepsilon$ (right panel, Δt = 0.01 ε): the Δt term dominates; going from $h_f=\varepsilon/32$ to $\varepsilon/8$ changes the speed by < 1 % and $E_s$ by < 2 %; at $\varepsilon/4$ the speed error increases by 6–11 % and $E_s$ by 7–12 %. **The $h/\varepsilon$ part cannot be resolved** at this Δt; consistent with the second-order behaviour of D2.

---

## F. The old "floor 1.6e-4" with three estimators

Old B1 set-up ($w_0=0$, datum $A(1-e^{-c_\varepsilon(x-2)})^+$, $T=12$, default mesh). Speed error by: cycle speed / least squares on $[T/2,T]$ / end points on $[T/2,T]$. Old code from git worktrees (`results/paper/F_floor_old_code.json`); current code `tw_error.py` group F.

| code | ε | Δt | cycle | lsq | end points |
|---|---|---|---|---|---|
| before anchor fix (d5cd9b8) | 0.01 | 0.02 | +1.66e-4 | +1.06e-4 | +4.4e-5 |
| | 0.01 | 0.005 | +1.42e-4 | +7.9e-5 | +2.0e-5 |
| | 0.001 | 0.02 | +1.83e-4 | +1.22e-4 | +6.9e-5 |
| | 0.001 | 0.005 | +1.55e-4 | +9.1e-5 | +3.8e-5 |
| anchor fixed (46cef8c) | 0.01 | 0.02 | +4.99e-5 | +4.97e-5 | +4.91e-5 |
| | 0.01 | 0.005 | +2.62e-5 | +2.62e-5 | +2.61e-5 |
| | 0.001 | 0.02 | +7.80e-5 | +7.77e-5 | +7.71e-5 |
| | 0.001 | 0.005 | +4.66e-5 | +4.65e-5 | +4.64e-5 |
| current | 0.01 | 0.02 | +1.66e-5 | +1.65e-5 | +1.62e-5 |
| | 0.01 | 0.005 | +1.77e-5 | +1.77e-5 | +1.77e-5 |
| | 0.001 | 0.02 | +3.74e-5 | +3.73e-5 | +3.70e-5 |
| | 0.001 | 0.005 | +3.73e-5 | +3.72e-5 | +3.72e-5 |
| current, rigid mesh + implicit v (no remeshing, no taper) | 0.01 | 0.02 / 0.005 | – | +1.15e-5 / +1.64e-5 | +1.12e-5 / +1.63e-5 |
| | 0.001 | 0.02 / 0.005 | – | +3.05e-5 / +3.53e-5 | +3.02e-5 / +3.53e-5 |

Cause of the old floor: (i) the band anchor (remeshing spikes; the cycle estimator, which starts and ends at remeshing steps, picked them up), (ii) the reference point on the auxiliary grid (sets the pinning phase; with the anchor offset of 0.35 $h_f$ per remeshing the phase also changed at every remeshing), (iii) the residual is the phase-dependent pinning error $g(\theta)$ of D3. Remeshing and the non-rigid zone now contribute < 1e-5 for $w_0=0$ (current vs rigid rows); the Δt dependence for ε = 0.01 (1.15e-5 vs 1.64e-5) is the transient from the non-ε-wave datum.

---

## Table `tab:simple` (`tab_simple.py`)

Eulerian scheme (S), $w_0=0$, Δt ∝ ε², exact nominal Δt. Datum and estimator as the old table (B1 datum, $T=12$, least squares on $[T/2,T]$; end points agree to 3 digits). Column "Φ_ε": datum $\Phi_\varepsilon$, $T=6$, slope of $s$.

| ε | Δt | Δt/ε | $c_{(S)}-c_\varepsilon$ (new) | before | from $\Phi_\varepsilon$ |
|---|---|---|---|---|---|
| 0.02 | 4ε² = 1.6e-3 | 0.08 | −1.863e-2 | −1.85e-2 | −2.00e-2 |
| 0.02 | ε² = 4e-4 | 0.02 | −4.85e-3 | −4.60e-3 | −5.21e-3 |
| 0.02 | ε²/4 = 1e-4 | 0.005 | −1.23e-3 | −9.30e-4 | −1.32e-3 |
| 0.01 | 4ε² = 4e-4 | 0.04 | −9.49e-3 | −9.41e-3 | −1.02e-2 |
| 0.01 | ε² = 1e-4 | 0.01 | −2.39e-3 | −2.21e-3 | −2.57e-3 |
| 0.01 | ε²/4 = 2.5e-5 | 0.0025 | −5.67e-4 | −3.69e-4 | −6.13e-4 |

The old values were smaller in magnitude because the actual step was smaller than nominal. The new values are linear in Δt at fixed ε (ratio 3.8–4.2 per factor 4) and all six satisfy $c_{(S)}-c_\varepsilon\approx-(0.39\text{–}0.42)\,c_\varepsilon\Delta t/\varepsilon$, i.e. they are a function of Δt/ε, not of Δt/ε² (B1 datum; −0.44 from $\Phi_\varepsilon$; D5 gives 0.46 on $T=4$).

---

## E. Transients, accumulated mesh-velocity error, phases (`transient.py`), and A5 (`mesh_snapshots.py`)

Added in the second revision (November 2026), in response to the second referee report.

- **A5 `A5_mesh_snapshots.png`** (`mesh_snapshots.py`, 10 s): the solution ($u$, $Y$) and the mesh (one tick per node) at four times, ALE, ε = 0.02, $h_f=\varepsilon/16$, $h_c=0.05$. Left: planar travelling wave ($w_0=0$, 1464 nodes at all times). Right: radial extinction of the self-similar ball ($N=2$, $w_0=0.5$, 1012 nodes). Insets: the reaction layer in units of ε. This figure and A1 are the mesh figures of the paper; A2–A4 are no longer in the paper.
- **E1** (figure `D6_transients.png`, left): Eulerian scheme (S) from $\Phi_\varepsilon$, $w_0=0$, ε ∈ {0.02, 0.01}, Δt ∈ {4ε², ε², ε²/4}, $T=60$. Speed error by least squares on successive windows of length 6, divided by Δt/ε. The speed relaxes slowly (time scale $4/c^2\approx12$) and is settled to 0.1 % on the last windows. On the window [54, 60] the speed error divided by Δt/ε is −0.2161, −0.2285, −0.2269 for ε = 0.01 and Δt/ε = 0.0025, 0.01, 0.04, and −0.2348, −0.2322, −0.2225 for ε = 0.02 and Δt/ε = 0.005, 0.02, 0.08; on the first window [0, 6] it is between −0.283 and −0.304. The averaged truncation error (`truncation_check.py`) gives −0.2345 and −0.2351. The values of `tab_simple.py` ($T=12$, window [6, 12]) and of D5 ($T=4$, window [2, 4]) are transient values (second referee report, N5).
- **E2** (figure `D6_transients.png`, right): ALE from $\Phi_\varepsilon$, ε = 0.01, Δt = 0.02, $T=24$: block means of $v_n-c_\varepsilon$ for the lattice phases ϑ = 1/32, 0.37, 17/32, 25/32.
- **E3**: the quantities of the ALE theorem, ALE from $\Phi_\varepsilon$, $T=8$, ϑ = 0.37:

| ε | Δt | $M_T=\Delta t\sum\lvert v_n-c_\varepsilon\rvert$ | $\max_n\lvert v_n-c_\varepsilon\rvert$ | $E$ | $E_s$ |
|---|---|---|---|---|---|
| 0.01 | 0.1 | 8.16e-4 | 2.74e-4 | 4.68e-4 | 6.10e-5 |
| 0.01 | 0.02 | 7.78e-4 | 5.65e-4 | 4.46e-4 | 4.29e-5 |
| 0.01 | 0.0025 | 7.70e-4 | 1.45e-3 | 4.42e-4 | 4.44e-5 |
| 0.001 | 0.1 | 8.31e-4 | 2.76e-4 | 4.80e-4 | 7.31e-5 |
| 0.001 | 0.02 | 7.80e-4 | 5.78e-4 | 4.50e-4 | 4.26e-5 |
| 0.001 | 0.0025 | 7.67e-4 | 1.59e-3 | 4.43e-4 | 4.45e-5 |

  $M_T$ does not depend on Δt or ε; the maximum of $\lvert v_n-c_\varepsilon\rvert$ is attained in the first steps and grows as Δt decreases.
- **E4**: over 16 lattice phases (ε = 0.01, Δt = 0.02, $T=8$): $E$ from 1.3e-5 to 9.5e-4, $E_s$ from 1.3e-5 to 9.6e-5, speed error from −1.7e-4 to +1.0e-4, $M_T$ from 2.5e-5 to 1.7e-3. For ε = 0.001 (8 phases): $E$ from 1.1e-4 to 7.7e-4, speed error from −1.3e-4 to +9.6e-5.

## G and the third revision (after the third referee report)

- **G** (`in_hypotheses.py`, 1 min): uniform mesh with $h=\varepsilon^2/8$ and $\varepsilon^2/16$, ε ∈ {0.1, 0.05}, translating rigidly with a prescribed velocity, from $\Phi_\varepsilon$, $T=3$ (the setting of the moving-mesh theorem). With $v=c_\varepsilon$: $E\le7.6\cdot10^{-6}$ for Δt = 0.1 and 0.01 ($\rho_hT$ = 0.19–0.38). With $v=c_\varepsilon+0.01$: $E$ = 6.6e-3 / 1.0e-3 (ε = 0.1, Δt = 0.1 / 0.01) and 9.8e-3 / 2.0e-3 (ε = 0.05), below $K\mu T\approx1.7\cdot10^{-2}$.

- **H** (`convergence_thm.py`, 1 min, `H_convergence_thm.json`): the setting of the convergence theorem (Theorem 5.1 of the paper). Radial self-similar solution $v$ of the limit problem ($N=2$, $w_0=0.5$, extinction time 1), data (eq:initdata), $0\le t\le 0.5$, no remeshing, $h_f=\varepsilon/16$, ε ∈ {0.02, 0.01, 0.005}, Δt/ε ∈ {1/64, 1/16, 1/4, 1, 4}. $E_v=\max_{n,j}|U^n_j-v(x^n_j,t_n)|$.
  - *Fixed uniform mesh* ($v_n=0$): $E_v$ = 0.022–0.023, 0.076–0.077, 0.177, 0.264, 0.296–0.302 for the five values of Δt/ε, the same for the three ε to two digits (a function of Δt/ε only); the front $\{u=\varepsilon\}$ at $t=0.5$ is behind by 0.002–0.018, 0.06–0.08, 0.21–0.23, 0.43–0.45, 0.57–0.62.
  - *Band mesh with prescribed velocity* $v_n=\dot R(t_n)$ (ALE with `fixed_velocity`, `remesh=False`, `ell0=0.5`, $L=5$; the cells outside the rigid zone stretch up to 0.044): $E_v=0.5\,\varepsilon=\varepsilon w_0$ for Δt ≤ ε (the value of $v-u^\varepsilon$ in the burnt region), and 0.64–0.67 ε for Δt = 4ε; the front is at −1.01…−1.07 ε from $R(0.5)$ for Δt ≤ ε (the shift $-\varepsilon/\sqrt{2M}=-0.94\,\varepsilon$ plus −0.07…−0.13 ε), and −0.33…−0.37 ε for Δt = 4ε.
  - Since this revision the figures A1 (node trajectories) and C3 (quenching) and the cost table A3 are no longer in the paper; their scripts and outputs are kept here as extra material.
  - New options of `simulate_ale`: `fixed_velocity` may be a function of $t$ (evaluated at $t_n$), and `remesh=False` disables the rebuilding of the mesh.
- **D3 recomputed on a late window.** The mesh velocity relaxes on the time scale $4/c^2$; D3 now uses $T=32$ ($w_0=0$) and $T=8$ ($w_0=1$, with $\ell_0=4$) and the mean of $v_n-c_\varepsilon$ on the last eighth of the run. Settled values (ε = 0.01, 16 phases), min / max / mean:

| $h_f/\varepsilon$ | $w_0=0$ | $w_0=1$ |
|---|---|---|
| 1/8 | −6.08e-4 / +3.84e-4 / +1.2e-5 | −3.52e-3 / +2.03e-3 / −5.1e-5 |
| 1/16 | −1.45e-4 / +1.02e-4 / +9.9e-6 | −9.24e-4 / +4.54e-4 / −6.3e-5 |
| 1/32 | −2.90e-5 / +3.26e-5 / +9.5e-6 | −2.81e-4 / +6.27e-5 / −6.6e-5 |

  The oscillating part scales exactly like $(h_f/\varepsilon)^2$ (half-ranges 4.96e-4, 1.24e-4, 3.08e-5 and 2.78e-3, 6.89e-4, 1.72e-4); the mean is an offset independent of $h_f$ (about +1e-5 and −6e-5). The values quoted in section D3 above (window $[T/2,T]$ with $T=4$ and 2) were transient values.
- **Eulerian residual.** The Δt-independent part of the Eulerian speed error (about +4.5e-5 at ε = 0.01) is a coarse-mesh contribution: with $h_c=0.005$ or with margins of 30ε the slope at ε = 0.01, Δt = 4ε², on [30, 36], changes from −9.089e-3 to −9.132e-3; at ε = 0.02 with $h_c=0.02$ the same shift appears.

## Issues in `combustion1d/`

1. **`simulate_ale` may stop when the front disappears** (kernel quenching, extinction). `common.simulate_ale_safe` keeps the partial result and records the remeshing steps. Not re-examined after the fixes.
2. **The mesh velocity is tapered to zero near $r=0$** (`chi` multiplied by `clip(z/ell_end)`, `ell_end = min(0.5, 0.1 L)`). For $r_f<$ `ell_end` the scheme is partly Eulerian and the $O(\Delta t\,v/\varepsilon)$ lag returns; $v\sim R^2/(2r)$ blows up near extinction. Effect: C2 dip (1.4–1.6 % at ε = 0.002). Unchanged.
3. **Rigid zone / default $\ell_0$ for fast fronts.** Now $R_1=\text{band}+\ell_0$ (was $+\ell_0/2$). For $w_0=1$ ($c=1.63$, $\ell_0=2$) the remaining Δt term is ≈ 0.007–0.012 Δt (was 0.04–0.07 Δt); $\ell_0=4$ removes it. Not fixed further.
4. **The rigid mode is only usable for moderate travel distances** (nodes in the tapered end zones collapse at $x=L$). B3, D4 limit the travel to about 5–6.
5. **Lattice pinning** (new, D3): with implicit velocity the speed error depends on the lattice phase of the front, amplitude $\propto(h_f/\varepsilon)^2$, zero mean over the phase (−1.9e-4…+1.1e-4 for $w_0=0$, −8.9e-4…+5.3e-4 for $w_0=1$ at $h_f=\varepsilon/16$). Speed estimators no longer matter (all agree).

## Runtimes (4 processes, `OMP_NUM_THREADS=1`)

A ~4 min (BDF runs 5.5, 12, 39, 122 s); B ~2 min; C ~12 min (C3 bisection ~7 min); D+F 6 min (largest: Eulerian ε = 0.001, Δt = 1e-5, 400k steps, 2 min); `tab_simple` ~3 min (ε = 0.01, Δt = 2.5e-5, 480k steps, 2.5 min).
