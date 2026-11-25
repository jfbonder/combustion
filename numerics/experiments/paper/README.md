# Numerical experiments of the paper

This directory holds the scripts that produce the figures, tables and in-line numbers
of the numerical section of the paper, and a few extra experiments. This file says what
each script computes, what it writes, and which item of the paper it corresponds to.
The package itself is documented in [`../../README.md`](../../README.md).

## Running

Everything is run from `numerics/`, with that directory on the Python path:

```bash
cd numerics
export PYTHONPATH=$PWD OMP_NUM_THREADS=1
python -m experiments.paper.convergence          # run and plot
python -m experiments.paper.convergence plot     # only redraw from the stored results
python -m experiments.paper.tw_error run D3 D4   # only some groups (tw_error, transient, limit_examples)
```

- The group scripts accept `run` or `plot`; with no argument they do both. The parallel
  scripts use 2 or 4 worker processes (`OMP_NUM_THREADS=1` avoids oversubscription).
- `common.py` holds the test problems, the I/O helpers, the mesh analysis and the plot
  helpers; `../util.py` holds the provenance record and the process pool.
- Results are written to `results/paper/*.json` (plus `A_mesh_*.npz`, the node histories of
  group A, which are not stored in git). Figures are written to `figures/paper/*.png`
  (200 dpi).
- **Provenance.** Every JSON file written by the scripts carries a record
  `provenance = {git, git_dirty, python, numpy, scipy, platform}` (`experiments.util.provenance`):
  the commit of the code (`git rev-parse HEAD`), whether the working tree of `numerics/`
  differed from it, and the library versions. A file that is a dictionary has it under the
  key `provenance`; in a file that is a list of records, every record has it (so that
  records kept from an earlier run keep their own). Files without the record were written
  before it was introduced.
- **Process pools.** The scripts create their pools with `experiments.util.worker_pool`,
  which installs the reaction rate of the parent process in every worker. (A rate set with
  `combustion1d.model.set_reaction` is state of one process and would otherwise be lost in
  workers started with the `spawn` method, the default on macOS and Windows.) All the
  experiments here use the default rate.
- **Failures.** A numerical failure of a run (violation of the Péclet condition, crossing
  mesh nodes, failure of the BDF integrator) is never used as a result. In the scripts that
  run many independent cases (`convergence`, `tw_error`) the failed case is stored as
  `{..., "failed": <message>}` and printed with the tag `FAILED`; in `limit_examples` the
  script stops. The disappearance of the front is not a failure: the solvers report it as
  `status = "extinct"`.

Model in all the experiments: $f(s)=(1-s)^+$, $M(w_0)=(1+w_0)^3/6$; $c_0=\sqrt{2M}/A$ is the
planar speed of the limit problem and $c_\varepsilon$ the exact speed of the travelling wave of
$(P_\varepsilon)$ (ODE shooting, `exact.tw_speed_eps`).

Default mesh, `MeshParams()`: a band of width $(K_b+K_f+2K_m)\varepsilon=55\varepsilon$ with
$h_f=\varepsilon/16$, a geometric transition with ratio 1.1, and coarse spacing $h_c=0.01$.

## Scripts, outputs and items of the paper

The labels in the first column are those of `paper/main.tex`. Times are the estimates
recorded when the results were generated (4 processes on one desktop machine); they are
**approximate** and machine dependent, and a laptop can be several times slower.

| item of the paper | script | solver | output | time (approx.) |
|---|---|---|---|---|
| Table `tab:speedcoef` (Proposition `prop:speed`) | `speed_coefficient.py` | ODE shooting | `speed_coefficient.json` | ≈ 10 s |
| eq. `an:eq:solv` (averaged truncation error) | `../truncation_check.py` | none (exact profile) | `truncation_check.json` | ≈ 15 s |
| Fig. `fig:A5` | `mesh_snapshots.py` | ALE | `A5_mesh_snapshots.png` | ≈ 10 s |
| Figs. `fig:D1`, `fig:D2`, `fig:D3`; drag ratios and Eulerian errors quoted in Section `sec:exp:tw` | `tw_error.py` (groups D1–D5, F) | ALE, Eulerian (S) | `D_tw_error.json`, `D1`–`D5_*.png` | ≈ 6–8 min |
| Fig. `fig:D6`; accumulated mesh-velocity error; errors over the lattice phase | `transient.py` (groups E1–E4) | ALE, Eulerian (S) | `E_transient.json`, `D6_transients.png` | ≈ 30 min (E1 ≈ 25 min) |
| Uniform mesh inside the hypotheses of Theorem `an:thm:ale` (Section `sec:exp:tw`) | `in_hypotheses.py` (group G) | ALE, rigid mesh | `G_in_hypotheses.json` | ≈ 1–2 min |
| Section `sec:exp:thm`: setting of Theorem `cv:thm:main` | `convergence_thm.py` (group H) | (S) on a fixed mesh; ALE with prescribed velocity | `H_convergence_thm.json` | ≈ 1–3 min |
| Fig. `fig:B4` | `convergence.py` (group B; B4 is the figure) | ALE, Eulerian (S) | `B_convergence.json`, `B1`–`B4_*.png` | ≈ 2 min |
| Figs. `fig:C1`, `fig:C2` | `limit_examples.py` (groups C1–C3) | ALE, BDF, limit solver | `C1_limit_comparison.json`, `C2_self_similar.json`, `C3_quenching.json`, `C1`–`C3_*.png` | ≈ 12 min (C3 ≈ 7 min) |
| Fig. `fig:instability` | `../traveling_wave.py`, then `../plots.py` | BDF | `results/traveling_wave.json`, `figures/tw_*.png` | ≈ 20 min |
| Table `tab:limit` | `../limit_comparison.py`, then `../plots.py` | BDF, limit solver | `results/limit_comparison.json` | ≈ 15 min |

Material that is produced by these scripts but is not in the paper:

| what | script | output |
|---|---|---|
| node trajectories, local mesh size, node counts and cost (A1–A4) | `meshes.py` (group A, ≈ 4 min) | `A_meshes.json`, `A_mesh_*.npz`, `A1`–`A4_*.png` |
| Eulerian vs ALE speed errors, mesh sizes, rigid mesh (B1–B3) | `convergence.py` | `B1`–`B3_*.png` |
| quenching of small burnt kernels (C3) | `limit_examples.py` | `C3_quenching.json`, `C3_quenching.png` |
| drag of the front by a rigid mesh (D4), Eulerian error (D5), three speed estimators (F) | `tw_error.py` | `D4_drag.png`, `D5_eulerian_error.png` |
| Eulerian scheme with $\Delta t\propto\varepsilon^2$ | `tab_simple.py` (≈ 3 min) | `tab_simple.json` |
| non-constant $w_0$ (group I) | `variable_w0.py` (≈ 10 min, see its docstring) | `I_variable_w0.json`, `I1_variable_w0.png` |
| work–precision comparison (group J) | `work_precision.py` (≈ 45 min, see its docstring) | `J_work_precision.json`, `J1_work_precision.png` |
| radial validation with the BDF solver | `../radial.py` (≈ 5 min) | `results/radial.json` |

`../../reproduce_paper.sh` runs the scripts of the first table, `meshes.py` and
`../radial.py`.

Which solver is behind each figure:

| figure / table | solver |
|---|---|
| A1 left, A2, A4, A5, A3 "ALE"; B1 "ALE", B2, B3, B4; C1–C3 coloured curves; D1–D4, D6 right, F, G, H "band" | scheme (S) on the moving mesh, `ale.simulate_ale` |
| B1 "Eulerian", D5, D6 left, `tab_simple` | scheme (S) on band meshes rebuilt around the front, `simple.simulate_simple`, fixed $\Delta t$ |
| H "fixed" | scheme (S) on a fixed uniform mesh (assembled in `convergence_thm.py` with `solver.laplacian`) |
| A1 right, A3 "Eulerian BDF"; C1 grey, C2 dotted, C3 "Eulerian BDF"; `fig:instability`, `tab:limit`, `radial.py` | adaptive BDF with PCHIP transfer, `solver.simulate` (not covered by the theorems) |

## How the runs are set up and measured

- **Travelling-wave runs, groups B and F.** Initial datum $A(1-e^{-c_\varepsilon(x-2)})^+$ (the
  fresh part of the ε-wave), $A=1$, $T=12$ ($w_0=0$) or $T=8$ ($w_0=1$),
  $L=2+c_0T+12/c_0+2$. The ALE speed (`err`) is the mean speed between the first and the
  last remeshing in $[T/2,T]$ (`common.cycle_speed`); also stored are the least-squares
  slope on $[T/2,T]$ (`err_lsq`) and the end-point quotient (`err_end`). The Eulerian speed
  is the least-squares slope on $[T/2,T]$.
- **Travelling-wave runs, groups D, E, G.** Initial datum $U^0_j=\Phi_\varepsilon(x_j-x_0)$, the
  exact ε-wave (`exact.tw_profile_eps`), `exact_initial=True`, $w=\varepsilon w_0$;
  $x_0=(k+\vartheta)h_f$ with $kh_f\approx2$ and lattice phase $\vartheta=0.37$ unless it is
  varied. A monitor records at every step (ALE) or every few steps (Eulerian):
  - $E=\max_n\max_j|U^n_j-\Phi_\varepsilon(x^n_j-x_0-c_\varepsilon t_n)|$;
  - the front position error $s_n=x_f^n-x_0-c_\varepsilon t_n$; the speed error is the
    least-squares slope of $s$ on $[T/2,T]$;
  - $E_s$: $E$ with the shift $s_n$ removed;
  - $\mu=\max_n|v_n-c_\varepsilon|$ and $\mu_1$, the same for $t_n\ge1$ ($v_n$ = mesh velocity).
- **Radial extinction (B4, A, H).** Exact self-similar solution of the limit problem
  (`exact.SelfSimilar`, $w_0=0.5$, extinction time 1) as initial datum, with the data of the
  paper (`solver.initial_data`); ALE with `ell0=1`, `remesh_distance=0.1`.
- **Front and shift.** The front is the level set $\{u=\varepsilon\}$, at distance
  $\approx\varepsilon/\sqrt{2M}$ from the free boundary of the limit problem; "shifted" errors
  have this shift removed.
- **Rigid-mesh runs (B3, D4, F "rigid", G).** The whole mesh translates (`ell0=1e6`,
  `remesh_distance=1e9`), with a prescribed velocity (`fixed_velocity`) or the implicit one.
  The cells in the tapered zones at the ends of $[0,L]$ are then compressed far below $h_f$,
  and these experiments rely on the mesh not being rebuilt: they pass
  `min_cell=common.RIGID_MIN_CELL` (= 0), which switches off the remeshing of compressed
  cells that `simulate_ale` does by default (`min_cell=0.5`). The travel distance is kept
  to about 5–6 by the choice of $T$.
- **Number of trial velocities.** Each trial velocity of the implicit mesh velocity costs
  one assembly and two tridiagonal solves; the counters `n_trials`, `n_solves`,
  `n_secant_fail` of the result record them. In the travelling-wave runs the iteration
  uses three trials per step (3.00–3.01 on average: the previous velocity, a perturbed
  one, one secant update that meets the tolerance $10^{-3}h_f$), i.e. six tridiagonal
  solves per step. For the accelerating fronts of B4 it uses more when the step is large:
  4.1–4.2 trials per step at $\Delta t=0.04$, 3.4–3.7 at 0.02, 3.07–3.13 at 0.01, and at most
  3.04 for $\Delta t\le0.005$ (all converged: `n_secant_fail = 0`).

## A. Meshes (`meshes.py`, `mesh_snapshots.py`)

**A5 `A5_mesh_snapshots.png` (Fig. `fig:A5`).** The solution ($u$, $Y$) and the mesh (one
tick per node) at four times, ALE, ε = 0.02, $h_f=\varepsilon/16$, $h_c=0.05$. Left: planar
travelling wave ($w_0=0$, 1464 nodes at all times). Right: radial extinction of the
self-similar ball ($N=2$, $w_0=0.5$, 1012 nodes). Insets: the reaction layer in units of ε.

**A1 `A1_spacetime_nodes_tw.png`.** Space–time diagram of the mesh nodes for the planar
travelling wave ($w_0=0$, $A=1$, ε = 0.01, $0\le t\le4$; $h_f=\varepsilon/16$, $h_c=0.04$, chosen so
that individual coarse nodes are visible). Each node is drawn as its trajectory $x_i(t)$,
coloured by the local mesh size; the orange line is the front and the ticks on the right
border mark remeshings. Left: ALE ($\Delta t=0.01$); the band, the transition and the
coarse nodes near the front translate with the flame, and one remeshing occurs. Right:
BDF solver; the nodes are fixed and the mesh is rebuilt 31 times. Bottom row: zoom on
the fresh edge of the band.

**A2 `A2_mesh_size.png`.** Local mesh size $h_{i+1/2}$ (log scale) for the ALE solver with
the default mesh (ε = 0.01) at four times; triangles mark the front. Left: planar
travelling wave ($\Delta t=0.01$): the plateau $h=h_f$ and the geometric ramps move with
the front; behind the flame the cells are stretched (decay zone of the mesh velocity),
ahead of it they are slightly compressed until the next remeshing. Right: radial
extinction ($N=2$, $w_0=0.5$, $\ell_0=1$, remeshing every 0.1, $\Delta t=0.005$).

**A3 `A3_node_count.png`.** Planar travelling wave, $w_0=0$, domain $[0,27.09]$, default
mesh, $t\le4$. (a) Nodes against time, ALE and BDF, ε = 0.01 and 0.001. (b) Nodes per
zone against ε, with the formula: band $=(K_b+K_f+2K_m)\,\mathrm{pts\_per\_eps}+1=881$;
transition $=2\lceil\log(h_c/h_f)/\log1.1\rceil$; coarse $\approx(L-\text{band}-\text{transition})/h_c$.
(c) Time steps to reach $t=4$.

| ε | total | band | transition | coarse | ALE: steps / tridiagonal solves / s | BDF: steps / LU / RHS evaluations / remeshings / s |
|---|---|---|---|---|---|---|
| 0.01 | 3573 | 882 | 58 | 2633 | 400 / 2404 / 1.1 | 932 / 1178 / 4512 / 31 / 5.5 |
| 0.005 | 3615 | 882 | 72 | 2661 | 400 / 2404 / 1.2 | 2018 / 2734 / 9913 / 62 / 11.7 |
| 0.002 | 3650 | 882 | 90 | 2678 | 400 / 2404 / 1.3 | 6533 / 8881 / 31773 / 155 / 39.4 |
| 0.001 | 3671 | 882 | 106 | 2683 | 400 / 2404 / 1.6 | 22687 / 30504 / 107801 / 309 / 122.1 |

The formula gives 3573, 3614, 3650 and 3670 nodes. The seconds are those of unoptimized
Python implementations on one machine; the steps, solves and factorizations are the
implementation-independent measures of cost. An ALE step costs six tridiagonal solves
($2404=6\cdot400+4$), an Eulerian step of (S) two, a BDF step on average 1.3 sparse LU
factorizations of the $2n\times2n$ Jacobian and 4.8 right-hand side evaluations.

**A4 `A4_radial_meshes.png`.** Node trajectories $r_i(t)$ for the ALE solver in radial
symmetry ($N=2$, $w_0=0.5$, ε = 0.01, $h_c=0.04$). Left: extinction of the self-similar
ball ($T=1$, $t\le0.9$, $\Delta t=0.005$, $\ell_0=1$, remeshing distance 0.1): 180 steps, 16
remeshings. Right: expanding flame from a burnt kernel of radius $r_0=1.5$
($u_0=\tanh(\sqrt{2M}(r-1.5))^+$, $L=12$, $\Delta t=0.01$, default $\ell_0$): 600 steps, 6
remeshings.

## B. Convergence (`convergence.py`)

**B1 `B1_tw_speed_eulerian_vs_ale.png`.** Speed error $|c_h-c_\varepsilon|$ against ε for the
planar travelling wave, $w_0=0$ and $w_0=1$, $\Delta t\in\{0.02,0.01,0.005\}$: Eulerian scheme
(S) (dashed squares) and ALE (solid circles). Right panel: relative error against
$\Delta t/\varepsilon$ (Eulerian also with $\Delta t=0.002$, 0.001). For $w_0=1$ the ALE runs use
`ell0=4` (kind `ale`) and the default `ell0=2` (kind `ale_default`).

Signed errors $c_h-c_\varepsilon$:

| scheme | $w_0$ | Δt | ε=0.01 | ε=0.005 | ε=0.002 | ε=0.001 |
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

- The Eulerian error is a function of $\Delta t/\varepsilon$ alone (identical values along the
  diagonals with equal $\Delta t/\varepsilon$).
- The ALE error does not depend on Δt. Its size and its variation with ε are set by the
  lattice phase of the front (see D3); the three speed estimators agree to 2–3 digits.
- With the default $\ell_0=2$ and $w_0=1$ ($c=1.63$) a term of about $0.007\,\Delta t$ from the
  non-rigid zone remains; `ell0=4` removes it (see D1).

**B2 `B2_ale_mesh_sizes.png`.** ALE speed error ($w_0=0$, $\Delta t=0.01$, B1 datum) against
$h_f/\varepsilon$ at $h_c=0.0025$ (left) and against $h_c$ at $h_f=\varepsilon/32$ (right).

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

Below $10^{-4}$ the values are dominated by the phase-dependent pinning error and do not
give a rate (sign changes); the measurement at a fixed phase is D2.

**B3 `B3_rigid_mesh_dt.png`.** Whole mesh translating rigidly with $v=c_\varepsilon$
(`fixed_velocity=c_eps`, `ell0=1e6`, `remesh_distance=1e9`, `min_cell=0`), default mesh;
$T=10$ ($w_0=0$) and 3.6 ($w_0=1$), front travel about 6.

| $w_0$ | ε | Δt=0.1 | 0.05 | 0.02 | 0.01 | 0.005 |
|---|---|---|---|---|---|---|
| 0 | 0.005 | -1.8e-06 | +3.0e-07 | +1.1e-06 | +1.0e-06 | +9.9e-07 |
| 0 | 0.001 | -5.1e-07 | -1.3e-07 | +3.1e-08 | +3.1e-08 | +2.9e-08 |
| 1 | 0.005 | -2.7e-05 | -7.7e-06 | +5.1e-06 | +7.4e-06 | +2.2e-06 |
| 1 | 0.001 | -3.5e-06 | -4.9e-06 | -4.0e-07 | +1.1e-07 | +1.0e-07 |

With the mesh velocity prescribed, the front is dragged by the mesh (D4): for
$\Delta t\gg\varepsilon$ the front speed is mostly the mesh speed, so this test shows locking to
the mesh rather than accuracy. D1 is the test of accuracy.

**B4 `B4_radial_extinction_dt.png` (Fig. `fig:B4`).** ALE, exact self-similar extinction
($w_0=0.5$, $T=1$, $t\le0.75$, $L=1.5R+1$, $\ell_0=1$, remeshing distance 0.1, default mesh).
Shifted front error $\max_{t\ge0.05}|r_h(t)+\varepsilon/\sqrt{2M}-R\sqrt{T-t}|$ against Δt.

| N | ε | Δt=0.04* | 0.02* | 0.01 | 0.005 | 0.0025 | error/Δt at 0.0025 |
|---|---|---|---|---|---|---|---|
| 2 | 0.01 | 6.41e-02 | 3.07e-02 | 1.55e-02 | 7.98e-03 | 4.28e-03 | 1.71 |
| 2 | 0.002 | 6.19e-02 | 2.95e-02 | 1.47e-02 | 7.42e-03 | 3.81e-03 | 1.52 |
| 2 | 0.001 | 6.19e-02 | 2.96e-02 | 1.49e-02 | 7.59e-03 | 3.98e-03 | 1.59 |
| 3 | 0.01 | 7.08e-02 | 3.40e-02 | 1.72e-02 | 8.88e-03 | 4.78e-03 | 1.91 |
| 3 | 0.002 | 6.90e-02 | 3.31e-02 | 1.67e-02 | 8.56e-03 | 4.55e-03 | 1.82 |
| 3 | 0.001 | 6.80e-02 | 3.21e-02 | 1.57e-02 | 7.57e-03 | 3.54e-03 | 1.42 |

\* The actual steps are $0.75/19=0.0395$ and $0.75/38=0.0197$ (`dt_actual` in the JSON).
First order in Δt, with a constant between 1.4 and 1.9.

## C. Limit problem and examples (`limit_examples.py`)

**C1 `C1_ale_vs_limit.png` (Fig. `fig:C1`).** Radial data $u_0=a(1-r^2/R_0^2)^+$,
$a=\sqrt{2M}R_0/2$, $R_0=2$, $w_0=0.5$, computed up to the time when the limit radius is
$R_0/2$ ($t=0.565$ for $N=2$, $0.393$ for $N=3$). The limit problem is solved by front
fixing (`limit.solve_radial`, $n=1600$). Left: fronts for $N=2$. Middle: shifted front
error. Right: $\max_r|u_h^+-u|$ at the final time. ALE with $\Delta t=0.002$, $0.0005$ and
$0.000125$ (`ell0=1`, `remesh_distance=0.1`); grey: BDF solver, read from
`results/limit_comparison.json`.

| N | ε | Δt | front error | shifted | sup | steps |
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

The unshifted front error is $\approx\varepsilon/\sqrt{2M}$ in all cases. For $\varepsilon\ge0.01$ the
ε part of the error dominates; for small ε the Δt part dominates at the larger steps.

**C2 `C2_self_similar_convergence.png` (Fig. `fig:C2`).** Same data as C1, followed until
the limit radius is $\max(0.15,30\varepsilon)$. Top row: $s(t)/\sqrt{T-t}$ against $T-t$; black:
limit problem with $T$ fitted from $s^2\sim R^2(T-t)$ on $s<0.3$ ($T=0.725090$ for $N=2$,
$0.505756$ for $N=3$); dashed: ALE, ε = 0.01 ($\Delta t=10^{-4}$) and 0.002
($\Delta t=2.5\cdot10^{-5}$), using $r_h+\varepsilon/\sqrt{2M}$ with $T_\varepsilon$ fitted in the same way;
dotted: ε = 0.002, BDF. Domain $L=R_0+0.5$. Bottom row: rescaled profiles and
$\kappa\phi(\eta)$.

| $s$ | $s/\sqrt{T-t}$, N=2 | N=3 |
|---|---|---|
| 1.5 | 2.4474 | 2.9202 |
| 1.0 | 2.4922 | 2.9716 |
| 0.5 | 2.5102 | 2.9971 |
| 0.25 | 2.5130 | 3.0020 |
| $R$ (exact) | 2.513927 | 3.003951 |

| | $T$ (limit) | BDF, ε=0.002 | ALE, ε=0.01 | ALE, ε=0.002 |
|---|---|---|---|---|
| N=2 | 0.725090 | 0.725024 | 0.72464 | 0.72533 |
| N=3 | 0.505756 | 0.505712 | 0.50543 | 0.50594 |

The ALE curve for ε = 0.002 dips by 1.4 % ($N=2$) and 1.6 % ($N=3$) below $R$ at
$T-t\approx10^{-2}$ (ε = 0.01: 0.3 %, 0.5 %): near $r=0$ the mesh velocity is tapered to
zero (see the limitations below). The JSON stores, for every ALE run, `status`, `reason`,
`n_secant_fail` and `n_mesh_remesh`; a failed run stops the script.

**C3 `C3_quenching.png` (not in the paper).** Expanding flame from a burnt kernel of
radius $r_0$: $u_0=\tanh(\sqrt{2M}(r-r_0))^+$, $A=1$, $w_0=0.5$, $N=2$, $c_0=1.0607$,
$L=12$. Left: front radius $r_f(t)$ for several $r_0$ (ε = 0.01, ALE, $\Delta t=0.001$).
Right: critical radius $r_*$ by bisection on $[0.6,0.8]$ with tolerance $10^{-3}$.

Classification of a kernel run (`run_kernel`), to $t=6$:

- *quenching*: the front disappeared (solver status `extinct`);
- *ignition*: the run ended with a front radius above $r_0+0.5$;
- *undecided*: the run ended with a front and a smaller radius; the run is then repeated
  to $t=12$;
- *unresolved*: as undecided, and the front had left the fine band of the ALE mesh.

A step of the bisection is decided only by *quenching* or *ignition*. A run that fails
numerically, is unresolved, or is still undecided at $t=12$ aborts the bisection: the
steps are written to `C3_quenching_aborted.json`, `C3_quenching.json` is left untouched
and the script stops with an error. Every bisection step is stored in `history` as a
record with `r0`, `t_end`, `outcome`, `reason` (the criterion that decided, with the time
of extinction or the final radius), `status`, `r_min`, `t_min`, `steps`, `seconds`,
`n_secant_fail`, `n_mesh_remesh`, `n_front_off_band`, `t_front_off_band`.

Two named constants of `limit_examples.py` control this:

- `REMESH_OFF_BAND = True`: the ALE runs use `simulate_ale(remesh_off_band=True)`. A
  collapsing kernel enters the zone $r<0.5$ where the mesh velocity is tapered; the band
  nodes then move more slowly than the front, and for small ε (narrow band) the front
  leaves the fine band. Without the option the front of the kernels with ε = 0.00125
  stalls at $r\approx0.26$ in cells that do not resolve the reaction layer and never
  disappears (for $\Delta t$ = 0.001, 0.0005 and 0.00025 alike); for ε = 0.005 and 0.0025
  the last part of the collapse (from $r\approx0.09$ and $0.26$) runs outside the band.
  With the option the mesh is rebuilt around the front and the collapse is followed until
  the front disappears. The option does not act in runs where the front stays in the band
  (all the igniting kernels).
- `QUENCH_BELOW = None`: with a value $q$, a run that ends with a front radius below
  $q\,r_0$ is also counted as quenching although the front is still present (the
  `reason` says so).

| ε | Δt | $r_*\in$ |
|---|---|---|
| 0.01 | 0.001 | [0.7148, 0.7156] |
| 0.005 | 0.001 | [0.7148, 0.7156] |
| 0.0025 | 0.001 | [0.7148, 0.7156] |
| 0.00125 | 0.001 | [0.7141, 0.7148] |
| 0.00125 | 0.00025 | [0.7156, 0.7164] |
| 0.01 | BDF | [0.7148, 0.7156] |

Kernels that ignite first shrink (down to $r\approx0.53$ for $r_0$ just above critical) and
then expand with speed tending to $c_0$. The heuristic value $1/c_0=0.943$ overestimates
$r_*$ by 30 %.

## D. Travelling wave from $\Phi_\varepsilon$ (`tw_error.py`)

$L=x_0+c_\varepsilon T+12/c_\varepsilon+2$, $T=8$ ($w_0=0$) or 4 ($w_0=1$) unless stated.

**D1 `D1_ale_error_vs_dt.png` (Fig. `fig:D1`).** ALE, error against Δt. Rows $w_0=0,1$;
columns $E$, $E_s$, $\mu_1$; ε = 0.01, 0.003, 0.001 with the default $\ell_0$; for $w_0=1$
also `ell0=4` (dashed) and a shorter rigid zone, `rigid_extra=0.5` (dotted).

| $w_0$ | ε | quantity | Δt=0.1 | 0.05 | 0.02 | 0.01 | 0.005 | 0.0025 |
|---|---|---|---|---|---|---|---|---|
| 0 | 0.01 | E | 4.7e-4 | 4.5e-4 | 4.5e-4 | 4.4e-4 | 4.4e-4 | 4.4e-4 |
| 0 | 0.001 | E | 4.8e-4 | 4.6e-4 | 4.5e-4 | 4.5e-4 | 4.4e-4 | 4.4e-4 |
| 0 | all | $E_s$ | 6.1–7.3e-5 | 4.1e-5 | 4.3e-5 | 4.4e-5 | 4.4e-5 | 4.4e-5 |
| 0 | all | $\mu_1$ | 1.17e-4 | 1.15e-4 | 1.13e-4 | 1.13e-4 | 1.12e-4 | 1.12e-4 |
| 0 | all | speed error | +9.0–9.4e-5 | +8.5e-5 | +8.1e-5 | +7.9e-5 | +7.9e-5 | +7.8e-5 |
| 0 | all | $\mu$ (all steps) | 2.7e-4 | 3.7e-4 | 5.7e-4 | 7.8–8.1e-4 | 1.1e-3 | 1.5–1.6e-3 |
| 1 | 0.01 | E | 5.3e-3 | 3.8e-3 | 3.1e-3 | 2.8e-3 | 2.7e-3 | 2.6e-3 |
| 1 | 0.001 | E | 7.4e-3 | 4.8e-3 | 3.4e-3 | 2.9e-3 | 2.7e-3 | 2.6e-3 |
| 1 | 0.01 | $E_s$ | 4.3e-4 | 2.0e-4 | 8.6e-5 | 8.7e-5 | 8.9e-5 | 9.0e-5 |
| 1 | 0.001 | $E_s$ | 7.1e-4 | 3.3e-4 | 1.2e-4 | 7.9e-5 | 8.0e-5 | 8.2e-5 |
| 1 | 0.01 | $\mu_1$ | 1.1e-3 | 7.1e-4 | 5.1e-4 | 4.3e-4 | 4.0e-4 | 3.8e-4 |
| 1 | 0.001 | $\mu_1$ | 1.6e-3 | 9.5e-4 | 6.0e-4 | 4.9e-4 | 4.4e-4 | 4.2e-4 |
| 1, ℓ0=4 | all | E / $E_s$ / $\mu_1$ | 2.4–2.5e-3 / 8.4–9.2e-5 / 3.6–3.9e-4 | | flat | | | |
| 1, `rigid_extra=0.5` | 0.001 | E / $\mu_1$ | 3.4e-2 / 7.1e-3 | 1.8e-2 / 3.6e-3 | 8.4e-3 / 1.6e-3 | 5.4e-3 / 1.0e-3 | 4.0e-3 / 7.0e-4 | 3.2e-3 / 5.5e-4 |

- $w_0=0$: $E$, $E_s$ and $\mu_1$ do not depend on Δt or ε from $\Delta t=0.0025$ to 0.1
  ($\Delta t/\varepsilon$ up to 100). $E$ is dominated by the drift of the front at the pinned
  speed (D3).
- $\mu$ over all steps is set by the first steps (the discrete profile relaxing from
  $\Phi_\varepsilon$ to the discrete wave) and grows as Δt decreases.
- $w_0=1$, default $\ell_0=2$: a Δt term from the non-rigid zone remains,
  $\mu_1-\mu_1(\ell_0{=}4)\approx0.007\,\Delta t$ (ε = 0.01) to $0.012\,\Delta t$ (ε = 0.001); with
  `rigid_extra=0.5` it is $0.04$–$0.07\,\Delta t$. With `ell0=4` the change from
  $\Delta t=0.0025$ to 0.1 is below 1 %.

**D2 `D2_ale_error_vs_mesh.png` (Fig. `fig:D2`).** ALE, $\Delta t=0.01$, $\vartheta=0.37$.
Left and middle: against $h_f/\varepsilon\in\{1/2,\dots,1/32\}$ at $h_c=0.0025$; right: against
$h_c$ at $h_f=\varepsilon/32$.

| $h_f/\varepsilon$ | 1/2 | 1/4 | 1/8 | 1/16 | 1/32 |
|---|---|---|---|---|---|
| $w_0=0$: $E_s$ | 3.3e-3 | 7.7e-4 | 1.9e-4 | 4.5e-5 | 1.1e-5 |
| $w_0=0$: speed error | +5.9e-3 | +1.37e-3 | +3.3e-4 | +8.3e-5 | +2.2e-5 |
| $w_0=1$: $E_s$ | 7.3e-3 | 1.6e-3 | 3.8e-4 | 8.9e-5 | 3.8–6.3e-5 |
| $w_0=1$: speed error | +3.3e-2 | +7.4e-3 | +1.8e-3 | +4.8–5.3e-4 | +1.6–2.1e-4 |

(The values for ε = 0.01, 0.003, 0.001 agree to 1–2 digits; ranges are given where they
differ.) Second order in $h_f/\varepsilon$ (ratio 4.0–4.4 per halving), independent of ε. The
dependence on $h_c$ is roughly $O(h_c^2)$ for $w_0=1$ and below $4\cdot10^{-5}$ in the speed
for $w_0=0$.

**D3 `D3_phase_pinning.png` (Fig. `fig:D3`).** ALE, ε = 0.01 (lines, 16 phases) and 0.001
(crosses, 8 phases), $\Delta t=0.02$, $T=32$ ($w_0=0$) and $T=8$ ($w_0=1$, with `ell0=4`);
mean of $v_n-c_\varepsilon$ over the last eighth of the run, as a function of the lattice phase
$\vartheta\in[0,1)$ of the initial front inside a band cell. Values for ε = 0.01,
min / max / mean over the phase:

| $h_f/\varepsilon$ | $w_0=0$ | $w_0=1$ |
|---|---|---|
| 1/8 | −6.08e-4 / +3.84e-4 / +1.2e-5 | −3.52e-3 / +2.03e-3 / −5.1e-5 |
| 1/16 | −1.45e-4 / +1.02e-4 / +9.9e-6 | −9.24e-4 / +4.54e-4 / −6.3e-5 |
| 1/32 | −2.90e-5 / +3.26e-5 / +9.5e-6 | −2.81e-4 / +6.27e-5 / −6.6e-5 |

With the implicit velocity the front keeps its phase, and the realised speed is
$c_\varepsilon+g(\vartheta)$. The oscillating part of $g$ scales like $(h_f/\varepsilon)^2$ (half-ranges
4.96e-4, 1.24e-4, 3.08e-5 for $w_0=0$ and 2.78e-3, 6.89e-4, 1.72e-4 for $w_0=1$) and does
not depend on ε; its mean over the phase is an offset that does not depend on $h_f$
(about $+10^{-5}$ and $-6\cdot10^{-5}$), i.e. small compared with the amplitude at
$h_f=\varepsilon/8$ and $\varepsilon/16$ but not zero. The mesh velocity relaxes on the time scale
$4/c^2$, which is why the late window is used.

**D4 `D4_drag.png`.** Rigid mesh with prescribed $v=c_\varepsilon+\delta v$, $\delta v=10^{-3}$, from
$\Phi_\varepsilon$; $c_h$ from the slope of $s$ on $[T/2,T]$ ($T=6$, 3). Ratio
$(c_h-c_\varepsilon)/\delta v$:

| Δt/ε | 0.03 | 0.1 | 0.3 | 1 | 3 | 10 | 30 | 100 |
|---|---|---|---|---|---|---|---|---|
| $w_0=0$ (ε = 0.01 / 0.003 / 0.001) | 0.04 | 0.07 / 0.08 / 0.09 | 0.14 / 0.15 / 0.15 | 0.33 | 0.60 | 0.83–0.85 | 0.94 | 0.98 |
| $w_0=1$ | 0.21 | 0.34 / 0.34 / 0.37 | 0.46–0.48 | 0.67–0.68 | 0.82–0.85 | 0.95 | 0.98 | 0.99 |

The ratio is a function of $\Delta t/\varepsilon$ alone and tends to 1 for $\Delta t\gg\varepsilon$. For
$w_0=1$ at small $\Delta t/\varepsilon$ the pinning amplitude is comparable to $\delta v$, so those
values are not conclusive.

**D5 `D5_eulerian_error.png`.** Eulerian scheme (S) from $\Phi_\varepsilon$, $T=4$ ($w_0=0$) and
2.5 ($w_0=1$), monitor every $\max(1,0.005/\Delta t)$ steps.

| Δt/ε | 0.01 | 0.03 | 0.1 | 0.3 | 1 |
|---|---|---|---|---|---|
| $w_0=0$: $\lvert c_h-c_\varepsilon\rvert/c_\varepsilon$ | 4.2–4.8e-3 | 1.4e-2 | 4.5e-2 | 0.122 | 0.30 |
| $w_0=0$: $E_s$ | 1.4e-3 | 4.0e-3 | 1.3e-2 | 3.4e-2 | 8.3e-2 |
| $w_0=1$: relative speed error | 1.4–1.6e-2 | 4.4–4.5e-2 | 0.128 | 0.277 | 0.50 |
| $w_0=1$: $E_s$ | 5.8e-3 | 1.7e-2 | 4.9e-2 | 0.11 | 0.22 |

(Ranges over ε = 0.01, 0.003, 0.001.) The scheme is slower than $c_\varepsilon$; the error is a
function of $\Delta t/\varepsilon$ alone and linear for $\Delta t/\varepsilon\lesssim0.1$. These are values on a
short window; the settled values are those of E1. Right panel ($\Delta t=0.01\,\varepsilon$): the
dependence on $h_f/\varepsilon$ is not resolved at this Δt.

**F (three speed estimators).** B1 set-up ($w_0=0$, $T=12$, default mesh); speed error by
cycle speed / least squares on $[T/2,T]$ / end points. Kind `floor`: default ALE; kind
`floor_rigid`: rigid mesh with implicit velocity (no remeshing).

| kind | ε | Δt | cycle | lsq | end points |
|---|---|---|---|---|---|
| `floor` | 0.01 | 0.02 | +1.66e-5 | +1.65e-5 | +1.62e-5 |
| | 0.01 | 0.005 | +1.77e-5 | +1.77e-5 | +1.77e-5 |
| | 0.001 | 0.02 | +3.74e-5 | +3.73e-5 | +3.70e-5 |
| | 0.001 | 0.005 | +3.73e-5 | +3.72e-5 | +3.72e-5 |
| `floor_rigid` | 0.01 | 0.02 / 0.005 | – | +1.15e-5 / +1.64e-5 | +1.12e-5 / +1.63e-5 |
| | 0.001 | 0.02 / 0.005 | – | +3.05e-5 / +3.53e-5 | +3.02e-5 / +3.53e-5 |

## E. Transients and accumulated mesh-velocity error (`transient.py`)

- **E1** (`D6_transients.png`, left; Fig. `fig:D6`): Eulerian scheme (S) from $\Phi_\varepsilon$,
  $w_0=0$, ε ∈ {0.02, 0.01}, Δt ∈ {4ε², ε², ε²/4}, $T=60$. Speed error by least squares on
  successive windows of length 6, divided by $\Delta t/\varepsilon$. The speed relaxes on the time
  scale $4/c^2\approx12$. On the window [54, 60] the quotient is −0.2161, −0.2285, −0.2269
  for ε = 0.01 and $\Delta t/\varepsilon$ = 0.0025, 0.01, 0.04, and −0.2348, −0.2322, −0.2225 for
  ε = 0.02 and $\Delta t/\varepsilon$ = 0.005, 0.02, 0.08; on the first window [0, 6] it is between
  −0.283 and −0.304. The prediction of `truncation_check.py` is −0.2345 (ε = 0.01) and
  −0.2351 (ε = 0.02).
- **E2** (`D6_transients.png`, right): ALE from $\Phi_\varepsilon$, ε = 0.01, $\Delta t=0.02$, $T=24$:
  block means of $v_n-c_\varepsilon$ for the lattice phases ϑ = 1/32, 0.37, 17/32, 25/32.
- **E3**: the quantities of the ALE theorem, ALE from $\Phi_\varepsilon$, $T=8$, ϑ = 0.37:

| ε | Δt | $M_T=\Delta t\sum\lvert v_n-c_\varepsilon\rvert$ | $\max_n\lvert v_n-c_\varepsilon\rvert$ | $E$ | $E_s$ |
|---|---|---|---|---|---|
| 0.01 | 0.1 | 8.16e-4 | 2.74e-4 | 4.68e-4 | 6.10e-5 |
| 0.01 | 0.02 | 7.78e-4 | 5.65e-4 | 4.46e-4 | 4.29e-5 |
| 0.01 | 0.0025 | 7.70e-4 | 1.45e-3 | 4.42e-4 | 4.44e-5 |
| 0.001 | 0.1 | 8.31e-4 | 2.76e-4 | 4.80e-4 | 7.31e-5 |
| 0.001 | 0.02 | 7.80e-4 | 5.78e-4 | 4.50e-4 | 4.26e-5 |
| 0.001 | 0.0025 | 7.67e-4 | 1.59e-3 | 4.43e-4 | 4.45e-5 |

- **E4**: over 16 lattice phases (ε = 0.01, $\Delta t=0.02$, $T=8$): $E$ from 1.3e-5 to 9.5e-4,
  $E_s$ from 1.3e-5 to 9.6e-5, speed error from −1.7e-4 to +1.0e-4, $M_T$ from 2.5e-5 to
  1.7e-3. For ε = 0.001 (8 phases): $E$ from 1.1e-4 to 7.7e-4, speed error from −1.3e-4 to
  +9.6e-5.

## G. Uniform mesh inside the hypotheses of Theorem `an:thm:ale` (`in_hypotheses.py`)

Uniform mesh with $h=\varepsilon^2/8$ and $\varepsilon^2/16$, ε ∈ {0.1, 0.05}, translating rigidly with
a prescribed velocity, from $\Phi_\varepsilon$, $T=3$. With $v=c_\varepsilon$: $E\le7.6\cdot10^{-6}$ for
$\Delta t=0.1$ and 0.01 ($\rho_hT$ = 0.19–0.38). With $v=c_\varepsilon+0.01$: $E$ = 6.6e-3 / 1.0e-3
(ε = 0.1, $\Delta t$ = 0.1 / 0.01) and 9.8e-3 / 2.0e-3 (ε = 0.05), below $K\mu T\approx1.7\cdot10^{-2}$.

## H. Setting of the convergence theorem (`convergence_thm.py`)

Radial self-similar solution $v$ of the limit problem ($N=2$, $w_0=0.5$, extinction time
1), data of the paper, $0\le t\le0.5$, no remeshing, $h_f=\varepsilon/16$, ε ∈ {0.02, 0.01, 0.005},
$\Delta t/\varepsilon$ ∈ {1/64, 1/16, 1/4, 1, 4}. $E_v=\max_{n,j}|U^n_j-v(x^n_j,t_n)|$.

- *Fixed uniform mesh* ($v_n=0$): $E_v$ = 0.022–0.023, 0.076–0.077, 0.177, 0.264,
  0.296–0.302 for the five values of $\Delta t/\varepsilon$, the same for the three ε to two digits;
  the front at $t=0.5$ is behind by 0.002–0.018, 0.06–0.08, 0.21–0.23, 0.43–0.45, 0.57–0.62.
- *Band mesh with prescribed velocity* $v_n=\dot R(t_n)$ (ALE with `fixed_velocity`,
  `remesh=False`, `ell0=0.5`, $L=5$; the cells outside the rigid zone stretch up to 0.044):
  $E_v=0.5\,\varepsilon=\varepsilon w_0$ for $\Delta t\le\varepsilon$ and 0.64–0.67 ε for $\Delta t=4\varepsilon$; the front is
  at −1.01…−1.07 ε from $R(0.5)$ for $\Delta t\le\varepsilon$ (the shift $-\varepsilon/\sqrt{2M}=-0.94\,\varepsilon$
  plus −0.07…−0.13 ε), and −0.33…−0.37 ε for $\Delta t=4\varepsilon$.

## `tab_simple.py`: Eulerian scheme with $\Delta t\propto\varepsilon^2$

Eulerian scheme (S), $w_0=0$. Set-up `B1`: B1 datum, $T=12$, least squares on $[T/2,T]$.
Set-up `Phi`: datum $\Phi_\varepsilon$, $T=6$, slope of $s$.

| ε | Δt | Δt/ε | $c_{(S)}-c_\varepsilon$ (`B1`) | from $\Phi_\varepsilon$ (`Phi`) |
|---|---|---|---|---|
| 0.02 | 4ε² = 1.6e-3 | 0.08 | −1.863e-2 | −2.00e-2 |
| 0.02 | ε² = 4e-4 | 0.02 | −4.85e-3 | −5.21e-3 |
| 0.02 | ε²/4 = 1e-4 | 0.005 | −1.23e-3 | −1.32e-3 |
| 0.01 | 4ε² = 4e-4 | 0.04 | −9.49e-3 | −1.02e-2 |
| 0.01 | ε² = 1e-4 | 0.01 | −2.39e-3 | −2.57e-3 |
| 0.01 | ε²/4 = 2.5e-5 | 0.0025 | −5.67e-4 | −6.13e-4 |

The values are linear in Δt at fixed ε and are a function of $\Delta t/\varepsilon$:
$c_{(S)}-c_\varepsilon\approx-(0.39\text{–}0.42)\,c_\varepsilon\Delta t/\varepsilon$ (`B1`), −0.44 (`Phi`). These are
values on the windows $[6,12]$ and $[3,6]$; the settled values are those of E1.

## `truncation_check.py` and `speed_coefficient.py`

- `../truncation_check.py` evaluates the truncation error of (S) on the exact travelling
  wave and the drift of the speed predicted by the solvability condition (eq.
  `an:eq:solv`), averaged over 16 offsets of the grid. It prints two tables and stores
  them in `truncation_check.json`: `time_step` ($h=\varepsilon/64$, $\Delta t=\kappa\varepsilon^2$; the
  quotient `dc_over_dt_eps` is −0.2351, −0.2345, −0.2342 for ε = 0.02, 0.01, 0.005 with
  κ = 1) and `mesh_size` ($\Delta t=10^{-3}\varepsilon^2$, $h=\varepsilon/8,\varepsilon/16$).
- `speed_coefficient.py` compares $(c_\varepsilon-c_0)/\varepsilon$, with $c_\varepsilon$ from the shooting
  code, with the prediction of Proposition `prop:speed`, for $w_0\in\{-0.5,0,1\}$ and
  $\varepsilon\in\{10^{-2},10^{-3},10^{-4}\}$ (`speed_coefficient.json`).

## Limitations of the code that affect the experiments

1. **Taper near the ends.** The mesh velocity is multiplied by `clip(z/ell_end)` and
   `clip((L-z)/ell_end)`, `ell_end = min(0.5, 0.1 L)`, so that no node crosses $r=0$ or
   $x=L$. When a collapsing front enters this zone the scheme is partly Eulerian and an
   $O(\Delta t\,v/\varepsilon)$ lag appears (the dip in C2); the band cells there are compressed,
   and `simulate_ale` rebuilds the mesh when one falls below `min_cell`·$h_f$ (this
   happens in the last steps of the collapsing kernels of C3; the counter `n_mesh_remesh`
   records it). In this zone the band nodes move more slowly than the front, so a narrow
   band (small ε) can be left behind by the front; `n_front_off_band` counts the steps
   that end with the front outside the band, and `remesh_off_band=True` rebuilds the
   mesh when it happens (used in C3; see there).
2. **Rigid zone and default $\ell_0$ for fast fronts.** For $w_0=1$ ($c=1.63$) and the
   default $\ell_0=2$ a term of about $0.007$–$0.012\,\Delta t$ remains in the mesh velocity;
   `ell0=4` removes it (D1).
3. **The rigid mode** (whole mesh translating) is only usable for moderate travel
   distances: the cells of the tapered end zone at $x=L$ are compressed. B3, D4, F and G
   limit the travel to about 5–6 and run with `min_cell=0`.
4. **Lattice pinning** (D3): with the implicit velocity the speed error depends on the
   lattice phase of the front, with amplitude $\propto(h_f/\varepsilon)^2$.
5. **Remeshing interval of the Eulerian solvers.** `simulate` and `simulate_simple` rebuild
   the mesh every $\min(7.5\,\varepsilon/v_{est},\,$`max_dt_remesh`$)$ with
   `max_dt_remesh = 0.05*t_end` by default, so the sequence of meshes depends on the final
   time unless `max_dt_remesh` is passed.

## Status of the stored results

The files in `results/` were generated before the provenance record existed and before
the following changes of the code, which alter some numbers slightly when the files are
regenerated (the tables above are those of the stored files):

- the Eulerian solvers (`simulate`, `simulate_simple`) now locate the initial front on the
  mesh of the run. This changes the first sample of `front` and, through the first
  remeshing interval, the BDF runs and some Eulerian (S) runs in the last digits
  (`A_meshes.json` BDF rows, `C2_self_similar.json` `bdf_run`, `C3_quenching.json` BDF
  bisection, `results/traveling_wave.json`, `results/radial.json`,
  `results/limit_comparison.json`; Eulerian rows of `B_convergence.json`,
  `D_tw_error.json`, `E_transient.json`, `tab_simple.json`);
- the mesh generator merges a remainder shorter than $h_f/2$ at an end of the domain into
  the neighbouring cell, and `simulate_ale` rebuilds the mesh when a cell is compressed
  below $h_f/2$. This acts only when the band reaches $r=0$: the rows ε = 0.05 of
  `C1_limit_comparison.json` and the collapsing kernels of `C3_quenching.json`.
