# The numerical experiments of the paper

Reference for the scripts that produce the figures, the tables and the numbers quoted in
the article and in its supplementary material: which script produces each item, how to
run it, how long it takes, what it writes. The package itself is documented in
[`../../README.md`](../../README.md).

**Numbering.** Numbers of figures, tables, sections and equations are those of the
drafts in `paper/` at the top of the repository (`main.pdf`; with the prefix S,
`supplement.pdf`). They may still change: the LaTeX label and the output file, given next
to each number, identify the item.

## Running

Everything is run from `numerics/`:

```bash
cd numerics
export OMP_NUM_THREADS=1                          # the scripts start 2 or 4 processes
python -m experiments.paper.convergence           # run, then plot
python -m experiments.paper.convergence plot      # only redraw the figures from the stored JSON
python -m experiments.paper.tw_error run D3 D4    # only some groups
```

- **Arguments.** `meshes`, `convergence`, `limit_examples`, `tw_error`, `transient`,
  `variable_w0` and `work_precision` accept `run` or `plot` and do both without argument.
  `limit_examples`, `tw_error` and `transient` also accept a list of groups after `run`;
  the records of the other groups are kept. The other scripts take no argument.
- **Output.** Data: `results/paper/*.json` (`results/*.json` for the three scripts in
  `experiments/`). Figures: `figures/paper/*.png` (`figures/*.png`). A script overwrites
  its files.
- **Provenance.** Every JSON file carries `provenance = {git, git_dirty, python, numpy,
  scipy, platform}`: under that key if the file is a dictionary, in every record if it
  is a list of records (`experiments/util.py`).
- **Failures.** A numerical failure of a run is never used as a result. `convergence` and
  `tw_error` store the failed case as `{..., "failed": message}` and print `FAILED`;
  `limit_examples` stops. The extinction of the flame is a regular outcome
  (`status = "extinct"`).
- **Times** are those of the last full run of `../../reproduce_paper.sh`, on a 4-core
  machine (92 min; `variable_w0` and `work_precision` were run separately, 9 and 26 min).
  They are approximate and machine dependent.
- **What was checked** when this file was last revised. All the `plot` modes,
  `experiments.plots` and `mesh_snapshots` were run, and reproduce the stored figures
  byte for byte. `speed_coefficient`, `truncation_check`, `in_hypotheses`,
  `convergence_thm` and `../radial.py` were run in full, and reproduce the stored numbers
  exactly (only the run times and the provenance differ). The other `run` modes, the
  selection of groups, `tab_simple.py` and the two movie scripts **were not run again**.

Common set-up: $f(s)=(1-s)^+$, $M(w_0)=(1+w_0)^3/6$; $c_0=\sqrt{2M}/A$ is the speed of the
travelling wave of the limit problem and $c_\varepsilon$ that of $(P_\varepsilon)$ (by shooting). Default
mesh: $h_f=\varepsilon/16$ in a band of width $55\,\varepsilon$, $h_c=0.01$. Solvers: **ALE**
(`ale.simulate_ale`), **Eulerian (S)** (`simple.simulate_simple`), **BDF**
(`solver.simulate`), **limit** (`limit.solve_radial`, `solve_planar`).

## Items of the paper

Scripts are modules of `experiments.paper` unless they start with `experiments.`; JSON
files are in `results/paper/` and figures in `figures/paper/` unless a directory is given.

### Article

| item | script, group | solver | data | figure | time |
|---|---|---|---|---|---|
| Table 1 (`tab:speedcoef`), Proposition 2.2 (`prop:speed`) | `speed_coefficient` | shooting | `speed_coefficient.json` | – | 5 s |
| Eq. (37) (`an:eq:solv`), Section 5.4 | `experiments.truncation_check` | none | `truncation_check.json` | – | 10 s |
| Figure 1 (`fig:A5`) | `mesh_snapshots` | ALE | – | `A5_mesh_snapshots.png` | 5 s |
| Section 6.1: remeshings, numbers of nodes | `meshes`, A | ALE, BDF | `A_meshes.json` | – | 3 min |
| Table 2 (`tab:limit`) | `experiments.limit_comparison` | BDF, limit | `results/limit_comparison.json` | – | 8 min |
| Figure 2 (`fig:C1`) and the ALE errors of Section 6.2 | `limit_examples`, C1 | ALE, limit | `C1_limit_comparison.json`; grey points: `results/limit_comparison.json` | `C1_ale_vs_limit.png` | 1 min |
| Figure 3 (`fig:I1`) and the numbers of Section 6.3 | `variable_w0`, I | BDF, ALE, limit | `I_variable_w0.json` | `I1_variable_w0.png` | 9 min |
| Figure 4 (`fig:instability`), Section 6.4 | `experiments.traveling_wave`, part B; then `experiments.plots` | BDF, limit | `results/traveling_wave.json` | `figures/tw_instability.png` | 21 min |
| Section 6.5: setting of Theorem 5.2 (`cv:thm:main`) | `convergence_thm`, H | (S) on a fixed mesh; ALE with prescribed velocity | `H_convergence_thm.json` | – | 25 s |
| Figure 5 (`fig:B4`) | `convergence`, B4 | ALE | `B_convergence.json` | `B4_radial_extinction_dt.png` | 3 min |
| Section 6.6 (summary of the supplement) | see the next table | | | | |
| Figure 6 (`fig:C2`) and the numbers of Section 6.7 | `limit_examples`, C2 | limit, ALE, BDF | `C2_self_similar.json` | `C2_self_similar_convergence.png` | 2 min |

### Supplementary material (and Section 6.6 of the article)

| item | script, group | solver | data | figure | time |
|---|---|---|---|---|---|
| Section S1, uniform meshes inside the hypotheses of Theorem 5.5 (`an:thm:ale`); the same with a $C^{1,1}$ rate; ε = 0.1 and 0.05 | `in_hypotheses`, G | ALE | `G_in_hypotheses.json`, `G_regular_f.json`, `G_large_eps.json` | – | 45 s |
| Figure S1 (`fig:D1`) | `tw_error`, D1 | ALE | `D_tw_error.json` | `D1_ale_error_vs_dt.png` | 8 min for all the groups of `tw_error` |
| Figure S2 (`fig:D3`) | `tw_error`, D3 | ALE | `D_tw_error.json` | `D3_phase_pinning.png` | |
| Figure S3 (`fig:D2`) | `tw_error`, D2 | ALE | `D_tw_error.json` | `D2_ale_error_vs_mesh.png` | |
| Section S1, "the mesh drags the front": the ratios α | `tw_error`, D4 | ALE, rigid mesh | `D_tw_error.json` | – | |
| Figure S4 (`fig:D6`); the Eulerian lag | `transient`, E1 (left), E2 (right) | Eulerian (S), ALE | `E_transient.json` | `D6_transients.png` | 33 min |
| Section S1: $\mathcal M_T$; range of $E$, $E_s$ over the position of the front | `transient`, E3, E4 | ALE | `E_transient.json` | – | |
| Figure S5 (`fig:J1`), Table S1 (`tab:work`) | `work_precision`, J | ALE, Eulerian (S), BDF, LSODA | `J_work_precision.json` | `J1_work_precision.png` | 26 min |

Three statements of the paper are not backed by a stored file: the number of trial
velocities per step for the accelerating fronts (Section 4.3; the records of group B4 do
not store it), the accuracy of the limit solver against the exact solutions (Section 6.2;
`tests/test_combustion1d.py` checks its order), and the contribution of the coarse mesh
to the Eulerian lag (Section S1, after Figure S4).

To read a number, load the JSON file. For instance, Table 2:

```python
import json
rows = json.load(open("results/limit_comparison.json"))
for eps in (0.05, 0.02, 0.01, 0.005, 0.002, 0.001):
    line = [f"{eps:<6g}"]
    for geom in ("planar", "radial2", "radial3"):
        r = next(r for r in rows if r["geom"] == geom and r["eps"] == eps)
        line += [f"{r['front_err_shifted'] / eps:.3f}", f"{r['sup_err'] / eps:.3f}"]
    print(" ".join(line))
```

## The scripts

To redraw a figure without recomputing: `python -m experiments.paper.<script> plot`, or
`python -m experiments.plots` for the figures in `figures/`. `mesh_snapshots` always
recomputes (5 s).

### `speed_coefficient.py`: Table 1

$(c_\varepsilon-c_0)/\varepsilon$ against the coefficient of Proposition 2.2, for $A=1$,
$w_0\in\{-0.5,0,1\}$, $\varepsilon\in\{0.1,0.05,10^{-2},10^{-3},10^{-4}\}$. Prints the table.
`speed_coefficient.json`: `A`, `kappa`, `rows` = `[{w0, c0, prediction, coef: {eps: value}}]`.

### `../truncation_check.py`: eq. (37)

Truncation error of (S) on the exact travelling wave and the drift of the speed
predicted by the solvability condition, averaged over 16 positions of the grid; no PDE
is solved. `truncation_check.json`: `time_step` ($w_0=0$, $h=\varepsilon/64$, $\Delta t=\kappa\varepsilon^2$;
keys `eps`, `kappa`, `dt`, `h`, `max_tau`, `dc`, `max_tau_eps2_over_dt`, `dc_over_dt_eps`;
the values of eq. (37) are `dc_over_dt_eps` for `kappa = 1`) and `mesh_size`
($\Delta t=10^{-3}\varepsilon^2$, $h=\varepsilon/8,\varepsilon/16$; not quoted in the paper).

### `mesh_snapshots.py`: Figure 1

$u$, $Y$ and the nodes at four times, ALE, ε = 0.02, $h_c=0.05$. Left: planar travelling
wave ($w_0=0$, Δt = 0.01). Right: radial extinction of the self-similar ball ($N=2$,
$w_0=0.5$, Δt = 0.0025, `ell0=1`, `remesh_distance=0.1`). No data file.

### `meshes.py` (group A): numbers of Section 6.1

Planar travelling wave ($w_0=0$, $0\le t\le4$): ALE with Δt = 0.01 and BDF (called
`euler` in the file), default mesh for ε ∈ {0.01, 0.005, 0.002, 0.001}, and $h_c=0.04$
(`mesh = "ill"`) for ε = 0.01. Also two radial ALE runs with ε = 0.01 (extinction;
expanding flame from a kernel of radius 1.5).

- `A_meshes.json`: `tw`, a list of records with `method`, `eps`, `mesh`, `L`, `c_eps`,
  `c_num`, `steps`, `n_remesh`, `parts` and `parts_pred` (nodes in the band, the
  transition and the coarse zone: counted, and predicted by the formula), `n_nodes`,
  `node_t`, `t`, `front`, `trials`, `solves`, `lu`, `fev`, `seconds`; and `radial`,
  `expanding`, `radial_def`, `expanding_def`.
- `A_mesh_*.npz`: the nodes at every step. They are not in git; `run` writes them and
  `plot` needs them (`run` takes 2 min, `plot` 1 min).
- Figures `A1`–`A4`: not in the paper, see below.

### `../limit_comparison.py`: Table 2

$w_0=0.5$; planar datum $u_0=\tanh(\sqrt{2M}(x-2))^+$ on $[0,16]$ up to $t=6$; radial
data $u_0=a(1-r^2/R_0^2)^+$, $a=\sqrt{2M}R_0/2$, $R_0=2$, $N=2,3$, until the radius of the
limit solution is $R_0/2$; ε from 0.05 to 0.001; BDF against the front-fixing solver
($n=1600$). `results/limit_comparison.json`: 18 records with `geom` (`planar`,
`radial2`, `radial3`), `eps`, `t_end`, `front_err`, `front_err_shifted` (the shift
$\varepsilon/\sqrt{2M}$ removed), `sup_err`, `steps`, `seconds`, `t`, `front`, `lim_t`, `lim_s`.
Table 2 shows `front_err_shifted/eps` and `sup_err/eps`.

### `limit_examples.py` (groups C1, C2, C3): Figures 2 and 6

- **C1** (Figure 2): the radial data of Table 2 with ALE, Δt ∈ {0.002, 0.0005, 0.000125},
  `ell0=1`, `remesh_distance=0.1`. `C1_limit_comparison.json`: 36 records with the keys of
  `limit_comparison.json` plus `N`, `dt`, `n_remesh`, `nodes`. The grey points of the
  figure are the BDF records of `results/limit_comparison.json`, which must exist.
- **C2** (Figure 6): the same data followed until the radius of the limit solution is
  $\max(0.15, 30\,\varepsilon)$, on $[0,2.5]$: ALE with ε = 0.01 (Δt = 10⁻⁴) and ε = 0.002
  (Δt = 2.5·10⁻⁵), BDF with ε = 0.002; the extinction times are fitted from
  $s^2\sim R^2(T-t)$. `C2_self_similar.json`: one record for each $N$ with `R`, `kappa`,
  `T`, `lim_t`, `lim_s`, `profiles` (the limit solution when the radius is 1.5, 1, 0.5,
  0.25), `eta`, `phi`, `eps_runs` (`eps`, `dt`, `T_eps`, `t`, `s`, `profiles`, `status`,
  `n_secant_fail`, `n_mesh_remesh`) and `bdf_run`.
- **C3**: not in the paper, see below.

### `variable_w0.py` (group I): Figure 3

A non-constant $w_0$: a planar front in a stratified mixture ($L=24$, $t\le6$, ALE with
Δt = 2.5·10⁻⁴) and a radial pocket of fresh gas ($N=3$, $L=6$, until the radius halves,
Δt = 1.25·10⁻⁴); the formulas are in the docstring. BDF and ALE for ε from 0.05 to 0.001,
against the front-fixing solver with $M$ evaluated at the front. `I_variable_w0.json`:
`eps`, `t_end`, `rows` (`geom`, `solver`, `eps`, `dt`, `status`, `front_err`,
`front_err_shifted`, `sup_err`, `W_err`, `t`, `front`, `steps`, `seconds`) and
`limits[geom]` (`t`, `s`: front of the limit problem; `t_const`, `s_const`: the same with
$w_0\equiv0.5$; `W_front`: $w_0$ at the free boundary).

### `../traveling_wave.py`, `../plots.py`: Figure 4

Part B: the same $u_0$ (the limit travelling wave for $w_0=0$) with
$w_0\in\{-0.5,0,1\}$, ε = 0.005, $t\le20$, BDF. `results/traveling_wave.json`: `B`, three
records with `w0`, `c0`, `t`, `front` and the solution at $t=20$ (`x`, `u`, `Y`); `A`,
see below. `plots.py` draws the figure and computes its dashed lines (the limit
problem) with `solve_planar`. Part A takes most of the 21 minutes.

### `convergence_thm.py` (group H): numbers of Section 6.5

Radial self-similar solution $v$ ($N=2$, $w_0=0.5$, extinction time 1), $0\le t\le0.5$,
no remeshing, $h_f=\varepsilon/16$, ε ∈ {0.02, 0.01, 0.005}, Δt/ε ∈ {1/64, 1/16, 1/4, 1, 4}, on
a uniform mesh that does not move (`mesh = "fixed"`) and on a band mesh whose band moves
with the velocity of the free boundary (`"band"`: `fixed_velocity`, `remesh=False`,
`ell0=0.5`, $L=5$). `H_convergence_thm.json`: 30 records with `mesh`, `eps`,
`dt_over_eps`, `E` ($=\max_{n,j}|U^n_j-v(x^n_j,t_n)|$), `front_error` (at $t=0.5$),
`hmax` (largest cell during the run), `nodes`, `steps`, `seconds`.

### `convergence.py` (groups B1–B4): Figure 5

**B4**: ALE for the extinction of the self-similar ball ($w_0=0.5$, $T=1$, $t\le0.75$,
$N=2,3$, ε ∈ {0.01, 0.002, 0.001}, Δt from 0.04 to 0.0025, `ell0=1`,
`remesh_distance=0.1`). Records of `B_convergence.json` with `group = "B4"`: `N`, `eps`,
`dt`, `dt_actual`, `front_err`, `front_err_shifted` (the quantity of the figure),
`front_err_shifted_end`, `sup_err`, `sup_err_int`, `steps`, `n_remesh`, `nodes`,
`seconds`. The script always runs B1–B3 too (not in the paper, see below).

### `in_hypotheses.py` (group G): numbers of Sections S1 and 6.6

From $\Phi_\varepsilon$, $w_0=0$, ε ∈ {0.1, 0.05}, Δt ∈ {0.1, 0.01}:

- `G_in_hypotheses.json`: uniform mesh $h=\varepsilon^2/k$, $k=8,16$, translating with the
  prescribed velocity $c_\varepsilon+\mu$, $\mu\in\{0,0.01\}$, $T=3$. Keys `eps`, `k`, `dt`, `mu`,
  `h`, `rho_h`, `E`, `E_shift`, `s_max`, `nodes`, `steps`, `seconds`.
- `G_regular_f.json`: the same with $f(s)=3((1-s)^+)^2$, $h=\varepsilon^{3/2}/k$, $k=2,\dots,16$,
  $\mu=0$ (`rho_h` is then $h^2/\varepsilon^3$).
- `G_large_eps.json`: the method itself (default mesh, implicit velocity), $T=8$; keys as
  in group D1 below.

### `tw_error.py` (groups D1–D5, F): Figures S1, S2, S3, drag ratios

Runs from $U^0_j=\Phi_\varepsilon(x_j-x_0)$, $w=\varepsilon w_0$, $w_0\in\{0,1\}$, ε ∈ {0.01, 0.003,
0.001}, $T=8$ ($w_0=0$) or 4 ($w_0=1$) unless stated; $x_0=(k+\vartheta)h_f\approx2$ with
$\vartheta=0.37$ unless it is varied. Every record of `D_tw_error.json` has `group`,
`kind`, `w0`, `eps`, `dt`, `T`, `c_eps` and

- `E` $=\max_n\max_j|U^n_j-\Phi_\varepsilon(x^n_j-x_0-c_\varepsilon t_n)|$, `E_shift` ($E_s$: the same
  with the error of the front removed), `s_max`, `s_end` (error of the front position);
- `speed_lsq` (slope of the front error on $[T/2,T]$), `speed_end`;
- ALE: `mu`, `mu1` ($\max|v_n-c_\varepsilon|$ over all steps, over $t_n\ge1$), `v_tail` (mean of
  $v_n-c_\varepsilon$ on the last eighth of the run), `trials_per_step`;
- `series` (`t`, `E`, `Es`, `s` at 160 times), `steps`, `n_remesh`, `nodes`, `seconds`.

| group | what varies | in the paper |
|---|---|---|
| D1 | Δt from 0.1 to 0.0025; for $w_0=1$ also `variant = "ell0=4"` and `"old rigid zone"` (`rigid_extra=0.5`) | Figure S1 |
| D2h, D2c | $h_f/\varepsilon$ from 1/2 to 1/32 at $h_c=0.0025$; $h_c$ from 0.04 to 0.0025 at $h_f=\varepsilon/32$ (key `mesh`); Δt = 0.01 | Figure S3 |
| D3 | 16 values of $\vartheta$ (key `theta`), $h_f/\varepsilon$ = 1/8, 1/16, 1/32, ε = 0.01 (8 values, $h_f=\varepsilon/16$ for ε = 0.001); Δt = 0.02, $T=32$ ($w_0=0$) or 8 ($w_0=1$, `ell0=4`); the figure shows `v_tail` | Figure S2 |
| D4 | rigid mesh with prescribed velocity $c_\varepsilon+10^{-3}$ (key `dv`), Δt/ε from 0.03 to 100, $T=6$ or 3; `ratio` = `speed_lsq/dv` | the ratios α |
| D5, D5h, F | | not in the paper, see below |

### `transient.py` (groups E1–E4): Figure S4

Records of `E_transient.json`, selected by `group`:

- **E1** (Figure S4, left): Eulerian (S) from $\Phi_\varepsilon$, $w_0=0$, ε ∈ {0.01, 0.02},
  Δt = $k\varepsilon^2$ with $k\in\{1/4,1,4\}$, $T=60$. `windows` = `[t_a, t_b, slope,
  slope/(dt/eps)]`: speed error by least squares on successive windows of length 6. Also
  `k`, `dt`, `steps`, `t`, `s`. One of these runs has 2.4·10⁶ steps and takes 33 min.
- **E2** (Figure S4, right): ALE, ε = 0.01, Δt = 0.02, $T=24$, four values of `theta`;
  `vt`, `vv`: means of $v_n-c_\varepsilon$ over blocks of 0.25 time units.
- **E3**: ALE, $T=8$, ε ∈ {0.01, 0.001}, Δt ∈ {0.1, 0.02, 0.0025}: `M_T`
  $=\Delta t\sum_n|v_n-c_\varepsilon|$, `X_T`, `mu`, `E`, `E_shift`, `s_max`.
- **E4**: ALE, $T=8$, Δt = 0.02, 16 values of `theta` for ε = 0.01 and 8 for ε = 0.001:
  the same keys.

### `work_precision.py` (group J): Figure S5, Table S1

Planar travelling wave from $\Phi_\varepsilon$, $w_0=0$, $T=4$, ε ∈ {0.01, 0.001}. Rows of
`J_work_precision.json`, selected by `kind`: `ale` (`pts` = 8, 16, 32, 64; `dt` from 0.2
to 0.0025), `eul` (`dt` = ε/4, …, ε/1024; up to ε/64 for ε = 0.001), `uni` (uniform mesh,
ε = 0.01), `bdf` and `lsoda` (`rtol` from 10⁻³ to 10⁻⁶). Keys: `front_err` (the error of
the figure and of the table), `sup_err`, `seconds`, `steps`, `nodes`, and the counters
`solves` (scheme (S)), `n_lu`, `n_fev` (BDF), `n_jev`, `n_fev`, `n_remesh` (LSODA). The
column "solves" of Table S1 is `solves`, `n_lu` and `n_jev`, respectively. The runs use 2
processes; `seconds` is the time of the solver call and depends on the machine and on
its load.

## In the repository but not in the paper

None of the following appears in the article or in the supplement.

| what | script | output | time |
|---|---|---|---|
| Node trajectories, local mesh size, node counts against time and ε (A1–A4) | `meshes.py` | `A1_spacetime_nodes_tw.png`, `A2_mesh_size.png`, `A3_node_count.png`, `A4_radial_meshes.png`; entries `radial*`, `expanding*` of `A_meshes.json` | part of the 3 min of `meshes` |
| Speed error of a travelling wave started from the limit profile: Eulerian (S) against ALE; mesh sizes; rigid mesh at speed $c_\varepsilon$ (B1–B3) | `convergence.py` | `B1_tw_speed_eulerian_vs_ale.png`, `B2_ale_mesh_sizes.png`, `B3_rigid_mesh_dt.png`; groups B1, B2h, B2c, B3 of `B_convergence.json` | part of the 3 min of `convergence` |
| Quenching of small burnt kernels: critical radius by bisection (C3) | `limit_examples.py` | `C3_quenching.json`, `C3_quenching.png` | 10 min |
| Figure of the drag ratios (the ratios themselves are quoted) | `tw_error.py`, D4 | `D4_drag.png` | |
| Eulerian (S) from $\Phi_\varepsilon$ on short windows (D5, D5h); three estimators of the speed (F) | `tw_error.py` | `D5_eulerian_error.png`; groups D5, D5h, F of `D_tw_error.json` | part of the 8 min of `tw_error` |
| Flame speed against ε with BDF (part A) | `../traveling_wave.py`, `../plots.py` | `A` in `results/traveling_wave.json`; `figures/tw_speed_convergence.png` | most of the 21 min of `traveling_wave` |
| Figures of the comparison of Table 2 | `../plots.py` | `figures/limit_fronts.png`, `figures/limit_convergence.png` | 5 s |
| Radial self-similar solution with BDF | `../radial.py`, `../plots.py` | `results/radial.json`; `figures/radial_convergence.png`, `figures/radial_snapshots.png` | 1 min |
| Eulerian (S) with Δt ∝ ε² (a table of an earlier draft) | `tab_simple.py` | `tab_simple.json` | not measured |
| Animations of a shrinking and of an expanding flame (BDF) | `../movie_radial.py`, `../movie_expanding.py` | `figures/flame_extinction.mp4`, `flame_extinction_still.png`, `flame_expanding.mp4`, `flame_expanding_still.png` | not measured |
| Rows `mesh_size` of the truncation check; rows with `mu = 0.01` of `G_in_hypotheses.json` | `../truncation_check.py`, `in_hypotheses.py` | | |

Notes:

- `reproduce_paper.sh` runs `../radial.py` and all the groups of the scripts above, but
  not `tab_simple.py` nor the movie scripts.
- `tab_simple.json` is the only stored file without a provenance record: it was written
  by an earlier version of the code and has not been regenerated.
- **C3.** A kernel run counts as quenching only if the front disappeared (status
  `extinct`), and as ignition if it ends with a front radius above $r_0+0.5$; an
  undecided run is repeated up to $t=12$, and a bisection that cannot decide a step is
  aborted (its steps go to `C3_quenching_aborted.json`, and `C3_quenching.json` is left
  untouched). Every step is stored in `critical[i].history`. The runs use
  `remesh_off_band=True` (constant `REMESH_OFF_BAND` of the script), because a collapsing
  kernel enters the zone near $r=0$ where the mesh velocity is tapered.
- **Rigid-mesh runs** (B3, D4, F, G) translate the whole mesh (`ell0=1e6`,
  `remesh_distance=1e9`) and pass `min_cell=0`, so that the cells compressed at the far
  end of the domain do not trigger a remeshing; the travel distance is limited instead.

The limitations of the code that affect the experiments (taper of the mesh velocity near
$r=0$, rigid zone for fast fronts, lattice pinning) are listed in
[`../../README.md`](../../README.md#known-limitations).
