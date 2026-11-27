# combustion1d: numerics for the high activation energy limit

Python code of the article

> J. Fernández Bonder, *The high activation energy limit of an equidiffusional flame
> model: rate of convergence and a convergent front-following scheme* (2026, draft):
> [`../paper/main.pdf`](../paper/main.pdf), with supplementary material
> [`../paper/supplement.pdf`](../paper/supplement.pdf).

The code solves the flame model, written with $u=T_f-T$ and $w=Y-u$,

$$(P_\varepsilon)\qquad u_t=\Delta u-(u+w)\,f_\varepsilon(u),\qquad w_t=\Delta w,\qquad f_\varepsilon(s)=\varepsilon^{-2}f(s/\varepsilon),$$

in one space dimension ($N=1$) and in radial symmetry ($N=2,3$), and the free boundary
problem obtained in the limit $\varepsilon\to0$ when $w/\varepsilon\to w_0$:

$$(P)\qquad u_t=\Delta u\ \text{ in }\{u>0\},\qquad |\nabla u|=\sqrt{2M(w_0)}\ \text{ on }\partial\{u>0\},\qquad M(w_0)=\int_{-w_0}^1(s+w_0)f(s)\,ds .$$

The default reaction rate is $f(s)=(1-s)^+$, for which $M(w_0)=(1+w_0)^3/6$.

| path | contents |
|---|---|
| `combustion1d/` | the package: the scheme of the article on a moving band mesh (ALE) and on meshes that do not move, a reference BDF solver, exact solutions, a solver for $(P)$ |
| `examples/` | three short examples |
| `experiments/`, `experiments/paper/` | the scripts that produce every figure and table of the article and of the supplement; see [`experiments/paper/README.md`](experiments/paper/README.md) |
| `results/`, `figures/` | the stored data (JSON) and figures (PNG) |
| `tests/`, `reproduce_paper.sh` | the test suite; one script that reruns everything |

## Installation

Python ≥ 3.10 with NumPy, SciPy, Matplotlib and pytest. Everything is run **from this
directory** (`numerics/`):

```bash
cd numerics
python -m pip install -r requirements.txt
```

`requirements.txt` gives minimal versions. `requirements-lock.txt` pins the versions with
which the stored results and figures were produced (Python 3.13, NumPy 2.5.3, SciPy 1.18.1,
Matplotlib 3.11.2, pytest 9.1.1); use it instead to reproduce them exactly.

There is nothing to build: `python -m ...` puts the current directory on the Python path.
To import `combustion1d` from a script elsewhere, add `numerics/` to `PYTHONPATH`. The
experiment scripts start several processes; `export OMP_NUM_THREADS=1` avoids
oversubscription (`reproduce_paper.sh` sets it).

## Quick start

```bash
python -m pytest -q                       # 39 tests, about 15 s
python -m examples.ex1_travelling_wave    # about 2 s
python -m examples.ex2_radial_vs_limit    # about 3 s
python -m examples.ex3_custom_reaction    # about 3 s
```

| example | what it shows |
|---|---|
| `ex1_travelling_wave.py` | ALE scheme started from the exact travelling wave $\Phi_\varepsilon$, with $\Delta t=10\,\varepsilon$: speed, error against the exact solution, cost |
| `ex2_radial_vs_limit.py` | a shrinking flame in 2D radial symmetry: $(P_\varepsilon)$ with the ALE scheme against $(P)$ with the front-fixing solver |
| `ex3_custom_reaction.py` | another reaction rate, $f(s)=((1-s)^+)^2$: travelling-wave speeds by shooting and with the ALE scheme |

The core of example 1:

```python
from combustion1d.ale import simulate_ale
from combustion1d.exact import tw_profile_eps, tw_speed_eps

A, w0, eps = 1.0, 0.0, 0.002
c_eps = tw_speed_eps(A, w0, eps)                  # exact speed of the eps-travelling wave
u0 = lambda x: tw_profile_eps(x - 2.0, A, w0, eps)
res = simulate_ale(u0, w0, eps, L=30.0, t_end=4.0, dt=0.02,
                   speed_guess=c_eps, exact_initial=True)
print(res.status, res.t[-1], res.front[-1])       # "ok", final time, front position
```

## The package `combustion1d`

| module | contents |
|---|---|
| `model.py` | the reaction rate `f`, `f_eps`, `M(w0)`, `fb_slope(w0)` $=\sqrt{2M(w_0)}$; `set_reaction` installs another $f$ |
| `solver.py` | band meshes (`MeshParams`, `build_mesh`, `find_front`), the finite volume Laplacian, the initial data of the paper (`initial_data`), the `Result` container, the exceptions, and the reference **BDF solver** `simulate` |
| `ale.py` | `simulate_ale`: the **ALE scheme**, i.e. the monotone linearly implicit scheme (S) of the paper on a mesh that moves with the flame |
| `simple.py` | `simulate_simple`: the **Eulerian scheme (S)**, the same scheme on band meshes that do not move |
| `exact.py` | travelling waves of $(P)$ and of $(P_\varepsilon)$, radial self-similar solution of $(P)$ |
| `limit.py` | `solve_radial`, `solve_planar`: the **limit problem** $(P)$ by front fixing |

### Meshes

All solvers use vertex-centred finite volumes on $[0,L]$ with zero flux at both ends; in
radial symmetry the origin needs no special treatment. The *front* is the point where
$u=\varepsilon$ (there is no reaction for $u\ge\varepsilon$). A **band mesh** has a band of uniform
spacing $h_f$ around the front, a geometric transition on each side, and coarse spacing
elsewhere:

| `MeshParams` field | default | meaning |
|---|---|---|
| `pts_per_eps` | 16 | $h_f=\varepsilon/$`pts_per_eps` in the band |
| `h_coarse` | 0.01 | $h_c$, mesh size away from the front |
| `K_fresh`, `K_burnt`, `K_margin` | 5, 20, 15 | extent of the band on each side of the front, plus a margin on both sides, in units of ε (total width 55ε) |
| `growth` | 1.1 | ratio of consecutive cells in the transition |

For a fixed ratio $h_f/\varepsilon$ the number of nodes grows only logarithmically in $1/\varepsilon$.

### The three solvers for $(P_\varepsilon)$

Common arguments: `u0`, a vectorized function of $x$ (or $r$); `w0`, a number or a
function of $x$ (then $w(\cdot,0)=\varepsilon w_0$); `eps`; `L`; `t_end`; `N` (1, 2 or 3);
`mesh` (a `MeshParams`); `speed_guess`; `snapshot_times`; `reaction`.

**`simulate_ale(u0, w0, eps, L, t_end, dt, ...)`**, the method of the paper. Each step is
linearly implicit (two tridiagonal solves) and satisfies a discrete comparison principle
for every Δt. The nodes move with velocity $v\,\chi(x)$, where $\chi=1$ on a *rigid zone*
that contains the band and extends `rigid_extra*ell0` beyond it, and $\chi$ decays to 0
over a further distance `ell0`; near the ends of $[0,L]$ the velocity is tapered to zero.
By default $v$ is *implicit*: at each step it is found by a secant iteration (normally
three trial steps) so that the front ends the step at the reference point of the band.
The mesh is rebuilt, on the lattice of the old band, when it has travelled
`remesh_distance`, when the front has drifted from the reference point, or when a cell
has been compressed below `min_cell`·$h_f$.

| option | default | meaning |
|---|---|---|
| `dt` | (required) | fixed time step; it need not be small with respect to ε |
| `speed_guess` | 1.0 | first trial for the mesh velocity; also sets the default `ell0` |
| `ell0`, `rigid_extra` | `max(2, 3/\|speed_guess\|)`, 1.0 | decay length of the mesh velocity and size of the rigid zone |
| `remesh_distance` | `ell0/4` | distance travelled by the mesh between remeshings |
| `exact_initial` | `False` | `False`: the initial data of the paper, $u_0^\varepsilon=u_0-\varepsilon w_0(1-u_0/\varepsilon)^+$, $w=\varepsilon w_0$. `True`: `u0` as given (for $\Phi_\varepsilon$; needs a constant `w0`) |
| `fixed_velocity` | `None` | prescribe the mesh velocity (a number or a function of $t$) |
| `implicit` | `True` | `False`: explicit estimate of the velocity from the measured front speed |
| `remesh`, `min_cell` | `True`, 0.5 | `remesh=False`: never rebuild the mesh; `min_cell=0`: do not rebuild compressed cells |
| `remesh_off_band` | `False` | also rebuild the mesh when the front has left the band |
| `record_mesh`, `monitor` | `False`, `None` | store the nodes of every step in `res.mesh_history`; call `monitor(t, x, u, v)` after every step |

**`simulate_simple(u0, w0, eps, L, t_end, N=1, dt=None, ...)`**: the same scheme on band
meshes that are rebuilt around the front and do not move. `dt` defaults to
`dt_factor*eps**2`. Its error is proportional to $\Delta t/\varepsilon$.

**`simulate(u0, w0, eps, L, t_end, N=1, ...)`**: method of lines on the same meshes with
SciPy's adaptive BDF and an analytic sparse Jacobian (`rtol=1e-5`, `atol=1e-4*eps`).
Accurate, but not covered by the analysis of the paper, and its number of steps grows
faster than $1/\varepsilon$. Both Eulerian solvers rebuild the mesh after every time interval
`min(7.5*eps/v_est, max_dt_remesh)`, where `v_est` is the measured speed of the front and
`max_dt_remesh = 0.05*t_end` by default.

All three return a `Result`:

| field | content |
|---|---|
| `t`, `front` | time levels and position of the front; `nan` once the flame is extinct |
| `snapshots` | `{t: (x, u, w)}` at the requested times |
| `status`, `reason`, `extinct` | see below |
| `n_steps`, `n_remesh`, `n_nodes` | number of steps and of remeshings; number of nodes after each step (ALE) or each remeshing |
| `n_solves`, `n_trials`; `n_lu`, `n_fev` | tridiagonal solves and trial velocities (scheme (S)); LU factorizations and right-hand sides (BDF) |
| `v`, `mesh_history` | ALE: mesh velocity of each step; `[(t, nodes)]` if `record_mesh=True` |
| `n_secant_fail`, `n_mesh_remesh`, `n_front_off_band` | ALE diagnostics, 0 in a well-resolved run: steps accepted without meeting the tolerance of the implicit velocity; remeshings forced by a compressed cell; steps that ended with the front outside the band |

### Statuses and exceptions

- `status == "ok"`: the run reached `t_end` with a front.
- `status == "extinct"` (and `extinct == True`): the front disappeared; `reason` says
  when. This is an outcome of the model, not an error, and nothing is raised. The ALE
  solver stops at that step; the Eulerian solvers continue on a uniform mesh.
- **Numerical failures raise** a subclass of `combustion1d.solver.NumericalFailure`
  (itself a `RuntimeError`): `PecletError` (ALE: no usable mesh velocity satisfies the
  cell Péclet condition $|\dot x_i|\,h\le2$, which keeps the scheme monotone),
  `MeshError` (ALE: nodes would cross), or `NumericalFailure` itself (non-finite values;
  the BDF integrator failed). The exception carries the partial result, with
  `status == "failed"`, in its attribute `result`:

```python
from combustion1d.ale import simulate_ale
from combustion1d.exact import tw_profile_limit
from combustion1d.solver import MeshParams, NumericalFailure

u0 = lambda x: tw_profile_limit(x, 0.0, 1.0, 0.577, x0=2.0)
try:                                    # a mesh velocity far too large for h_c = 0.05
    res = simulate_ale(u0, 0.0, 0.02, 12.0, 0.1, 0.01, mesh=MeshParams(h_coarse=0.05),
                       fixed_velocity=60.0)
except NumericalFailure as ex:
    res = ex.result
print(res.status, "|", res.reason)      # failed | step 1 (t = 0.01): PecletError: ...
```

- Invalid arguments (nonpositive `dt` or `ell0`, an initial datum without a front) raise
  `ValueError`.

Always look at `res.status` before using a result: the Eulerian solvers return
`status == "failed"`, without raising, if they find non-finite values when the front
disappears.

### Exact solutions and the limit problem

- `tw_speed_limit(A, w0)` $=\sqrt{2M}/A$ and `tw_profile_limit(x, t, A, c, x0)`
  $=A(1-e^{-c(x-x_0-ct)})^+$: the travelling wave of $(P)$.
- `tw_speed_eps(A, w0, eps)`: the speed $c_\varepsilon$ of the travelling wave of $(P_\varepsilon)$ with
  $w=\varepsilon w_0$, by shooting; `tw_profile_eps(z, A, w0, eps)`: its profile $\Phi_\varepsilon$,
  normalized by $\Phi_\varepsilon(0)=\varepsilon$.
- `SelfSimilar(N, w0, T)`: the radial self-similar solution of $(P)$ that becomes extinct
  at time `T`; `V(r, t)` evaluates it, `V.radius(t)` is the radius $R\sqrt{T-t}$ of its
  support, `V.R` and `V.kappa` its constants.
- `solve_radial(u0, w0, N, s0, t_end, n=400, t_eval=None)` solves $(P)$ for a fresh
  region $\{r<s(t)\}$, $s(0)=s_0$ (it stops when $s<0.05\,s_0$);
  `solve_planar(u0, w0, L, s0, t_end, n=800, t_eval=None)` for a fresh region
  $\{s(t)<x<L\}$. Second order in space, BDF in time. `w0` is a number or a function
  `w0(x, t)`, evaluated at the front. The result has `t`, `s` (the free boundary) and
  `u_on(x, k)` (the solution at time `t[k]`).

## Reproducing the paper

```bash
./reproduce_paper.sh          # tests, all experiments, all figures: about 2 h 10 min
./reproduce_paper.sh plots    # only redraw the figures from the stored results: about 1 min
```

The data are written to `results/` and `results/paper/` (JSON) and the figures to
`figures/` and `figures/paper/` (PNG). Every JSON file records the commit of the code and
the versions of Python, NumPy and SciPy (key `provenance`).

**Times** are those of the last full run, on a 4-core machine; they are approximate and
machine dependent. That run took 92 min for the script without its steps `variable_w0`
(9 min) and `work_precision` (26 min), which were run separately. Checked when this file
was last revised: `./reproduce_paper.sh plots` (1 min 16 s) reproduces the stored figures
byte for byte, and the scripts marked ✓ below reproduce the stored numbers exactly; the
other scripts, and the full `./reproduce_paper.sh`, were not run again.

Numbers of figures, tables and sections are those of the draft in `../paper/` and may
still change; the output file identifies the item. Each line is run with
`python -m <script>`.

| item | output | script | time | |
|---|---|---|---|---|
| Table 1 (flame speed, first-order correction) | `results/paper/speed_coefficient.json` | `experiments.paper.speed_coefficient` | 5 s | ✓ |
| Eq. (37) (averaged truncation error) | `results/paper/truncation_check.json` | `experiments.truncation_check` | 10 s | ✓ |
| Figure 1 | `figures/paper/A5_mesh_snapshots.png` | `experiments.paper.mesh_snapshots` | 5 s | ✓ |
| Section 6.1, nodes and remeshings | `results/paper/A_meshes.json` | `experiments.paper.meshes` | 3 min | |
| Table 2 | `results/limit_comparison.json` | `experiments.limit_comparison` | 8 min | |
| Figures 2 and 6 | `figures/paper/C1_ale_vs_limit.png`, `C2_self_similar_convergence.png` | `experiments.paper.limit_examples` | 13 min | |
| Figure 3 | `figures/paper/I1_variable_w0.png` | `experiments.paper.variable_w0` | 9 min | |
| Figure 4 | `figures/tw_instability.png` | `experiments.traveling_wave`, then `experiments.plots` | 21 min | |
| Section 6.5, setting of the convergence theorem | `results/paper/H_convergence_thm.json` | `experiments.paper.convergence_thm` | 25 s | ✓ |
| Figure 5 | `figures/paper/B4_radial_extinction_dt.png` | `experiments.paper.convergence` | 3 min | |
| Sections 6.6 and S1, uniform meshes inside the hypotheses | `results/paper/G_*.json` | `experiments.paper.in_hypotheses` | 45 s | ✓ |
| Figures S1, S2, S3; drag ratios | `figures/paper/D1_ale_error_vs_dt.png`, `D3_phase_pinning.png`, `D2_ale_error_vs_mesh.png` | `experiments.paper.tw_error` | 8 min | |
| Figure S4; Eulerian lag | `figures/paper/D6_transients.png` | `experiments.paper.transient` | 33 min | |
| Figure S5, Table S1 | `figures/paper/J1_work_precision.png`, `results/paper/J_work_precision.json` | `experiments.paper.work_precision` | 26 min | |

[`experiments/paper/README.md`](experiments/paper/README.md) is the reference for these
scripts: parameters, contents of each JSON file, how to run only a part or redraw only a
figure, and the experiments in the repository that are not used in the paper.

## Writing your own experiment

Choose the data ($u_0$, $w_0$, ε), the geometry ($N$, $L$), a solver, and what to
measure from the `Result`; the three examples are templates.

**Data.** `u0` is a vectorized function with $u_0\ge0$; the fresh region is $\{u_0>0\}$
and the code assumes **a single front** in $[0,L]$. To compare with $(P)$, take data with
$|\nabla u_0|=\sqrt{2M(w_0)}$ on the free boundary (`fb_slope(w0)`). In radial symmetry
the burnt region can be inside (an expanding flame) or outside (a shrinking one). For a
non-constant $w_0$ pass a function of $x$ to the solver, and the corresponding caloric
function `w0(x, t)` to the limit solver (`experiments/paper/variable_w0.py`).

**Parameters of the ALE scheme.**

- *Band.* The error is of second order in $h_f/\varepsilon$ and does not depend on ε for a fixed
  ratio; with the default `pts_per_eps=16` the error of the speed of a travelling wave is
  $10^{-4}$ to $7\cdot10^{-4}$ (supplement, Section S1). The coarse mesh adds $O(h_c^2)$.
- *Time step.* For a travelling wave the error does not depend on Δt (tested up to
  Δt = 100ε). For general data it is of first order in Δt, with a constant that does not
  grow as ε → 0 (between 1.4 and 1.9 in the extinction test of Section 6.5): choose Δt
  for the accuracy you need in the position of the front, not according to ε.
- *Rigid zone.* It must contain several decay lengths $1/c$ of the profile ahead of the
  front: pass `speed_guess` ≈ the expected speed, and for fast fronts ($c\gtrsim1.5$) use
  `ell0=4`.
- *Accelerating fronts and small domains.* Use a smaller `ell0` and `remesh_distance`
  (`ell0=1`, `remesh_distance=0.1` in example 2), and smaller time steps when the front
  gets close to $r=0$.
- *`PecletError`.* Reduce `h_coarse` or `ell0`, or check `speed_guess`.

**Measuring.** `res.front` is the level set $u=\varepsilon$, at distance
$\approx\varepsilon/\sqrt{2M}$ from the free boundary of $(P)$: add this shift before comparing
(example 2). Speeds are least-squares slopes of `res.front`. For profiles use
`snapshot_times`; to record something at every step, pass a `monitor`.

<a id="another-reaction-rate"></a>
**Another reaction rate.** `model.set_reaction(f, df=None, M=None)` installs a vectorized
$f$ for every solver and exact solution of the package (`df` is used by the Jacobian of
the BDF solver, `M` is computed by quadrature if omitted); `model.set_reaction()`
restores the default. The theory and the monotone schemes need $f$ Lipschitz and
nonincreasing, positive on $(-\infty,1)$ and zero on $[1,\infty)$. `set_reaction` changes
the state of the current process only: it is not thread safe and it does not reach
worker processes started with the `spawn` method (the default on macOS and Windows). For
parallel runs use `experiments.util.worker_pool`, which installs the rate of the parent
in every worker (the functions must then be defined at module level, not as lambdas):

```python
from experiments.util import worker_pool      # ProcessPoolExecutor + initializer
with worker_pool(4) as ex:
    results = list(ex.map(job, jobs))
```

Alternatively the three solvers accept `reaction=f` (or `reaction=(f, df)` for the BDF
solver) for one call, without touching the module state; the exact solutions and `M`
always use the installed rate.

## Tests

```bash
python -m pytest -q
```

39 tests in about 15 s (`tests/test_combustion1d.py`). They cover the exact solutions,
the mesh generator and the transfer at remeshing, the three solvers against the exact
travelling wave, the bounds and the discrete comparison principle of the scheme, the
failure modes (Péclet condition, extinction, compressed cells, front outside the band),
the limit-problem solvers (second order), custom reaction rates (including `spawn`
workers), and one stored result of `results/paper/B_convergence.json`.

## Known limitations

- **Single front, one dimension or radial symmetry** only.
- **Near $r=0$** (on $r<\min(0.5,\,0.1L)$) the ALE mesh velocity is tapered to zero, so
  that no node crosses the origin. A collapsing front that enters this zone is computed
  on an almost fixed mesh, and a lag $O(\Delta t\,v/\varepsilon)$ appears (about 1.5 % in
  $s/\sqrt{T-t}$ for ε = 0.002 in Section 6.7, `C2_self_similar_convergence.png`). The
  band nodes also move more slowly than the front there; for small ε the front can leave
  the band and stall in cells that do not resolve the layer. The run counts these steps
  in `res.n_front_off_band`; `remesh_off_band=True` rebuilds the mesh when it happens.
- **Rigid zone for fast fronts.** With the default `ell0` and $c\approx1.6$ an error of
  about $0.01\,\Delta t$ remains in the speed (supplement, Section S1); `ell0=4` removes it.
- **Lattice pinning.** The front keeps its position inside a band cell and its speed
  depends on that position, with an amplitude proportional to $(h_f/\varepsilon)^2$ (supplement,
  Section S1, `D3_phase_pinning.png`).
- **A mesh that translates as a whole** (`ell0=1e6`) is usable only for moderate travel
  distances, because the cells at the far end are compressed; the experiments that do
  this pass `min_cell=0` and limit the travel distance.
- **Sequence of meshes of the Eulerian solvers.** With the default `max_dt_remesh` it
  depends on `t_end`; pass `max_dt_remesh` to make runs with different final times
  comparable.
- **Reaction rate.** The rate installed with `set_reaction` is state of one process (see
  [above](#another-reaction-rate)).

## Citing

If you use this code, please cite the article above (see `../CITATION.cff`). The analysis
of the singular limit builds on J. Fernández Bonder and N. Wolanski, *Interfaces Free
Bound.* 2 (2000) 381–411, and *SIAM J. Math. Anal.* 36 (2004) 172–185.
