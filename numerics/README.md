# combustion1d: numerics for the high activation energy limit

Python code for the paper

> J. Fernández Bonder, *The high activation energy limit of an equidiffusional flame
> model: rate of convergence and a front-following numerical method* (preprint, 2026).
> Source in [`../paper/main.tex`](../paper/main.tex).

The code solves the flame model

$$T_t=\Delta T+Y f_\varepsilon(T_f-T),\qquad Y_t=\Delta Y-Y f_\varepsilon(T_f-T),\qquad f_\varepsilon(s)=\varepsilon^{-2}f(s/\varepsilon),$$

written, with $u=T_f-T$ and $w=Y-u$ (a caloric function), as

$$(P_\varepsilon)\qquad u_t=\Delta u-(u+w)\,f_\varepsilon(u),\qquad w_t=\Delta w,$$

in one space dimension ($N=1$) and in radial symmetry ($N=2,3$). It also solves the
limit free boundary problem obtained when $w/\varepsilon\to w_0$:

$$(P)\qquad u_t=\Delta u\ \text{ in }\{u>0\},\qquad |\nabla u|=\sqrt{2M(w_0)}\ \text{ on }\partial\{u>0\},\qquad M(w_0)=\int_{-w_0}^1(s+w_0)f(s)\,ds .$$

The default reaction rate is $f(s)=(1-s)^+$, for which $M(w_0)=(1+w_0)^3/6$; any other
rate can be installed (see [Another reaction rate](#another-reaction-rate)).

Contents:

1. [Installation](#installation)
2. [Quick start](#quick-start)
3. [The package `combustion1d`](#the-package-combustion1d)
4. [Reproducing the paper](#reproducing-the-paper)
5. [Writing your own experiments](#writing-your-own-experiments)
6. [Tests](#tests)
7. [Known limitations](#known-limitations)

---

## Installation

Python ≥ 3.10 with `numpy`, `scipy`, `matplotlib` and `pytest`:

```bash
python -m venv .venv
```

```bash
.venv/bin/pip install -r numerics/requirements.txt
```

All commands below are run **from the `numerics/` directory**, with that directory on the
Python path (`pytest.ini` and `reproduce_paper.sh` take care of it):

```bash
cd numerics
```

```bash
export PYTHONPATH=$PWD OMP_NUM_THREADS=1
```

`OMP_NUM_THREADS=1` avoids oversubscription: the experiment scripts run several
cases in parallel processes.

## Quick start

Three short examples are in [`examples/`](examples). Each one runs in a few seconds:

```bash
python -m examples.ex1_travelling_wave
```

```bash
python -m examples.ex2_radial_vs_limit
```

```bash
python -m examples.ex3_custom_reaction
```

| example | what it shows |
|---|---|
| `ex1_travelling_wave.py` | ALE scheme started from the exact travelling wave $\Phi_\varepsilon$; error against the exact solution, speed, cost |
| `ex2_radial_vs_limit.py` | a shrinking flame in 2D radial symmetry: $(P_\varepsilon)$ with the ALE scheme against $(P)$ with the front-fixing solver |
| `ex3_custom_reaction.py` | installing another reaction rate, $f(s)=((1-s)^+)^2$, and computing travelling-wave speeds |

The core of example 1:

```python
from combustion1d.ale import simulate_ale
from combustion1d.exact import tw_profile_eps, tw_speed_eps

A, w0, eps = 1.0, 0.0, 0.002
c_eps = tw_speed_eps(A, w0, eps)                  # exact speed of the eps-travelling wave
u0 = lambda x: tw_profile_eps(x - 2.0, A, w0, eps)
res = simulate_ale(u0, w0, eps, L=30.0, t_end=4.0, dt=0.02,
                   speed_guess=c_eps, exact_initial=True)
print(res.t[-1], res.front[-1])                   # final time and front position
```

---

## The package `combustion1d`

| module | contents |
|---|---|
| `model.py` | the reaction rate `f`, `f_eps`, their derivatives, `M(w0)`, `fb_slope(w0)` $=\sqrt{2M(w_0)}$, and `set_reaction` to install another $f$ |
| `exact.py` | exact solutions: limit travelling wave (`tw_speed_limit`, `tw_profile_limit`), travelling wave of $(P_\varepsilon)$ (`tw_speed_eps`, `tw_profile_eps`, by ODE shooting), and the radial self-similar solution of $(P)$ (`SelfSimilar`, via Kummer's function) |
| `solver.py` | band meshes (`MeshParams`, `build_mesh`, `find_front`), the finite volume Laplacian (`laplacian`), the initial data of the paper (`initial_data`), the `Result` container, and the reference **BDF solver** (`simulate`) |
| `ale.py` | the **ALE scheme** (`simulate_ale`): the monotone linearly implicit scheme (S) on a mesh that moves with the flame |
| `simple.py` | the **Eulerian scheme** (S) (`simulate_simple`): the same scheme on band meshes rebuilt around the front, fixed time step |
| `limit.py` | solvers for the **limit problem** (P) by front fixing: `solve_radial`, `solve_planar` |

### Spatial discretization and meshes

All solvers use vertex-centred finite volumes on non-uniform meshes. In radial
symmetry the flux weight $r^{N-1}$ vanishes at $r=0$, so the origin needs no special
treatment. Zero-flux conditions are imposed at both ends of $[0,L]$.

A **band mesh** (`MeshParams`, `build_mesh`) has three zones:

- a band around the front with uniform spacing $h_f=\varepsilon/$`pts_per_eps`;
- a geometric transition with ratio `growth`;
- coarse spacing `h_coarse` elsewhere.

The band covers `K_burnt`·ε on the burnt side and `K_fresh`·ε on the fresh side,
plus `K_margin`·ε on each side. The front is the point where $u=\varepsilon$; the
reaction vanishes for $u\ge\varepsilon$.

| `MeshParams` field | default | meaning |
|---|---|---|
| `pts_per_eps` | 16 | $h_f=\varepsilon/16$ in the band |
| `h_coarse` | 0.01 | $h_c$, mesh size away from the front |
| `K_fresh`, `K_burnt`, `K_margin` | 5, 20, 15 | band extent in units of ε (total width 55ε) |
| `growth` | 1.1 | ratio of consecutive cells in the transition |

The number of nodes is $O(L/h_c)+O(\text{band}/h_f)+O(\log(h_c/h_f))$; for a fixed ratio
$h_f/\varepsilon$ it grows only logarithmically in $1/\varepsilon$.

### The three solvers for $(P_\varepsilon)$

**`simulate_ale(u0, w0, eps, L, t_end, dt, ...)`: ALE scheme (the method of the paper).**

- *Time stepping.* The scheme is linearly implicit: the reaction rate is lagged and the
  reactant is implicit, $(Y^{n+1}-Y^n)/\Delta t=\Delta_hY^{n+1}-f_\varepsilon(U^n)Y^{n+1}$, and the
  same for $U$. Each step costs two tridiagonal solves, and the scheme satisfies a
  discrete comparison principle for every Δt.
- *Moving mesh.* The nodes move with velocity $v\chi(x)$. Here $\chi=1$ on the band and
  on a distance `rigid_extra*ell0` on each side of it, and then decays linearly to 0 over
  `ell0`, with `ell0` $=\max(2,3/c)$ by default. The equation is written along the node
  trajectories, which adds a convection term. The code checks that the matrices stay
  M-matrices (cell Péclet number $|v|h\le2$).
- *Implicit mesh velocity.* At each step, $v$ is chosen by a secant iteration so that the
  computed front ends the step at the reference point of the band. This takes three trial
  velocities, i.e. six tridiagonal solves per step.
- *Remeshing.* The mesh is rebuilt when it has travelled `remesh_distance` (default
  `ell0/4`). The new band lies on the lattice of the old band, so the transfer by linear
  interpolation is exact there.

Main options:

| option | meaning |
|---|---|
| `N` | 1 (planar), 2 or 3 (radial) |
| `mesh` | a `MeshParams` instance |
| `speed_guess` | initial guess for the mesh velocity; also sets the default `ell0` |
| `ell0`, `rigid_extra` | size of the rigid zone and of the decay zone (see the limitations below) |
| `remesh_distance` | travel distance between remeshings |
| `exact_initial` | `True`: use `u0` as given and $w=\varepsilon w_0$. `False` (default): the initial data of the paper, $u_0^\varepsilon=u_0-\varepsilon w_0(1-u_0/\varepsilon)^+$, which have $Y_0^\varepsilon=0$ in the burnt region |
| `snapshot_times` | times at which $(x,u,w)$ is stored in `res.snapshots` |
| `fixed_velocity` | prescribe the mesh velocity instead of computing it (rigid-mesh tests) |
| `record_mesh` | store the node positions at every step in `res.mesh_history` |
| `monitor` | callable `monitor(t, x, u, v)` called after every step |

**`simulate_simple(u0, w0, eps, L, t_end, N=1, dt=None, ...)`: Eulerian scheme (S).**
Same time stepping; the band mesh is rebuilt around the front at regular intervals.
`dt` defaults to `dt_factor*eps**2`. The error of this scheme grows like
$\Delta t/\varepsilon$.

**`simulate(u0, w0, eps, L, t_end, N=1, ...)`: reference BDF solver.** Adaptive BDF
(scipy) with an analytic sparse Jacobian, on band meshes rebuilt around the front;
the solution is transferred by monotone cubic interpolation. It is accurate but not
covered by the analysis of the paper, and its number of steps grows like
$\varepsilon^{-1.4}$.

All three solvers return a `Result` with:

| field | content |
|---|---|
| `t`, `front` | time levels and front positions ($u=\varepsilon$); `nan` once the flame is extinct |
| `snapshots` | `{t: (x, u, w)}` at the requested times |
| `extinct` | `True` if the front disappeared (ALE stops there) |
| `n_steps`, `n_remesh`, `n_nodes` | cost counters |
| `n_solves`, `n_trials` | tridiagonal solves and trial velocities (scheme (S)) |
| `n_lu`, `n_fev` | LU factorizations and RHS evaluations (BDF) |
| `v` | ALE mesh velocity of each step |
| `mesh_history` | `[(t, nodes)]` if `record_mesh=True` |

The arguments are:

- `u0`: a vectorized function of $x$ (or $r$);
- `w0`: a number, the limit of $w/\varepsilon$ (`simulate` also accepts a function);
- `L`: the length of the domain $[0,L]$.

### Exact solutions and the limit problem

- **Limit travelling wave.** `tw_speed_limit(A, w0)` $=\sqrt{2M}/A$, and
  `tw_profile_limit(x, t, A, c, x0)` $=A(1-e^{-c(x-x_0-ct)})^+$.
- **Travelling wave of $(P_\varepsilon)$.** `tw_speed_eps(A, w0, eps)` is the exact speed
  $c_\varepsilon$ of the travelling wave with $w=\varepsilon w_0$, computed by shooting with bisection
  on the energy. `tw_profile_eps(z, A, w0, eps)` is the profile, normalized by
  $\Phi_\varepsilon(0)=\varepsilon$.
- **Self-similar solution.** `SelfSimilar(N, w0, T)` is the radial self-similar solution of
  $(P)$ that becomes extinct at time `T`; `V(r, t)` evaluates it and `V.radius(t)` gives
  the radius $R\sqrt{T-t}$ of its support.
- **Limit problem.** `solve_radial(u0, w0, N, s0, t_end, n=400, t_eval=None)` solves $(P)$
  for a radial support $\{r<s(t)\}$ with $s(0)=s_0$. `solve_planar(u0, w0, L, s0, t_end, ...)`
  solves it for a planar support $\{s(t)<x<L\}$. Both use front fixing (Landau
  transformation), are second order, and integrate in time with BDF. The result has `t`,
  `s` (the free boundary) and `u_on(x, k)` (the solution at time `t[k]`).

---

## Reproducing the paper

One command reruns every experiment and redraws every figure (about one hour on
4 cores):

```bash
./reproduce_paper.sh
```

To only redraw the figures from the stored results (about one minute):

```bash
./reproduce_paper.sh plots
```

The results are stored as JSON in `results/` and `results/paper/`, and the figures as
PNG in `figures/` and `figures/paper/`. Detailed captions, every number quoted in the
paper, and the history of the changes made after the referee report are in
[`experiments/paper/README.md`](experiments/paper/README.md).

| in the paper | script (run with `python -m ...`) | solver | time |
|---|---|---|---|
| Table `tab:speedcoef` (Prop. 2.1) | `experiments.paper.speed_coefficient` | shooting | 10 s |
| Section 5.5, averaged truncation error (`an:eq:solv`) | `experiments.truncation_check` | — | 15 s |
| Fig. `fig:A5` (solution and mesh at successive times) | `experiments.paper.mesh_snapshots` | ALE | 10 s |
| Fig. `fig:A1` and the cost figures of Section 6.1 (node trajectories, nodes per zone, steps) | `experiments.paper.meshes` | ALE, BDF | 4 min |
| Figs. `fig:D1`, `fig:D2`, `fig:D3`; drag ratios (error from $\Phi_\varepsilon$, lattice pinning) | `experiments.paper.tw_error` | ALE, Eulerian (S) | 6 min |
| Section 6.2, uniform mesh inside the hypotheses of the travelling wave theorem | `experiments.paper.in_hypotheses` | ALE (rigid mesh) | 1 min |
| Section 6.3, setting of the convergence theorem (fixed uniform mesh; band with prescribed velocity, no remeshing) | `experiments.paper.convergence_thm` | (S) on a fixed mesh, ALE with `fixed_velocity`, `remesh=False` | 1 min |
| Fig. `fig:D6`; accumulated mesh-velocity error $M_T$; errors over the lattice phase | `experiments.paper.transient` | ALE, Eulerian (S) | 30 min |
| Fig. `fig:B4` (radial extinction, error against Δt) | `experiments.paper.convergence` | ALE | 2 min |
| Figs. `fig:C1`, `fig:C2` (limit problem, self-similar extinction); `C3_quenching.png` (quenching of small kernels, an extra example that is not in the paper) | `experiments.paper.limit_examples` | ALE, BDF, limit | 12 min |
| Fig. `fig:instability` | `experiments.traveling_wave`, then `experiments.plots` | BDF | 20 min |
| Table `tab:limit` | `experiments.limit_comparison`, then `experiments.plots` | BDF, limit | 15 min |

Some scripts also produce figures and tables that were in earlier drafts and are not in
the paper (`A2`–`A4`, `B1`–`B3`, `D4`, `D5`, `tab_simple.py`); they are described in
[`experiments/paper/README.md`](experiments/paper/README.md).

The scripts in `experiments/paper/` accept `run` or `plot` as an argument; with no
argument they do both. Some also accept a list of groups, for example
`python -m experiments.paper.tw_error run D3 D4`. The node histories behind
Figures A1 and A4 (`results/paper/A_mesh_*.npz`) are not stored in git and are
regenerated by `experiments.paper.meshes run`.

Other scripts in `experiments/`, not used in the paper:

- `radial.py`: convergence to the self-similar solution with the BDF solver;
- `movie_radial.py`, `movie_expanding.py`: animations of a shrinking and an expanding
  flame (`figures/*.mp4`).

---

## Writing your own experiments

A new experiment has four ingredients:

1. the **data**: $u_0$, $w_0$ and $\varepsilon$;
2. the **geometry**: $N$ and the domain $[0,L]$;
3. a **solver** with its parameters;
4. what you **measure** from the `Result`.

The examples in `examples/` are templates.

**Initial data.** `u0` is any vectorized function of $x$ (or $r$) with $u_0\ge0$; the
fresh region is $\{u_0>0\}$. If the flame is to be compared with the limit problem, take
data compatible with the free boundary condition: $|\nabla u_0|=\sqrt{2M(w_0)}$ on
$\partial\{u_0>0\}$ (use `fb_slope(w0)`). The code assumes **a single front**: one
interface in $[0,L]$. In radial symmetry the burnt region can be either inside
(an expanding flame, as in `experiments/paper/common.py: expanding_u0`) or outside
(a shrinking flame). For the travelling wave of $(P_\varepsilon)$ itself, use
`tw_profile_eps` with `exact_initial=True`.

**Choosing the parameters of the ALE scheme.**

- **Mesh sizes.** The error is of second order in $h_f/\varepsilon$ and essentially
  independent of $\varepsilon$ for a fixed ratio; `pts_per_eps=16` gives errors of order
  $10^{-4}$ in the speed. The coarse size `h_coarse` contributes $O(h_c^2)$.
- **Time step.** For a travelling wave the error does not depend on Δt (tested up to
  Δt = 100ε). For general data the error is $O(\Delta t)$ with a constant of order 1 to 2 that
  does not grow as $\varepsilon\to0$; choose Δt for the accuracy you need in the front
  position, not according to ε.
- **Rigid zone.** It must contain several decay lengths $1/c$ of the profile ahead of the
  front. Pass `speed_guess` ≈ the expected speed (it sets `ell0`=max(2, 3/c)); for fast
  fronts ($c\gtrsim1.5$) use `ell0=4` or larger.
- **Accelerating fronts** (extinction) and **small domains.** Use a smaller `ell0` and
  `remesh_distance`, as in `examples/ex2_radial_vs_limit.py` (`ell0=1`,
  `remesh_distance=0.1`). Near $r=0$ the mesh velocity is tapered to zero (see the
  limitations), so take smaller time steps when the front gets close to the origin.
- **Péclet bound.** A `RuntimeError: cell Peclet number too large` means $|v|h>2$ on
  some cell. Reduce `h_coarse` or `ell0`, or check `speed_guess`.

**Measuring.**

- *Front.* `res.front` is the position of the level set $u=\varepsilon$. It lies at a
  distance $\approx\varepsilon/\sqrt{2M}$ from the free boundary of the limit problem; add this
  shift before comparing with $(P)$, as in example 2.
- *Speeds.* Least-squares slopes of `res.front` over a time window.
- *Profiles.* Use `snapshot_times`. To record anything at every step (errors against an
  exact solution, the mesh, the mesh velocity) without storing all the states, pass a
  `monitor(t, x, u, v)` callback.

<a id="another-reaction-rate"></a>
**Another reaction rate.** `model.set_reaction(f, df=None, M=None)` installs a
vectorized $f$, which every solver and exact solution in the package then uses:

- `df` is used by the Jacobian of the BDF solver; if it is omitted, a finite difference is
  used.
- `M(w0)` is computed by quadrature if it is not given.
- `model.set_reaction()` with no arguments restores $(1-s)^+$.

The theory and the monotone schemes need $f$ Lipschitz and **nonincreasing**, positive on
$(-\infty,1)$ and zero on $[1,\infty)$. Example 3 does this for
$f(s)=((1-s)^+)^2$. Some scripts in `experiments/` and the self-similar solution assume
the default $f$ only through `M`, which they read from `model`.

**Non-constant $w_0$.** `solver.simulate` accepts `w0` as a function of $x$; then $w$
solves the heat equation with $w(\cdot,0)=\varepsilon w_0$. The ALE and Eulerian (S) schemes
take a constant $w_0$.

---

## Tests

```bash
python -m pytest -q
```

The tests (`tests/test_combustion1d.py`, about 10 s) check:

- the closed form of $M$ and the finite volume Laplacian;
- the exact solutions (including that the self-similar profile solves the PDE);
- the travelling-wave speed by shooting and the profile $\Phi_\varepsilon$ against its ODE;
- the band meshes;
- the BDF solver and the ALE scheme against the exact travelling wave (the ALE error does not change with Δt);
- the limit-problem solvers against the exact solutions (second order);
- a custom reaction rate.

## Known limitations

- **Single front, one dimension or radial symmetry.** Several interfaces and genuinely
  2D/3D problems are not implemented.
- **Near $r=0$** the ALE mesh velocity is tapered to zero so that no node crosses the
  origin. When a collapsing front enters this zone the scheme becomes partly Eulerian, and
  an $O(\Delta t\,v/\varepsilon)$ lag appears (about 1.5% in $s/\sqrt{T-t}$ at ε = 0.002 near
  extinction, Figure C2).
- **Rigid zone for fast fronts.** With the default `ell0` and $c\approx1.6$ a term
  $\approx0.01\,\Delta t$ remains in the speed; `ell0=4` removes it.
- **Lattice pinning.** With the implicit mesh velocity the front keeps its position
  inside a band cell, and its speed depends on that position by $O((h_f/\varepsilon)^2)$, with
  zero mean (Figure D3).
- **The rigid mode** (`fixed_velocity`) is only usable for moderate travel distances.

## Citing

If you use this code, please cite the paper above. The analysis of the singular limit
builds on J. Fernández Bonder and N. Wolanski, *Interfaces Free Bound.* 2 (2000) 381–411,
and *SIAM J. Math. Anal.* 36 (2004) 172–185.
