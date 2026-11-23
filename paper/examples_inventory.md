# Inventory of explicit examples in the reference papers

This file lists the concrete examples in the PDFs in `paper/pdfs/`: explicit solutions, special solutions, counterexamples, asymptotic profiles and numerical illustrations. For each one it gives where it is, the formula, the claimed behaviour, the normalization the paper uses, and whether our code can reproduce it.

Conventions used below:

* "p." is the printed page number. "pdf p." is the page index in the PDF file.
* Every formula is copied from the paper, from the text layer or, where that was garbled, from a rendered image of the page.
* Lines marked **[ours]** are not in the paper. They are our own remarks or numerical values, computed with `scipy` in the project venv. Each one is labelled so it can be checked or discarded.
* Verdict codes:
  * **(a)** ε-problem, 1D planar or radial, one front;
  * **(b)** direct limit solver, 1D planar or radial, one front. `combustion1d/limit.py` currently handles only two cases: radial with the fresh region a ball `{r < s(t)}`, and planar with the fresh region `{s < x < L}` and zero flux at `L`;
  * **(c)** an exact solution that is already implemented;
  * **"needs multiple fronts/annulus"**, **"needs 2D"**, or **"not a numerical target"**.

## Normalization (important for every comparison)

All the papers read here set the free boundary constant to 1. Most write `|∇u| = 1`, or `u_ν = -1`, or `-∂u/∂ν = 1`. The exceptions are Lederman–Vázquez–Wolanski and Vázquez 1996 §3, which use `√(2M)`.

* In the papers that use 1, the ε-problem is `u_t = Δu − β_ε(u)`, with `β_ε(s) = ε^{-1} β(s/ε)`, `supp β = [0,1]` and `∫β = 1/2`. Then `|∇u| = √(2·½) = 1`.
* Our problem fits this form. Write `(u+εw0) f_ε(u) = β_ε(u)` with `β(σ) = (σ+w0) f(σ)`. Then `∫β = M`, and the free boundary constant is `√(2M)`. **[ours]**
* **Scaling rule [ours]:** the problem is linear in `u` and both free boundary conditions are homogeneous. So if `u` solves a paper's problem with `|∇u| = 1`, then `√(2M)·u` solves ours with `|∇u| = √(2M)`.
  * The free boundary, the extinction time and every similarity variable `η` stay the same.
  * Only the amplitudes (the values of `u`) are multiplied by `√(2M)`.
  * Mass identities scale the same way. For example, `d/dt ∫u = −√(2M)·H^{N−1}(Γ(t))`.

---

## 1. Caffarelli & Vázquez, Trans. AMS 347 (1995), 411–441 — `CaffarelliVazquez1995.pdf`

**Normalization.**
* Free boundary conditions (0.2): `u = 0`, `|∇u| = c`. The paper then sets `c = 1` ("without loss of generality").
* (0.7): `β_ε(s) = ε^{-1}β(s/ε)`.
* Assumptions on β, p. 413: (i) `β > 0` exactly on (0,1); (ii) smooth; (iii) increasing on (0,½), decreasing on (½,1); (iv) `∫β = 1/2`.
* In this paper `u` is minus the temperature: `u = λ(T_c − T)` (p. 412). This sign convention holds in all the papers read here.

### CV-1. Radial self-similar extinction solution (§1, eqs. (1.1)–(1.8), pp. 415–416)

**Formula.**
* (1.1): `U(x,t) = (T−t)^α f(x/(T−t)^β)`. Plugging in forces `α = β = 1/2`.
* (1.2): `Δf = −½ f + ½ ξ·∇f`.
* Radial problem (E):
  * (1.3a) `f'' + ((n−1)/r − r/2) f' + ½ f = 0` for `0 < r < R`;
  * (1.3b) `f'(0) = 0`, `f > 0` on `[0, R)`;
  * (1.3c) `f(R) = 0`, `f'(R) = −1`.
* Proposition 1.1: problem (E) has a unique solution. The proof shoots from `f(0) = 1`, `f'(0) = 0`, finds the first zero `R`, and divides by `c = −f'(R)`.

**Remarks in the paper (p. 416).**
1. In 1D, (1.6): `f''(r) = −(a/2) e^{r²/4}`, with `f(0) = a`, `f'(0) = 0` and `a` chosen properly. Also `f'' ≤ −a/2`.
2. `sup U(·,0) = a T^{1/2}`.
   * (1.7): the support is `{0 ≤ t ≤ T, |x| ≤ R (T−t)^{1/2}}`.
   * (1.8): the free boundary is `|x| = R (T−t)^{1/2}`.
   * The solution is only `C^{1/2}` in t at the extinction point `(0, T)`. It decays like `O((T−t)^{1/2})`.
3. `|f'| < 1` and `|∇U| < 1` inside the support.
   * At the extinction point the nontangential limit of `|∇U|` is 0.
   * The limit along the free surface is 1.

**Behaviour (Proposition 1.2, p. 416).** `U` is the unique limit of (P_ε) from the data `U(x,0;T)`. This holds for any admissible `β_ε` and any approximating data `u0ε` with `u0ε → u0` uniformly and `supp u0ε → supp u0`.
* Corollary 1.4: the comparison principle holds between `U` and any limit solution.
* Proposition 1.5: `U` is stable.
* Theorem 1.6: compactly supported data have compact support for all t and vanish in finite time.

Rescaling (1.9)–(1.10): `U^λ(x,t) = λ^{-1} U(λx, λ²t)` has extinction time `λ²T`.

Inside the proof (Lemma 1.3, eq. (1.18), p. 418), the paper gives `μ0(ε) → μ* = 1/(2Mλ²)` and an estimate `R < R1 < R + Cε` for the sub/supersolution construction. This is a construction estimate, **not** a stated convergence rate for `u_ε`.

**Verdict.** Reproducible now with **(c)**, which is already implemented, and **(a)/(b)**.
* Possible ε-test: start from `U(x,0;T)` and measure `T_ε − T` and the front position. The paper proves convergence but gives no rate.

### CV-2. Planar travelling waves (§2.1, eqs. (2.1)–(2.6), pp. 420–421)

**Formula.**
* (2.1): `u = f(x1 − ct)`, `ξ = x1 − ct`.
* Advancing wave (2.2): `f_c^a(ξ) = (1/c){e^{−cξ} − 1}` on `{x1 < ct}`. The fresh domain expands, and (2.5) `f_c^a ≈ (1/c) e^{c|ξ|}` as `ξ → −∞`.
* Receding wave (2.3): `f_c^r(ξ) = (1/c){1 − e^{−cξ}}` on `{x1 > ct}`. The fresh domain shrinks, so the flame advances. Then (2.6) `f_c^r → 1/c` as `ξ → ∞`.
* Stationary waves (2.4): `[x1]_+` and `[−x1]_+` (the case `c = 0`).

**Behaviour.** "The methods of Section 1 can be used to show the uniqueness and stability of these solutions" (p. 421). No proof is given.

**Normalization `|∇u| = 1`.** [ours] In our normalization the receding wave is `√(2M)(1/c)(1 − e^{−cξ})`, so `A = √(2M)/c`. This is exactly `exact.py` (`c0 = √(2M)/A`).

**Verdict.**
* Receding wave: reproducible now with **(c)** (implemented) and **(a)/(b)**.
* Advancing wave: unbounded, exponentially growing data. It would need Dirichlet data at the far boundary, which (b) does not have (zero flux at L). It is linearly unstable in practice. Low priority.

### CV-3. Stationary solutions with a hole (§2.2, eqs. (2.7)–(2.10), pp. 421–422)

**Formula.**
* (2.7): `Δu = 0` for `|x| > R`, `u(R) = 0`, `u_ν(R) = −1`, with `n > 1`.
* `n = 2`, (2.8): `u(r) = R log(r/R)`.
  * If `u(1) = A` is prescribed, then (2.9) gives `A = −R log R`.
  * There is no solution if `A > A* = 1/e`, one solution `R = 1/e` if `A = A*`, and two solutions (big hole and small hole) if `0 < A < A*`.
  * As `A → 0` the small hole radius tends to 0 and the large one tends to 1.
* `n ≥ 3`, (2.10): `u(r) = (R^{n−1}/(n−2)) {1/R^{n−2} − 1/r^{n−2}}`. This is bounded on `{|x| > R}`.
  * For `n = 3` it is "the famous Zel'dovich flame …, which is unstable" (p. 422). The instability is stated without proof.
  * [ours] For `n = 3`, (2.10) reads `u = R − R²/r`. Prescribing `u(1) = A` gives `A = R − R²`, so two holes exist for `0 < A < 1/4`. This is derived by us from (2.10), not stated in the paper.

**Verdict.** The fresh region is the *exterior* of a ball.
* **(b)** cannot do this now (radial = fresh ball only).
* **(a)**: `solver.find_front` returns the burnt side, so an exterior configuration may work. This is untested; it needs a Dirichlet condition at `r = L`.
* Good test of the instability of radial stationary flames (see also GHV Remark 2 below). Medium priority.

### CV-4. Exterior initial data (§8, eq. (8.13), p. 433)

**Formula.** "Another interesting form is (for n ≥ 3)", (8.13): `u0(x) = A − b/|x|^{n−2}` for `|x| > a`.
* `a` is fixed by `u0(a) = 0`, so `A a^{n−2} = b`.
* The boundary condition forces `b ≤ a^{n−1}/(n−2)`.

It is given as typical data satisfying the paper's hypotheses (H1)–(H2) (`u_t ≤ 0`). The "typical example" is compactly supported concave data like the self-similar solutions.

**Verdict.** Data example only, with an exterior fresh region. Same status as CV-3.

### CV-5. Mass law (eq. (3.8), p. 423)

**Formula.** For classical solutions, `d/dt ∫_{Ω(t)} u dx = −λ_{n−1}(Γ(t))` (the area of the front). The paper says this also holds for the weak solutions of Theorem 1.6.

**Verdict.** Not an example. It is an exact diagnostic for every computation: with `√(2M)` the right side becomes `−√(2M)·|Γ(t)|`.

### CV-6. "Trivial" nonuniqueness of approximation (p. 413)

If `u0ε ≥ ε`, the absorption never acts and the limit is just the positive solution of the heat equation, with no free boundary.

**Verdict.** Not a numerical target. It explains why our initial data must contain a burnt region with `u0ε < ε`. LVW Remark 6.1 says the same.

---

## 2. Galaktionov, Hulshof & Vázquez, J. Math. Pures Appl. 76 (1997), 563–608 — `GalaktionovHulshofVazquez1997.pdf`

**Normalization.**
* (0.2): `u = 0`, `∂u/∂ν = −1`.
* In Part I, `u` is replaced by `−u`, giving (1.1): `u(ζ,t) = 0`, `u_r(ζ,t) = 1`.
* §2 (p. 574): `β_ε` is scaled from `β` with `∫β = 1/2`, "needed to obtain … slope precisely 1".
* The pdf page of printed p. `k` is `k − 562`.

### GHV-1. Rescaled problem and self-similar extinction profile f0 (§4.1–4.2, eqs. (4.1)–(4.8), pp. 580–581; §5, (5.1)–(5.3), p. 582)

**Problem.** (P), p. 580: radial, `u_t = r^{1−N}(r^{N−1}u_r)_r` for `0 < r < ζ(t)`, `u_r(0,t) = 0`, `u(ζ,t) = 0`, `u_r(ζ,t) = −1`.

**Self-similar solution.**
* (4.1): `u(r,t) = √(T−t) f0(η)`, `η = r/√(T−t)`.
* (4.2): `A(f0) = 0` on `(0, η0)`, `f0'(0) = 0`, `f0 > 0`, `f0(η0) = 0`, `f0'(η0) = −1`.
* (4.3)/(5.1): `A(f) ≡ η^{1−N}(η^{N−1} f')' − ½ η f' + ½ f`.

**Rescaled variables.**
* (4.4): `w(η,τ) = (T−t)^{−1/2} u(r,t)`, `τ = −log(T−t)`.
* (4.5): `ψ(τ) = e^{τ/2} ζ(t)`.
* (4.6): `w_τ = A(w)`, with `w = 0`, `w_η = −1` on the lateral boundary.

**Interface equations.**
* (4.7): `ζ'(t) = u_t(ζ,t) = u_rr(ζ,t) − (N−1)/ζ(t)`.
* (4.8): `ψ'(τ) = ½ψ + w_ηη − (N−1)/ψ`.

**Kummer series (5.2).** `f0(η) = c0[1 − η²/(4N) − … − (2k−1)!!/(4^{k+1}(k+1)!) · η^{2k+2}/(N(2+N)(4+N)…(2k+N)) − …]`.
* `f0` is decreasing and concave for `η > 0`.
* It has a unique positive root `η0(N) < 2√N`.
* (5.3): `f0(η) ≈ −a0 η^{−(N+1)} e^{η²/4}` for `η ≫ 1`.

[ours] (5.2) is `c0·M(−½, N/2, η²/4)`, which is the form used in `exact.py`.

**Mass equation (6.8), p. 586.** `dE/dτ = ((N+1)/2) E − ψ^{N−1}`, with `E(τ) = ∫_0^{ψ} η^{N−1} w dη`.
* The Remark on p. 586 says that `((N+1)/2)E − ψ^{N−1} → 0`.
* [ours] We checked numerically that `f0` satisfies `((N+1)/2)∫_0^{η0} η^{N−1} f0 = η0^{N−1}` to 1e-8 for N = 1, …, 5.

**Numerical constants [ours], computed with `scipy`.** `f0` is normalized by `f0'(η0) = −1`, i.e. `|∇u| = 1`. Multiply `f0(0)` by `√(2M)` for our normalization.

| N | η0 (rescaled front radius) | f0(0) = lim u(0,t)/√(T−t) |
|---|---|---|
| 1 | 1.8482777460 | 0.7868021977 |
| 2 | 2.5139266836 | 0.9817148738 |
| 3 | 3.0039505365 | 1.1053303017 |
| 4 | 3.4084194168 | 1.1981718577 |
| 5 | 3.7602307821 | 1.2735715438 |

For N = 1, solving the paper's (7.5) gives the same `η0`, and the paper's formula for `c0` gives the same `f0(0) = c0`.

**Verdict.** This is the same profile as **(c)**, already implemented.

### GHV-2. Asymptotic convergence: single-point extinction for general radial ball data (§6, Theorem 6.6, p. 588; Propositions 6.1–6.2, Corollary 6.3, pp. 584–586)

**Setting.** Bounded radial data `u0 ≥ 0` supported in a ball (§4, p. 579), with N > 1 (§4) and regular, eventually monotone data.

**Result.** Theorem 6.6: as `τ → ∞`, `w(·,τ) → f0(·)` uniformly, and `ψ(τ) → η0`.
* The convergence of the interface is completed by a semiconvexity argument in §7.
* In original variables (our reformulation of (4.4)–(4.5)):
  * `ζ(t) = η0 √(T−t)(1+o(1))`;
  * `u(0,t) = f0(0) √(T−t)(1+o(1))`.
* No rate of convergence for the `o(1)` is given.

**Intermediate bounds.** These hold for all `τ` and are not explicit numbers:
* (4.9) `|w_η| ≤ M`;
* (4.10) `w ≤ M1`;
* (6.1) `ψ ≤ M2`;
* (6.5) `ψ ≥ m2`;
* (6.6) `w(0,τ) ≥ m1`.

**Verdict.** Reproducible now with **(a)/(b)**, radial N = 2, 3, from *non-self-similar* data such as a parabola-like hump with slope 1 at the edge.
* Compare `ζ(t)/√(T−t) → η0` and `u(0,t)/√(T−t) → √(2M) f0(0)`.
* `T` is not explicit. Estimate it from the numerical extinction time, or by a fit `ζ² ≈ η0²(T−t)`.

### GHV-3. One-dimensional extinction, general (non-symmetric) data (§7, Theorem 7.1, eqs. (7.1)–(7.9), pp. 588–589)

**Setting.** N = 1. Any data with connected support. The extinction point is moved to 0. Interfaces (7.2): `w = 0`, `w_η = −1` at `η = ψ(τ)`; `w = 0`, `w_η = 1` at `η = φ(τ)`.

**Explicit profile (7.4).** `f0(η) = c0( e^{η²/4} − (η/2) ∫_0^η e^{ζ²/4} dζ )` on `(−η0, η0)`.
* (7.5): `η0` is the unique positive root of `η ∫_0^η e^{ζ²/4} dζ = 2 e^{η²/4}`.
* `c0 = 2(∫_0^{η0} e^{ζ²/4} dζ)^{−1}`.
* The paper notes `η0 > √2`.
* [ours] Numerically `η0 = 1.8482777460` and `c0 = 0.7868021977`.

**Theorem 7.1.** The convergence (6.20) `w → f0` holds, with `φ(τ) → −η0` and `ψ(τ) → η0`. The limit is **symmetric even for non-symmetric data**.

**Exact mass identity (7.6), p. 589.** Integrating `u_t = u_xx` gives `(d/dt)∫u = −2`, i.e. `∫u(t) dx = 2(T − t)`. In rescaled form, `∫ w dη = 2`.
* [ours] With `√(2M)`: `∫u(t) dx = 2√(2M)(T−t)`.
* So the extinction time is known exactly from the data: `T = ∫u0 dx / (2√(2M))`, as long as the support stays one interval.
* For symmetric data on the half-line with zero flux at the symmetry point, this becomes `T = ∫_0^∞ u0 dx / √(2M)`.

**Semiconvexity (7.9).** `w_ηη ≤ M2 e^{−τ/2}`.

**Verdict.**
* Symmetric data: reproducible now with **(b)** planar. `solve_planar` has zero flux at `L`, so `(s, L)` is half of a symmetric interval. It is also reproducible with **(a)** planar.
* The exact `T` gives a sharp test of both the limit solver and of `T_ε → T`.
* Non-symmetric data (the symmetrization claim in Theorem 7.1): **needs multiple fronts** (two interfaces).

### GHV-4. Self-similar profiles with distant or small interfaces (§5.1–5.2, eqs. (5.4)–(5.13), pp. 582–583)

**Second solution `f1`.**
* (5.4): `f1(η) ≈ a1 η` as `η → ∞`.
* (5.5): as `η → 0`, `f1 ≈ c1 log η` (N = 2) or `f1 ≈ −c1 η^{2−N}` (N ≥ 3).
* `f1` is increasing, with a unique zero `η1(N) < η0` where `f1'(η1) = 1`.
* (5.6): the general solution is `f = A f0 + B f1`.

**Family `f(η; l)`.**
* (5.7): `f(l; l) = 0`, `f'(l; l) = −1`.
* (5.8)–(5.9a): `f(η; l) = A(l) f0 + B(l) f1`, with `A(l) = f1(l)/d(l)` and `B(l) = −f0(l)/d(l)`.
* (5.9b): `d(l) = ½ a0 a1 l^{1−N} e^{l²/4}`.
* Limits:
  * (5.10): `η*(l) → 0` and `f'(η*; l) → ∞` as `l → η0+`;
  * (5.11)–(5.12): `f(η; l) → 0`, `η*(l) → η1−` and `f'(η*; l) → 0` as `l → ∞`.

**Verdict.** These are auxiliary objects for comparison arguments. **Not a numerical target**, except through GHV-5 and GHV-8.

### GHV-5. Annular self-similar solution: single-point extinction on an annulus (§5.3, Proposition 5.1, eq. (5.15), p. 584)

**Formula.** There exists `l* > η0` such that `f*(η) = f(η; l*)` satisfies (5.14): `f*'(η*(l*)) = 1`.
* (5.15): `u*(x,t) = (T−t)^{1/2} f*(η)`, `η = |x|/(T−t)^{1/2}`.
* Outer interface: `|x| = (T−t)^{1/2} l*`. Inner interface: `|x| = (T−t)^{1/2} η*`.
* Both reach the origin at the same time `T`.

**Behaviour.** This is a "new type of single-point extinction". It "must be unstable", since it needs both interfaces to arrive exactly together. It is the borderline between GHV-2 (the hole closes first) and GHV-6 (extinction on a sphere).

**Numerical values [ours].** Shooting backwards from `f(l) = 0`, `f'(l) = −1` and solving `f'(η*) = 1`:
* N = 2: `l* = 3.70060`, `η* = 0.281527`;
* N = 3: `l* = 4.07471`, `η* = 0.724703`.

**Verdict.** **Needs multiple fronts/annulus.**
* A good target once two radial fronts are available: start from `u*(·,0)` and observe the instability. Perturbing the data should give either a closing hole or extinction on a sphere.

### GHV-6. Extinction on a sphere (§8, Lemma 8.1, Theorem 8.2, eqs. (8.1)–(8.28), pp. 595–600)

**Data.** Annular support, e.g. (8.1) `supp u0 = {1 < r < 2}`. Condition (8.2):
* N = 2: `u0(r) ≤ min{log r, 2 log(2/r)}`;
* N ≥ 3: `u0(r) ≤ (1/(N−2)) min{1 − r^{−(N−2)}, 2^{N−1} r^{−(N−2)} − 2}`.

Under (8.2), Lemma 8.1 gives (8.3): `supp u(·,t) ⊂ [1,2]` for all t, so `0 ∉ E`.

**Rescaled problem.**
* (8.4): `w = (T−t)^{−1/2} u`, `η = (r − r*)/(T−t)^{1/2}`.
* (8.5): `w_τ = A(w) + e^{−τ/2} B(w)`.
  * (8.6): `A(w) = w_ηη − ½ w_η η + ½ w`;
  * (8.7): `B(w) = ((N−1)/(r* + e^{−τ/2} η)) w_η`.

**Result.** Theorem 8.2: there is a unique `r* > 0` such that `w → f0` of the **one-dimensional** problem (7.4), and the interfaces converge.
* So `ζ(t) − r* ≈ η0^{(1D)} √(T−t)` and `ξ(t) − r* ≈ −η0^{(1D)} √(T−t)`, with `η0^{(1D)} = 1.84828` [ours].
* Mass (8.19): `E(τ) = ∫ w dη = 2 + O(e^{−τ})`.
* (8.27)–(8.28): `∫f = 2` and `f'' ≤ 0`.
* Figure, p. 595: a scheme of the three extinction types for annular data (single point, annular type, sphere).

**Normalization `|∇u| = 1`.** [ours] With `√(2M)`, the bound (8.2) is multiplied by `√(2M)`.

**Verdict.** **Needs multiple fronts/annulus.**
* Natural first target after the two-front extension. Measure `r*`, the symmetric √ collapse, and `∫u dr` near `T`.

### GHV-7. Focusing at the origin (§9, Theorems 9.1–9.2, eqs. (9.1)–(9.18), pp. 600–603; Remark 1, (9.35), p. 605)

**Data.** (9.1): N ≥ 2, `u0 > 0` and `u0' > 0` on `(1, ∞)`, `u0(1) = 0`, with `|u0'| ≤ M` and `|u0''| ≤ M`.
* The fresh region is the exterior of the hole `B_1`.
* (9.3): `u0(x) ≥ u1(x,0,T)` for some T guarantees that the hole closes in finite time.

**Self-similar model (9.2).** `u1(r,t;T) = (T−t)^{1/2}[f1(r/(T−t)^{1/2})]_+`, with `f1` from GHV-4.

**Results.**
* Theorem 9.1:
  * (9.4) `w(η,t) = (T−t)^{−1/2} u(η(T−t)^{1/2}, t) → [f1(η)]_+` uniformly on compact sets;
  * (9.5) `ζ(t) = η1 (T−t)^{1/2}(1 + o(1))`.
* Theorem 9.2, the profile at the focusing time: (9.7) `u(r,T) = a1 r(1 + o(1))` as `r → 0`.

**Explicit cases.**
* N = 3, (9.8): `f1(η) = ½(η − 2/η)`, `η1 = √2`, `a1 = ½`.
  * [ours] In original variables, `u1(r,t) = ( r/2 − (T−t)/r )_+`. This is an elementary closed-form exact solution, with free boundary `r = √(2(T−t))`.
* N = 5, (9.9): `f1(η) = a1(η − 4/η − 4/η³)`, with `η1 = [2(1+√2)]^{1/2}` and `a1 = (3+2√2)/(4(2+√2))`.

**After focusing.** Remark 1, (9.35): "For t > T the hole and the free boundary disappear", and `u(0,t) ≈ c2 (t−T)^{1/2}`. `c2` is not given.

**Expanding waves.** Remark 2: the interface can instead go to infinity, behaving like a 1D travelling wave. "Radially symmetric stationary solutions … are known to be unstable" (no proof).

**Comparison subsolution (9.10).** `ū(r,t) = (T0−t)^{1/2} c_l f(η; l)`, `η = r/(T0−t)^{1/2}`, with `c_l f'(η*(l); l) ≥ 1`. Lemma 9.3: if `u0 ≥ ū(·,0)` for `r ≥ 1`, focusing occurs with `T(u0) ≤ T0`.

**Numerical values [ours].** `f1 ∝ U(−½, N/2, η²/4)`, the Tricomi function. This reproduces the paper's `η1`, `a1` for N = 3 and N = 5 to 10 digits.
* N = 2: `η1 = 0.8838633896`, `a1 = 0.5798256453`. Not in the paper.
* N = 3: `η1 = 1.4142135624`, `a1 = 0.5`.
* N = 5: `η1 = 2.1973682269`, `a1 = 0.4267766953`.

**1D contrast.** p. 600, and Vázquez 1996 (89): in 1D the interface arrives with finite speed, `ζ(t) = O(t1 − t)`.

**Verdict.** The fresh region is outside a hole, with linearly growing data.
* **(b)** cannot do this now. It needs an exterior radial limit solver with Dirichlet data `u1(L,t)` at the outer boundary. This is a moderate extension of `limit.py`.
* **(a)** may work if the burnt-inside orientation is supported (untested).
* Single front, so no multiple-front machinery is needed. High priority: there is an explicit exact solution for N = 3 and a universal rate.

### GHV-8. Classification of annular extinction modes (Theorem 3.2(ii), p. 579)

**Statement.** For annular data, both interfaces converge as `t → T_m`.
* If the smaller limit is positive, the two limits coincide: extinction on a sphere.
* If the smaller limit is zero, the larger may be zero (annular extinction at a point) or positive. In the latter case the solution continues as a ball solution.
* Remark (p. 579): when only the inner sphere collapses, the outer free boundary is `C^∞` at that time.

**Verdict.** **Needs multiple fronts/annulus.** Gives the qualitative phase diagram for annulus experiments: hole closes first, sphere extinction, or the borderline GHV-5.

---

## 3. Vázquez 1996, "The free boundary problem for the heat equation with fixed gradient condition" (lecture notes, Zakopane 1995) — `Vazquez1996.pdf`

The text layer is garbled. Everything below was read from rendered page images; page numbers are the printed ones, which equal the pdf pages.

**Normalization.**
* (5), p. 2: `∂u/∂n = −1`. "The constant in the second member of (5) is put to 1 as a normalization."
* §3 derives it from the travelling wave analysis as `|u'(0)| ≈ √(2M)/√D` (34), with `M = ∫β(u)du` (29).
* §4, p. 11: "u is rescaled so that the slope √(2M) becomes 1 and D is also put to 1".
* (41), p. 12: `β_ε(s) = ε^{−1}β(s/ε)`, `∫β = 1/2`. The heat production `M` is "normalized to 1/2".

### V-1. Travelling wave derivation of the jump condition (§3, eqs. (29)–(35), pp. 9–10)

**Formulas.**
* Fresh zone, (30): `Du'' = cu'`, `u(0) = h`, `u(−∞) = A`. Its solution is (31): `u(ξ) = A(1 − e^{λξ}) + h e^{λξ}`, with `λ = c/D`.
* Reaction zone, (32)–(33): `(D/2)(u'(ξ))² + c∫_ξ^∞ (u')² = ∫_0^h β(u)du` and `Du'(0)² = 2M − 2c∫(u')²`.
* As `β` concentrates: (34) `|u'(0)| ≈ √(2M/D)`, and (35) `c = √(2M)/(AD)`, with `A = H0 − T0`.
* (37): the KPP minimal speed for comparison, `c = √(2β'(0))`.

**Verdict.** The ε-problem TW speed against the limit is already done with **(c)** and **(a)**: `exact.py` has `c0 = √(2M)/A` and `tw_speed_eps`. Formula (35) is the same as ours with `D = 1`.

### V-2. Explicit classical solutions (§6, eqs. (44)–(47), pp. 14–16)

* **Travelling waves**, N = 1:
  * (44) receding: `u(x,t) = (1/c)(1 − e^{c(x+ct)})_+`, moving toward the negative x-axis;
  * (45) progressing: `u(x,t) = (1/c)(e^{c(ct−x)} − 1)_+`;
  * `c = 0` gives `u = (−x)_+`.
  * These are CV (2.2)/(2.3) reflected.
* **Self-similar solution** (46): `u = (t1 − t)^{1/2} f(|x|/(t1 − t)^{1/2})`, with free boundary (47) `|x| = R√(t1 − t)`. Remarks:
  * extended by zero it is globally Lipschitz for `0 ≤ t < t1`, with maximum gradient 1 at the free boundary;
  * it is only `C^{1/2}` at extinction;
  * near extinction the curvature radii become comparable to the reaction-zone thickness, so the model is valid only as intermediate asymptotics.
* **Stationary holes**, p. 16: identical to CV (2.8)–(2.10), including `A = −R log R`, `A* = 1/e`, and the Zeldovich flame for N = 3.

**Verdict.** Same as CV-1, CV-2 and CV-3.

### V-3. Kinetic equation of the front (§7, eqs. (48)–(50), p. 17)

**Formulas.**
* (48): `v = n u_t = n Δu`.
* Radial, (49): `v = n(∂_r²u|_Γ − (N−1)/r(t))`.
* General, (50): `v = n(∂_n²u|_Γ − K(x,t))`, with K "the Gauss curvature" (sic), linking the problem to motion by curvature.

**Verdict.** A diagnostic for **(b)**: `limit.py` already uses `s' = Δu/λ`.

### V-4. Relation to the (undercooled) Stefan problem (§11, eqs. (65)–(70), pp. 21–23)

**Formulas.** In 1D, on `{0 < x < ζ(t)}` with `u_x(0,t) = g(t)`:
* (66): `w = −1 − u_x` solves the heat equation;
* (67): `w(0,t) = g(t)`;
* on the front, `w = 0` and (68) `v = u_xx = −w_x`, i.e. the Stefan condition.
* If `u_x ≥ −1` on the parabolic boundary, the result is an **undercooled** Stefan problem with `w < 0`. "This is precisely what happens for the self-similar solutions (46)."
* Two-interface case, (70): `w = −2`, `w_x = v_l` on the left interface.

**Verdict.** **Not a numerical target.** It is a possible independent check of the 1D limit solver through a Stefan code.

### V-5. Nonuniqueness: two humps (§13, eq. (79), pp. 27–28)

**Construction.**
* A radial hump expands until time `t1` to its maximum radius `r1`, then collapses.
* Two such humps centred `x1` apart:
  * if `x1 > 2r1` they never interact;
  * at `x1 = 2r1` one solution is the superposition, touching at `t = t1` and then separating.
* A different continuation merges the supports at `t1` and never separates. This is seen first in 1D with `x1 < 2r1`: the inner pair of fronts is cancelled and the solution continues with support the union of the two intervals.
  * "we will lose twice the production rate (normalized to 1)".
  * Taking `x1 → 2r1` gives a second solution for `x1 = 2r1`.

**Mass rule.** (79), as printed: `∫u2(x,t)dx − ∫u1(x,t)dx = 2`, "so that even their extinction times are different".
* [ours] By CV (3.8) in 1D, the difference of masses grows at rate 2 per unit time after `t1`. So (79) presumably means `d/dt` of the difference, or `2(t − t1)`; the printed formula has no `t`. Flagged, not corrected.
* "Clearly, one of them is the maximal and another one the minimal solution."

**Pulse splitting**, p. 28: with `x1 > 2r1` and a thin sliver of data joining the supports along `{|x'| < ε}`, the connected support splits in two after a short time that depends continuously on the sliver size. "Detailed proofs will appear elsewhere."

**Verdict.**
* 1D minimal solution: **needs multiple fronts** (two intervals = four interfaces, or two on the symmetric half-line).
* 1D merged (maximal) solution: after `t1` it is a single symmetric interval, so one front on the half-line. It is reproducible with **(a)/(b)** from data at time `t1`.
* The interesting experiment is which solution **the ε-problem selects** for touching data. That needs (a) to handle a burnt *point or interval* inside the fresh region.
* Pulse splitting: **needs 2D**.

### V-6. Asymptotics: extinction and focusing (§14, eqs. (80)–(89), pp. 28–31)

A summary of GHV, with the same formulas:
* (82): focusing self-similar solution;
* (83): N = 3 values `η1 = √2`, `a1 = ½`;
* (84): data conditions;
* (86)–(88): Theorems 5–6, equal to GHV 9.1–9.2.

New statement here: (89), in 1D focusing "generally" has `ζ(t) = O(t1 − t)`, i.e. finite arrival speed, "in sharp contrast with (87)".

The text also lists the three extinction modes: (i) single point, (ii) sphere, (iii) annular self-similar with support size `O(τ − t)` (sic, presumably `O((τ−t)^{1/2})`; the claim is copied as printed).

**Verdict.** See GHV-2, GHV-3, GHV-5, GHV-6 and GHV-7. The 1D focusing rate (89) is reproducible with **(a)** in planar 1D: the burnt interval sits at the symmetry wall and fresh fronts close in. The (b) planar solver has the burnt region on the left and zero flux at L, so it would also need Dirichlet data at L to drive the front to x = 0.

### V-7. Other models (§§8–10, 12, pp. 17–24)

* Florin filtration, (51)–(54).
* Deflagration–detonation transition (Burgers two-phase), (56)–(59).
* Elliptic–parabolic model, (60)–(64), and its use for radial solutions, (74)–(78) = GHV §1.
  * The radial harmonic extension is `u(r,t) = (ζ/(N−2))[1 − (ζ/r)^{N−2}]` (p. 25).
* Meirmanov's and Andreucci–Gianni's existence results.

**Verdict.** **Not a numerical target.** The elliptic–parabolic reformulation (78) is a possible alternative limit solver.

---

## 4. Petrosyan & Yip, J. Geom. Anal. 18 (2008), 1098–1126 — `PetrosyanYip2008.pdf`

**Normalization.**
* (P): `Δu − ∂t u = 0` in `{u>0}`, `|∇u| = 1`.
* (P_ε): `Δu − ∂t u = β_ε(u)`, `supp β_ε = [0,ε]`, `∫β_ε = 1/2`.
* (2.1)–(2.3): β is Lipschitz, positive on (0,1), with `∫β = 1/2`, and `β_ε(s) = ε^{−1}β(s/ε)`.

### PY-1. Two touching humps: nonuniqueness of limit solutions (Example 1.1, eqs. (1.1)–(1.2), pp. 1099–1100; Theorem 1.2, p. 1101; Fig. 1)

**Data.**
* `u0` is radial and smooth in `B̄_{r0}`, with `{u0 > 0} = B_{r0}` and `|∇u0| = 1` on `∂B_{r0}`.
* The classical radial solution has radius `r(t)`, which tends to 0 as `t → T_ext`.
* (1.1): `r* = r(t*) = max r(t)`. Both `t* = 0` and `t* > 0` are possible, depending on the shape of `u0`.
* (1.2): `w0(x) = u0(x − r* e_n) + u0(x + r* e_n)`.

**Two solutions.**
* Minimal: approximate `u0` from inside. The superposition `w = u(x − r*e_n, t) + u(x + r*e_n, t)` touches at `t*`, then separates and collapses at `T_ext`.
* Maximal: approximate from outside. The supports overlap at `t*` and may merge.

**Results.**
* Theorem 1.2: for `n ≥ 2` the maximal and minimal limit solutions differ.
* In `n = 1` this "can be validated by a simple use of the maximum principle". Fig. 1 shows the minimal and maximal 1D solutions.
* Theorem 2.1: maximal and minimal limit solutions always exist.

**Verdict.**
* `n ≥ 2`: **needs 2D**.
* `n = 1`: **needs multiple fronts**.
* [ours] Simplest version: choose `u0` with `t* = 0`, i.e. the hump shrinks from the start. Then `w0` is a single interval `[−2r0, 2r0]` with an interior zero at `x = 0`, where `w0 ≈ |x|`.
  * The maximal solution is then a single symmetric interval. It is reproducible with (b), using zero flux at the centre after a reflection.
  * The minimal solution needs two fronts on the half-line.
  * Running **(a)** with this `w0` and various `u0ε` near `x = 0` tests the selection mechanism directly. Compare masses: the merged and separated solutions lose mass at rates 2 and 4 (times `√(2M)`), by CV (3.8).

### PY-2. Catenoidal Alt–Caffarelli minimizer and radial barrier (§3, Proposition 3.2, eqs. (3.2)–(3.4), p. 1108)

**Explicit barrier.** `ψ(x) = (1/(n−2))(|x|^{2−n} − 1)` for `n > 2`, and `ψ = log(1/|x|)` for `n = 2`, on `B1 \ B_{1/2}`, with `|∇ψ| = 1` on `∂B1`.
* Extended by 0 outside `B1`, it minimizes (AC) in `B_ρ \ B_{1/2}`.
* Proposition 3.2: in a flat cylinder `Q_{a,1}` with boundary data `1 − κ` in the middle, the minimizer has a connected "catenoidal" positivity set containing `Q_{1,1}`, when `a > a0` and `κ < κ0`. No explicit values of `a0`, `κ0` are given.

**Verdict.** Elliptic. **Needs 2D** (axisymmetric 3D).

### PY-3. Self-similar "neck" subsolution (§4, eqs. (4.8)–(4.17), Proposition 4.1, pp. 1112–1115)

**Formula.** `υ(x,t) = √(αt) v̂_c(x/√(αt))`, where `Δv̂_c = c` is a perturbed catenoidal minimizer. It is joined to two shrinking radial solutions `θ±(x,t) = θ(x ∓ ρ0 e_n, t)` (Fig. 4).

**Verdict.** A proof device. **Not a numerical target.**

### PY-4. Analogy with motion by mean curvature (§6, pp. 1119–1121)

* §6.1 restates the GHV annulus trichotomy (focusing, extinction on a circle, extinction on an annulus) in terms of `a*`, citing GHV p. 595. It is analogous to the shrinking torus in MMC.
* §6.2: two humps joined by a neck. If the neck is too narrow it pinches; if it is thick the support becomes convex; a critical thickness is anticipated. This is the analogue of Hamilton's "incredibly shrinking dumbbell".
* Both are stated as conjectures or directions.

**Verdict.**
* §6.1: **needs multiple fronts/annulus**. This is the same experiment as GHV-6.
* §6.2: **needs 2D**.

### PY-5. Explicit pinching supersolution (Appendix A.3, pp. 1123–1125, Fig. 5)

**Formula.**
* `u(x1,x2,t) = (h² − x2²)_+/h = (h − x2²/h)_+`;
* `h(x1,t) = min{g(x1), 2} − at`, with `g(x1) = 1 + (bx1)²`.

**Properties.**
* `supp u(·,t) = {|x2| ≤ h}`, and `|∂_{x2} u| = 2` on the boundary.
* Pinching at `t1 = 1/a`; the function vanishes at `t2 = 2/a`.
* It is a supersolution if `a + b² ≤ 1/2` and `0 ≤ a ≤ 1/2`, after multiplying by a small constant so that `|∇u| ≤ 1`.
* ε-version: `u_ε = φ_ε(u)`, where `φ'' = β(φ)`, `φ(−∞) = 0`, `φ'(+∞) = 1`, `φ(0) = 1`, and `φ_ε(s) = εφ(s/ε)`.
* The unburnt zones `{u_ε > ε} = {u > 0}` pinch in finite time, independently of ε.

**Verdict.** **Needs 2D.** It is a supersolution, not a solution, so it gives no quantitative target.

---

## 5. Weiss & Zhang, arXiv:0904.1261 (2009) — `WeissZhang2009.pdf`

**Normalization.** Elliptic, two-phase.
* `Δu_ε = β_ε(u_ε)`, with `β_ε(s) = ε^{−1}β(s/ε)`, β Lipschitz, positive on (0,1), `∫_0^1 β = ½`.
* Limit (1.2): `(u_ν^+)² − (u_ν^−)² = 1`.

### WZ-1. Cross-shaped domain-variation solution (Introduction, p. 2, Fig. 1)

`u(x1,x2) = |x1² − x2²|` is a domain-variation solution with a cross-shaped free boundary. It is cited from Andersson–Weiss. It is not known whether it arises as an ε-limit.

### WZ-2. Triple junction (p. 2, Fig. 2)

Andersson's example, homogeneous of degree 3/2. "Domain variation solutions with a degree of homogeneity ∈ (1, 3/2) are not possible."

### WZ-3. Degenerate singularity as an ε-limit (Theorem A, p. 3; §3, eq. (3.1); Theorem 6.4; Remark 6.5)

**Construction.**
* (3.1): solve `Δu = β_ε(u(x) − u(0) + ε)` in the sector `K = {0 < r < 1, 0 < θ < π/m}`, `m ≥ 3`, with `u = g` on the arc and `∂_ν u = 0` on the straight sides.
* Reflect evenly, giving `2m` symmetry sectors.
* `v_ε = u_ε − u_ε(0) + ε` satisfies `v_ε(0) = ε`.

**Result.** The limit has a degenerate free boundary point at the origin: `u(x0 + rx)/r → 0` (Theorem A).
* Remark 6.5: arbitrarily many symmetry lines are possible, suggesting growth `≤ |x|^m`.
* Fig. 3 shows an "imaginable ε-solution" with 4 symmetry lines.
* The singularity is expected to be unstable. §3 "suggests a numerical approach" for constructing it: the fixed-point problem (3.1) on one sector.

**Verdict.**
* **Needs 2D** (elliptic). It is a possible future elliptic 2D benchmark: compute (3.1) in a sector and check `sup_{B_r} u = o(r)`.
* Note that the problem is two-phase (`u` can be negative).

---

## 6. Lederman, Vázquez & Wolanski, Trans. AMS 353 (2001), 655–692 — `LedermanVazquezWolanski2001.pdf`

**Normalization.**
* `Δu + Σ a_i u_{x_i} − u_t = 0` in `{u > 0}`, `u = 0`, `−∂u/∂ν = √(2M)` (p. 656).
* `β_ε(s) = ε^{−1}β(s/ε)` (1.1), with `∫β = M`.
* **This is exactly our normalization.**

**Examples.** There are no explicit solutions. Two remarks act as counterexamples:
* **Remark 3.1** (p. 664) and **Remark 6.1** (p. 685): the bounded solution of the heat equation with the same Dirichlet/Neumann data is a viscosity solution, a "counterexample" to uniqueness if the free boundaries are not required to start correctly.
  * If `u_ε ≥ ε` on `∂_D Q`, the limit is the heat equation solution (same as CV-6).

**Main results.** In a semi-cylinder `(0,∞) × Σ`, or the full cylinder `R × Σ`, with Neumann data on the lateral boundary:
* Theorems 3.1 and 3.2: a bounded classical solution with `u_{x1} < 0` on `{u>0} ∩ ∂_D Q` is the unique viscosity solution. This means monotone in the cylinder axis, with a bounded free boundary.
* Theorem 6.1: any family of ε-solutions with data converging on `∂_D Q`, and positivity sets converging, converges uniformly to it. The limit is independent of β.

[ours] The 1D version (`Σ` trivial) is the half-line `x1 > 0` with Dirichlet data `u(0,t) > 0` and a single receding front. This is the geometry of a flame entering a tube.

**Verdict.** **Not an explicit example.** The theorem guarantees that (a) → (b) convergence in 1D planar with a monotone profile is well-posed and independent of `f`. A test with Dirichlet data at `x = 0` would need a small change in (b), which currently has zero flux at `L`.

---

## 7. Petrosyan, Comm. PDE 27 (2002), 763–789 — `Petrosyan2002.pdf`

**Normalization.** `|∇u| = 1`, with `supp β = [0,1]` and `∫β = 1/2` (1.1).

**Examples.** No explicit solutions. The relevant results are uniqueness and regularity statements for classes of data:
* Theorem 2.7: if `u0` is starshaped in the sense (S), `u0(λx) ≥ u0(x)` for `λ ∈ (0,1)` (p. 766), and `u0ε → u0` uniformly with supports converging, then the limit solution is unique and equals the minimal supersolution.
* Theorems 8.1 and 9.3: if `u0` satisfies (C1), i.e. superharmonic and log-concave on a bounded convex `Ω0` (p. 770), plus the technical conditions (C2)–(C3), then the limit solution is classical up to the extinction time. The sections `Ω(t)` are convex, shrinking and `C^{2,α}`.

**Verdict.** **Not a numerical target.** It certifies that our radial hump experiments (concave, log-concave data in a ball) have a unique, classical limit up to extinction, so the GHV-2 comparison is meaningful. It also says non-starshaped data (two humps) is where nonuniqueness can occur (PY).

---

## 8. Weiss, Calc. Var. 17 (2003), 311–340 — `Weiss2003.pdf`

**Normalization.**
* (1): `∂t u_ε − Δu_ε = −β_ε(u_ε)`, with `β_ε(z) = ε^{−1}β(z/ε)`, `β ∈ C0¹([0,1])`, `∫β = ½`.
* (2): `|∇u| = 1`.

### W03-1. Two colliding travelling waves (eq. (3), p. 312, Fig. 1)

**Formula.** `u(t,x) = χ_{x+t>1}(exp(x + t − 1) − 1) + χ_{−x+t>1}(exp(−x + t − 1) − 1)` for `t ∈ [0,1)`.
* Two advancing (expanding) waves with unit slope at the fronts `x = ±(1 − t)`.
* The burnt interval `(t−1, 1−t)` closes at `t = 1`.

**Claimed behaviour.** "Even for flawless initial data, classical solutions of (2) develop singularities after a finite time span". The parabolic blow-up limit at `(t,x) = (1,0)` is `|x|` (p. 313). This is the 1D analogue of focusing, with a finite arrival speed (cf. Vázquez (89)).

**Verdict.** Symmetric, so one front on the half-line `x ≥ 0`, with the burnt region at the symmetry wall. The data grow exponentially, so a Dirichlet condition `u(L,t) = e^{L+t−1} − 1` is needed at `x = L`.
* (a): possibly reproducible now, if the burnt-side orientation and Dirichlet data at L are supported.
* (b): needs Dirichlet data at L.
* Interesting quantities:
  * the ε-dynamics near the collision;
  * the slope at `x = 0` just after `t = 1`, where `|x|` is the blow-up limit;
  * the post-collision mass.

### W03-2. 1D stationary ε-solutions with arbitrary slope (eq. (5), p. 313; Remark 4.1, p. 324)

**Formula.**
* `u_ε'' = β_ε(u_ε)` in R, `u_ε'(0) = 0`.
* `u_ε(0) = B_ε^{−1}(δ)` for `δ ∈ (0, ½]` (case `δ = 0` as printed: `B_ε^{−1}(ε)`).
* Here `B(z) = ∫_0^z β`, `B_ε(z) = ∫_0^z β_ε`, and `B_ε^{−1}` is the inverse on `[0,ε]`.

**Claim.** `u_ε → √(1 − 2δ)|x|` as `ε → 0`.
* So limits can have any slope in `[0,1)` at a two-sided zero, and need not satisfy (4).
* Remark 4.1: with `δ = ½` this "suggests that χ can be different from `χ_{u>0}`". The limit is then identically 0 while the reaction is concentrated.

**Verdict.** Stationary 1D. Reproducible now with an ODE shooting check, or with (a) as a stationary state.
* It is a good **unit test of the ε-limit of the two-sided slope**, directly relevant to merged fronts (V-5, W03-1, PY-1).
* [ours] In our normalization the analogue is `√(2M − 2B(u_ε(0)/ε))|x|`, by the same first integral. This is not stated in the paper.

### W03-3. Stationary distributional solution in 2D (p. 313, Fig. 2)

**Formula.** `u(x) = max(−log|x − e1|, 0) + max(−log|x + e1|, 0)` in `R² \ (B_{1/2}(e1) ∪ B_{1/2}(−e1))`. The two annuli `{1/2 < |x ∓ e1| < 1}` touch at `x = 0`, where the blow-up is `|x_n|`.

[ours] Each term is the CV hole solution (2.8) with `R = 1`, reflected inside.

**Verdict.** **Needs 2D.**

### W03-4. Citations of Alt–Caffarelli examples (p. 312)

The Alt–Caffarelli singular cone for `n = 3` (their example 2.7) and nonuniqueness for the 2D Dirichlet problem (their example 2.6). See AC below.

---

## 9. Alt & Caffarelli, J. Reine Angew. Math. 325 (1981), 105–144 — `AltCaffarelli1981.pdf` (elliptic Bernoulli problem)

**Normalization.**
* Minimize `J(v) = ∫_Ω (|∇v|² + χ({v>0}) Q²)` (1.2), with `v = u0` on `S ⊂ ∂Ω`.
* The free boundary condition is `∂_{−ν}u = Q`, i.e. `|∇u| = Q`.
* `Q = 1` in all examples.
* [ours] Our elliptic analogue `Δu = β_ε(u)` with `∫β = M` corresponds to `Q = √(2M)`.

**Notation (1.5), p. 107.** Fundamental-solution-type function:
* `φ_R^m(r) = (R/(m−2))((r/R)^{2−m} − 1)` for `m ≥ 3`;
* `φ_R^m(r) = R log(R/r)` for `m = 2`;
* `φ_R^m(r) = R − r` for `m = 1`.

### AC-1. Linear and radially symmetric solutions; nonuniqueness of the Dirichlet problem (§2.6, pp. 109–110)

**Exterior (hole) case.** `Ω = {x = (x¹,x²) ∈ R^m × R^{n−m} : |x¹| ≥ ρ, |x²| ≤ 1}`, `S = {|x¹| = ρ}`, `u0 = const > 0`, `Q = 1`.
* For `R > ρ`, `u = max(0, φ_R^m(|x¹|))` is harmonic in `{u>0}`, with `∂_{−ν}u = 1` on the free boundary.
* It matches the data iff `φ_R^m(ρ) = u0`. For given `u0, ρ` there is exactly one such `R > ρ`, and `u` is the unique absolute minimum.

**Interior case.** `|x¹| ≤ ρ`, `u = max(0, −φ_R^m(|x¹|))` with `0 < R < ρ`, condition `−φ_R^m(ρ) = u0`.
* A solution exists iff `0 < u0/ρ ≤ c_m`, where `c_m = (m−1)^{−(m−1)/(m−2)}` for `m ≥ 3` and `c_m = 1/e` for `m = 2`.
* If `u0/ρ < c_m` and `m ≥ 2`, there are exactly two solutions `0 < R1 < R2 < ρ`, called `u¹` and `u²`.
* `J(u¹) > J(u0)`, and `J(u0) ≥ J(u²)` iff `R2/ρ ≥ (2/m)^{1/(m−2)}` (`m ≥ 3`) or `R2/ρ ≥ 1/√e` (`m = 2`).
* `u¹` is a saddle point ("hyperbolic solution"). `u0` (the trivial solution `≡ u0`, no free boundary) and `u²` are local minima.

[ours] This is the elliptic counterpart of CV (2.9): `A* = 1/e`, big and small hole.

**Verdict.** Radially symmetric elliptic problem.
* Reproducible with a 1D radial elliptic or stationary solver, i.e. the stationary limit of (a)/(b) with Dirichlet data. Not covered by the current time-dependent code without a Dirichlet boundary.
* For the parabolic code, the two radial solutions are stationary states of the evolution. Their (in)stability is a target, with fresh region outside or inside.

### AC-2. Singular cone in n = 3 (§2.7, p. 110)

**Formula.**
* `u(x) = r max(f(θ)/f'(θ0), 0)` in polar coordinates `x = r(cos φ sin θ, sin φ sin θ, cos θ)`;
* `f(θ) = 2 + cos θ log((1 − cos θ)/(1 + cos θ))`, which solves `(sin θ f')' + 2 sin θ f = 0` with `f'(π/2) = 0`;
* `θ0 ≈ 33.534°` is the unique zero of `f` between 0 and `π/2`.

**Properties.**
* Harmonic in `{u>0}`, with `∂_{−ν}u = 1` on `∂{u>0} \ {0}`.
* Singular free boundary (a double cone) at the origin.
* `|∇u| ≤ 1 − cε` for `θ0 + ε ≤ θ ≤ π − θ0 − ε`.
* It is a weak solution but **not a minimum** (stated).

**Verdict.** **Needs 2D.** It is axisymmetric, so an axisymmetric (r,z) elliptic code could do it. It is a classic benchmark for the elliptic study: check `θ0` and the instability of the cone.

### AC-3. Two touching balls (§5.7, p. 125)

**Formula.**
* `Ω = {x ∈ Rⁿ : |x − e_n| > ½ and |x + e_n| > ½}`;
* `u = φ_1^n(|x − e_n|)` in `|x − e_n| < 1`, `u = φ_1^n(|x + e_n|)` in `|x + e_n| < 1`, and `u = 0` otherwise.

**Properties.** A weak solution for `Q = 1`, with a singular point: `∂{u>0} = {0} ∪ ∂_red{u>0}`. `{u = 0}` has density zero at the origin, so it is **not a local minimum**.

[ours] For n = 2 this is exactly Weiss's example W03-3.

**Verdict.** **Needs 2D.**

### AC-4. Two parallel planes, ε-gap (§5.8, p. 125)

**Formula.** `u = q₊(x_n − ε)` for `x_n ≥ ε`, `0` for `|x_n| ≤ ε`, and `q₋(ε − x_n)` for `x_n ≤ −ε`.
* For `ε > 0`, it is a weak solution if `Q = q±` on `x_n = ±ε`.
* For `ε = 0` it is not, since `q_u = q₋ + q₊`: the reduced boundary is empty.

**Verdict.** 1D, elliptic. It is the stationary picture of two fronts meeting, and relates to W03-2. Reproducible in principle with (a) planar as a stationary configuration, but low priority.

### AC-5. Discontinuous Q (§5.9, p. 125)

**Formula.** `Ω = (0,∞)`, `S = {0}`, `Q = Q1` on (0,1) and `Q2` on (1,∞), with `Q1 < u0 < Q2`. The absolute minimum is `u = u0 max(1 − x, 0)`, so `∂_{−ν}u(1) = u0`. The free boundary condition is set by the data, not by Q.

**Verdict.** Not relevant to our constant-M problem. **Not a numerical target.**

---

## Papers in `pdfs/` not inventoried

`HuangRussell2011.pdf` (book, 446 pp.) and `ShenDaiHuang2024.pdf` were not in the list and were not read.

---

## Recommended comparisons

Ranked by how much each comparison tells us, and whether it fits the code now.

### 1. 1D extinction with an exact extinction time (GHV-3, (7.4)–(7.6), Theorem 7.1)

**Fits now:** (b) `solve_planar` (zero flux at L means a symmetric interval) and (a) planar.

**Data.** Any symmetric hump on `(−a, a)`, e.g. concave with edge slope ≤ `√(2M)`.

**Compare:**
1. Extinction time. The exact value from the data is `T = ∫_0^a u0 dx / √(2M)` (half-line), from (7.6) rescaled [ours]. Check `T_limit` against it, then `T_ε → T`. For the ε-problem, use `∫u0ε` in the same formula as a consistency check of the mass law.
2. Rescaled front: `(L − s(t))/√(T−t) → η0 = 1.848278`.
3. Peak: `u(center,t)/√(T−t) → √(2M)·0.786802`.
4. Mass linearity: `∫u dx` should be exactly linear in t for the limit solution.

This is the cleanest available check: an exact number for general, non-self-similar data.

### 2. Radial extinction from non-self-similar data (GHV-2, Theorem 6.6)

**Fits now:** (a)/(b) radial, N = 2, 3.

**Compare:**
1. `s(t)/√(T−t) → η0` (2.513927 for N = 2, 3.003951 for N = 3).
2. `u(0,t)/√(T−t) → √(2M)·f0(0)` (0.981715 and 1.105330).
3. The rescaled mass identity `((N+1)/2)E(τ) − ψ^{N−1} → 0`.
4. `T_ε − T_0` as ε → 0.

**Point of the test.** It shows the self-similar solution (c) is an *attractor*, not just an exact solution. Use T from the limit run; it is not explicit.

### 3. Focusing (GHV-7, Theorems 9.1–9.2), with the explicit N = 3 solution

**Exact solution.** `u1 = √(2M)( r/2 − (T−t)/r )_+`, from (9.8) after scaling [ours].

**Fit.** Single front, but the fresh region is *outside* the hole.
* (b) needs an exterior radial version with Dirichlet data `u1(L,t)` at `r = L`: a moderate change to `limit.py`.
* (a) may work if the burnt-inside orientation is supported and a Dirichlet condition is added at L.

**Compare:**
1. `ζ(t) = η1√(T−t)`: `η1 = √2` (N = 3), 0.883863 (N = 2).
2. Slope at the focusing time, `u(r,T) ≈ √(2M) a1 r`: `a1 = ½` (N = 3), 0.579826 (N = 2).
3. General monotone data (9.1): convergence to `f1`.
4. Post-focusing growth `u(0,t) ∝ (t−T)^{1/2}`, (9.35).
5. 1D contrast: linear arrival `ζ = O(t1 − t)`, from Vázquez (89) and Weiss (3).

### 4. Colliding planar fronts and the two-sided slope (W03-1, W03-2; also the merged branch of V-5 / PY-1)

**Fit.** 1D symmetric, one front on the half-line, with the burnt region against the symmetry wall. Fits (a) if the orientation is supported. (b) needs Dirichlet data at L for (3), or an equivalent bounded setup.

**Compare:**
1. Arrival time `t = 1` in (3).
2. The ε-solution's behaviour when the fronts meet. After the collision it is a smooth heat solution, but the ε-reaction layer has a stationary analogue (5) with limit slope `√(1−2δ)`.
3. Mass loss: rate 2, then 0, times `√(2M)`.

**Payoff.** A direct numerical look at how the ε-problem resolves front merging. This is the mechanism behind nonuniqueness in PY-1.

### 5. Annulus experiments (GHV-5, GHV-6, GHV-8; PY §6.1) — needs two fronts

**Fit.** Needs the two-front (annulus) extension. This is the first target once it exists.

**Compare:**
1. The three modes and the threshold between them: hole closes first, sphere extinction, or the borderline.
2. For sphere extinction: `r*`, and the 1D profile `f0` with `η0^{(1D)} = 1.848278` and `∫w dη → 2`.
3. The unstable annular self-similar solution: `l* = 3.70060`, `η* = 0.281527` (N = 2); `l* = 4.07471`, `η* = 0.724703` (N = 3). These values are [ours]. Start from `u*(·,0)` and watch it drift to one of the two stable modes.

### 6. Stationary radial holes and their instability (CV-3, AC-1)

**Fit.** Exterior fresh region with Dirichlet data. Same requirements as item 3.

**Compare:** perturbations of `u = √(2M) R log(r/R)` (N = 2) or `√(2M)(R − R²/r)` (N = 3) should either collapse (focusing) or expand (TW-like, GHV Remark 2). Also check the big/small-hole branches against `A = −R log R`.

**Priority.** Lower than items 1–3, because the instability is only asserted in the papers, not quantified.
