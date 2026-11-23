# Referee report (SINUM)

**Manuscript:** "The high activation energy limit of an equidiffusional flame model: rate of convergence and ε-uniform numerical approximation"

**Version reviewed:** `paper/main.tex` at commit `ba7224b` (1854 lines). Line numbers refer to that commit; the file was edited while I was reading (abstract and eq. (1) changed between `53b8e76` and `ba7224b`), so section/label references are the stable ones.

**Code reviewed:** `numerics/combustion1d/{ale,simple,solver,exact,model}.py`, `numerics/experiments/paper/*`, `numerics/results/paper/B_convergence.json`.

---

## Summary and recommendation

The paper has two parts: (i) a rate O(ε) for the singular limit, obtained by composing an inner profile with a rescaled classical solution (Section 3), and (ii) a monotone, linearly implicit scheme on a front-following band mesh, with error estimates for travelling waves (Sections 4–6) and experiments (Section 7).

**Recommendation: major revision (not acceptable in its present form).**

What I checked and found correct:

- Proposition `prop:speed` and Table `tab:speedcoef` (reproduced to all printed digits).
- Lemmas `lem:bounds`, `lem:comparison`, `an:lem:comp` (discrete comparison).
- Lemma `an:lem:profile`, Lemma `an:lem:residual` and Theorem `an:thm:main` (Eulerian estimate): I went through the inequalities line by line and found no error.
- Lemma `an:lem:ale`: the sign argument (σ ≥ μ makes the lag favourable) is valid for a uniform grid translating rigidly with a spatially constant velocity.
- There is no logical circularity between Sections 3 and 6: Lemma `an:lem:profile` does not use Section 3. The forward references in the proof of Lemma `q:lem:inner` are a presentation problem only.

What is wrong or unproved, in decreasing order of importance:

1. The lower bound of Theorem `q:thm:rate` is not proved: the subsolution construction has constants that degenerate when the strictness δ is of order ε, which is exactly how it is used.
2. The upper bound relies on constants that the lemma itself declares γ-dependent, and on regularity that is not assumed.
3. Theorem `an:thm:band` (the headline numerical result) has a gap: the barrier drifts out of the band unless σT ≲ ε, which forces h_c ≲ ε.
4. The proof of Theorem `an:thm:ale` contains an arithmetic error (a lost factor 1/ε).
5. The "ε-uniform", "time step plays no role" and "logarithmic number of nodes" claims do not follow from the theorems, and every experiment is run outside the theorems' hypotheses.
6. The implementation corollary is not proved, and three of its factual statements about the code are false (measured below).
7. The speed-error metric used in Section 7 does not test what the theorems bound; the "floor 1.6e-4 = band term" interpretation is not supported.

---

## Major issues

### M1. Subsolution construction: constants blow up when δ ~ ε. Lower bound of Theorem `q:thm:rate` unproved
**Severity: fatal for the two-sided rate as written.**

**Location.** Lemma `q:lem:inner`(b) and its proof (l. 309–314, 352–355); Proposition `q:prop:sub` (l. 467–496); "Lower bound" in the proof of Theorem `q:thm:rate` (l. 546–554); Remark item 1 (l. 559–562).

**Problem.** In the proof of Lemma `q:lem:inner`(b) the lower bound on the slope is explicit: `P ≥ sqrt(2γ)/g_2`. In Proposition `q:prop:sub`, γ = cδ, and in Theorem `q:thm:rate` δ = c·a·ε'. So:

- `Ψ_s ≥ c` holds only with c = O(√ε). At Ψ = −W the energy identity gives P² = 2(G²/2 − M − κ∫P)/G² ≤ 2γ/G² when G² = 2M + 2γ.
- The orbit passes at distance O(√γ) from the saddle (−W, 0), so `|s_*| ≤ C` is false uniformly: s_* ~ −log(1/γ)/λ.
- The sensitivities Ψ_G, Ψ_{sG} solve the linearised equation along this near-saddle passage and grow like γ^(−1/2) (equivalently ∂_G P ~ 1/P at fixed Ψ). So `|T_p| ≤ C` is not uniform either.
- The proof then takes `K ≥ K_0 + C/c`, i.e. K ≳ γ^(−1/2) (γ^(−1) with the sensitivity growth). But the profile reaches Ψ = −W − q̄ only if κ∫P < G²/2 − M ≈ γ, i.e. Kε ≲ γ.
- Together these require ε ≲ γ^(3/2) (or γ²), i.e. δ ≳ ε^(2/3) (or ε^(1/2)). With δ ~ ε, as used in Theorem `q:thm:rate`, the two requirements are incompatible for small ε.

**Why it matters.** The claim "the constant does not depend on δ ... gives the rate ε rather than √ε" is the stated novelty of Section 3. As written, the argument gives the lower bound only with rate between ε^(1/2) and ε^(2/3).

**Also.** The energy identity in l. 352 has the wrong sign. From G²Ψ'' − κΨ' = (Ψ+W)f, integrating back from P(1) = 1 gives G²P²/2 = B(Ψ;W) + (G²/2 − M) **−** κ∫_Ψ^1 P. This is why κ must be ≲ γ; with the printed "+" the constraint is hidden.

**Suggested fix.** Either a different subsolution near Ψ = −W, in which the terms to absorb are bounded by C·Ψ_s with a δ-uniform C (this needs an actual argument through the saddle region), or state the lower bound with the rate the proof gives. State the dependence of every constant on γ explicitly in Lemma `q:lem:inner`.

### M2. Supersolution construction: γ-uniformity not shown; regularity hypotheses insufficient; clipping step incorrect
**Severity: major (probably repairable).**

**Location.** Lemma `q:lem:inner`(a) (l. 291–308) and proof (l. 318–350); Definition of regular strict supersolution (l. 360–377); Proposition `q:prop:super` (l. 379–465).

(a) **Constants.** Lemma `q:lem:inner` says κ_0, C, λ_0 depend on γ. Proposition `q:prop:super` asserts constants "*not* on δ" while γ = cδ → 0. The proof uses `1−Ψ_m ≥ c(γ)`, `(Ψ+W)f(Ψ) ≥ c(γ)`, and an implicit-function bound `(Ψ_m+W)f(Ψ_m) − κ ≥ c > 0`, all of which degenerate as γ → 0. Since θ → 0 and f(Ψ_m) → 0, ∂_pθ ~ 1/f(Ψ_m) is unbounded, and so is Ψ_pp. I believe this can be repaired: Ψ_pp enters T_p only multiplied by ε, and ε/f(Ψ_m) ≲ ε/γ is bounded. But that argument is absent, and K, S, c_*, a are chosen in a circular order (K depends on S and C; κ = Kε ≤ κ_0(γ); γ depends on a; a depends on c_*).

(b) **Regularity of f.** Item (iv) ("Ψ is C² in p") requires f ∈ C² on (−∞, 1) (C^{1,1} at least). Hypothesis (F) only gives Lipschitz. It holds for f = (1−s)^+, not under (F).

(c) **Regularity of ṽ.** T_p contains ε Ψ_G (ΔG − G_t) with G = |∇ṽ|. This involves D³ṽ and ∇ṽ_t. The definition only gives ṽ ∈ C^{2,1}, and Theorem `q:thm:rate` assumes "v is C^{2,1} up to its free boundary". The claim `|∇p| + |ΔG − G_t| ≤ C in N_r` (l. 444) does not follow. The extension "by reflection in normal coordinates" loses a further derivative.

(d) **Clipping.** G = min{|∇ṽ|, G_γ(W)} is only Lipschitz, so V is not C¹_x across {|∇ṽ| = G_γ} and "it suffices to check the inequality a.e." (l. 417–419) is not justified there. The kink has the right sign only if Ψ_G ≥ 0, which is not shown.

(e) **"Clipped parameters" (l. 458–461).** The statement "if G ≠ |∇ṽ| then s ≤ −cδ_0 r/ε" is false. The margin between |∇ṽ|² on Γ and G_γ² is O(δ) = O(ε), so clipping can start at distance ~δ/K_0 ~ ε from Γ, where s = O(1). The extra term (1/ε)(|∇ṽ|² − G²)Ψ_ss is then of size C|s|e^{λ_0 s}, not (C/ε)e^{−c/ε}. It can be absorbed in the tail estimate only after showing that no clipping occurs for s ≥ −S, which couples a, c_* and S.

(f) **Far burnt region.** Outside N_r there is no lower bound on |∇ṽ|, so G may fall below g_1 and "in all cases p ∈ P_γ" fails.

(g) **Generality of the statement.** P_γ allows G as small as g_1. For W > 0 and G² < 2B(0;W) one gets 1−θ < 0; then the region {1−θ ≤ s < 0} lies outside Ω, where "V = ũ + ε²q̃" and Δṽ − ṽ_t ≤ 0 are both false. This is harmless in Theorem `q:thm:rate` (there |∇ṽ| ≈ √(2M)) but the proposition is stated in a generality in which its proof fails.

**Suggested fix.** Assume f ∈ C² on (−∞,1] and v ∈ C^{3,1} (or avoid G-dependence, e.g. by using the signed distance as in the Allen–Cahn literature). Use a smooth clipping bounded below. Restate Lemma `q:lem:inner` with γ-uniform constants and prove it.

### M3. Theorem `q:thm:rate`: smaller gaps
**Severity: minor to major.**

- l. 543–544: `v_λ − v ≤ C(1−λ)` is justified by "v is Lipschitz in x and 1/2-Hölder in t". Hölder-1/2 in t gives only (1−λ)^(1/2) = O(√ε). The conclusion is true because v_t is bounded up to Γ and v is continuous, so v is Lipschitz in t; say that.
- "C^{2,1} up to the free boundary" on [0, T*] including t = 0 requires compatibility of u_0 (C² up to ∂Ω_0, |∇u_0| = √(2M(w_0)) there). Not stated. Existence of such classical solutions under (S') is not discussed; cite the classical-solvability literature or give examples beyond the self-similar one.
- Initial ordering for the subsolution (l. 549–553) is not "as above". In the strip Ω_0 \ λ⁻¹Ω_0 one needs εΨ(ṽ/ε) ≤ u_0^ε, which uses s_* > −C (see M1) and Y_0^ε = u_0^ε + w^ε ≥ 0. The latter is not among the hypotheses.
- q̃ is fixed to 1 in the tail estimate but must be ≥ C_0 for the ordering outside Ω_0.
- Remark item 2: the "sharpness" example (travelling waves) does not satisfy the hypotheses (compact support, starshaped).
- l. 778: "Combined with Theorem `q:thm:rate` they give the convergence of the discrete solutions to the limit problem." For travelling waves Theorem `q:thm:rate` does not apply; the link is Lemma `an:lem:profile`(5).

### M4. Theorem `an:thm:band`: the barrier leaves the band
**Severity: major (the stated estimate is not proved; a corrected hypothesis destroys the headline claim).**

**Location.** Section `an:sec:band`, mesh definition (l. 1218–1229), "Truncation" (l. 1231–1244), proof (l. 1264–1317).

**Problem.** The band is fixed in the node frame, centred at y_*. The barrier is Φ(z) with z = y + X_n − c t_n + σ t_n, so its layer sits at y − y_* = −(X_n − c t_n) − σ t_n and moves relative to the nodes by up to 2σT. For the subsolution it moves the same amount towards the fresh side. The proof identifies regions in ζ = z'/ε with mesh regions ("(II_+) 2 ≤ ζ ≤ K_f/2: fresh side of the band"; "on the burnt side of the band ζ ≤ −K_b log(1/ε)"), which is legitimate only while

    2σT ≤ min(K_f − 2, K_b log(1/ε) − Z)·ε,   σ = μ + C_4(ρ_f + T h_c).

Otherwise region (II) of the barrier lies in the transition or coarse mesh, where the truncation error is O(h_c/ε²). The theorem's hypotheses (ρ_f + μ + h_c ≤ γ) do not imply this. The missing condition is μ + ρ_f + T h_c ≲ ε/T; in particular **h_c ≲ ε/T²**.

**Why it matters.** The claim "outside the band the error is O(h_c), without ε" (abstract, contribution 2, conclusions) is exactly what this condition forbids: with h_c ≲ ε the mesh is quasi-uniform at scale ε and has O(1/ε) nodes.

**Suggested fix.** Add the condition, or take a band of width ≳ σT (no longer O(ε log 1/ε)), or re-centre the barrier at each remeshing and pay for it (which again brings in the realised front error, i.e. an a posteriori statement).

### M5. Theorem `an:thm:ale`: arithmetic error in region (II)
**Severity: major (fixable).**

**Location.** l. 1173–1186.

**Problem.** q_n = q_* + C_2(h^a + σ²Δt)t_n, and the proof bounds q_n ≤ C ε²ρ_h + C εσγ. In region (II), B ≤ C ε⁻² q_n, hence B ≤ C(ρ_h + γσ/ε), not `C(ρ_h + γσ)` as printed at l. 1184. The term C_2 T σ²Δt/ε² can be absorbed by σθ only if σΔt ≲ ε², i.e. Δt(ρ_h + μ) ≲ ε², one power of ε stronger than the hypothesis of the theorem. C_3 is also redefined mid-proof ("here we also take q_* = C_3ε²(ρ_h + γμ)") with C_3 and C_4 each depending on the other; the order of choice (C_3, then C_4, then γ) must be made explicit.

**Suggested fix.** In region (I), Φ'' = −cΦ', so r_1 ≤ Cσ²Δt·Φ' can be absorbed by σΦ' when σΔt is small, and r_2 + r_3 ≤ C h^a Φ' by σ ≥ C_4ρ_h. Then q need not grow like σ²Δt and the stated hypothesis suffices. (Theorem `an:thm:band` avoids the problem through Qχ, but inherits M4.)

### M6. "ε-uniform", "the time step plays no role", "logarithmic number of nodes": not supported
**Severity: major (title, abstract, introduction, conclusions).**

**Location.** Title; abstract (l. 32–40); contribution 2 (l. 123–134); l. 1162–1165; l. 1486–1492; conclusions (l. 1828–1837).

- **Meaning of ε-uniform.** The estimates are C(ε + μ + ρ_f + h_c) with ρ_f = h_f/ε² under (F) and h_f²/ε³ under (F'), and require ρ_f ≤ γ. So h_f ≲ ε² (resp. ε^(3/2)), and the band has ≳ ε⁻¹ log(1/ε) (resp. ε^(−1/2) log) nodes. This is not parameter-uniform in the accepted sense (Miller–O'Riordan–Shishkin, which is in the .bib but not cited), and "the number of nodes grows only logarithmically in 1/ε" is true only for h_f ∝ ε, where the theorems do not apply.
- **Experiments are outside the theorems.** All runs use f = (1−s)^+ (only (F)) and h_f = ε/16, so ρ_f = 1/(16ε) = 6.25 to 62.5. The hypothesis ρ_f ≤ γ fails, and the bound C·ρ_f is larger than ‖u‖_∞. "The experiments confirm the estimates / confirm each term" (abstract, l. 1834) is therefore not accurate: they are consistent with the conjectured sharp estimate of Section `an:sec:sharp`, which is not proved.
- **Time step.** Theorem `an:thm:ale` requires Δt(ρ_h + μ) ≤ γε. At h = ε/16 this is Δt ≲ ε², the usual parabolic restriction. "Δt ≤ γ independent of ε" needs ρ_h + μ ≲ ε, i.e. h ≲ ε³ under (F). See also M5.
- **Scope.** All discrete estimates are for one-dimensional travelling waves with constant w, exact ε-wave initial data U⁰ = Φ_ε(x_j), an infinite uniform (or rigid band) mesh and spatially constant mesh velocity. The new abstract ("prove error estimates that are uniform in ε: the time step plays no role") drops all of these qualifiers. For general data the paper's own experiments show an error bΔt with b ≈ 1.35–2.6 (Fig. `fig:B4`, Fig. `fig:C1`), so reaching the O(ε) regime needs Δt ≲ ε/3 and O(1/ε) steps. The cost advantage over the adaptive Eulerian solver (ε^(−1.4) steps) is then modest and should be stated as such.

### M7. Implementation corollary (Section `an:sec:impl`) is not proved, and three of its premises are false in the code
**Severity: major.**

**Location.** l. 1333–1379; Section 4, l. 624–647.

(a) **"Experimentally μ ≲ 10⁻⁴" (l. 1345) is false.** I recorded the realised band velocity v_n = (displacement of a band node)/Δt in the default ALE run (w_0 = 0, Δt = 0.02, T = 12):

| ε | first step | t ∈ (0.1, 0.5] | t ∈ (1, 2] | after each remeshing (t > 2) |
|---|---|---|---|---|
| 0.01 | 0.288 | 1.5e-2 | 3.0e-4 | 3.3e-3 |
| 0.001 | 0.297 | 1.5e-2 | 2.7e-4 | 3.8e-3 |

(values are max |v_n − c_ε| over the interval). So μ = max_n |v_n − c_ε| ≈ 0.3 over the run, and ≈ 4e-3 even after the transient, larger than ε at ε = 0.001. With this μ the bound C(ε + μ + ...) is dominated by μ, and the hypothesis Δt·μ ≤ γε fails (0.02 × 0.3 ≫ 0.001).

(b) **μ is not controlled a priori.** "For the travelling wave the implicit choice gives v_n = c_ε + O(ρ_f + h_c)" is asserted, not proved; a bootstrap from the error bound gives only |v_n − c| ≲ (position error)/Δt. Existence of a root of the secant equation is not shown (the code falls back to the previous velocity). The corollary is an a posteriori statement and should be labelled so.

(c) **"The new band lies on the same lattice as the old one, so the transfer is exact in the band" (l. 642–646, 1351–1353) is false in the code.** `ale.py` takes `anchor = band_nodes[0]` with `band_nodes = x[diff < 1.5*h_f]`. With growth factor 1.1 the four transition cells next to the band (1.1, 1.21, 1.331, 1.4641 h_f) pass the 1.5 h_f test, so the anchor is a transition node 3.64 h_f to the left of the band. I measured the new band lattice offset from the old one: 0.349 h_f (ε = 0.01) and 0.359 h_f (ε = 0.001), at every remeshing. Linear interpolation is therefore applied inside the reaction layer at each remeshing, with error ~h_f²|Φ''|/8 ~ h_f²/ε per remeshing where the barrier has no room (a constant added in the layer costs q/ε²). The velocity spikes in (a) occur exactly at the remeshing steps.

(d) **Absorbing the remeshing error by "adding C h_c² to Q" does not work near the band edge.** Q acts through Qχ, and χ(z) ≈ C_χ(z − K_fε/2) is O(ε) at the fresh edge of the band. The interpolation error at distance d into the transition is ~((g−1)d)²|Φ''|/8, the available increment is C h_c² C_χ d, and the inequality fails for d ~ h_c/(g−1) as h_c → 0.

(e) **Rigid zone.** Section 4 says "we take R_1 ≥ max(2, 3/c)" and l. 1370–1371 "the rigid zone extends a distance ℓ_0 ≥ max(2, 3/c) beyond the band". In `ale.py`, `R1 = (max(K_burnt, K_fresh) + K_margin)*eps + 0.5*ell0` with `ell0 = max(2, 3/c)`, so the rigid radius is ≈ max(1, 1.5/c) and e^(−cR_1) ≈ e^(−1.5) ≈ 0.22: the "additional term O(Δt e^(−cR_1))" is not small. The README confirms that the default fails for w_0 = 1 (error ≈ 0.06Δt, growing as ε decreases) and that ℓ_0 ≈ 6.5/c "seems needed"; the paper reports only the tuned value ℓ_0 = 4.

(f) The far-field argument does not mention that the support of χ' must lie inside the rigid zone (it does in the code, but only because R_1 ≥ 1 + 35ε).

**Suggested fix.** Demote the corollary to a remark. Fix the anchor (use the nodes with |h − h_f| < tol). Report the realised μ honestly, or modify the algorithm (start from the ε-wave, damp v_n) so that μ is small, and prove something about it.

### M8. The speed-error experiments do not test the theorems; the "band term" interpretation is unsupported
**Severity: major.**

**Location.** Section 7.2 (l. 1516–1582), Table `tab:B1`, Figs `fig:B1`–`fig:B3`.

(a) **The rigid-mesh test (Fig. `fig:B3`) is tautological.** On a mesh translating at exactly c_ε, any discrete profile that is steady in the node frame has front speed exactly c_ε. I ran it with coarser bands (ε = 0.005, Δt = 0.02, T = 10):

| h_f | c_h − c_ε |
|---|---|
| ε/16 | +1.1e-6 |
| ε/4 | −2.2e-7 |
| ε/2 | +2.2e-7 |

The "error" does not see h_f at all. The test shows that the discrete solution locks to the mesh, not that it is accurate.

(b) **The front is dragged by the mesh, by a Δt/ε-dependent amount.** Prescribing v = c_ε + 10⁻³ on the rigid mesh (w_0 = 0):

| ε | Δt | Δt/ε | (c_h − c_ε)/(v − c_ε) |
|---|---|---|---|
| 0.001 | 0.02 | 20 | 0.91 |
| 0.001 | 0.002 | 2 | 0.50 |
| 0.01 | 0.02 | 2 | 0.49 |
| 0.01 | 0.002 | 0.2 | 0.10 |
| 0.01 | 0.0002 | 0.02 | 0.03 |

So the front-speed error is ≈ α(Δt/ε)·(v − c_ε) with α → 1 for Δt ≫ ε. This is consistent with the μ-term of Theorem `an:thm:ale` (and shows it is sharp), but it means: the ALE error does depend on Δt/ε, through μ; for Δt ≫ ε the implicit-velocity equation is nearly degenerate (d(gap)/dv ≈ −(1−α)Δt); and "front speed" in ALE mode measures mostly the mesh speed. This should be reported.

(c) **The quoted ALE error depends on the estimator and on T.** For the Table `tab:B1` entry w_0 = 0, ε = 0.01, Δt = 0.02 I reproduce +1.66e-4 with `cycle_speed`, but the least-squares slope over the same window gives +1.06e-4 (also stored as `err_lsq` in `B_convergence.json`; the ratio is ≈ 0.6 for all w_0 = 0 entries). Only 2 remeshing cycles enter the average. With T = 10 instead of 12 (ε = 0.005) I get +2.6e-4 (cycle), +9.9e-5 (lsq on [T/2, T]), +1.3e-4 (lsq on [3T/4, T]), +3.8e-5 (end points). For h_f = ε/4 the estimators give −2.2e-3, −6.6e-4, +1.1e-3, +6e-5: not even the sign is determined. Statements such as "roughly first order in h_f/ε" and "the curves for different ε coincide" (the ε = 0.01, h_f = ε/4 point is 1.54e-3 against 8.7e-4) are not robust at this level.

(d) **"The floor ≈ 1.6e-4 is the band term for h_f = ε/16" (l. 1560–1562) contradicts the paper's own computation.** Section `an:sec:sharp` reports a solvability drift "below 4e-5 for h ≤ ε/8"; `truncation_check.py` gives −2.5e-6 at h = ε/16, w_0 = 0 (opposite sign, 60 times smaller, and not decreasing between ε/8 and ε/16). The rigid runs in (a) agree with ~10⁻⁶. With implicit velocity on a rigid mesh without remeshing I get +2e-5 to +2.7e-5. The 1.6e-4 is therefore produced by the remeshing cycle / non-rigid zone / estimator, not by the band discretisation. The "O(h/ε) from the kink of f" explanation in Section `an:sec:sharp` is likewise unsupported by the numbers printed there.

(e) **No experiment reports the quantity the theorems bound**, max_j |U_j^n − Φ_ε(x_j^n − c_ε t_n)|, although Φ_ε is available from the shooting code. Report it against Δt, h_f/ε, h_f/ε² and h_c, starting from U⁰ = Φ_ε as in the theorems.

(f) There is a small but systematic Δt-trend in the ALE rows for w_0 = 0 (1.83, 1.67, 1.55 ×10⁻⁴ for Δt = 0.02, 0.01, 0.005 at ε = 0.001; the Δt = 0.01 row is omitted from Table `tab:B1`). "Does not depend on Δt" should be "depends weakly".

### M9. Code versus the scheme that is analysed
**Severity: major for reproducibility, minor individually.**

- **Which solver produced what.** Table `tab:tw`, Figs `fig:tw`, `fig:instability`, `fig:limit`, Table `tab:limit` and half of Section 7.5–7.6 come from the adaptive BDF solver with PCHIP transfer (`solver.py`): neither linearly implicit nor monotone, and covered by no result of Sections 5–6. Say so in each caption.
- **Nominal versus actual Δt in the Eulerian monotone runs.** `simple.py` sets `h = interval/nsteps` with `interval = 7.5ε/v_est`, so the actual step is ≤ the nominal one and varies with the measured speed. In my reruns the mean step was 0.0091 for nominal 0.01 (ε = 0.001) and 0.0181 for nominal 0.02 (ε = 0.002). The "function of Δt/ε alone" collapse in Fig. `fig:B1` is affected by this.
- **Cost of the implicit velocity.** "Two to four linear solves" (l. 636): each trial velocity costs two tridiagonal solves and a mesh/matrix rebuild, at least two trials are made, up to ten. "ALE: 400 steps" in Table `tab:A3` is not comparable to BDF steps without this factor.
- **Band width.** Abstract/intro/Theorem `an:thm:band`: O(ε log(1/ε)) with K_b log(1/ε) on the burnt side. Code and Section 7: fixed 55ε. Conclusions: "band of width O(ε)".
- **Initial data.** Section 7.2 does not say that the convergence runs start from A(1 − e^(−c_ε(x−2)))^+ with c_ε from the shooting code (README), i.e. neither the limit wave nor the ε-wave of the theorems.
- **Table `tab:simple`** has no generating script in the repository or in the README (see Reproduction).
- **What does match.** Reaction lagged as f_ε(Uⁿ)Y^{n+1}; two tridiagonal solves per step; M-matrix/Péclet check on the off-diagonals (correct condition |ẋ_i|h ≤ 2 on non-uniform meshes); centred three-point convection; `np.interp` (linear) at remeshing in `ale.py` and `simple.py`; W = Y − U discrete-caloric.

### M10. Novelty claims and missing literature
**Severity: major for the introduction.** (Titles from memory; please check details.)

- **Proposition `prop:speed`.** Two-term expansions of the flame speed in the inverse activation energy are classical: Bush and Fendell, *Combust. Sci. Technol.* 1 (1970), and the monographs of Williams and of Zeldovich et al. The specific formula with w_0 may be new; "the first-order correction ... in closed form" should be positioned against this.
- **"First rate for this singular limit" / technique.** Barriers of the form (inner profile)∘(modified distance) giving optimal O(ε) interface estimates are standard for Allen–Cahn-type limits: X. Chen (JDE 1992), de Mottoni–Schatzman, Bellettini–Paolini (1995), Alfaro–Hilhorst–Matano (JDE 2008). The claim may be true for this particular problem, but the method is not new and these must be cited.
- **Numerical analysis of singular limits with ε–h–Δt coupling.** Nochetto–Paolini–Verdi (optimal interface error estimates for mean curvature flow via Allen–Cahn/double obstacle), Feng–Prohl (Numer. Math. 2003). This is the natural SINUM context and is absent.
- **Layer-adapted meshes and parameter-uniform convergence.** Miller–O'Riordan–Shishkin and Kopteva–O'Riordan are in the .bib but not cited; Roos–Stynes–Tobiska, Linß.
- **Implicit mesh velocity so that the front stays at a reference point** is the freezing method / phase condition of Beyn–Thümmler (SIADS 2004) and Rowley et al.; moving-mesh combustion computations (Dorfi–Drury, Petzold, Larrouturou) are also relevant.
- **"No numerical study of the limit ε → 0 is available"** is a strong negative claim; the .bib itself notes an unverified numerical appendix in Caffarelli–Lee–Mellet. Several .bib entries still carry "unverified" check notes (authors, journal versions) that must be resolved before submission.

---

## Minor issues

1. **Proposition `prop:speed`** has only a "sketch of proof" (l. 186). The missing step (φ' → √(2G) uniformly in φ) follows from the two-sided energy bound in the proof of Lemma `an:lem:profile`(1); give it, and move Lemma `an:lem:profile` before its first use (Sections 2 and 3 both rely on it).
2. **Forward references** in the proof of Lemma `q:lem:inner` ("as in the proof of Lemma `an:lem:profile`", l. 325, 344) to a lemma about a different ODE (G = 1, constant W, damping σ). Not circular, but the reader has to redo the proof with parameters.
3. **Hypotheses.** (F) is stated three times with small variations (l. 59–60 without monotonicity, l. 268–269, l. 681–684). Theorem `q:thm:rate` uses "the setting of [FBW04]" and "(2.2)" without restating them.
4. l. 607: "O(L/h_c) + O(1/ε · ε/h_f) + O(log(h_c/h_f))" — the middle term reads as O(1/h_f); it should be (band width)/h_f.
5. l. 100–102: "a uniform mesh needs O(ε^(−N)) cells" — in the regime the theorems need (h ≲ ε²) it is ε^(−2N).
6. Remark after Theorem `an:thm:band`, item 2: "it should be h_c², which is what the experiments show" rests on two points per ε (5e-5 and 1e-5 at h_c = 0.04, 0.02, ratio 5; the ε = 0.001 column gives 6.7e-5, 2.6e-5, ratio 2.5) at the noise level of the measurement.
7. Fig. `fig:B3`: "no trend for Δt from 0.1 to 0.005" — for w_0 = 1, ε = 0.005 the stored values go from −2.7e-5 (Δt = 0.1) to +2.2e-6 (Δt = 0.005). See M8(a) for the relevance of this test.
8. Section 7.1: cells behind the flame "stretched up to 6 h_c": then h_c is not the maximal mesh size of Theorem `an:thm:band`, and the ratio condition g ≤ 2 should be checked for the distorted mesh.
9. Section 7.6, quenching: "r_* = 0.715 ± 0.001 for every ε" — state Δt and that both the threshold and "shrinks to r ≈ 0.52" are for this particular profile.
10. Section `sec:limitations` item 2 admits an O(Δt v/ε) lag near r = 0; this qualifies the "extinction times agree to four digits" statement (ALE, ε = 0.01: 0.72459 against 0.725090 is three digits).
11. Sharp-estimate subsection: "ψ* is bounded" ahead of the front, but it does not decay; the solvability integrals converge only because of Φ'. Say which weighted space is used.
12. Notation: μ is both the mesh-velocity error and f(−w_0)/4 (renamed μ_0 only later); κ is the damping in Section 3 and the constant (3√3−2)/5 in Section 2; a is f(−w_0), an exponent (1 or 2), and the rescaling constant; C_3 is redefined mid-proof in Theorem `an:thm:ale`.
13. Table `tab:B1` omits the Eulerian w_0 = 1 rows and the default-ℓ_0 ALE rows (errors up to 1.2e-3, O(Δt), growing as ε ↓) that are in the README. Include them: they are the honest picture of the method with its default parameters.

---

## Reproduction

Environment: `/home/user/combustion/.venv/bin/python`, `PYTHONPATH=/home/user/combustion/numerics`, `OMP_NUM_THREADS=1`, run from `numerics/`. Total compute about 3 minutes. Scratch scripts are in `/tmp/claude-0/-home-user-combustion/c5917af1-f076-5e7f-b209-4cc1d3f058ea/scratchpad/referee/` (`t1.py` … `t8.py`).

| # | What | Command | Claimed | Obtained | Verdict |
|---|---|---|---|---|---|
| 1 | Prop. `prop:speed`: I/√(2M) = κ(1+w_0) | `python t1.py` | κ = 0.6392 | 0.63923 (quadrature) for w_0 = −0.5, 0, 1 | reproduced |
| 2 | Table `tab:speedcoef`, (c_ε − c_0)/ε | `python t1.py` | +0.139876, +0.209515, −0.448230 (ε = 1e-2); +0.138982, +0.208412, −0.454072 (1e-3) | identical to 6 digits | reproduced |
| 3 | Table `tab:B1`, ALE, w_0 = 0, ε = 0.01, Δt = 0.02 / 0.005 | `python t2.py 0 0.01 0.02 None`, `... 0.005 None` | +1.66e-4 / +1.42e-4 | +1.663e-4 / +1.418e-4 (cycle); +1.06e-4 / +7.9e-5 (lsq) | reproduced with the paper's estimator only |
| 4 | Table `tab:B1`, Eulerian, w_0 = 0 | `python t5.py simple 0 0.001 0.01`, `simple 0 0.01 0.01`, `simple 0 0.001 0.001`, `simple 0 0.002 0.02` | −3.88e-1, −1.44e-1, −0.144, (README) −3.87e-1 | −3.884e-1, −1.443e-1, −1.442e-1, −3.875e-1 | reproduced; actual mean Δt 0.0091, 0.0097, 0.00097, 0.0181 |
| 5 | Fig. `fig:B3`, rigid mesh, ε = 0.001 | `python t5.py rigid 0 0.001 0.1` (and 0.005; w_0 = 1) | < 3e-5, no Δt trend | −5.1e-7, +2.9e-8; −3.5e-6, +1.0e-7 | reproduced, but see M8(a) |
| 6 | Rigid mesh, coarser band, and wrong velocity | `python t6.py`, `python t7.py` | (not in paper) | tables in M8(a), M8(b) | new |
| 7 | Table `tab:simple` (scheme (S), Δt ∝ ε²) | `python t5.py simple 0 0.02 0.0004`, `... 0.02 0.0016`, `... 0.01 0.0001` | −4.60e-3, −1.85e-2, −2.21e-3 | −4.84e-3, −1.857e-2, −2.39e-3 | **not reproduced exactly** (5–8 % off in two of three); no script; my setup: T = 12, c_ε-data, lsq on [T/2, T] |
| 8 | Section `an:sec:sharp`, solvability drift | `python -m experiments.truncation_check` (13 s) | −0.2351, −0.2345, −0.2342; b ≈ 0.085, 0.47; spatial part < 4e-5 | same; spatial drift −2.5e-6 (w_0 = 0, h = ε/16), −3.9e-5 (w_0 = 1, ε = 0.01, h = ε/8) | reproduced; contradicts the "1.6e-4 band term", see M8(d) |
| 9 | "experimentally μ ≲ 1e-4" | `python t3.py 0 0.01 0.02 12`, `python t3.py 0 0.001 0.02 12` | ≲ 1e-4 | 0.29 (first step), 1.5e-2 (t < 0.5), 3–4e-3 after each remeshing | **not reproduced** |
| 10 | "new band on the same lattice" | `python t4.py` | exact transfer in the band | offset 0.349 h_f (ε = 0.01), 0.359 h_f (ε = 0.001) | **false in the code** |
| 11 | ALE floor "1.6e-4" | `python t8.py` (T = 10) | 1.6e-4, independent of Δt, ε | cycle +2.6e-4 to +2.9e-4; lsq +0.8e-4 to +1.3e-4; end points +2e-5 to +4e-5; rigid + implicit v: +2e-5 to +2.7e-5 | **estimator- and T-dependent** |
| 12 | Fig. `fig:B4` constants | read `results/paper/B_convergence.json` | 1.35Δt (N = 2), 1.55Δt (N = 3) at ε = 0.001 | 1.34–1.44, 1.53–1.64; ε = 0.01: 1.5–2.0 | consistent with stored data (not rerun) |
| 13 | Table `tab:tw` | compared with `numerics/results_tw.log` | e.g. +7.4e-5 (w_0 = 0, ε = 0.01) | 0.579519 − 0.579445 = 7.4e-5 | consistent with stored log (not rerun; BDF solver) |

Not rerun (cost or scope): group A (meshes), group C (limit comparison, self-similar convergence, quenching), the BDF travelling-wave table, the instability figure. All figure files referenced in `main.tex` exist; all `\ref`/`\cite` keys resolve.
