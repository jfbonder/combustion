# Referee report, second round (SINUM)

**Manuscript:** "The high activation energy limit of an equidiffusional flame model: rate of convergence and a front-following numerical method"

**Version reviewed:** `paper/main.tex` at commit `664e555` (2133 lines; line numbers refer to that commit), `paper/references.bib`.

**Code reviewed:** `numerics/combustion1d/{ale,simple,solver}.py`, `numerics/experiments/paper/{tw_error,tab_simple}.py`, `numerics/experiments/paper/README.md`, `numerics/results/paper/D_tw_error.json`.

**First-round report:** `paper/referee_report.md` (issues M1–M10).

---

## Summary and recommendation

**Recommendation: major revision, of much narrower scope than in round 1. I found no fatal error.**

The revision is substantial and mostly successful.

- **Section 3.** I went through Lemma `q:lem:inner`, Lemma `q:lem:ext`, Propositions `q:prop:super` and `q:prop:sub` and Theorem `q:thm:rate` line by line. The fatal problem of round 1 (M1) is removed: with G_±(W)² − 2M(W) = ±2γ independent of W, the W-derivative of the energy has no O(1) term at the saddle, and the bounds |P_W| + |P_WW| ≤ CP are γ-uniform. I also checked the three bounds of (a)(ii) and (b)(ii) numerically for f = (1−s)^+, γ from 10⁻² to 10⁻⁹ (see Reproduction). I believe the theorem as stated is true. What remains are gaps in the write-up (N4) and one new unproved claim (N3).
- **Proposition `prop:speed`.** The proof is complete and correct.
- **Section 6.** The estimates are true as far as I can tell, but the two rewritten proofs do not close as written: a circular choice of constants in Theorem `an:thm:ale` (N1), and a too-small q_b plus a sign error in Theorem `an:thm:band` (N2). Both are repairable; I give the repair.
- **Section 7.** The code fixes are real and every number I reran is reproduced. Two things are not reported correctly: the Eulerian "constants" depend on the time window by 20–30 % (N5), and the realised μ grows without bound as Δt → 0 (N6).

Issues in decreasing order of importance:

1. Theorem `an:thm:ale`: C_3 and C_4 each need to be large in terms of the other; the displayed inequality in region (III) is false for C_4 > 6 (N1).
2. Theorem `an:thm:band`: region (III) fails with q_b = C_3ε²ρ_f; the transport of χ has the unfavourable sign for the subsolution; the hypotheses imply Δt ≲ ε/h_c and a band of width ≳ C_T h_c, which the introduction does not say (N2).
3. "The logarithm disappears for well-prepared data" (introduction, Remark `q:rem:rate`(2)) is not proved, and the data described violate a hypothesis of the theorem (N3).
4. Section 3 bookkeeping: order of constants (Λ, K_1, Λ_1), dependence of c_* on q̃, a false regularity statement, a reversed inequality, the extension lemma (N4).
5. Section 7: speeds measured during a slow transient are presented as constants (N5); μ over all steps diverges as Δt → 0 (N6); pinning is measured at one lattice phase (N7).

---

## Status of the round-1 major issues

| # | Issue | Status | Reason |
|---|---|---|---|
| M1 | Subsolution: constants blow up for δ ~ ε | **Resolved** | New profile (b): P² ≥ c(γ+Θ²), so Θ_s ≥ c√γ and s_* ≥ −C log(1/γ). ∂_W of the energy is O(Θ²) + O(κ), hence \|P_W\| ≤ CP uniformly (needs κ ≤ cγ, which is assumed). The (1+\|s\|) factors are absorbed by the margin Λ_1ε\|s\| from Lemma `q:lem:ext`. For Θ < 0 the proof uses Θ_s ≥ c√γ against O(δ + ε log(1/ε)); correct. Sign of the energy identity corrected. |
| M2 | Supersolution: γ-uniformity, regularity, clipping | **Resolved in substance; write-up gaps** | (a) uniformity proved, but see N4 for the order of constants. (b) f ∈ C²((−∞,1]) now assumed. (c),(d),(e) gone: no G-dependence, no clipping. (f) far region handled by the region s < −ε^(−1/2). (g) see N4(f). |
| M3 | Theorem `q:thm:rate`: smaller gaps | **Mostly resolved** | Lipschitz in t: fixed. Y_0^ε ≥ 0: now assumed. q̃ = max{1, C_0}: fixed. Initial ordering for the lower bound: now a correct case analysis. Not resolved: the sharpness example is still a travelling wave, outside the class of the theorem (Remark item 3, and the unqualified "The order ε cannot be improved" at l. 131); compatibility of u_0 at t = 0 is only implicit in "C^{2,1} up to the free boundary on [0,T*]". New: N3. |
| M4 | Band theorem: barrier leaves the band | **Partially resolved** | The margin hypothesis ℓ ≥ C_T(μ+ρ_f+h_c) is the correct repair of the stated estimate. Its consequences are not drawn, and the proof has two further defects (N2). |
| M5 | Theorem `an:thm:ale`: lost factor 1/ε | **Partially resolved** | The 1/ε error is gone (constant q, r_1 absorbed by σΦ'). The new proof has a circular choice of constants (N1). |
| M6 | "ε-uniform", "time step plays no role", "logarithmic nodes" | **Mostly resolved** | Title changed; h_f ≲ ε² and "experiments outside the hypotheses" are stated in the introduction and conclusions. Remaining: the Δt restriction is not stated next to "the time step does not appear" (N2(d)); l. 839–840 states "for general data it is O(Δt) with a constant independent of ε" as a fact (it is an observation on two examples; the README says it "holds only roughly", 1.4–1.9, and that aε + bΔt "no longer describes the data well" for C1). |
| M7 | Implementation corollary | **Resolved** | Corollary removed; Section `an:sec:impl` says what is not proved. Anchor fix verified: lattice offset at remeshing < 3·10⁻¹¹ h_f in 3 runs (was 0.35 h_f). First step from Φ_ε: \|v_1 − c_ε\| = 5.3·10⁻⁴ at Δt = 0.02 (was 0.29). Rigid zone is band + ℓ_0 in the code, as in the text. Remaining: N6. |
| M8 | Speed-error experiments | **Resolved, with a new problem** | max\|U − Φ_ε\| from U⁰ = Φ_ε is now reported (D1–D5); the tautological rigid-mesh test is acknowledged; drag is reported (reproduced); the three speed estimators agree (reproduced); the "floor" is explained by lattice pinning, which I reproduced. New: N5, N7. |
| M9 | Code vs analysed scheme | **Mostly resolved** | Solver named per figure; `simple.py` uses the nominal Δt (verified: `dts` = 1.0·10⁻⁴); cost of the implicit velocity stated (3 trials = 6 solves per step, verified); `tab_simple.py` exists and reproduces the table. Remaining: band width (minor 6); Section 4 still presents the BDF/PCHIP solver as "the" method. |
| M10 | Novelty and literature | **Mostly resolved** | Allen–Cahn barriers, NPV, Feng–Prohl, Beyn–Thümmler, Bush–Fendell, layer-adapted meshes, Baconneau–Lunardi are cited and the claims are positioned. Remaining: more than 30 entries of `references.bib` still carry `checknote = {unverified: ...}` (one with "Sun, [first name unverified] and others"); literature for the new pinning claims is missing (N7). |

---

## New major issues

### N1. Theorem `an:thm:ale`: circular choice of C_3 and C_4; the inequality in region (III) is false
**Severity: major (fixable; the statement is true).**

**Location.** l. 1358–1379.

**Problem.** The proof takes q_* = C_3ε²(ρ_h + γμ), σ = μ + C_4ρ_h, "with C_3, then C_4, then γ". In region (III) it claims

    (7μ_0/4)·C_3(ρ_h + γμ) − Cρ_h − Cγ(1+C_4)(μ + C_4ρ_h) ≥ 0   if C_3 ≥ 4C/μ_0 and γ(1+C_4)C_4 ≤ 1.

The coefficient of γμ is (7μ_0/4)C_3 − C(1+C_4). With C_3 = 4C/μ_0 this is C(6 − C_4), negative for C_4 > 6. Take ρ_h → 0 with μ fixed: the left side is negative. So (III) needs C_3 ≥ (4C/7μ_0)(1+C_4), while (II) needs C_4 ≥ 4C(1+C_3)/θ. Both hold only if 16C²/(7μ_0θ) < 1, which is not available (θ = θ(Z) is small).

**Why it matters.** This is the theorem that carries the claim "no Δt in the error". In round 1 the same proof was criticised for "C_3 and C_4 each depending on the other"; the dependence is still there.

**Suggested fix.** Let q absorb r_1 directly: q_* = C_3(ε²ρ_h + εσ²Δt). Then

- (III): (7μ_0/4)C_3(ρ_h + σ²Δt/ε) ≥ C(ρ_h + σ²Δt/ε) for C_3 ≥ 4C/(7μ_0), with no reference to C_4;
- (II): B ≤ Cq_*/ε² = CC_3ρ_h + CC_3(σΔt/ε)σ, absorbed by θσ if C_4θ ≥ 4C(1+C_3) and then (σΔt/ε) ≤ (1+C_4)γ ≤ θ/(4CC_3);
- the order C_3, C_4, γ is then consistent, and q_* ≤ C(ρ_h+μ)ε.

**Smaller points in the same proof.**

- Lemma `an:lem:ale` uses "sup_I" with I of Lemma `an:lem:residual`, I = [z'−h, z'+max(h, cΔt)]. In the ALE case the Taylor points lie in [z'−2σΔt, z'], which is not in I when 2σΔt > h. Redefine I = [z'−max(h, 2σΔt), z'+h].
- For the subsolution in region (I) one needs V^n ≥ ε, i.e. q_* ≤ Φ(2ε) − ε ≈ 2√(2M) ε. This is said in Theorem `an:thm:band` but not here.

### N2. Theorem `an:thm:band`: region (III), sign for the subsolution, and what the hypotheses imply
**Severity: major (items (a),(b) fixable; (c),(d) concern the claims).**

**Location.** l. 1437–1516; introduction l. 132–148; abstract.

(a) **Region (III) does not close with q_b = C_3ε²ρ_f.** The proof says "as in Theorem `an:thm:ale`: the stability term dominates B and r if C_3 is large". The stability term is 2μ_0C_3ρ_f. The residual r in (III) contains r_1 ≤ Cσ²Δt/ε ≤ Cγ'σ with σ = μ + C_4(ρ_f + Th_c), and, in the coarse burnt part, Ch_cε⁴. Neither μ nor h_c is controlled by ρ_f (take h_f → 0 with μ, h_c fixed). In Theorem `an:thm:ale` the term γμ was included in q_* for this reason; here it was dropped. Fix: q_b = C_3(ε²(ρ_f + h_c) + εσ²Δt), as in N1; then (II) needs C_4 large in terms of C_3 and the subsolution condition q_b ≤ (Ac_0/16)K_fε still holds.

(b) **Subsolution: the transport of χ has the unfavourable sign.** For W = Φ(z−2σt) − q_b − Q_nχ(z−2σt) the argument of χ decreases at a fixed node (z' ≤ z_old), so

    −[Q_{n+1}χ(z') − Q_nχ(z_old)]/Δt = −Q̇_nχ(z') + Q_n(χ(z_old) − χ(z'))/Δt,

and the last term is ≥ 0, up to 2σC_χQ_n on supp χ'. For a subsolution this is the wrong sign. "All signs are reversed as in Theorem `an:thm:ale`" is therefore not correct. The term is absorbed by σθ_1 once 4C_χ sup Q ≤ θ_1, which holds for γ small; say so.

(c) **The margin hypothesis changes what the theorem says about the mesh.** ℓ ≥ C_T(μ+ρ_f+h_c) with C_T ≥ 2T·C_4·T makes the band at least C_T h_c wide on each side. With h_c fixed this is not O(ε log(1/ε)), and the band has ≳ C_T h_c/h_f nodes. A band of width O(ε log(1/ε)) requires h_c + μ ≲ ε log(1/ε), and then the coarse mesh has ≳ 1/(ε log(1/ε)) nodes. The introduction (l. 133–143) still describes "a fine band of width O(ε log(1/ε))", a node count that "grows only logarithmically", and the estimate C_T(ε+μ+ρ_f+h_c) in the same paragraph. These three statements do not hold together under the hypotheses of the theorem. The code uses a margin of 15ε (0.015 at ε = 0.001) against h_c = 0.01–0.04, so it does not satisfy the margin hypothesis either; Section `an:sec:impl` lists three differences between code and theorem and this is a fourth (together with the burnt-side width 35ε instead of K_b ε log(1/ε) with K_b ≥ 12/κ_0, i.e. 55ε at ε = 0.01 and 83ε at ε = 0.001 for w_0 = 0).

(d) **The Δt restriction.** The hypothesis Δt(ρ_f+μ+h_c) ≤ γε gives Δt ≤ γε/h_c for a fixed coarse mesh: Δt = O(ε), hence O(1/ε) steps. The abstract ("error estimates in which the time step does not appear") and l. 142 should state this restriction; "the number of ALE steps does not depend on ε" (l. 1686) is an experimental statement at Δt = 0.01, outside this hypothesis for small ε.

### N3. "The logarithm disappears for well-prepared data" is not proved
**Severity: major for the claim in the introduction (l. 125); minor if the claim is withdrawn or restricted to the lower bound.**

**Location.** Remark `q:rem:rate`(2), l. 733–742.

**Problem.** The remark assumes Y_0^ε ≥ cε·exp(−β dist(x,Ω_0)/ε) with β small and concludes "the rate is O(ε+η_ε) on both sides".

- The lower bound with A = 0 does follow (I checked: Θ ≤ (1+W)e^{c(s−1)}, s ≤ −c'(ρ + caε')/ε, a large, β ≤ cc').
- The upper bound of Theorem `q:thm:rate` assumes u_0^ε ≤ −w^ε + C_0ε² outside Ω_0, i.e. Y_0^ε ≤ C_0ε² there. Data with Y_0^ε ≥ cε just outside Ω_0 violate it. So the theorem does not apply to the data of the remark.
- Repairing the upper bound needs εΘ^+(s) ≥ Y_0^ε outside Ω_0, i.e. an *upper* bound Y_0^ε ≤ Cε·exp(−β'ρ/ε) with β' at least the decay rate of Θ^+ (in the variable ρ/ε). With the constants of Lemma `q:lem:inner` (generic c, C, λ_0) one needs β small for the lower bound and β' large for the upper bound. For data with exactly the travelling-wave tail this would require decay rates of Θ^± sharp to within O(δ), which the lemma does not provide.

**Suggested fix.** Either prove a version with sharp tail exponents, or state only: "for such data the lower bound holds without the logarithm", and remove "the logarithm disappears for well-prepared data" from the introduction.

### N4. Section 3: gaps in the write-up
**Severity: minor individually; together they must be fixed, because the non-circular order of constants is the point of the section.**

(a) **Λ, K_1 and Λ_1 (l. 545–556, 566–578, 595–596).** The proof needs Λ_1 ≥ C_T, where C_T bounds 2|∇ṽ||∇W| in T and so depends on sup|∇ṽ|. It obtains Λ_1 = (Λ−C)/K_1 from d ≥ ε|s|/K_1, where K_1 is the C^{1,1} norm of ṽ and depends on Λ (Lemma `q:lem:ext`: "K_1 depending only on δ_0, r, K_0, Λ"). As written, "Λ_1 as large as we wish by taking Λ large" is circular. It is true, because in the region d ≤ r_1 one has |ṽ^±| ≤ (5/4)|Eũ| and |∇ṽ^±| ≤ (3/2)|∇Eũ|, both independent of Λ (only D²ṽ and r_1 depend on Λ). State these two bounds in Lemma `q:lem:ext` and use them.

(b) **c_* depends on q̃.** In the region Θ > θ_* the proof takes K ≥ K_1 + C + Cq̃/c_0, and then c_* ≤ c·c_1/K. Proposition `q:prop:super` says c_*, δ_* depend only on f, δ_1, W̄, C_W, δ_0, r, K_0, C_1 and allows any q̃ ≥ 1. Add q̃ to the list (harmless in the theorem, where q̃ = max{1, C_0}).

(c) **"Θ is C² in W" (Lemma `q:lem:inner`(a)(ii)) is false across s = s_m(W).** Θ_W is continuous there (I agree with the computation −s_W P = 1), but Θ_WW and Θ_sW jump (Θ_sW = P_Θ + P_W(Θ_m) on the left, 0 on the right). Θ^+ is C^{1,1} in (s,W), which is what the a.e. argument needs; one must add that {ṽ/ε = s_m(W)} has measure zero (true: |∇ṽ|/ε ≥ c/ε against |s_m'||∇W| ≤ C; I find Θ_{m,W} bounded, numerically d s_m/dW ∈ [0.001, 0.5]). Note also that the bound Θ_{m,W} ≤ C/γ used for the γ^(−1) term is pessimistic: F_W(Θ_m) = ∫_{Θ_m}^{1+W} t f'(t−W) dt − κJ_W and F_Θ = b(Θ_m) − κ give |Θ_{m,W}| ≤ C, so the γ^(−1) in the Θ_WW bound is not needed (numerically the ratio |Θ_WW|/((1+|s|)²Θ_s) stays below 17 for γ → 0).

(d) **Reversed inequality (l. 449).** "Θ(0) ≥ c (since s_m = ∫ dθ/P ≤ (1/c) log(Θ_m/Θ(0)))" gives an upper bound for Θ(0). The lower bound follows from P ≤ CΘ: s_m ≥ (1/C) log(Θ_m/Θ(0)) and s_m < 1.

(e) **Lemma `q:lem:ext`.** "Reflection in normal coordinates" does not give a C^{2,1} extension: for Γ_t of class C² the normal-coordinate map is only C¹, so the reflected function is C¹, not C^{1,1} (the first-round report made the same point). Use a higher-order reflection after flattening Γ as a graph x_N = g(x',t) with g ∈ C^{2,1}; this requires stating the regularity of Γ in Definition `q:def:strict`. Also "d ≤ K_1ε^{1/2}" at l. 573 should be d ≤ ε^{1/2}/c.

(f) **Definition `q:def:strict` is too weak for the generality of the propositions.** Nothing prevents two components of Γ_t from being closer than r_1 across the burnt region; then the closest point y and the extension are not well defined. Assume a uniform exterior tubular neighbourhood (true for the starshaped Ω of the theorem). Similarly "d ≤ Cε/δ_0 in Ω where s ≤ s_m" needs ũ ≥ cδ_0 d for δ_0/K_0 ≤ d ≤ r, which is not in the definition.

(g) **Comparison paragraph (l. 352–357).** The class is "C^{1,1} in x, Lipschitz in t, inequality a.e."; the maximum in Proposition `q:prop:sub` is not in this class (kink at s = s_*). The conclusion is right (at s = s_* the maximum is touched from above only where the smooth subsolution −w−ε²q̃ is), but say in which sense (viscosity, or comparison on {s > s_*} with boundary values).

### N5. Section 7: the Eulerian "constants" are values taken during a slow transient
**Severity: major for the quantitative statements of Section 7.2 and Section `an:sec:sharp` (the qualitative conclusion, error ∝ Δt/ε, stands).**

**Location.** l. 943–945 and Table `tab:simple`; l. 1603–1608; l. 1799–1803 and Fig. `fig:D5`; Table `tab:B1`.

**Evidence.** Scheme (S), w_0 = 0, ε = 0.02, Δt = 4ε² (Δt/ε = 0.08), T = 12, speed error by least squares on successive windows:

| window | [1,3] | [3,6] | [6,9] | [9,12] |
|---|---|---|---|---|
| from Φ_ε | −2.354e-2 | −2.012e-2 | −1.891e-2 | −1.841e-2 |
| from A(1−e^{−c_ε(x−2)})^+ | −2.355e-2 | −2.012e-2 | −1.891e-2 | −1.841e-2 |

The two data give the same numbers; the window changes the result by 28 %, and it is still moving by 2.6 % between the last two windows. The paper's own numbers show the same thing:

- Table `tab:simple` (T = 12, window [6,12]) gives −(0.39–0.42)c_εΔt/ε; D5 (T = 4, window [2,4]) gives "≈ 0.46 Δt/ε". The README attributes the difference to the datum ("−0.44 from Φ_ε"); it is the window.
- At Δt/ε = 1: −0.149 in Table `tab:B1` (T = 12) and 0.30·c_ε = −0.174 in D5 (T = 4).
- "δc/(Δt/ε) = −0.2351 ... in agreement with the speed error ≈ −0.23Δt/ε measured" (l. 1605–1607): the measured coefficient is −0.294, −0.251, −0.236, −0.230 on the four windows above, so it passes through the predicted value and continues.

The ALE runs have the same slow relaxation: from Φ_ε, w_0 = 0, ε = 0.01, Δt = 0.02, the mean of v_n − c_ε is 1.06e-4 on (1,1.5], 8.35e-5 on (3.5,4], 7.97e-5 on (7.5,8], then 8.4e-5 at t = 24. So "μ_1 ≈ 1.1·10⁻⁴" is the value at t = 1 of a decaying quantity, and D3 (T = 4 and 2, mean over [T/2,T]) overestimates the pinned speed by about 10 % (8.7e-5 against 8.0e-5 at θ = 0.37).

**Why it matters.** Three different "constants" for the same scheme appear in the paper (0.23, 0.46·c_ε = 0.27, and 0.149 at Δt = ε), and the agreement with the solvability prediction is quoted to three digits.

**Suggested fix.** Report the time window with every speed; show v(t) or s(t) for one case; quote the solvability comparison to the accuracy the data support (about 5 %), or run long enough for the speed to settle (the relaxation time is several units of 1/c²).

### N6. The realised μ = max_n |v_n − c_ε| grows as Δt → 0
**Severity: major for Section `an:sec:impl` ("Theorem `an:thm:band` holds with μ = max_n|v_n − c_ε|"), minor for the experiments.**

**Evidence.** From Φ_ε, w_0 = 0, I obtain μ = 2.7e-4, 5.7e-4, 1.45e-3 (ε = 0.01) and 2.8e-4, 5.8e-4, 1.59e-3 (ε = 0.001) for Δt = 0.1, 0.02, 0.0025; the maximum is at the first step. The README states the law (≈ Δt^(−0.45)). The paper reports only the range "3·10⁻⁴ to 1.6·10⁻³" (l. 1736) and not the trend.

**Why it matters.** The a posteriori bound C(ε+μ+ρ_f+h_c) with this μ gets worse as Δt decreases, while the measured error does not change. The max norm is the wrong quantity.

**Suggested fix.** The barrier argument works with a time-dependent drift σ_n ≥ |v_n − c_ε| + C_4ρ (Lemma `an:lem:residual` is already written with σ_n). Then the error bound contains Δt·Σ_n|v_n − c_ε| instead of T·max_n|v_n − c_ε|, and the restriction becomes σ_nΔt ≲ ε for each n. In the runs above Δt·Σ|v_n − c_ε| = max|s_n| = 7.7–8.3e-4 for every Δt. State the theorem in that form.

### N7. Lattice pinning: scope of the evidence and missing literature
**Severity: minor to major (it is now the main explanation of the ALE error).**

- I reproduce the phenomenon: the front stays at phase 0.3700 for 1200 steps and 10 remeshings, and the speed error depends on the phase (−1.9e-4 at θ = 1/32, +4.6e-5 at 9/32, +1.1e-4 at 17/32, +5e-7 at 25/32; ε = 0.01, h_f = ε/16).
- Consequently E "≈ 4.5·10⁻⁴, independent of Δt and ε" (l. 1731) is the value at the single phase θ = 0.37 and at T = 8. At θ = 25/32 I get E = 8·10⁻⁶; at θ = 1/32, E = 5.6·10⁻⁴ already at T = 4; and E grows linearly in T (s = 2.1·10⁻³ at T = 24). Report max over θ, or E_s and the speed error separately, and say that the phase is set by where the initial front falls in the lattice.
- "g has zero mean over θ" (l. 1765) is a numerical observation: the means are 1–3 % of the amplitude at three resolutions with 16 phases. Say so.
- "Second order in h_f/ε" (Fig. `fig:D2`) is measured at one phase; it is the scaling of the amplitude only if the shape of g does not change. D3 supports this at three resolutions; say that the two figures must be read together.
- "as in Theorem `an:thm:band`" (l. 1750): in the theorem ε enters through h_f/ε² or h_f²/ε³, not through h_f/ε.
- For Δt = 0.1, E_s is 6.1–7.3·10⁻⁵, 40–65 % above the 4.4·10⁻⁵ quoted for "none of them depends on Δt" (README table; reproduced: 6.10e-5 and 7.31e-5).
- A remedy is suggested by the authors' own data and should be mentioned: g vanishes near θ ≈ 0.2 and 0.8, so the reference point can be chosen there, or the phase can be dithered at remeshing so that g averages out.
- Literature (from memory; please check details): the dependence of the speed of a discretised travelling wave on its position relative to the grid, and propagation failure, are classical (Keener, SIAM J. Appl. Math. 1987; Elmer and Van Vleck on travelling waves of spatial discretisations of bistable equations; Hupkes–Sandstede and co-authors on lattice travelling waves). The effect of discretisation on frozen waves is analysed by Thümmler (freezing and discretisation of relative equilibria / travelling waves, 2008). These should be cited where l. 1769–1772 states the phenomenon without reference.

---

## Minor issues

1. **Proposition `prop:speed`.** Correct. Table `tab:speedcoef` shows a remainder ≈ 0.1ε² with no visible logarithm; say whether the log is believed to be an artefact of Step 3.
2. **Remark `q:rem:rate`(4).** Baconneau–Lunardi is cited "for the problem with constant M"; for non-constant caloric w_0 the existence of classical solutions is not covered. The self-similar example does satisfy all hypotheses (I checked (S′): h is decreasing, the ball is strictly starshaped).
3. **Theorem `q:thm:rate`.** State explicitly that u_0 ∈ C²(Ω̄_0) with |∇u_0| = √(2M(w_0)) on ∂Ω_0 (implied by the regularity of v at t = 0). (W) requires W^ε(·,0) ∈ C²; say that this restricts the data (it holds for w_0 constant).
4. **Abstract and l. 56.** "f is Lipschitz" in the introduction; (F) of Section 3 adds f ∈ C²((−∞,1]). Say in the introduction that the rate needs it.
5. **l. 100–101.** "a uniform mesh needs at least O(ε^(−N)) cells": unchanged from round 1 (in the regime of the theorems it is ε^(−2N)). **l. 785.** "O(1/ε · ε/h_f)": unchanged; write (band width)/h_f.
6. **Band width.** Introduction and Theorem `an:thm:band`: O(ε log(1/ε)); Sections 4 and 7 and the code: fixed 55ε. The logarithmic growth of the node count in Table `tab:A3` comes from the transition, not from the band. Say which is meant.
7. **Section 4** still describes the adaptive BDF solver with PCHIP transfer as the method, and scheme (S)/ALE as versions of it. The theorems and the main experiments concern (S). Reverse the presentation.
8. **l. 86–91.** "no numerical study of the limit ε → 0 ... is available": still a strong negative claim; the `.bib` still records an unverified numerical appendix in one of the cited papers.
9. **l. 964, l. 1194.** "error estimates, uniform in ε": the constants are uniform but the hypotheses (ρ ≤ γ) tie h and Δt to ε; this is not parameter-uniform convergence in the sense of the literature cited at l. 101–103. Use "with constants independent of ε".
10. **Section `an:sec:sharp`.** "its projection ... averages out" (l. 1584–1587) is true only when the front slides through the lattice; for the pinned ALE front nothing averages, which is the point of Fig. `fig:D3`. "ψ* is bounded" but does not decay; state the weighted space in which the solvability condition is used (round-1 minor 11).
11. **Table `tab:simple`.** The caption does not give the datum (A(1−e^{−c_ε(x−2)})^+, not Φ_ε), T = 12, or the window; see N5.
12. **Secant tolerance.** The implicit velocity is solved to 10⁻³h_f (code); the text does not say. The reference point after remeshing is reset to the current front (`xc = xf`); harmless but not described.
13. **Default ℓ_0.** The default max(2, 3/c) leaves a term ≈ 0.01Δt for w_0 = 1 (reproduced: E = 7.4e-3 at Δt = 0.1 against 2.4e-3 with ℓ_0 = 4, ε = 0.001), growing as ε decreases. It is disclosed; a default that works (≈ 6.5/c according to the README) would be better than a disclosed defect.
14. **Notation.** Still overloaded: κ (constant (3√3−2)/5, damping in Section 3, κ_0 = √a, κ_N); γ (profile parameter in Section 3, smallness constant in Section 6); μ (mesh-velocity error, and f(−w_0)/4 in the proof of Theorem `an:thm:main`); a (f(−w_0), exponent 1 or 2, rescaling constant, a_0); θ (θ(Z), θ_*, θ_1, lattice phase); T (time, remainder in (`q:eq:basic`), temperature); W (W^ε, discrete caloric W^n, and the subsolution W^n_j at l. 1505); G (G(φ) in Section 2, G_±(W) in Section 3).
15. **Bibliography.** Resolve the `checknote = {unverified ...}` entries before submission.

---

## Reproduction

Environment: `/home/user/combustion/.venv/bin/python`, run from `numerics/` with `PYTHONPATH=$PWD OMP_NUM_THREADS=1`. Total compute about 3 minutes. Scratch scripts are in `/tmp/claude-0/-home-user-combustion/c5917af1-f076-5e7f-b209-4cc1d3f058ea/scratchpad/referee2/` (`r1.py` … `r4.py`, `inner.py`, `inner2.py`); below `$S` is that directory.

| # | What | Command | Claimed | Obtained | Verdict |
|---|---|---|---|---|---|
| 1 | D1, ALE from Φ_ε, w_0 = 0, ε ∈ {0.01, 0.001}, Δt ∈ {0.1, 0.02, 0.0025} | `python $S/r1.py D1` | E ≈ 4.5e-4, E_s ≈ 4.4e-5, μ_1 ≈ 1.1e-4, independent of Δt, ε | E = 4.68/4.46/4.42e-4 (ε = 0.01), 4.80/4.50/4.43e-4 (0.001); E_s = 6.1/4.3/4.4e-5 and 7.3/4.3/4.4e-5; μ_1 = 1.17/1.13/1.12e-4 | reproduced; E_s at Δt = 0.1 is 40–65 % larger |
| 2 | D1, μ over all steps | same | 3e-4 to 1.6e-3 | 2.7e-4, 5.7e-4, 1.45e-3 (ε = 0.01); 2.8e-4, 5.8e-4, 1.59e-3 (0.001); increasing as Δt ↓ | reproduced; trend not in the paper (N6) |
| 3 | D1, w_0 = 1, ε = 0.001, default ℓ_0 and ℓ_0 = 4 | `python $S/r1.py D1w1` | default: Δt term ≈ 0.01Δt; ℓ_0 = 4: E ≈ 2.5e-3, E_s ≈ 9e-5 | default: E = 7.4e-3 (Δt = 0.1), 2.9e-3 (0.01); ℓ_0 = 4: E = 2.42e-3, 2.43e-3, E_s = 8.4e-5, 8.6e-5 | reproduced |
| 4 | D3, phase points, w_0 = 0, ε = 0.01, Δt = 0.02, h_f = ε/16 | `python $S/r1.py D3` | −1.9e-4 ≤ g ≤ 1.1e-4 | lsq slope: −1.88e-4 (θ = 1/32), +4.6e-5 (9/32), +8.7e-5 (0.37), +1.08e-4 (17/32), +5e-7 (25/32) | reproduced |
| 5 | Front stays at its phase; pinned speed | `python $S/r2.py 0 0.01 0.02 8`, `... 24` | fixed phase; speed c_ε + g(θ) | phase 0.3700 throughout (10 remeshings); v − c_ε: 1.06e-4 at t ≈ 1.25, 7.97e-5 at t ≈ 8, 8.43e-5 at t ≈ 24 | reproduced; slow relaxation (N5) |
| 6 | D4, drag, w_0 = 0 | `python $S/r1.py D4` | α = 0.33 (Δt/ε = 1), 0.84 (10), 0.98 (100) | 0.328 (ε = 0.01, Δt = 0.01), 0.845 (0.01, 0.1), 0.333 (0.001, 0.001), 0.981 (0.001, 0.1) | reproduced |
| 7 | Table `tab:simple`, ε = 0.01, Δt = ε² | `python $S/r1.py tab` | −2.39e-3 | −2.390e-3 (B1 datum, T = 12; end points −2.393e-3); −2.573e-3 (Φ_ε, T = 6); all steps 1.0e-4 | reproduced (round 1: not reproduced) |
| 8 | Eulerian speed vs time window, ε = 0.02, Δt = 4ε² | `python $S/r4.py` | −1.86e-2 (Table `tab:simple`); coefficient −0.23 | −2.35e-2, −2.01e-2, −1.89e-2, −1.84e-2 on [1,3], [3,6], [6,9], [9,12]; identical for both data; [6,12]: −1.864e-2 | table entry reproduced; window-dependent (N5) |
| 9 | Band transfer exact at remeshing | `python $S/r3.py` | Π is the identity on band nodes | lattice offset ≤ 2.8e-11 h_f at all remeshings (3 runs, 3–12 remeshings) | reproduced (round 1: 0.35 h_f) |
| 10 | Table `tab:B1`, ALE, w_0 = 0 | `python $S/r3.py` | +1.66e-5 (ε = 0.01, Δt = 0.02), +3.73e-5 (0.001, 0.005); three estimators agree | cycle/lsq/end: 1.66/1.65/1.62e-5 and 3.73/3.72/3.72e-5 | reproduced |
| 11 | D2, h_f/ε = 1/8, 1/32 (w_0 = 0, ε = 0.01) | `python $S/r3.py` | E_s = 1.9e-4, 1.1e-5; speed +3.3e-4, +2.2e-5 | E_s = 1.85e-4, 1.06e-5; speed +3.34e-4, +2.13e-5 | reproduced |
| 12 | Cost of the implicit velocity | all ALE runs above | 3 trials per step | 3.00 in every run | reproduced |
| 13 | Lemma `q:lem:inner`, γ-uniform bounds, f = (1−s)^+, κ = 0.05γ, W ∈ {−0.5, 0, 1} | `python $S/inner.py`, `python $S/inner2.py` | ratios bounded uniformly in γ | (a): \|Θ_W\|/((1+\|s\|)Θ_s) ≤ 2.7, \|Θ_sW\|/(…) ≤ 2.1, \|Θ_WW\|/((1+\|s\|)²Θ_s) ≤ 16.1 for γ ≤ 1e-3; (b): ≤ 2.8, 3.1, 25.2 for γ down to 1e-9, with Θ_s ≥ 0.83√γ and s_* ≥ −0.9 log(1/γ) | consistent with the lemma (finite differences in W, two step sizes) |
| 14 | References and citations | grep on `main.tex`, `references.bib` | – | all `\ref`, `\eqref`, `\cite` keys resolve; more than 30 `unverified` checknotes | – |

Not rerun (cost or scope): groups A, B4, C (meshes, radial extinction, limit comparison, quenching), D5 (read from `D_tw_error.json` only), the BDF tables. Checked by hand, not numerically: Proposition `prop:speed`, Lemma `q:lem:ext`, Propositions `q:prop:super`/`q:prop:sub`, Theorem `q:thm:rate`, Lemma `an:lem:ale`, Theorems `an:thm:ale` and `an:thm:band`.
