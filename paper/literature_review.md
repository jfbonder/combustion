# Literature review: numerics for the high activation energy limit (Lewis number 1, w_0 problem)

Compiled 2026-11-22. **How this was checked:** WebFetch and direct HTTP access were blocked in this environment (arxiv.org, link.springer.com, crossref, zbMATH and the other bibliographic hosts all failed to resolve). Every entry below was therefore checked against web-search result pages (publisher pages, arXiv listings, MaRDI/zbMATH/Numdam records, and citing papers' reference lists). The URL recorded for each entry is a page that appeared in those results. Details that no result page confirmed are marked **[unverified]**. Before submission, open each DOI once, and in particular check page ranges.

---

## Bottom line

**I found no numerical work on the high activation energy limit of the equidiffusional (Le = 1) system, whether with the O(ε) offset w_0 or for the reduced case u = Y.** Concretely, I found nothing that:
(a) solves u_t = Δu − β_ε(u) (or the Fernández Bonder–Wolanski system/w_0 formulation) for a sequence ε → 0 and measures convergence to the free boundary problem u_t = Δu in {u>0}, |∇u| = √(2M) on ∂{u>0};
(b) computes the Caffarelli–Vázquez limit problem (planar waves, radial/self-similar solutions) numerically; or
(c) cites Fernández Bonder–Wolanski (2000, 2004) and contains numerics.

I searched for numerics on the singular perturbation problem, the parabolic Bernoulli problem with a fixed-gradient condition, the NEF/free-interface models, and works citing Caffarelli–Vázquez 1995 and Fernández Bonder–Wolanski. The PDE papers on this limit (Caffarelli–Vázquez; CLW; LVW; Petrosyan; Weiss; Kim; Petrosyan–Yip; Weiss–Zhang; Furukawa–Giga–Kajiwara 2024–26) are purely analytical. Weiss–Zhang (arXiv 0904.1261) say explicitly that the degenerate singularity they build "would be otherwise hard to capture" numerically and that their construction "suggests a numerical approach". They give no computations. Search-engine coverage of citing works is incomplete, so a forward-citation check in MathSciNet/zbMATH/Google Scholar before submission is still recommended. The claim "to our knowledge, no numerical study of this limit exists" is defensible.

**Closest works:**
1. **Numerics for combustion fronts at large activation energy, without convergence to the free boundary problem.**
   - Denet–Haldenwang (1992) solve the 2-D thermo-diffusive model. Snippets of the paper report that the gap between numerical growth rates and the large-activation-energy theory behaves like O(1/β). This is the nearest example of measuring convergence in the activation energy, but for a linear-stability quantity, not front position.
   - Bayliss–Matkowsky and coauthors (1987–1995) built adaptive pseudo-spectral methods for combustion fronts with large activation energy.
   - Smooke (1982, 1986) used equidistribution-based adaptive grids for premixed flames.
2. **Numerics for the elliptic Bernoulli problem.** This is the stationary analogue, |∇u| = const on the free boundary. Trial methods, shape-Newton, level sets, CutFEM, and, closest in method, the **moving-mesh FEM of Shen–Dai–Huang (CiCP 2024)**. Its pseudo-time evolution is parabolic-like, but it targets the steady problem and has no singular-perturbation parameter.
3. **Numerics in neighbouring flame free-boundary problems.** One search snippet reports that a Caffarelli–Lee–Mellet homogenization paper has an appendix computing 1-D pulsating fronts of the limit free boundary problem in periodic media [unverified which of the two papers; not read]. Free-interface (ignition-temperature/NEF) models by Brauner and coauthors are studied mostly analytically.

---

## 1. Numerics on this problem or close relatives

**Result: nothing directly on (P_ε) → free boundary problem with convergence in ε.** The items below are the nearest.

- **Weiss, G. S., Zhang, G.** *Existence of a degenerate singularity in the high activation energy limit of a reaction-diffusion equation.* arXiv:0904.1261 (2009); journal version [unverified]. https://arxiv.org/abs/0904.1261
  Constructs u_{ε_k} → u with a degenerate free-boundary point for u_t − Δu = −β_ε(u). The abstract and Section 3 state that the singularity is unstable and hard to capture numerically, and that the construction suggests a numerical approach. There are no computations. Directly relevant to our "instability of the limit" discussion.

- **Denet, B., Haldenwang, P.** *Numerical study of thermal-diffusive instability of premixed flames.* Combust. Sci. Technol. 86 (1992) 199–221. DOI 10.1080/00102209208947195. https://scilit.com/publications/0bb9fc95b2f8c20c0911f73f73495e65
  2-D numerical solution of the thermo-diffusive model, comparing growth rates with large-activation-energy linear theory. Search snippets report a discrepancy that scales as O(1/β). This is a precedent for checking convergence in the activation energy numerically, at finite β, for a different quantity.

- **Bayliss, A., Gottlieb, D., Matkowsky, B. J., Minkoff, M.** *An adaptive pseudo-spectral method for reaction diffusion problems.* J. Comput. Phys. 81 (1989) 421–443. DOI 10.1016/0021-9991(89)90215-5. https://ntrs.nasa.gov/api/citations/19880002032/downloads/19880002032.pdf
  Adaptive (mapped) Chebyshev method for thin combustion fronts in 1-D/2-D. A precedent for adaptive resolution of high-activation-energy fronts. It does not study the ε → 0 limit.

- **Bayliss, A., Matkowsky, B. J.** *Fronts, relaxation oscillations, and period doubling in solid fuel combustion.* J. Comput. Phys. 71 (1987) 147–168. https://portal.mardi4nfdi.de/wiki/Publication:1821540
  Adaptive pseudo-spectral computations of pulsating fronts at large activation energy. These are condensed-phase (gasless) fronts, not the Le = 1 gas model.

- **Bayliss, A., Kuske, R., Matkowsky, B. J.** *A two-dimensional adaptive pseudo-spectral method.* J. Comput. Phys. 91 (1990) 174–196. https://portal.mardi4nfdi.de/wiki/Publication:2638742
  2-D extension of the above (Chebyshev × Fourier), applied to pulsating gasless combustion.

- **Caffarelli, L. A., Lee, K.-A., Mellet, A.** *Singular limit and homogenization for flame propagation in periodic excitable media.* Arch. Ration. Mech. Anal. 172 (2004) 153–190. DOI 10.1007/s00205-003-0299-9. https://link.springer.com/article/10.1007/s00205-003-0299-9
  Same singular perturbation (ε, high activation energy) in periodic media. Shows that pulsating fronts converge to travelling waves when ε ≪ |L| ≪ 1.
- **Caffarelli, L. A., Lee, K.-A., Mellet, A.** *Homogenization and flame propagation in periodic excitable media: the asymptotic speed of propagation.* Comm. Pure Appl. Math. 59 (2006) 501–525. DOI 10.1002/cpa.20094. https://zendy.io/title/10.1002/cpa.20094
  One search snippet (from https://www.math.umd.edu/~mellet/publi/flame2.pdf) reports an appendix with numerical computations of 1-D pulsating travelling fronts of the free boundary problem [unverified: not read; which paper has the appendix is unconfirmed]. If confirmed, it is the only numerical computation of the limit free boundary problem that I located. It does not do ε-convergence of (P_ε).

- **Shen, J., Dai, H., Huang, W.** *A moving mesh finite element method for Bernoulli free boundary problems.* Commun. Comput. Phys. 36 (2024) 248–273. DOI 10.4208/cicp.OA-2023-0214. arXiv:2404.04418. https://ojs.global-sci.org/index.php/cicp/article/view/17190
  Moving-mesh (MMPDE) FEM. A pseudo-transient moving-boundary problem whose steady state solves the Bernoulli problem. The closest methodological relative, since it combines a moving mesh with a |∇u| = const free-boundary condition. It is elliptic/steady, front-tracking and has no ε.

- **Sun, ... et al.** *An adaptive moving finite element method for steady low Mach number compressible combustion problems.* Int. J. Numer. Methods Fluids (2020). DOI 10.1002/fld.4818. arXiv:1903.11030 [full author list unverified]. https://arxiv.org/abs/1903.11030
  r-adaptive MMPDE FEM for premixed laminar flames (2-D, steady). Precedent for moving meshes in flame computations.

- **Analytical papers in the same family; checked, no numerics found:**
  Kim (JDE 2003); Petrosyan–Yip (JGA 2008); Furukawa–Giga–Kajiwara (arXiv 2410.15662, 2608.19479, 2611.15089); the parabolic-Bernoulli self-shrinker paper arXiv 2608.22574. See Sections 4–5.

## 2. Elliptic Bernoulli problem: numerics

- **Flucher, M., Rumpf, M.** *Bernoulli's free-boundary problem, qualitative theory and numerical approximation.* J. Reine Angew. Math. 486 (1997) 165–204. zbl 0909.35154. https://ins.uni-bonn.de/publication/b-ernoulli-s-free-boundary-problem-qualitative-theory-and-numerical-approximation-927.bib?pk=927
  Classical reference. Variational theory plus numerical schemes based on local parametrization of the free boundary.
- **Tiihonen, T.** *Shape optimization and trial methods for free boundary problems.* M2AN 31 (1997) 805–825. https://numdam.org/item/M2AN_1997__31_7_805_0/
  Compares formulations (shape optimization vs. fixed-point "trial" methods) and the conditioning of the shape Hessian.
- **Tiihonen, T., Järvinen, J.** *On fixed point (trial) methods for free boundary problems.* In: Free Boundary Problems in Continuum Mechanics, ISNM 106, Birkhäuser, Basel, 1992, pp. 339–350 [editors per citing paper: Antontsev, Hoffmann, Khludnev]. Cited in Tiihonen 1997: https://www.esaim-m2an.org/articles/m2an/pdf/1997/07/m2an1997310708051.pdf
  Original analysis of trial methods.
- **Kärkkäinen, K. T., Tiihonen, T.** *Free surfaces: shape sensitivity analysis and numerical methods.* Int. J. Numer. Methods Engrg. 44 (1999) 1079–1098. https://jyx.jyu.fi/bitstream/123456789/13270/1/9513923959.pdf
  Shape linearization, the "total linearization method". Fast convergence for Bernoulli-type problems.
- **Bouchon, F., Clain, S., Touzani, R.** *Numerical solution of the free boundary Bernoulli problem using a level set formulation.* Comput. Methods Appl. Mech. Engrg. 194 (2005) 3934–3948. DOI 10.1016/j.cma.2004.09.008. https://hal.archives-ouvertes.fr/hal-00086529
  Level-set evolution in pseudo-time. Handles topology changes and multiple solutions.
- **Bouchon, F., Clain, S., Touzani, R.** *A perturbation method for the numerical solution of the Bernoulli problem.* J. Comput. Math. 26 (2008) 23–36. https://ojs.global-sci.org/index.php/JCM/article/view/11865
- **Ito, K., Kunisch, K., Peichl, G. H.** *Variational approach to shape derivatives for a class of Bernoulli problems.* J. Math. Anal. Appl. 314 (2006) 126–149. DOI 10.1016/j.jmaa.2005.03.100. https://portal.mardi4nfdi.de/wiki/Publication:2577457
  Shape gradient plus level-set steepest descent.
- **Eppler, K., Harbrecht, H.** *Efficient treatment of stationary free boundary problems.* Appl. Numer. Math. 56 (2006) 1326–1339. DOI 10.1016/j.apnum.2006.03.017. https://portal.mardi4nfdi.de/wiki/Publication:2507725
  First- and second-order shape optimization with fast boundary elements.
- **Harbrecht, H., Mitrou, G.** *Improved trial methods for a class of generalized Bernoulli problems.* J. Math. Anal. Appl. 420 (2014) 177–194 [DOI unverified]. https://edoc-vmtest.ub.unibas.ch/34137/4/20140713103342_53c24466cd73d.pdf
  New update rules that speed up trial methods.
- **Bugeanu, M., Harbrecht, H.** *A second order convergent trial method for a free boundary problem in three dimensions.* Interfaces Free Bound. (2015) [volume/pages unverified]. https://ems.press/journals/ifb/articles/13635
- **Burman, E., Elfverson, D., Hansbo, P., Larson, M. G., Larsson, K.** *A cut finite element method for the Bernoulli free boundary value problem.* Comput. Methods Appl. Mech. Engrg. 317 (2017) 598–618. DOI 10.1016/j.cma.2016.12.021. arXiv:1609.02836. https://discovery-pp.ucl.ac.uk/id/eprint/1518358
  CutFEM on a fixed mesh with a level-set/shape-optimization update.
- **Fan, Y., Billingham, J., van der Zee, K.** *A shape-Newton method for free-boundary problems subject to the Bernoulli boundary condition.* SIAM J. Sci. Comput. (2024) [volume/pages unverified]. DOI 10.1137/23M1590263. arXiv:2305.14254. https://nottingham-repository.worktribe.com/index.php/output/39715821/a-shape-newton-method-for-free-boundary-problems-subject-to-the-bernoulli-boundary-condition
  Superlinearly convergent shape-Newton.
- **Shen–Dai–Huang (2024)**: see Section 1.
- **Regularization approach (theory only):** Berestycki–Caffarelli–Nirenberg (1990) and Caffarelli (1995) give uniform estimates for the elliptic singular perturbation Δu = β_ε(u) (Section 4). I found **no** paper that uses the β_ε regularization as a *numerical* method for Bernoulli with a convergence study in ε. Our parabolic study would also cover this gap.

## 3. Moving-mesh and adaptive methods for reaction fronts

- **Huang, W., Russell, R. D.** *Adaptive Moving Mesh Methods.* Applied Mathematical Sciences 174, Springer, New York, 2011. DOI 10.1007/978-1-4419-7916-2. https://www.springer.com/gp/book/9781441979155
  Standard monograph covering equidistribution, monitor functions and MMPDEs. Our 1-D MMPDE/equidistribution choices should cite it.
- **Huang, W., Ren, Y., Russell, R. D.** *Moving mesh partial differential equations (MMPDEs) based upon the equidistribution principle.* SIAM J. Numer. Anal. 31 (1994) 709–730. DOI 10.1137/0731038. https://kuscholarworks.ku.edu/entities/publication/37b63a2a-55ba-4003-a477-b855eae2aa72
  The MMPDE family (MMPDE4/5/6) in 1-D.
- **Huang, W., Russell, R. D.** *Analysis of moving mesh partial differential equations with spatial smoothing.* SIAM J. Numer. Anal. 34 (1997) 1106–1126. DOI 10.1137/S0036142993256441. https://whuang.ku.edu/research/paper/HR97.pdf
  Spatial smoothing of the monitor. No node crossing.
- **Budd, C. J., Huang, W., Russell, R. D.** *Adaptivity with moving grids.* Acta Numerica 18 (2009) 111–241. DOI 10.1017/S0962492906400015. https://www.cambridge.org/core/journals/acta-numerica/article/adaptivity-with-moving-grids/3F7AA4F21C502125FBC832CF9EF84820
  Survey of r-adaptivity. The author-hosted PDF is paginated 1–131. The page range 111–241 comes from citing papers [verify].
- **de Boor, C.** *Good approximation by splines with variable knots.* In: Meir, A., Sharma, A. (eds.), Spline Functions and Approximation Theory, Birkhäuser, 1973, pp. 57–72 [pages unverified]. https://rd.springer.com/chapter/10.1007/978-3-0348-5979-0_3 ; and **de Boor, C.** *Good approximation by splines with variable knots. II.* In: Numerical Solution of Differential Equations (Dundee 1973), Lecture Notes in Math. 363, Springer, 1974, pp. 12–20. DOI 10.1007/BFb0069121.
  The usual source for equidistribution. That these papers contain the equidistribution principle is the standard attribution in the literature; the search snippets did not show it.
- **Stockie, J. M., MacKenzie, J. A., Russell, R. D.** *A moving mesh method for one-dimensional hyperbolic conservation laws.* SIAM J. Sci. Comput. 22 (2001) 1791–1813. DOI 10.1137/S1064827599364428. https://epubs.siam.org/doi/10.1137/S1064827599364428
  Conservative finite-volume solution update on a moving mesh, using a decoupled (predictor–corrector) mesh/PDE solve. Methodologically close to our FV scheme.
- **Smooke, M. D.** *Solution of burner-stabilized premixed laminar flames by boundary value methods.* J. Comput. Phys. 48 (1982) 72–105. DOI 10.1016/0021-9991(82)90036-5. https://portal.mardi4nfdi.de/wiki/Publication:1168062
- **Smooke, M. D.** *Use of adaptive methods in premixed combustion.* AIChE J. 32(8) (1986) [pages unverified]. DOI 10.1002/aic.690320802. https://scilit.com/publications/4830f40a12864647956dfa98334f8813
  Grid selection by equidistributing positive weight functions. The classical adaptive-gridding approach in 1-D flame codes.
- **Day, M. S., Bell, J. B.** *Numerical simulation of laminar reacting flows with complex chemistry.* Combust. Theory Model. 4 (2000) 535–556. DOI 10.1088/1364-7830/4/4/309. https://ccse.lbl.gov/Publications/day/db99/db99_abs.html
  Block-structured AMR low-Mach method for premixed and non-premixed flames.
- **Esclapez, L., et al.** *PeleLMeX: an AMR low Mach number reactive flow simulation code without level sub-cycling.* J. Open Source Softw. 8(90) (2023) 5450. DOI 10.21105/joss.05450. https://joss.theoj.org/papers/10.21105/joss.05450
  The AMReX-based Pele suite, the current exascale AMR combustion code. For context on h-adaptive (AMR) versus r-adaptive (moving mesh) approaches.
- **Miller, J. J. H., O'Riordan, E., Shishkin, G. I.** *Fitted Numerical Methods for Singular Perturbation Problems*, revised ed., World Scientific, Singapore, 2012. ISBN 978-981-4390-73-6. https://worldscientific.com/worldscibooks/10.1142/8410
  Shishkin (piecewise-uniform, layer-adapted) meshes. Parameter-uniform convergence.
- **Kopteva, N., O'Riordan, E.** *Shishkin meshes in the numerical solution of singularly perturbed differential equations.* Int. J. Numer. Anal. Model. 7 (2010) 393–415. https://ojs.global-sci.org/index.php/ijnam/article/view/9961
  Review. Relevant contrast: Shishkin meshes need a priori layer locations, but our interior layer moves, which motivates r-adaptivity.
- **Clavero, C., Gracia, J. L.** — Shishkin-mesh schemes for 1-D parabolic singularly perturbed reaction-diffusion problems. A 2010 Appl. Math. Comput. paper is cited in https://link.springer.com/content/pdf/10.1007%2F978-3-642-19665-2_9.pdf [exact title/volume/pages unverified; not included in the .bib].

## 4. Theory relevant to the limit problem (mainly after 2000, plus classical background)

**Classical background**
- **Zeldovich, Ya. B., Frank-Kamenetskii, D. A.** *K teorii ravnomernogo rasprostraneniya plameni* [On the theory of uniform flame propagation]. Doklady Akad. Nauk SSSR 19(9) (1938) 693–697. https://elib.gstu.by/handle/220612/31930
  Origin of the high activation energy (ZFK) approach. A parallel Acta Physicochimica URSS version is [unverified].
- **Buckmaster, J. D., Ludford, G. S. S.** *Theory of Laminar Flames.* Cambridge Univ. Press, 1982. https://www.cambridge.org/core/books/theory-of-laminar-flames/4AB085396CB9DA40C0207266DCAFDDDF
  Activation-energy asymptotics, including reaction-sheet jump conditions such as |∇u| = √(2M).
- **Matkowsky, B. J., Sivashinsky, G. I.** *An asymptotic derivation of two models in flame theory associated with the constant density approximation.* SIAM J. Appl. Math. 37 (1979) 686–699. DOI 10.1137/0137051. https://epubs.siam.org/doi/10.1137/0137051
  Derives the NEF free-interface model, in which β(Le−1) = O(1).
- **Berestycki, H., Nicolaenko, B., Scheurer, B.** *Traveling wave solutions to combustion models and their singular limits.* SIAM J. Math. Anal. 16 (1985) 1207–1242. DOI 10.1137/0516088. https://epubs.siam.org/doi/10.1137/0516088
  Existence of travelling waves and convergence of the speed in the high activation energy limit. This is the theoretical basis for our planar-wave validation.
- **Berestycki, H., Caffarelli, L. A., Nirenberg, L.** *Uniform estimates for regularization of free boundary problems.* In: Analysis and Partial Differential Equations, Lecture Notes in Pure and Appl. Math. 122, Dekker, New York, 1990, pp. 567–619. Cited in https://bibliotecadigital.exactas.uba.ar/download/paper/paper_00217824_v86_n6_p552_Lederman.pdf
- **Caffarelli, L. A.** *Uniform Lipschitz regularity of a singular perturbation problem.* Differential Integral Equations 8 (1995) 1585–1590. https://portal.mardi4nfdi.de/wiki/Publication:1894156
- **Caffarelli, L. A., Vázquez, J. L.** *A free-boundary problem for the heat equation arising in flame propagation.* Trans. Amer. Math. Soc. 347 (1995) 411–441. DOI 10.1090/S0002-9947-1995-1260199-7. https://www.ams.org/tran/1995-347-02/S0002-9947-1995-1260199-7/
  Uniform estimates for u_t = Δu − β_ε(u) and limit solutions of the free boundary problem. The source of our radial self-similar test.
- **Vázquez, J. L.** *The free boundary problem for the heat equation with fixed gradient condition.* In: Free Boundary Problems, Theory and Applications (Zakopane 1995), Pitman Res. Notes Math. 363, Longman, Harlow, 1996, pp. 277–302. https://chenle02.github.io/SPDEs-Bib/bib_entries/vazquez:96:free.html
  Lecture notes. Contains Vázquez's nonuniqueness example, later made rigorous by Petrosyan–Yip.
- **Galaktionov, V. A., Hulshof, J., Vázquez, J. L.** *Extinction and focusing behaviour of spherical and annular flames described by a free boundary problem.* J. Math. Pures Appl. 76 (1997) 563–608. DOI 10.1016/S0021-7824(97)89963-1. https://portal.mardi4nfdi.de/wiki/Publication:1385313
  Radial solutions: asymptotically self-similar near extinction. Relevant to our radial tests.
- **Caffarelli, L. A., Lederman, C., Wolanski, N.** *Uniform estimates and limits for a two phase parabolic singular perturbation problem.* Indiana Univ. Math. J. 46 (1997) 453–489 (some citations give 453–490). Cited in https://numdam.org/item/ASNSP_1998_4_27_2_253_0/ and https://bibliotecadigital.exactas.uba.ar/download/paper/paper_00217824_v86_n6_p552_Lederman.pdf
- **Caffarelli, L. A., Lederman, C., Wolanski, N.** *Pointwise and viscosity solutions for the limit of a two phase parabolic singular perturbation problem.* Indiana Univ. Math. J. 46 (1997) 719–740. (Same sources.)
- **Lederman, C., Wolanski, N.** *Viscosity solutions and regularity of the free boundary for the limit of an elliptic two phase singular perturbation problem.* Ann. Sc. Norm. Super. Pisa Cl. Sci. (4) 27 (1998) 253–288. https://numdam.org/item/ASNSP_1998_4_27_2_253_0/

**After 2000**
- **Lederman, C., Vázquez, J. L., Wolanski, N.** *Uniqueness of solution to a free boundary problem from combustion.* Trans. Amer. Math. Soc. 353 (2001) 655–692. DOI 10.1090/S0002-9947-00-02663-5. https://www.ams.org/journals/tran/2001-353-02/S0002-9947-00-02663-5
  Classical, limit and viscosity solutions agree and are unique under geometric assumptions.
- **Lederman, C., Vázquez, J. L., Wolanski, N.** *Uniqueness in a two-phase free-boundary problem.* Adv. Differential Equations 6 (2001) 1409–1442. https://projecteuclid.org/journals/advances-in-differential-equations/volume-6/issue-12/Uniqueness-in-a-two-phase-free-boundary-problem/ade/1357139953.full
- **Petrosyan, A.** *On existence and uniqueness in a free boundary problem from combustion.* Comm. Partial Differential Equations 27 (2002) 763–789. MR1900562. Cited in https://www.math.purdue.edu/~yip/papers/nonuniqueFLAME.pdf [DOI unverified].
- **Petrosyan, A.** *Convexity and uniqueness in a free boundary problem arising in combustion theory.* Rev. Mat. Iberoam. [vol./year/pages unverified; EMS lists it]. https://ems.press/journals/rmi/articles/5060
- **Weiss, G. S.** *Self-similar blow-up and Hausdorff dimension estimates for a class of parabolic free boundary problems.* SIAM J. Math. Anal. 30 (1999) [vol./pages unverified]. DOI 10.1137/S0036141097327409. https://epubs.siam.org/doi/10.1137/S0036141097327409
- **Weiss, G. S.** *A singular limit arising in combustion theory: fine properties of the free boundary.* Calc. Var. Partial Differential Equations 17 (2003) 311–340. DOI 10.1007/s00526-002-0171-z. https://portal.mardi4nfdi.de/wiki/Item:Q1404476
  Cauchy problem for u_t − Δu = −β_ε(u). The limit free boundary is rectifiable up to a small set. The limit is a solution in a weak "domain variation" sense.
- **Kim, I. C.** *A free boundary problem arising in flame propagation.* J. Differential Equations 191 (2003) 470–489 [volume unverified]. DOI 10.1016/S0022-0396(02)00195-X. https://www.sciencedirect.com/science/article/pii/S002203960200195X
  Viscosity solutions, uniform convergence of u_ε and uniqueness under geometric conditions on the data.
- **Danielli, D., Petrosyan, A., Shahgholian, H.** *A singular perturbation problem for the p-Laplace operator.* Indiana Univ. Math. J. 52 (2003) 457–476. https://asu.elsevierpure.com/en/publications/a-singular-perturbation-problem-for-the-p-laplace-operator
  Elliptic, uniform Lipschitz estimates.
- **Lederman, C., Wolanski, N.** *Singular perturbation in a nonlocal diffusion problem.* Comm. Partial Differential Equations 31 (2006) 195– [end page unverified]. DOI 10.1080/03605300500358111. https://www.tandfonline.com/doi/full/10.1080/03605300500358111
  Treats both a single equation and a two-equation system (fuel/temperature) with nonlocal diffusion, a direct generalization of the Fernández Bonder–Wolanski setting.
- **Lederman, C., Wolanski, N.** *A two phase elliptic singular perturbation problem with a forcing term.* J. Math. Pures Appl. 86 (2006) 552–589. DOI 10.1016/j.matpur.2006.10.008. https://portal.mardi4nfdi.de/wiki/Publication:864197
- **Lederman, C., Oelz, D.** *A quasilinear parabolic singular perturbation problem.* Interfaces Free Bound. 10 (2008) 447–482. DOI 10.4171/IFB/197. https://ems.press/journals/ifb/articles/1893
- **Petrosyan, A., Yip, N. K.** *Nonuniqueness in a free boundary problem from combustion.* J. Geom. Anal. 18 (2008) 1098–1126. DOI 10.1007/s12220-008-9044-9. https://link.springer.com/article/10.1007/s12220-008-9044-9
  Rigorous two-hump example with at least two limit solutions in n ≥ 2, from minimal and maximal approximations. Strongly relevant to the "instability of the limit" theme. Uniqueness needs geometric conditions such as starshapedness.
- **Andersson, J., Weiss, G. S.** *A parabolic free boundary problem with Bernoulli type condition on the free boundary.* J. Reine Angew. Math. 627 (2009) 213–235. DOI 10.1515/CRELLE.2009.016. https://arxiv.org/abs/math/0608746
  One-sided flatness implies regularity.
- **Monneau, R., Weiss, G. S.** *Pulsating traveling waves in the singular limit of a reaction-diffusion system in solid combustion.* Ann. Inst. H. Poincaré Anal. Non Linéaire 26 (2009) 1207–1222. DOI 10.1016/j.anihpc.2008.09.002. https://www.numdam.org/articles/10.1016/j.anihpc.2008.09.002
- **Weiss–Zhang (2009)**: see Section 1.
- **Kamburov, N.** *Nondegeneracy and stability in the limit of a one-phase singular perturbation problem.* arXiv:2302.13422 (2023). https://arxiv.org/abs/2302.13422v1 [journal version possibly Math. Ann., DOI 10.1007/s00208-023-02591-0, unverified].
- **Furukawa, K., Giga, Y., Kajiwara, N.** *No formation of a new phase for a free boundary problem in combustion theory.* arXiv:2410.15662; journal version DOI 10.1007/s41808-025-00333-8 [journal name unverified]. https://arxiv.org/html/2410.15662v2
  Parabolic Bernoulli problem with a source, 1-D, Laplace transform.
- **Furukawa, K., Giga, Y., Kajiwara, N.** *Solutions to a one-dimensional combustion-type free boundary problem via maximal regularity.* arXiv:2608.19479 (2026). https://arxiv.org/pdf/2608.19479
  1-D local well-posedness of the limit problem. Relevant to our 1-D/radial setting.
- *Well-posedness for the Cauchy problem with free boundary in half-space arising from combustion problem.* arXiv:2611.15089 (2026) [authors attributed to Furukawa–Giga–Kajiwara by a citing reference; unverified]. https://arxiv.org/pdf/2611.15089
- *Singularity formation and self-shrinkers for the parabolic Bernoulli problem.* arXiv:2608.22574 (2026) [authors unverified]. https://arxiv.org/pdf/2608.22574
  Self-similar shrinking solutions of the limit problem.
- *Regularity of the free boundary in an unstable parabolic problem.* arXiv:2508.15997 (2025) [authors unverified]. https://arxiv.org/html/2508.15997v1

**NEF and free-interface models**
- **Brauner, C.-M., Lunardi, A.** *Instabilities in a two-dimensional combustion model with free boundary.* Arch. Ration. Mech. Anal. 154 (2000) 157–182. DOI 10.1007/s002050000099. https://link.springer.com/doi/10.1007/s002050000099
  Instability of the planar NEF travelling wave.
- **Brauner, C.-M., Hulshof, J., Lorenzi, L.** *Stability of the travelling wave in a 2D weakly nonlinear Stefan problem.* Kinet. Relat. Models 2 (2009) 109–134. DOI 10.3934/krm.2009.2.109. https://portal.mardi4nfdi.de/wiki/Publication:1042712
- **Ducrot, A., Marion, M.** *Two-dimensional travelling wave solutions of a system modelling near equidiffusional flames.* Nonlinear Anal. (2005) 1105–1134 [volume unverified]. https://oskar-bordeaux.fr/handle/20.500.12278/190935
  Travelling waves in a strip for high activation energy (NEF scaling).
- **Lederman, C., Roquejoffre, J.-M., Wolanski, N.** *Mathematical justification of a nonlinear integro-differential equation for the propagation of spherical flames.* C. R. Math. Acad. Sci. Paris 334 (2002) 569–574. DOI 10.1016/S1631-073X(02)02299-9. https://www.numdam.org/item/CRMATH_2002__334_7_569_0 ; full version in Ann. Mat. Pura Appl. 183 (2004) 173– [end page unverified]. DOI 10.1007/s10231-003-0085-1. https://link.springer.com/article/10.1007/s10231-003-0085-1
  Lewis number close to 1. Cites Fernández Bonder–Wolanski 2000.
- **Brauner, C.-M., Gordon, P., Zhang, W.** *An ignition-temperature model with two free interfaces in premixed flames.* Combust. Theory Model. 20 (2016) 976–994. DOI 10.1080/13647830.2016.1220625. https://portal.mardi4nfdi.de/wiki/Item:Q5073667
- **Langlois thesis:** not found by any search. Omitted.

## 5. Works citing Fernández Bonder–Wolanski (2000, 2004)

- **The two papers themselves:**
  - Fernández Bonder, J., Wolanski, N. *A free-boundary problem in combustion theory.* Interfaces Free Bound. 2 (2000) 381–411. DOI 10.4171/IFB/26. https://ems.press/journals/ifb/articles/40 . Note: the EMS page appears to show a wrong abstract (one about the modified Stefan problem). The UBA record https://bibliotecadigital.exactas.uba.ar/collection/paper/document/paper_14639963_v2_n4_p381_Bonder has the correct abstract.
  - Fernández Bonder, J., Wolanski, N. *Uniqueness of limit solutions to a free boundary problem from combustion.* SIAM J. Math. Anal. 36 (2004) 172–185 [end page as given by author]. DOI 10.1137/S0036141002412938. https://epubs.siam.org/doi/10.1137/S0036141002412938
  - Fernández Bonder, J. *Un problema de frontera libre en teoría de combustión.* PhD thesis, Universidad de Buenos Aires (advisor N. Wolanski), 2002 [year per repository record]. https://repositoriosdigitales.mincyt.gob.ar/vufind/Record/BDUBAFCEN_c450590507fa30794538c4aab2e917cb
- **Confirmed citing work:** Lederman–Roquejoffre–Wolanski, Ann. Mat. Pura Appl. 2004. Its reference list includes FB–W 2000.
- **Likely citing works, not confirmed from reference lists:** Lederman–Wolanski CPDE 2006 (nonlocal two-equation system), Lederman–Oelz IFB 2008, Petrosyan–Yip JGA 2008 (on uniqueness of limit solutions). Search results placed FB–W 2004 next to these, but their reference lists were not visible.
- **No citing paper with numerics was found.** A full forward-citation list needs MathSciNet/zbMATH/Google Scholar, which could not be reached from this environment.
