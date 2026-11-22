# Papers para conseguir

No se pudieron bajar desde el entorno, porque los sitios están bloqueados. Están ordenados por prioridad. Subirlos a `paper/pdfs/`; esa carpeta conviene no versionarla si son copias con derechos.

## Prioridad 1: ejemplos para comparar con los experimentos

1. **Caffarelli, L. A.; Vázquez, J. L.** *A free-boundary problem for the heat equation arising in flame propagation.* Trans. Amer. Math. Soc. 347 (1995), 411–441. doi:10.1090/S0002-9947-1995-1260199-7. Es de acceso libre en ams.org.
   - Origen del problema límite.
   - Ondas viajeras y soluciones autosemejantes.
   - Ejemplos.
2. **Galaktionov, V. A.; Hulshof, J.; Vázquez, J. L.** *Extinction and focusing behaviour of spherical and annular flames described by a free boundary problem.* J. Math. Pures Appl. 76 (1997), 563–608. doi:10.1016/S0021-7824(97)89963-1.
   - Llamas esféricas y anulares: extinción y focalización.
   - Es exactamente la geometría radial de nuestros experimentos. Probablemente tenga perfiles y tasas para comparar.
3. **Alt, H. W.; Caffarelli, L. A.** *Existence and regularity for a minimum problem with free boundary.* J. Reine Angew. Math. 325 (1981), 105–144. doi:10.1515/crll.1981.325.105. Hay un escaneo en eudml.org/doc/152360.
   - Problema de Bernoulli elíptico, con ejemplos. Sirve como banco de prueba para el elíptico 2D.
4. **Berestycki, H.; Caffarelli, L. A.; Nirenberg, L.** *Uniform estimates for regularization of free boundary problems.* En: Analysis and Partial Differential Equations, Lecture Notes in Pure and Appl. Math. 122, Dekker, New York, 1990, pp. 567–619.
   - Regularización de β_ε en el caso elíptico.
5. **Vázquez, J. L.** *The free boundary problem for the heat equation with fixed gradient condition.* En: Free Boundary Problems, Theory and Applications (Zakopane, 1995), Pitman Res. Notes Math. 363, Longman, 1996, pp. 277–302.
   - Panorama del problema límite, con ejemplos.

## Prioridad 2: casos difíciles y posibles desafíos numéricos

6. **Petrosyan, A.; Yip, N. K.** *Nonuniqueness in a free boundary problem from combustion.* J. Geom. Anal. 18 (2008), 1098–1126. doi:10.1007/s12220-008-9044-9.
   - Ejemplos de no unicidad. Hay que ver si son reproducibles numéricamente.
7. **Weiss, G. S.; Zhang, G.** *Existence of a degenerate singularity in the high activation energy limit of a reaction-diffusion equation.* arXiv:0904.1261 (2009).
   - Singularidad degenerada. Los autores dicen que es "difícil de capturar" numéricamente.
8. **Weiss, G. S.** *A singular limit arising in combustion theory: fine properties of the free boundary.* Calc. Var. PDE 17 (2003), 311–340. doi:10.1007/s00526-002-0171-z.
9. **Caffarelli, L. A.; Lederman, C.; Wolanski, N.** *Uniform estimates and limits for a two phase parabolic singular perturbation problem.* Indiana Univ. Math. J. 46 (1997), 453–490.
   - Choque de frentes: u ≈ α x₁⁺ + ᾱ x₁⁻. Es el análogo de nuestro Lema 2.3 de 2000.
10. **Caffarelli, L. A.; Lederman, C.; Wolanski, N.** *Pointwise and viscosity solutions for the limit of a two phase parabolic singular perturbation problem.* Indiana Univ. Math. J. 46 (1997), 719–740.
11. **Lederman, C.; Vázquez, J. L.; Wolanski, N.** *Uniqueness of solution to a free boundary problem from combustion.* Trans. Amer. Math. Soc. 353 (2001), 655–692. doi:10.1090/S0002-9947-00-02663-5. Es de acceso libre en ams.org.
12. **Petrosyan, A.** *On existence and uniqueness in a free boundary problem from combustion.* Comm. PDE 27 (2002), 763–789.

## Prioridad 3: para el artículo (numérica y ondas viajeras)

13. **Berestycki, H.; Nicolaenko, B.; Scheurer, B.** *Traveling wave solutions to combustion models and their singular limits.* SIAM J. Math. Anal. 16 (1985), 1207–1242. doi:10.1137/0516088.
    - Convergencia de la velocidad de las ondas viajeras cuando ε → 0. Hay que comparar con nuestra Proposición sobre la corrección O(ε).
14. **Shen, J.; Dai, H.; Huang, W.** *A moving mesh finite element method for Bernoulli free boundary problems.* Commun. Comput. Phys. 36 (2024), 248–273. doi:10.4208/cicp.OA-2023-0214. También en arXiv:2404.04418.
    - El trabajo numérico más cercano al nuestro, aunque es elíptico.
15. **Caffarelli, L. A.; Lee, K.-A.; Mellet, A.** *Homogenization and flame propagation in periodic excitable media: the asymptotic speed of propagation.* Comm. Pure Appl. Math. 59 (2006), 501–525. doi:10.1002/cpa.20094.
    - Según un resultado de búsqueda, tendría un apéndice numérico sobre frentes 1D del problema límite. Hay que verificarlo.
16. **Caffarelli, L. A.** *Uniform Lipschitz regularity of a singular perturbation problem.* Differential Integral Equations 8 (1995), 1585–1590.
17. **Huang, W.; Russell, R. D.** *Adaptive Moving Mesh Methods.* Springer, 2011. doi:10.1007/978-1-4419-7916-2.
