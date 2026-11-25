# Referato interno: decisión editorial y recomendaciones

**Manuscrito:** *The high activation energy limit of an equidiffusional flame model: rate of convergence and a front-following numerical method* (J. Fernández Bonder), `paper/main.tex`, borrador de noviembre de 2026.
**Revistas consideradas:** SIAM J. Numer. Anal. (SINUM), Math. Comp., Numer. Math.
**Fecha:** 25 de noviembre de 2026.

El referato es interno y adversarial. Lo hicieron cuatro referís independientes:

| Referí | Foco | Recomendación |
|---|---|---|
| R1, físico de combustión | Modelo, casos estudiados, interpretación física de los resultados | Revisión mayor |
| R2, experto en EDPs | Estimaciones cuantitativas (Secs. 2–3), pruebas de la Sec. 5, nivel de detalle | Revisión mayor |
| R3, analista numérico | Método, originalidad, pruebas de convergencia (Secs. 4–6), adecuación a cada revista | Revisión mayor (SINUM); rechazo con invitación a reenviar (Math. Comp., Numer. Math.) |
| R4, ingeniero de software científico | Código `numerics/`, correspondencia con el artículo, rendimiento, reproducibilidad | Revisiones menores a moderadas |

Los cuatro informes completos están en el Anexo.

**Convenciones de este documento**
- Los números de línea remiten a `paper/main.tex` en el commit `52c3129`.
- Los resultados se numeran según el PDF compilado: Lema 2.1 y Prop. 2.2 (Sec. 2); Lema 3.1, Def. 3.2, Lema 3.3, Props. 3.4–3.5, Teorema 3.6 (tasa) y Obs. 3.7 (Sec. 3); Lemas 4.1–4.2 (Sec. 4); Teorema 5.1 (malla fija), Teorema 5.2 (malla móvil), Obs. 5.3, Teorema 5.4 (malla de banda) y Obs. 5.5 (Sec. 5).
- En el Anexo algunos referís numeran distinto. R1 llama "Teorema 3.7" al de tasa. R2 llama "Teoremas 4.3–4.5" a los de la Sec. 5. R3 llama "Thm 5.3" al de banda. En todos los casos se refieren a los resultados de la lista anterior.

---

## 1. Decisión editorial

**Revisión mayor.** Tal como está, no recomiendo enviarlo a ninguna de las tres revistas.

Los cuatro referís coinciden en que **la matemática es correcta en lo esencial**:
- R2 revisó línea por línea las Secciones 2, 3 y 5 y no encontró errores fatales.
- R3 verificó las pruebas de los Lemas 4.1–4.2 y de los Teoremas 5.1, 5.2 y 5.4, también contra el código.
- R4 confirmó que el código implementa fórmula por fórmula el esquema descripto y reprodujo parte de los resultados publicados.

Los problemas no están en la corrección sino en tres puntos:

1. **El paper son dos artículos pegados.** La Sección 3 (tasa del límite singular) es EDP pura y **no alimenta** el análisis numérico; sólo se usa para interpretar la Tabla 3. La Sección 5 (análisis numérico) cubre sólo la onda viajera exacta en 1D. R2 y R3, por separado, sugieren dividir el trabajo.
2. **El resultado numérico es estrecho para SINUM, Math. Comp. o Numer. Math.** Las cotas de error no dicen nada cuantitativo sobre ningún experimento del paper (ver M5). Además, no hay teorema alguno sobre el esquema que realmente se implementa.
3. **Varias afirmaciones exceden lo probado:**
   - la optimalidad del orden ε;
   - el alcance del Teorema 3.6;
   - "Δt no aparece en el error";
   - el "quenching" y la "extinción";
   - la ventaja de costo frente a BDF.

### Recomendación estratégica del editor

Hay dos caminos. Recomiendo intentar primero el **Camino A**, porque resuelve de un golpe las objeciones 1 y 2.

**Camino A: un solo artículo para SINUM, unificado por un teorema nuevo.**
Probar la convergencia del esquema totalmente discreto al problema límite (P), para datos generales (estrellados, como en el Teorema 3.6) y con tasa combinada, del tipo

  ‖U^n_h − v(t_n)‖ ≤ C(ε log(1/ε) + η_ε + h_f/ε + …).

R3 señala que la prueba saldría de combinar el principio de comparación discreto (Lema 4.2) con las barreras explícitas de la Sección 3, al estilo Nochetto–Verdi. Así las barreras de la Sección 3 pasan a ser la herramienta del análisis numérico, y la unidad del artículo queda justificada. Con este teorema, SINUM es un destino realista. Math. Comp. y Numer. Math. requerirían además que la Sec. 5.5 pase a ser un teorema riguroso (ver M8).

**Camino B: dividir el trabajo** si el teorema del Camino A no sale en un tiempo razonable.
- **Artículo EDP** (Secs. 2–3 más mejoras): para SIMA, IFB, JDE o AIHP-C. Las mejoras son eliminar el log(1/ε) si se puede (M2), delimitar la novedad frente a Petrosyan 2002 y FBW 2004 (M3), reformular la optimalidad (M1) y precisar la clase de soluciones (M4).
- **Artículo numérico** (Sec. 2.1 y Secs. 4–6 más mejoras): para SINUM, o para JCP o SISC si el énfasis queda en el método. Las mejoras son un teorema sobre el esquema implementado (M6), la reformulación de "Δt no aparece" (M7), la comparación trabajo-precisión (M11) y la discusión de originalidad (M10).

No recomiendo Math. Comp. ni Numer. Math. en ninguno de los dos caminos, salvo que el resultado numérico se amplíe sustancialmente.

---

## 2. Fortalezas (consenso de los referís)

1. **La idea central de la Sección 3 es sólida y nueva** (R2). Lo nuevo es obtener el umbral ε₀(δ) *lineal* en la "strictness" δ mediante una barrera explícita, en lugar del blow-up de FBW 2004. Las herramientas son los perfiles internos parametrizados sólo por W, con energía G±² = 2M(W) ± 2γ y cotas uniformes en γ, y la extensión monótona ṽ± = Ê ± Λ'(Ê⁻)². R2 verificó las constantes; por ejemplo, 8Λ'δ₀³ − 2K₀² = Λ + δ₀K₀. No encontró ningún trabajo que contradiga que sea la primera tasa para el problema de evolución.
2. **El Lema 2.1 y la Prop. 2.2 son limpios** (R1, R2).
   - κ* = (3√3−2)/5 y la columna teórica de la Tabla 1 están verificados.
   - La constante |∇u| = √(2M) coincide con la fórmula ZFK; con f = ½e^{−s} reproduce la normalización de Buckmaster–Ludford.
   - La Tabla 1 coincide con `speed_coefficient.json`, que se regenera bit a bit (R4).
3. **El principio de comparación discreto (Lema 4.2) vale para todo Δt y toda malla.** Es la pieza analítica más valiosa de la parte numérica (R2, R3): reduce el análisis de error a construir barreras.
4. **El método es robusto.** R3 corrió el ALE desde la onda exacta hasta Δt/ε ≈ 333: el error queda por debajo de 10⁻³ con la velocidad implícita.
5. **La relación entre el código y el artículo es buena** (R4).
   - El Laplaciano, el operador D, los signos ALE, la condición de Péclet, la secante, el remallado y el solver del problema límite coinciden con el texto, sin errores de signo ni de upwinding.
   - Las cotas del Lema 4.1 se cumplen numéricamente.
   - `pytest` pasa (19/19).
   - `B_convergence.json` se reproduce con una diferencia máxima de 1.2·10⁻⁸.
   - `reproduce_paper.sh` cubre todas las figuras y tablas.
6. **Las explicaciones de fenómenos numéricos finos son correctas** (R1, R3): el residuo 0.3ε de la Tabla 4 a partir de la Prop. 2.2, el retraso Δt/ε del esquema euleriano por la condición de solubilidad, y la relajación lenta en escala 4/c².
7. **El texto es honesto sobre lo que no prueba** (Sec. 5.4, el log, el alcance limitado a ondas viajeras).

---

## 3. Problemas mayores, por prioridad

Cada ítem indica quién lo marcó, dónde está, qué falla y qué cambio se propone.

### Parte EDP (Secciones 2–3)

**M1. La optimalidad del orden ε está sobredimensionada.** (R2; verificado por el editor)
- **Dónde:** l. 109 ("and the order ε cannot be improved"), abstract, Obs. 3.7 (l. 1005–1014), Conclusiones.
- **Qué falla:**
  - La introducción lo afirma sin condiciones. La Obs. 3.7 sólo lo justifica en la región quemada, donde u^ε ≈ −εW^ε, y sólo si w₀ ≢ 0.
  - No dice nada sobre el caso w₀ ≡ 0 (u = Y, el de Caffarelli–Vázquez), ni sobre (u^ε)⁺ − v (lo que mide la Tabla 3), ni sobre la posición del frente.
  - El ejemplo de la onda viajera queda fuera del teorema, y su desplazamiento ε/√(2M) depende del nivel que se elija como "frente".
- **Cambio:** (a) reformular como "óptimo en L^∞ para u^ε − v cuando w₀ ≢ 0" en el abstract, la introducción y las conclusiones; o (b) probar una cota inferior genuina para el frente, por ejemplo usando c_ε − c₀ = βε con β ≠ 0 de la Prop. 2.2.

**M2. El log(1/ε) de la cota inferior del Teorema 3.6 probablemente es un artefacto.** (R2)
- **Dónde:** Teorema 3.6, l. 625, Obs. 3.7(1) (l. 995–1004).
- **Qué falla:** el log viene de comparar en t = 0 con Y₀ = 0 en la región quemada. El perfil Θ⁻ necesita una longitud log(1/γ) para llegar a 0, con γ ≲ δ ~ ε. El propio texto admite (l. 1002–1004) que los experimentos no muestran el log.
- **Ruta propuesta** (R2 no la verificó en todos los detalles):
  1. Usar que Y_t ≥ ΔY − Cε⁻²Y.
  2. Tomar una subsolución de tipo "frente exponencial" εe^{−β(d−Vt)/ε}, con V ~ 1/ε y β² > C. Genera una cola Y ≥ cεe^{−βd/ε} en un tiempo t₁ = Cε² log(1/ε).
  3. Desde t₁, aplicar la Obs. 3.7(2) a v(·, · + t₁), a costo O(t₁). Requiere que v(·, t₁) siga estrictamente estrellada, lo que vale por continuidad.
- **Resultado esperado:** −C(ε + η) ≤ u^ε − v ≤ C(ε + η), simétrico.
- **Prioridad:** alta. Un resultado simétrico y limpio es mucho más vendible. Si no sale, dejar el log y decir con precisión de dónde viene.

**M3. La novedad frente a los trabajos previos no está delimitada.** (R2)
- **Qué falla:**
  - FBW 2004 (`final.tex`, l. 242–246) atribuye a **Petrosyan 2002** una construcción explícita de supersoluciones aproximantes para w = 0. El manuscrito no la discute: hay que explicarla y decir si ya da ε₀(δ) ≥ cδ.
  - La cercanía O(ε) ya estaba en el Step II de la prueba del Thm 4.1 de FBW 2004 (`final.tex`, ~l. 690). Lo que no estaba cuantificado era ε₀(δ), escondido en el blow-up del Step III.
  - Falta explicitar la diferencia con Chen 1992, Alfaro–Hilhorst–Matano 2008 y de Mottoni–Schatzman:
    - el perfil está acoplado al gradiente exterior;
    - la energía es G(W) con margen cδ, no la del gradiente local;
    - el reescalado sustituye al corrimiento temporal ±εp(t);
    - la barrera es "gruesa" y no captura el término de orden ε.
- **Cambio:**
  - Agregar un párrafo de "relación con trabajos previos" en la Sección 3.
  - Revisar Glangetas 1992 (Asymptotic Anal. 5), Glangetas–Roquejoffre 1996, Langlois–Marion y Brauner–Lunardi.
  - Formular la contribución como "una contribución real pero acotada", sin inflarla.

**M4. La clase de soluciones a la que se aplica el Teorema 3.6 es más chica de lo que sugiere el texto.** (R1, R2; verificado por el editor)
- **Qué falla:**
  - El teorema exige u₀ con soporte compacto (l. 881) y datos estrictamente estrellados: la región fresca {u > 0} es una bolsa acotada rodeada de gas quemado. Quedan fuera la llama plana, las llamas en expansión y la ignición.
  - Por eso "in agreement with Theorem 3.6" (l. 1821) excede el alcance en el caso plano.
  - La existencia de soluciones clásicas regulares sólo está asegurada por Baconneau–Lunardi 2004 para M constante, con datos lisos y compatibles y en tiempos cortos. Para w₀ variable no hay referencia.
  - La regularidad uniforme hasta t = 0 que se exige ((R1)) es una hipótesis de compatibilidad.
- **Cambio:**
  - Enunciar explícitamente la clase no vacía conocida: w₀ constante, datos C^{2+α} compatibles y estrellados, soluciones autosemejantes.
  - Aclarar que (S) sólo se usa en t = 0.
  - Corregir la l. 1821 y cualquier frase que sugiera aplicabilidad a frentes planos.

**M5. Las hipótesis están mal atribuidas.** (R2)
- l. 114–115: "for the rate we need f nonincreasing and of class C² up to 1". Con w₀ constante, el término 𝒫 se anula y basta (F); (F₂) sólo hace falta para ∂_W y ∂²_W. La frase es falsa en el caso principal.
- La monotonía de f se usa en el Lema 2.1, en el Lema 3.1(a) y en el Lema 4.2, pero **no** en la Prop. 3.4 (l. 763), aunque el texto diga que sí. En la parte EDP podría reemplazarse por una condición de no degeneración en 1.
- l. 178: "it holds … for Arrhenius-type rates" se presta a confusión (R1). e^{−s} es monótona pero no cumple el resto de (F), porque no se anula en [1, ∞). Decir "Arrhenius truncada" o explicar el corte.

### Parte numérica (Secciones 4–6)

**M6. No hay teorema sobre el esquema implementado.** (R3, R2)
- **Qué falla:** los Teoremas 5.1, 5.2 y 5.4 son para una versión idealizada. En particular:
  - Velocidad implícita M_n: no se prueba existencia, monotonía ni una cota a priori. Los Teoremas 5.2 y 5.4 son enunciados *a posteriori* (dependen de M_n y de la condición ℓ ≥ 2ξ_N).
  - Remallado: no está cubierto.
  - Campo lejano no rígido: la cota afirmada debe probarse o corregirse.
- **Cambio:**
  - Probar existencia y unicidad de M_n y una cota a priori.
  - Cerrar el remallado. Hay un argumento casi gratis: Φ es cóncava del lado fresco, así que la interpolación lineal favorece a la barrera superior y la inferior pierde sólo O(h_c²) por remallado.
  - Si no se puede, decir explícitamente cuál es la brecha.

**M7. Las cotas no dicen nada sobre los experimentos.** (R3)
- **Qué falla:**
  - Todos los experimentos usan f = (1−s)⁺, que no es C^{1,1} en ℝ, así que el factor de la cota es ρ_f = h_f/ε². Con h_f = ε/16 eso da 1/(16ε), unos 62 para ε = 10⁻³.
  - La hipótesis σΔt ≤ νε termina exigiendo Δt ≲ ε².
  - Ninguna tabla ilustra un teorema en un régimen donde éste sea informativo.
- **Cambio:**
  - Agregar al menos un experimento con una f regular (C^{1,1} en ℝ) y Δt dentro de las hipótesis, que muestre que la cota es ajustada.
  - Discutir explícitamente qué régimen cubren los teoremas y cuál los experimentos.

**M8. "Δt no aparece en el error" tiene que reformularse.** (R3)
- **Dónde:** abstract, introducción, Teorema 5.2, "The mesh drags the front", Sec. 5.5.
- **Qué falla:**
  - En esencia, la afirmación se reduce a que el esquema linealmente implícito deja invariante un estado estacionario en el sistema comóvil.
  - En el esquema implementado, Δt vuelve a entrar por M_n, y el propio experimento de arrastre lo muestra.
  - La corrida independiente de R3 muestra un error residual que crece con Δt: 4.5·10⁻⁴ con Δt = 0.02 y 9.7·10⁻⁴ con Δt = 1, para ε = 0.003. R3 lo atribuye al remallado y al campo lejano.
- **Cambio:** reformular vía la onda discreta comóvil (c_h(ϑ), Φ_h). El efecto de pinning es exactamente c_h(ϑ) − c_ε. El enunciado honesto sería: "el error no depende de Δt *para la onda viajera con malla trasladada rígidamente*".
- **Además:** la afirmación de la l. 1542 ("of order Δt e^{−cR}") no se deduce de barreras de traslación, que sólo dan CΔt·T. Corregirla o probarla.

**M9. La Sec. 5.5 ("What is lost: averaged truncation errors") es heurística.** (R3, R2)
- **Qué falla:** ψ* = e^{cz}Φ' no decae en +∞, y el espectro esencial del linealizado toca el 0. Eso explica la relajación lenta en escala 4/c². Una prueba requiere espacios con peso.
- **Cambio:** rotularla claramente como heurística y acortarla (para SINUM), o convertirla en teorema con espacios con peso (para Math. Comp. o Numer. Math.).

**M10. Originalidad: falta posicionarse frente a la literatura cercana.** (R3, R1)

El manuscrito debe discutir, no sólo citar:

| Trabajo | Por qué es relevante |
|---|---|
| **Nochetto–Paolini–Verdi 1994**, EJAM, ondas viajeras continuas y semidiscretas | El análogo más cercano a las Secs. 5.1 y 5.5. **Es el más importante.** |
| Nochetto–Verdi 1996/1997 | Convergencia al límite con tasa combinada; modelo para M6 y para el Camino A |
| Fife–McLeod 1977 | Origen de las barreras de traslación usadas en la Sec. 5 |
| Patankar; Larrouturou 1991 | Reacción implícita que preserva positividad |
| Smooke–Miller–Kee 1983 (PREMIX) | Fijar el frente para determinar la velocidad |
| Furzeland–Verwer–Zegeling 1990 | Mallas móviles aplicadas a un problema de llama |
| Gracia–O'Riordan 2012 | Capas interiores móviles con métodos ε-uniformes |
| Dwyer–Sanders 1978 y literatura de mallas móviles para llamas 1D desde fines de los 70 (R1) | Muestra que la idea de "seguir el frente" no es nueva en combustión |

R3 verificó en la web NPV 1994, Nochetto–Verdi 1996, Gracia–O'Riordan 2012 y Furzeland–Verwer–Zegeling 1990. **Todas las demás se citaron de memoria y hay que verificarlas antes de incluirlas.** Lo mismo vale para las referencias de combustión que propone R1 (Ferziger–Echekki 1993, Brailovsky et al. 2015, Chen–Ju 2007).

La originalidad defendible no es "mallas que siguen el frente", sino la **combinación de un esquema monótono para todo Δt con mallas ALE y análisis de error en el régimen ε → 0**. El texto debe decirlo así.

**M11. La comparación de costo no es justa.** (R1, R3, R4)
- **Dónde:** Tabla `tab:A3` y el texto asociado.
- **Qué falla:**
  - Se comparan costos sin igualar el error.
  - En los segundos (1.6 s contra 122 s), el BDF de referencia gasta el 41 % del tiempo en factorizar con SuperLU una matriz de 2n × 2n y el 13 % en rearmar el Jacobiano, y además se reinicia en cada remallado.
  - El Jacobiano es triangular por bloques con bloques tridiagonales: la LU (4.4 ms) se puede reemplazar por dos solves tridiagonales (0.32 ms) (R4).
  - La comparación con el esquema euleriano usa Δt/ε hasta 20, con errores de velocidad de hasta 80 %.
- **Cambio:**
  - Reemplazar la tabla por **diagramas trabajo-precisión** (error contra tiempo de CPU o contra número de pasos) frente a malla uniforme fina, BDF con Jacobiano tridiagonal por bloques y un método MMPDE estándar.
  - Reportar pasos y factorizaciones, que es la comparación justa, y no sólo segundos.

**M12. La condición de Péclet.** (R3)
"Para todo Δt y toda malla" depende de la condición de Péclet, que se viola cerca de la extinción. Hay que enunciarla explícitamente, incluida su forma radial, y considerar upwinding híbrido fuera de la banda.

### Interpretación física (afecta credibilidad ante referís de SINUM con formación aplicada)

**M13. No se explica qué es w₀ y "unstable" se presta a confusión.** (R1)
- w₀ es la perturbación de la temperatura de llama en unidades de Zeldovich: es el escalado NEF de Matkowsky–Sivashinsky 1979.
- Con Arrhenius, √(2M) ∝ e^{w₀/2}: la "inestabilidad" del límite es la conocida sensibilidad de la velocidad ZFK a la temperatura de llama.
- La frase de l. 54 ("T_f is the adiabatic flame temperature") sólo vale para w₀ = 0.
- ZFK 1938 y Matkowsky–Sivashinsky 1979 están en el `.bib` pero no se citan.
- **Cambio:** agregar un párrafo de interpretación física en la introducción; usar "sensibilidad" o "no robustez del límite" en lugar de "unstable", o al menos explicar el término.

**M14. La corrección de velocidad (Prop. 2.2) mezcla un término trivial con el cinético.** (R1)
- Se puede escribir 1 − I/√(2M) = −w₀ + ∫(1 − √(B/M)). El término −w₀ viene sólo de que el salto real de temperatura es A + εw₀ y no A.
- El cambio de signo en w₀* ≈ 0.564 (l. 360) es esencialmente ese efecto de normalización.
- **Cambio:** separar los dos términos y dar la interpretación correcta.

**M15. El ejemplo de "quenching" está mal interpretado y, además, no es auditable.** (R1, R4)
- **Física (R1):**
  - La heurística c₀ − 1/r es la relación de Markstein lineal y no vale en r ≈ 0.7δ.
  - En N = 3 existe la bola de llama de Zeldovich con R_Z = 1/c₀ exacto. En 2D no hay bola estacionaria, así que el umbral depende necesariamente del perfil inicial.
  - No se discute por qué r* coincide en 4 cifras para ε entre 0.01 y 0.0025.
- **Código (R4):**
  - `common.py:87-99` captura **cualquier** excepción, y `limit_examples.py:180` la cuenta como "quenching" en la bisección de r*.
  - El JSON no guarda el error de cada paso, así que **r* = 0.715 no se puede auditar** sin volver a correr C3.
  - Según R1, r* depende de Δt, y el texto no lo menciona. El único archivo con otro Δt (`C3_quenching_dt0.005.json`) proviene del código viejo y contiene errores contados como quenching, así que tampoco sirve para evaluarlo.
- **Cambio:**
  - Volver a correr C3 distinguiendo extinción real de falla numérica.
  - Reportar r* en función de Δt.
  - Reinterpretar físicamente el resultado y, si es posible, validar contra la bola de llama 3D exacta.

**M16. Terminología y régimen de parámetros.** (R1)
- "Extinction" designa en el paper el consumo total de una bolsa de gas fresco, no que la llama se apague. Cambiar por "burn-out" o "consumption of the fresh region", o definirlo.
- "Accelerates without bound" (l. 1745) describe la solución exacta, pero el experimento se detiene en t ≤ 0.75.
- Los parámetros son matemáticos y no físicos:
  - número de Zeldovich 100–800, cuando lo realista es 5–20;
  - con ε ≈ 0.05–0.1 la banda de 55ε cubre toda la llama;
  - h_c = 0.01 equivale a unos 170 nodos por espesor térmico, así que la malla "gruesa" no es gruesa.

  Para SINUM no es un problema en sí mismo, pero hay que decirlo y agregar ε = 0.05–0.1 a la Tabla 1 y al menos a un experimento.

---

## 4. Problemas menores, por sección

### Sección 1 (introducción)
- Acortar el párrafo l. 413–424 (comienzo de la Sec. 3), que repite la introducción (R2).
- Agregar un párrafo de motivación física (M13) y mencionar las limitaciones del modelo: Le = 1 (la comparación depende de esa estructura), sin pérdidas de calor, 1D/radial (R1).

### Sección 2
- **Prop. 2.2:** el resto es O(ε²), no O(ε² log(1/ε)) (R2, arreglo de una línea).
  - El Paso 1 usa J ≤ I. En la prueba del Lema 2.1 (l. 255) ya está J(φ) ≤ ½α₀(φ + w₀)².
  - De ahí √(2B) − p = 2σJ/(√(2B) + p) ≤ σα₀t²/(√(2b₀)t) = Cσt, luego I − J(1) ≤ Cσ y el resto es O(ε²), lo que coincide con la Tabla 1.
  - **Borrar l. 363–364** ("the factor log(1/ε) … may be an artefact of Step 1") y corregir el enunciado.
- l. 216 y l. 273: los intervalos [−Zε, 2ε] y "[−Z, 0] en ζ" no coinciden; uniformizar.
- l. 217: Φ'' es continua, así que la cota vale para todo z y no a.e.
- l. 285: Φ ∈ C^{2,1}; puede ir al enunciado.
- l. 393: la reducción a Kummer está verificada (sin cambio).

### Sección 3
- **Marco de comparación (l. 452–461):**
  - El cap. VII de Lieberman es para dominios acotados y Q_T no lo es. Hace falta un argumento tipo Phragmén–Lindelöf para funciones acotadas con coeficiente de orden cero |c| ≤ Cε⁻², por ejemplo con la función auxiliar μ(|x|² + 2Nt + 1)e^{Kt}. Decir que D es acotado en la Prop. 3.5.
  - La regla de la cadena a.e. depende de que {ṽ = εs_m(W)} tenga medida nula; decirlo.
  - Pedir u₀^ε continua, o precisar el sentido de la comparación con datos medibles.
  - Dar una referencia precisa para la frontera parabólica de dominios no cilíndricos.
- **Lema 3.1(ii)** (l. 575–586, 636–637):
  - Citar el teorema de variedad inestable C² con parámetros, que también da la finitud de m(Θ).
  - Escribir las cotas de P_WW, que hoy se omiten con "similarly"; de ellas sale |Θ_WW| ≤ C(1 + |s|)²Θ_s.
  - La fórmula de Θ_WW (l. 607) está verificada.
- **Constantes del Teorema 3.6:**
  - Explicitar que C depende de (δ₀, r, K₀), de la constante de Lipschitz de v en (x, t), de ϱ₀, C_w, C₀ y de T* − T.
  - Destacar que **no crece exponencialmente en T**, porque no hay Gronwall.
  - Decir que v es Lipschitz en t y que su soporte está uniformemente acotado.
- **Orden de elección de las constantes** (Λ → K → c₁ → c* → δ*, ε):
  - Escribir que c* ≤ c₁ (para γ ≥ ε, l. 864), que C_{Λ,K} no depende de c₁ y que c* depende de q̃.
  - Una tabla de dependencias ayudaría mucho.
- l. 443–450: que w₀ sea Lipschitz requiere η_ε → 0; además hace falta η_ε ≤ min{δ₁/2, 1}.
- l. 699: la cota es Eũ ≥ −(3/2)K₀d, no −K₀d. No afecta al −2K₀d de (i).
- l. 739: incluir Λ, K, c₁ en el enunciado de la Prop. 3.4.
- l. 887: formular (q:eq:mono) en t = 0 y deducirla por comparación.
- l. 893–895: destacar que la regularidad hasta t = 0 es una hipótesis de compatibilidad.
- l. 1011: z_ε = log(1 − ε/A)/c_ε = −ε/√(2M) + O(ε²); puede escribirse así.

### Sección 4
- **Lema 4.2 en mallas infinitas:** escribir el argumento, que son dos líneas (R2). Con S = sup Z > 0 y Z_j ≥ S − η, se obtiene (S − η)/Δt ≤ ηC/h_f², y luego η → 0. R3 pide además justificar existencia en ℓ^∞.
- Explicitar la condición de Péclet, también en forma radial (M12).
- **Discrepancias código ↔ texto (R4):**
  - La zona rígida es simétrica (`ale.py:102`): del lado fresco mide 15ε más de lo que dice §4.3. Corregir el texto o el código.
  - La secante acepta en silencio cuando no converge. Cerca de la extinción con Δt grande hace 3.3–4.5 evaluaciones por paso, no 3 (en las corridas guardadas sí son ≤ 3.013). La tolerancia no está escalada con Δt (R3).
  - El intervalo de remallado euleriano está acotado por 0.05·t_end, así que la secuencia de mallas depende del tiempo final. Decirlo o cambiarlo.

### Sección 5
- **Teorema 5.4 (malla de banda):**
  - Exigir que los bordes de la banda sean nodos (R2, R3).
  - En el caso (b) falta un factor e^{2ch_c} (R2).
  - ℓ ≥ 2ξ_N es una condición a posteriori; decirlo.
- **Teorema 5.1:** la asimetría es real en la prueba (Φ'' ≥ −cΦ' en la cota inferior). Aclarar que C_E no depende de v_max.
- **Teorema 5.2:** es correcto y no circular, porque la comparación se aplica a la sucesión realizada v_n. Decirlo.
- Las barreras no satisfacen flujo nulo en 0 ni en L, lo que da errores exponencialmente pequeños en L/ε. Mencionarlo en la Sec. 5.4 (R2).
- Acortar "The mesh drags the front" y parte de la Sec. 5.5 (R2), o reescribirlas según M8 y M9.

### Sección 6 (experimentos)
- La figura `fig:instability` sale de código anterior que no se volvió a correr (R3). Regenerarla con el código actual o quitarla.
- El error de pinning se reduciría trivialmente variando ϑ en cada remallado (R3). Probarlo, o explicar por qué no se hace.
- Faltan tablas de **órdenes de convergencia observados** para datos generales y varios ε (R3).
- Hay choques de notación: (S), M, A, E, n y N se usan con más de un significado (R3). Hacer una pasada de notación.

---

## 5. Qué agregar y qué quitar

### Agregar

**Resultados** (en orden de impacto):
1. El teorema de convergencia del esquema discreto al problema límite con tasa combinada (Camino A; M6 y M10).
2. La eliminación del log(1/ε) mediante un lema de capa inicial (M2), o al menos su bosquejo.
3. Existencia y cota a priori de la velocidad implícita M_n, y cierre del remallado (M6).
4. La Prop. 2.2 con resto O(ε²) (sección 4).
5. La reformulación de la optimalidad (M1) y la clase explícita de soluciones regulares (M4).

**Experimentos:**
1. **w₀(x) variable** (mezcla estratificada). Es justamente lo que distingue este problema del de Caffarelli–Vázquez, y hoy no aparece en ningún experimento (R1). Prioridad alta y costo bajo.
2. Diagramas **trabajo-precisión** frente a malla uniforme, BDF con Jacobiano tridiagonal por bloques y MMPDE (M11).
3. Un experimento con f regular y Δt dentro de las hipótesis, para mostrar que las cotas son ajustadas (M7).
4. Tablas de órdenes observados para datos generales (R3).
5. ε = 0.05–0.1 en la Tabla 1 y en al menos un experimento (R1).
6. *(Deseable)* La bola de llama 3D (Zeldovich, R_Z = 1/c₀) como solución exacta y prueba del radio crítico (R1, M15).
7. *(Deseable)* Una cinética tipo Arrhenius truncada, por ejemplo f(s) = (e^{−s} − e^{−1})⁺ normalizada, para mostrar que nada depende de la elección f = (1−s)⁺ (R1).
8. *(Opcional)* Ignición desde un núcleo con exceso de entalpía (w no constante) (R1).

**Texto:**
- Un párrafo de motivación e interpretación física (M13).
- La relación con Petrosyan 2002 y con el Step II de FBW 2004 (M3).
- El posicionamiento frente a NPV 1994 y el resto de la tabla de M10.
- El marco de comparación en dominios no acotados (sección 4).
- Una tabla de dependencia de las constantes.

### Quitar o acortar
- Las l. 363–364 (el log de la Prop. 2.2 ya no hace falta).
- El párrafo l. 413–424, que repite la introducción.
- "The mesh drags the front" y la parte heurística de la Sec. 5.5, salvo que se conviertan en teorema.
- Si se sigue el Camino A: pasar a un apéndice las pruebas de los Lemas 3.1 y 3.3 (R2).
- La Tabla `tab:A3` en su forma actual (reemplazarla por trabajo-precisión).
- Toda frase que sugiera aplicabilidad general (multi-D, frentes planos para el Teorema 3.6, "quenching" en sentido físico).

### Arbitraje entre referís

R1 pide más física y más experimentos; R2 y R3 piden un paper más corto y enfocado. Mi criterio como editor:
- Los experimentos 1–5 de la lista son **imprescindibles**: son baratos y responden objeciones matemáticas.
- Los experimentos 6–7 son **deseables** si entran sin alargar el paper; si no, mencionarlos en la discusión.
- El experimento 8 y cualquier extensión a Le ≠ 1, pérdidas de calor o multi-D van **sólo como trabajo futuro** en las conclusiones.

El objetivo es un paper de análisis numérico, no de combustión. Pero la física tiene que estar bien dicha, porque un referí de SINUM con formación aplicada va a notar los problemas de M13–M16.

---

## 6. Código y reproducibilidad

R4 corrió todo en una copia aislada, con un venv armado desde `requirements.txt`.

### Lo que funciona
- `pytest`: 19 tests pasan. Tardan 50 s de reloj (17 s de CPU); el README dice "about 10 s".
- `speed_coefficient`: regenera el JSON bit a bit idéntico.
- `convergence` (grupo B, 156 corridas): reproduce `B_convergence.json` con una diferencia máxima de 1.2·10⁻⁸. Tarda unos 10 min; el README dice "~2 min".
- `truncation_check` reproduce −0.2351 / −0.2345 / −0.2342, pero sólo imprime y no guarda en `results/`.
- `reproduce_paper.sh` cubre todas las figuras y tablas. R4 estima 2–3 h en un portátil de 4 núcleos; el README dice 1 h y el script 1.5 h. **Corregir las estimaciones de tiempo en el README y en el script.**

### Problemas a corregir

| # | Problema | Dónde | Cambio |
|---|---|---|---|
| C1 | Excepciones contadas como "quenching"; r* = 0.715 no se puede auditar | `common.py:87-99`, `limit_examples.py:180` | Capturar sólo las excepciones esperadas; guardar en el JSON el motivo de cada paso de la bisección; volver a correr C3 (ver M15) |
| C2 | El fallback puede lanzar el error de Péclet sin capturar; `extinct` no distingue extinción de falla | `ale.py:166` | Separar los estados "extinción", "falla numérica" y "Péclet violado" |
| C3 | Celdas diminutas en los bordes (hasta 1.6·10⁻¹⁰·h_f, cociente entre vecinas 6·10⁹) cuando la banda toca 0 o L con un ancla fuera del retículo | `solver.py:103-129` | Fusionar las celdas menores que un umbral con su vecina |
| C4 | Nada impide que los nodos se crucen; en la extinción radial la celda mínima cae a 0.50·h_f con Δt = 0.05 | `ale.py:125-130` | Validar la malla en cada paso y remallar cuando una celda se comprima por debajo de un umbral |
| C5 | El error del frente inicial corregido en ALE sigue en BDF y en el esquema euleriano | `solver.py:171-178`, `simple.py:60-67` | Aplicar el mismo arreglo |
| C6 | `find_front` devuelve `None` si un nodo vale exactamente ε, y ALE lo interpreta como extinción (comprobado) | `find_front` | Usar una desigualdad no estricta o interpolar; hoy los experimentos lo esquivan desplazando la fase en 1/32 |
| C7 | `set_reaction` es estado global y con el arranque *spawn* de macOS no llega a los procesos hijos | `set_reaction` | Pasar la reacción como parámetro |
| C8 | Defaults con `x or ...` pisan un 0 pasado explícitamente | varios | Usar `x if x is not None else ...` |
| C9 | La docstring de `_KeepResult` (`common.py:63-67`) habla de un `TypeError` que ya no existe. `C3_quenching_dt0.005.json` (que contiene esos `TypeError` del código viejo) y `F_floor_old_code.json` no los produce ningún script actual | `common.py`, `results/paper/` | Corregir la docstring; borrar los JSON o documentar su origen |

### Tests que faltan
- `simple.py`.
- Las cotas del Lema 4.1 y el principio de comparación (Lema 4.2): ya se cumplen numéricamente, falta el test.
- El error de Péclet.
- ALE radial.
- Regresión contra los JSON publicados.

### Rendimiento (medido por R4)
- **ALE** (caso de `tab:A3`, ε = 0.001, 1.85 s): el 68 % se va en armar matrices y sólo el 20 % en resolver. Armar las diagonales con NumPy y resolver con `lapack.dgtsv` da **1.5×** con frentes idénticos.
- **Euleriano:** el 65 % es `solve_banded`, que refactoriza en cada paso una matriz de calor constante. Con el mismo cambio da **3.9×**.
- **Numba:** R4 estima 5–10×. Importa sobre todo para E1 (2.4·10⁶ pasos, unos 25 min).
- **BDF de referencia:** usar el Jacobiano tridiagonal por bloques (ver M11). Es imprescindible para que la comparación de costo sea justa.

### Publicación como material suplementario
- **Falta:** `LICENSE`, `CITATION.cff`, DOI de Zenodo, `pyproject.toml` y versiones fijadas (`requirements.txt` con `==` o un lockfile).
- **Los JSON** deberían guardar la versión de las librerías y el commit con que se generaron.
- **`experiments/paper/README.md`** cuenta la historia de cambios tras referatos internos; reescribirlo como documentación neutral.
- **Derechos de distribución:** revisar los de `ifbj.pdf` y `final.pdf` (artículos publicados en IFB y SIMA) incluidos en el repo, y no subir la carpeta `pdfs/` (otros autores).

**Veredicto de R4:** aceptable para publicar como suplemento, con revisiones menores a moderadas. Lo imprescindible es C1–C5, matizar en el texto M11 y las discrepancias de la sección 4, la licencia, las versiones fijadas y los tests del Lema 4.1.

---

## 7. Plan de trabajo

**Prioridad A:** imprescindible antes de enviar a cualquier revista.
**Prioridad B:** muy recomendable; un referí probablemente lo pida.
**Prioridad C:** deseable.

### Decisión previa
- [ ] **A**: Decidir entre el Camino A (un artículo unificado por el teorema discreto → límite) y el Camino B (dividir). Darse un plazo para intentar el teorema del Camino A antes de optar por dividir.

### Teoría EDP (Secs. 2–3)
- [ ] **A**: M1, reformular la optimalidad en el abstract, la introducción, la Obs. 3.7 y las conclusiones.
- [ ] **A**: M3, discutir Petrosyan 2002 y el Step II de FBW 2004 y delimitar la novedad.
- [ ] **A**: M4, enunciar la clase de soluciones a la que se aplica el teorema y corregir la l. 1821.
- [ ] **A**: M5, corregir la atribución de hipótesis (F₂, monotonía, Arrhenius).
- [ ] **A**: Prop. 2.2 con resto O(ε²) y borrar las l. 363–364.
- [ ] **A**: Marco de comparación en dominios no acotados (Phragmén–Lindelöf).
- [ ] **A**: Cotas de P_WW y cita del teorema de variedad inestable en el Lema 3.1(ii).
- [ ] **B**: M2, intentar eliminar el log(1/ε) con el lema de capa inicial.
- [ ] **B**: Dependencia de las constantes (tabla) y orden de elección.
- [ ] **B**: Revisar Glangetas 1992, Glangetas–Roquejoffre 1996, Langlois–Marion y Brauner–Lunardi.
- [ ] **C**: Los menores de la Sección 3.

### Teoría numérica (Secs. 4–5)
- [ ] **A**: Teorema 5.4, bordes de la banda como nodos y el factor e^{2ch_c}.
- [ ] **A**: Lema 4.2, argumento en mallas infinitas y existencia en ℓ^∞.
- [ ] **A**: M8, reformular "Δt no aparece"; corregir la l. 1542.
- [ ] **A**: M9, rotular la Sec. 5.5 como heurística y acortarla.
- [ ] **A**: M12, enunciar la condición de Péclet, también en forma radial.
- [ ] **A** (Camino A) / **B** (Camino B): teorema de convergencia al problema límite con tasa combinada.
- [ ] **B**: M6, existencia y cota a priori de M_n; cierre del remallado.
- [ ] **C**: Sec. 5.5 como teorema en espacios con peso (sólo si se apunta a Math. Comp. o Numer. Math.).

### Literatura y texto
- [ ] **A**: M10, discutir NPV 1994, Nochetto–Verdi, Fife–McLeod, Furzeland–Verwer–Zegeling y Gracia–O'Riordan. **Verificar cada referencia** antes de citarla (las que no verificó R3 se citaron de memoria).
- [ ] **A**: M13, párrafo de interpretación física de w₀; citar ZFK 1938 y Matkowsky–Sivashinsky 1979.
- [ ] **A**: M14, separar el término −w₀ en la Prop. 2.2 y reinterpretar w₀*.
- [ ] **A**: M16, terminología ("extinction", "accelerates without bound").
- [ ] **B**: Pasada de notación ((S), M, A, E, n, N).
- [ ] **B**: Acortar el párrafo l. 413–424 y "The mesh drags the front".
- [ ] **C**: Mencionar Le ≠ 1, pérdidas de calor y multi-D como trabajo futuro.

### Experimentos
- [ ] **A**: Volver a correr C3 (quenching) con el manejo de fallas corregido; r* en función de Δt; reinterpretarlo (M15, C1).
- [ ] **A**: Diagramas trabajo-precisión con BDF tridiagonal por bloques, malla uniforme y MMPDE (M11).
- [ ] **A**: Experimento con w₀(x) variable.
- [ ] **A**: Regenerar `fig:instability` con el código actual, o quitarla.
- [ ] **B**: Experimento con f regular y Δt dentro de las hipótesis (M7).
- [ ] **B**: Tablas de órdenes observados para datos generales.
- [ ] **B**: ε = 0.05–0.1 en la Tabla 1 y en al menos un experimento.
- [ ] **B**: Variar ϑ en cada remallado para eliminar el pinning, o explicar por qué no.
- [ ] **C**: Bola de llama 3D; cinética Arrhenius truncada.

### Código
- [ ] **A**: C1–C5 (manejo de fallas, celdas diminutas, cruce de nodos, frente inicial en BDF y euleriano).
- [ ] **A**: Alinear código y texto: zona rígida, secante, intervalo de remallado euleriano.
- [ ] **A**: `LICENSE`, versiones fijadas, `CITATION.cff`, DOI de Zenodo.
- [ ] **A**: Tests del Lema 4.1 y del principio de comparación.
- [ ] **B**: C6–C9; tests de `simple.py`, Péclet, ALE radial y regresión contra los JSON.
- [ ] **B**: Guardar en los JSON la versión de las librerías y el commit; `truncation_check` debe escribir en `results/`.
- [ ] **B**: Corregir las estimaciones de tiempo del README y de `reproduce_paper.sh`.
- [ ] **B**: Revisar los derechos de distribución de los PDFs del repo; reescribir `experiments/paper/README.md`.
- [ ] **C**: Optimización (`dgtsv` / factorización reutilizada; Numba).

---

## Anexo: informes completos de los referís

Los informes se reproducen sin cambios, salvo el nivel de los encabezados. Cada referí usa su propia etiqueta M1, M2, …, que **no** coincide con la numeración M1–M16 de la sección 3 de este documento. Recordatorio de numeración: R1 llama "Teorema 3.7" al de tasa, que en el PDF es el **Teorema 3.6**. R2 llama "Teoremas 4.3–4.5" a los de malla fija, móvil y de banda, que en el PDF son los **Teoremas 5.1, 5.2 y 5.4**. R3 llama "Thm 5.3" al de banda, que es el **Teorema 5.4**.


---

---

### A.1. R1, físico de combustión

#### Informe de referí: modelo y aplicaciones

**Manuscrito:** "The high activation energy limit of an equidiffusional flame model: rate of convergence and a front-following numerical method" (J. Fernández Bonder)
**Revista objetivo:** SINUM / Math. Comp. / Numer. Math.
**Perfil del referí:** teoría de combustión (llamas premezcladas, asintótica de alta energía de activación). El informe se concentra en el modelo, su interpretación física y las aplicaciones numéricas. No evalúo en detalle las demostraciones de las Secciones 3–4.

---

##### 1. Qué afirma el paper

1. Para el modelo térmico-difusivo con número de Lewis igual a 1, (1)/(P_ε), con reacción escalada f_ε(s)=ε⁻²f(s/ε) y exceso de entalpía w^ε=Y−u de orden ε (w^ε/ε→w₀>−1), el autor prueba una tasa de convergencia hacia el problema de frontera libre (P), u_t=Δu en {u>0} con |∇u|=√(2M(w₀)). La tasa es −C(ε log(1/ε)+η_ε) ≤ u^ε−v ≤ C(ε+η_ε). Vale si (P) tiene una solución clásica regular y los datos son estrictamente estrellados (Teorema 3.7).
2. Para la onda viajera da una corrección de primer orden de la velocidad, con resto riguroso O(ε² log 1/ε): c_ε = c₀[1+(ε/A)(1−I/√(2M))]. En el caso f=(1−s)⁺ obtiene la constante explícita κ*=(3√3−2)/5 (Prop. 2.2, Tabla 1).
3. Construye las soluciones radiales autosemejantes de (P) a partir de la función de Kummer ₁F₁(−½; N/2; r²/4) (Sec. 2.2).
4. Propone un esquema linealmente implícito con principio de comparación discreto para todo Δt, sobre una malla con banda fina que sigue al frente en formulación ALE con velocidad de malla implícita. Para ondas viajeras prueba cotas de error sin Δt (Teoremas 5.1–5.3).
5. Los experimentos, todos 1D o radiales y con f=(1−s)⁺, muestran:
   - errores (h_f/ε)², independientes de Δt y de ε, para ondas viajeras;
   - un error O(Δt) con constante uniforme en ε para un frente acelerado;
   - convergencia O(ε) hacia un solver de seguimiento de frente para (P);
   - la "inestabilidad" respecto de w₀;
   - la aproximación al perfil autosemejante de Galaktionov–Hulshof–Vázquez;
   - un radio crítico de "quenching" r*≈0.715, independiente de ε.

---

##### 2. Fortalezas

- **El modelo está bien planteado y las constantes son correctas.** La condición |∇u|=√(2M) con M=∫(s+w₀)f(s)ds es exactamente la fórmula de Zeldovich–Frank-Kamenetskii: c₀²A²=2∫ω dT, con A = salto de temperatura. Con f(s)=½e^{−s} y w₀=0 se obtiene M=½, c₀=1 cuando A=1, que es la normalización de Buckmaster–Ludford (prefactor β²/2). Verifiqué la Tabla 1 contra `speed_coefficient.json` y contra c₀=√(2M)/A con M=(1+w₀)³/6: las tres filas cierran. Además c_ε−c₀−(predicción)·ε decae como ε², sin rastro visible del logaritmo, tal como dice el texto (l. 361–364).
- **La Prop. 2.2 aporta algo concreto.** Da un resto riguroso para la corrección de segundo orden de la velocidad y su dependencia explícita en w₀. Se usa bien después: explica el residuo planar ≈0.3ε de la Tabla 4 (l. 1825–1828). Este es el tipo de conexión cuantitativa entre teoría y cálculo que un lector de combustión valora.
- **La interpretación del error en tiempo de los esquemas eulerianos es correcta y elegante.** El retraso del frente proporcional a Δt/ε (y no a Δt/ε², como da la cota) se explica con la condición de solvencia (proyección sobre el núcleo adjunto e^{cz}Φ′), y el coeficiente −0.235 coincide con el medido (l. 1549–1572, 1722–1731). La escala de relajación 4/c² también es la correcta: es el borde del espectro esencial en espacios con peso.
- **Los ejemplos son exactos y no triviales.** La solución autosemejante de Kummer y la comparación con un solver independiente de (P), de segundo orden y validado hasta 6·10⁻⁶ (l. 1790–1792), dan una referencia confiable.
- **El autor es honesto con los límites de su análisis.** Dice que los experimentos están fuera de las hipótesis (l. 1518–1521, 1584–1585) y que las cotas no son óptimas (l. 1404–1407).
- **El código y los scripts de reproducción son públicos** (l. 165–167).

---

##### 3. Debilidades mayores

###### M1. Falta la interpretación física de w₀; la palabra "inestabilidad" es engañosa (l. 75–84, 160–163, 1760–1777, 1761)

Este es el punto más importante para un lector de combustión.

- **Qué es w = Y+T−T_f.** Es el exceso de entalpía total respecto de la referencia T_f. La temperatura adiabática real de la mezcla es T_f+εw₀, no T_f. Por eso la frase de l. 54 ("T_f is the adiabatic flame temperature") sólo es cierta para w₀=0.
- **Qué es w₀.** Es la perturbación de la temperatura de llama medida en unidades de Zeldovich (RT_b²/E). El escalado w^ε=O(ε) es exactamente el escalado NEF (near-equidiffusional flames) de Matkowsky–Sivashinsky (1979) y de Buckmaster–Ludford.
- **Qué dice el resultado de FBW2000 en lenguaje físico.** Con cinética Arrhenius sin corte, f(s)=Ce^{−s}, se obtiene M(w₀)=M(0)e^{w₀}, es decir √(2M) ∝ e^{w₀/2}. Esta es precisamente la condición de salto NEF, [∂T/∂n] ∝ exp(H/2). El resultado es la sensibilidad exponencial ZFK de la velocidad de llama a la temperatura de llama: una perturbación O(1/Ze) de la entalpía cambia la velocidad en O(1). Es material de libro de texto. Explica, por ejemplo, la fuerte dependencia de S_L con la riqueza.
- **Por qué "unstable" es un problema.** En combustión, "inestabilidad" designa casi siempre las inestabilidades del frente: difusivo-térmicas (celulares o pulsantes), Darrieus–Landau, etc. Un lector físico entenderá algo que el paper no hace. Además, con Le=1 el frente plano del modelo térmico-difusivo es incondicionalmente estable.

**Pido:**
1. Reemplazar "instability" por "sensitivity", "non-uniformity of the limit" o "non-commutation of limits", o explicitar que se trata de inestabilidad respecto de los datos en el sentido de Hadamard.
2. Agregar en la introducción un párrafo físico con la interpretación de w₀, el vínculo con NEF y la fórmula ZFK.
3. Citar ZFK 1938, Matkowsky–Sivashinsky 1979 y Zeldovich–Barenblatt–Librovich–Makhviladze 1985. Las dos primeras ya están en `references.bib` pero no se citan.

###### M2. La condición w₀>−1 y la cinética con corte: artefactos del modelo que se presentan como física (l. 56–58, 177–180, 1577)

- **Qué significa w₀>−1.** La temperatura quemada real, T_f+εw₀, tiene que superar la temperatura de corte T_f−ε. Si w₀≤−1, la mezcla es inerte a su propia temperatura adiabática. Esto es un "límite de inflamabilidad" producido por el corte de f en s=1. Con Arrhenius no existe: en el modelo adiabático con Le=1 no hay límite de inflamabilidad sin pérdidas de calor. Hay que decirlo explícitamente.
- **(F) no admite Arrhenius.** La l. 178 afirma que (F) "holds … for Arrhenius-type rates", pero e^{−s} no se anula en [1,∞). Lo correcto es "Arrhenius truncado" o "cinética con temperatura de ignición". Corresponde citar la literatura de cinética con temperatura de ignición o escalonada, por ejemplo Ferziger–Echekki (1993) y Brailovsky–Gordon–Kagan–Sivashinsky (Combust. Flame 2015). Esa literatura muestra que esta cinética reproduce el mismo salto que el límite de alta energía de activación.
- **Los resultados cuantitativos dependen de la cinética elegida.** Todos los experimentos usan f=(1−s)⁺, que da M ∝ (1+w₀)³, una dependencia polinomial muy distinta de la exponencial. Por ejemplo, de w₀=0 a w₀=1 la velocidad crece un factor 2.83 con f=(1−s)⁺ y un factor 1.65 con Arrhenius. Las magnitudes observadas son propias de esta cinética: las posiciones de frente de la Fig. 6, el valor w₀*≈0.564 y r*≈0.715.
- **Sugerencia.** Repetir al menos la Tabla 1 y la Fig. 6 con una cinética más cercana a Arrhenius que cumpla (F₂). Por ejemplo f(s)=(e^{−s}−e^{−1})⁺, convenientemente acotada para s≪0, que es no creciente, Lipschitz y C² hasta 1.

###### M3. La corrección de velocidad mezcla un efecto trivial con uno cinético; el cambio de signo en w₀* es casi todo trivial (l. 301–308, 359–360)

Se puede reescribir 1−I/√(2M) = −w₀ + ∫_{−w₀}^{1}(1−√(B(φ)/M))dφ. Esto separa dos efectos:

- **El término −w₀ es trivial.** Viene sólo de que el salto real de temperatura (= Y fresco) es A+εw₀ y no A. La fórmula ZFK con el salto correcto da c=√(2M)/(A+εw₀) ≈ c₀(1−εw₀/A).
- **El segundo término es la corrección cinética genuina.** Es siempre positivo. Para f=(1−s)⁺ vale (1−κ*)(1+w₀) ≈ 0.361(1+w₀). Su integral converge también para Arrhenius sin corte en s=1, lo que permite comparar con las expansiones clásicas de segundo orden (Bush–Fendell, Buckmaster–Ludford). En la forma actual, I diverge para Arrhenius.

El "cambio de signo" de l. 360 sale entonces de −w₀+0.361(1+w₀)=0, es decir w₀*=0.564. Es esencialmente un efecto de normalización (la elección de A como T_f−T_u), no un fenómeno cinético. **Pido** presentar la Prop. 2.2 en esta forma, comentar el sentido físico de cada término y moderar l. 359–360.

###### M4. El alcance físico del Teorema 3.7 es mucho más estrecho de lo que sugieren la introducción y las conclusiones (l. 101–109, 650–651, 881, 1821, 1888–1893)

- La Definición 3.3 exige Ω_t={ũ>0} acotado y el Teorema 3.7 exige u₀ de soporte compacto. Como {u>0} es la región **fresca**, el teorema sólo cubre **bolsas de mezcla fresca rodeadas de gas quemado** que se consumen en tiempo finito.
- Quedan fuera las configuraciones físicamente más comunes: la llama plana (que el propio autor excluye en l. 1010–1011), la llama en expansión desde un núcleo y el problema de ignición. Todas tienen región fresca no acotada.
- Esto hay que decirlo en la introducción (l. 101–109) y en las conclusiones.
- En l. 1821 se afirma "All the errors are O(ε), in agreement with Theorem 3.7", pero la columna planar de la Tabla 4 no está cubierta por el teorema (soporte no compacto). Hay que reformular.

###### M5. Terminología: "extinction" y "quenching" (l. 1590, 1600, 1744–1758, 1844–1856, 1867–1883)

- **"Extinction".** En combustión significa que la llama se apaga, por pérdidas de calor, estiramiento, dilución, etc. Aquí designa algo distinto: la llama consume completamente una bolsa de gas fresco, es decir el agotamiento de una bolsa de gas fresco ("burnout of a fresh-gas pocket", "collapse of an unburnt pocket"). Es la terminología matemática de FBW2004 y GHV1997. Si se conserva, hay que explicarla. En los pies de figura conviene cambiarla.
- **Por qué el frente acelera.** Corresponde explicar físicamente la aceleración: es la curvatura de un frente cóncavo hacia el gas fresco. Y hay que acotar que (P) deja de ser válido cuando el radio es comparable al espesor térmico o al de reacción.
- **"Accelerates without bound" (l. 1745).** El experimento se detiene en t≤0.75 con T=1, donde la velocidad del frente está acotada (≈R/(2√0.25)). La frase describe la solución exacta, no lo que se testea.

###### M6. El ejemplo de "quenching de núcleos pequeños" está mal interpretado físicamente (l. 1867–1883)

- **La heurística c₀−1/r es la relación de Markstein lineal.** Para Le=1 en el modelo térmico-difusivo, la longitud de Markstein es el espesor térmico δ=1/c₀, y la relación vale para curvatura débil, r≫δ. Aplicarla en r≈0.7δ y decir que "overestimates it by 30%" no tiene contenido: la heurística está fuera de su rango. Hay que citar Markstein, Clavin (1985) o Matalon–Matkowsky (1982).
- **En N=3 el radio crítico es exacto.** El problema (P) tiene la bola de llama de Zeldovich u=A(1−R/r)⁺ con R_Z=A/√(2M)=1/c₀ exactamente, y esa bola es inestable. Es decir, 1/c₀ no es una heurística sino el radio de la bola de llama en 3D. En N=2 no existe bola de llama estacionaria en un dominio no acotado (el potencial logarítmico no tiende a A). El umbral en 2D depende entonces necesariamente del perfil inicial y del dominio (aquí L=12 con flujo nulo, cf. `limit_examples.py`). Esto es justamente lo que observa el autor ("presumably depends on the whole initial profile").
- **Pido.**
  1. Rehacer el ejemplo en N=3 y comparar r* con R_Z=1/c₀. Es un test físico limpio, con referencia exacta, y además un test exigente del método cerca de un equilibrio inestable.
  2. Comparar r* directamente con el solver de (P) de la Sec. 6.5.
  3. Citar la literatura de ignición y bolas de llama: Zeldovich (1944), Buckmaster–Joulin–Ronney (1990), Chen–Ju (Combust. Theory Model. 2007).
- **La independencia de ε es sorprendente y no se discute.** r* es idéntico a 4 cifras para ε=0.01, 0.005 y 0.0025 (`C3_quenching.json`: 0.715234 en los tres casos). Uno esperaría una corrección O(ε), de tipo Markstein de orden superior. Hay que discutir por qué está ausente o es menor que la tolerancia de la bisección (≈8·10⁻⁴), lo que implicaría un coeficiente de la corrección O(ε) por debajo de ≈0.1.
- **Detalle de reproducibilidad.** En `run_kernel` una simulación que termina con error se clasifica como "quenching" (`if err is not None ... outcome = "quenching"`). Hay que justificar ese criterio o separar los fallos numéricos. También hay que mencionar que con Δt=0.005 aparecen corridas "undecided" y un r* que depende de Δt (`C3_quenching_dt0.005.json`: 0.7129 contra 0.7105).
- **"Quenching" no es el término adecuado.** En la literatura esto es "ignition failure" o "critical kernel radius". Además, en este modelo el núcleo no tiene exceso de entalpía (w constante), así que no se modela la energía de ignición. Ver la sugerencia A3.

###### M7. El régimen de parámetros es matemático, no físico, y el paper no lo dice (l. 1057–1066, 1064, 1599, 1630–1635, 1655)

- Los números de Zeldovich realistas son Ze≈5–20, es decir ε≈0.05–0.2. Los experimentos usan ε=0.01–0.00125 (Ze=100–800). El único ε "físico" aparece en la Tabla 4 (ε=0.05).
- **A ε realista la banda cubre toda la llama.** Con ε=0.05 la banda de 55ε mide 2.75, más que el espesor térmico 1/c₀≈1.73 (w₀=0). A ε≈0.1 la banda abarca la llama entera y el método se reduce a una malla fina uniforme sobre la llama. La ventaja del método (costo que crece sólo logarítmicamente en 1/ε) sólo aparece a Ze≫50.
- **La malla "gruesa" no es gruesa.** h_c=0.01 equivale a ≈170 nodos por espesor térmico, y para ε=0.01 coincide con ε. La comparación de costo de la Tabla 2 sería más convincente con h_c≈0.1–0.2 (≈10–20 nodos por espesor de precalentamiento, lo habitual en combustión) y comparando a **igual precisión** con BDF, que corre con rtol=10⁻⁵ y cuyo error no se informa.
- **Esto no invalida el paper,** que es de análisis numérico y apunta al límite ε→0 como herramienta para testear teorías asintóticas. Pero conviene decirlo explícitamente y usarlo como motivación: un método uniforme en ε permite verificar numéricamente predicciones asintóticas (velocidad, Markstein, NEF) en el régimen donde los códigos de β fijo se vuelven caros.
- **Pido** agregar ε=0.1 y 0.05 a la Tabla 1, para mostrar qué tan útil es la corrección de primer orden a Ze realista, y al menos un experimento ALE a ε=0.05–0.1.

###### M8. Cobertura de casos: el modelo con Le=1 es el caso dinámicamente más pobre, y la parte "nueva" del modelo (w no constante) no se testea (l. 441–447, 1577, 1909–1921)

- **Ningún experimento usa w₀ no constante.** Todos usan w^ε≡εw₀, así que w es trivial y M es una constante. Lo que distingue este problema del de Caffarelli–Vázquez (u=Y) es el acoplamiento con un w calórico no constante (hipótesis (W), barreras con Θ(s;W)), y eso nunca se ejercita numéricamente. Físicamente, w₀(x) no constante corresponde a una mezcla estratificada en entalpía o en composición (combustión estratificada) o a un núcleo caliente. Ver A1 y A3.
- **Faltan los fenómenos físicamente relevantes**, y conviene al menos discutirlos con honestidad en las conclusiones:
  1. **Le≠1 (NEF).** Aparecen inestabilidades celulares y pulsantes; las pulsantes existen ya en 1D y serían un test exigente en la geometría del paper.
  2. **Pérdidas de calor.** Dan la extinción verdadera (curva en C, Joulin–Clavin, Buckmaster).
  3. **Multi-D.** Frentes celulares, Kuramoto–Sivashinsky.
  4. **Ignición y bolas de llama.**
- **El principio de comparación discreto depende de la estructura Le=1.** Con Le≠1 o con pérdidas de calor, w deja de ser calórico y el sistema deja de ser cooperativo. Hay que decir explícitamente que la pieza central del método no se extiende de manera directa. Esto es importante para lectores que piensen en aplicar el método.
- **¿Es aceptable que todo sea 1D/radial?** Para SINUM, con teoremas para ondas viajeras y experimentos 1D, sí, siempre que el paper no sugiera aplicabilidad general a llamas. Las mallas móviles en 1D para llamas son estándar desde fines de los 70: Dwyer–Sanders 1978, y los trabajos de malla móvil y adaptativa de Bayliss–Matkowsky, Hyman–Larrouturou y otros. Lo novedoso es el análisis uniforme en ε y el principio de comparación, no seguir el frente en sí. La introducción (l. 117–133) debería reconocer esa literatura de combustión, no sólo Smooke y BGMM.

---

##### 4. Observaciones menores (por línea de main.tex)

- **l. 29:** "DRAFT, \today", a eliminar.
- **l. 47–58:**
  - Declarar el modelo: densidad constante, térmico-difusivo (Zeldovich–Barenblatt), difusividades iguales, temperatura adimensionalizada con el calor de reacción para que T+Y se conserve. Hoy el lector no-matemático no sabe qué se desprecia.
  - "ε is essentially the inverse of the activation energy" → ε=1/Ze, Ze=E(T_b−T_u)/(RT_b²), con longitudes en unidades del espesor térmico D/S_L.
  - Explicar por qué el prefactor es ε⁻²: es el límite distinguido de BL/ZFK que mantiene la velocidad de llama O(1).
- **l. 54:** "T_f is the adiabatic flame temperature" sólo vale para w₀=0 (ver M1).
- **l. 56:** f positiva en (−∞,0) significa reacción a temperaturas superadiabáticas (w₀>0). Conviene comentarlo.
- **l. 65–70:** aclarar que u=Y es el caso adiabático exacto (w=0).
- **l. 77:** decir que (P) es el modelo de "flame sheet" de Le=1, que u≡0 en el quemado (temperatura uniforme T_f) y que por eso es un problema de una fase.
- **l. 82–84:** "unstable" (ver M1).
- **l. 96–99:** el único antecedente "cuantitativo" citado es BNS/JKS. Agregar Bush–Fendell (ya citado más abajo) y ZFK.
- **l. 119:** "We are not aware of numerical studies of the limit ε→0" es demasiado fuerte. Hay estudios numéricos a Ze creciente comparados con asintótica, por ejemplo para llamas pulsantes en Bayliss–Matkowsky. Reformular como "ε-uniform numerical analysis".
- **l. 178:** "and for Arrhenius-type rates": Arrhenius no cumple (F) (ver M2).
- **l. 186–188:** identificar A=T_f−T_u (salto de temperatura en la referencia) y c₀=√(2M)/A como la fórmula ZFK, con cita a ZFK 1938, que está en el .bib pero no se cita.
- **l. 196:** la normalización Φ(0)=ε (frente = isoterma T_f−ε) es una convención. En combustión se usa también el máximo de la tasa de reacción. El corrimiento ε/√(2M) de l. 1011–1013 depende de esa convención y conviene decirlo.
- **l. 359–360:** ver M3.
- **l. 361–364 y Tabla 1:** agregar ε=0.1 y 0.05 (ver M7).
- **l. 381–408:** explicar físicamente la solución autosemejante (bolsa fresca que colapsa, velocidad del frente ∝(T−t)^{−1/2}) y su validez física (radio ≫ ε).
- **l. 650–651, 881:** destacar que {u>0} acotado significa región fresca acotada (ver M4).
- **l. 1005–1009 (orden ε alcanzado):** el ejemplo usa u^ε−v=−εW en el quemado, que es simplemente la diferencia entre la temperatura quemada real T_f+εw₀ y la de referencia, un offset trivial. Para un lector físico, la medida relevante es la posición del frente o la velocidad, no el sup de u. Sugiero reforzar la afirmación con el corrimiento del frente (que ya se menciona) y aclarar el carácter del offset.
- **l. 1054:** definir el frente en términos físicos, como la isoterma T=T_f−ε, el borde frío de la zona de reacción.
- **l. 1059:** la "cola de 20ε del lado quemado" tiene longitud de decaimiento 2/α₀. Decir que depende de f y de w₀.
- **l. 1187:** ℓ₀=max{2,3/c}: c pequeña (w₀ cerca de −1) da zonas rígidas muy largas, lo que es otra manifestación del límite de inflamabilidad artificial.
- **l. 1590, 1600, 1744, 1756, 1844–1856:** "extinction" (ver M5).
- **l. 1623, Tabla 2:** "the number of steps of the adaptive BDF solver grows like ε^{−1.4}". Comparar a igual error (ver M7).
- **l. 1745:** "accelerates without bound": en el rango testeado (t≤0.75) no es así (ver M5).
- **l. 1761:** "instability proved in [FBW2000]" (ver M1).
- **l. 1762–1769:** explicar físicamente: mezclas con Y fresco distinto en O(1/Ze) (algo más o menos de reactivo limitante) y por lo tanto temperaturas adiabáticas distintas en O(RT_b²/E). Mencionar que w₀=−0.5 está a mitad de camino del límite artificial w₀=−1, lo que explica el retroceso inicial y la velocidad 8 veces menor.
- **l. 1773:** la figura sale de `numerics/figures/` y no de `figures/paper/`, con un estilo distinto (título dentro de la figura). Unificar.
- **l. 1789:** cluster en el caso plano; indicar L y el tratamiento del borde x=L.
- **l. 1821–1822:** "in agreement with Theorem 3.7" excede el alcance para el caso plano (ver M4). Además, front*/ε en radial todavía decrece (0.187→0.053), lo que indica que el término O(ε) no está asintóticamente resuelto. Comentarlo.
- **l. 1868–1870:** ver M6. "extinguished by curvature" debería ser "fails to ignite because of stretch/curvature".
- **l. 1871:** r*=0.715±0.001 "for every ε" (ver M6 sobre la sensibilidad a Δt y la clasificación de fallos).
- **l. 1886–1921 (conclusiones):** agregar el alcance físico (Le=1, adiabático, densidad constante), la estructura que hace funcionar el principio de comparación y qué fenómenos quedan fuera (M8).
- **Bibliografía:** no se citan ZFK 1938 ni Matkowsky–Sivashinsky 1979, aunque están en `references.bib`. Sugiero agregar también Sivashinsky (1977, 1983), Matalon–Matkowsky (1982), Clavin (1985), Zeldovich et al. (1985), Williams (1985), Buckmaster (1993) y Matalon (2007, Annu. Rev.). Del lado matemático, Brauner–Lunardi 2000 y Brauner–Hulshof–Lorenzi 2009 (en el .bib) estudian el problema NEF con Le≠1, del que (P) es el caso degenerado. Citarlos situaría bien el trabajo.

---

##### 5. Qué agregar y qué quitar

###### Agregar
- **A1. Un experimento con w₀(x) no constante**, por ejemplo una llama plana que entra en una mezcla con un escalón suavizado de entalpía w₀(x). Se compara con el solver de (P) con M(w₀(x,t)), donde w resuelve la ecuación del calor. Testea la parte del modelo que lo distingue de Caffarelli–Vázquez y tiene lectura física directa (combustión estratificada). Hay que revisar la compatibilidad con (q:eq:mono), que sólo hace falta para el teorema y no para el experimento.
- **A2. La bola de llama de Zeldovich en N=3**, R_Z=1/c₀, como solución exacta adicional (Sec. 2) y como test del radio crítico (M6).
- **A3. Ignición por un núcleo con exceso de entalpía**, w(x,0)/ε = w₀+E·1_{|x|<r₀}, que entra dentro de (W). Daría un análogo genuino de la energía mínima de ignición dentro del modelo.
- **A4. La Prop. 2.2 reescrita** separando el término de normalización −w₀ de la corrección cinética (M3), y una cinética tipo Arrhenius truncada (M2).
- **A5. Un párrafo de motivación física en la introducción:**
  - (P) es el modelo de "flame sheet" para Le=1;
  - la tasa cuantifica el error de reemplazar llamas de Ze finito por láminas;
  - la Prop. 2.2 da la primera corrección en 1/Ze;
  - w₀ es la perturbación de entalpía del NEF;
  - un método uniforme en ε sirve para testear teorías asintóticas.
- **A6. Una tabla de parámetros físicos:** relación de ε con Ze, de h_f y h_c con el espesor térmico, y de las escalas de tiempo de las corridas en tiempos de tránsito δ/S_L.
- **A7. Una discusión de por qué el principio de comparación requiere Le=1** y de qué se perdería con Le≠1 o con pérdidas de calor.

###### Quitar o acortar
- Las afirmaciones sobre "inestabilidad", o redefinirlas.
- La comparación con la heurística 1/c₀ en 2D, o reemplazarla por la versión 3D.
- La Fig. 1 (A5) y la Fig. 2 (A1) se solapan bastante en lo que muestran (la malla sigue al frente). Una sola podría alcanzar y liberaría espacio para A1–A2.
- Los detalles de pinning de red (l. 1674–1691) son interesantes numéricamente pero largos. Se pueden condensar.

---

##### 6. Recomendación

**Revisión mayor.**

**Justificación.** La contribución matemática y numérica parece seria. La tasa de convergencia es nueva para este límite. La corrección de velocidad con resto riguroso es útil y está verificada con precisión. El esquema monótono con banda móvil tiene un análisis original y los experimentos están cuidadosamente hechos. Desde el punto de vista del modelo y de las aplicaciones, en cambio, el manuscrito tiene carencias que deben corregirse antes de publicar:

1. No explica qué significan físicamente w₀, la condición w₀>−1 ni la "inestabilidad". Son, respectivamente, la perturbación de entalpía del NEF, un artefacto de la cinética con corte y la sensibilidad ZFK. Usa además "instability" y "extinction" en sentidos que un lector de combustión interpretará mal.
2. Omite las referencias físicas fundamentales: ZFK, Matkowsky–Sivashinsky, Markstein/Clavin, bolas de llama.
3. Presenta como física efectos que son de normalización (el cambio de signo en w₀*) o de rango de validez (la heurística 1/c₀ en 2D).
4. Sugiere un alcance mayor del que tiene. El teorema de tasa sólo cubre bolsas frescas acotadas, y los experimentos están en un régimen de Ze≫100 con una cinética polinomial.
5. No ejercita numéricamente la parte del modelo con w no constante.

Ninguno de estos puntos requiere nueva teoría. Son cambios de presentación y de interpretación, más dos o tres experimentos adicionales: w₀(x) variable, la bola de llama 3D y un ε realista. Con esos cambios, el paper sería una contribución sólida y atractiva también para lectores con formación en combustión. Que sea 1D/radial es aceptable para SINUM o Numer. Math., siempre que se acoten las afirmaciones como se indica arriba.

---

### A.2. R2, experto en EDPs

#### Informe de referí — "The high activation energy limit of an equidiffusional flame model: rate of convergence and a front-following numerical method"

Foco del informe: Secciones 2 y 3 (estimaciones cuantitativas, barreras, tasa de convergencia) y las demostraciones de la Sección 5 (comparación discreta, teoremas de error), desde el punto de vista del análisis. Los números de línea remiten a `paper/main.tex`.

---

##### 1. Resumen

El manuscrito estudia el límite ε→0 de (P_ε), u_t = Δu − (u+w^ε) f_ε(u) con w^ε calórica y w^ε/ε→w_0, cuyo límite es el problema de frontera libre (P) con |∇u| = √(2M(w_0)). Las contribuciones analíticas son las siguientes.

- **Sección 2.** Lema 2.1: existencia y cotas finas de la onda viajera Φ_ε, con |c_ε−c_0|≤Kε, cotas |Φ^{(k)}|≤Kε^{1−k}Φ', cola exponencial y |Φ_ε−u_0|≤Kε. Proposición 2.2: corrección de primer orden de la velocidad, c_ε = c_0[1+(ε/A)(1−I/√(2M))] + O(ε² log(1/ε)), con el valor explícito κ_* para f=(1−s)^+.
- **Sección 3.** Versión cuantitativa de los Teoremas 4.1/4.2 de FBW 2004. Una barrera explícita V = −w + εΘ(ṽ/ε; W) + ε²q̃ se construye con perfiles internos Θ^± (Lema 3.1), parametrizados solo por W y con energía G_±² = 2M(W) ± 2γ, y con extensiones ṽ^± cuyo gradiente es monótono en la dirección normal (Lema 3.3). La barrera aproxima a la super/subsolución estricta con error O(ε) bajo la única condición ε+η_ε ≤ c_*δ, donde δ es la "strictness". Combinado con el reescalado parabólico v_λ de datos estrellados, esto da el Teorema 3.6:
  −C(ε log(1/ε)+η_ε) ≤ u^ε − v ≤ C(ε+η_ε).
- **Sección 5.** El esquema linealmente implícito satisface un principio de comparación discreto para todo Δt (Lema 4.2). Para la onda viajera hay tres teoremas de error vía barreras trasladadas: malla fija (Teorema 4.3), malla que se traslada (Teorema 4.4) y malla de banda (Teorema 4.5).

Revisé línea por línea las demostraciones de las Secciones 2, 3 y 5. **No encontré errores fatales.** Las cuentas centrales son correctas: la identidad básica (q:eq:basic), las cotas del Lema 3.1, las extensiones monótonas, los márgenes (3.x) y las tres regiones (A)–(C), la verificación de los datos iniciales en el Teorema 3.6 y las barreras discretas. Sí encontré:
- una mejora inmediata (el log de la Prop. 2.2 es un artefacto removible en una línea);
- una afirmación de optimalidad que está probada solo en un sentido trivial;
- huecos de marco funcional (comparación en dominios no acotados, existencia de soluciones clásicas regulares);
- una cuestión de novedad que hay que discutir mejor (Petrosyan 2002, y el propio Step II de FBW 2004).

---

##### 2. Fortalezas

1. **La idea central es correcta y está bien identificada.** Para obtener una tasa hay que conseguir que el umbral ε_0(δ) de las aproximaciones de FBW 2004 sea lineal en δ. Notemos que la cercanía O(ε) de la aproximación ya estaba en FBW 2004, Thm 4.1, Step II (`final.tex`, l. ~690: "sup|u^ε − ũ| ≤ Cε"); lo que no estaba cuantificado era ε_0(δ), escondido en el blow-up por contradicción del Step III. El manuscrito reemplaza el blow-up por una verificación directa con constantes independientes de δ, y esa es la contribución genuina. La Observación 3.7(1) lo explica bien.
2. **La elección técnica de la energía es acertada.** Se parametrizan los perfiles solo por W, con G_±²−2M(W) = ±2γ constante en W. Así se obtienen cotas |Θ_W|, |Θ_sW|, |Θ_WW| uniformes en γ (Lema 3.1(ii)). Las cotas se apoyan en:
   - la acotación de ∫_{Θ_m}^{1+W} t f'(t−W)dt por 2b(Θ_m), que usa la monotonía de f (l. 590–593);
   - el Gronwall "desde arriba" en (b).

   Ambos argumentos son elegantes y correctos.
3. **La extensión con gradiente monótono (Lema 3.3) es simple y eficaz.** Con ṽ^± = Ê ± Λ'(Ê^−)², el margen Λ ε|s| absorbe los términos (1+|s|)Θ_ss que produce la variación de W a través de la capa. Verifiqué las constantes: 8Λ'δ_0³ − 2K_0² = Λ + δ_0K_0, y 4Λ'δ_0K_0² r_1 ≤ δ_0K_0/4.
4. **El Lema 2.1 y la Prop. 2.2 son limpios.** El uso de la energía E = ½p² − B, con dE/dζ = −σp², da en pocas líneas existencia, monotonía, cotas inferiores de p y la expansión de la velocidad. Verifiqué el cálculo de κ_* = (3√3−2)/5 y los valores de la Tabla 1 (columna "Prop. 2.2": 0.13888, 0.20829, −0.45473; todos coinciden).
5. **Sección 5: el principio de comparación discreto para todo Δt es correcto y relevante.** La observación clave es que congelar la tasa f_ε(U^n) con f no creciente hace monótono el esquema (Lema 4.2). Las barreras trasladadas de los Teoremas 4.3–4.5 están bien hechas. En particular, la asimetría del Teorema 4.3 es real en la prueba (para la cota inferior se usa Φ'' ≥ −cΦ', por eso aparece Δt y no Δt/ε²). También es correcto que en el Teorema 4.4 el retraso de la tasa tiene el signo favorable.
6. **Honestidad.** El autor declara qué no está probado: la Sección 5.4 dice que el código está fuera de las hipótesis, el log es posiblemente innecesario, y los teoremas numéricos solo cubren ondas viajeras.

---

##### 3. Errores y huecos

Clasificación: **[F] fatal**, **[S] serio**, **[R] reparable**, **[m] menor**.

###### [S1] La optimalidad del orden ε (l. 109, l. 1005–1014, Conclusiones) está probada solo en un sentido trivial

La introducción afirma sin calificación "the order ε cannot be improved". La Obs. 3.7(3) lo justifica así: en {v=0}, lejos del frente, u^ε = −εW^ε + O(ε²), y por lo tanto ‖u^ε − v‖_∞ ≥ cε si w_0 ≢ 0. El argumento es correcto, pero es una trivialidad: mide la diferencia entre el estado quemado de (P_ε), u = −w^ε, y el de (P), u = 0. No dice nada sobre los siguientes casos.
- **w_0 ≡ 0**, que es el caso u^ε = Y^ε de Caffarelli–Vázquez y el caso más estudiado. Ahí el argumento no da nada.
- **La cantidad relevante (u^ε)^+ − v**, que es la que se mide en la Tabla 3, o la posición de la frontera libre.

El argumento de la onda viajera (desplazamiento ε/√(2M) del conjunto de nivel {u^ε=ε}) depende de qué nivel se elija para definir el "frente". Además, como reconoce el propio texto, la onda viajera no está cubierta por el teorema.

**Corrección propuesta.** Hay dos opciones.
- (a) Reformular: "the order ε is optimal for u^ε − v in L^∞ when w_0 ≢ 0, because of the burnt state", y dejar abierta la optimalidad para (u^ε)^+ y para la frontera.
- (b) Probar una cota inferior genuina. Por ejemplo, para w_0 constante con corrección de velocidad no nula, comparar con barreras construidas con Φ_ε. La Prop. 2.2 da c_ε − c_0 = βε con β ≠ 0 explícito (β = 0 solo si w_0 = w_0^*), y eso daría un error de posición ≥ cε t en geometría plana local.

Como mínimo hay que corregir la frase del abstract/introducción.

###### [S2] El log(1/ε) de la cota inferior: muy probablemente un artefacto de comparar en t=0. Hay una ruta concreta para eliminarlo

El origen del log está correctamente identificado en la Obs. 3.7(2). Como los datos tienen Y_0^ε = 0 en la región quemada, la subsolución debe cumplir Θ ≤ 0 allí. El perfil Θ^−, con P(0) ~ √γ, tarda una longitud s ~ log(1/γ) en llegar de Θ ~ 1 a Θ = 0 (l. 625). Con γ ~ δ ~ ε, esto fuerza un desplazamiento ε log(1/ε).

La asimetría es por lo tanto un artefacto de la construcción, no un fenómeno; los experimentos (Tabla 3) no la muestran. La propia Obs. 3.7(2) dice que con datos Y_0^ε ≥ cε e^{−β dist/ε} la cota inferior es O(ε). Esos datos no cumplen la tercera condición de (q:eq:data), pero eso no es un obstáculo real: la cota superior se establece desde t=0 y se propaga sola. Lo que falta es un **lema de capa inicial**:

- (i) Y = u^ε + w^ε satisface Y_t ≥ ΔY − (C/ε²)Y. En Ω_0, a distancia ≥ ε de ∂Ω_0, Y_0 ≥ cε.
- (ii) Una subsolución explícita de tipo "frente exponencial", εe^{−β(d − V t)/ε} con V ~ 1/ε (que es subsolución si β² > C), muestra que tras un tiempo t_1 = Cε² log(1/ε) vale Y(·,t_1) ≥ cε e^{−β d/ε} para d ≤ Cε log(1/ε). Esa es la única franja que importa, porque más lejos Θ^− < 0.
- (iii) Desde t_1 se aplica el argumento de la Obs. 3.7(2) con v(·, · + t_1). Esto cuesta O(t_1) por la Lipschitz-continuidad en t, y requiere que v(·,t_1) siga siendo estrictamente estrellada, lo que es cierto por continuidad para t_1 pequeño si v es regular hasta t=0.

Esto debería dar −C(ε+η_ε) ≤ u^ε − v. No lo verifiqué en todos los detalles, pero no veo obstrucción. Si el autor no lo hace, la frase "We do not know whether the logarithm is needed" (l. 1002) es aceptable. Aun así, recomiendo fuertemente intentarlo: un resultado simétrico O(ε) es mucho más limpio y elimina la principal debilidad del teorema principal.

Lo que **no** resuelve el log es cambiar γ, porque la condición cδ − 2γ > 0 fuerza γ ≲ δ ~ ε.

###### [S3] Novedad respecto de Petrosyan 2002 y de FBW 2004: hay que precisarla

- FBW 2004 (`final.tex`, l. 242–246) dice que "the simple construction in [CP] of supersolutions of (P_ε) that approximate a strict classical supersolution, when w^ε=0, does not work in the general case". [CP] es Petrosyan 2002, CPDE. Es decir, para w^ε = 0 ya existía una construcción explícita de perfil-compuesto-con-supersolución. El manuscrito cita a Petrosyan solo por unicidad (l. 74, l. 89).

  El autor debe explicar:
  - (a) en qué consiste la construcción de Petrosyan;
  - (b) si ella ya da ε_0(δ) ≥ cδ para w ≡ 0, en cuyo caso la tasa para el caso u=Y sería esencialmente consecuencia de trabajo previo más el argumento de reescalado;
  - (c) qué es exactamente lo nuevo.

  Mi lectura es que lo nuevo es (1) la dependencia en W variable con cotas uniformes en γ y (2) la extensión monótona. Eso es una contribución técnica real pero modesta.
- La afirmación "to our knowledge this is the first rate for the present singular limit beyond travelling waves" (l. 113–114) es plausible. No encontré trabajos que la contradigan para el problema de evolución. Sin embargo, el autor debería revisar y citar o descartar explícitamente:
  - Glangetas (Asymptotic Anal. 5, 1992), sobre el límite singular de alta energía de activación con estimaciones para frentes planos;
  - Glangetas–Roquejoffre (ARMA 1996), sobre bifurcación de ondas en el modelo termo-difusivo;
  - Langlois–Marion (Adv. Differential Equations), con estimaciones uniformes y límite de alta energía de activación con Lewis ≈ 1;
  - Brauner–Lunardi y colaboradores (estabilidad para el problema de frontera libre "NEF").

  Varios de estos tratan estimaciones cuantitativas en contextos de ondas o casi-ondas.
- Respecto de la técnica Allen–Cahn (Chen 1992; Alfaro–Hilhorst–Matano 2008; de Mottoni–Schatzman 1995): allí la "strictness" la provee el corrimiento ±εp(t) de la distancia con signo, y el orden O(ε) es estándar. El manuscrito debería decir con claridad qué diferencia estructural hay aquí:
  - el perfil interno está acoplado al gradiente exterior (energía G²), no a la curvatura;
  - el perfil **no** se ajusta al gradiente local |∇ũ(y)|, sino a una energía G_±(W) con margen cδ;
  - el reescalado de datos estrellados sustituye al corrimiento temporal.

  Una consecuencia, que conviene decir explícitamente, es que la barrera es "gruesa": sirve para O(ε) pero no puede capturar el término de orden ε (el corrimiento ε/√(2M)), a diferencia de las expansiones de orden superior de de Mottoni–Schatzman.

###### [S4] Existencia de soluciones clásicas regulares: el teorema puede ser casi vacuo fuera de w_0 constante

La Obs. 3.7(4) cita a Baconneau–Lunardi 2004 para M constante, datos lisos y compatibles, y tiempos cortos. El Teorema 3.6 está formulado para (W) general, con w^ε no constante y M = M(w_0(x,t)) variable, y para eso no se da ninguna referencia de existencia. Además, "regular classical solution" exige (R1) con |D²Ev| ≤ K_0 y |∂_t Ev| ≤ K_0 **uniformes hasta t=0**. Eso requiere u_0 ∈ C^{1,1}(Ω̄_0) y que se cumpla la condición de frontera libre en t=0 (compatibilidad de orden 1), y quizás algo más para la acotación de v_t cerca de (∂Ω_0, 0).

**Corrección.**
- (a) Enunciar explícitamente la clase no vacía conocida: w_0 constante, datos C^{2+α} compatibles y estrictamente estrellados, intervalo [0,T^*] de existencia de Baconneau–Lunardi, y las soluciones autosemejantes.
- (b) Para w_0 variable, indicar que la existencia local es una extensión esperable pero no probada, o limitar el enunciado.
- (c) Comentar si (S) se preserva en el tiempo. No hace falta para la prueba, solo en t=0, y conviene decirlo, porque en FBW 2004 la geometría se usa de otro modo.

###### [R1] Proposición 2.2: el log(1/ε) en el resto es un artefacto, removible en una línea (l. 313–326, 362–364)

En el Paso 1 se usa J ≤ I, lo que da 2B − p² = 2σJ ≤ 2σI y por eso el corte en t ~ √σ y el logaritmo. Pero en la demostración del Lema 2.1 (l. 255) ya se probó J(φ) ≤ ½α_0 (φ+w_0)² = ½α_0 t². Entonces

  √(2B) − p = 2σJ/(√(2B)+p) ≤ σα_0 t² / (√(2b_0) t) = Cσ t,

de donde 0 ≤ I − J(1) ≤ Cσ, r_ε = O(ε), y el resto de la Prop. 2.2 es **O(ε²)**. Esto coincide con la Tabla 1 (sin log visible). Hay que corregir el enunciado y eliminar la frase de l. 363–364 ("may be an artefact").

###### [R2] Principio de comparación: marco funcional incompleto (l. 452–461)

- (a) El Cap. VII de Lieberman (ABP–Krylov–Tso para W^{2,1}_{N+1}) está formulado en dominios **acotados**. En la Prop. 3.4 el dominio es Q_T = R^N × (0,T). Hace falta un argumento de tipo Phragmén–Lindelöf para funciones acotadas, con coeficiente de orden cero acotado c = [(V+w)f_ε(V) − (u+w)f_ε(u)]/(V−u), |c| ≤ Cε^{−2} para ε fijo. Por ejemplo, se puede usar la función auxiliar μ(|x|² + 2Nt + 1)e^{Kt} con μ→0. Es estándar, pero hay que escribirlo o citarlo bien.

  En la Prop. 3.5 el dominio D sí es acotado, porque fuera de un compacto ṽ ≡ −m(1+Λ'm) < εs_*, y conviene decirlo.
- (b) La regla de la cadena a.e. para D²(Θ∘(ṽ, W)), con Θ ∈ C^{1,1} y C² fuera de {s = s_m(W)}, requiere que {ṽ = εs_m(W)} tenga medida nula. Eso se prueba (l. 805), pero solo "cerca" de ese conjunto; hay que decir dónde se usa.
- (c) La continuidad de u^ε en t=0: (q:eq:data) solo pide desigualdades puntuales. Hay que pedir u_0^ε continua (o, al menos, la comparación vale con datos acotados medibles en sentido L^∞ con la solución mild); hay que decir cuál.
- (d) La definición de frontera parabólica para dominios no cilíndricos (l. 459) está bien como conjunto suficiente, pero conviene citar un enunciado preciso, por ejemplo Lieberman, Lema 2.1/Cor. 7.4, para dominios generales.

###### [R3] Lema 3.1(ii): pasos insuficientemente justificados (l. 575–586, 636–637)

- La afirmación de que "la variedad inestable depende de forma C² de W con P_W = O(Θ), P_WW = O(Θ)" es cierta: el campo (Θ,P) ↦ (P, (b−κP)/G²) es C² cerca de (0,0) porque f ∈ C² en (−∞,1] y M'' = f(−W) es continua. Pero la finitud a priori de m(Θ) = sup|P_W|/θ se apoya en eso, y conviene citar el teorema de variedad estable con parámetros.
- "Differentiating once more ... |P_WW| ≤ CP" y "similarly |P_WW| ≤ CP" en (b) se omiten. Son las cotas que dan |Θ_WW| ≤ C(1+|s|)²Θ_s, que a su vez es el término εΘ_WW|∇W|² en 𝒫. Hay que escribirlas al menos en (a), o moverlas a un apéndice.
- Verifiqué la fórmula de Θ_WW (l. 607): es correcta.

###### [R4] Hipótesis: (F_2) y la monotonía de f se usan menos de lo que se dice (l. 114–115, 440)

- **Si w^ε ≡ εw_0 con w_0 constante** (η_ε = 0, ∇W ≡ 0), entonces 𝒫 ≡ 0 y no se necesitan las cotas en W del Lema 3.1(ii). En ese caso basta (F) (f Lipschitz, no creciente). La frase "For the rate we need f nonincreasing and of class C² up to 1" es entonces falsa para el caso principal, que es el de los experimentos y el de las soluciones autosemejantes. Conviene separar: (F) basta para w_0 constante, y (F_2) se usa solo para controlar ∂_W y ∂²_W.
- **La monotonía de f** se usa en el Lema 2.1 (f ≤ α_0², B ≥ ½f(φ)(φ+w_0)²), en el Lema 3.1(a) (b ≥ ½b(Θ_m) en [(1+W)/2,Θ_m] y la cota de |Θ_m'|) y, sobre todo, en el esquema (Lema 4.2). En la Prop. 3.4 (l. 763–767) **no** se usa, aunque el texto dice "Since f is nonincreasing and Lipschitz". Para la parte de EDP, la monotonía podría reemplazarse por una condición de no degeneración cerca de 1 (por ejemplo, f'(1^-) < 0, o f ~ (1−s)^p). Al menos hay que indicar dónde se usa.

###### [R5] Teorema 3.6: dependencia de las constantes

Hay que explicitar de qué depende C:
- de los márgenes (δ_0, r, K_0) de v en [0,T^*];
- de la constante de Lipschitz de v en (x,t);
- de ϱ_0, C_w, C_0;
- de T^*−T, a través de λ²T ≤ T^*.

Hay que decir que **no** crece exponencialmente en T, porque no hay Gronwall: es comparación pura. Es una ventaja frente a estimaciones de tipo energía y merece mencionarse.

También: en (q:eq:vlam) se necesita que v sea Lipschitz en t globalmente, lo que se deduce de v_t acotada en Ω y v=0 fuera (l. 933), y que el soporte esté uniformemente acotado. Ambas cosas están implícitas; conviene escribirlas.

###### [R6] Sobre la "strictness" y el orden de las constantes (Props. 3.4–3.5)

Verifiqué que el orden Λ → K → c_1 → c_* → δ_*, ε es consistente. Hay que escribir explícitamente:
- en la Prop. 3.5 se necesita c_* ≤ c_1 para tener γ = c_1δ ≥ ε (l. 864);
- C_{Λ,K} no depende de c_1 (l. 860–863);
- c_* depende de q̃ a través de K (región (B) de la Prop. 3.4), algo que sí se dice.

Sugiero un cuadro o lista con "constante: depende de".

###### [m] Menores (Secciones 2–3)

- l. 196–202 y l. 222: en el Lema 2.1(4) el exponente α_0/(2ε) y la constante K absorben e^{α_0Z_*/2}. Correcto, pero conviene decirlo.
- l. 216: "Φ' ≥ θ_Z on [−Zε, 2ε]". En la demostración (l. 273) se escribe "[−Z,0] in the variable ζ"; uniformizar.
- l. 217–218: "for a.e. z, |Φ''| ≤ ..." — Φ'' es continua, así que vale para todo z.
- l. 285–286: Φ ∈ C^{2,1}(R) en general (Φ''' salta en z=0 si f'(1^−) ≠ 0). Correcto; podría ir al enunciado.
- l. 443, 448–450: "Note that w_0 is Lipschitz in x with constant C_W" requiere η_ε → 0 a lo largo de alguna sucesión; además 𝒲 ∋ W requiere η_ε ≤ min{δ_1/2, 1}, no solo η_ε ≤ δ_1/2.
- l. 699: en (q:eq:Eu), la cota inferior correcta es E ũ ≥ −K_0 d − ½K_0 d² ≥ −(3/2)K_0 d para d ≤ 1, no −K_0 d. No afecta (i) (−2K_0d), porque (5/4)(3/2) < 2.
- l. 739: "(for suitable Λ, γ = c_1δ and κ = Kε)" — Λ, K, c_1 son constantes del enunciado; mejor incluirlas en él.
- l. 887: (q:eq:mono) para la función calórica w_0(x,t) se deduce de la condición en t=0 por comparación, como en FBW 2004, l. 430–432. Conviene formularla en t=0 y deducirla.
- l. 893–895: "In particular v(·,0)=u_0 is C^{1,1} up to ∂Ω_0 and satisfies the free boundary condition there": es una **hipótesis de compatibilidad** y conviene destacarla como tal.
- l. 917–919: verificar que (q:eq:data) se cumple para los datos (eq:initdata). Se cumple, con C_0 = |w_0| y Y_0^ε = 0 fuera de Ω_0.
- Obs. 3.7(3), l. 1011–1012: z_ε = log(1−ε/A)/c_ε = −ε/(Ac_ε) + O(ε²) = −ε/√(2M) + O(ε²). Correcto; podría escribirse así.
- l. 393: la corrección del signo en FBW 2004 (4.6) está bien. Verifiqué la reducción a Kummer con a = −1/2, b = N/2.

###### Sección 5 (análisis)

- **[R7] Lema 4.2 en mallas infinitas (l. 1262–1264).** El argumento "evaluar en un punto donde Z está a η del supremo" funciona si los coeficientes off-diagonal están acotados (h ≥ h_f > 0) y F_j Z_j ≥ 0 en ese punto. Hay que escribir las dos líneas: con S = sup Z > 0 y Z_j ≥ S−η, se obtiene (S−η)/Δt ≤ ηC/h_f², y luego η → 0.
- **[R8] Teoremas 4.3–4.5 en R vs. el código en [0,L] con flujo nulo.** Las barreras Φ(x − ct ± ξ) **no** satisfacen flujo nulo en x=0 ni en x=L. En x=L el error es Φ'(L) ~ e^{−cL}; en x=0 es exponencialmente pequeño en L/ε. Hay que mencionarlo en la Sección 5.4, junto con las otras diferencias.
- **Teorema 4.3.** Verifiqué ambas desigualdades. La condición c_0Δt ≤ ε/2 asegura que los puntos de Taylor están a distancia ≤ ε, y que vale la "Harnack" e^K. Nota: C_E no depende de v_max en este teorema (l. 1266); conviene aclararlo.
- **Teorema 4.4.** Correcto. Es un enunciado a posteriori (M_n depende de la secuencia realizada v_n, que a su vez depende de U). Como la comparación se aplica a la sucesión realizada, no hay circularidad. El autor lo dice (Sección 5.4); conviene mencionarlo también en el enunciado.
- **Teorema 4.5.**
  - La clasificación de nodos (a)/(b)/(c) cubre todos los nodos gracias a la holgura h_f en las definiciones, lo que exige que los extremos de la banda sean nodos de la malla. Hay que decirlo.
  - La condición ℓ ≥ 2ξ_N es de nuevo a posteriori (la malla se fija antes que M_N).
  - En (b), la cota |Δ_hΦ − Φ''| ≤ (2/3)h_c sup|Φ'''| sobre el stencil requiere comparar sup Φ' en el stencil con Φ'(z'), lo que da un factor e^{2c h_c}. Hay que escribirlo.
  - q = C_Bε³h_f² es inofensivo.
- **Obs. 4.6(4).** El conteo de nodos y la restricción Δt ≲ 1 frente a Δt ≲ ε³ son correctos.
- **Alcance.** Los resultados de la Sección 5 son esencialmente estabilidad de la onda viajera exacta bajo el esquema. Son correctos pero estrechos: no cubren datos generales, y el propio autor lo reconoce como el problema abierto (ii). Desde el punto de vista del análisis, el resultado más interesante de la Sección 5 es el Lema 4.2.

---

##### 4. Nivel de detalle

- **Sección 2:** adecuado.
- **Sección 3:** densa y en algunos puntos insuficiente: la dependencia C² de la variedad inestable, |P_WW|, "similarly", la regla de la cadena a.e. y la comparación en dominios no acotados. Las cuentas de márgenes y regiones están bien presentadas. Para una revista de EDP el nivel es casi adecuado, salvo [R2]–[R3]. Para SINUM/Math. Comp./Numer. Math. es excesivo para lo que el resto del artículo usa: la Sección 3 no interviene en el análisis del esquema. Solo se invoca para interpretar la Tabla 3 y en el problema abierto (ii).
- **Sección 5:** el nivel es adecuado. Las Secciones 5.4–5.5 (implementación y truncamiento promediado) son honestas y útiles.

**¿Encaja?** La Sección 3 es un resultado de EDP puro (singular perturbation / free boundary) y encajaría mejor en SIMA, Interfaces Free Bound., JDE, Calc. Var. o Ann. IHP (C). Las Secciones 4–6 son de análisis numérico. Hoy el manuscrito son dos artículos cosidos por las soluciones exactas de la Sección 2. Mis recomendaciones, por orden de preferencia:
1. separar en (a) un artículo de EDP (Secciones 2 y 3, con el log eliminado según [S2] y la discusión de novedad [S3]) y (b) un artículo numérico (Secciones 2.1, 4–6, que cita a (a));
2. si se mantiene uno solo, mover las pruebas de los Lemas 3.1 y 3.3 a un apéndice, y reforzar el vínculo con la numérica. Por ejemplo, usar el Teorema 3.6 más una barrera discreta para datos generales (problema abierto (ii)), aunque sea en 1D/radial con malla fija.

---

##### 5. Qué agregar / quitar

**Agregar**
- Discusión precisa de Petrosyan 2002 y del Step II de FBW 2004 (la cercanía O(ε) ya conocida), y de qué es nuevo [S3].
- Prop. 2.2 con resto O(ε²) [R1].
- Lema de capa inicial para eliminar el log, o al menos el bosquejo de [S2] como observación.
- Reformulación de la optimalidad [S1].
- Clase no vacía explícita de soluciones regulares, y la observación de que la tasa no crece exponencialmente en T [S4], [R5].
- Comparación en dominios no acotados [R2].
- Separación (F) vs (F_2) según w_0 sea constante o no [R4].

**Quitar o acortar**
- La frase "the factor log(1/ε) ... may be an artefact of Step 1" (tras [R1]).
- Repeticiones entre la introducción, la Obs. 3.7(1) y la Sección 3 (el párrafo l. 413–424 repite la introducción).
- En la Sección 5, el párrafo "The mesh drags the front" y parte de la Sección 5.5 podrían condensarse.

---

##### 6. Recomendación

**Revisión mayor** (major revision).

**Justificación.** Las demostraciones de las Secciones 2, 3 y 5 son, hasta donde pude verificar línea por línea, correctas en lo esencial; no hay errores fatales. La idea de hacer el umbral ε_0(δ) lineal en δ mediante una barrera explícita, con perfiles parametrizados por W y extensiones de gradiente monótono, es sólida y produce la primera tasa conocida (a mi entender) para el problema de evolución. Pero:
1. la afirmación de optimalidad está sobredimensionada [S1];
2. la asimetría log(1/ε) es muy probablemente evitable y su eliminación haría el teorema principal mucho más limpio [S2];
3. la novedad respecto de Petrosyan 2002 y FBW 2004 no está bien delimitada [S3];
4. la existencia de las soluciones a las que se aplica el teorema solo está asegurada en un caso particular [S4];
5. hay huecos reparables de marco funcional [R2]–[R3].

Hay además un problema editorial de fondo: la Sección 3 es un resultado de EDP que no alimenta el análisis numérico, y una revista de análisis numérico puede considerar que esa mitad está fuera de su ámbito. Sugiero al editor pedir al autor que justifique la unidad del artículo, o que lo divida.

---

### A.3. R3, analista numérico

#### Informe de referí: parte numérica (Secciones 4–6)

**Manuscrito:** "The high activation energy limit of an equidiffusional flame model: rate of convergence and a front-following numerical method" (J. Fernández Bonder), `paper/main.tex`, 1926 líneas.
**Perfil del referí:** análisis numérico de EDPs, esquemas monótonos, mallas móviles y ALE, problemas singularmente perturbados, aproximación de fronteras libres.
**Foco:** el método (Sec. 4), los teoremas de error (Sec. 5), los experimentos (Sec. 6), la originalidad y la adecuación a SINUM, Math. Comp. y Numer. Math. Leí el manuscrito completo, incluida la Sec. 3, que comento sólo donde toca la parte numérica.

---

##### 1. Resumen

El trabajo estudia el sistema equidifusional (P_ε), reducido a `u_t = Δu − (u+w)f_ε(u)` con `w` calórica, y su límite de frontera libre (P). Tiene dos partes:

1. **Parte analítica (Sec. 2–3).** Hay una tasa O(ε) por arriba y O(ε log 1/ε) por abajo para la convergencia a soluciones clásicas regulares con datos estrictamente estrellados (Thm 3.5). También hay una corrección de primer orden de la velocidad de la onda (Prop. 2.2).
2. **Parte numérica (Sec. 4–6).** El esquema (S) (líneas 1083–1086) es linealmente implícito: la tasa de reacción `f_ε(U^n)` se congela y el reactivo `Y^{n+1}` es implícito. Bajo una condición de Péclet de celda cumple un principio de comparación discreto para todo Δt (Lemas 4.1–4.2). Se usa sobre una malla con banda fina (h_f) que se traslada rígidamente con el frente (ALE). La velocidad de malla se determina implícitamente, imponiendo que el frente termine cada paso en un punto de referencia.
   - Para la onda viajera 1D exacta se prueban tres cotas de error: malla fija (Thm 5.1), malla uniforme trasladada (Thm 5.2) y malla de banda no uniforme (Thm 5.3).
   - En las dos últimas Δt no aparece en la cota; aparecen M_n (error acumulado de la velocidad de malla) y ρ_f = h_f/ε² (o h_f²/ε³ si f ∈ C^{1,1}).
   - Los experimentos usan h_f = ε/16, fuera de las hipótesis. Muestran errores O((h_f/ε)²) para ondas viajeras y O(Δt), uniforme en ε, para un frente acelerado.

El manuscrito es honesto sobre sus límites: la Sec. 5.4 enumera las brechas entre teoría e implementación y la 5.5 dice qué se pierde. Las demostraciones de la Sec. 4–5 que verifiqué son, salvo detalles, correctas. El problema central es otro. **El contenido de análisis numérico es delgado** frente a lo que esperan las tres revistas:
- los teoremas cubren sólo la onda viajera exacta en 1D;
- con la f de todos los experimentos, f(s) = (1−s)^+, que no es C^{1,1}, las cotas son vacías en todos los experimentos;
- el resultado más publicitado ("Δt no aparece") es en esencia la invariancia de estados estacionarios bajo Euler linealmente implícito en el sistema comóvil;
- no se analiza el esquema realmente implementado (velocidad implícita, remallado, campo lejano no rígido).

---

##### 2. Fortalezas

1. **El principio de comparación discreto (Lemas 4.1–4.2) es limpio y correcto.** Verifiqué los signos.
   - La matriz `A_F = Δt^{-1}I − 𝒟^{n+1} + diag F` es una Z-matriz con sumas de fila `Δt^{-1}+F_i > 0`, luego M-matriz. Esto vale porque Δ_h y D tienen filas de suma cero, lo que verifiqué también en `ale.py::_matrix` y `_derivative_coeffs`.
   - El argumento para supersoluciones usa `F_U ≥ F_V` (f no creciente) y `V+W ≥ 0`. El de subsoluciones usa `Y^{n+1} ≥ 0`, del Lema 4.1. La asimetría está bien identificada: las subsoluciones no necesitan condición de signo.
   - La observación de la l. 1163 es correcta y útil: el tratamiento totalmente implícito de `(u+w)f_ε(u)` exige Δt ≲ ε², porque `g'` puede valer −Cε^{-2}.
2. **Los Teoremas 5.1–5.3 son correctos en lo esencial.** Revisé línea por línea los desarrollos de Taylor, la cota de Harnack (Lema 2.1(3)) y el signo del término de reacción retrasado.
   - En el Thm 5.1, la supersolución hereda el término desfavorable `𝓑 ≤ C(Δt/ε²)Φ'`, mientras que la subsolución sólo ve `(Φ'')^- ≤ cΦ'`, de donde sale la asimetría Δt/ε² contra Δt. Esta asimetría es una observación fina.
   - En el Thm 5.2, el corrimiento `d_n ≥ 0` de la barrera respecto de los nodos hace que el retraso de la reacción tenga el signo favorable. Es la idea clave y está bien explotada.
   - En el Thm 5.3, la clasificación en (a) banda, (b) lado fresco con Φ''' = c²Φ' sin ε y (c) lado quemado con cola exponencial, más el término q = C_B ε³h_f², es correcta.
3. **La identificación de la mesa de vibraciones del frente (lattice pinning) es valiosa.** La Fig. D3 muestra la función g(ϑ) de la velocidad respecto de la posición del frente en la celda, con amplitud ∝ (h_f/ε)². También la medición del retraso −0.235 Δt/ε en malla fija, contrastada con el argumento de resolubilidad de la Sec. 5.5 (ec. (5.10)). Es un trabajo cuidadoso.
4. **El código está disponible y es legible.** Coincide con lo descrito en las Sec. 4.1–4.3: la condición de Péclet se verifica en cada paso (`ale.py`, l. 54–55), el remallado se ancla a la red de la banda y la secante tiene tolerancia 10^{-3}h_f. Reproduje parcialmente el experimento D1 y amplié su rango:

   | ε | Δt | Δt/ε | E (velocidad implícita) | E (v = c_ε prescripta) |
   |---|---|---|---|---|
   | 0.01 | 0.02 | 2 | 4.5·10^{-4} | 1.4·10^{-4} |
   | 0.01 | 1.0 | 100 | 8.6·10^{-4} | 7.2·10^{-4} |
   | 0.003 | 0.02 | ≈ 6.7 | 4.5·10^{-4} | 4.4·10^{-5} |
   | 0.003 | 1.0 | ≈ 333 | 9.7·10^{-4} | 8.3·10^{-4} |

   La robustez del método ALE para la onda viajera es real. Con Δt grande se ve un crecimiento suave del error con Δt, atribuible al remallado y al campo lejano no rígido.
5. **El texto es honesto.** Las Secciones 5.4–5.5 y las Conclusiones reconocen las brechas. Eso facilita el trabajo del referí, aunque no las cierra.

---

##### 3. Problemas mayores

###### M1. Las cotas de los teoremas son vacías en todos los experimentos, no sólo "fuera de hipótesis" (l. 1257–1260, 1423, 1523–1527, 1584–1585)

- Todos los experimentos usan f(s) = (1−s)^+ (l. 1577). Esta f no está en C^{1,1}(ℝ), porque f'(1^-) = −1 ≠ 0 = f'(1^+). Lo relevante es entonces ρ_h = h/ε², no h²/ε³.
- Con h_f = ε/16 queda ρ_f = 1/(16ε). Para ε = 10^{-3}, ρ_f ≈ 62. La cota del Thm 5.3 es entonces K·62·T, mientras que el error observado es ~5·10^{-4}.
- La hipótesis `σ_nΔt ≤ νε` del Thm 5.3 exige Δt ≲ ε², exactamente la restricción que el método pretende evitar. Lo mismo vale para `ℓ ≥ 2ξ_N`, que pide bandas de ancho O(1).
- Hay además una incongruencia con la introducción. Las líneas 141–149 presentan la cota C(ε + M_T + (ρ_f+h_c)T) como el resultado y dicen que Δt "no aparece". Pero, para la f usada, garantizar un error O(ε) exige h_f ≲ ε³. Con esa malla, la banda de ancho O(ε log(1/ε)) tiene O(ε^{-2} log) nodos: no es el método que se usa.

**Corrección sugerida:**
- (a) Declarar en la introducción y el resumen que los teoremas no cubren ningún experimento y por qué.
- (b) Mejor aún, sustituirlos por la versión aguda (ver M3). Como mínimo, añadir un experimento en el régimen del teorema, por ejemplo ε = 0.05 con h_f = ε³, y una f ∈ C^{1,1} como f(s) = ((1−s)^+)², para ilustrar cuán pesimista es la cota.

###### M2. "Δt no aparece en el error" es esencialmente la invariancia de estados estacionarios, y la hipótesis lo oculta parcialmente (l. 37–39, 141–149, 1334–1409, 1209–1214)

**La explicación elemental.** En el sistema que se mueve con la banda a velocidad v ≡ c_h, la onda viajera discreta (si existe) es un estado estacionario U* de `−𝒟U* + f_ε(U*)(U*+εw_0) = 0`. Euler linealmente implícito con coeficiente congelado deja invariante a U* para todo Δt, porque `f_ε(U^n)U^{n+1} = f_ε(U*)U*`. Por lo tanto el error se reduce a:
- el error de la onda discreta, |U* − Φ| y |c_h − c_ε|, que es un problema estacionario no lineal de autovalor, como en PREMIX o en el método de congelamiento;
- la estabilidad, uniforme en Δt, de U* bajo el esquema.

La demostración por barreras de traslación del Thm 5.2 es correcta, pero no explicita esta estructura. El resumen y la introducción venden como fenómeno de discretización temporal lo que es una propiedad trivial de cualquier método de un paso aplicado a un estado estacionario.

**Dónde está el contenido no trivial.** Está en (i) la condición `σ_nΔt ≤ νε` (estabilidad uniforme en Δt cuando v_n ≠ c) y (ii) el término M_n.

**Dónde vuelve a entrar Δt.** En el esquema implementado, M_n depende de Δt: el experimento de arrastre (l. 1714–1720) muestra que la velocidad del frente es c + α(Δt/ε)δv. La frase "the error of the ALE scheme does not contain Δt (Theorem 5.2)" (l. 1211–1212) es entonces engañosa: lo que se prueba es una cota a posteriori en términos de una cantidad que sí depende de Δt.

**Corrección sugerida:**
- (a) Reformular el Thm 5.2 en dos pasos: existencia y error de la onda discreta comóvil (c_h(ϑ), Φ_h), y después la estabilidad uniforme en Δt.
- (b) Moderar el resumen y la introducción.
- (c) Discutir la relación con la literatura de ondas discretas, que ya se cita parcialmente (Keener 1987; Thümmler 2008). Agregar Elmer–Van Vleck 1999, que ya está en el .bib, y la literatura sobre propagación en redes. El "pinning" g(ϑ) es exactamente c_h(ϑ) − c_ε.

###### M3. No hay teorema sobre el esquema realmente implementado: velocidad implícita, remallado, campo lejano (Sec. 5.4, l. 1516–1545)

La Sec. 5.4 lo reconoce, pero las brechas son sustanciales.

1. **Velocidad implícita.** No se prueba que la ecuación escalar `x_f^{n+1}(v) = x_c + Δt v` (l. 1192–1194) tenga solución, ni que sea única, ni que la secante converja. No se prueba tampoco ninguna cota a priori para M_n.
   - El bootstrap obvio no funciona. Las barreras dan |x_f^n − x_0 − c t_n| ≤ ξ_n, y como el punto de referencia sigue al frente, |v_n − c|Δt ≤ ξ_n + ξ_{n+1}. Eso da una cota de |v_n − c| de orden ξ/Δt, que explota.
   - Pero el teorema necesita la variación total M_n = Δt Σ|v_k − c|, no el desplazamiento neto. Ese es el punto técnico a resolver.
   - Sugerencia: probar la monotonía de v ↦ x_f^{n+1}(v), lo que parece accesible con el Lema 4.2 aplicado a dos mallas comparables o mediante un argumento de traslación. Junto con un argumento de contracción, eso daría |v_n − c_h| ≤ C ρ^{n} + (error de la onda discreta).
2. **Remallado.** El argumento se puede cerrar casi gratis, y conviene hacerlo.
   - En el lado fresco (z > 0), Φ'' = −cΦ' < 0: Φ es cóncava, así que ΠΦ ≤ Φ y la interpolación lineal es favorable para la barrera superior.
   - Para la barrera inferior se pierde |ΠΦ − Φ| ≤ (h_c²/8)|Φ''| ≤ (h_c²c/8)Φ'. Eso se absorbe corriendo la barrera C h_c² en cada remallado. Con O(T/ℓ_0) remallados, la pérdida total es O(h_c² T/ℓ_0).
   - En el lado quemado, Φ + εw_0 ≤ Kεh_f^4 y la pérdida se absorbe en q.
   - Esto convertiría la frase de la l. 1535–1538 ("not covered by the proof") en un lema.
3. **Campo lejano no rígido (l. 1540–1545).** La afirmación de que el error es "of order Δt e^{-cR}" no se sigue de barreras de traslación.
   - Para un nodo con velocidad ẋ_i ≠ v_n, el residuo es ½(ẋ_i − c + σ)²Δt|Φ''| ≤ CΔtΦ', relativo a Φ'. Una barrera de traslación debe absorberlo con un corrimiento uniforme σ ≥ CΔt, lo que da un error CΔt·T y no CΔt e^{-cR}.
   - Obtener e^{-cR} requiere una barrera no de traslación, por ejemplo Φ(z+ξ) más un término de la ecuación del calor localizado en el lado fresco.
   - Corregir la afirmación o probarla.

**Corrección sugerida:** al menos los puntos 2 y 3 son factibles con las herramientas del manuscrito y deberían incluirse. El punto 1 es el que haría del Thm 5.3 un resultado sobre el método.

###### M4. Alcance: sólo ondas viajeras 1D. Falta un teorema para datos generales, y la parte analítica ya da la herramienta (l. 1244, 1909–1912)

Para SINUM, Math. Comp. o Numer. Math., una estimación de error sólo para la solución exacta de onda viajera de un problema 1D es insuficiente como contribución central de análisis numérico. Lo es aún más cuando las cotas no son agudas y no cubren los experimentos.

Lo paradójico es que el propio manuscrito tiene el ingrediente para algo mucho más fuerte:
- el principio de comparación discreto (Lema 4.2) vale en cualquier dimensión, sobre mallas que satisfagan la condición de M-matriz;
- la Sec. 3 construye barreras explícitas, perfil interior compuesto con extensiones monótonas, con márgenes del orden de KΘ_s independientes de la estrictez δ.

El teorema natural sería, al menos para el esquema euleriano (S) en malla uniforme de ℝ^N o en simetría radial:

> Bajo las hipótesis del Thm 3.5, si h, Δt ≤ c ε², entonces |U^n_j − v(x_j, t_n)| ≤ C(ε log(1/ε) + η_ε + (h + Δt)/ε^α) para algún α explícito.

La demostración consistiría en verificar que las barreras V^ε de las Props. 3.3–3.4, evaluadas en la malla, son super/subsoluciones discretas:
- los errores de truncamiento serían ~ (h/ε² + Δt/ε²)Θ_s, más el término de reacción retrasado, con el signo controlado como en el Thm 5.1;
- esos errores se absorben en el margen KΘ_s o reduciendo la estrictez δ.

Esto daría, por primera vez, **convergencia de un esquema totalmente discreto al problema de frontera libre (P), con tasa combinada en ε, h y Δt**. Es el análogo exacto del programa de Nochetto–Paolini–Verdi y Nochetto–Verdi para Allen–Cahn → flujo por curvatura media, con barreras discretas y acoplamientos h, Δt = o(ε^k). Sería un resultado claramente de nivel SINUM o Numer. Math.

Una segunda opción, más ambiciosa, es probar la versión ALE con datos generales. El esquema no es estacionario en el sistema de la banda, y se espera O(Δt) con una constante uniforme en ε, como sugiere la Fig. B4. Probar esa uniformidad sería el verdadero resultado distintivo del método.

###### M5. La "versión aguda" (Sec. 5.5) es heurística, y su marco tiene un obstáculo espectral que no se menciona (l. 1547–1572, 1722–1731)

- La fórmula λ = −∫τψ*/∫Φ'ψ* con ψ* = e^{cz}Φ' es correcta: verifiqué que ψ* está en el núcleo del adjunto ∂² − c∂ − g'_ε.
- Pero ψ* no decae en z → +∞: tiende a (A−ε)c. Además, el operador linealizado 𝒜 = ∂² + c∂ − g'_ε tiene, en el lado fresco donde g'_ε = 0, espectro esencial {−k² + ick} que toca el origen. No hay brecha espectral en L^∞ ni en L² sin peso.
- Esto explica la relajación lenta en "the time scale 4/c²" (l. 1724) y el transitorio de la Fig. D6. También implica que el corrector Ψ con 𝒜Ψ = τ + λΦ' no está acotado uniformemente en espacios sin peso.
- Una prueba debe trabajar en espacios con peso e^{az} (Sattinger) o usar la estructura de tipo ignición.
- El manuscrito afirma que "the experiments agree" con λ = −0.235 Δt/ε, pero el valor medido está entre −0.216 y −0.235 y la diferencia se atribuye a términos de orden superior sin cuantificarlos.

**Sugerencia:** o bien convertir la Sec. 5.5 en un teorema, por ejemplo para la velocidad de la onda discreta (M2(a)): |c_h(ϑ) − c_ε| ≤ C(h/ε)² en malla comóvil y, en malla fija, c_h − c_ε = λ_1 Δt/ε + O(...); o bien reducirla a una observación breve.

###### M6. Originalidad: hay solapamientos y omisiones bibliográficas relevantes

La introducción (l. 117–133) y la Sec. 4.3 posicionan el método de manera incompleta. En particular:

1. **Reacción linealmente implícita que preserva positividad.** Tratar implícitamente el consumo con la tasa congelada es la linealización de fuentes de Patankar (Patankar, *Numerical Heat Transfer and Fluid Flow*, 1980). Es estándar para preservar positividad en química y combustión: Larrouturou, *J. Comput. Phys.* 95 (1991), "How to preserve the mass fractions positivity…"; los esquemas Patankar modificados de Burchard–Deleersnijder–Meister, *Appl. Numer. Math.* 47 (2003). La novedad del Lema 4.2 es el principio de comparación para el sistema acoplado con f monótona, que es modesta (un argumento de M-matriz de tres líneas). Hay que decirlo y citar.
2. **Barreras de traslación Φ(x − ct ± ξ(t)) ± q.** Son las super/subsoluciones de Fife–McLeod, *Arch. Rational Mech. Anal.* 65 (1977). Los Thms 5.1–5.3 son su versión discreta, y debe citarse.
3. **Condición de fijación del frente para determinar la velocidad.** Fijar la temperatura en un punto para que la velocidad de llama sea incógnita es el método clásico de PREMIX (Smooke–Miller–Kee, *Combust. Sci. Tech.* 34 (1983), "Determination of adiabatic flame speeds by boundary value methods"; Kee–Grcar–Smooke–Miller 1985). Se cita Beyn–Thümmler 2004 (freezing), pero no la tradición en combustión.
4. **Mallas móviles para propagación de llamas.**
   - Dwyer–Sanders (1978), problema de llama usado como test estándar.
   - Furzeland–Verwer–Zegeling, *J. Comput. Phys.* 89 (1990), que compara tres métodos de malla móvil (incluidos Petzold y MFE) precisamente sobre un problema de llama.
   - Miller (MFE), Coyle–Flaherty–Ludwig (1986), Petzold (1987).
   - Estos son los competidores naturales del método ALE propuesto y deberían citarse y, idealmente, compararse (ver M7).
5. **Capas interiores móviles con métodos uniformes en el parámetro.** Gracia–O'Riordan, *IJNAM* 9 (2012), "A singularly perturbed convection-diffusion problem with a moving interior layer", y trabajos siguientes, que construyen mallas de Shishkin que siguen una capa móvil con análisis ε-uniforme. Ahí la posición de la capa la da la característica del problema reducido. Aquí es incógnita, pero la conexión es directa y la afirmación de la l. 127–129 ("layers at known positions") debe matizarse. También Shishkin–Shishkina (2009) y la monografía de Linß (LNM 1985, 2010).
6. **Límite singular con h, Δt y ε acoplados.**
   - Además de Nochetto–Paolini–Verdi 1994 y Feng–Prohl 2003, son directamente relevantes Nochetto–Verdi, *Numer. Math.* 74 (1996), "Combined effect of explicit time-stepping and quadrature for curvature driven flows", y Nochetto–Verdi, *SINUM* 34 (1997), "Convergence past singularities…". Usan barreras discretas y principios de comparación discretos con escalas h², Δt = o(ε^4) y tasas de interfaz explícitas: el programa que sugiero en M4.
   - Sobre todo, Nochetto–Paolini–Verdi, *Eur. J. Appl. Math.* 5 (1994), "Continuous and semidiscrete travelling waves for a phase relaxation model": estima asintóticamente, en términos del parámetro de relajación y del paso Δt, las ondas viajeras de un problema de Stefan con relajación de fase y de su contraparte semidiscreta en tiempo. Es el análogo más cercano a la Sec. 5.1 y 5.5 (retraso ~ Δt/ε). Debe citarse y compararse.
7. **Método front-fixing para (P) (Sec. 6.5, l. 1781–1792).** Es la transformación de Landau (1950) y el método clásico para Stefan 1D (Crank, *Free and Moving Boundary Problems*, 1984). Citar. Sobre mallas comóviles para problemas de frontera móvil hay también trabajos recientes ("comoving mesh method", arXiv 2106.00202) que conviene revisar.

Con estas referencias, la originalidad numérica se reduce a la combinación de cuatro elementos y a su análisis para la onda viajera: (i) esquema Patankar-monótono; (ii) banda rígida ALE; (iii) velocidad por fijación del frente; (iv) el estudio cuidadoso del límite ε → 0. La combinación es razonable y el estudio de uniformidad en ε es valioso, pero el texto debe posicionarlo con precisión.

###### M7. Los experimentos no comparan costo contra error con métodos estándar (Tabla A3, l. 1619–1640; Sec. 6.2)

1. **La Tabla A3 compara costos sin igualar errores.** Compara ALE con Δt = 0.01 (error ~5·10^{-4}) contra BDF con rtol = 10^{-5}, pero no da el error de BDF. Hace falta un diagrama trabajo-precisión (error contra tiempo de CPU o contra número de resoluciones lineales) para:
   - ALE;
   - Euleriano (S);
   - BDF sobre la banda;
   - un método de malla móvil estándar (MMPDE5/6 con monitor de arco o curvatura, Huang–Russell);
   - una malla uniforme fina con BDF o Euler implícito con Newton (no monótono).
2. **La comparación con Euleriano (S) en B1 favorece artificialmente al método propuesto.** Los errores relativos de velocidad llegan al 80% con Δt/ε = 20 (README, tabla B1). Ningún usuario correría así. La comparación honesta es a error fijo.
3. **Faltan órdenes observados para datos generales.** Para datos generales sólo está la Fig. B4 (frente acelerado, O(Δt)). Para los datos de la Sec. 6.5 debería darse una tabla de convergencia en Δt y h_f con EOC, a ε fijo y para varios ε, que muestre la uniformidad en ε.
4. **Falta la convergencia en el error combinado ε + h + Δt contra (P).** La Fig. C1 lo insinúa, pero no hay tabla de EOC en una norma y una escala de refinamiento simultáneo, por ejemplo h_f = ε/16 y Δt ∝ ε.
5. **El error de pinning es fácil de reducir y no se reduce.** El error dominante del ALE para ondas viajeras es g(ϑ) (l. 1674–1691), de amplitud (h_f/ε)² y media de 1–3% de la amplitud. Cambiar ϑ en cada remallado, de forma aleatoria o equidistribuida, o elegir ϑ en un cero de g, lo reduciría en un orden de magnitud. Se dice "we have not pursued this"; debería hacerse, porque es trivial e impacta directamente las tablas.

###### M8. La condición de Péclet es restrictiva cerca de la extinción y no se discute el upwinding (l. 1096–1105, 1209)

- "Para todo Δt y toda malla" (l. 1114, 1142, 1895) es condicional: requiere (4.2), es decir |ẋ_i| max h ≤ 2.
- Con convección centrada, cerca de la extinción radial la velocidad del frente diverge como R/(2√(T−t)). En la zona rígida, los nodos de transición y gruesos (h hasta h_c) violan la condición.
- El código lo maneja así: si una velocidad de prueba viola Péclet, la secante devuelve ∞ y se cae a la velocidad previa (`ale.py`, l. 149–175). El texto lo atribuye al "cut-off near r=0" (l. 1854–1856), pero no informa si hubo violaciones de Péclet.

**Corrección sugerida:**
- Informar si hubo violaciones de Péclet en los experimentos.
- Considerar un esquema híbrido: upwind donde |ẋ|h > 2, centrado en la banda. Mantiene la M-matriz incondicionalmente, al precio de difusión numérica O(|ẋ|h) fuera de la banda, que el análisis de la zona (b) del Thm 5.3 tolera (truncamiento O(h_c) sin ε).
- Así "para toda malla" sería literalmente cierto.

---

##### 4. Problemas menores

1. **Conflictos de notación.**
   - (S) designa el esquema (l. 1082) y la hipótesis de estrellado (l. 900).
   - M es M(w_0) y también M_n, M_T (l. 1346).
   - A es la amplitud y también A_F (l. 1118).
   - E es la energía (l. 242), el error (l. 1647) y C_E.
   - n es el número de nodos (x_n = L, l. 1036) y el índice temporal.
   - N es la dimensión y también el número de pasos (l. 1357).
2. **l. 1100–1102.** La condición de Péclet se escribe sólo para N = 1. Dar la forma radial explícita, con los pesos x_{i±1/2}^{N−1}/|K_i|.
3. **l. 1261–1264 (mallas infinitas).** Hay que justificar la existencia y unicidad de la solución acotada de (S) en ℓ^∞: invertibilidad por dominancia diagonal y serie de Neumann. El argumento de comparación "evaluate at a point where Z is within η of its supremum" requiere Σ_j|a_ij| acotado uniformemente, cierto en malla uniforme. Escribirlo.
4. **Thm 5.3, clasificación de nodos.** Requiere que los extremos de la banda (5.7) sean nodos. Si no lo son, un nodo con y_{j−1} ∈ (3ε+ℓ−2h_f, 3ε+ℓ−h_f) no cae en (a) ni en (b). Añadir la hipótesis.
5. **l. 1299 y l. 1259.** Decir explícitamente que la opción ρ_h = h²/ε³ excluye f = (1−s)^+, porque f ∈ C^{1,1}(ℝ) obliga a f'(1) = 0.
6. **l. 1213–1214 y l. 1748–1750.** "with a constant that does not grow as ε → 0": es un único ejemplo (autosemejante, N = 2, 3) con ε ∈ {0.01, 0.002, 0.001}. Moderar.
7. **l. 1230.** "ε^{-1} to ε^{-3/2}" frente a "ε^{-1.4}" en la l. 1623: unificar.
8. **l. 1255.** El factor (1 + Ac_0t) en (5.1) merece una línea de explicación: viene de |c_ε − c_0|t·sup Φ'.
9. **Fig. `fig:instability` (l. 1773).** Usa `numerics/figures/tw_instability.png`, generado por código anterior y no recalculado (según `experiments/paper/README.md`). Regenerarla con los scripts de `experiments/paper/` por reproducibilidad.
10. **Tabla A3.** Los segundos de CPU en Python no son una medida de costo comparable entre un solver tridiagonal propio y el BDF de SciPy. Reportar también resoluciones lineales y evaluaciones del lado derecho, como en el README, y los errores.
11. **l. 1191–1196.** La tolerancia 10^{-3}h_f de la secante introduce un error de velocidad ~10^{-3}h_f/Δt por paso. Para Δt pequeño (0.0025, ε = 10^{-3}) eso es ~2.5·10^{-5}, comparable al error de pinning. Comentar o usar una tolerancia relativa a Δt.
12. **l. 1092–1094.** El ALE se escribe en forma no conservativa. Δ_h es conservativo, pero el término ẋD no lo es: la masa de Y no se conserva exactamente. Mencionarlo, y mencionar que no se requiere la ley de conservación geométrica (GCL).
13. **l. 1658–1660.** "M_T = 7.7–8.3·10^{-4}" es mayor que E. Dado que M_T entra en la cota con constante K, comentar la relación cuantitativa.
14. **Lema 2.1(3) (l. 217–220).** Precisar que Φ''' existe a.e. y que tiene un salto en z = 0 para f = (1−s)^+. El Thm 5.1 usa sup|Φ'''| en un entorno de 0.
15. **Prop. 2.2 y Tabla 1.** El comentario "no visible logarithm" (l. 363–364) es razonable; podría añadirse ε = 10^{-5} para reforzarlo.

---

##### 5. Qué agregar y qué quitar

**Agregar (por prioridad):**
1. Un teorema de convergencia del esquema totalmente discreto para datos generales al problema límite (P), con tasa combinada ε + h + Δt, combinando el Lema 4.2 con las barreras de la Sec. 3 (M4). Es lo que más cambiaría el nivel del trabajo.
2. La reformulación del Thm 5.2 vía la onda discreta comóvil, con una estimación |c_h(ϑ) − c_ε| ≤ C(h/ε)² (M2, M5), es decir, hacer rigurosa la Sec. 5.5 al menos para la velocidad.
3. Los lemas de remallado y del campo lejano (M3.2–3), que son casi gratuitos, y un resultado sobre la velocidad implícita: existencia, monotonía y cota de M_n (M3.1).
4. Diagramas trabajo-precisión frente a MMPDE (Huang–Russell), BDF en banda y malla uniforme (M7.1).
5. Tablas de EOC en Δt y h_f para datos generales, para varios ε (M7.3–4).
6. La aleatorización o equidistribución de ϑ en el remallado (M7.5).
7. Las referencias de M6.

**Quitar o acortar:**
- La Sec. 5.5, si no se convierte en teorema: reducirla a un párrafo.
- La Sec. 6.6 ("Two examples from the literature"): es interesante, pero no aporta al análisis numérico; trasladarla a material suplementario o a un artículo aplicado.
- La Sec. 6.3 (instabilidad con respecto a w_0), calculada con BDF, ilustra un resultado de 2000; puede ir a material suplementario.
- Si se busca una revista de análisis numérico, considerar **separar el trabajo**: la Sec. 3 (tasa del límite singular, con la Prop. 2.2) es un artículo de EDP autocontenido y probablemente la parte más sólida; la parte numérica, fortalecida según M4, sería otro artículo.

---

##### 6. Adecuación a cada revista

**SINUM.** Es la más natural para el método y las Sec. 4–5. En la forma actual:
- la parte numérica es demasiado estrecha: sólo la onda viajera exacta, cotas no agudas y vacías para los experimentos, y sin teorema sobre el esquema implementado;
- casi la mitad del manuscrito (Sec. 2–3) es análisis de EDP sin contenido numérico, lo que SINUM tolera mal en un artículo cuyo título promete numérica.

Con M4 (convergencia totalmente discreta para datos generales), M2 y M3, sería un candidato serio. Veredicto: **revisión mayor**, con un margen apreciable de rechazo si no se agrega un resultado para datos generales.

**Math. Comp.** Exige profundidad técnica en el análisis numérico. Los Teoremas 5.1–5.3 son correctos pero elementales: barreras de Fife–McLeod discretas y un principio de comparación de M-matriz. La versión aguda, que sería lo técnicamente interesante, queda heurística. Veredicto: **no adecuado en la forma actual**. Lo sería si se probaran la Sec. 5.5 (velocidad de la onda discreta con tasa aguda (h/ε)² y Δt/ε) y la convergencia para datos generales.

**Numer. Math.** Valora el análisis agudo y las estimaciones explícitas en los parámetros, en la tradición de Nochetto–Verdi. El paralelismo con el límite Allen–Cahn → curvatura media es atractivo, pero hoy el artículo no tiene el análogo de esos teoremas (convergencia de esquemas totalmente discretos con acoplamientos explícitos h, Δt = o(ε^k) para datos generales). Veredicto: **no adecuado en la forma actual**; podría serlo con M4 y M5.

**Alternativas.** La Sec. 3 sola es publicable en una revista de EDP o fronteras libres: *Interfaces Free Bound.*, *J. Differential Equations*, *Calc. Var.*, *Ann. IHP (C)*, *Nonlinearity*. La parte numérica, tal como está, encajaría mejor en *IMA J. Numer. Anal.*, *J. Comput. Phys.* (con más comparaciones) o *J. Sci. Comput.*

---

##### 7. Recomendación

**Revisión mayor (SINUM); rechazo en la forma actual, con invitación a reenviar tras ampliar el análisis (Math. Comp., Numer. Math.).**

- **Lo que está bien.** Las demostraciones numéricas son correctas, salvo los detalles de M3.3 y de los menores 3–4. El principio de comparación es limpio. La observación sobre el signo favorable del retraso de la reacción en una malla que sigue al frente es la idea central y es buena. Los experimentos son cuidadosos y el código es reproducible.
- **Lo que falta.** Primero, que el análisis cubra algo más que la onda viajera exacta en 1D y que se refiera al esquema implementado. Segundo, que las cotas tengan contenido para las mallas que se usan; hoy, con f = (1−s)^+ y h_f = ε/16, son vacías. Tercero, que la originalidad se posicione con precisión frente a Patankar/Larrouturou, Fife–McLeod, PREMIX/freezing, las mallas móviles para llamas (Furzeland–Verwer–Zegeling), los métodos para capas interiores móviles (Gracia–O'Riordan) y Nochetto–Paolini–Verdi 1994 (ondas viajeras semidiscretas).
- **El camino más corto a un artículo de nivel SINUM** es probar la convergencia de (S) al problema (P) para datos generales con las barreras de la Sec. 3 (M4), y separar o condensar el resto.

---

###### Verificaciones realizadas por el referí
- Signos de los Lemas 4.1–4.2, de la condición de Péclet (l. 1101) y de las entradas de D y Δ_h (contra `ale.py::_matrix`, `_derivative_coeffs` y `solver.py::laplacian`).
- Desarrollos de Taylor y cotas de Harnack de los Thms 5.1–5.3. Las desigualdades de los casos (a), (b) y (c) del Thm 5.3 cierran, salvo el detalle de los extremos de la banda.
- El núcleo del adjunto ψ* = e^{cz}Φ', que no decae en +∞, y el espectro esencial que toca el 0 (M5).
- Una corrida propia del esquema ALE desde Φ_ε (f = (1−s)^+, w_0 = 0, T = 8, h_f = ε/16, h_c = 0.01), con Δt ∈ {0.02, 0.2, 0.5, 1}, ε ∈ {0.01, 0.003}, velocidad implícita y v = c_ε prescripta. Resultados en la tabla de la Fortaleza 4. Es consistente con el argumento de estado estacionario (M2) y muestra un efecto residual O(Δt) del remallado y del campo lejano (M3). Script: `/private/tmp/claude-501/-Users-julianfernandezbonder-Repos-combustion/f1d4d1b7-7ef1-4a81-98c9-d895295d116b/scratchpad/t1.py`.

###### Fuentes bibliográficas verificadas en la web
- [Nochetto–Paolini–Verdi, Continuous and semidiscrete travelling waves for a phase relaxation model, EJAM 5 (1994)](https://cambridge.org/core/journals/european-journal-of-applied-mathematics/article/continuous-and-semidiscrete-travelling-waves-for-a-phase-relaxation-model/20E5B8BD9508B07DA07886C9125EB559) ([registro](https://air.unimi.it/handle/2434/187298))
- [Nochetto–Verdi, Combined effect of explicit time stepping and quadrature for curvature driven flows, Numer. Math. 74 (1996)](https://math.umd.edu/~rhn/papers/abstracts/explicit.html); [Convergence past singularities…, SINUM 34 (1997)](https://math.umd.edu/~rhn/papers/abstracts/implicit.html)
- [Gracia–O'Riordan, A singularly perturbed convection-diffusion problem with a moving interior layer, IJNAM (2012)](https://ojs.global-sci.org/index.php/ijnam/article/view/10079)
- [Furzeland–Verwer–Zegeling, A numerical study of three moving-grid methods for 1D PDEs (CWI / JCP 1990)](https://ir.cwi.nl/pub/1640/)
- [Comoving mesh method for certain classes of moving boundary problems (arXiv 2106.00202)](https://arxiv.org/pdf/2106.00202)

---

### A.4. R4, ingeniero de software científico

#### Informe de referato: código `numerics/` (paquete `combustion1d`) que acompaña a *The high activation energy limit of an equidiffusional flame model: rate of convergence and a front-following numerical method*

Revisor: ingeniería de software científico (Python/NumPy/SciPy, álgebra lineal rala, reproducibilidad).
Revisé el commit `52c3129` en una copia del repositorio, sin modificar el repo original. Entorno: venv con `numpy 2.5.5`, `scipy 1.18.1`, `matplotlib 3.11.2`, `pytest 9.2.0` (lo que instala hoy `requirements.txt`), Python 3.12.5, macOS, Intel i5-1030NG7 (4 núcleos), `OMP_NUM_THREADS=1`.

---

##### 1. Resumen

El código es chico (unas 1000 líneas en el paquete), está escrito con cuidado y en general se corresponde bien con la Sección 4 del artículo. Verifiqué fórmula por fórmula el Laplaciano de volúmenes finitos, el operador `D` de tres puntos, los signos del término ALE (`+ẋ D`), la condición de Péclet, la estructura linealmente implícita de (S), el dato inicial (4.x) `initdata`, la velocidad de malla implícita por secante (tolerancia `1e-3 h_f`), los criterios de remallado (`ℓ0/4`, `15ε/4`), el anclaje de la banda en el retículo, la transferencia por interpolación lineal y el solver de frente fijo del problema límite (ecuaciones transformadas, nodos fantasma, velocidad del frente). No encontré errores de signo ni de *upwinding*. Las cotas del Lema 4.1 se cumplen numéricamente: `Y ≥ 0`, `W ≡ εw0` exacto y `max U` no crece, incluso con Δt = 0.05 en la extinción radial.

Los tests pasan: **19 passed in 50.01 s** (17 s de CPU). Ejecuté dos scripts del artículo: `speed_coefficient` regenera `speed_coefficient.json` **bit a bit idéntico** con NumPy/SciPy actuales, y `truncation_check` reproduce exactamente los −0.2351/−0.2345/−0.2342 de (5.x). También corrí `convergence` (grupo B completo, 156 corridas), que reproduce `B_convergence.json` (detalle en §6). Los tres ejemplos funcionan.

Los problemas que encontré son de **robustez** (casos borde del mallador y del detector de frente), de **manejo de fallas** en los scripts (una excepción cualquiera se clasifica como "quenching"), algunas **discrepancias menores con el texto** y la **preparación para publicación**: no hay licencia, ni DOI, ni versiones fijadas, y quedan archivos huérfanos y referencias al proceso de referato. En rendimiento, el esquema ALE ya es barato, pero la mayor parte de su tiempo se va en ensamblar matrices, no en resolverlas. El solver BDF de referencia está penalizado por la implementación (LU rala genérica de 2n×2n), y eso afecta la comparación en segundos de la Tabla `tab:A3`.

---

##### 2. Fortalezas

- **Correspondencia con el texto.** `ale.py:33-57` implementa exactamente `D` y las entradas extradiagonales de `-𝒟^{n+1}` de la Sección 4.2: comprobé a mano `-(2+ẋh_{i-1})/(h_i(h_{i-1}+h_i))` y `-(2-ẋh_i)/(...)`. La verificación M-matriz (`ale.py:55`) es la condición (4.peclet) tal cual, chequeada en cada ensamblado. Como las filas suman `1/Δt + F > 0`, el chequeo de signo de las extradiagonales alcanza para garantizar la M-matriz.
- **Diseño del esquema.** La reacción con tasa congelada y reactivo implícito, dos solves tridiagonales por paso y `W = Y − U` resolviendo el calor discreto: todo está implementado de manera que el Lema 4.1 vale *por construcción*. Lo comprobé: `W/ε ∈ [0.500000, 0.500000]` en todas las instantáneas.
- **Bandas en retículo.** El anclaje de la banda nueva en el retículo de la anterior (`ale.py:194-196`, `solver.py:100-107`) hace exacta la transferencia en la banda. Es una idea simple y eficaz, y está bien documentada.
- **Instrumentación.** Están los contadores `n_trials`, `n_solves`, `n_lu`, `n_fev` y `v`, más el callback `monitor`. Con eso los números de costo del artículo se pueden auditar: en el caso de `tab:A3`, ε = 0.001, obtuve 400 pasos, 1202 pruebas y 2404 solves, igual que la tabla.
- **Soluciones exactas independientes.** Shooting con bisección por energía, perfil Φ_ε con cola exponencial exacta y Kummer para la solución autosemejante. Están testeadas contra sus EDOs.
- **Determinismo.** No hay aleatoriedad, así que no hacen falta semillas. `ProcessPoolExecutor.map` preserva el orden.
- **Documentación de usuario.** `numerics/README.md` es muy bueno: incluye tabla script ↔ figura, guía de parámetros y limitaciones conocidas.

---

##### 3. Problemas mayores

###### M1. El mallador puede generar celdas arbitrariamente pequeñas en los extremos del dominio
`solver.py:103-113` y `solver.py:117-129`.
Cuando la banda toca `x = 0` o `x = L` y el ancla del retículo no es múltiplo de `h_f`, queda un resto `a = band[0] > 0` (o `L − band[-1]`) menor que `h_f`. `_graded(length, …)` con `length < h0` devuelve `[0, length]` (la fusión de la última celda sólo actúa si `len(d) > 2`), es decir, una celda del tamaño del resto. Lo medí:

| caso | min h / h_f | cociente máximo entre celdas vecinas |
|---|---|---|
| ancla = 0.3 h_f, frente cerca de 0 | 0.30 | 3.3 |
| ancla = 1e-13 | **1.6e-10** | **6e9** |
| frente cerca de L, ancla = 0.3 h_f | 0.70 | — |

En el ALE el ancla es `x[exact[0]]`, arbitraria respecto de 0, así que el caso es alcanzable. Pasa cuando una banda se acerca a `r = 0` (núcleos que colapsan en C3, extinción hasta el final). Las consecuencias son pérdida de cuadratura, un número de condición enorme y que el chequeo de Péclet salte espuriamente.
*Probablemente no afecta a los números del artículo*: en C2 la corrida termina en `s ≈ max(0.15, 30ε)` y la banda del lado fresco mide 20ε. Aun así es un bug del paquete publicado.
**Fix:** en `_graded`, si `length < h0` (o `< 0.5 h0`), devolver `[0]` y absorber el resto en la celda de banda adyacente (`band[0] = 0` o `band[-1] = L`). En general, imponer `min(h) ≥ 0.5 h_f` y agregar un test con ancla aleatoria y frente cerca de cada extremo.

###### M2. Las fallas se tragan y se clasifican como resultados físicos
- `experiments/paper/common.py:87-99`: `simulate_ale_safe` captura **cualquier** `Exception` y devuelve el resultado parcial más `repr(ex)`.
- `experiments/paper/limit_examples.py:180-181`: `if err is not None or … : outcome = "quenching"`. Cualquier error numérico entra en la bisección del radio crítico r_* como "quenching": un Péclet violado, un `TypeError`, un fallo de `solve_banded`. En el JSON almacenado `C3_quenching.json` no quedan registrados los errores de cada paso de bisección (`history` sólo guarda el outcome), así que no se puede auditar si alguno de los r_* = 0.715 se decidió por una excepción.
- La docstring de `_KeepResult` (`common.py:63-67`) dice que `simulate_ale` "raises a TypeError when the front disappears". Eso ya no es cierto: `ale.py:170-174` hace `extinct = True; break`. Además, el archivo huérfano `results/paper/C3_quenching_dt0.005.json` contiene exactamente esos `TypeError` del código viejo.
- `ale.py:166-167`: si las dos pruebas de la secante fallan, el *fallback* `advance(v)` puede lanzar `RuntimeError` (Péclet) sin capturar. Y `res.extinct` (`ale.py:170`) no distingue "el frente desapareció" de "no se encontró el frente porque falló la iteración".

**Fix:** una excepción propia (`PecletError`), un campo `res.status ∈ {ok, extinct, failed}` y que C3 clasifique como quenching **sólo** cuando `extinct` es verdadero o el radio cae por debajo de un umbral, abortando ante cualquier otra falla. Hay que guardar `error` en cada entrada de `history` y re-correr C3 para confirmar que ningún paso de bisección se decidió por excepción.

###### M3. No se verifica la validez de la malla móvil
`ale.py:125-130`.
Los nodos se mueven con `v χ(x)`. En las zonas de decaimiento de χ y en el *taper* de los extremos (`ell_end = min(0.5, 0.1L)`, pendiente 2) las celdas se comprimen. Si `Δt |v| |χ'| ≥ 1` (en el *taper*, `Δt|v| ≥ 0.25`) los nodos **se cruzan**. En la extinción radial `v ~ R²/(2r)` crece sin cota cerca de `T`. El propio `experiments/paper/README.md` (issue 4) y el comentario `convergence.py:28-29` admiten que "the cells at the tapered end x = L collapse" en modo rígido, y por eso limitan la distancia recorrida en B3/D4. En mis corridas de extinción (N = 2, ε = 0.01) la celda mínima fue 0.80·h_f (Δt = 0.002), 0.76·h_f (Δt = 0.01) y 0.50·h_f (Δt = 0.05). La compresión de la banda crece con Δt. Hoy la única defensa es indirecta: una malla invertida suele dar extradiagonales positivas y el chequeo de Péclet salta, pero eso no está garantizado (el término `ẋD` puede compensar).
**Fix:** chequear `np.diff(xn).min() > c·h_f` en `advance` y **remallar** (no fallar) cuando una celda se comprime por debajo de, digamos, `h_f/2`, o cuando el cociente entre celdas vecinas supera `2·growth`. Documentar el criterio en el texto.

###### M4. La comparación de costo en segundos de la Tabla `tab:A3` mezcla algoritmo e implementación
`solver.py:195-207`.
Para ε = 0.005 (BDF, T = 4, 3615 nodos) medí 34.8 s, repartidos así:
- **41 %** en la factorización SuperLU (`gstrf`, 2695 LU de 7230×7230, unos 5 ms cada una);
- 8 % en `SuperLU.solve`;
- **13 %** en reconstruir el Jacobiano con `sp.bmat` en cada llamada;
- el resto, *overhead* de `scipy.sparse` y de `solve_ivp`, que además se **reinicia en cada remallado** (62 reinicios, cada uno vuelve a orden 1).

El Jacobiano es triangular superior por bloques con bloques tridiagonales: `[[Lap − diag d, −diag f], [0, Lap]]`. Dos solves en bandas cuestan **0.32 ms**, contra **4.4 ms** de `splu` (unas 14× menos).
Las columnas "steps" y "LU" de `tab:A3` son una métrica justa. La columna de segundos (ALE 1.6 s contra BDF 122 s) exagera la ventaja por razones de implementación.
**Fix:** aclarar en el texto que los tiempos son de implementaciones Python no optimizadas y que la comparación relevante son pasos y factorizaciones. Si se quiere mantener los segundos, usar `solve_ivp(method="LSODA", lband=2, uband=2)` con incógnitas intercaladas `(u_i, w_i)`, que dan un Jacobiano pentadiagonal, o un BDF propio con solves por bloques.

###### M5. Discrepancias concretas entre el texto y el código
1. **Zona rígida** (`ale.py:102`). `R1 = (max(K_burnt, K_fresh) + K_margin)·ε + rigid_extra·ℓ0` se mide **simétricamente** desde `x_c`. Del lado fresco la zona rígida es entonces (20+15)ε + ℓ0 y no "banda + ℓ0" = (5+15)ε + ℓ0 como dice §4.3: son 15ε de más de ese lado. Es inocuo, pero conviene escribir lo que el código hace, o usar `R1` asimétrico.
2. **"In all our runs it takes three evaluations"** (§4.3). En los JSON almacenados es cierto (pruebas/paso ≤ 3.013 en B y D). Pero si la secante no converge en 8 iteraciones, `ale.py:152-163` acepta `v1` **en silencio**. Con Δt mayores en la extinción medí 3.3 pruebas/paso (Δt = 0.01) y 4.5 (Δt = 0.05), con un fallo de Péclet dentro de la iteración. Sugiero contar las no-convergencias (`res.n_secant_fail`), reportarlas y matizar la frase. Detalle de eficiencia: la tolerancia nunca se chequea sobre `g0`, así que siempre se hacen al menos 3 evaluaciones aunque `v_{n-1}` ya sea solución.
3. **Remallado euleriano** (§4.4: "every time the front has travelled 7.5ε"). En `simple.py:75` y `solver.py:202` el intervalo es además `≤ max_dt_remesh = 0.05·t_end` por defecto. **La secuencia de mallas depende entonces del tiempo final** (para ε = 0.02, `7.5ε/c ≈ 0.37 > 0.05·T`). Hay que documentarlo o quitar el tope implícito.
4. **Frente inicial en BDF y en (S) euleriano.** El bug que se corrigió en ALE (`ale.py:93-96`) sigue presente en `solver.py:171-178` y `simple.py:60-67`: `fronts[0]` y el primer `v_est` usan el frente localizado en la grilla auxiliar `linspace(0, L, 20001)` (paso L/20000, comparable o mayor que ε para ε ≤ 10⁻³). Afecta al primer punto de las series `s_n` (p. ej. la primera ventana de E1) y al primer intervalo de remallado. **Fix:** `xf = find_front(x, u, eps)[0]` después de construir la malla, como en ALE.

---

##### 4. Problemas menores

1. **`find_front` pierde cruces exactos** (`solver.py:80-81`). Usa `sign(u−ε)` y producto `< 0`: si un nodo vale exactamente ε devuelve `None`. Verificado: `find_front(linspace(0,1,11), x, 0.5) → (None, 0)`. En ALE eso dispara `extinct = True`. Los experimentos lo esquivan con fases `θ = k/16 + 1/32` (no con `θ = 0`). **Fix:** detectar el cambio de `u ≥ level`, o tratar `s == 0` explícitamente.
2. **No se verifica que el frente siga dentro de la banda** entre remallados (`simple.py:73-107`, `solver.py:187-238`); `v_est` es retrasado. En una prueba de extinción radial (ε = 0.01, BDF y Eulerian) la distancia mínima del frente al borde de la banda fue de 11.4ε, o sea que está bien. Igual vale la pena un `assert` barato, por ejemplo para frentes que aceleran más.
3. **`set_reaction` es estado global de módulo** (`model.py:40-84`). En macOS `ProcessPoolExecutor` usa *spawn*, así que una reacción instalada en el proceso padre **no** llega a los *workers*, y las lambdas no se serializan. Hay que documentarlo o pasar `f` como parámetro de los solvers.
4. **`laplacian` construye una matriz `csc`** (`solver.py:43-61`) que ALE y Eulerian usan sólo para extraer tres diagonales: es 9× más caro que calcular las diagonales con NumPy (380 µs contra 43 µs, n = 3600).
5. **`_matrix`** (`ale.py:43-57`) se llama dos veces por prueba con la misma `(x, ẋ)`, cambiando sólo la diagonal de reacción. Además usa indexado elegante `ab[1, i] -= …` con `i = arange`, cuando alcanzan rebanadas.
6. **Matriz de calor constante refactorizada en cada paso** (`simple.py:84-88`). `heat` no cambia entre remallados, pero `solve_banded` la factoriza en cada paso. Ver §5.
7. **Los tests no cubren lo que más importa del artículo** (`tests/test_combustion1d.py`):
   - no hay test de `simple.py` (el esquema de los Teoremas 5.x);
   - no hay test de las propiedades del Lema 4.1 (positividad de Y, conservación exacta de W, `max U` no creciente) ni del principio de comparación discreto (dos datos ordenados, soluciones ordenadas);
   - no se testea que se lance el error de Péclet;
   - no hay test radial de ALE ni del remallado (exactitud de la transferencia en la banda);
   - no hay test de regresión contra un JSON almacenado.
   Todos son baratos de agregar. El README dice "about 10 s"; medí 50 s de reloj y 17 s de CPU.
8. **Duplicación**: `init()` aparece en `ale.py:85-88` y en `simple.py:56-59`; `_banded` y `_matrix` hacen casi lo mismo; la matriz de esparsidad está dos veces en `limit.py`, con bucles sobre `lil_matrix`.
9. **`truncation_check.py` sólo imprime**: los números (5.x) no quedan guardados en `results/`.
10. **Docstrings desactualizadas**: `common.py:63-67` (TypeError); `paper/README` dice "zero mean" para g(θ), mientras el artículo dice "media 1–3 % de la amplitud". `README` dice "about one hour" y `reproduce_paper.sh` dice "~1.5 h".
11. **`df` por defecto** (`model.py:31-32`) toma −1 en `s < 1` y 0 en `s ≥ 1`. Es aceptable para el Jacobiano de BDF, pero conviene un comentario sobre el *kink*.
12. **Falsy defaults**: `dt = dt or …` (`simple.py:53`), `ell0 = ell0 or …` y `max_dt_remesh = … or …` hacen que un valor 0 pasado explícitamente se reemplace en silencio. Mejor usar `is None`.

---

##### 5. Rendimiento (medido)

Configuración: i5-1030NG7, un hilo; tiempos de cProfile. La máquina tenía frecuencia variable, así que los absolutos fluctúan en torno a ±30 %; los repartos porcentuales son estables.

**ALE** (caso de `tab:A3`: ε = 0.001, L = 27, T = 4, Δt = 0.01; 3671 nodos; 400 pasos, 1202 pruebas, 2404 solves): **1.85 s**, unos 4.6 ms por paso.

| función | tiempo acumulado | % |
|---|---|---|
| `laplacian` (sparse `diags` → csc, + `.diagonal()` ×3) | 0.65 s | 35 % |
| `_matrix` (ensamblado en bandas + `_derivative_coeffs`) | 0.61 s | 33 % |
| `solve_banded` (con `check_finite`) | 0.37 s | 20 % |
| `chi`, `find_front`, varios | ~0.2 s | 12 % |

O sea, sólo el 20 % es álgebra lineal.

**Eulerian (S)** (ε = 0.01, Δt = 4ε², 2500 pasos, unos 3574 nodos): 1.22 s, de los cuales **65 % es `solve_banded`** (0.12–0.15 ms por solve).

**BDF** (ε = 0.005, T = 4): 34.8 s; el reparto está en M4.

Microbenchmarks con n = 3600:

| operación | µs |
|---|---|
| `solve_banded` (por defecto) | 151 |
| `solve_banded(check_finite=False)` | 126 |
| `lapack.dgtsv` | 68 |
| `dgttrs` prefactorizado (`dgttrf` una vez) | **38** |
| `laplacian()` actual | 380 |
| tres diagonales con NumPy puro | **43** (idénticas a `allclose`) |
| `splu` del Jacobiano BDF 2n×2n | 4380 |
| dos solves en bandas (sustitución por bloques) | 321 |

**Prototipo medido** (por *monkeypatch* en la copia, sin tocar el repo): Laplaciano como tres diagonales NumPy más `dgtsv` en lugar de `solve_banded`. Los frentes salen **idénticos** (diferencia máxima 0.0).

| solver | antes (mín. de 3) | después | aceleración |
|---|---|---|---|
| ALE | 2.79 s | 1.90 s | **1.5×** |
| Eulerian | 2.74 s | 0.70 s | **3.9×** |

En ALE la aceleración es menor porque `_matrix` pasa a dominar.

**Aceleraciones adicionales concretas** (estimadas):
- ALE: ensamblar una vez por prueba, usar rebanadas, saltar la tercera evaluación cuando `|g0| < tol` y escribir el algoritmo de Thomas en Numba (fusionando ensamblado, solve y búsqueda del frente). Estimo **5–10×**, es decir menos de 0.5 ms por paso. En absoluto no importa (ALE ya corre en segundos), pero sí para barridos.
- Eulerian: `dgttrf` una vez por remallado más `dgttrs` para U, y `dgtsv` (o Thomas en Numba) para Y. Estimo **4–6×**. Esto **sí** importa: E1 (2.4·10⁶ pasos, unos 25 min) y `tab_simple`/D5 dominan el tiempo de `reproduce_paper.sh`.
- BDF: interleaving + LSODA en bandas, o BDF propio con solves por bloques; Jacobiano con patrón fijo (actualizar sólo `.data`); no reiniciar el integrador en cada remallado. Estimo **≥5×**, con la salvedad de M4.

---

##### 6. Reproducibilidad

- **Cobertura.** Verifiqué la correspondencia figura/tabla → script: `reproduce_paper.sh` regenera **todas** las figuras y tablas citadas en `main.tex` (A5, A1, `tab:A3`, D1, D2, D3, D6, B4, `fig:instability`, `tab:limit`, C1, C2, C3, `tab:speedcoef`) y los números en línea de E1–E4 y de (5.x). También genera material que no está en el artículo (A2–A4, B1–B3, D4, D5, `tab_simple`, `radial`), que conviene separar o marcar.
- **Ejecuté:**
  - `pytest`: **19 passed in 50.01 s**.
  - `speed_coefficient`: unos 27 s; JSON **bit a bit idéntico** al almacenado y a la tabla del artículo.
  - `truncation_check`: unos 29 s; reproduce −0.2351, −0.2345, −0.2342. No guarda JSON.
  - `convergence` (grupo B, 156 corridas): 596 s de reloj (unos 10 min, 1922 s de CPU) con 4 procesos. Las 156 corridas reproducen el JSON almacenado: la diferencia absoluta máxima en el error medido es 1.2·10⁻⁸, la relativa mediana 2·10⁻⁸ y la máxima 8·10⁻⁴, en un error de 4·10⁻⁶. No son bit a bit idénticas (diferencias de redondeo entre versiones de BLAS/SciPy), pero son indistinguibles a los dígitos que publica el artículo. Las figuras B1–B4 se regeneraron sin errores.
  - los tres ejemplos: OK (9 s, 24 s, 22 s).
- **Tiempo total.** La suma de `seconds` registrada en los JSON almacenados es de unos 7 650 s de CPU (unas 2.1 h). Con 4 procesos son 35–60 min en la máquina del autor; el README dice "≈1 h" y el script "≈1.5 h". En mi máquina el grupo B tardó unos 10 min (el README dice "~2 min"). Estimo **2–3 h** en un portátil de 4 núcleos.
- **Versiones.** `requirements.txt` sólo tiene cotas inferiores (`numpy>=2.0`, etc.). Por ahora funciona con las versiones más nuevas (numpy 2.5.5 / scipy 1.18.1). Pero no hay *lock file*, y los JSON no registran versión de NumPy/SciPy, commit ni máquina. **Fix:** `requirements-lock.txt` (o `environment.yml`) con las versiones exactas usadas, y guardar `{"git": sha, "numpy": …, "scipy": …, "platform": …}` en cada JSON.
- **Semillas.** No aplica: todo es determinista.
- **Artefactos.** `results/paper/C3_quenching_dt0.005.json` y `F_floor_old_code.json` no los produce ningún script actual (este último refiere a commits `d5cd9b8`/`46cef8c`). Los `.npz` de A1/A4 no están en git (se regeneran en unos 4 min); eso está bien documentado.
- **Empaquetado.** No hay `pyproject.toml`: el paquete funciona sólo con `PYTHONPATH`. Sugiero `pip install -e .`.
- **Listo para material suplementario: todavía no.**
  - **No hay LICENSE.** Sin licencia el código no es legalmente reutilizable; sugiero MIT o BSD-3.
  - Falta `CITATION.cff` y un DOI de Zenodo (release etiquetada con el commit de las figuras).
  - `experiments/paper/README.md` contiene la historia "after the referee report" (M7a, M8(c), N5…). No corresponde en material suplementario público; conviene moverla a un CHANGELOG interno o borrarla.
  - El README raíz está parcialmente en castellano y el repo incluye los PDFs de los artículos de tesis (IFB 2000, SIAM 2004). Conviene revisar los derechos de distribución de esos PDFs, o publicar sólo `numerics/` y `paper/` en el depósito de Zenodo.

---

##### 7. Recomendación

**Aceptable con revisiones menores-moderadas del código antes de publicarlo como suplemento.** El núcleo numérico es correcto y corresponde al esquema descrito; no encontré nada que invalide los resultados del artículo.

Lo imprescindible:
1. M2: no clasificar errores como quenching; re-verificar C3.
2. M1 y M3: robustez del mallador y de la malla móvil.
3. M5.4: frente inicial en BDF y Eulerian.
4. Matizar en el texto la comparación en segundos de `tab:A3` (M4), la zona rígida (M5.1) y las evaluaciones de la secante (M5.2).
5. LICENSE, versiones fijadas con metadatos en los JSON, DOI de Zenodo y limpieza de archivos huérfanos y de las notas del referato.
6. Tests del Lema 4.1, del principio de comparación y de `simple.py`.

Las optimizaciones de §5 son opcionales: con un esfuerzo chico (unas 50 líneas) reducirían el tiempo de `reproduce_paper.sh` probablemente a la mitad.
