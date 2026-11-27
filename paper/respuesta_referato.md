# Respuesta al referato interno

Responde a `paper/referato_editorial.md` (25 de noviembre de 2026). El referato evaluó `main.tex` en el commit `52c3129`; esta versión parte del borrador posterior (`78bd186`), que ya contenía el teorema de convergencia del esquema al problema límite, y agrega lo que se detalla abajo.

**Decisión previa.** Se siguió el Camino A del editor: un solo artículo, unificado por el teorema de convergencia del esquema totalmente discreto al problema de frontera libre (Teorema 5.1), que usa las barreras de la Sección 3 evaluadas en la malla. Ese teorema es posterior al referato y **ningún referí lo revisó todavía**: es lo primero que debería mirar la nueva evaluación.

**Estado general.** El artículo compila en 38 páginas (antes 31). Las pruebas de los Lemas 3.1 y 3.3 pasaron a un apéndice. Todos los experimentos se volvieron a correr con el código corregido (92 minutos en 4 núcleos, más 35 de los dos grupos nuevos).

Convención de estados: **Hecho**, **Parcial** (se indica qué falta), **No hecho** (con el motivo), **Decisión del autor** (no corresponde que la tome yo).

---

## 1. Problemas mayores

| Ítem | Estado | Qué se hizo |
|---|---|---|
| **M1** Optimalidad del orden ε | Hecho | La introducción y la Obs. 3.7(3) dicen ahora que el orden ε no se mejora en norma del supremo cuando w₀ ≢ 0, y que esa diferencia es el offset de temperatura del gas quemado. Para la posición del frente se agrega un argumento que no depende del nivel elegido: por la Prop. 2.2, c_ε − c₀ = βε + O(ε²), así que cualquier conjunto de nivel se separa de la frontera libre en βεt. Vale para ondas viajeras, que quedan fuera del teorema, y así se dice. |
| **M2** El log(1/ε) | No hecho | Analicé la ruta que propone R2 y no cierra directamente. Hacen falta dos estimaciones en la capa inicial que no pude probar: una cota inferior para la cola de Y^ε con el exponente exacto (la tasa de decaimiento de la barrera difiere de la natural solo en un factor 1 + o(1), así que una cota gruesa no alcanza) y una cota inferior para u^ε en la región fresca que pierda solo O(ε). La Obs. 3.7(2) explica ahora con precisión de dónde viene el logaritmo y qué faltaría para sacarlo. |
| **M3** Novedad frente a trabajos previos | Parcial | Párrafo nuevo "Relation with previous constructions" al comienzo de la Sección 3: qué hace FBW 2004 (resolver (P_ε) en una capa; cercanía O(ε) inmediata; umbral ε₀(δ) por blow-up, sin cuantificar), qué es nuevo (umbral lineal en δ con error independiente de δ) y las tres diferencias con las barreras tipo Allen–Cahn, incluida la observación de que la barrera es "gruesa". **Falta:** no pude leer Petrosyan 2002, así que el texto solo dice que sus resultados son cualitativos y que su construcción no se extiende a w ≠ 0 (esto último lo afirma FBW 2004). La pregunta de R2, si esa construcción ya da ε₀(δ) ≥ cδ para w ≡ 0, sigue sin respuesta. |
| **M4** Clase de soluciones | Hecho | La introducción, la definición de solución clásica regular, la Obs. 3.7(4) y las conclusiones dicen que el teorema cubre bolsas acotadas de gas fresco, que la regularidad hasta t = 0 es una condición de compatibilidad, que para w₀ variable no hay resultado de existencia conocido y cuál es la clase a la que se sabe que se aplica. Se corrigió la frase "in agreement with Theorem" para el caso plano. |
| **M5** Atribución de hipótesis | Hecho | (F₂) solo se usa para derivar los perfiles respecto de W; se dice en la introducción y en la Sección 3 que con W constante alcanza (F). En la Prop. 3.4 se corrigió "nonincreasing and Lipschitz" por "Lipschitz" y se indica dónde entra la monotonía. "Arrhenius-type rates" pasó a "Arrhenius rates with a cut-off", con la aclaración de que e^{−s} no cumple (F). |
| **M6** Teorema sobre el esquema implementado | No hecho | Siguen abiertos la cota a priori de la velocidad implícita y el remallado; la Sección 5.4 lo dice explícitamente. El argumento "casi gratis" para el remallado no lo es: la interpolación de la barrera en la cola quemada no tiene holgura relativa (la barrera varía en un factor e^{h_c/ε} dentro de una celda gruesa), y cerrarlo exige rehacer parte de la prueba. |
| **M7** Las cotas no dicen nada sobre los experimentos | Hecho | Experimento nuevo con una tasa C^{1,1}(ℝ), f(s) = 3((1−s)⁺)², dentro de las hipótesis del teorema de ondas viajeras: el error es aproximadamente proporcional a ρ_h = h²/ε³ y casi no depende de Δt. La constante de la cota sigue lejos de ser ajustada, y el texto lo dice. |
| **M8** "Δt no aparece en el error" | Hecho | Reformulado como pide R3: es la invariancia de los estados estacionarios bajo un esquema linealmente implícito, vale para la onda viajera en mallas que se trasladan rígidamente, y no es una afirmación sobre el método en general. La frase sobre "Δt e^{−cR}" ya no está. |
| **M9** Sección 5.5 heurística | Hecho | La expansión formal queda rotulada como no justificada, con la obstrucción espectral que señala R3 (ψ* no decae, el espectro esencial toca el origen) y su relación con la relajación lenta. |
| **M10** Originalidad | Parcial | Párrafo de la introducción reescrito: mallas móviles para llamas (Dwyer–Sanders, Furzeland–Verwer–Zegeling), velocidad de llama como incógnita (Smooke–Miller–Kee), capas interiores móviles (Gracia–O'Riordan), Nochetto–Verdi 1996/97, Nochetto–Paolini–Verdi (EJAM 1994), Fife–McLeod para las barreras de traslación y Patankar para la reacción linealmente implícita. Se dice que lo nuevo no es seguir el frente sino combinar un esquema monótono para todo Δt con una banda móvil y analizarlos cuando ε → 0. **Falta:** ver la sección 5 de este documento sobre el grado de verificación de cada referencia. |
| **M11** Comparación de costo | Parcial | La tabla de costo se reemplazó por una subsección "Work and accuracy" con un diagrama trabajo–precisión y una tabla de pasos y resoluciones lineales, para ALE, (S) euleriano en banda y en malla uniforme, BDF y un método de líneas con Jacobiano en banda (LSODA). **Falta** la comparación con MMPDE, que no implementé; el texto lo dice. |
| **M12** Condición de Péclet | Hecho | Forma radial explícita; se aclara que "para todo Δt y toda malla" está sujeto a ella, qué hace el código cuando se viola y que no se usa upwinding. Se corrigió además que la condición plana es unilateral. |
| **M13** Interpretación física de w₀ | Hecho | La introducción declara el modelo, explica que w es el exceso de entalpía, que T_f es la temperatura adiabática solo si w ≡ 0, que w₀ es la perturbación de la temperatura de llama en el escalado de Matkowsky–Sivashinsky, y que w₀ > −1 es un efecto del corte de la cinética. "Unstable" se reemplazó en todo el texto. Se citan ZFK 1938 y Matkowsky–Sivashinsky 1979. |
| **M14** Término −w₀ en la Prop. 2.2 | Hecho | La proposición enuncia la descomposición 1 − I/√(2M) = −w₀ + ∫(1 − √(B/M)) y el texto interpreta cada término; el cambio de signo en w₀* se presenta como efecto de normalización. |
| **M15** "Quenching" | Hecho (por eliminación) | El experimento ya no figura en el artículo. En el repositorio queda con el criterio corregido; ver la sección 3. |
| **M16** Terminología y parámetros | Hecho | "Extinction" se define al introducir la solución autosemejante (se extingue la bolsa fresca, no la llama). Se corrigió "accelerates without bound". Párrafo nuevo sobre el rango de parámetros (números de Zeldovich, tamaño de la banda y de h_c frente al espesor térmico). ε = 0.1 y 0.05 en la Tabla 1 y en un experimento con ALE. |

## 2. Problemas menores

**Hechos**

- Prop. 2.2 con resto O(ε²): la prueba usa ahora J(φ) ≤ ½α₀t², como indica R2, y se eliminó el comentario sobre el logaritmo. Los datos lo confirman: el resto dividido por ε² vale 0.10, 0.12 y 0.65 para los tres valores de w₀.
- Lema 2.1: intervalos uniformizados, Φ'' acotada para todo z, Φ ∈ C^{2,1} en el enunciado.
- Marco de comparación: dominios no acotados (Phragmén–Lindelöf), datos continuos, regla de la cadena a.e., D acotado en la Prop. 3.5.
- Lema 3.1(ii): la dependencia C² en W se prueba ahora con el principio de contracción uniforme en lugar de invocar la variedad inestable, y se escribe la identidad para P_WW.
- Lema 3.3: cota corregida a −(3/2)K₀d. **Esto obligó a cambiar r₁ de m/(2K₀) a m/(3K₀)**, cosa que el referí no había notado; revisé que las cotas que usan r₁ siguen valiendo.
- Constantes: Λ, K, c₁ en los enunciados de las Props. 3.4–3.5; orden de elección explicitado; observación nueva sobre de qué depende C y que no crece exponencialmente en T.
- (W): corregida la nota sobre w₀ Lipschitz y la condición η_ε ≤ min{δ₁/2, 1}. La condición de monotonía se puede pedir en t = 0.
- Lema 4.2 en mallas infinitas: existencia en ℓ^∞ por serie de Neumann.
- Sección 5: carácter a posteriori del teorema de ondas viajeras, barreras que no cumplen flujo nulo en los extremos, ρ_h = h²/ε³ excluye f = (1−s)⁺.
- Sección 4: forma no conservativa del ALE; descripción de la secante, del chequeo de malla y del remallado euleriano según lo que hace el código.
- Notación: la hipótesis de estrellado es (H); la energía del Lema 2.1, las matrices A_F y el error acumulado M_n cambiaron de símbolo para no chocar con E, A y M.
- La figura de dependencia en w₀ se regeneró con el código actual.
- Por qué no se varía ϑ en cada remallado: explicado en el texto.

**No hechos**

- Tabla de dependencia de las constantes (prioridad B): hay una observación en prosa, no una tabla.
- Tablas de órdenes observados para datos generales (prioridad B): sigue solo la figura de extinción radial.
- Sección 5.5 como teorema en espacios con peso (prioridad C).
- Bola de llama 3D, cinética Arrhenius truncada, ignición por núcleo (prioridad C / opcional).

## 3. Código

| Ítem | Estado | Detalle |
|---|---|---|
| C1–C2 Fallas | Hecho | Excepciones propias (`NumericalFailure`, `PecletError`, `MeshError`), estado `ok` / `extinct` / `failed` en los resultados. La bisección de quenching cuenta solo extinciones reales y aborta ante cualquier otra falla. |
| C3 Celdas diminutas | Hecho | El mallador fusiona el resto con la celda vecina; test con anclas aleatorias. |
| C4 Validez de la malla | Hecho | Se remalla si una celda cae por debajo de h_f/2; un ensayo con nodos cruzados se descarta. |
| C5 Frente inicial | Hecho | Corregido en BDF y en el euleriano. |
| C6–C9 | Hecho | `find_front` con un nodo exactamente en el nivel; reacción por llamada y en los procesos hijos; defaults con `is None`; docstrings y JSON huérfanos eliminados. |
| Tests | Hecho | 39 tests (antes 19), 14 s: cotas del Lema 4.1, principio de comparación, error de Péclet, ALE radial, transferencia exacta en la banda, regresión contra un JSON. |
| Procedencia | Hecho | Los JSON nuevos guardan commit y versiones de las librerías. |
| README de experimentos | Hecho | Reescrito como documentación neutral. |
| Rendimiento | No hecho | Fuera de alcance (prioridad C). |

**Hallazgos de la auditoría que el referato no tenía**

- El código anterior reproducía bit a bit los resultados guardados en los 95 casos baratos que se probaron. Con el código nuevo cambian 9, todos de BDF o eulerianos, en los últimos dígitos.
- El temor de R4, que el radio crítico del experimento de quenching se hubiera decidido por una excepción, **no se confirmó**. El defecto real era otro: para ε = 0.00125 todos los "quenching" guardados se contaron con el frente todavía presente, porque se salía de la banda fina y quedaba frenado. Con un remallado cuando el frente sale de la banda esos núcleos se extinguen de verdad, y en la corrida completa los seis intervalos guardados para el radio crítico se reproducen sin cambios.
- "In all our runs it takes three evaluations" era falso para la figura de extinción radial: hasta 4.2 evaluaciones por paso en promedio con Δt grande. Corregido en el texto.
- Las cifras de R4 sobre compresión de celdas (0.80 / 0.76 / 0.50 h_f) no se reprodujeron: en las 30 corridas de esa figura ninguna celda bajó de h_f/2.

## 4. Experimentos nuevos y números que cambiaron

- **w₀ no constante** (subsección nueva, prioridad A): un frente plano en una mezcla estratificada y una bolsa radial con N = 3. Los errores contra el problema límite calculado son proporcionales a ε sin logaritmo visible: front*/ε ≈ 0.73 (plano) y 0.40 (radial). El caso radial cumple la condición de monotonía del teorema de tasa.
- **Trabajo–precisión** (subsección nueva): para la onda viajera el error del ALE lo fija la malla de la banda (4.5·10⁻⁴ con h_f = ε/16, 10⁻⁴ con ε/32, 10⁻⁵ con ε/64), no depende de Δt ni de ε, y cuesta décimas de segundo. BDF y LSODA llegan a errores entre 10⁻⁵ y 6·10⁻⁵ en 5–18 s con ε = 0.01 y en 43–190 s con ε = 0.001. El texto aclara que la onda viajera es el caso más favorable para el ALE.
- **ε = 0.1 y 0.05 con ALE**: mismos errores que con ε chico.
- **Tasa regular**: ver M7.
- **Tabla de comparación con el problema límite**: cambiaron las dos filas planas de ε = 0.002 (0.321 / 0.333) y ε = 0.001 (0.296 / 0.302). El resto de los números del texto se verificó contra las corridas nuevas y no cambió.

## 5. Bibliografía: qué está verificado y qué no

La verificación se hizo solo con resultados de búsqueda: no se pudo abrir ninguna página de editorial ni resolver ningún DOI desde esta sesión. Por eso las afirmaciones sobre el contenido de los trabajos citados son deliberadamente neutras.

- **Datos bibliográficos confirmados:** Nochetto–Paolini–Verdi (EJAM 1994), Fife–McLeod 1977, Gracia–O'Riordan 2012, Dwyer–Sanders 1978, Matkowsky–Sivashinsky 1979, Langlois–Marion 2009.
- **Confirmados con algún campo dudoso:** Nochetto–Verdi 1996 y 1997 (exponentes de los acoplamientos no vistos en la fuente), Furzeland–Verwer–Zegeling 1990, Smooke–Miller–Kee 1983, Patankar 1980, Chow–Hale 1982 (la ubicación del principio de contracción uniforme, Sección 2.2, viene de fuentes secundarias).
- **Correcciones a lo que proponía el referato:**
  - El artículo de EJAM es sobre un modelo de relajación de fase, no sobre el doble obstáculo.
  - Larrouturou 1991 trata la positividad bajo la discretización convectiva, no términos de reacción; no se cita.
  - Que Smooke–Miller–Kee fijen la temperatura en un punto no se pudo confirmar en ese artículo; el texto solo dice que calculan la velocidad de llama como incógnita.
- **Sin poder determinar:** si Glangetas 1992, Glangetas–Roquejoffre 1996 o Langlois–Marion 2009 contienen alguna tasa. Glangetas 1992 y Langlois–Marion se citan de forma neutra; conviene leerlos antes de sostener "primera tasa".

## 6. Decisiones que quedan para el autor

- **Licencia** del código (`LICENSE`), `CITATION.cff` y DOI de Zenodo.
- **Derechos de distribución** de `ifbj.pdf` y `final.pdf`, que están en un repositorio público.
- **Largo:** 38 páginas. No verifiqué el límite de SINUM. Si hay que acortar, los candidatos son la subsección de dependencia en w₀, el detalle del pinning y parte de la figura de transitorios.
- **Petrosyan 2002** y las referencias de la sección 5: leerlas y ajustar las dos frases que dependen de ellas.
- **Experimento de quenching:** ya no está en el artículo; decidir si se conserva en el repositorio.

## 7. Reorientación posterior: *Interfaces and Free Boundaries*

Después de escrita esta respuesta se decidió apuntar a *Interfaces and Free Boundaries* en lugar de SINUM, porque la parte numérica es esencialmente unidimensional. Cambios respecto de lo descripto arriba:

- **Eje del artículo:** el resultado principal es la tasa del límite singular; el segundo, la convergencia de un esquema totalmente discreto al problema de frontera libre. Título, abstract e introducción se reordenaron en ese sentido.
- **Pruebas:** las de los Lemas 3.1 y 3.3 volvieron al cuerpo del artículo (ya no hay apéndice).
- **Sección 6:** se reordenó (comparación con el problema límite, w₀ no constante, dependencia en w₀, marco del teorema de convergencia, ondas viajeras, perfil autosemejante). El detalle de los experimentos con ondas viajeras (pinning, tamaños de malla, arrastre, retraso en mallas fijas) y la comparación de trabajo–precisión pasaron a `paper/supplement.tex`; en el artículo queda un resumen.
- **Largo:** 34 páginas el artículo y 5 el suplemento.

Con este destino pesan más los puntos de EDP que quedaron abiertos: M2 (el logaritmo), M3 (Petrosyan 2002) y la existencia de soluciones clásicas para w₀ variable.

