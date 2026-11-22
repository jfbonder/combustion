# Numérica del límite de alta energía de activación: etapa 1 (1D y radial)

Este código resuelve el problema regularizado

$$u_t=\Delta u-(u+w)\,f_\varepsilon(u),\qquad w_t=\Delta w,\qquad f_\varepsilon(s)=\varepsilon^{-2}f(s/\varepsilon),$$

donde $w=Y-u$ y $w(\cdot,0)=\varepsilon w_0$. Lo compara con el problema de frontera libre límite

$$u_t=\Delta u\ \text{ en }\{u>0\},\qquad |\nabla u|=\sqrt{2M(w_0)}\ \text{ en }\partial\{u>0\},\qquad M(w_0)=\int_{-w_0}^1(s+w_0)f(s)\,ds,$$

en geometría plana ($N=1$) y radial ($N=2,3$).

Tomamos $f(s)=(1-s)^+$, que da $M(w_0)=(1+w_0)^3/6$ en forma cerrada.

## Método

- **Espacio:** volúmenes finitos centrados en los nodos, sobre una malla no uniforme. El esquema es conservativo y en el caso radial trata $r=0$ sin casos especiales.
- **Tiempo:** BDF de `scipy`, de orden y paso variables, con jacobiano analítico disperso.
- **Malla que persigue el frente.** En cada intervalo de re-mallado:
  1. se ubica el frente $\{u=\varepsilon\}$; para $u\ge\varepsilon$ la reacción es nula;
  2. se arma una banda que cubre la zona de reacción, más un margen que el frente no alcanza a cruzar antes del próximo re-mallado;
  3. se construye una malla con tamaño $h_f=\varepsilon/16$ en la banda, celdas que crecen geométricamente (factor 1,1) hacia afuera, y tamaño $h_c$ lejos;
  4. $u$ y $w$ se transfieren por interpolación cúbica monótona (PCHIP).

  El refinamiento delante del frente y el des-refinamiento detrás ocurren en el mismo paso. El intervalo entre re-mallados se adapta a la velocidad medida del frente. La cantidad de nodos es $O(1/h_c)+O(\text{banda}/h_f)+O(\log(h_c/h_f))$, independiente del recorrido del frente: entre 500 y 900 nodos en todos los casos.
- **Datos iniciales:** $u_0^\varepsilon=u_0-\varepsilon w_0\,(1-u_0/\varepsilon)^+$ e $Y_0=u_0^\varepsilon+\varepsilon w_0$. Así, $(Y_0-u_0^\varepsilon)/\varepsilon=w_0$ exactamente, la región quemada tiene $Y_0=0$ y $u_0^\varepsilon\to u_0$ uniformemente.

## Soluciones exactas usadas para validar

1. **Ondas viajeras planas.**
   - En el límite: $u=A\,(1-e^{-c_0(x-c_0t)})^+$ con $c_0=\sqrt{2M(w_0)}/A$.
   - Para $(P_\varepsilon)$ con $w=\varepsilon w_0$: la velocidad $c_\varepsilon$ se calcula por tiro sobre la EDO de la onda, con bisección guiada por la monotonía de la energía. Es independiente del resolvedor de la EDP.
2. **Solución autosemejante radial** (Caffarelli–Vázquez; Prop. 4.3 del paper de SIAM):

   $$V(x,t)=\kappa\,(T-t)^{1/2}\,\phi\big(|x|(T-t)^{-1/2}\big),\qquad \phi(r)=M\big(-\tfrac12,\tfrac N2,\tfrac{r^2}{4}\big),$$

   donde $M(a,b,z)$ es la función de Kummer, $R$ es el primer cero de $\phi$ y $\kappa=\sqrt{2M}/|\phi'(R)|$. El soporte es la bola de radio $R\sqrt{T-t}$, que se extingue en $t=T$.

   **Errata en el paper de SIAM:** en (4.6), el término debe ser $\big(\frac{N-1}{r}-\frac r2\big)h'$, con signo **menos**. Con $+\frac r2$, la función $V$ no es calórica. El test `test_self_similar_is_exact_solution` lo verifica.

## Resultados

Las corridas completas están en curso; las tablas y figuras se agregan al terminar.

## Uso

```bash
python -m venv .venv && .venv/bin/pip install -r requirements.txt
```

```bash
.venv/bin/pytest -q
```

```bash
.venv/bin/python -m experiments.traveling_wave
```

```bash
.venv/bin/python -m experiments.radial
```

```bash
.venv/bin/python -m experiments.plots
```

Todos los comandos se corren desde `numerics/`. Los resultados quedan en `results/*.json` y las figuras en `figures/`.

## Limitaciones y próximos pasos

- El código supone **un solo frente**, que es lo que pasa en 1D y en el caso radial. Para varios frentes, o en 2D, el indicador tiene que marcar celdas en lugar de un intervalo; esa es la etapa 2.
- El costo por corrida crece como $\varepsilon^{-1}$ a $\varepsilon^{-1.5}$, por la cantidad de pasos de tiempo: el frente cruza $O(1/\varepsilon)$ celdas finas. Una malla genuinamente móvil (ALE o MMPDE), con nodos que se desplazan con el frente, permitiría pasos de tiempo independientes de $\varepsilon$ en régimen de onda viajera.
- No hay todavía un resolvedor directo del problema límite para datos generales, por ejemplo por transformación de Landau en 1D. Por ahora la comparación usa las dos soluciones exactas.
