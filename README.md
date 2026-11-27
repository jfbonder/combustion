# Combustion

This repository contains:

- **`paper/`** — the article *The high activation energy limit of an equidiffusional flame
  model: rate of convergence and a front-following numerical method* (J. Fernández Bonder,
  2026): LaTeX source `main.tex`, bibliography and the compiled draft.
- **`numerics/`** — the Python code of the article (package `combustion1d`), the scripts
  that reproduce every figure and table, and examples. See
  [`numerics/README.md`](numerics/README.md) for installation, usage and reproduction.
- The two thesis papers (below), on which the article builds.

`paper/referato_editorial.md` is an internal adversarial referee report on an earlier
draft, and `paper/respuesta_referato.md` the point-by-point response (both in Spanish).
`paper/supplement.tex` is the supplementary material (experiments on travelling waves and
on the cost of the method); compile `main.tex` first, since it uses its labels.

## Tesis doctoral

Papers de la tesis doctoral de Julián Fernández Bonder, en colaboración con Noemí Wolanski.

| Archivo | Referencia |
|---|---|
| `ifbj.tex`, `ifbj.pdf` | J. Fernández Bonder and N. Wolanski. *A free boundary problem in combustion theory*. Interfaces Free Bound. **2** (2000), 381–411. |
| `final.tex`, `final.pdf` | J. Fernández Bonder and N. Wolanski. *Uniqueness of limit solutions to a free boundary problem from combustion*. SIAM J. Math. Anal. **36** (2004), no. 1, 172–185. |

## Compilación

- `final.tex` usa la clase `siamltex.cls` de SIAM, que no viene con TeX Live y hay que tenerla
  aparte. Necesita `amsmath`/`amssymb`, que ahora se cargan en el preámbulo.
- `ifbj.tex` está escrito en LaTeX 2.09 (`\documentstyle{amsart}`) y compila en modo de
  compatibilidad. Los dos avisos sobre `\negmedspace`/`\negthickspace` son inocuos.

Los PDFs del repositorio se generaron en noviembre de 2026 a partir de estos fuentes:
`ifbj.pdf` con `amsart` y `final.pdf` con la clase original de SIAM (`siamltex.cls` y
`siam10.clo`, versión 1.0 de 1995). Esos dos archivos no están en el repositorio, porque su
licencia no permite distribuirlos por separado del paquete completo de SIAM.

Los comentarios `% NOTE (2026 review)` en los fuentes marcan puntos que conviene revisar.
No modifican el texto compilado.
