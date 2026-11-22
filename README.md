# Combustion: tesis doctoral

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
