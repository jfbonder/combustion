# The high activation energy limit of an equidiffusional flame model

Source of the article

> J. Fernández Bonder, *The high activation energy limit of an equidiffusional flame model:
> rate of convergence and a convergent front-following scheme* (2026, draft),

and the code that produces all its numerical results.

The article studies the limit of large activation energy in a model of premixed flames
with equal diffusivities, whose limit is a free boundary problem for the heat equation.
It proves a rate of convergence for this singular limit and the convergence, with an
explicit error bound, of a fully discrete scheme whose mesh has a fine band that follows
the flame front.

## Contents

| Path | What it is |
|---|---|
| [`paper/main.pdf`](paper/main.pdf), `paper/main.tex` | The article. |
| [`paper/supplement.pdf`](paper/supplement.pdf), `paper/supplement.tex` | Supplementary material: experiments on travelling waves and on the cost of the method. |
| `paper/references.bib` | Bibliography. |
| [`numerics/`](numerics/README.md) | The Python package `combustion1d`, examples, tests, and the scripts, data and figures of every numerical experiment of the article. |

## The code

`numerics/` contains a small Python package (NumPy, SciPy, Matplotlib) with

- the monotone, linearly implicit scheme of the article on a moving band mesh (ALE) and
  on meshes that do not move;
- a reference solver (method of lines with BDF);
- the exact travelling waves and self-similar solutions;
- a front-fixing solver for the limit free boundary problem.

Quick start:

```bash
cd numerics
python -m pip install -r requirements.txt
python -m pytest -q
python -m examples.ex1_travelling_wave
```

Installation, usage, the description of each solver and the instructions to reproduce
each figure and table are in [`numerics/README.md`](numerics/README.md); the experiments
are documented one by one in
[`numerics/experiments/paper/README.md`](numerics/experiments/paper/README.md).

## Building the article

The figures are read from `numerics/figures/`, so the article must be compiled inside
`paper/`. With a TeX distribution that provides `latexmk`:

```bash
cd paper
latexmk -pdf main.tex
latexmk -pdf supplement.tex
```

Compile `main.tex` first: the supplement takes the numbers of the sections, theorems
and equations of the article from `main.aux`.

## Reproducing the numerical results

```bash
cd numerics
./reproduce_paper.sh
```

reruns every experiment and redraws every figure (about two hours on four cores).
`./reproduce_paper.sh plots` only redraws the figures from the stored results in
`numerics/results/`. See [`numerics/README.md`](numerics/README.md) for the individual
experiments.

## Use of AI tools

The code, the numerical experiments, the figures and the tables were produced with
Claude Code (Anthropic), which also assisted in the proofs and in the writing of the
article; see the statement at the end of the article.

## License

The code in `numerics/` is distributed under the MIT License (see [`LICENSE`](LICENSE)).
The article and its supplement in `paper/` are © the author and are not covered by that
license.

## Citing

Until the article is published, please cite it as a preprint, with the address of this
repository; see [`CITATION.cff`](CITATION.cff).
