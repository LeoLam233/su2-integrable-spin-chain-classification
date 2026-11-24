# All-spin classification of isotropic nearest-neighbour integrable chains

This repository contains the manuscript **All-spin classification of isotropic nearest-neighbour integrable chains** (24 November 2026), its LaTeX source, and supplementary exact checks. The manuscript is the expanded reading edition (13 pages), with a proof roadmap, detailed orbit-coefficient derivation, and explicit construction of the singlet test vectors.

The manuscript classifies Hermitian SU(2)-invariant nearest-neighbour Hamiltonian densities on irreducible spin-S sites satisfying the Reshetikhin condition, up to affine changes. It derives the four classical families and the exceptional spin-three family using angular-momentum recoupling, rank-one and reflection identities, and three-site permutation symmetry. A corollary uses the local-conservation-law dichotomy of Shiraishi and Sanatani.

## Manuscript

- [PDF](classification.pdf)
- [LaTeX source](classification.tex)

To compile the manuscript, run `pdflatex -halt-on-error classification.tex` twice.

## Supplementary verification

- [`verify_reconstruction.py`](verify_reconstruction.py): exact symbolic identities and finite-sector regression checks.
- [`verification.log`](verification.log): recorded output (SymPy 1.14.0).

```sh
python -m pip install sympy==1.14.0
python verify_reconstruction.py
```

The general argument is in the manuscript; the finite checks are supplementary and do not establish the classification by extrapolation.

The manuscript is a research preprint. The current version has not undergone independent external peer review.
