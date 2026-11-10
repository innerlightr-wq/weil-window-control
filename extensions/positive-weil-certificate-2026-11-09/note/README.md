# Technical note — reproducibility package

**Certified Weil Positivity at `L = 4/5`: Two-Sided Bounds, Parity Separation,
and Cutoff Diagnostics**
Elias De Jesús · Independent Researcher · ORCID
[0009-0007-0190-9143](https://orcid.org/0009-0007-0190-9143) · version 1.0.0.

This directory holds the manuscript and the submission package. The computation
it reports lives one level up, in
`extensions/positive-weil-certificate-2026-11-09/`.

> **Status.** Technical note / preprint. Not peer-reviewed. **No DOI has been
> assigned to this note.** It is not a DOI-bearing publication, and the
> repository's root `CITATION.cff` — which refers to the separate sensitivity
> paper — is unchanged and remains the repository's preferred citation. Cite
> *this note* with `CITATION.cff` / `CITATION.bib` in this directory.

## The result, in one line

```
   1.41e-17  <=  lambda_Q(4/5)  <=  1.9014e-17
```

for the window Weil form at `L = 4/5` in the normalization of `../DERIVATION.md`
§1, with the odd sector separately bounded below by `2.3e-17`, so the **even**
infimum is the one that determines `lambda_Q`.

**Quote the upper bound as `1.9014e-17`.** The decimal
`1.90138291114546440e-17` that appears in earlier internal documents is the
enclosure's *midpoint* and is not a valid upper bound for its own ball; the
strict endpoint is `U = 1.901382911145464406030000e-17`. See
`CLAIM_INDEX.md` §2 and the note's Remark 8.

## Contents

| file | what it is |
|---|---|
| `weil_L45_certified_positivity.tex` | manuscript source (LaTeX, self-contained, inline bibliography) |
| `weil_L45_certified_positivity.pdf` | typeset manuscript, 14 pages |
| `CLAIM_INDEX.md` | every claim → artifact file, field, precision, enclosure, verifier |
| `consistency_check.py` | re-derives every quoted decimal from the artifacts and checks rounding direction, scope and parity labels (47 checks) |
| `zenodo.json` | deposition metadata, prepared locally — **nothing is uploaded** |
| `ZENODO_CHECKLIST.md` | manual upload checklist |
| `CITATION.cff`, `CITATION.bib` | how to cite this note |
| `make_archive.sh` | builds the reproducibility archive, deterministically |
| `weil-L45-certificate-1.0.0.tar.gz` | the archive (build product; not committed — regenerate with `./make_archive.sh`) |
| `SHA256SUMS.txt` | checksums for every file in this package and every artifact it cites |

## Reproducing

All commands are run from
`extensions/positive-weil-certificate-2026-11-09/`. The interpreter is the
repository's virtualenv.

### Fast check — verify everything from the frozen artifacts

This re-derives every certified implication without reassembling any matrix. It
is the check to run first.

```sh
V=../../.venv/bin/python

PYTHONPATH=src:../../src                    $V verify.py              # 34 checks
PYTHONPATH=phase2/src:src                   $V phase2/verify_phase2.py # 19 checks
PYTHONPATH=phase2/src:src                   $V phase2/verify_phase3.py # 20 checks
PYTHONPATH=phase2/src:src                   $V phase3/verify_audit.py  # 16 checks
PYTHONPATH=phase2/src:phase4/src            $V phase4/verify_t400.py   # 18 checks
PYTHONPATH=phase2/src:phase4/src:phase5/src $V phase5/test_prec.py     # 12 checks
PYTHONPATH=phase2/src                       $V note/consistency_check.py # 47 checks
```

Each exits non-zero on any failure. Measured wall time for the whole sequence in
a fresh process on the machine described below: **97 s**, split as

| verifier | checks | time |
|---|---|---|
| `verify.py` | 34 | 20 s |
| `phase2/verify_phase2.py` | 19 | 1 s |
| `phase2/verify_phase3.py` | 20 | 31 s |
| `phase3/verify_audit.py` | 16 | 36 s |
| `phase4/verify_t400.py` | 18 | 9 s |
| `phase5/test_prec.py` | 12 | 0 s |
| `note/consistency_check.py` | 47 | 0 s |

for **166 checks** in total, all passing.

### Full regeneration — rebuild the certificates from scratch

Hours, not minutes. The two matrix assemblies dominate.

```sh
PH="PYTHONPATH=phase2/src:src ../../.venv/bin/python"
eval $PH phase2/src/baseline.py        # gate / complement / coupling      ~45 s
eval $PH phase2/src/run_block.py       # assemble + certify at T=200       ~22 min
eval $PH phase3/src/audit.py           # beta monotonicity, K_T identity   ~2 min
eval $PH phase3/src/minimizer.py       # freeze f_200, r_f and e_f         ~7 min

PH4="PYTHONPATH=phase2/src:phase4/src ../../.venv/bin/python"
eval $PH4 phase4/src/cutoff.py         # cutoff identity, beta_400         ~1 min
eval $PH4 phase4/src/eval_f200.py      # compact-route Q(f_200)            ~7 min
eval $PH4 phase4/src/run_block400.py   # assemble + certify at T=400       ~60 min
eval $PH4 phase4/src/cross.py          # the cross table, Q(f_400)         ~26 min

PH5="PYTHONPATH=phase2/src:phase4/src:phase5/src ../../.venv/bin/python"
eval $PH5 phase4/src/emit_certificate_T400.py  # rewrite certificate_T400  ~1 min
eval $PH5 phase5/src/conventions.py    # degree-2 convention checks        ~1 min
eval $PH5 phase5/src/parity.py         # replay + parity separation        ~2 min
eval $PH5 phase5/src/pilot800.py       # R_200(f_400) + T=800 scalar pilot ~3.5 h
```

Total for a cold rebuild is roughly **six hours**, of which `pilot800.py` is
about `3.5 h` and the two assemblies about `1.4 h`. These are measured figures
from the run that produced the artifacts, not estimates; see `CLAIM_INDEX.md`
§6.

### Building the PDF

```sh
cd note && pdflatex -interaction=nonstopmode weil_L45_certified_positivity.tex
         && pdflatex -interaction=nonstopmode weil_L45_certified_positivity.tex
```

No bibliography pass is needed: the reference list is an inline
`thebibliography`, following the convention of the repository's existing paper.

## Pinned environment

| component | version |
|---|---|
| OS | Linux 7.2.5 (x86_64) |
| Python | 3.14.7 |
| `python-flint` | 0.9.0 (FLINT / Arb) |
| `mpmath` | 1.3.0 |
| `numpy`, `gmpy2`, `sympy` | **absent** — not used |

`mpmath` is used only to *locate* candidate test vectors, which are then frozen
as exact rationals and certified in Arb. It appears in no pass/fail decision.
Arb working precisions on sign-critical paths are `256`–`2000` bits and are
recorded per artifact in `CLAIM_INDEX.md`.

## What is in the archive

`weil-L45-certificate-v1.0.0.tar.gz` (built by `make_archive.sh`) contains only
redistributable files authored in this repository: the manuscript source and
PDF, this package, and the extension's documentation, source and certificates.
It contains **no third-party material** — in particular no copy of
arXiv:2608.24827, whose reproducibility materials were never obtained and whose
matrix step is not reproduced here (`../SOURCE_AUDIT.md` §3).

Licensing follows the repository split: code (`*.py`) under **MIT**
(`../../LICENSE`), the manuscript and documentation under **CC BY 4.0**
(`../../LICENSE-CC-BY-4.0`).

## Scope, stated once

This is a computation at a single window. It establishes no arithmetic estimate,
no statement about the location of any zero of `zeta`, no verified range of
zeros, no result at any larger support, and nothing about the Riemann
hypothesis. The finite-reduction strategy is from arXiv:2608.24827v2 and is
credited there; parity decoupling is likewise not new. What is offered here is a
certified computation inside that framework, with its artifacts and its errors
both published.
