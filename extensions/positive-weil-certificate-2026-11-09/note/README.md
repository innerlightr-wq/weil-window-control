# Technical note — reproducibility package

**Certified Weil Positivity in Two Windows: Spectral Bounds, Conditioning, and
Reproducible Certificates**
Elias De Jesús · Independent Researcher · ORCID
[0009-0007-0190-9143](https://orcid.org/0009-0007-0190-9143) · **version 2.0.0**

> **Status.** Version 2.0.0 is the **first archived release** of this note.
> Versions 1.0.x–1.2.0 were internal drafts and were **never deposited**.
> Technical note / preprint; **not peer-reviewed**.
> **No DOI has been assigned** at the time of writing — there is deliberately no
> DOI field in `zenodo.json`, `CITATION.cff` or `CITATION.bib`. The repository's
> root `CITATION.cff` refers to the separate sensitivity paper and is unchanged.

This directory holds the manuscript and the submission package. The computation
it reports lives one level up, in
`extensions/positive-weil-certificate-2026-11-09/`, with the two-window and
conditioning work in its `phase7/` subdirectory.

## The results, in one place

```
    1.43300048827856e-17  <=  lambda_Q(4/5)  <=  1.674e-17
    8.27343868254787e-18  <=  lambda_Q(L~)   <=  9.475e-18
                                      L~ = 80471895621705018 / 10^17  <  (log 5)/2
```

Monotonicity of `lambda_Q` in `L` extends the first lower bound to all
`0 < L <= 4/5`, the second to all `0 < L <= L~`, and the second **upper** bound
to every `L >= L~`, including exactly `(log 5)/2`. **No positive lower bound at
`(log 5)/2` follows** — monotonicity runs the other way.

Pairing the four endpoints bounds the change in the true minima:

```
    0.494  <=  lambda_Q(L~) / lambda_Q(4/5)  <=  0.662
```

a fall of between `33.8 %` and `50.6 %`. At `L = 4/5`, per parity sector:

| conditions | certified constant |
|---|---|
| 4 explicit even directions | `Q >= 0.617553710936882 ||f||^2` |
| 5 explicit even directions | `Q >= 0.853492736815556 ||f||^2` |
| 3 explicit odd directions | `Q >= 0.106805801391493 ||f||^2` |
| 4 explicit odd directions | `Q >= 0.784130096434765 ||f||^2` |

and **exactly four** even eigenvalues of `Q` lie below `10^-2` (at most four from
the four-condition certificate, at least four from a certified four-dimensional
test space). In the odd sector only **at most three** is certified. The parity
separation at `L = 4/5` is `lambda_Q,odd >= 1.38569999999861e-14`.

Short decimals are **directed**: an upper bound is rounded up, a lower bound
rounded down. `1.674e-17` is a valid short form for the `L = 4/5` upper
endpoint and `1.6731e-17` is **not**; `0.6175` is valid for the four-condition
even constant and `0.6176` is **not**. The release gate checks this.

## Excluded from this release

An earlier compact time-domain evaluator (`phase4/src/compact.py`) omitted the
truncation error of its outer panelled Gauss–Legendre rule. Ball arithmetic
cannot supply that term: it encloses nodes, weights and rounding, not the
difference between a quadrature sum and its integral. That route is **excluded**
from every claim above, and the quantities it produced — `q(f_200)`, `q(f_400)`,
the endpoint once quoted from the latter, and the discarded energies
`e_T = q - r_T` with their ratios — remain **uncertified**, retained only as
labelled measurements. See `../RETRACTIONS.md` R14. The open outer remainder is
off the critical path: closing it would restore a second route to a bracket this
release already has.

## Reproducing

### The release gate — run this first

```sh
cd extensions/positive-weil-certificate-2026-11-09
python phase7/acceptance.py
```

**61 checks, 0 failed, about 210 s.** Nonzero exit on any failure, so it is
usable as a release gate. It verifies the **integrated** candidate from the
packaged files: it applies and reverts no patch, reads no value held in memory
from the run that produced it, and can be run twice. Last run from a clean
extraction of the archive with `HOME=/nonexistent`, so it reads no
home-directory path and no hidden cache.

Its stages: the repository suite below; the closed-form upper bounds
re-certified from their stored exact rationals; a Rump replay of **both**
improved unconditioned `L = 4/5` shifts (`phase4/verify_t400.py` does *not*
cover these — it tests the earlier committed shifts, and the strict Schur test
assumes `A >= aI` rather than checking it); all eight Schur checks from their own
stored `a, d, kappa, mu`, with rejection required at `mu = a`; the four
conditioned rank-`K` bumps rebuilt from the stored rational vectors; both
second-window blocks; the count certificate with its Gram matrix; the
directed-decimal checks and the two-window ratio; and the complement, coupling
and Bernstein-remainder provenance fields.

### The repository suite, on its own

```sh
EXT=extensions/positive-weil-certificate-2026-11-09
PP=$EXT/phase2/src:$EXT/phase4/src:$EXT/phase5/src:$EXT/src:src
cd $EXT
PYTHONPATH=$PP python verify.py                      # 34 checks
PYTHONPATH=$PP python phase2/verify_phase2.py        # 19 checks
PYTHONPATH=$PP python phase2/verify_phase3.py        # 20 checks
PYTHONPATH=$PP python phase3/verify_audit.py         # 16 checks
PYTHONPATH=$PP python phase4/verify_t400.py          # 18 checks
PYTHONPATH=$PP python phase5/test_prec.py            # 17 checks
PYTHONPATH=$PP python note/consistency_check.py      # 52 checks
PYTHONPATH=$PP python phase5/src/upper_highN.py --verify
```

**176 checks**, all passing, every exit code 0.

### What is replay and what is a cold rebuild

All of the above is **frozen-artifact replay**. The three matrix assemblies are
**not** re-run; they are cold rebuilds, measured at `1343 s` (`T = 200`,
`L = 4/5`), `3532 s` (`T = 400`, `L = 4/5`) and `3678 s` (`T = 400`, `L~`), and
are regenerated by `phase2/src/run_block.py`, `phase4/src/run_block400.py` and
`phase7/src/run_block_lam5.py`. Checksums establish that the packaged files are
the ones produced; they are **not** evidence of numerical correctness, which is
what the replay is for. Full detail in `../phase7/REPRODUCE.md`.

The verifiers are **not** independent implementations: they share the assembly
and constants modules with the computation they check, and no independence is
claimed for them.

### Building the PDF

```sh
cd note
export SOURCE_DATE_EPOCH=1794528000 FORCE_SOURCE_DATE=1   # reproducible build
pdflatex -interaction=nonstopmode weil_two_window_positivity.tex
pdflatex -interaction=nonstopmode weil_two_window_positivity.tex
```

No bibliography pass is needed: the reference list is an inline
`thebibliography`. The committed PDF is 25 pages and builds with no warnings and
no overfull boxes. With `SOURCE_DATE_EPOCH` set as above the build is
**byte-reproducible**, which in turn makes `make_archive.sh` byte-reproducible,
so the archive checksum in `SHA256SUMS.txt` can be re-derived from a clean
checkout. Without it, `pdflatex` embeds a wall-clock `CreationDate` and both
hashes change.

### Building the archive

```sh
cd note && ./make_archive.sh      # -> weil-two-window-certificate-2.0.0.tar.gz
```

Deterministic: sorted entries, fixed mtime, numeric owner `0:0`, gzip without a
timestamp, so two builds of an unchanged tree are byte-identical. Entries carry
their **repository-relative** paths, so extracting into an empty directory
reproduces the layout the gate expects, including the one repo-level module the
code imports (`src/zeta_window.py`) and both licence files.

## What is in this directory

| file | role |
|---|---|
| `weil_two_window_positivity.tex` / `.pdf` | the manuscript (25 pp.) |
| `CLAIM_INDEX.md` | every claim → artifact file, field, precision, enclosure, verifier |
| `consistency_check.py` | re-derives every quoted decimal from the artifacts and checks rounding direction, scope and parity labels (52 checks) |
| `zenodo.json` | deposition metadata, prepared locally — **nothing is uploaded, no DOI** |
| `ZENODO_METADATA.txt` | the same metadata as plain text, for pasting into the upload form |
| `ZENODO_CHECKLIST.md` | manual upload checklist |
| `CITATION.cff`, `CITATION.bib` | how to cite this note (no DOI field yet) |
| `make_archive.sh` | builds the reproducibility archive, deterministically |
| `weil-two-window-certificate-2.0.0.tar.gz` | the archive (build product) |
| `SHA256SUMS.txt` | checksums for every file in this package and every artifact it cites |

## Pinned environment

| component | version |
|---|---|
| OS | Linux 7.2.5 (x86_64) |
| Python | 3.14.7 |
| `python-flint` | 0.9.0 (FLINT / Arb) |
| `mpmath` | 1.3.0 |
| `numpy`, `gmpy2`, `sympy` | **absent** — not used |

`mpmath` is used only to *locate* candidate test vectors and conditioning
subspaces, which are then frozen as exact rationals and certified in Arb. It
appears in no pass/fail decision of the gate. Arb working precisions on
sign-critical paths are `256`–`2920` bits and are recorded per artifact in
`CLAIM_INDEX.md`.

## Licensing

Code (`*.py`, `*.sh`) under **MIT** (`../../LICENSE`); the manuscript and
documentation under **CC BY 4.0** (`../../LICENSE-CC-BY-4.0`). Both licence
files are included in the archive. No third-party material is redistributed —
in particular no copy of arXiv:2608.24827 or arXiv:2611.13129, which are cited.

## Scope, stated once

This is a computation at two windows. It establishes no arithmetic estimate, no
zero-free region, no verified range of zeros, no result at any window above
`L~`, and nothing about the Riemann hypothesis. It settles no open question of
arXiv:2611.13129v1, and the conditioning spaces used here are **not** identified
with the prolate spaces of that paper. The finite-reduction strategy is from
arXiv:2608.24827v2 and is credited there; parity decoupling is likewise not new,
and that paper's Theorem 6.2 establishes the simplicity and evenness of the
ground state at `L = 4/5`, which this note does not. No priority is claimed.
