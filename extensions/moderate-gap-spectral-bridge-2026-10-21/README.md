# Moderate zero gaps, signed Gram matrices, and finite-window Weil positivity

An independent investigation of whether the techniques of

> César Angón, Erick Ross, Ze Sen Tang, Boyang Zhao,
> *Unconditional Bounds for Two, Three and Four Consecutive Moderate Gaps between
> Zeros of the Riemann Zeta Function*, arXiv:2610.22927v1 (16 Oct 2026)

yield new results for my finite-window Weil-positivity program
(`doi:10.5281/zenodo.23250932`, *A δ²–L³ Sensitivity Law for Finite-Window Weil
Positivity*).

**Date:** 2026-10-21.  **Author:** Elias De Jesús.
**Status:** investigation + verification package, kept in this repository as an
extension. **Not a manuscript, and it carries no DOI of its own.**

---

## Verdict

**E — the proposed transfer is obstructed**, with one **C**-level new inequality
and one **D**-level exact correspondence preserved. Per-route verdicts in
`11_FINAL_VERDICT.md`. Verdict **A is not awarded**; no new constraint on an
individual hypothetical off-critical-line zero was obtained.

**The three things worth keeping:**

1. **A new, costless strengthening of the paper's budget inequality** — the slack
   discarded in their Lemma 5 is the *only* `δ`-sensitive quantity in the whole
   Alpöge–Furman / Lamzouri / ARTZ framework, and it equals `(κ−1)²` with
   `κ = ∫η_τ²(u)cosh(2δu log T)du`. Its `δ²` coefficient is **twice the centred
   second moment of the window** — the same constant, and the same description, as
   in my published law.
2. **Anti-localisation** — every off-line zero has a partner at ordinate distance
   **exactly 0**, so it is structurally excluded from every run the moderate-gap
   theorem certifies. A 4-gap run gives **three** simple critical-line zeros, not
   five. No gap theorem can ever be localised at a hypothetical off-line zero.
3. **The non-additivity law** — the `δ`-sensitive functional is log-additive, not
   additive, over off-line zeros at ordinary spacing (`0.98 + 1.33 ln g`), because
   the relevant directions *align* as the depth grows. Had it been additive, the
   refined inequality would have proved the **Density Hypothesis**; that is why it
   cannot be, and why the literature extracts only the trace.

---

## Layout

```
01_SOURCES_AND_PROVENANCE.md      sources, hashes, chronology, notation fix, dependencies
02_NEW_PAPER_RECONSTRUCTION.md    ARTZ reconstructed: Z_T, kernel, G=P+Q, inertia, Theorem 1
03_GRAM_WEIL_COMPARISON.md        the exact dictionary, and the four structures it breaks
04_SPECTRAL_DEFECT_BUDGET.md      the refined budget; Routes A-D; the measured non-additivity
05_GAP_CONDITIONING.md            anti-localisation + the sub-Nyquist conditioning ceiling
06_CONFLUENT_BLIND_SPOTS.md       against the 2026-10-15 no-go: cases A-E
07_HORIZONTAL_DISPLACEMENT.md     the second-order expansion, sign, coefficient, uniformity
08_PARAMETER_SCALING.md           the exact conversion table and the five failure modes
09_POINTWISE_VS_AVERAGE.md        four explicit configurations; scope discipline
10_CANDIDATE_THEOREMS.md          T1-T7: proofs, counterexamples, missing lemmas
11_FINAL_VERDICT.md               verdicts, R1-R7, three next steps
THEOREM_STATUS.md                 every statement, its class, its tier
ASSUMPTIONS_AND_DEPENDENCIES.md   what is assumed and what is not verified
SHA256SUMS                        checksums for everything in the package

scripts/iv.py                     rigorous rational interval arithmetic (stdlib only)
scripts/linalg.py                 Jacobi eigensolver, PSD square root (diagnostic only)
scripts/verify.py                 the main suite: blocks A-J
scripts/cluster_scan.py           cluster cost vs size, eta^2=1 window
scripts/cluster_smooth.py         cluster cost vs size, admissible C^infty window
scripts/additivity.py             the additivity law and the alignment mechanism

results/verify.txt                captured output of the main suite
results/cluster_scan.txt          "
results/cluster_smooth.txt        "
results/additivity.txt            "
results/predecessor_hashes.txt    hashes of the antecedent PDFs
sources/SOURCE_HASHES.txt         fingerprints of the PDFs read (not redistributed)
```

## Reproducing

```
cd scripts
python3 verify.py            # main suite, exact/interval + diagnostic; ~5 min
python3 cluster_smooth.py    # admissible-window cluster costs;         ~3 min
python3 additivity.py        # the log-additivity law and alignment;    ~4 min
```
Python 3 standard library only. **No numpy, sympy, mpmath or pip on this machine** —
every sign-critical conclusion is therefore established in exact-rational interval
arithmetic (`scripts/iv.py`), and everything double-precision is labelled
`[diagnostic]`.

## Notation warning

ARTZ and Lamzouri use **δ** for the *smoothing cutoff*; my program uses **δ** for the
*horizontal displacement* `β − 1/2`. Throughout this package the cutoff is **τ**, the
displacement is **δ**, and `a = δ·log T` is the only dimensionless displacement.
See `01_SOURCES_AND_PROVENANCE.md` §3.

## Boundaries observed

The investigation was carried out entirely in a scratch directory. **No existing
manuscript, certificate, Zenodo record or `CITATION.cff` was modified**, no version
number was bumped, and no DOI was created; `src/`, `paper/`, `data/` and the other
extensions are untouched. This directory was added to `main` afterwards, at the
author's instruction, as an update of ongoing work — it is new material only.

No EOC/Collatz repository was touched; no researcher was contacted; no publication
manuscript was written. No Lean build was run and no formal verification is claimed,
for this work or for the Lean certificates reported by Alpöge–Furman and by
Lamzouri's Appendix A. **No statement in this package depends on the location of any
actual zero of ζ, and nothing here proves, advances, or assumes the Riemann
hypothesis.**
