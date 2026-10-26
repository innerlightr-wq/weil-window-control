# Suzuki spectral-function control

> **Verdict.** Phase 0 **PASS** · Phase 1 **PASS** · Phase 2 **CONDITIONAL** (on three
> unverified cited inputs; `X = 5, 7` additionally conditional on an unrefereed source) ·
> Phase 3 **PASS** · Phase 4 **MIXED/NULL** (preregistered prediction partly confirmed,
> partly refuted; its interpretation **not** supported).

Extension of [`weil-window-control`](../../README.md). **Carries no DOI, is not part of the
paper, and has not been peer-reviewed.** Preregistration (written before any computation):
[`PREREG.md`](PREREG.md). Full record with tiers and retractions: [`REPORT.md`](REPORT.md).

**Source:** M. Suzuki, *Spectral convergence associated with Weil's quadratic form*,
[arXiv:2610.32704v1](https://arxiv.org/abs/2610.32704).

---

## Non-claims — read these first

1. **Nothing here proves or advances the Riemann hypothesis.** Suzuki's Theorems 1.1, 1.2
   and Prop 6.1 all assume RH; the global space `H_W` requires it.
2. **Roots of `W_N` at a fixed window are local spectral values of a finite-dimensional
   Riesz approximation.** They are **not** zeros of `W`, and **not** zeros of `ζ`. Their
   numerical closeness to `γ_j` is an observation, not an identification.
3. **Suzuki's `a → ∞` limit stays RH-conditional** whatever is found here. No fixed-window
   result can relax it, because Theorem 1.1 runs through Theorem 3.3 and the global space.
4. **Planted spectra are fabricated diagnostics**, not `L`-functions of anything.
5. **Liu's half-width-1 and 17/16 certificates are UNREFEREED** (alphaXiv). Everything
   resting on them is labelled conditional on an unrefereed source, every time.

## What was established

**Phase 0 (GATE) — PASS.** `a = L` and `λ_a = λ*(L)` **exactly, constant 1**. Pinned
algebraically, then confirmed to **49 digits**: Suzuki's exponential-basis Gram, assembled
here independently, predicts this repository's own cosine-basis window form through
`cos(kπx/L) = (√(2L)/2)(χ_k + χ_{−k})`, with worst relative deviation `6.8e−49` against
`src/zeta_window.build_parts_idx` — code not written in this extension. `G` is exactly
real, symmetric, index-reversal symmetric and positive definite, all **derived** from the
arithmetic assembly with no reference to any zero of `ξ`.

**Phase 1 — PASS.** Suzuki's §6 reproduced: Table 1 to **all ten printed decimals** at
every `X`, Table 2's `d_6^{(N)}` and `κ₂(G_{a,18})` to their printed figures, the
forced-zero lattice `W_N(kπ/a) ≡ 0` for `|k| > N`, and his extra zero of `W_14` at
`30π/log 11 ≈ 39.3044` identified as lattice point `k = 15`.

**Phase 2 — CONDITIONAL.** Every inspectable step of Prop 6.1's proof factors through
coercivity `(2.4)` alone, so on a window where `λ_a > 0` is *certified* the conclusions
follow without assuming RH — **conditional on three cited inputs in papers not obtained**.
See [`notes/prop61_dependency.md`](notes/prop61_dependency.md).

**Phase 3 — PASS.** `λ_min(G_{a,N}) ≥ λ*(a)` wherever `λ*` is certified, and `ln κ₂` tracks
`−ln λ_min` as the Landau–Widom picture predicts.

**Phase 4 — MIXED, and the preregistered reading is NOT supported.** The function-field
analogue exists and is exact (`φ_i(w) = w^i`; Gram = the Hallouin–Perret Toeplitz matrix;
`J` = the Christoffel–Darboux kernel; `W(a,θ;·)` ↔ the para-orthogonal polynomials). The
preregistered prediction — that the `δ` where roots lose reality coincides with the `δ`
where positivity is lost — holds **exactly** at `R = 3, 4` and **fails by a factor 2.4** at
`R = 5`. The departure rate is **square-root**, `p ≈ 0.54–0.71`, not `δ²`. And where the
thresholds do coincide, the coincidence is the classical Szegő implication, **not** the
`δ²`–`L³` law — so "the `δ²`–`L³` law governs a second pipeline" is not supported either way.

## Six errors of my own, all logged

Four were caught by the Phase 0 gate before any result (a wrong pole-term conjugation; a
factor 2 in the prime term; a factor `2L` in the gate's own test code; a catastrophic
cancellation in the `d = 0` archimedean integrand). Two more came later: an overstated
"equivalence" in the function-field note, and an off-by-one in which Toeplitz matrix
governs the para-orthogonal zeros. **Two preregistration defects** are recorded rather than
silently repaired: an unreachable `1e−8` threshold in P0-T3, and — the one worth reading —
P1-T5, which asserted a discrepancy between the brief's restatement of (6.4) and the paper.
**There is no such discrepancy**: the typeset (6.4) carries an overline that both
`pdftotext` and `pypdf` silently drop. Rendering the page as an image settled it. Full list:
[`REPORT.md`](REPORT.md), Retractions 21–26.

## Layout

```
PREREG.md       preregistration, committed before any computation
REPORT.md       phase results, a tier on every claim, retractions
notes/
  normalization.md        the three conventions, and the closed forms
  prop61_dependency.md    where RH enters Prop 6.1, step by step
  ff_construction.md      the function-field analogue, decided before computing
  prior_2607.24830.md     what the earlier numerical realization already covers
src/            suzuki.py, phase0.py, phase1.py, phase23.py, phase4.py (standalone)
data/           JSON results and cached Gram matrices
results/        console transcripts, one per phase
```

## Licence

As the parent repository: **MIT** for code under `src/`, **CC BY 4.0** for the Markdown and
`data/`, `results/`.
