# Suzuki spectral-function control

> **Verdict.** Phase 0 **PASS** · Phase 1 **PASS** · Phase 2 **PASS**, scoped to `X = 3`
> (upgraded from CONDITIONAL once the three cited inputs were located and each verified
> RH-free; `X = 5, 7` stay conditional on an unrefereed source, `X = 11` numerical only) ·
> Phase 3 **PASS** · Phase 4 **MIXED/NULL** (preregistered prediction partly confirmed,
> partly refuted; its interpretation **not** supported).

Extension of [`weil-window-control`](../../README.md). **Carries no DOI, is not part of the
paper, and has not been peer-reviewed.** Preregistration (written before any computation):
[`PREREG.md`](PREREG.md). Full record with tiers and retractions: [`REPORT.md`](REPORT.md).

**Source:** M. Suzuki, *Spectral convergence associated with Weil's quadratic form*,
[arXiv:2610.32704v1](https://arxiv.org/abs/2610.32704).

---

## Non-claims — read these first

1. **Nothing here proves or advances the Riemann hypothesis.** Suzuki's Theorems 1.1 and
   1.2 assume RH irreducibly; the global space `H_W` requires it. Phase 2 does remove RH
   from **Prop 6.1** at the single window `a = ½log 3`, where `λ_a > 0` is certified rather
   than inferred — and that changes nothing about this non-claim. RH is *equivalent* to
   `λ_a > 0` for **every** `a` (Weil's criterion); positivity at finitely many windows was
   never the obstacle, and `λ_a > 0` for *some* window is already an unconditional theorem
   (arXiv:2606.09096 Thm 1.4). The certificate makes that theorem **effective** at a named
   `a`; it does not extend it.
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

**Phase 2 — PASS, scoped to `X = 3`.** Every step of Prop 6.1's proof factors through
coercivity `(2.4)` alone, so on a window where `λ_a > 0` is *certified* the conclusions
follow without assuming RH. The three cited inputs that were previously unverified —
`[2, Lemma 2.2]` and `[10, §3.2]` (the two form cores) and `[10, §3]` (the form-domain
characterization) — have been located in Connes–Consani (arXiv:2106.01715) and Suzuki
(arXiv:2606.09096), **read from page renders rather than text extraction**, and each is
**RH-free**. So at `a = ½log 3 = 0.5493`, inside Zhu's certified range, Prop 6.1 is
unconditional. Recorded and *not* erased by the upgrade: Prop 2.1's third de Branges axiom
genuinely needs `Γ ⊂ ℝ` — it simply is not invoked by Prop 6.1's chain.
See [`notes/prop61_dependency.md`](notes/prop61_dependency.md) and
[`notes/prop61_dependency_chain.md`](notes/prop61_dependency_chain.md).

**Phase 3 — PASS.** The result is a **matched-window comparison**: this extension's
`λ_min(G_{a,18}) = 1.7419e−17` at `a = 0.8047` against the repository's Gate B
`λ*(0.8) ≤ 2.2702e−17` and Zhu's certified `λ*(0.8) ≥ 8.9e−18`. The 0.59 % window offset is
stated, not smoothed over, and **Zhu's lower bound does not transfer** to the larger window
(monotonicity runs the wrong way), so nothing here certifies `λ_a > 0` at `X = 5`. The
finite-basis `λ_min` is a **variational upper bound** decreasing to `λ*` — that the limit is
exactly `λ*` is Connes–Consani Cor. 2.4, itself RH-free. The earlier `ln κ₂ / (−ln λ_min)`
agreement is **demoted to a sanity check** (Retraction 28): it is forced to within 7 % of 1
by arithmetic, since `κ₂ = λ_max/λ_min` with `λ_max ≈ 3`.

**Phase 4 — MIXED, and the preregistered reading is NOT supported.** The function-field
analogue exists and is exact (`φ_i(w) = w^i`; Gram = the Hallouin–Perret Toeplitz matrix;
`J` = the Christoffel–Darboux kernel; `W(a,θ;·)` ↔ the para-orthogonal polynomials). The
preregistered prediction — that the `δ` where roots lose reality coincides with the `δ`
where positivity is lost — holds at `R = 3, 4` **to within the bisection resolution**
(relative gap `≤ 1.5e−35`; the two bisections return bit-identical values, but they share a
bracket and a halving sequence, so that is a bound and not an identity — Retraction 30) and
**fails by a factor 2.4** at `R = 5`. The departure rate is **square-root**, `p ≈ 0.54–0.71`, not `δ²`. And where the
thresholds do coincide, the coincidence is the classical Szegő implication, **not** the
`δ²`–`L³` law — so "the `δ²`–`L³` law governs a second pipeline" is not supported either way.

## Ten errors of my own, all logged

Four were caught by the Phase 0 gate before any result (a wrong pole-term conjugation; a
factor 2 in the prime term; a factor `2L` in the gate's own test code; a catastrophic
cancellation in the `d = 0` archimedean integrand). Two more came later: an overstated
"equivalence" in the function-field note, and an off-by-one in which Toeplitz matrix
governs the para-orthogonal zeros.

**Four of the ten are the same error — a number that was really a floor, a cast or a
resolution, reported as if it were a measurement — and that pattern is itself the finding.**
Two precision floors (a `γ₆` transcribed at 10 decimal places, then roots serialized at 14
significant figures) silently capped the `d₆` column, and were caught not by reading the
numbers but by *a value repeating identically across rows that should differ*. Third: Phase
4b's claim that the POPUC zeros recover the Frobenius spectrum **"exactly"**, `0.000e+00`,
was a `float64` cast rounding `7.78e−62` down — the one place `float64` had survived here,
in a quantity carrying a headline; it is now stated as agreement **to `10^{−61}`**. Fourth:
Phase 4c's thresholds, reported as coinciding **"exactly"** with relative difference `0.0`
and genuinely bit-identical at 120 digits — because both come from one bisection helper
sharing a bracket and a halving sequence, so identical returns are built in; it is now a
**bound**, `≤ 1.5e−35`. None of the four changed a verdict. All four turned an artefact
into a claim, and the fourth was found only by going looking for the third's signature.

**Four preregistration defects** are recorded rather than silently repaired: an unreachable
`1e−8` threshold in P0-T3; a Phase 4c interpretation weakened on a claim that was itself too
strong; a P1-T3 non-vanishing test implemented as an arbitrary **relative-to-peak** cut that
failed 9 of 12 rows while the dichotomy it exists to measure held by a factor `10^110` at
every row (both the old cut and the preregistered one are now reported per row, so the
change is auditable); and — the one worth reading — P1-T5, which asserted a discrepancy
between the brief's restatement of (6.4) and the paper. **There is no such discrepancy**:
the typeset (6.4) carries an overline that both `pdftotext` and `pypdf` silently drop.
Rendering the page as an image settled it, and that lesson set the method for Phase 2, where
every cited statement was read from a page render rather than extracted text.

A further preregistered "test" is retracted as **not a test at all**: Phase 3's
`ln κ₂ / (−ln λ_min)` agreement is forced to within 7 % of 1 by arithmetic, since
`κ₂ = λ_max/λ_min` with `λ_max ≈ 3`. Its factor-3 threshold could not have failed. It is
demoted to a sanity check. Full list: [`REPORT.md`](REPORT.md), Retractions 21–30.

## Layout

```
PREREG.md       preregistration, committed before any computation
REPORT.md       phase results, a tier on every claim, retractions
notes/
  normalization.md        the three conventions, and the closed forms
  prop61_dependency.md       where RH enters Prop 6.1, step by step
  prop61_dependency_chain.md the three cited inputs, located and each verified RH-free
  ff_construction.md      the function-field analogue, decided before computing
  prior_2607.24830.md     what the earlier numerical realization already covers
src/            suzuki.py, phase0.py, phase1.py, phase23.py, phase4.py (standalone)
data/           JSON results and cached Gram matrices
results/        console transcripts, one per phase
```

## Licence

As the parent repository: **MIT** for code under `src/`, **CC BY 4.0** for the Markdown and
`data/`, `results/`.
