# REPORT.md — Suzuki spectral-function control

Lab record. Extension of `weil-window-control`; **no DOI, not part of the paper, not
peer-reviewed.** Preregistration, committed before any computation:
[`PREREG.md`](PREREG.md) (commit `e8b6a9d`). Derivations:
[`notes/`](notes/). Transcripts: [`results/`](results/).

Source: M. Suzuki, *Spectral convergence associated with Weil's quadratic form*,
[arXiv:2610.32704v1](https://arxiv.org/abs/2610.32704), §§1–3 and §6 read in full before
any code was written.

---

## 0. Non-claims

Fixed in PREREG before any number existed.

1. **Nothing here proves or advances RH.** Theorems 1.1, 1.2 and Prop 6.1 all assume it.
2. Roots of `W_N` at a **fixed** window are local spectral values of a finite-dimensional
   Riesz approximation — **not** zeros of `W`, **not** zeros of `ζ`.
3. Suzuki's `a → ∞` limit stays **RH-conditional** whatever is found here.
4. Planted spectra are **fabricated diagnostics**.
5. Liu's certificates are **UNREFEREED**; anything resting on them says so every time.

### Evidence tiers

**T1** exact/certified (exact over `ℚ`, or mpmath with precision-doubling agreement at a
stated tolerance) · **T2** numerical · **T3** interpretation. `float64` appears nowhere.
Non-detections certify nothing. Upper and lower bounds are never interchanged.

---

## 1. Verdicts

| phase | verdict | one line |
|---|---|---|
| 0 (GATE) | **PASS** | `a = L`, `λ_a = λ*(L)`, constant 1 — confirmed to 49 digits against the repository's own code |
| 1 | **PASS** | Suzuki's Tables 1 and 2, the forced-zero lattice, and the `39.3044` artifact all reproduced |
| 2 | **CONDITIONAL** | Prop 6.1 reduces to coercivity for every inspectable step; three cited inputs unverified |
| 3 | **PASS** | `λ_min(G) ≥ λ*(a)` where certified; `ln κ₂` tracks `−ln λ_min` |
| 4 | **MIXED / NULL** | prediction exact at `R = 3,4`, fails at `R = 5`; rate is square-root, not `δ²`; interpretation unsupported |

---

## 2. Phase 0 — normalization reconciliation (GATE) · **PASS**

Full derivation: [`notes/normalization.md`](notes/normalization.md). Transcript:
`results/phase0.txt`.

**Q0-a: `a = L`. Q0-b: `λ_a = λ*(L)`, constant `1`.** The Fourier conventions are
identical (`û(z) = ∫u e^{izx}dx` on both sides); the spectral sides are identical (sum over
distinct ordinates with multiplicity weights = sum over zeros with multiplicity); the
supports are identical; and Suzuki's `λ_a := inf σ(A_a)` is literally Zhu's
`inf Q(f)/‖f‖²`.

| test | result | tier |
|---|---|---|
| **P0-T1** algebraic identification | pinned exactly, constant `1` | T1 |
| **P0-T2** cross-basis vs `src/zeta_window.build_parts_idx` | worst relative deviation **`6.8e−49`** at `L = 0.5493`, **`8.5e−49`** at `L = 0.8047`, full `\|j\|,\|k\| ≤ 8` block (threshold was `1e−25`) | T1 |
| **P0-T3** spectral side (assumes RH) | deviation sign `= (−1)^{m+n}` in 5/5 pairs; `\|dev\|` ratios `1.551, 1.589` vs predicted `T`-ratios `1.699, 1.733` — the `1/T` truncation tail, to 8.7% | T2 |
| **P0-T4** structural | `G` real to `0.0`, symmetric to `0.0`, index-reversal symmetric to `0.0`, positive definite (`λ_min = 1.04e−7` at `a = 0.5493`, `2.00e−13` at `a = 0.8047`) | T1 |

**P0-T2 is the decisive one.** Suzuki's exponential-basis Gram, assembled here from an
independent closed form, predicts the repository's cosine-basis window form through
`cos(kπx/L) = (√(2L)/2)(χ_k + χ_{−k})`, `f_0 = √(2L)χ_0`. The repository's code was not
written in this extension and was not modified.

**One derivation worth recording.** Suzuki truncates the archimedean series and bounds the
tail at `10^{−188}`. That is **impossible by truncation** here: the summand bracket is
`O(1/k²)`, so `10^{−188}` would need `k ~ 10^{188}` terms. The closed form derived in
`notes/normalization.md` §4 replaces it with **one analytic, non-oscillatory 1-D integral
per entry** — no series tail to defend. Yoshida (5.15)/(5.16), which Suzuki cites for the
same object, was **not obtained**; the assembly here is independent and reconciled against
the repository instead.

---

## 3. Phase 1 — reproduce Suzuki §6 · **PASS**

Transcript: `results/phase1.txt`. Data: `data/phase1.json`. Gram matrices cached in
`data/gram_X*_N18_dps*.json`, built once at `N = 18` per `(X, dps)` and sliced for
`N = 14, 16` (the index sets are sub-blocks at the same `a`).

PHASE1_TABLES_PLACEHOLDER

---

## 4. Phase 2 — unconditional fixed-window statement · **CONDITIONAL**

Dependency analysis: [`notes/prop61_dependency.md`](notes/prop61_dependency.md).
Transcript: `results/phase23.txt`.

**The reduction.** Every step of Prop 6.1's proof that can be inspected in
arXiv:2610.32704v1 factors through coercivity `(2.4)` alone. RH's role is located
precisely: §2.2 takes `λ = 0` *"under RH"*, against the RH-free construction with `λ` below
`inf σ(A_a)` — so **RH is used to supply `λ_a > 0`**, exactly what a certificate replaces.
Prop 6.1 does **not** cite Theorem 3.3, the one irreducibly RH-dependent result (it builds
the global `H_W`, needing `Q_W ≥ 0` on all of `C_c^∞(ℝ)`); Thm 3.3 serves **Theorem 1.1**,
the `a → ∞` limit.

**Why CONDITIONAL and not PASS.** Three cited inputs are load-bearing and were **not
verified**, because the papers were not obtained: `[2, Lemma 2.2]` and `[10, §3.2]` (the two
form cores) and `[10, §3]` (the form-domain characterization). None is *about* zeros, so
none is expected to need RH — but expected is not checked.

**Recorded though it does not bite.** Prop 2.1's **third** de Branges axiom (norm
preservation under the nonreal-zero transformation) argues from `|F_w| = |F|` on the real
axis and therefore uses `Γ ⊂ ℝ` **essentially**, not via `(2.4)`. Prop 6.1's chain never
invokes that axiom — but its existence shows the reduction to `(2.4)` is **not automatic**,
which is the failure mode PREREG named in advance.

PHASE2_TABLES_PLACEHOLDER

---

## 5. Phase 3 — Landau–Widom consistency · **PASS**

PHASE3_TABLES_PLACEHOLDER

---

## 6. Phase 4 — function-field and planted-spectrum control · **MIXED / NULL**

Construction fixed before computing: [`notes/ff_construction.md`](notes/ff_construction.md).
Prior art recorded before Phase 4: [`notes/prior_2607.24830.md`](notes/prior_2607.24830.md).
Transcript: `results/phase4.txt`.

### 4a — the analogue exists; the preregistered KILL does not fire

`φ_i(w) = w^i` on the unit circle generates the Hallouin–Perret Toeplitz Gram from the
Frobenius angles exactly as `φ_{n,a}(γ)` generates `G_{a,N}` from the ordinates. The
singularity of `T_R` for `R > d` is a **range restriction** (the construction runs on
`R < d`), not an obstruction.

Where the analogue genuinely stops: `W_N` comes from the deficiency structure of `i d/dx`
— indices `(1,1)`, deficiency points `±i`, a `θ`-circle of self-adjoint extensions — and a
finite Toeplitz matrix has none of that. What replaces it is classical and already contains
the Gram: `T_R` is the moment matrix of the Frobenius measure on the circle, `J` is the
Christoffel–Darboux kernel **verbatim**, and the `θ`-family is the **para-orthogonal
polynomials** `B_R(z;β) = Φ_R(z) − βΦ*_R(z)`. So `W_N ↔ B_R` and *zeros real* ↔ *zeros on
the unit circle*.

### 4b — genuine spectra

PHASE4B_PLACEHOLDER

### 4c — planted off-circle spectra

PHASE4C_PLACEHOLDER

---

## 7. Retraction log

Continues the **main** series of the parent `REPORT.md`, which reaches 20. (The
Davenport–Heilbronn extension keeps a separate `R1–R6` series.)

**Retraction 21 — the pole term's conjugation.** `src/suzuki.py` first wrote the pole term
as `φ_m(i/2)conj(φ_n(i/2)) + c.c.` The correct term is `2 Re[φ_m(i/2)φ_n(i/2)]`: the
continuation of `û(t)conj(f̂(t))` is `û(z)conj(f̂(z̄))`, so
`h(i/2) = φ_m(i/2)conj(φ_n(−i/2)) = φ_m(i/2)φ_n(i/2)`. Both versions are real and they
**coincide for a real basis**, so the repository's real-basis code cannot distinguish them
— only the derivation does. Cross-check: the correct form reproduces `zeta_window`'s
`pole_sign = −1` for the odd sector automatically; the wrong one does not. *Found before any
Phase 1 result; no reported number was affected.*

**Retraction 22 — a factor 2 in the prime term.** First written as
`2Λ(n)n^{−1/2}S(log n)`; the correct coefficient is `Λ(n)n^{−1/2}S(log n)` with
`S = C(v)+C(−v)`. The repository's `2Λ/√n · Csym` equals this because `Csym = S/2` for a
real basis. *Found before any Phase 1 result.*

**Retraction 23 — a factor `2L` in the gate's own test code.** P0-T2 used `1/√(2L)` for the
`j = 0` basis coefficient instead of `√(2L)`: `χ_0` is the **constant** `1/√(2L)` and
`f_0 = 1`, so `f_0 = √(2L)χ_0`. This produced a **1.2–3.0 % cross-basis disagreement and an
initial Phase 0 FAIL**. The gate caught a bug in the gate. After the fix the agreement is
`6.8e−49`.

**Retraction 24 — catastrophic cancellation in the `d = 0` archimedean integrand.**
Computed as `S(v) − 2`, a difference of two numbers near `2`, which loses `~|log v|`
digits; tanh–sinh quadrature samples `v ~ 10^{−2·dps}`, and at `v = 10^{−55}` with dps 60
the naive form is wrong by a **relative `2.65e−7`**. Replaced by the half-angle form
`−4sin²(θ/2) − (v/a)cos θ`. Verified **not** to change the dps 40/60/90 values (which agree
to 28 digits), but necessary at the 140/180 digits Phase 1 uses.

**Retraction 25 — a discrepancy I reported in the source that does not exist.** PREREG
asserted that the brief's restatement `J(w,z) = φ(z)ᵀ conj(G^{-1}φ(w))` and "the paper's
(6.4)" `J(w,z) = φ(z)ᵀ G^{-1}φ(w)` were different functions. **They are — but the paper's
is the first one.** The typeset (6.4) carries an **overline over the whole group**
`G^{-1}φ(w)`, which **both `pdftotext` and `pypdf` silently drop**. Rendering page 15 as an
image shows it plainly. The brief transcribed the paper correctly; my text extraction did
not. Measured consequence: the conjugated form reproduces **6/6** of Suzuki's Table-1
values at `X = 5`; the overline-dropped misreading reproduces **2/6**, agreeing on the first
two roots and then diverging (`25.0108575815` against `25.0108575802`, and
`37.5862119696` against `37.5861810892`). *Lesson, recorded: do not infer a discrepancy in a
source from a text-extraction artifact — render the page.*

**Retraction 26 — an overstated equivalence of my own, and an off-by-one.**
(a) `notes/ff_construction.md` §4 asserted that positive definiteness of the moment matrix
is **equivalent** to all zeros of `B_R` lying on the circle, and used that to pre-emptively
weaken Phase 4c's interpretation. The classical theorem is **one-directional**
(positive definite ⟹ zeros on the circle); the converse is false, and `R = 5` exhibits it —
the zeros stay on the circle to `δ = 1.611` while positivity is lost at `0.662`.
(b) The first 4c run compared the zero-departure threshold against `λ_min(T_R)`, but `B_R`
has degree `R` and is built from `T_{R−1}`; against the correct matrix the agreement at
`R = 3, 4` is **exact** rather than `28.7 %` and `187 %` off. Both the wrong and the right
comparison are reported in `results/phase4.txt`.

## 8. Preregistration defects

Recorded rather than silently repaired, because a preregistration that is quietly edited is
not a preregistration.

**D1 — P0-T3's threshold was unreachable.** PREREG fixed `1e−8` absolute on the
arithmetic-vs-zero-side comparison. The zero sum truncated at `NZ` zeros has a tail
`~ 2(2a)³log(T/2π)/(2πT)`, so `1e−8` needs `T ~ 10^8`, i.e. `~10^7` zeros. Replaced —
**before** looking at whether the replacement passes — by the test that actually
discriminates a normalization error from a truncation tail: the deviation must carry the
sign `(−1)^{m+n}` and must fall like `1/T`. A fixed tolerance would not have distinguished
a constant offset from a shrinking tail; this does.

**D2 — P1-T5 rested on a false premise.** See Retraction 25. The test was rewritten, before
the data was examined, as a measurement of how far the overline-dropped misreading
diverges.

**D3 — Phase 4c's interpretation was weakened before the data existed, on a claim that
itself turned out to be too strong.** See Retraction 26(a). The amendment in
`notes/ff_construction.md` §4 stands as written and is corrected here.

## 9. Out of scope this round

Rigorous `W − W_N` error bounds via Prop 6.1's residual identity. Certified coercivity
supplies only the **lower-bound** ingredient and the constants are numerically vacuous at
these windows: Liu's `2^{−151}` gives `λ^{−1/2} ≈ 2^{75.5}`, so the Cauchy–Schwarz step in
Prop 6.1's proof carries no usable information here. Recorded as future work.

## 10. Repository effect

**None outside this directory.** No file in `paper/`, no core module in `src/`, and no part
of the v1.0 note was modified. The repository's `src/zeta_window.py`, `src/planted.py` and
`src/weilform.py` are **imported, never edited**. One line is added to the root README's
extension list, in a separate commit so it can be dropped independently.
