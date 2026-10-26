# PREREG.md — preregistration

**Written before any computation in this extension.** Every kill condition below is fixed
here; results go in `REPORT.md`. Where a later phase needs a threshold, it is stated here
even if the phase is never reached.

Source: M. Suzuki, *Spectral convergence associated with Weil's quadratic form*,
[arXiv:2610.32704v1](https://arxiv.org/abs/2610.32704), read in full for §§1–3 and §6
before any code was written. Tiering follows the repository: **T1** exact/certified (exact
over ℚ, or mpmath with precision-doubling agreement at a stated tolerance), **T2**
numerical, **T3** interpretation. `float64` only as a deliberate contrast. Non-detections
certify nothing. Upper and lower bounds are never interchanged.

---

## Standing non-claims (fixed now, repeated in README/REPORT)

1. **Nothing in this extension proves or advances RH.** Suzuki's Theorems 1.1, 1.2 and
   Prop 6.1 all assume RH; the global space `H_W` requires it.
2. Roots of `W_N` at a **fixed** window are local spectral values of a finite-dimensional
   Riesz approximation. They are **not** zeros of `W`, and **not** zeros of ζ. Any
   numerical closeness to `γ_j` is an observation, not an identification.
3. Suzuki's `a → ∞` limit (Thm 1.1, Thm 1.2) stays **RH-conditional** whatever this
   extension finds. No fixed-window result can relax it.
4. Planted spectra are **fabricated diagnostics**, not `L`-functions of anything.
5. Liu's half-width-1 and 17/16 certificates are **UNREFEREED** (alphaXiv). Anything
   resting on them is labelled conditional on an unrefereed source, every time.

---

## Phase 0 — Normalization reconciliation (**GATE**)

**What must be pinned.** An exact dictionary between (i) Suzuki's `Q_W` with his Fourier
convention `û(z) = ∫u e^{izx}dx`, `‖Fu‖² = 2π‖u‖²`,
`Q_W(u,f) = Σ_{γ∈Γ} m_γ û(γ) conj(f̂(γ))`; (ii) Zhu's
`λ*(L) = inf{Q(f)/‖f‖² : supp f ⊂ [−L,L]}`; (iii) this repository's window form
(`src/zeta_window.py`, `notes/stage3_assembly.md`) and its Gate B normalization.

**Two questions to answer with a yes/no and an explicit constant:**
- Q0-a: does Suzuki's `a` equal Zhu's `L`?
- Q0-b: does `λ_a` (Suzuki, `= inf σ(A_a)`) equal `λ*(L)` (Zhu), and if so up to what
  constant?

**Tests, all preregistered:**

- **P0-T1 (algebraic).** Write the three conventions side by side and identify
  `Q_W ↔ Q`. PASS iff the identification is exact with an explicitly named constant
  (`1` allowed).
- **P0-T2 (cross-basis numerical).** The repository's `EVEN` basis `cos(kπx/L)` is an
  exact real combination of Suzuki's `χ_n` at `a = L`:
  `cos(kπx/L) = (√(2L)/2)(χ_k + χ_{−k})`, and `χ_0 = 1/√(2L)`. So the repo's
  Gram entry must satisfy
  ```
  Q(f_j, f_k)/(‖f_j‖‖f_k‖) = (2L/4)·[G_{j,k}+G_{j,−k}+G_{−j,k}+G_{−j,−k}]/(‖f_j‖‖f_k‖)
  ```
  with `‖f_0‖² = 2L`, `‖f_k‖² = L` for `k ≥ 1`. Computed from
  `src/zeta_window.py:build_parts_idx` (code **not** written here) against `G` computed
  here. **PASS: max relative deviation ≤ 1e−25 at dps 60, on the full index block
  |j|,|k| ≤ 8, at a = L ∈ {0.5493, 0.8047}.** FAIL at > 1e−20.
- **P0-T3 (spectral-side numerical).** Under RH (assumed only for this check, and said
  so), `G_{m,n} = Σ_ρ φ_m(γ_ρ)φ_n(γ_ρ)` over the actual ζ zeros. Compare the arithmetic
  `G` against a zero-side partial sum with a tail estimate. **PASS: agreement to ≤ 1e−8
  absolute on the (0,0) entry and the `|m|,|n| ≤ 3` block at a = 0.8047**, with the tail
  estimated and reported, not hidden. This is a *consistency* check of the normalization,
  **T2**, and assumes RH — it cannot be used to support an unconditional claim.
- **P0-T4 (structural).** `G` must be real, symmetric, positive definite, and invariant
  under simultaneous index reversal `G_{−m,−n} = G_{m,n}` (which follows from `Γ = −Γ`,
  `m_{−γ} = m_γ`, `φ_{−n}(−γ) = φ_n(γ)` — and must therefore hold for the *arithmetic*
  formula too, with no reference to zeros). **PASS: all four hold, reversal symmetry to
  ≤ 1e−30 relative at dps 60.** A failure here means the arithmetic assembly is wrong.

> **KILL (Phase 0).** If P0-T1 cannot be pinned exactly, **or** P0-T2 fails at > 1e−20,
> **or** P0-T4 fails: the unconditional claim of Phase 2 is **STOPPED** and Phases 1, 3, 4
> continue only as conditional/numerical work. This is Retraction-18 territory — the
> failure mode being guarded against is substituting one normalization's certified bound
> into another's inequality.

---

## Phase 1 — Reproduce Suzuki §6 [T2, then T1 where possible]

Implement `G_{a,N}`, `κ₂(G_{a,N})`, `W_N` for `X ∈ {3,5,7,11}`, `N ∈ {14,16,18}`,
`a = ½ log X`, at **140 and 180 digits**, archimedean tail `< 1e−188` entrywise.

- **P1-T1 (Table 1).** First six positive roots of `W_18`, all four `X`.
  **PASS: agreement with every printed digit of Suzuki's Table 1 (10 d.p.).**
- **P1-T2 (Table 2).** `d₆^{(N)} = t_{a,N,6} − γ₆` for `N = 14,16,18` and `κ₂(G_{a,18})`.
  **PASS: same leading digits and same exponent as printed** (Table 2 gives 2–3
  significant figures, so the match criterion is those figures).
- **P1-T3 (forced-zero lattice).** Derived in chat and to be verified here: at
  `z = kπ/a` with `|k| > N`, every `φ_{n,a}(z) = 0` (since
  `φ_n(kπ/a) = √(2a)·sin(kπ+nπ)/(kπ+nπ) = 0` for `k+n ≠ 0`), hence `W_N(kπ/a) = 0`
  identically. **PASS: `|W_N(kπ/a)|` at the precision floor for all `|k| > N` tested, and
  bounded away from zero at sampled `|k| ≤ N`, at every `(a,N)`.**
- **P1-T4 (the 39.3044 artifact).** Suzuki reports an extra zero of `W_14` at
  `30π/log 11 ≈ 39.3044` for `X = 11`, i.e. `k = 15 > N = 14`.
  **PASS: reproduced, and identified as lattice point `k = 15`.** Preregistered
  corollary to check: for `N = 18` the first lattice point is `k = 19`, i.e.
  `19π/a = 38π/log 11 ≈ 49.79`, **outside** Suzuki's `[0,40]` search range — so no such
  artifact should appear in his `N = 18` table. Verify and state.
- **P1-T5 (`J` convention).** The brief's restatement
  `J(w,z) = φ(z)^T conj(G^{-1}φ(w))` and the paper's `(6.4)`
  `J(w,z) = φ(z)^T G^{-1} φ(w)` are **not the same function** when `w = ±i` (they differ
  by `w ↦ w̄`, which swaps the two terms of `W_N`). Both make `−iW_N` real on ℝ.
  **Resolution is empirical: whichever reproduces Table 1 is adopted, and the other is
  recorded as a discrepancy in the brief's restatement, not in the paper.**

**Mismatch policy.** Any disagreement with Suzuki's printed digits is investigated and
logged in `REPORT.md` with its cause. **Nothing is tuned to make a number match.** If a
mismatch cannot be explained, it is reported as an open discrepancy.

---

## Phase 2 — Unconditional fixed-window statement [T1] (only if Phase 0 passes)

**The claim under test.** Prop 6.1's proof uses RH only through:
(i) coercivity `(2.4)`: `Q_W(u) ≥ λ_a‖u‖²`, `λ_a > 0`, on `C_c^∞(−a,a)`;
(ii) whatever else the proof invokes.
**Deliverable:** `notes/prop61_dependency.md`, a line-by-line pass over Suzuki's proof of
Prop 6.1 **and every result it cites**, marking each use of RH and whether it factors
through `(2.4)` alone.

> **PREREGISTERED HONESTY CONDITION.** The claim "Prop 6.1 is unconditional given certified
> `λ_a > 0`" is to be *tested*, not assumed. If the dependency analysis finds **any** use of
> RH not reducible to `(2.4)` — including inside `[10]`'s cited lemmas, the definition of
> `H(A_a)`, the de Branges property (Prop 2.1), Lemma 3.2, Thm 3.3, or the cores argument
> — the unconditional claim is **withdrawn for that part** and reported as such. The
> expected failure mode is that `H(A_a)` itself, or the identification with
> `D(A_a^{1/2})`, needs more than positivity.

- **P2-T1 (dependency).** Verdict per proof step: `RH-free` / `via (2.4) only` /
  `needs RH`. No threshold; the output is the classification.
- **P2-T2 (certified roots at X = 3).** `a = ½log 3 = 0.5493` lies inside Zhu's certified
  range `L ≤ 0.8` (`λ*(0.8) ≥ 8.9e−18`, and `λ*` is non-increasing in `L`, so
  `λ*(0.5493) ≥ 8.9e−18 > 0`). Certify the first several real roots of `W_N` by:
  sign changes on a grid, **precision-doubling stability** (dps 140 vs 180), **exclusion
  of the forced lattice** `kπ/a`, and a root count on a contour if feasible.
  **PASS: roots stable to ≤ 1e−30 under precision doubling, none on the lattice, and the
  contour count (if run) matching the grid count.**
  **Reported as "certified roots of the finite approximation `W_N`". NOT as zeros of `W`.
  NOT as zeta zeros.**
- **P2-T3 (X = 5, 7).** Same computation, labelled **conditional on Liu's unrefereed
  certificate**. Note explicitly that `a(X=5) = 0.8047 > 0.8` is **outside Zhu's range**.
- **P2-T4 (X = 11).** `a = 1.1990`. No certified window. **Numerical only**, labelled so.
- **P2-T5 (prime-free window control).** For `2a < log 2` (`a < 0.3466`, i.e. `X < 2`) the
  prime sum is empty and `Q_W` is archimedean + pole only. Preregistered check: `G` is
  then computable in closed form with no prime input, and `λ_min(G) > 0` should hold with
  room. Used as a positive control on the assembly.

---

## Phase 3 — Landau–Widom consistency [T2]

Compare `λ_min(G_{a,N})` and `κ₂(G_{a,N})` from Phase 1 with this repository's `λ*(L)`
data and Gate B values at matched windows (after the Phase 0 mapping).

- **P3-T1.** `λ_min(G_{a,N}) ≥ λ*(a)` must hold for every `(a,N)`, since `G` is the form
  restricted to a subspace and `λ*` is the infimum over the whole admissible class — *but
  only after accounting for the basis normalization*: `G` is the Gram in the `Q_W` inner
  product on an **orthonormal** `L²` basis, so `λ_min(G) = min_{c} Q_W(Σc_nχ_n)/‖Σc_nχ_n‖²`
  is a Rayleigh quotient over a subspace, hence `≥ λ*(a)`. **PASS: the inequality holds at
  every tested `(a,N)`. A violation means the Phase 0 mapping is wrong** and feeds back
  into the Phase 0 gate.
- **P3-T2 (preregistered growth check).** Zhu's law `−ln λ*(L) ≈ C·N(T*)/ln N(T*)`,
  `T* = 2πe^{2L}`, `C ≈ 2π²`, predicts `λ*` falling super-exponentially in `L`. Since
  `κ₂ = λ_max/λ_min` and `λ_max(G_{a,N})` is `O(1)`-ish, the prediction is that
  `ln κ₂(G_{a,18})` tracks `−ln λ_min(G_{a,18})` and grows with `X` at a rate comparable
  to the repo's `−ln λ*(L)` column, **up to basis-dependent factors**.
  **Report: the measured `ln κ₂` against the repo's `−ln λ*` at matched `L`, with the
  ratio. Agreement within a factor 3 → PASS. Otherwise report the discrepancy as a
  finding; do NOT tune.** This is a consistency check between two finite bases, not a
  test of Zhu's law — the repo already records (projection-formula Stage 4) that its own
  `−ln λ*` tracks Zhu's law only to a steadily worsening ratio.

---

## Phase 4 — Function-field and planted-spectrum control (**main research test**)

Before any code: read [arXiv:2607.24830](https://arxiv.org/abs/2607.24830) (numerical
realization of Suzuki's earlier operator, including an injected off-line zero) and record
in `notes/prior_2607.24830.md` what it already covers, so this phase adds what it does not.

### 4a — the function-field analogue
Build the `W_N` analogue from the Hallouin–Perret Frobenius Gram matrix already in
`src/stage1.py`: same construction (Gram solve against basis transforms at `z = ±i`, then
the `W_N` formula), with the basis transforms defined consistently with that Gram matrix.
**Document in `notes/ff_construction.md` BEFORE computing.**
> **KILL (4a).** If no faithful analogue of `φ_{n,a}` exists for that basis — i.e. if there
> is no function `φ_n(z)` whose values at the Frobenius angles generate the Gram matrix in
> the way `φ_{n,a}(γ)` generates `G_{a,N}` — **STOP 4a and record why.** The specific risk:
> the function-field Gram is a *finite* Toeplitz matrix in `t(n) = Σ_j 2cos(nθ_j)` with a
> kernel for `R > d`, so `G` is **singular**, and `G^{-1}` in (6.4) does not exist.

### 4b — genuine spectra (RH = Weil's theorem here)
Check that the roots of the `W_N` analogue recover the Frobenius angles.
**PASS: every angle `θ_j` matched to ≤ 1e−20 (exact over ℚ where available), with no
spurious roots other than identified lattice artifacts.**

### 4c — planted off-circle spectra (reuse `src/stage2.py`)
Vary the displacement `δ` and record, per window:
- (i) the `δ` at which `λ_min` turns negative — predicted by the repo's `δ²`–`L³` law with
  `C ≍ (window length)³/6`;
- (ii) the `δ` at which `W_N` first acquires non-real roots.

> **PREREGISTERED PREDICTION.** (ii) coincides with (i). Reality of `W`'s zeros comes from
> positivity — they are eigenvalues of a self-adjoint operator, and the construction of
> `H(A_a)` needs `Q_W > 0` — so loss of reality should track loss of positivity at the
> `δ²`–`L³` threshold.
>
> **Tolerance fixed now:** "coincide" means `|δ_(ii) − δ_(i)| / δ_(i) ≤ 0.10`.
>
> **OUTCOMES.**
> - coincide within 10% → **"the `δ²`–`L³` law governs a second pipeline"** (T2; T1 where
>   exact over ℚ);
> - a systematic offset reproducible across windows and curves → **report the offset law**
>   as the finding, with its fitted form and residuals;
> - no clean relation → **report negative.**
>
> **A ratio measured on one curve at one `R` is not an offset law.** At least three
> `(curve, R)` pairs must show the same ratio to ±20% before the word "law" is used.

### 4d — optional ζ-side
Planted spectra via the spectral side `Q_W = Σ m_γ û f̂‾` with fabricated zero sets (repo
Part II convention), same comparison. **Clearly labelled fabricated diagnostics.** Run only
if 4a–4c complete.

---

## Out of scope this round (recorded, not attempted)

Rigorous `W − W_N` error bounds via Prop 6.1's residual identity. Certified coercivity
supplies only the **lower-bound** ingredient, and the constants are tiny: Liu's `2^{−151}`
gives `λ^{−1/2} ≈ 2^{75.5}`, so the Cauchy–Schwarz step in Prop 6.1's proof is numerically
vacuous at these windows. Noted in `REPORT.md` as future work.

---

## Verdict line format

Per phase: **PASS** / **CONDITIONAL** / **NULL** / **STOP**, with the condition named when
CONDITIONAL.

## Retractions

Logged in `REPORT.md` as they happen, continuing this repository's numbering (the parent
`REPORT.md` reaches Retraction 20; the Davenport–Heilbronn extension keeps a separate R1–R6
series, so this extension continues the **main** series from **Retraction 21**).
