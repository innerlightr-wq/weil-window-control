# Weil window control — Stages 1 and 2

Control experiments for the Connes / Connes–van Suijlekom (CvS) "finite-prime Weil form"
strategy, run in the function-field setting where the answer is known
(Weil's theorem = RH for curves), plus planted off-line configurations.

**Tiering.** T1 = exact / proved. T2 = numerical, with stated working precision and a
convergence check. T3 = interpretation.

All eigenvalue computations use `mpmath` at 50 decimal digits (sometimes 80) with the
sign of `λ_min` asserted only when `|λ_min| > 100·|λ_min(dps) − λ_min(2·dps)|`.
Where the data is rational, the inertia is additionally certified **exactly over Q**
by the Sylvester/Jacobi leading-minor rule (`src/exact_inertia.py`) — those rows are T1.
Backend: pure-Python mpmath (no gmpy2 on this machine).

---

## Stage 1 — function-field control

### Setup, derived not assumed (T1, `notes/stage1_derivation.md`)

From `Z_C(T) = exp(Σ N_n T^n/n) = P(T)/((1−T)(1−qT))`, `P(T) = Π(1 − α_j T)`:

```
N_n = 1 + q^n − p_n,        p_n = Σ_j α_j^n
t(n) := Σ_j cos(n θ_j) = p_n q^{−n/2} = q^{n/2} + q^{−n/2} − N_n q^{−n/2},   t(0) = 2g
```

which is the identity in the brief. `(T_R)_{ij} = t(|i−j|)` for `i,j = 0..R`.

**Weil positivity is automatic given `|α_j| = √q` (T1).** With `(w_m)_i = e^{i i θ_m}`,

```
T_R = Σ_{m=1}^{2g} w_m w_m^H,   so   Q(c) = Σ_m |ĉ(θ_m)|² ≥ 0,   rank T_R ≤ min(R+1, d)
```

where `d = #distinct θ_m ≤ 2g`. Confirmed: `inertia_neg = 0` in **every** row of
`data/stage1_spectra.csv` (40 rows, 4 curves, R = 0..4g).

### Curves (T1 — exact integer point counts)

| tag | curve | q | g | N₁..N₃ | P(T) coefficients |
|---|---|---|---|---|---|
| E5 | y²=x³+x | 5 | 1 | 4, 32, 148 | 1, −2, 5 |
| H2F3 | y²=x⁵+x²+1 | 3 | 2 | 6, 18, 18 | 1, 2, 6, 6, 9 |
| G3F5 | y²=x⁷+x+1 | 5 | 3 | 9, 35, 123 | 1, 3, 9, 17, 45, 75, 125 |
| G3F3 | y²=x⁷+x³+x+1 | 3 | 3 | 7, 11, 37 | 1, 3, 5, 9, 15, 27, 27 |

H2F3 reproduces the note's `N₁ = 6, N₂ = 18`. No Sage/PARI was available, so the
L-polynomials are cross-checked two independent ways instead (see H3).

`max | |α_j|/√q − 1 | ≤ 1.34e−51` on all four curves (T2, 50 dps) — Weil's theorem,
numerically confirmed.

### Hypotheses

**H1 — CONFIRMED, and it is RH-free (T1 reasoning + T2 evidence).**
*If `λ_min(T_R)` is simple, the ground-state polynomial has all `R` zeros on `|z| = 1`.*
`T_R − λ_min I` is Toeplitz, PSD **by construction** (λ_min is the smallest eigenvalue —
the sign of λ_min is irrelevant), singular, and of rank exactly `R` when λ_min is simple.
Carathéodory–Fejér / Pisarenko then forces `R` simple zeros on the unit circle.
Measured `max | |z| − 1 |` ≤ 2.4e−50 at every applicable row (50 dps).
Adversarially re-tested on 400 random real symmetric Toeplitz matrices with simple
`λ_min`: worst deviation 7.5e−38 at 40 dps (`src/adversarial.py`, A1). **H1 survives.**

**H2 — CONFIRMED as stated for these curves, but REFUTED as a general statement; restated.**
At `R = 2g`: `λ_min = 0` (|λ| ≤ 3e−50), **simple**, with a healthy gap (1.38 – 3.30),
ground state **even**, and zero angles matching `{θ_j}` to ≤ 6e−50. Exact finite-window
zero recovery.
Moreover the kernel vector is predicted in closed form (T1): the kernel polynomial is
`Π_j (z − α_j/√q) = Σ_k A_{2g−k} q^{−(2g−k)/2} z^k`, which the functional equation
`A_{2g−i} = q^{g−i} A_i` collapses to

```
c_k ∝ A_k q^{−k/2} ,
```

visibly **palindromic — so evenness of the R = 2g ground state is forced, not observed.**
Measured agreement `‖v − c_pred‖_∞ ≤ 1.3e−50` on all four curves.

*The refutation:* H2 needs the `θ_j` to be **distinct**. Attack A2 found real
counterexamples — e.g. **y² = x⁵ + x over F₅** has `P(T) = 1 + 10T² + 25T⁴ = (1+5T²)²`,
so `θ_j = ±π/2` each twice, `d = 2 < 2g = 4`. There `rank T_R = 2`, exact recovery
happens at `R = 2 = d` (λ_min = 0 simple, zeros exactly `e^{±iπ/2}`), and at `R = 2g = 4`
λ_min = 0 has multiplicity 3. **Restated H2 (T1):** at `R = d = #distinct θ_j`, `T_R` has
rank `d`, `λ_min = 0` is simple, and the kernel polynomial's zeros are exactly the
distinct `e^{iθ_j}`. The original form is the `d = 2g` case.

**H3 — CONFIRMED (T1, exact integers).** `P` built from `N_1..N_g` alone (Newton +
functional equation) predicted every further brute-forced point count exactly:
6, 9, 4 and 8 independent integer predictions for E5, H2F3, G3F5, G3F3. The L-polynomial
rebuilt from all of `p_1..p_{2g}` without the functional equation was identical in all
four cases.

**H4 — CONFIRMED but NOT monotone; the right monotone quantity is `λ_min`, not the zeros.**
Coverage error `max_j min_z |θ_j − arg z|` for `R = 0..2g`:

| curve | R=1 | R=2 | R=3 | R=4 | R=5 | R=6 |
|---|---|---|---|---|---|---|
| H2F3 (g=2) | 1.571 | 0.342 | 0.373 | **1.1e−50** | — | — |
| G3F5 (g=3) | 2.097 | 0.990 | 0.671 | 0.420 | 0.344 | **2.7e−51** |
| G3F3 (g=3) | 2.130 | 1.281 | 0.559 | 0.646 | 0.196 | **5.3e−50** |

Generally decreasing, strictly so only for G3F5 — so H4 as a monotone-error statement is
**not** supported. Two structural caveats (T1): (a) with `R < 2g` there are only `R`
zeros for `2g` angles, so the error cannot vanish; (b) parity forces a spurious zero —
an odd ground state always has a zero at `z = 1`, an even one at odd `R` always has one
at `z = −1`, which is what dominates the reverse (`zero → θ`) error at odd `R`
(e.g. exactly `π − 2.533 = 0.608` for G3F5 at R = 1, 3, 5).

**The monotone quantity is `λ_min(R)`,** which decreases to 0 by Cauchy interlacing
(T1: `T_R` is a principal submatrix of `T_{R+1}`). Verified with no violations at
tolerance 1e−40, e.g. G3F5: 6, 4.658, 3.000, 2.220, 0.889, 0.335, 0.

### Recorder split `T = A − 2M`

With `A = 2g·I + pole part (q^{n/2}+q^{−n/2})`, `2M = prime part (N_n q^{−n/2})`
(split residual ≤ 3.1e−46 at 50 dps).

**Block F holds**: for `R` large enough, `A` and `M` are both indefinite and
`S = T_R` is PSD, on all four curves.

**The split is convention-dependent**, exactly as in the 389a1 appendix. The only
ambiguity is the `n = 0` bookkeeping, and it changes `A` completely:

- brief convention (`A(0) = 2g+2`, `2M(0) = 2`): `A` inertia `(1, 0, R)` — full rank.
- uniform convention (`A(0) = 2`, `2M(0) = N_0 := 2−2g`): `A` inertia `(1, R−1, 1)` — **rank 2**.

The rank-2 fact is exact (T1): `q^{|n|/2} + q^{−|n|/2} = u_i w_j + w_i u_j` with
`u_i = q^{i/2}`, `w_i = q^{−i/2}`, so the pole part is a single hyperbolic plane.
Residual 0.0 – 1.1e−47. **The pole part is indefinite precisely because the poles of
`ζ_C` sit off the critical circle** — it is the `|β| ≠ 1` signature of Stage 2, appearing
already in the honest case.

### Task 1e — `T_R` is a Gram matrix in the trace form on `H¹` (T1 + T2)

Let `F` be Frobenius on `H¹` (dim 2g), with characteristic polynomial the reversed
L-polynomial `Σ A_i T^{2g−i}` — an integer companion matrix. The polarization supplies a
nondegenerate Hermitian `S` with `F* S F = q S`; the adjoint is `φ† = S^{-1} φ* S`, which
gives `F† = q F^{-1}`, so `U := F/√q` has `U† = U^{-1}` and

```
Gram_ij = Tr( U^i (U^j)† ) = Tr( U^{i−j} ) = Σ_m β_m^{i−j} = t(|i−j|),
```

i.e. **`T_R` is exactly the Gram matrix of `{1, F/√q, …, F^R/q^{R/2}}` in `Tr(φψ†)`**.
Cor. 14 is the `g = 1, R = 1` case. Verified for all four curves at `R ≤ 2g+2`
(`src/gram_trace.py`, `data/stage1e_gram_trace.csv`, 40 dps):

| curve | `F*SF − qS` | `UU† − I` | inertia of S | `max |Gram − T|` |
|---|---|---|---|---|
| E5 (g=1) | 2.5e−40 | 9.2e−41 | (0,0,2) pos.def. | 1.8e−40 |
| H2F3 (g=2) | 9.8e−40 | 3.6e−39 | (0,0,4) pos.def. | 2.8e−40 |
| G3F5 (g=3) | 7.8e−40 | 2.0e−38 | (0,0,6) pos.def. | 2.8e−40 |
| G3F3 (g=3) | 2.4e−40 | 1.1e−38 | (0,0,6) pos.def. | 1.1e−39 |

**It does not fail — and locating exactly where RH enters is the useful part.** Writing
`F = V diag(α) V^{-1}`, the constraint in the eigenbasis reads `conj(β_j) β_k S_jk = q S_jk`,
so `S_jk = 0` unless `β_k = q/conj(β_j)`. The functional equation makes that pairing an
involution, so a Hermitian `S` always exists — but `S_jj ≠ 0` requires `|β_j|² = q`.
**If any eigenvalue is off the circle its diagonal entry of `S` is forced to zero, and a
Hermitian matrix with a zero diagonal entry is never positive definite.** Hence (T1)

```
all |α_j| = √q   ⟺   a positive-definite q-isometry form exists
                 ⟺   φ ↦ Tr(φφ†) is a positive form   ⟺   T_R is PSD for all R.
```

**That equivalence is a classical linear-algebra fact, not a finding of this work.**
"An automorphism admits an invariant positive-definite Hermitian form iff its eigenvalues
all have the same modulus" is standard, and it is precisely why Weil's proof of RH for
curves runs through positivity of the Rosati involution. I state it here only to locate
where RH enters the window form; no novelty is claimed for it.

**What this section does contribute is the bridge**: that the *window Toeplitz matrix*
`T_R` of the Connes/CvS pipeline is exactly the Rosati Gram matrix of
`{1, F/√q, …, F^R/q^{R/2}}`, so the pipeline's finite-window object and Weil's
Rosati-positivity object are the same matrix; and that on planted spectra the signature of
the polarization `S` reproduces the Stage 2 Toeplitz signature term for term.

Run on the Stage 2 planted spectra, the **Gram identity still holds** (residuals
≤ 7.3e−39) while `S` goes indefinite, with inertia equal to the Stage 2 Toeplitz
signature in every case:

| config | forced-zero diagonals | inertia of S | stable inertia of `T_R` |
|---|---|---|---|
| C0-online | 0/4 | (0,0,4) pos.def. | (0,·,4) |
| C1-g1real | 2/2 | **(1,0,1)** | (1,·,1) |
| C2-quartet | 4/4 | **(2,0,2)** | (2,·,2) |
| C4-mixed | 4/6 | **(2,0,4)** | (2,·,4) |

So the `H¹` trace-form picture and the Toeplitz picture are the same object: the
signature of the polarization *is* the stable signature of the window form.

---

## Stage 2 — teeth (planted off-circle configurations)

Writing `β_j = α_j/√q`, the functional equation is `β → 1/β` and RH-for-curves is
`|β| = 1`. **`t(n) = Σ_j β_j^n` depends only on the normalised `β_j`, so `q` drops out of
Stage 2 entirely** (T1). Blocks: on-line pair `{e^{±iφ}}`; off-line quartet
`{ρ e^{±iφ}, ρ^{−1} e^{±iφ}}` = `{β, 1/β, β̄, 1/β̄}`; off-line real pair `{ρ, 1/ρ}`
(the `g = 1`, `|a| > 2√q` tooth, here `a = 2.069√q`).

Signature per block (T1, same `u/w` identity): on-line → rank 2, **PSD**;
off-line quartet → rank 4, signature **(2,·,2)**; off-line real pair → rank 2,
signature **(1,·,1)**. Predicted total signatures matched the computed inertia in
**every** row of `data/stage2_spectra.csv` once `R+1 ≥ rank`.

### Main table

| config | planted | first R with λ_min<0 | exact over Q | H1 holds? | λ_min simple? | definite parity? |
|---|---|---|---|---|---|---|
| C0-online (control) | all on-line, φ = ±0.7, ±2.0 | **never** | yes, and PSD by T1 proof | yes, always | fails for R ≥ 5 (kernel) | fails where non-simple |
| C1-g1real | one real pair ρ=1.3 | **1** | **R = 1** | yes, always | never fails | never fails (odd) |
| C2-quartet | quartet ρ=1.3, φ=0.7 | **2** | **R = 2** | yes, always | never fails | never fails |
| C3-twoquartets | quartets ρ=1.3 at 0.7 and 2.0 | **4** | **R = 4** | yes, always | never fails | never fails |
| C4-mixed | quartet (1.3, 0.7) + on-line pair at 2.0 | **4** | **R = 4** | yes, always | never fails | never fails |

The "exact over Q" column uses rational surrogate angles (`cos φ = 19/25 → φ = 0.70137`,
`cos φ = −21/50 → φ = 2.00459`) so every `t(n)` is rational — those first-negative
windows are **T1**, not T2. By Cauchy interlacing, `λ_min < 0` at `R₀` implies
`λ_min < 0` for every `R ≥ R₀`, so the whole tail is certified too.

### Findings

**F1. The witness is free (T1 + T2) — and CvS already say so implicitly.** Their theorem
is stated for *any* real even distribution with the stated spectral hypothesis; no
arithmetic enters, so nothing in step (i) can distinguish an arithmetic spectrum from a
fabricated one. What the planted control adds is a concrete demonstration with numbers
attached, not the observation itself. H1 holds at **every** applicable row of **every**
planted configuration, including rows where `λ_min = −642.8`. Measured
`max||z|−1| ≤ 5.1e−48`. *(Corrected in Round 3: the figure 1.6e−48 first reported here was the maximum over a subset of rows; over all applicable planted rows it is 5.03e−48.)* All zeros of the approximant lie exactly on the critical circle
whether or not the underlying spectrum does. *"The approximant's zeros are real" carries
no RH content whatsoever.*

**F2a. The kernel is an ideal; the CvS hypothesis fails for every `R > d` (T1 + T2).**
Let `P*(z) = Π (z − β_j)` over the **distinct** `β_j`, `d = deg P*`. Then
`Q(c) = 0 ⟺ ĉ(β_j) = 0 ∀j ⟺ P* | ĉ`, so

```
ker T_R = { P*(z) q(z) : deg q ≤ R − d },     dim ker T_R = max(0, R + 1 − d).
```

Predicted nullity matched the computed nullity in **every** row of
`data/stage2_function_level.csv` (all four curves, the repeated-angle curve, and the
on-line control, `R` up to 14). At `R = d` the kernel is 1-dimensional — `λ_min = 0` is
simple and the ground state *is* `P*`. **At `R = d+1` it is 2-dimensional, so `λ_min = 0`
stops being simple and the CvS hypothesis (simple isolated lowest eigenvalue) fails for
every `R > d`.** A finite spectrum runs out of content past its own degree; for an honest
curve with distinct angles that threshold is exactly `R = 2g`.

**This explains the evenness failure, which is not a failure (T1).** `P*` is
self-reciprocal up to sign (the spectrum is closed under `β → 1/β`), so `J(P*q) = ±P*q^rev`:
the kernel is `J`-invariant and splits into even and odd parts of dimensions
`⌈(m+1)/2⌉` and `⌊(m+1)/2⌋` with `m = R − d`. Measured even/odd kernel dimensions:
`1/0, 1/1, 2/1, 2/2, 3/2, 3/3, …` for `m = 0,1,2,3,4,5` — exactly the prediction.
For `m ≥ 1` **both** parts are nonzero, so a generic kernel vector — which is what any
eigensolver returns — has no definite parity. The "mixed parity" rows at `R > 2g` in
Stage 1 are an artefact of degeneracy, not a property of the curve.

**F2. Simplicity and evenness are not diagnostics (T2).** In C1–C4, `λ_min` is simple and
the ground state has definite parity at **every** `R` from 0 to 20 — the planted
violation never disturbs either. Meanwhile parity alternates even/odd with `R` in the
honest control C0 too. So failure of simplicity/evenness tracks window bookkeeping and
spectral degeneracy, not the position of the zeros.

**F3. Zero ANGLES converge; the minimizer FUNCTION does not (T1 + T2).**
*(This finding supersedes an earlier version of F3; see retraction 6.)*

*Angles — free.* In C2 (planted `φ = 0.7`) the ground-state zero angles go
0.69036 → 0.70012 → 0.69804 → 0.69834 → 0.69920 for R = 2, 6, 12, 16, 20; C3 and C4 lock
onto both 0.7 and 2.0, reaching ~1e−3 by R = 30. The approximant recovers the angular
(ordinate) data of an off-line spectrum perfectly.

*Functions — not free, and provably so.*

> **Theorem (T1).** Suppose `λ_min(T_R)` is simple for infinitely many `R` and the
> normalized ground-state polynomials `ĉ_R` converge locally uniformly on `C` to some
> `F ≢ 0`. By H1 every zero of every `ĉ_R` lies on `|z| = 1`; by Hurwitz every zero of
> `F` is a limit of zeros of the `ĉ_R`, hence lies on `|z| = 1`. So if the target has
> **any** zero off the unit circle, the minimizer functions cannot converge to it.

**H1 — the free witness — is precisely the obstruction to function-level convergence in
the off-line case. (i) and (iv) are two sides of one coin.** Three measurements:

- **M1, coefficient distance at `R = d`** between the unit-normalized ground state and
  the unit-normalized target. Honest curves and the on-line control: ≤ 1.5e−50 (the
  minimizer *is* the target). All planted configs: **exactly `√2 = 1.41421`** — the
  maximum possible, i.e. the minimizer is *orthogonal* to the target.
- **M1′, distance to the ideal `(P*)`** truncated at degree `R`, for `R ≥ d`.
  Honest: ≤ 1.5e−50 for every `R ≥ d` — the ground state always lies in the ideal.
  Planted: **exactly 1.0**. (This one is T1 but downstream of (ii): once `λ_min ≠ 0` the
  ideal is a different eigenspace, so orthogonality is automatic.)
- **M2, the independent measure — `|ĉ_R(β)| / (‖c_R‖₂ (Σ_k |β|^{2k})^{1/2}) ∈ [0,1]`** at
  the off-circle planted `β`; it is 0 iff `β` is a zero of the minimizer, so function
  convergence requires M2 → 0. Measured to R = 36 it does the opposite:

  | config | R=2 | R=6 | R=10 | R=20 | R=36 | min over R≥4 |
  |---|---|---|---|---|---|---|
  | C1-g1real | 0.207 | 0.443 | 0.580 | 0.691 | **0.707** | 0.339 |
  | C2-quartet | 0.091 | 0.391 | 0.502 | 0.550 | **0.578** | 0.303 |
  | C4-mixed | 0.126 | 0.388 | 0.502 | 0.550 | **0.578** | 0.244 |

  C1 increases monotonically to `1/√2`; the quartets oscillate about ≈ 0.55. The
  **radial gap `| |β| − 1 | = 0.3` is constant in `R` and never closes** — it cannot, by
  H1. See `figures/stage2_iv_rescored.png`.

**F4. The sign of `λ_min` is the discriminator (T1).** It is the only one of the four
steps that separates C0 from C1–C4, it does so at tiny windows (R = 1, 2, 4), and
attack A4 found no off-line configuration anywhere in a scan over
ρ ∈ {1.0001, 1.01, 1.1, 1.3, 2, 5} × 6 angles that stays PSD out to R = 30.

**F5. Exact detection law (T1).** For the real off-line pair, `T_R = u w^T + w u^T` with
`u_i = ρ^i`, `w_i = ρ^{−i}`, so the two nonzero eigenvalues are `(u·w) ± |u||w|` and

```
λ_min(T_R) = (R+1) − |u||w|,      |u|² = (ρ^{2R+2}−1)/(ρ²−1), |w|² = (ρ^{−2R−2}−1)/(ρ^{−2}−1)
```

verified against direct eigensolve to relative 4.5e−80. Cauchy–Schwarz gives
`|u||w| ≥ u·w = R+1` with equality iff `ρ = 1`, so **`λ_min < 0` for every `R ≥ 1` the
instant `ρ ≠ 1`**. Expanding at `ρ = 1 + ε`:

```
λ_min = − C(R+2,3) · ε² + O(ε³)
```

with `C(R+2,3) = R(R+1)(R+2)/6` matching the fitted constant to 1e−6 for
R = 1,2,3,4,5,8,16,32 and fitted exponent 1.9999996. Asymptotically
`λ_min ~ −ρ^{R+2}/(ρ²−1)` — exponential divergence at rate ρ.

**F6. Window size is limited by precision, not by ε (T1 for the law; T3/toy for the
ζ-side numbers).** Because detection happens at
R = 1 for any ε > 0, the only obstruction is resolving `ε²` above the arithmetic noise
floor. At a floor of `10^{−17}`, R = 1 suffices down to ε ≈ 3e−9, and R = 18 is needed at
ε = 1e−10; in float64 (`10^{−16}`) nothing below ε ≈ 1e−8 is visible at small R.
**T3 / toy.** The sentence "a `λ_min` certified at the 1.7e−17 level excludes off-line
deviations only down to ε ≈ 4e−9" is a statement about **this rank-2 toy Toeplitz form
and nothing else**. It is *not* a sensitivity claim about Connes's or Zhu's computation:
the ζ window form has a different kernel, an archimedean term, a continuum of
frequencies and infinitely many zeros, and its ε↦λ_min constant has not been computed.
Quoting the 4e−9 figure as a bound on the ζ side would be unjustified. Stage 3e is where
the corresponding ζ-side constant gets measured; until then this is an analogy.
See `figures/stage2_detection_law.png`.

**F7. Masking by on-line zeros is weak — in the toy (T1 arithmetic, T3/toy transfer).**
Adding `k` honest on-line pairs around one off-line quartet pushes the detection window
from `R = 2` (k = 0) only to `R = 12` (k = 24) — a window of size 13 still certifies a
violation hidden among 26 on-line zeros, while the total rank is 52. Growth is sublinear
and steps in plateaus. The *arithmetic* is exact over Q; the **extrapolation is T3/toy**.
The masking angles are 24 generic rational cosines on a finite spectrum, nothing like the
density, spacing statistics or infinitude of the ζ zeros, and the ζ form is not of finite
rank at all. This experiment shows masking is weak **for finite generic spectra**; it is
not evidence about ζ.

**F8. float64 reproduces Zhu's warning here (T2).** At `R ≥ 2g`, where `λ_min = 0`
exactly, float64 `eigvalsh` returns spurious values of size `10^{−16}` to `10^{−15}`,
negative in most rows (e.g. E5 at R = 2: `−5.1e−17`; C0-online at R = 20: `−6.1e−15`).
Every sign conclusion in this report is from mpmath or exact rational arithmetic.

### Deliverable question

> In the CvS/Connes pipeline, which step carries the RH content — (i) real zeros of the
> approximant, (ii) sign of λ_min, (iii) simplicity/evenness, (iv) convergence of the
> approximant zeros to the true zeros?

**Answer in this setting (T3, resting on F1–F7).** The answer depends on how (iv) is
read, and that turns out to be the crux:

| step | reading | verdict |
|---|---|---|
| (i) real zeros of the approximant | — | **free** — Carathéodory–Fejér gives it from Toeplitz structure plus simplicity. True verbatim for planted off-line spectra, even at `λ_min = −642` (F1). This is *implicit in CvS's own statement*: their theorem is about an arbitrary real even distribution, with no arithmetic input, so RH-freeness of step (i) is built into the hypothesis class. The planted control makes it concrete and quantitative; it does not discover it. |
| (ii) sign of `λ_min` | — | **carries the RH content** (F4, F5). The only step separating on-line from off-line, at windows R = 1, 2, 4. No off-line config tested survives it. |
| (iii) simplicity / evenness | — | **free, and worse than free** — never fails under planting, while it fails routinely in the honest control for structural reasons (F2, F2a). |
| (iv) convergence to the true zeros | as **zero ordinates/angles** | **free** — converges to ~1e−3 for off-line spectra too (F3). |
| (iv) convergence to the true zeros | as **Hurwitz convergence of the normalized minimizer FUNCTION** | **carries RH content** (F3). Honest: the minimizer *is* the target exactly at `R = d`. Planted: orthogonal to the target (M1 = √2), M2 rises to ≈ 0.6–0.7 instead of falling to 0, radial gap constant. |

So the brief's expectation — "(i) free; RH sits in (ii) and (iv)" — is **correct once (iv)
is read at the function level**, which is what CvS actually assert. My earlier scoring of
(iv) as free was an artefact of measuring zero angles instead of functions; see
retraction 6.

The sharper statement the control yields is that (i) and (iv) are **not independent**:
H1 forces every approximant zero onto the critical circle, and that is exactly what makes
Hurwitz convergence to an off-line target impossible. The free witness is the obstruction.
What distinguishes the two is *which* object you track — the zero set (free) or the
function (not free) — and a numerical table of approximant-zero ordinates against true
zero ordinates, of the kind Connes reports for the first 50 zeta zeros, lives on the free
side of that line.

Two corollaries to carry into Stage 3 (both T3):

1. The RH content reachable from *geometric-side data alone* sits in a **sign** that
   scales as `ε²` in the toy, so the question is one of **certified numerical precision** —
   which is why Zhu's float64 warning is not a technicality but the centre of the problem.
2. The function-level reading of (iv) is the one to instrument on the ζ side. Reproducing
   Connes's ordinate table is reproducing the free half; the content is whether the
   normalized ground state converges **as a function**, which Stage 3b should measure
   directly rather than infer from zero locations.

---

## Stage 3 — number-field replication (zeta)

### Assembly (T1 derivation, `notes/stage3_assembly.md`)

Zhu's geometric side with `f` real supported in `[-L,L]`, `F = f-hat`:

```
Q(f) = 2 F(i/2) F(-i/2) + (1/2pi) \int_R |F(t)|^2 Psi_L(t) dt
Psi_L(t) = Re psi(1/4 + it/2) - log pi - sum_{log n < 2L} 2 Lambda(n) n^{-1/2} cos(t log n)
```

The brief's `2F(i/2)^2` is the even case of the pole term `h(i/2)+h(-i/2) = 2F(i/2)F(-i/2)`;
for odd `f`, `F(-i/2) = -F(i/2)` and **the pole term flips sign to `-2 P_j P_k`**. I need
the odd sector for 3c, so both are implemented.

Three bases, all with closed-form transforms, all complete in their sector:

| tag | basis | `F_k(t)` | `P_k = F_k(i/2)` | vanishes at `±L`? |
|---|---|---|---|---|
| `even_d` | `cos((2k+1)pi x/2L)` | `L[sinc((w-t)L)+sinc((w+t)L)]` | `(-1)^k 2w cosh(L/2)/(w^2+1/4)` | yes |
| `even` | `cos(k pi x/L)` | same | `(-1)^k sinh(L/2)/(w^2+1/4)` | no |
| `odd` | `sin(k pi x/L)` | `L[sinc((w-t)L)-sinc((w+t)L)]`, times `i` | `(-1)^k 2w sinh(L/2)/(w^2+1/4)` | yes |

**No oscillatory quadrature is used anywhere.** Everything moves to `u`-space by Parseval
with the cross-correlation `C_jk(u) = \int f_j(x) f_k(x-u) dx`:

```
(1/2pi)\int F_j conj(F_k) dt = C_jk(0),    (1/2pi)\int F_j conj(F_k) cos(tu) dt = C_jk^sym(u)
```

so the `-log pi` term is diagonal and the prime term is a **finite sum** of closed-form
values. `C_jk` is an explicit finite combination of `sin(alpha u + beta)` and
`(2L-u)cos(gamma u)`, verified against direct numerical convolution to 1e-26.

For the archimedean term, `psi(z) = -gamma + sum_m [1/(m+1) - 1/(m+z)]` at `z = 1/4+it/2`
gives the inverse Fourier transform `e^{-c_m|u|}`, `c_m = 2m+1/2`, hence

```
Arch_jk = -gamma C(0) + sum_{m>=0} [ C(0)/(m+1) - 2 \int_0^{2L} C_jk(u) e^{-c_m u} du ].
```

The `m`-sum is then done **analytically** (digamma + Lerch transcendent + Hurwitz zeta),
turning a double quadrature into a closed form — the difference between ~17 s and ~0.01 s
per matrix entry, which is what made the high-precision runs feasible. Checked three ways:
closed form vs `mp.nsum` of the same series (agree to 1e-30), and both against direct
numerical integration of the original `t`-integral truncated at `T`, where the discrepancy
falls as `log T / T^3` exactly as predicted (1.6e-7 → 3.2e-9 for `T` = 300 → 1200).

### 3a — GATE A: is the assembly right? (T2, independent of Zhu)

Under RH the explicit formula says `Q(f) = sum_rho |F(gamma_rho)|^2`. Comparing the
assembled geometric-side matrix with the zeros-side matrix built from `mpmath.zetazero`:

| sector | K=50 | K=200 | K=800 | predicted decay |
|---|---|---|---|---|
| `even_d` | 6.07e−5 | 3.58e−6 | **1.63e−7** | `gamma_K^{-3}` (f continuous) |
| `odd` | 7.28e−5 | 4.67e−6 | **2.26e−7** | `gamma_K^{-3}` |
| `even` | 2.18e−2 | 1.04e−2 | **4.23e−3** | `gamma_K^{-1}` (f jumps at ±L) |

Every discrepancy is the truncation tail of the zeros sum, decaying at exactly the rate
the basis's smoothness predicts (for `even_d` the measured discrepancy sat at ≈0.5× a
crude analytic tail bound at all five `K`). **Control: flipping the sign of the pole term
blows the discrepancy up to 5.35 (even) and 0.212 (odd)** — so the test has teeth, and it
is what validates my self-derived odd-sector pole sign.

### 3a — GATE B: `lambda*(0.8)` vs Zhu's certified interval (T2, 50 digits)

`lambda_min(N)` is a **variational upper bound** that must decrease to `lambda*(L)`.
Two independent complete even bases:

| N | `even_d` | `even` |
|---|---|---|
| 8 | 4.968e−15 | 4.005e−14 |
| 12 | 4.314e−17 | 1.386e−16 |
| 16 | 2.541e−17 | 3.199e−17 |
| 20 | 2.307e−17 | 2.753e−17 |
| 24 | **2.2702e−17** | 2.3357e−17 |

Monotone decreasing at every step in both bases (variational monotonicity holds), positive
at every `N` and both precisions. The best upper bound obtained is

```
lambda*(0.8) <= 2.2702e-17        (Zhu's certified interval: [8.9e-18, 2.27e-17])
```

**Gate B passes** in the sense that matters: an independently derived assembly, two
independent bases and independent arithmetic give a positive `lambda_min` that descends
monotonically onto Zhu's certified range. Two things I will *not* claim:

- The agreement with Zhu's upper end to four figures (`1.00009x`) is **largely a
  coincidence of where the sequence happens to sit at `N = 24`** — the sequence is still
  falling, so it will pass *through* 2.27e−17 and keep going. The honest content is the
  inequality, not the digits.
- I have **not** reproduced the quoted central value ≈1.66e−17. My result is consistent
  with the whole interval but has not converged tightly enough to locate `lambda*` within
  it, and I have no rigorous lower bound. Pinning it down needs a faster-converging basis
  (Legendre) or `N` well beyond 24.

### 3b — Connes's experiment, support [1,13] (T2, 110 digits)

`2L = log 13`, so the geometric side carries exactly the prime powers `<= 13`. At `N = 24`:

```
lambda_min = 3.654e-43,  simple,  even sector,  lam_2/lam_min ~ 1.3e6
```

Zeros of the ground state's Fourier transform against `mpmath.zetazero`:

| n | 1 | 2 | 3 | 5 | 7 | 9 | 10 | 11 | 12 | 15 | 25 | 50 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| abs err | **4.2e−39** | 1.5e−35 | 2.4e−33 | 1.8e−28 | 2.5e−23 | 2.4e−18 | 6.7e−17 | 1.0e−14 | 2.6e−11 | 0.197 | 0.603 | 0.192 |

The error profile has Connes's shape — astronomically accurate for the leading zeros,
degrading fast — but **the ground state stops tracking the zeros past `n ≈ 12`**: beyond
that the "nearest zero of `F`" is simply whichever of `F`'s 57 real zeros happens to be
closest, and its position is basis-determined, not precision-determined (those entries are
bit-identical between the `N=20`/50-digit and `N=24`/110-digit runs). Empirically the
ground state annihilates roughly `N/2` leading zeros, so Connes's 50 would need `N ~ 100`.
Raising precision alone does not help: going from (N=20, 50 digits) to (N=24, 110 digits)
improved `gamma_1` from 5.6e−35 to 4.2e−39 but moved nothing beyond `n = 12`.

### 3c — simplicity, evenness, and the even/odd gap (T2, 40 digits, N=12)

| L | `lambda_min` even | `lambda_2` even | `lambda_min` odd | odd − even | global ground state |
|---|---|---|---|---|---|
| 0.5 | 1.0796e−6 | 2.265e−2 | 2.4409e−4 | +2.430e−4 | **even** |
| 0.8 | 4.3141e−17 | 1.749e−11 | 4.9808e−14 | +4.976e−14 | **even** |
| 1.0 | 9.9817e−23 | 2.192e−17 | 1.4256e−20 | +1.416e−20 | **even** |
| log13/2 | 4.8493e−29 | 1.366e−23 | 1.4589e−26 | +1.459e−26 | **even** |

The CvS hypothesis holds at every `L` tested: `lambda_min` is **simple** with a relative
gap `lambda_2/lambda_min` of 2e4 to 4e5, and the global ground state is **even**, the odd
sector sitting 2 to 5 orders of magnitude above it. The even/odd gap *widens* in relative
terms as `L` grows (odd/even = 226, 1154, 1.4e5, 3.0e5).

### 3d — recorder split (T2)

`A = Pole + Arch - logpi`, `2M = Prime`, `S = A - 2M` (the Weil form), N=12:

| L | sector | inertia A | inertia M | inertia S |
|---|---|---|---|---|
| 0.5 | even / odd | (1,0,11) / (1,0,11) | (6,0,6) / (6,0,6) | (0,0,12) / (0,0,12) |
| 0.8 | even / odd | (1,0,11) / (1,0,11) | (6,0,6) / (7,0,5) | (0,0,12) / (0,0,12) |
| 1.0 | even / odd | (1,0,11) / (2,0,10) | (5,0,7) / (7,0,5) | (0,0,12) / (0,0,12) |
| log13/2 | even / odd | (2,0,10) / (2,0,10) | (5,0,7) / (6,0,6) | (0,0,12) / (0,0,12) |

**Block F holds on the zeta side**: `A` and `M` are both indefinite at every `L` and in
both sectors, while `S` is PSD throughout. And the split is **parameter-dependent exactly
as in the 389a1 appendix and as in Stage 1**: `A` carries 1 negative direction at `L=0.5`
but 2 by `L=log13/2`, and `M`'s inertia moves with both `L` and the sector, while `S`
stays `(0,0,N)`. The decomposition is bookkeeping; only `S` is invariant.

Worth recording: in Stage 1 the brief-convention `A` had inertia `(1,0,R)` for every curve,
and here it is `(1,0,N-1)` at small `L`. Same shape, independently.

### 3e — teeth on the zeta side (T2 — a DIAGNOSTIC, not a certificate)

Built from **zero data**: the first `K=400` true zeros plus one planted off-line quartet
`{rho, 1-rho, conj rho, 1-conj rho}`. With `rho = 1/2 + i z`, the quartet's ordinates are
`{+-gamma +- i delta}`, `delta = Re rho - 1/2` — *exactly* the Stage 2 off-line quartet,
with `e^delta` in the role of the radius and `L` in the role of the window `R`. Its
contribution is `4 Re[F(w)F(-w)]`, `w = gamma_* + i delta`. Planted at `gamma_* = 14.1347`:

| L | control (true zeros only) | δ=0.2 (Re ρ=0.7) | δ=0.1 | δ=0.05 | δ=0.02 |
|---|---|---|---|---|---|
| 0.4 | +2.09e−4 | +1.45e−4 | +1.94e−4 | +2.07e−4 | +2.10e−4 |
| **0.5** | +1.08e−6 | **−2.70e−5** | **−5.22e−6** | **−4.58e−7** | +8.32e−7 |
| **0.6** | +1.95e−9 | −1.47e−3 | −2.61e−4 | −2.81e−6 | **−3.38e−7** |
| 1.0 | +1.53e−24 | −3.64e−2 | −8.82e−3 | −2.13e−3 | −3.31e−4 |
| 2.0 | +5.11e−56 | −4.32e−1 | −1.04e−1 | −2.52e−2 | −3.79e−3 |

The control is a sum of rank-1 PSD terms, so it is **exactly PSD by construction** — I use
it as the measured noise floor, and at 80 digits it stays positive (down to 5e−56) at every
`L`, so no sign below is roundoff. (At 40 digits it went spuriously negative past `L≈3.5`;
that run was discarded.)

**Three findings, all echoing Stage 2:**

1. **At a fixed height, the window needed is almost independent of the DEPTH** — planted
   at `gamma_* = 14.13`, the form turns indefinite at `L=0.5` for `delta` = 0.2, 0.1, 0.05
   and at `L=0.6` for `delta = 0.02`. A 10× shallower zero costs essentially nothing in
   window size, only in signal size. That is the zeta-side version of Stage 2's F6.

   **But the window needed depends strongly on the HEIGHT, and that was not anticipated by
   the function-field toy** (which has no height — its spectrum is a finite set of angles
   on one circle). Sweeping the planted height at fixed `delta = 0.1`
   (`data/stage3e_height_sweep.csv`, K=300, 70 digits):

   | planted height `gamma_*` | 3.0 | 14.13 (`=gamma_1`) | 18.0 | 50.0 (`~gamma_10`) | 100.0 (`~gamma_29`) | 200.0 (`~gamma_79`) |
   |---|---|---|---|---|---|---|
   | first `L` with `lambda_min < 0` | 0.6 | **0.5** | 0.8 | **2.5** | not by `L=3` | not by `L=3` |

   A window of half-width `L` is band-limited, so it simply cannot see an off-line zero
   high in the critical strip: at `gamma_* = 100` and `200` the planted form is
   indistinguishable from the control out to `L = 3` (at `L = 3` both sit at the 1e−73
   noise floor, so "not detected" there means *unresolved*, not *certainly positive*).
   **Short windows are a probe of the low-lying zeros only.**
2. **`lambda_min ~ -C(L) delta^2`.** Fitted exponents 2.098, 2.033, 2.068, 2.050, 2.069 at
   `L` = 0.8, 1.0, 1.3, 1.6, 2.0 — the same quadratic law as the function-field toy.
3. **`C(L)` grows like `L^3`**: measured `C` = 0.310, 0.827, 2.157, 4.628, 9.472 at
   `L` = 0.8, 1.0, 1.3, 1.6, 2.0, against `L^3` = 0.512, 1, 2.197, 4.096, 8 — the same
   cubic-in-window law as Stage 2's exact `binom(R+2,3) ~ R^3/6`, with a constant about
   4–7× larger.

**This retires the T3/toy caveat on F6 with a measured zeta-side number.** Combining
`lambda*(0.8) <= 2.27e-17` with `C(0.8) = 0.310`:

> an off-line zero **at height `gamma ~ 14.13`** of depth `delta` contributes about
> `-0.31 delta^2` to the `L=0.8` window form, so a certified `lambda*(0.8) > 0` at the
> `2.27e-17` level is consistent with such a zero only if `delta < 8.6e-9`
> (and `delta < 5.2e-9` at `L=1.0`).

Stage 2's toy estimate was `eps ~ 4e-9`; the measured zeta-side value is `8.6e-9` — the toy
was right to within a factor of about two **for a low-lying zero**.

**The height restriction is not a footnote, it is the main limitation.** By finding 1, the
same `L = 0.8` window has essentially *no* sensitivity to an off-line zero at height 50 or
above at any depth — the corresponding `C(L)` is numerically zero there. So the statement
above constrains `delta` only for off-line zeros among the first few; it says nothing
whatever about the rest of the critical strip. And it remains a diagnostic built from zero
data with `K = 400`, not a theorem about zeta.

### 3f — Davenport–Heilbronn: NOT ATTEMPTED (difficulty flagged, as the brief asks)

`-f'/f` for the Davenport–Heilbronn function has no Euler product, so there is no
prime-power "geometric side" of the form used throughout Stage 3; its Dirichlet series has
a different abscissa of convergence, and the relevant explicit formula would have to be
derived and numerically validated from scratch (the analogue of Gate A) before any
eigenvalue computed from it would mean anything. Doing that properly is a project of its
own, and doing it improperly would produce numbers that look like the Stage 3 tables but
certify nothing. Not attempted.

---

## Stage 4 — literature check (done BEFORE any novelty claim)

Searched arXiv and the surrounding literature for prior work on each of H1–H4, the
Stage 2 diagnosis, and the Stage 3 replications. Findings, stated as what I verified:

**H1 and H2 are NOT novel, and I should have said so from the start.** The CvS paper
(arXiv:2511.23257, *Quadratic Forms, Real Zeros and Echoes of the Spectral Action*) is
organised as a five-step proof whose **Step 1 is explicitly "a C*-algebraic proof of a
corollary of Carathéodory–Fejér's 1911 structure theorem for Toeplitz matrices"** and whose
**Step 5 is Hurwitz's theorem on zeros of uniform limits of holomorphic functions**. So:

- H1 *is* CvS Step 1. My contribution is only to have measured that it survives planted
  off-line data, which is the "witness is free" point — not the theorem itself.
- H2 (exact recovery of the atoms from a rank-deficient PSD Toeplitz moment matrix) is the
  classical Carathéodory–Fejér / **Pisarenko harmonic decomposition** fact, standard in
  signal processing since 1973. The brief anticipated this; it is correct.
- The Hurwitz step confirms, from the source, that reading (iv) as function-level
  convergence (not zero-angle convergence) is the right reading — which is what the
  Stage 2 correction already concluded.

**The Stage 2 signature count appears to be known.** The statement that *"the negative
index of finite truncations of Weil's form equals the number of off-line zero pairs seen by
the truncation"* surfaced in the search as an existing claim in this literature. That is
precisely my Stage 2 prediction (off-line quartet → signature `(2,·,2)`, off-line real pair
→ `(1,·,1)`). I claim no novelty for it; what I add is the exact-over-Q certification and
the first-detection-window measurements.

**Stage 3b is an independent small-scale reproduction of published work.** Groskin,
arXiv:2605.20224 (*High-Precision Approximation of Riemann Zeros via the Truncated Weil
Form*) implements the CvS Galerkin matrix at cutoffs `c = 13 … 67` and `c = 100`. At
`c = 13`, `N = 100` the reported first-zero error is `~2e-55`; at `c = 100`, `N = 250` the
smallest even-sector eigenvalue reaches `~1e-334` and recovers the first ten zeros to
307–329 digits. My `N = 24` run (first-zero error `4.2e-39`, `lambda_min = 3.65e-43`) sits
on the same trajectory at much smaller `N`, and that paper's scale independently supports
my explanation that the depth of the error profile is basis-size-limited.

Also relevant and consulted: Groskin arXiv:2607.02828 (*A finite Guinand–Weil dictionary
and archimedean tail order*), which treats the archimedean tail I handle in closed form;
Suzuki arXiv:2606.09096 (*Weil's quadratic form via the screw function*); Zhu
arXiv:2608.24827, the source of the certified `lambda*(0.8)` interval; and CCM
arXiv:2511.22755 (*Zeta Spectral Triples*).

**What I did not find.** No function-field (curve over a finite field) control experiment
for the CvS/CCM pipeline, and no planted-off-line-configuration test of it, turned up in
these searches. **This is weak evidence.** Absence of search hits is not absence of prior
work, I did not read these papers in full, and the area is moving fast. I make **no
novelty claim** for Stages 1–2; the most I will say is that I did not find this particular
control experiment already done, and anyone building on it should check properly.

Sources consulted: [CvS 2511.23257](https://arxiv.org/abs/2511.23257),
[CCM 2511.22755](https://arxiv.org/pdf/2511.22755),
[Zhu 2608.24827](https://arxiv.org/pdf/2608.24827),
[Groskin 2605.20224](https://arxiv.org/abs/2605.20224),
[Groskin 2607.02828](https://arxiv.org/abs/2607.02828),
[Suzuki 2606.09096](https://arxiv.org/pdf/2606.09096).

---

# Round 2 — tightened Stage 3, settled novelty

## Task 1 — the 3e certification asymmetry

**The user's objection was right, and the earlier result was wrong.** A negative Rayleigh
quotient from any trial function certifies indefiniteness; a finite basis that finds none
certifies nothing. The Round-1 entries "not detected by L = 3" at `gamma_* = 100, 200`
were upper estimates of the onset reported as if they were limits. Working through it
turned up **two** errors, not one.

### Error A — the on-line side must come from the geometric assembly

Near the onset the negative Rayleigh quotient is ~1e-27, while truncating the zeros sum at
`K = 400` omits positive mass ~1e-6. No reachable `K` closes that gap (the tail falls like
`log(gamma_K)/gamma_K^3`; one would need `gamma_K ~ 1e11`). So the truncated assembly
cannot certify anything near the onset, and because it omits *positive* mass it is biased
toward false detections. Round 2 takes the on-line contribution from the **geometric
side**, which has no truncation error.

### Error B — I had conflated two different counterfactuals

For real even `f`, `F(conj z) = conj F(z)`, so with `w = gamma + i delta`

```
4 Re F(w)^2 = 4 F(gamma)^2 - 4 delta^2 (F'(gamma)^2 + F(gamma) F''(gamma)) + O(delta^4).
```

Simply **adding** a quartet therefore contributes `+4 F(gamma)^2 > 0` unless `F` nearly
vanishes at `gamma`. Confirmed directly: at `gamma_* = 100`, `L = 1.38`, the added-quartet
contribution at the ground state is `+1.8e-47`, not negative.

**This does not make the added-quartet experiment illegitimate.** Adding a zero and moving
a zero are two different, equally meaningful counterfactuals, and **Bombieri (2000) §13
studies the first**: his `rho_0 = 0.52 + 3.14i` plants at an ordinate `3.14` which is not a
zeta zero, and he still finds a critical window `t_c`. The honest reading of the Round-1
rows at heights 18, 50, 100, 200 is therefore *not* that they were meaningless, but that
they were subject to **the same certification asymmetry as everything else in Task 1**:
a finite basis and a finite window that find no negative direction certify nothing, and the
truncated on-line side made the positive entries unreliable in the other direction too.

What this round does is study the **second** counterfactual, which is the one matched to the
functional equation: remove the on-line pair at `+-gamma_n` and insert the quartet in its
place:

```
Delta(f) = -2 F(gamma_n)^2 + 4 Re F(gamma_n + i delta)^2
         =  2 F(gamma_n)^2 - 4 delta^2 (F'^2 + F F'') + O(delta^4),
```

negative exactly when the trial function nearly annihilates `gamma_n` — which is what the
window minimiser does for the leading zeros. So Round 2 plants **at actual zeta zero
ordinates** and sweeps `n`.

### (a) Does the onset move earlier as the basis grows? YES — decisively

Basis = low block `{0..N0-1}` plus a resonant block `{k : w_k ~ gamma_n}`; any subspace
gives a valid one-sided certificate, and the resonant block makes high heights reachable
without paying for every intermediate mode. `B1 = (8,0)`, `B3 = (24,3)`, `B5 = (44,7)`.

| plant at | `gamma_n` | `L_pred` | onset, B1 (dim 9) | onset, B3 (dim ~31) | onset, B5 (dim ~52) |
|---|---|---|---|---|---|
| `gamma_1` | 14.13 | 0.405 | 1.3× | 1.3× | 1.3× |
| `gamma_5` | 32.94 | 0.828 | 1.0× | 1.0× | 1.0× |
| `gamma_10` | 49.77 | 1.035 | never ≤2.2× | 1.0× | 1.0× |
| `gamma_30` | 101.32 | 1.390 | never ≤1.7× | 1.3× | **0.85×** |
| `gamma_80` | 201.26 | 1.733 | never | not at ≤1.15× | **not at ≤1.0×, dim 59** |

**`gamma_30 ~ 101`, which Round 1 reported as "not detected by L = 3", is detected at
`L = 1.18` — that is 0.85 x L_pred, and about 2.5x earlier than the window Round 1 gave
up at.** The limitation was basis size, not the window.

### (b) Onset vs the bandwidth prediction

`L_pred(gamma_*) = (1/2) log(gamma_*/2pi)` is the `L` at which Zhu's window bandwidth
`T*(L) = 2 pi e^{2L}` first reaches `gamma_*` — and Zhu's own Landau–Widom decay law is
stated in terms of that same `T*`, so this is a sourced scale, not a fitted one. Measured
onsets at the largest basis run: **1.3×, 1.0×, 1.0×, 0.85×** for `n` = 1, 5, 10, 30. The
onset tracks `L_pred` to within about ±30% and drifts *below* it as the basis grows.

### (c) Interval certification — partially achieved, and I will not overstate it

I implemented mpmath-interval certification (`src/iv_certify.py`), including working around
a bug in mpmath (`iv.sincpi` calls a non-existent `ctx.sinpi`) by using the exact closed
form `F_k(t) = (-1)^k 2 w_k cos(tL)/(w_k^2 - t^2)`. It certifies the **zeros-side**
assembly fine. But Error A means the zeros-side assembly is the wrong object here, and
mpmath's interval context has **no digamma and no Lerch transcendent**, so the geometric
side cannot be interval-evaluated with this library. What is certified instead:

- **Exact identity (T1):** the LDL^T witness satisfies `v^T Z v = D_j` identically —
  verified to 1e-30 against direct evaluation on 60 random symmetric matrices, along with
  the inertia of every leading block against `eigsy`.
- **Measured noise floor (T2):** the geometric control is PSD by construction, so whenever
  it comes out *negative* the magnitude is the arithmetic noise floor. Detections are
  reported as trustworthy only if they exceed `100 x` that floor. This filter is in the
  code and in the CSV (`trustworthy` column).

**Full interval certification of the geometric-side Rayleigh quotient was NOT achieved**
and would need an interval implementation of digamma/Lerch.

### (d) Retraction 10, reworded; and the delta-dependence

See retraction 11 below for the replacement text. On `delta`: the Round-2 sweep ran at
`delta = 0.1`; the Round-1 sweep covered `delta` = 0.2, 0.1, 0.05, 0.02 and moved the
onset only from `L = 0.5` to `L = 0.6` across a 10x change in depth. There is a reason:
detection is `C(L) delta^2 > lambda*(L)`, and by Zhu's Landau–Widom law `lambda*(L)` falls
super-exponentially in `L` while `C(L) ~ L^3` grows polynomially, so a factor in `delta^2`
shifts the crossing only logarithmically. **The onset is essentially delta-independent,
and that is a consequence of the decay law, not a coincidence.**

### (e) C(L) against the planted height

`C(L, gamma_*) = -lambda_min/delta^2` at fixed basis (24 low + 13 resonant), `delta = 0.02`:

| | `gamma_1`=14.1 | `gamma_5`=32.9 | `gamma_10`=49.8 | `gamma_20`=77.1 | `gamma_30`=101.3 | `gamma_50`=143.1 |
|---|---|---|---|---|---|---|
| **L=1.0** (`T*`=46.4) | 0.862 | 0.0137 | 2.9e−12 | none | none | none |
| **L=1.4** (`T*`=103.3) | 3.024 | 1.792 | 0.478 | 3.0e−13 | none | none |
| `gamma/T*` at L=1.4 | 0.14 | 0.32 | 0.48 | 0.75 | 0.98 | 1.39 |

`C` collapses by ~12 orders of magnitude as `gamma_*/T*` rises from 0.3 to 0.75, and above
`gamma_*/T* ~ 1` there is no detection at all at this basis. **But part of that fall-off
is basis size, not bandwidth**: `gamma_30` at `L = 1.4` shows "none" here at dim 37, yet
*is* detected at `L = 1.39` with dim 52. The honest statement is that detection needs
**both**

1. `L >~ L_pred(gamma_n)` — the bandwidth must reach the height; and
2. `N >~ 2n` — the basis must be able to annihilate the first `n` zeros,

and (2) *(superseded — see retraction 17)*. The `N ~ 2n` reading of (2) was imported from
3b, where it governs how many leading zeros a **contiguous** basis can annihilate. That is
not the detection criterion: detection only needs some `f` in the subspace with
`Q_geom(f) + 2F(gamma)^2 < 4 delta^2 (F'^2 + F F'')`. In the data `gamma_30` is detected at
dimension 46 (top index 45, against `2n = 60`) and `gamma_80` at dimension 59 (top index
~134, against `2n = 160`). The correct second condition is that the basis must **contain
modes resonant with** `gamma_n`, i.e. indices near `k* ~ gamma_n L/pi`, which the resonant
block supplies at a small fraction of the dimension a contiguous basis would need.

The `gamma_80` rows were re-run at a flat 220 digits so that non-detection could not be
blamed on precision: at `L = 1.47` (0.85x) and `L = 1.73` (1.0x) the geometric control
stays **positive** (1.47e-69, 2.01e-57), so the arithmetic is sound and the absence of a
negative direction is a genuine statement about the dim-59 subspace — and about nothing
larger. `gamma_30` at 220 digits likewise gives a clean margin: at `L = 1.60` the control
is `+1.36e-89` against a detection of `-3.41e-80`, a factor of `2.5e9`.

## Task 2 — the delta^2 coefficient, derived and checked

**(a) Expansion and normalisation (T1, verified T2).** Derivation above. Normalisation
against the code's zeros-side convention (the one Gate A validated): an on-line pair
`{+-gamma}` contributes `2 F(gamma)^2`; an off-line quartet contributes `4 Re F(w)^2`. The
brief's `2 Re F(gamma - i delta)^2` is the **pair-normalised half** of the quartet, and
`F(conj w) = conj F(w)` makes the two signs of `delta` agree. Checked numerically against
the exact `4 Re F(w)^2` at the window minimiser: relative error 3.6e−7 to 9.9e−7 at
`delta = 1e-3`, consistent with the `O(delta^4)` remainder.

**(b) Does the fitted `C(L)` match the predicted coefficient? Only at the right minimiser.**

| L | `C` measured from `lambda_min` | `4(F'^2+FF'')` at the **perturbed** min | same at the **unperturbed** ground state |
|---|---|---|---|
| 0.8 | 0.3108 | 0.3186 | 1.22e−4 |
| 1.0 | 0.8530 | 0.8668 | 5.64e−5 |
| 1.3 | 2.4987 | 2.6405 | 1.87e−5 |
| 1.6 | 4.9218 | 4.9356 | 8.38e−7 |
| 2.0 | 10.083 | 10.199 | 6.21e−5 |

Agreement with the perturbed minimiser is 2–6%; the unperturbed ground state is off by
**3 to 6 orders of magnitude**. The reason is structural: `Q_geom` has an enormous
near-null space (`lambda_2/lambda_1 ~ 1e5–1e6`), so the perturbed minimiser moves far from
the ground state to maximise `F'^2` at almost no cost in `Q_geom`. Evaluating the
perturbation at the *unperturbed* minimiser gives only a weak lower bound on `C`.

**(c) Not the Laguerre expression.** The coefficient here is `F'^2 + F F''`; the Laguerre
quantity is `L_1 = F'^2 - F F''`, opposite in sign on `F F''`. At the window minimiser
`F(gamma_1)` is 3.6e−16 down to 2.2e−53, so the two agree **to every printed digit** in the
table above — but that is a property of the minimiser, not an identity. **Any link to the
Partition note's Remark 1 is leading-order only** and fails as soon as `F(gamma)` is not
negligible — which is exactly the regime (planting away from a zero) where the sign of the
perturbation flips, as Error B above shows.

## Task 3 — the L^3 law as a Cauchy–Schwarz bound

**(a) Derivation (T1).** For real even `f` supported in `[-L,L]` with `||f||_2 = 1`,
`F'(gamma) = -int u f(u) sin(gamma u) du`, so

```
|F'(gamma)|^2 <= int_{-L}^{L} u^2 sin^2(gamma u) du  ->  L^3/3   as gamma -> infinity,
```

hence `C(L) = 4|F'|^2 <= 4L^3/3`. (Dropping the `sin` weight gives `2L^3/3` for `|F'|^2`,
i.e. a bound twice as weak.)

**(b) Measured ratio — the bound is saturated.**

| L | 0.6 | 0.8 | 1.0 | 1.3 | 1.6 | 2.0 |
|---|---|---|---|---|---|---|
| `C` measured | 3.4e−4 | 0.3108 | 0.8530 | 2.4987 | 4.9218 | 10.083 |
| `C / (4L^3/3)` | 0.0012 | 0.455 | 0.640 | 0.853 | 0.901 | **0.945** |
| `C / L^3` | 0.0016 | 0.607 | 0.853 | 1.137 | 1.202 | **1.260** |

The ratio increases monotonically toward 1 and `C/L^3` toward `4/3 = 1.333`: **the window
minimiser asymptotically saturates the Cauchy–Schwarz bound.** (`L = 0.6` is anomalous
because `lambda*(0.6)` is still large there, so the minimiser cannot move freely.)

**(c) Function-field analogue, and the two laws are the same bound.** For the off-line real
pair, expanding the Stage 2 exact formula `lambda_min = (R+1) - |u||w|` gives

```
lambda_min = -2 eps^2 [ S2 - S1^2/(R+1) ] + O(eps^3),  S1 = sum i, S2 = sum i^2,
S2 - S1^2/(R+1) = sum_{i=0}^{R} (i - R/2)^2 = R(R+1)(R+2)/12,
```

so `lambda_min = -C(R+2,3) eps^2`, reproducing the Stage 2 law **exactly** and exhibiting
it as `2 eps^2 x` the **centred** second moment of the window. The centring is not a
convention — it falls out of the optimisation, which removes the mean. The naive uncentred
bound `sum k^2 = R(R+1)(2R+1)/6` overshoots by `(2R+1)/(R+2) -> 2` (measured: 1.00, 1.25,
1.50, 1.70, 1.83, 1.91, 1.95 at `R` = 1..64).

**Same law in both settings.** Identifying the window *length* (`R+1` lattice points ↔ an
interval of length `2L`, i.e. `R <-> 2L`) and the pair-vs-quartet factor 2:

```
function field :  C_R  = R^3/6 + O(R^2)
zeta           :  C(L) = 4L^3/3 = (2L)^3/6 .
```

They are literally the same expression, `(window length)^3 / 6`.

## Task 4 — Stage 3 loose ends

**(a) 3c coverage.** The even/odd study covered exactly `L in {0.5, 0.8, 1.0, log(13)/2 =
1.2824746787}`, `N = 12`, 40 digits. The `L = 0.8` entry
(`lambda_min^even = 4.3141e-17` at `N=12`, improving to `2.2702e-17` at `N=24`) is a
variational upper bound **consistent with Zhu's certified `[8.9e-18, 2.27e-17]`**; Zhu also
certifies the ground state there is simple and even, which my computation reproduces. The
other three `L` values carry **no external certificate and are T2 only.**

**(b) 3d inertia as N grows — the stated expectation is REFUTED.**

| L | basis | N=6 | N=10 | N=14 | N=18 | N=24 |
|---|---|---|---|---|---|---|
| 0.5 | both | (1,0,5) | (1,0,9) | (1,0,13) | (1,0,17) | (1,0,23) |
| 0.8 | both | (1,0,5) | (1,0,9) | (1,0,13) | (1,0,17) | (1,0,23) |
| log13/2 | both | (2,0,4) | (2,0,8) | (2,0,12) | (2,0,16) | (2,0,22) |

`A` does gain positive directions as higher frequencies enter — but **the number of
negative directions does not grow**: it is 1 for `L <= 0.8` and 2 at `L = log(13)/2`, and
it is **identical in both even bases**. So the expectation that `(1,0,N-1)` would be
basis-dependent is not borne out: the negative index of `A` is a function of `L` alone,
stable in `N` and in the choice of basis, jumping from 1 to 2 somewhere in
`0.8 < L < 1.28`. `M` stays near-balanced, `(~N/2, 0, ~N/2)`.

*Bookkeeping convention, stated as asked.* On the zeta side I group the `-log pi` term with
`A`, so `A = Pole + Arch - log(pi)`, `2M = Prime`. There is no `n = 0` ambiguity here (no
discrete index); the zeta-side analogue of the Stage 1 convention choice is precisely where
`-log pi` is placed. With the Stage 1 **brief convention** (`A = 2g I + pole`) the
function-field `A` had inertia `(1, 0, R)` — the same shape as the zeta-side `(1, 0, N-1)`
at small `L`. Under the Stage 1 **uniform convention** the function-field `A` was rank 2,
`(1, R-1, 1)`; any comparison must name which convention is in use.

*Caveat (T2):* at `L = log(13)/2`, `N >= 14`, `S` is reported with nullity 1–3 at 40
digits. That is precision exhaustion (`lambda*` there is ~1e-29 and falling with `N`), not
true nullity.

**(c) Larger-N Connes window — done.** At `L = log(13)/2`, `N = 36`, 140 digits:

```
lambda_min = 8.977e-52   (was 3.654e-43 at N=24, 110 digits),  simple, even
```

and the ground state now tracks far more zeros — a direct look at step (iv) on zeta:

| n | 1 | 5 | 10 | 12 | 15 | 20 |
|---|---|---|---|---|---|---|
| abs err, N=24 | 4.2e−39 | 1.8e−28 | 6.7e−17 | 2.6e−11 | 0.197 | 0.0186 |
| abs err, N=36 | **7.2e−48** | 4.2e−38 | 8.7e−29 | 1.5e−25 | **9.8e−21** | **2.2e−13** |

At `N = 24` the ground state stopped tracking past `n ~ 12`; at `N = 36` it tracks past
`n ~ 20`. That is the `N/2` rule again, and it confirms that the Round-1 cutoff was basis
size and not precision. Both `lambda_min` and the per-zero errors are moving in the
direction of Groskin's `N = 100` figures (first-zero error ~2e−55 at this cutoff).

## Task 6 — explicit novelty check

Checked against CvS 2511.23257, CCM 2511.22755, Groskin 2605.20224, Groskin 2607.02828,
Suzuki 2606.09096, Zhu 2608.24827, Connes 2602.04022.

> **HOW THIS WAS CHECKED, AND ITS LIMIT.** Each paper was queried through an automated
> fetch of its arXiv landing page, i.e. **abstract and metadata level, not a full read**.
> Several are 30–40 pages. One fetch said so explicitly ("to answer these questions
> accurately would require examining the full paper text (33 pages)"). So every
> "not found" below means *not found at abstract level*, which is weak evidence. Nothing
> here should be treated as a literature clearance; anyone building on this must read the
> papers properly.

| # | Item | Verdict | Reference / note |
|---|---|---|---|
| 1 | Function-field control of the CvS/Connes pipeline with planted spectra | **not found** (all 7) | No paper mentions curves over finite fields or Frobenius angles. CCM is spectral-triple-theoretic; CvS is stated for arbitrary real even distributions. |
| 2 | Function-level re-scoring of step (iv) + the Hurwitz obstruction argument | **partially found** | The *Hurwitz step itself is CvS Step 5* ("a classical theorem of Hurwitz concerning the zeros of uniform limits of holomorphic functions"), and CCM state spectra converge to the zeros. The **obstruction direction** — that all-real approximant zeros *forbid* convergence to an off-line target — was not found. |
| 3 | `ker T_R = {P* q}`, CvS hypothesis failing for all `R > 2g` | **partially found** | CvS Steps 1–2 treat Toeplitz/convolution kernels and the eigenfunction–polynomial link; the ideal structure and the simplicity breakdown past the window degree were not found. |
| 4 | Toeplitz–Rosati bridge (`T_R` as a Rosati Gram matrix) | **not found** | No paper mentions Rosati involutions or `H^1` of a curve. *But see Task 5a: the underlying positivity criterion is classical linear algebra and is why Weil's proof works; only the identification with the window Toeplitz matrix is not found.* |
| 5 | `delta^2` sensitivity law with `C ~ L^3`, both sides | **not found** | Groskin 2605.20224 fits an empirical cutoff law `|log10 lambda_min| ~ 13.24 c^0.634` and then reports it falsified at larger `N`; that is decay in the cutoff, not sensitivity in off-line depth. **Zhu's Landau–Widom law** `-ln lambda*(L) ~ 2 pi^2 N(T*)/ln N(T*)`, `T* = 2 pi e^{2L}`, is adjacent and is the source of the bandwidth scale used in Task 1 — but it too is a decay law in `L`, not a `delta`-sensitivity law. |
| 6 | Block F: `A` and `M` indefinite on zeta windows, `S` PSD | **not found** | Groskin 2607.02828 proves the omitted archimedean tail is "a totally positive Cauchy–Stieltjes increment" with budget `B_T ~ (2N+1) rho log T/(pi^2 T)` — related in spirit, but it does not split the form into `A` and `M` and report inertia. |
| 7 | Height dependence of the detection onset | **not found** | No paper was found to discuss detectability as a function of position in the critical strip. |

**What may therefore be described as new (subject to the caveat above):** items 1, 4, 5, 7,
and the obstruction half of item 2. **Everything else is confirmation or reformulation:**
H1 is CvS Step 1 (Carathéodory–Fejér 1911) and Pisarenko harmonic decomposition; H2 is the
classical exact-recovery fact for rank-deficient PSD Toeplitz moment matrices; the Hurwitz
step is CvS Step 5; the Rosati positivity criterion is classical; Stage 3b reproduces
Groskin 2605.20224 at smaller scale; and the Stage 2 signature count (negative index of a
truncation = number of off-line pairs it sees) surfaced in the Round-1 search as an
existing claim.

**Two external numbers now confirmed against source**, which matter for Gate B:
Zhu normalises by `Q(f)/||f||^2` — **the same normalisation I used**, so my upper bound is
on the same quantity — and certifies `8.9e-18 <= lambda*(0.8) <= 2.27e-17`. Connes reports
the first 50 zeros to between `2.6e-55` and `1e-3` using primes below 13, which is the
profile shape Stage 3b reproduces at smaller `N`.

---

# Round 3 — prior art, inertia reconciliation, proofs (Tasks 1–3)

## Task 1 — prior art, verified against full texts

Full findings in `PRIOR_ART.md`. **The prior art is materially stronger than the Round-2
abstract-level check found, and it removes three claimed contributions.**

- **Hallouin–Perret, Trans. AMS 372 (2019), 5409–5451** (full text obtained from the
  authors' copy of the AMS offprint; not on arXiv under this title). Their
  **Proposition 5** (§1.2, pp. 5415–5416) gives `Gram(γ^0,…,γ^n)` as a Toeplitz matrix with
  `2g` on the diagonal and `x_n = ((q^n+1) − |X(F_{q^n})|)/q^{n/2}` off it — **my `T_R`,
  entry for entry**. Their **Theorem 6** (§1.3, pp. 5417–5418) states that `d`, the rank of
  the Frobenius space, is "the minimal integer such that `Gram(γ^0,…,γ^d)` is singular", and
  identifies a kernel generator — **my H2, restatement included**. Their **Theorem 36(ii)**
  (Appendix A.2, p. 5448) is **my H1** verbatim; **Lemma 33** is the leading-minor rank rule
  I used in Round 2; **Lemma 34** is the parity split.
  *Citation note:* the authors' own later bibliography prints the volume as "312"; the
  published header and CvS's reference [12] both give **372**. 372 is correct.
- **CvS 2511.23257, Introduction**: cites Hallouin–Perret as [12] for the function-field
  connection, and states "The key difficulty in this context, then, becomes the verification
  that zero is indeed the (simple) minimal eigenvalue of `T`." Their **Remark 2.3** states my
  Round-2 §3 kernel result in words.
- **Bombieri (2000)**, Rend. Lincei (9) Mat. Appl. 11, 183–233 — **obtained and read**, not
  "not checked". **Theorem 8** (§8) gives the negative-eigenvalue count = number of off-line
  conjugate pairs. And **§13 already runs the planted experiment**: "a fictitious zero `ρ_0`
  off the critical line, together with their images by complex conjugation and reflection" —
  my quartet — at `ρ_0 = 0.52 + 3.14i`, window `[−t,t]`, `N` up to 160, with a critical value
  `t_c` beyond which the negative eigenvalue stays bounded away from 0. **No `δ`-expansion**:
  the depth is fixed, never varied.

## Task 2 — the inertia of `A`, reconciled; and Block F withdrawn on ζ

**There was no contradiction between Rounds 1 and 2.** Round 1's columns headed "even / odd"
were `(EVEN_D, ODD)`; Round 2's "both bases" were `(EVEN_D, EVEN)` — two *even* bases. Both
give negative index 1 at `L ≤ 0.8` and 2 at `L = log13/2` in the even sector. Verified by
running all three bases explicitly (`src/task2_inertia_reconcile.py`,
`data/stage3d_inertia_reconciled.csv`).

**But the exercise broke the Block F claim.** Two bookkeeping conventions, both giving the
same `S = A − 2M`:

```
I  (used throughout):  A = Pole + Arch − log(π)·I,   2M = Prime
II (alternative)     :  A = Pole + Arch,              2M = Prime + log(π)·I
```

| L | conv I: A / M | conv II: A / M |
|---|---|---|
| 0.5 | (1,0,11) / (6,0,6) | **(0,0,12) / (0,0,12)** |
| 0.8 | (1,0,11) / (6,0,6) | **(0,0,12) / (0,0,12)** |
| log13/2 | (2,0,10) / (5,0,7) | (0,0,12) / (2,0,10) |

identical in all three bases. **Under Convention II both `A` and `M` are positive definite at
`L ≤ 0.8`.** So "A and M are both indefinite" on the ζ side is an artefact of where the
`−log π` term is placed, not a property of the form.

**Why the function-field case is different, and robust.** There `A` contains the pole part
`q^{|n|/2} + q^{−|n|/2} = u_i w_j + w_i u_j`, an exact rank-2 **hyperbolic** plane whose scale
grows like `q^{R/2}` against a fixed `2g·I`; so `A` is *eventually* indefinite under both
conventions tested (brief: from `R = 2`–3; uniform: from `R = 1`). On the ζ side the pole
term for even `f` is `2F(i/2)² ≥ 0` — **rank one and definite, not hyperbolic** (and rank-one
negative in the odd sector). There is no structural hyperbolic plane to find.

**Final convention, stated once:** Convention I, `A = Pole + Arch − log(π)·I`, `2M = Prime`.
All inertia tables in this report use it. The Block F claim is retained **only** for the
function field, where it is structural, and withdrawn on the ζ side.

## Task 3 — proofs

Full proofs in `notes/proofs.md`; every step numerically checked by `src/proofs_check.py`.

**(a) Function field — now T1.** For an off-circle real pair `{ρ, ρ^{-1}}`, `a = log ρ`:
`T_R = uw^T + wu^T` with `u_i = ρ^i`, `w_i = ρ^{-i}` (residual ≤ 6.2e−61); the nonzero
eigenvalues are `(R+1) ± |u||w|` and all others vanish (≤ 1.0e−59); so
`λ_min = (R+1) − |u||w| < 0` for every `R ≥ 1` whenever `ρ ≠ 1`, by strict Cauchy–Schwarz.
**Two refinements over Round 2:** `λ_min` is an **even** function of `a` (because `ρ ↦ ρ^{-1}`
swaps `u, w` and fixes `T_R`), so

```
λ_min(T_R) = −C(R+2,3) a² + O(a⁴),    a = log ρ,
```

with remainder `O(a⁴)` rather than `O(a³)`, and the constant exact at printed precision in
`a` (`1.0, 4.0, 20.0, …`) versus `0.99999999, …` in `ε = ρ−1`. The coefficient is `2 ×` the
**centred** second moment `Σ(i − R/2)² = R(R+1)(R+2)/12` (identity verified exactly for
`R = 1…199`). For `R = 1`, exactly `λ_min = −4 sinh²(a/2)`. *Lean: not attempted.*

**(b) ζ side — T1 bound.** The task's normalisation is confirmed to be the **count-preserving**
one (a *double* on-line pair at `γ` → one off-line quartet); the other two conventions leave a
`+2F²` or `+4F²` term. With that convention `Q_δ = Q_0 − 4δ²(F'² + FF'') + R_4`, verified to
ratio 0.99999999 at `δ = 10^{-4}`. Re-derived Cauchy–Schwarz bounds, valid for **all** `γ`:
`|F|² ≤ 2L`, `|F'|² ≤ 2L³/3`, `|F''|² ≤ 2L⁵/5`. Since `conj G(δ) = G(−δ)`, `Re G(δ)²` is even
in `δ`, giving the explicit remainder `|R_4| ≤ (16/3)L⁵δ⁴e^{2Lδ}`. Hence

> `λ_min(Q_δ) ≥ λ*(L) − c(L)δ² − (16/3)L⁵δ⁴e^{2Lδ}`,  `c(L) = 8L³(1/3 + 1/√5) ≈ 6.244 L³`.

**Unconditional exactly where positivity is certified:** for `L ≤ 0.8`, since Zhu certifies
`λ*(0.8) ≥ 8.9e−18` and `λ*` is non-increasing in `L`. For `L > 0.8` it holds with
`λ*(L) → 0`, which uses no positivity input at all. The bound is valid but loose by a factor
5–10 (`C/c` = 0.097 → 0.202), entirely because of the `F F''` allowance; at the minimiser
`F(γ) ≈ 0`, so the sharp behaviour is `4|F'|²`, which saturates `4L³/3` to within 5%.

**(c) Regime, in the theorem statement.** As `δ → 0`, `λ_min(Q_δ) → λ*(L) > 0`. The
`−C(L)δ²` law describes only `C(L)δ² ≫ λ*(L)`, and `C(L)` is a property of the
**re-optimised** minimiser — `4(F'²+FF'')` at the *unperturbed* ground state is 3–6 orders
too small.

---

## Non-claims

- **Nothing here proves, advances, or provides evidence for RH for ζ.** Stage 1 runs in a
  setting where the Riemann hypothesis is a theorem (Weil); Stage 2 uses *fabricated*
  spectral data chosen by hand to be off-line.
- The Stage 2 configurations are not L-functions of anything. They are
  functional-equation-closed multisets built to violate `|β| = 1`.
- F6 and F7 are quantitative statements **about the function-field window form only**.
  Their transfer to the ζ window form (different kernel, continuum of "frequencies",
  archimedean term, infinitely many zeros) is T3 conjecture until Stage 3 measures it.
- No claim of novelty is made for H1/H2: they are Carathéodory–Fejér and Pisarenko
  harmonic decomposition in Toeplitz form — and H1 is CvS's own Step 1 (Stage 4).
- **Stage 3 proves nothing about RH either.** Gate A *assumes* RH (it compares the
  geometric side to a sum over zeros written as `1/2 + i gamma` with real `gamma`); it is a
  check that my assembly is correct, not evidence for anything. Stage 3e is built from zero
  data and is a diagnostic only. Stage 3b reproduces, at smaller scale, a computation
  already published by Groskin.
- `lambda*(0.8) <= 2.2702e-17` is a **variational upper bound** from a finite basis. It is
  consistent with Zhu's certified interval but is not itself a certified enclosure: I have
  no rigorous lower bound, and the discretisation error is not bounded, only observed to be
  monotone.
- The Stage 3e sensitivity statement (`delta < 8.6e-9` at `L = 0.8`) is conditional on one
  planted height, `K = 400` zeros, and Zhu's number; it is not a statement that no off-line
  zero exists.
- Stage 3f (Davenport–Heilbronn) was not attempted.
- No Lean formalisation was attempted (Stage 1d, optional).

## Retraction log

1. **Frobenius radius check was wrong in the first run.** I computed the roots of `P(T)`
   rather than of its reverse, so I was measuring `|1/α_j|` and reported
   `max||α|/√q − 1| = 0.8` for q = 5. Fixed in `frobenius_angles`; the correct value is
   ≤ 1.34e−51. The `θ_j` were unaffected (the root set is negation-symmetric in angle),
   so no downstream result changed.
2. **The float64 column was truncated in the console display**, making `−5.148e−17` read
   as `−5.148`. Display-only; the stored values were always correct.
3. **The first sign test was tolerance-based and misclassified exact zeros.** A fixed
   `10^{−(dps−10)}` threshold called the on-line control "negative" at large `R`, where
   `λ_min` is exactly 0. Replaced by an error-based test (precision doubling, require
   `|λ| > 100·err`) and, where the data is rational, by an exact rational certificate.
   The on-line control is now correctly reported as never negative.
4. **H2 as originally stated was refuted by my own adversarial pass** and has been
   restated with a distinctness hypothesis; see Stage 1, H2.
5. **H4 as a monotone-error claim is not supported** and has been downgraded to
   "generally decreasing", with `λ_min(R)` identified as the correct monotone quantity.
6. **"(iv) is free" is RETRACTED.** I scored step (iv) by measuring convergence of the
   ground state's *zero angles* and found it free. That is the wrong object: in CvS, (iv)
   is Hurwitz transfer — convergence of the normalized minimizer **function**, with zero
   convergence as a corollary, not a substitute. Measured at the function level (M1, M1′,
   M2 in F3), (iv) **fails for every planted configuration and carries RH content**, and
   there is a short proof (H1 + Hurwitz) that it must. The verdict table now scores both
   readings separately. The earlier claim that this "refutes the brief's expectation" is
   withdrawn: the expectation was right.
7. **The ε ≈ 4e−9 sensitivity figure and the masking result are retagged T3/toy.** Both
   were stated with enough hedging to be defensible but not enough to stop them being
   quoted as ζ-side claims. They are properties of a finite-rank toy Toeplitz form; the
   ζ-side constants are not computed until Stage 3e.

8. **A first Stage 3e run was discarded, not reported as data.** At 40 digits the control
   (true zeros only, which is PSD *by construction*) came out negative for `L >= 3.5`, and
   the planted curves were non-monotone in `L`. Both were artefacts: the working precision
   had been exhausted (`lambda_min` there is ~1e−44 against entries of order 1), and the
   basis size was held fixed while `L` grew, so the basis's top frequency `N pi / 2L` was
   *falling*. Rerun at 80 digits with `N` scaling as `L`, using the control as an explicit
   noise floor. Only the corrected run is in `data/stage3e_teeth.csv`.
9. **Three bugs found and fixed during the Stage 3 refactor**, all caught by regression
   against previously verified numbers: when I renamed the even sector to add a second
   even basis, `pole_sign`, `F` and `F_at_half` were left dispatching on the old name, so
   the original basis silently picked up the odd-sector pole term (`lambda_min` went from
   `+1.73e-10` to `-6.34`). Caught because I had kept a known-good value to regress
   against; the lesson is that the regression check, not the new result, is what found it.

10. **(superseded by 11.)** A Stage 3e conclusion was corrected once in Round 1, on the
    basis of a height sweep that Round 2 then showed was itself unsound. Kept here only so
    the sequence of corrections is legible.

11. **"Low-lying zeros only" is RETRACTED as a structural limitation.** Round 1 concluded
    that an off-line zero at height `gamma_* = 100` or `200` is undetectable out to
    `L = 3`, and scoped the sensitivity statement to low-lying zeros as a property of the
    method. The conclusion was wrong, for two reasons:
      * the on-line side was a `K`-truncated zeros sum, which omits positive mass ~1e-6
        while the signal near the onset is ~1e-27 — so it was biased toward DETECTION and
        could certify nothing;
      * and the non-detections were limited by basis size and window, i.e. by the same
        **certification asymmetry** that this round is about: a finite basis finding no
        negative direction is not evidence of positivity.
    *Amended (Round 3):* Round 1 also planted by **adding** a quartet rather than moving a
    zero. That is a legitimate counterfactual in its own right — it is the one Bombieri
    (2000) §13 studies, at `rho_0 = 0.52 + 3.14i`, where `3.14` is likewise not an ordinate —
    and the added-quartet contribution `+4F(gamma)^2` simply makes detection harder away
    from a zero of `F`. The present work studies the **count-preserving** (zero-moving)
    counterfactual instead, because it is the one matched to the functional equation and
    the only one continuous as `delta -> 0`. Neither choice invalidates the other.
    Redone with the geometric on-line side and a zero-moving perturbation planted at actual
    ordinates: `gamma_30 ~ 101.3` is detected at `L = 1.18`, i.e. `0.85 x L_pred`.
    **The correct statement is (T2): detection requires `L >~ L_pred(gamma_n)` AND a basis
    with `N >~ 2n`; the Round-1 "height limitation" was the second condition misread as the
    first.** It is a limitation of that computation, not of the method.
12. **A Stage 3e run at 40 digits was discarded for a second, different reason.** Beyond
    the Round-1 precision problem, the Task-1 sweep at `dps = 30 + 40L` produced
    "detections" at `gamma_80` whose magnitude (6.1e-81) sat *below* the measured noise
    floor (the PSD-by-construction control had gone to −3.7e-76). Those rows are discarded
    and a `trustworthy` column now records the test. `gamma_80` is reported as
    **under-resolved, not undetectable** (`N >~ 160` would be needed; dim 59 was reached).
13. **Three framing corrections, applied rather than defended.** (i) The Rosati positivity
    criterion is classical linear algebra and is why Weil's proof works — only the
    identification of `T_R` with the Rosati Gram matrix is offered as a contribution.
    (ii) CvS's theorem is stated for arbitrary real even distributions, so RH-freeness of
    step (i) is implicit in their hypothesis class; the planted control makes it concrete,
    it does not discover it. (iii) The `epsilon ~ 4e-9` figure and the masking result keep
    their T3/toy tags from Round 1.
14. **The Task-4b expectation was refuted, and I report that rather than the expectation.**
    The brief anticipated that `A`'s inertia `(1, 0, N-1)` would prove basis-dependent as
    higher frequencies entered. It is not: the negative index of `A` is 1 for `L <= 0.8`
    and 2 at `L = log(13)/2`, stable across `N = 6..24` and **identical in both even
    bases**. `A` gains only positive directions.

15. **"Block F holds on the zeta side" is RETRACTED.** Rounds 1 and 2 reported `A` and `M`
    both indefinite on ζ windows with `S` PSD, noting only that the inertia *values* were
    convention-dependent. Task 2 shows the indefiniteness itself is: moving the `−log π`
    term from `A` to the prime side makes **both `A` and `M` positive definite** at
    `L ≤ 0.8`, in all three bases. Block F is retained only for the function field, where
    the pole part is an exact rank-2 hyperbolic plane; the ζ pole term `2F(i/2)²` is rank-one
    and definite, so there is nothing structural to find. (There was no Round-1/Round-2
    contradiction: the two tables had different second columns, ODD versus a second even
    basis.)
16. **Three claimed contributions are published theorems.** The Round-2 novelty check was
    abstract-level, and said so, but it was wrong in substance. Reading the full texts:
    H1 is **Hallouin–Perret Theorem 36(ii)**; H2 including the kernel generator is their
    **Theorem 6**; the window Toeplitz matrix is their **Proposition 5** entry for entry, so
    the "Toeplitz–Rosati bridge" is **not** a new identification; the Round-2 kernel-ideal
    result is **CvS Remark 2.3** in words; the Stage 2 negative-index count is **Bombieri
    Theorem 8**; and the planted-quartet windowed experiment with a critical window size is
    **Bombieri §13**, from 2000. What survives as plausibly new is the `δ²`–`L³` sensitivity
    law and the end-to-end control-experiment framing. See `PRIOR_ART.md`.

17. **The `N >~ 2n` detection rule is RETRACTED.** Round 2 stated that detecting an off-line
    zero at the `n`-th ordinate requires a basis with `N >~ 2n`, and Round 3 carried it into
    the paper, where it contradicted both the limitation ("largest subspace had dimension
    59") and the table row reporting `gamma_80` detected. Checking which basis the
    `gamma_80` detection actually used settles it: dimension 59, top mode index ~134, with
    `2n = 160`; and `gamma_30` was detected at dimension 46, top index 45, with `2n = 60`.
    **Neither reaches `2n`, so the rule is false as a necessary condition.** It was imported
    from 3b, where the `N/2` rule governs how many leading zeros a *contiguous* basis can
    annihilate — a strictly stronger demand than detection, which only needs
    `Q_geom(f) + 2F(gamma)^2 < 4 delta^2 (F'^2 + F F'')` for some `f` in the subspace. The
    sparse "low block + resonant block" subspace meets that at a small fraction of the
    dimension. The rule is replaced by: the window must reach the ordinate
    (`L >~ L_pred`), and the basis must contain modes resonant with `gamma_n`
    (`k* ~ gamma_n L/pi`). The table row and the limitation were both correct; only the rule
    was wrong.
