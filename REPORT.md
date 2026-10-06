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
`results/stage1_spectra.csv` (40 rows, 4 curves, R = 0..4g).

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
(`src/gram_trace.py`, `results/stage1e_gram_trace.csv`, 40 dps):

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
**every** row of `results/stage2_spectra.csv` once `R+1 ≥ rank`.

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

**F1. The witness is free (T1 + T2).** H1 holds at **every** applicable row of **every**
planted configuration, including rows where `λ_min = −642.8`. Measured
`max||z|−1| ≤ 1.6e−48`. All zeros of the approximant lie exactly on the critical circle
whether or not the underlying spectrum does. *"The approximant's zeros are real" carries
no RH content whatsoever.*

**F2a. The kernel is an ideal; the CvS hypothesis fails for every `R > d` (T1 + T2).**
Let `P*(z) = Π (z − β_j)` over the **distinct** `β_j`, `d = deg P*`. Then
`Q(c) = 0 ⟺ ĉ(β_j) = 0 ∀j ⟺ P* | ĉ`, so

```
ker T_R = { P*(z) q(z) : deg q ≤ R − d },     dim ker T_R = max(0, R + 1 − d).
```

Predicted nullity matched the computed nullity in **every** row of
`results/stage2_function_level.csv` (all four curves, the repeated-angle curve, and the
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
| (i) real zeros of the approximant | — | **free** — Carathéodory–Fejér gives it from Toeplitz structure plus simplicity. True verbatim for planted off-line spectra, even at `λ_min = −642` (F1). |
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
  harmonic decomposition in Toeplitz form. Stage 4 (literature check) has not been run.
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
