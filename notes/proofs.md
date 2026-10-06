# Proofs for the δ² law (Round 3, Task 3)

Every step below is checked numerically by `src/proofs_check.py`; the residuals quoted are
from that script at 60–80 digits.

---

## A. Function field (T1)

**Setting.** A functional-equation-closed off-circle *real* pair of normalised Frobenius
eigenvalues `{ρ, ρ^{-1}}`, `ρ > 1`. The window data is `t(n) = ρ^n + ρ^{-n}` and
`T_R = (t(|i−j|))_{0≤i,j≤R}`. Put `a = log ρ > 0`.

### A.1 Rank-2 factorisation

**Lemma A1.** With `u = (ρ^i)_{i=0}^R` and `w = (ρ^{-i})_{i=0}^R`,
`T_R = u w^T + w u^T`.

*Proof.* `(u w^T + w u^T)_{ij} = ρ^i ρ^{-j} + ρ^{-i} ρ^j = ρ^{i-j} + ρ^{-(i-j)}`, and
`m ↦ ρ^m + ρ^{-m}` is even, so this equals `ρ^{|i-j|} + ρ^{-|i-j|} = t(|i−j|)`. ∎
*Checked: residual ≤ 6.2e−61.*

### A.2 The spectrum

**Proposition A2.** `T_R` has rank ≤ 2. Its nonzero eigenvalues are
`(R+1) ± |u||w|`, and all other eigenvalues are 0. In particular

> **λ_min(T_R) = (R+1) − |u||w|,**
> `|u|² = (ρ^{2R+2}−1)/(ρ²−1)`, `|w|² = (ρ^{−2R−2}−1)/(ρ^{−2}−1)`.

*Proof.* If `x ⟂ u, w` then `T_R x = u(w·x) + w(u·x) = 0`, so `ran T_R ⊆ span{u,w}` and the
rank is ≤ 2. On `span{u,w}`, for `x = αu + γw`,
`T_R x = (αp + γ|w|²) u + (α|u|² + γp) w` where `p = u·w`. In the (generally non-orthogonal)
basis `{u,w}` the operator is therefore represented by `[[p, |w|²],[|u|², p]]`, whose
eigenvalues — being those of the operator, hence basis-independent — are `p ± |u||w|`.
Finally `p = u·w = Σ_{i=0}^R ρ^i ρ^{-i} = R+1`. Since `T_R` is symmetric these are real,
and `|u||w| ≥ |u·w| = R+1` by Cauchy–Schwarz, so `(R+1) − |u||w| ≤ 0` is the least. ∎
*Checked: eigenvalue residual ≤ 1.0e−59; middle eigenvalues ≤ 8.1e−60.*

**Corollary A3 (strictness).** `λ_min(T_R) < 0` for every `R ≥ 1` as soon as `ρ ≠ 1`, since
Cauchy–Schwarz is an equality iff `u ∥ w` iff `ρ = 1`. *A single off-circle pair is seen by
the very first nontrivial window.*

### A.3 Parity in `a`, and the expansion

**Lemma A4.** `λ_min` is an **even** function of `a = log ρ`.

*Proof.* `ρ ↦ ρ^{-1}` interchanges `u` and `w`, and `T_R = uw^T + wu^T` is symmetric in the
pair, hence unchanged; but `a ↦ −a`. ∎
*Checked: `|λ(ρ) − λ(1/ρ)| ≤ 3.4e−80`.*

This forces the remainder below to be `O(a⁴)`, not `O(a³)`.

**Theorem A5.** With `S_1 = Σ_{i=0}^R i`, `S_2 = Σ_{i=0}^R i²`,

```
λ_min(T_R) = −2a² [ S_2 − S_1²/(R+1) ] + O(a⁴) = −C(R+2,3) a² + O(a⁴),
```

where `C(R+2,3) = R(R+1)(R+2)/6` is the binomial coefficient.

*Proof.* `|u|² = Σ e^{2ia} = (R+1) + 2aS_1 + 2a²S_2 + O(a³)` and `|w|²` is the same with
`a → −a`. Hence
`|u|²|w|² = [(R+1) + 2a²S_2]² − (2aS_1)² + O(a³) = (R+1)² + 4a²[S_2(R+1) − S_1²] + O(a³)`,
so
`|u||w| = (R+1)·(1 + 4a²[S_2(R+1)−S_1²]/(R+1)²)^{1/2} + O(a³)
       = (R+1) + 2a²[S_2 − S_1²/(R+1)] + O(a³)`,
and `λ_min = (R+1) − |u||w| = −2a²[S_2 − S_1²/(R+1)] + O(a³)`. By Lemma A4 the `O(a³)` term
vanishes and the remainder is `O(a⁴)`.

For the constant, `S_2 − S_1²/(R+1) = Σ_{i=0}^R (i − R/2)²` is the **centred** second moment
of the window, and with `S_1 = R(R+1)/2`, `S_2 = R(R+1)(2R+1)/6`,

```
S_2 − S_1²/(R+1) = R(R+1)[(2R+1)/6 − R/4] = R(R+1)(R+2)/12,
```

whence `λ_min = −2a²·R(R+1)(R+2)/12 = −C(R+2,3) a²`. ∎
*Checked: the identity holds exactly for R = 1…199; the ratio `−λ/a²` equals
`C(R+2,3)` to 12 printed digits for R = 1, 2, 4, 8, 16, 32; `(λ + C a²)/a⁴` converges
(−1/12, −4/3, −68/3, −472 for R = 1, 2, 4, 8), confirming the `O(a⁴)` remainder.*

**Remark A6 (`a` versus `ε`).** In `ε = ρ − 1` the same statement reads
`λ_min = −C(R+2,3)ε² + O(ε³)` — this is the Round-2 form — but `a = log ρ` is the natural
variable: it makes the coefficient exact at the printed precision (`1.0, 4.0, 20.0, …`
versus `0.99999999, 3.99999996, …` in `ε`) and the remainder one order better.

**Remark A7 (the case R = 1, exactly).**
`λ_min(T_1) = 2 − (ρ + ρ^{-1}) = −(e^{a/2} − e^{−a/2})² = −4 sinh²(a/2) = −a² − a⁴/12 − …`
— manifestly even, with `C(3,3) = 1`. *Checked to 7.7e−80.*
**Lean formalisation was not attempted**; R = 1 is recorded here as the hand-checkable case.

---

## B. ζ side (T1 for the bound; T2 for saturation)

**Setting.** `f` real, even, supported in `[−L, L]`, `‖f‖₂ = 1`; `F = f̂`. For real even `f`,
`F` is real on `R` with real Taylor coefficients there, so `F(z̄) = conj F(z)`.

### B.1 Normalisation — which perturbation gives the stated formula

The on-line pair `{±γ}` contributes `2F(γ)²`; an off-line quartet `{±γ ± iδ}` contributes
`4 Re F(γ+iδ)²`. Expanding `F(γ+iδ) = F + iδF' − (δ²/2)F''+ O(δ³)` (all at `γ`):

```
4 Re F(γ+iδ)² = 4F² − 4δ²(F'² + F F'') + O(δ⁴),
```

the odd powers being purely imaginary. Hence **three different conventions**:

| convention | perturbation `Q_δ − Q_0` | `δ → 0` limit |
|---|---|---|
| add a quartet to all true zeros | `+4F² − 4δ²(F'²+FF'')` | `+4F²` (zeros added) |
| one pair → one quartet | `+2F² − 4δ²(F'²+FF'')` | `+2F²` (count doubles) |
| **double pair → one quartet** | **`−4δ²(F'² + F F'')`** | **0 (count preserved)** |

**The task's formula is the third**, i.e. a *double* on-line zero at `γ` splitting into an
off-line quartet. It is the only one continuous as `δ → 0` and the only one preserving the
zero count, so it is the convention adopted. *Checked: the ratio of the exact
`−4F² + 4ReF(γ+iδ)²` to `−4δ²(F'²+FF'')` is 0.999964, 0.9999996, 0.99999999 at
`δ = 10^{-2}, 10^{-3}, 10^{-4}`.*

### B.2 Cauchy–Schwarz bounds

For real even `f` with `‖f‖₂ = 1`, `F(t) = ∫ f(u)cos(tu)du`, so

```
F'(t)  = −∫ u f(u) sin(tu) du,     F''(t) = −∫ u² f(u) cos(tu) du.
```

Cauchy–Schwarz, then dropping the oscillatory weight (valid for **every** `γ`):

```
|F(γ)|²   ≤ ∫_{−L}^{L} cos²(γu) du     ≤ 2L,
|F'(γ)|²  ≤ ∫_{−L}^{L} u² sin²(γu) du  ≤ 2L³/3,
|F''(γ)|² ≤ ∫_{−L}^{L} u⁴ cos²(γu) du  ≤ 2L⁵/5.
```

*Checked at L = 0.5, 0.8, 1.0, 1.3: all three hold with room.* Keeping the weight gives
`∫u²sin² → L³/3` and `∫u⁴cos² → L⁵/5` as `γ → ∞`, i.e. a factor 2 better asymptotically.

Hence `|F'² + F F''| ≤ 2L³/3 + √(2L)·√(2L⁵/5) = 2L³(1/3 + 1/√5)`.

### B.3 The remainder

With `G(δ) = F(γ+iδ) = ∫ f(u)e^{iγu}e^{−δu}du = Σ_k (−δ)^k m_k/k!`, `m_k = ∫ f u^k e^{iγu}`,
Cauchy–Schwarz gives `|m_k| ≤ (∫_{−L}^{L}u^{2k})^{1/2} = (2L^{2k+1}/(2k+1))^{1/2} ≤ √(2L) L^k`.
So `G(δ)² = Σ_n c_n δ^n` with `|c_n| ≤ 2L (2L)^n/n!`. Since `conj G(δ) = G(−δ)`,
`Re G(δ)²` is **even** in `δ`, so all odd terms vanish and

```
|R_4| := |4 Re G(δ)² − 4F² + 4δ²(F'²+FF'')| ≤ 8L Σ_{n≥4} (2Lδ)^n/n! ≤ (16/3) L⁵ δ⁴ e^{2Lδ}.
```

*Checked at L = 0.8: actual remainder 4.4e−9, 4.4e−13, 4.4e−17 against the bound
2.1e−4, 1.8e−8, 1.8e−12 at `δ = 10^{-1}, 10^{-2}, 10^{-3}`.*

### B.4 The theorem

> **Theorem B.** Let `f` be real, even, supported in `[−L,L]`, with `‖f‖₂ = 1`. Let `Q_0` be
> the Weil window form (geometric side) and `Q_δ` the form obtained by moving a double
> on-line zero at ordinate `γ` off the critical line to `½ + δ ± iγ` (§B.1, third
> convention). Then
>
> ```
> Q_δ(f) ≥ λ*(L) − c(L) δ² − (16/3) L⁵ δ⁴ e^{2Lδ},     c(L) = 8L³(1/3 + 1/√5) ≈ 6.244 L³,
> ```
>
> and hence `λ_min(Q_δ) ≥ λ*(L) − c(L)δ² − (16/3)L⁵δ⁴e^{2Lδ} ≥ −c(L)δ² − (16/3)L⁵δ⁴e^{2Lδ}`.
> In particular `c(L) = O(L³)`.

*Proof.* `Q_δ(f) = Q_0(f) − 4δ²(F'(γ)² + F(γ)F''(γ)) + R_4` by §B.1, with `|R_4|` bounded in
§B.3. `Q_0(f) ≥ λ*(L)‖f‖² = λ*(L)`, and `4|F'²+FF''| ≤ 8L³(1/3+1/√5)` by §B.2. ∎

**Which windows is this unconditional on?** The inequality `Q_0 ≥ λ*(L) > 0` is *not* free —
it is Weil positivity on the window. It is **unconditional exactly where positivity has been
certified**: Zhu certifies `λ*(0.8) ∈ [8.9×10^{-18}, 2.27×10^{-17}]`, and `λ*` is
non-increasing in `L` (larger window = larger test-function space), so

> **for `L ≤ 0.8` the bound holds unconditionally with `λ*(L) ≥ 8.9×10^{-18} > 0`.**

For `L > 0.8` the bound still holds with `λ*(L)` replaced by `0`, but then it says only
`λ_min(Q_δ) ≥ −c(L)δ² − …`, which is unconditional (it uses no positivity at all).

*Checked: `c(L)` exceeds the measured `C(L)` at every `L` tested
(ratios `C/c` = 0.097, 0.137, 0.182, 0.192, 0.202 at `L` = 0.8, 1.0, 1.3, 1.6, 2.0), so the
bound is valid but loose by a factor 5–10. The looseness is entirely the `F F''` allowance:
at the minimiser `F(γ) ≈ 0` (3.6e−16 down to 2.2e−53), so the sharp behaviour is governed by
`4|F'|²` alone, for which the refined bound `4L³/3` is approached to within 5% — see §B.5.*

### B.5 Regime of validity — part of the statement, not a footnote

> **The `−C(L)δ²` law is not a statement about `λ_min(Q_δ)` for all `δ`.** As `δ → 0`,
> `λ_min(Q_δ) → λ*(L) > 0`: the form returns to the honest positive floor. The law
> `λ_min ≈ −C(L)δ²` describes only the regime
>
> ```
> C(L) δ² ≫ λ*(L).
> ```
>
> Moreover `C(L)` is **not** `4(F'²+FF'')` evaluated at the unperturbed ground state — that
> value is 3 to 6 orders of magnitude too small (Round 2, Task 2b). `C(L)` is a property of
> the **re-optimised** minimiser of `Q_δ`, which moves far from the ground state because
> `Q_0` has a near-null space with `λ₂/λ₁ ∼ 10⁶`. The measured `C(L)` agrees with
> `4(F'²+FF'')` at the *perturbed* minimiser to 2–6%, and `C(L)/(4L³/3) → 0.95`, i.e. the
> re-optimised minimiser asymptotically saturates the Cauchy–Schwarz bound.

---

## B′. A sharp *conditional* bound: completing the square

The unconditional Theorem B is loose by a factor 5–10, and the looseness is entirely the
`F F''` allowance. In the **zero-moving** setting that term can be absorbed, at the cost of
one hypothesis.

**The hypothesis.** Write `Q̃ := Q_0 − 4F(γ)²`, the contribution of every zero *other* than
the one being moved. If all of those lie on the critical line then `Q̃ ≥ 0`. This is RH for
every zero but one — strictly weaker than RH, but not free; hence "conditional".

**The identity.** For real `F, F''` and any `δ`,

```
4F² − 4δ² F F'' = (2F − δ² F'')² − δ⁴ F''² .
```

*Checked: residual ≤ 1.2e−60 over twelve sign/magnitude combinations.*

**Theorem B′ (conditional).** Let `f` be real, even, supported in `[−L,L]`. Suppose every
nontrivial zero of `ζ` other than the one at ordinate `γ` lies on the critical line, so that
`Q̃ = Q_0 − 4F(γ)² ≥ 0`. Move the double on-line zero at `γ` to `½ + δ ± iγ`. Then

```
λ_min(Q_δ) ≥ − 4δ² · sup_γ |F'(γ)|²/‖f‖²  −  (2/5 + (16/3) e^{2Lδ}) L⁵ δ⁴
           ≥ − (8L³/3) δ²                  −  (2/5 + (16/3) e^{2Lδ}) L⁵ δ⁴ .
```

*Proof.* `Q_δ = Q̃ + 4 Re F(γ+iδ)²` by construction. Expanding (§B.1) and completing the
square,

```
4 Re F(γ+iδ)² = 4F² − 4δ²(F'² + F F'') + R_4
              = (2F − δ² F'')² − δ⁴ F''² − 4δ² F'² + R_4 ,
```

so `Q_δ ≥ 0 + 0 − 4δ²F'² − δ⁴F''² − |R_4|`. Now `|F''|² ≤ 2L⁵/5` and
`|R_4| ≤ (16/3)L⁵δ⁴e^{2Lδ}` by §B.2–B.3, and `|F'|² ≤ 2L³/3` for every `γ`. ∎

**The constant is essentially sharp, and the minimiser nearly attains it.** Keeping the
`sin²` weight gives `|F'(γ)|² ≤ ∫u² sin²(γu)du → L³/3` as `γ → ∞`, i.e. a sharp constant
`4L³/3` in the large-`γ` regime. Measured against both:

| L | measured `C` | `8L³/3` (all γ) | ratio | `4L³/3` (γ→∞) | ratio |
|---|---|---|---|---|---|
| 0.8 | 0.3108 | 1.3653 | 0.228 | 0.6827 | 0.455 |
| 1.0 | 0.8530 | 2.6667 | 0.320 | 1.3333 | 0.640 |
| 1.3 | 2.4987 | 5.8587 | 0.427 | 2.9293 | 0.853 |
| 1.6 | 4.9218 | 10.923 | 0.451 | 5.4613 | 0.901 |
| 2.0 | 10.083 | 21.333 | 0.473 | 10.667 | **0.945** |

So `C(L)` reaches **94.5%** of the sharp large-`γ` constant `4L³/3` at `L = 2` and is still
rising: the re-optimised minimiser asymptotically saturates Cauchy–Schwarz.

*Checked: the bound holds against the directly computed `λ_min(Q_δ)` at
(L, δ) = (0.8, 0.1), (0.8, 0.02), (1.0, 0.1), (1.0, 0.02), (1.3, 0.1), (1.3, 0.02), with
2.3×–4.4× of slack.*

**Status.** Theorem B is unconditional for `L ≤ 0.8` (via Zhu's certificate) and carries a
constant `≈ 6.244 L³`. Theorem B′ is conditional on all other zeros being on the line and
carries `8L³/3 ≈ 2.667 L³`, sharp to `4L³/3` asymptotically. Both give `c(L) = O(L³)`.
