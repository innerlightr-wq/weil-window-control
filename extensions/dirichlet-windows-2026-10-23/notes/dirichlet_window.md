# The finite-window Weil form for L(s, χ)

Phase 1. Derivation, in the normalisation of `notes/stage3_assembly.md` of the parent
repository, with the same u-space (Parseval) treatment so that **no oscillatory quadrature
enters**. Every step carries a tier. Gate A (Phase 2) validates the assembled result
numerically; nothing below is adjusted to make Gate A pass (see the working rules in
`REPORT.md`).

Notation throughout: `f` real, even, supported in `[−L, L]`, `F(t) = ∫ f(x) e^{itx} dx`,
`C_{jk}(u) = ∫ f_j(x) f_k(x − u) dx`, `U = 2L`, `γ_E` = Euler–Mascheroni.

---

## 1.0 What is being generalised, and from where [T1, cited]

The ζ-side form used by the parent repository is Zhu's geometric side,

```
Q_ζ(f) = 2 F(i/2)² + (1/2π) ∫_R |F(t)|² Ψ_L(t) dt,
Ψ_L(t) = Re ψ(1/4 + it/2) − log π − Σ_{log n < 2L} 2 Λ(n) n^{−1/2} cos(t log n),
```

which is the finite-window specialisation of Weil's explicit formula. The finite-window Weil
quadratic form goes back to **Yoshida (1992)** and **Bombieri (2000)**
(`PRIOR_ART.md` of the parent repository, entries verified there in full text); the ζ
assembly above is reproduced, not claimed. The Dirichlet explicit formula is equally
classical — Weil's original statement is already for Hecke characters; textbook forms are
**Iwaniec–Kowalski**, *Analytic Number Theory*, Ch. 5 (Thm 5.12, general L-functions, of
which Dirichlet is the degree-1 case) and **Montgomery–Vaughan**, *Multiplicative Number
Theory I*, Ch. 12. **What is generalised here is only the parent repository's u-space
assembly and its basis-level matrix elements**, not the explicit formula itself.

The derivation below is self-contained so that every sign and every factor of ½ can be
checked against Gate A independently of how a textbook normalises its test function.

---

## 1.1 Completed L-function and functional equation [T1]

Let χ be a primitive character mod q, q > 1, with parity

```
a ∈ {0, 1},      χ(−1) = (−1)^a        (a = 0 even, a = 1 odd).
```

Put

```
Λ(s, χ) = (q/π)^{(s+a)/2} Γ((s+a)/2) L(s, χ),
Λ(s, χ) = ε(χ) Λ(1−s, χ̄),     ε(χ) = τ(χ) / (i^a √q),   τ(χ) = Σ_{r mod q} χ(r) e^{2πir/q}.
```

For q > 1 and χ primitive, `Λ(s, χ)` is **entire**: there is no pole, which is the first
structural difference from ζ. Its zeros are exactly the non-trivial zeros of `L(s, χ)`; the
trivial zeros of `L(s, χ)` at `s = −a, −a−2, −a−4, …` are the poles of `Γ((s+a)/2)` and are
absorbed by the gamma factor, exactly as for ζ.

For **real** primitive χ (the five characters of this pass) `τ(χ) = √q` if `a = 0` and
`τ(χ) = i√q` if `a = 1`, hence

```
ε(χ) = 1        for every real primitive χ.                                    [T1]
```

Consequently `Λ(1/2 + it, χ)` is **real** for real `t`: `Λ(s̄) = conj Λ(s)` because the
Dirichlet coefficients are real, and `Λ(s) = Λ(1−s)` then gives
`Λ(1/2+it) = Λ(1/2−it) = conj Λ(1/2+it)`. This is what makes a Hardy-type `Z` function
available in Phase 2, and it is **specific to real characters** — the reason complex χ are
out of scope for this pass.

## 1.2 The logarithmic derivative [T1]

From §1.1, for `Re s > 1`,

```
(Λ'/Λ)(s, χ) = ½ log(q/π) + ½ ψ((s+a)/2) + (L'/L)(s, χ),
(L'/L)(s, χ) = − Σ_{n≥2} Λ(n) χ(n) n^{−s}.
```

Compare ζ, where `Λ_ζ(s) = π^{−s/2}Γ(s/2)ζ(s)` gives
`½ log(1/π) + ½ψ(s/2) + (ζ'/ζ)(s)`. **Three changes, and only three:**

| term | ζ | L(s, χ) |
|---|---|---|
| conductor | `½ log(1/π)` | `½ log(q/π)` |
| archimedean | `½ ψ(s/2)` | `½ ψ((s+a)/2)` |
| arithmetic | `Λ(n) n^{−s}` | `Λ(n) χ(n) n^{−s}` (so `p ∣ q` drops out, since `χ(p) = 0`) |
| pole | `s = 1` present | **absent** |

## 1.3 The explicit formula in the repository's normalisation [T1]

Let `h` be even, analytic in `|Im z| ≤ ½ + ε` with enough decay, and
`g(u) = (1/2π) ∫ h(t) e^{−itu} dt`, so `h(t) = ∫ g(u) e^{itu} du`. Writing the non-trivial
zeros as `ρ = ½ + iγ_ρ` and integrating `(Λ'/Λ)(s,χ) h(−i(s−½))` over
`Re s = 1 + ε` minus `Re s = −ε`, then folding the left line onto the right with the
functional equation (which flips the sign of `Λ'/Λ` and, by evenness of `h`, leaves the test
factor unchanged), one gets

```
Σ_ρ h(γ_ρ) = (1/2π) ∫_R h(t) [ Re ψ((1/4 + a/2) + it/2) + log(q/π) ] dt
             − Σ_{n≥2} (Λ(n)/√n) [ χ(n) g(log n) + χ̄(n) g(−log n) ]        … (EF-χ)
```

with **no pole term**.

**The ζ consistency check, which fixes every constant.** Set `q = 1`, `a = 0`, `χ ≡ 1`,
and reinstate the pole of `ζ` (residue 1 at `s = 1`, contributing `h(i/2) + h(−i/2)`):

```
Σ_ρ h(γ_ρ) = h(i/2) + h(−i/2)
             + (1/2π) ∫ h(t) [ Re ψ(1/4 + it/2) − log π ] dt
             − Σ_n (Λ(n)/√n)[g(log n) + g(−log n)].
```

For real even `f` with `h(t) = F(t)F(−t) = F(t)²` and `g = f ⋆ f` even, this is **exactly**
`Q_ζ` as assembled in `notes/stage3_assembly.md`. So (EF-χ) is normalised the same way, and
the two risky constants are pinned:

* the conductor enters as **`+ log q`, with coefficient 1, not ½**. The ½ in
  `½ log(q/π)` of §1.2 is cancelled by the two contour lines, exactly as the ½ in
  `½ ψ(s/2)` is cancelled by the symmetrisation `½[ψ(z) + ψ(z̄)] = Re ψ(z)`. Its sign is
  **positive**: `−log π` becomes `log(q/π) = log q − log π`, so the conductor term pushes
  the form **up**, in the same direction as the archimedean term and opposite to the
  (subtracted) prime term.
* the archimedean digamma argument shifts `1/4 → 1/4 + a/2`, i.e. **`1/4` for even χ and
  `3/4` for odd χ**.

**Parity of the prime side for real χ.** `χ̄ = χ`, and for real even `f` the autocorrelation
`g = f ⋆ f` is even, so the two prime terms coincide and (EF-χ) becomes

```
Σ_ρ h(γ_ρ) = (1/2π) ∫ |F(t)|² [ Re ψ((1/4+a/2) + it/2) + log(q/π) ] dt
             − 2 Σ_{log n < 2L} Λ(n) χ(n) n^{−1/2} g(log n),
```

the sum being **finite** because `g = f ⋆ f` is supported in `[−2L, 2L]`.

> **(D1) [T1] The window form.** For `f` real, even, supported in `[−L, L]`, and χ real
> primitive mod q with parity `a`,
> ```
> Q_χ(f) = (1/2π) ∫_R |F(t)|² Ψ_{L,χ}(t) dt,
> Ψ_{L,χ}(t) = Re ψ((1/4 + a/2) + it/2) + log(q/π)
>              − Σ_{log n < 2L} 2 Λ(n) χ(n) n^{−1/2} cos(t log n),
> ```
> and, if every non-trivial zero of `L(s, χ)` lies on the critical line,
> `Q_χ(f) = Σ_ρ |F(γ_ρ)|²`. The repository's normalisation is `Q_χ(f)/‖f‖₂²`.

**Note on the missing pole, and why it is not a free gain.** On the ζ side the pole term
`2F(i/2)²` is a positive rank-one contribution — it is the term that makes ζ-window
positivity plausible at small `L` at all. It is **absent** for χ, and in its place there is
the positive diagonal `log q · ‖f‖₂²`. The two are structurally different: one is rank one
and `L`-dependent through `F(i/2)`, the other is a multiple of the identity. This is the
single most important difference for Phase 3 and is the reason the conductor cannot simply
be read as "more height".

---

## 2. u-space form of every term [T1]

Exactly as in `notes/stage3_assembly.md`: with `C_{jk}(u) = ∫ f_j(x) f_k(x−u) dx`,

```
(1/2π) ∫ F_j conj(F_k) dt          = C_{jk}(0) = L δ_{jk}
(1/2π) ∫ F_j conj(F_k) cos(tu) dt  = ½ (C_{jk}(u) + C_{jk}(−u))
```

so the whole matrix is a **finite** combination of closed-form evaluations of `C`. The
`C`-machinery is **unchanged** from the ζ case — it depends only on the basis, not on the
L-function — so `src/zeta_window.py`'s `Csym_terms` / `eval_terms` are reused verbatim.
The three terms:

### 2.1 Conductor term (diagonal) [T1]

```
Cond_{jk} = log(q/π) · C_{jk}(0) / (‖f_j‖ ‖f_k‖) = log(q/π) · δ_{jk}.
```

In code this is the ζ `LogPi` matrix with `log π` replaced by `log(π/q)`, equivalently
`M_χ = Arch_a − LogPi + LogQ − Prime_χ` with `LogQ = log q · C(0)/(‖·‖‖·‖)`. Keeping the two
pieces separate is deliberate: the sign-flip control of Gate A flips `LogQ` alone.

### 2.2 Prime term [T1]

```
Prime^χ_{jk} = Σ_{n: log n < 2L} 2 Λ(n) χ(n) n^{−1/2} · C^sym_{jk}(log n) / (‖f_j‖ ‖f_k‖).
```

Identical to ζ with the **real** weight `Λ(n) χ(n)`, which is `±log p` on `p ∤ q` and `0` on
`p ∣ q`. Note this term can now have **either sign**, where for ζ it is a sum of positive
weights — a second structural difference, and the reason `Q_χ` is not a sign-perturbation of
`Q_ζ`.

### 2.3 Archimedean term [T1]

From `ψ(z) = −γ_E + Σ_{m≥0} [1/(m+1) − 1/(m+z)]` at `z = (1/4 + a/2) + it/2`: writing
`m + z = ½(c_m + it)` with

```
c_m = 2m + ½ + a            (ζ / even χ: 2m + ½;   odd χ: 2m + 3/2),
```

one has `Re 1/(m+z) = 1/(c_m + it) + 1/(c_m − it)`, hence

```
Re ψ((1/4+a/2) + it/2) = −γ_E + Σ_{m≥0} [ 1/(m+1) − 1/(c_m + it) − 1/(c_m − it) ].
```

Since `1/(c+it) + 1/(c−it)` pairs with `C` as `2∫_0^{2L} C(u) e^{−cu} du` (the parent note's
step, unchanged), the archimedean matrix element is the convergent series

```
Arch_{jk} = −γ_E C(0) + Σ_{m≥0} [ C(0)/(m+1) − 2 ∫_0^{2L} C(u) e^{−c_m u} du ],
```

**summed in closed form** term by term over the `(c, α, β)` representation of `C`. For a
term `c · sin(αu + β)`, with `z = (1/4 + a/2) − iα/2` and `w = e^{−2U}`:

```
Σ_m [ sin β/(m+1) − 2∫_0^U sin(αu+β) e^{−c_m u} du ]
   = Im{ e^{iβ} (ψ(z) + γ_E) } + e^{−(1/2+a)U} Im{ e^{i(β+αU)} Φ(w, 1, z) },
```

and for a term `c · (U − u) cos(βu)`, with `z₀ = (1/4 + a/2) − iβ/2`:

```
Σ_m [ U/(m+1) − 2∫_0^U (U−u)cos(βu) e^{−c_m u} du ]
   = U · Re(ψ(z₀) + γ_E) + ½ Re{ ζ(2, z₀) − e^{iβU} e^{−(1/2+a)U} Φ(w, 2, z₀) },
```

`Φ` the Lerch transcendent, `ζ(2, ·)` Hurwitz. **Derivation of the first line** (the second
is the same computation with one integration by parts): `2∫_0^U e^{i(αu+β)} e^{−c_m u} du
= e^{iβ}·2[1 − e^{(iα−c_m)U}]/(c_m − iα)` and `c_m − iα = 2(m+z)`, so
`2/(c_m − iα) = 1/(m+z)` and `e^{−c_m U} = e^{−(1/2+a)U} w^m`; take imaginary parts and use
`Σ_m[1/(m+1) − 1/(m+z)] = ψ(z) + γ_E`, `Σ_m w^m/(m+z)^s = Φ(w,s,z)`.

> **(D2) [T1] The archimedean generalisation is exactly two substitutions** in the parent
> repository's `arch_closed`:
> ```
> z  : 1/4 − iα/2        →   (1/4 + a/2) − iα/2
> prefactor : e^{−U/2}   →   e^{−(1/2+a)U}
> ```
> and nothing else. At `a = 0` it reduces to the parent code character for character; this
> reduction is the regression test (`src/dirichlet_window.py::arch_closed_a`, asserted
> against `zeta_window.arch_closed` at `a = 0`).

## 3. Both parity sectors of the test function [T1]

The **test-function** parity sectors are those of the parent note and are *independent of*
the **character** parity `a`; the two must not be conflated, and the code names them
`sector` and `a_par` respectively.

```
EVEN_D : φ_k(x) = cos((2k+1)πx/2L),   k = 0,1,2,…      (vanishes at ±L)
EVEN   : φ_k(x) = cos(kπx/L),         k = 0,1,2,…      (contains the constant)
ODD    : ψ_k(x) = sin(kπx/L),         k = 1,2,…
```

All three give `∫ f_j f_k = L δ_{jk}` (except `EVEN, k=0`, where it is `2L`), and the
cross-correlation machinery is shared. **The pole-sign bookkeeping of the ζ case disappears**:
`pole_sign(sector)` multiplied the pole term, and for χ there is no pole term, so the odd
sector differs from the even one only through `C` — i.e. through the basis. In particular the
odd-sector χ matrix is

```
M^odd_χ = Arch_a − LogPi + LogQ − Prime_χ
```

with the odd-basis `C`, and there is no sector-dependent sign anywhere in the assembly. This
is a genuine simplification relative to ζ and it removes one of the two places where the ζ
code could carry a sector/sign bug — which is also why the ζ control must be kept in the
suite: a sector bug on the χ side would no longer be caught by an internal inconsistency.

## 4. Assembled matrix [T1]

In the orthonormal basis `f_k/√(norm2_k)`:

```
M^χ_{jk} = [ Arch^{(a)}_{jk} + log(q/π) C_{jk}(0)
             − Σ_{log n < 2L} 2 Λ(n) χ(n) n^{−1/2} C^sym_{jk}(log n) ] / (‖f_j‖‖f_k‖)
```

and `λ_min(M^χ)` over an `N`-dimensional nested basis is a **variational upper bound** for
`λ*_χ(L) = inf { Q_χ(f)/‖f‖₂² }`, monotone non-increasing in `N`. It is **never** used as a
lower bound (`REPORT.md` non-claim 8).

## 5. Planting: the three counterfactuals [T1]

Let `ρ = ½ + δ ± iγ` with `δ ≥ 0`. For real even `f`, `F(γ ± iδ) = A_δ ∓ i B_δ` with

```
A_δ(f) = ∫ f(u) cos(γu) cosh(δu) du,      B_δ(f) = ∫ f(u) sin(γu) sinh(δu) du,
```

both real, and the four z-ordinates of the quartet `{½ ± δ ± iγ}` are `±γ ± iδ`, giving the
quartet contribution

```
Σ_{quartet} F(z)² = 4 Re F(γ + iδ)² = 4 (A_δ² − B_δ²).
```

So the three conventions of `notes/proofs.md` §B.1 carry over verbatim, with `Q_ζ → Q_χ`:

| convention | perturbation | used in |
|---|---|---|
| **add** a quartet, true zeros left in place | `+4(A_δ² − B_δ²)` | Phase 4 (added-quartet), Phase 5 |
| **move** a double on-line zero at `γ` to the quartet | `−4F(γ)² + 4(A_δ² − B_δ²)` | Phase 4 (moved-zero) |
| move one on-line **pair** at `γ` to the quartet | `−2F(γ)² + 4(A_δ² − B_δ²)` | Phase 4 (moved-zero, count-preserving) |

The two Phase-4 counterfactuals are kept distinct, as in §10 of the paper.

### 5.1 The γ = 0 degeneration — the Siegel configuration [T1]

A **real** zero `β = ½ + δ` of `L(s, χ)` forces, by the functional equation, the partner
`1 − β = ½ − δ`. In the z-variable `ρ = ½ + iz` these are `z = ∓iδ`: a pair of **purely
imaginary** ordinates, which is precisely `γ = 0` in the parametrisation above. Its
contribution is

```
Σ_{pair} F(z)² = 2 F(iδ)² = 2 (∫ f(u) cosh(δu) du)² = 2 A_δ(f)²   ≥ 0,
```

since at `γ = 0` one has `B_δ ≡ ∫ f(u) sin(0) sinh(δu) du = 0` identically. Equivalently: the
`δ²` coefficient of Theorem 8, `4K(L, γ)` with `K(L,0) = ∫ u² sin²(0·u) du = 0`, vanishes —
but the exact statement is stronger than a vanishing leading coefficient, and is recorded as
Phase 5's main result rather than here.

---

## 6. What Gate A must discriminate [T1, design]

Gate A compares the assembled geometric side of (D1) against `Σ_{ρ} |F(γ_ρ)|²` over the first
`K` computed zeros. The tail is `Σ_{|γ|>γ_K} |F(γ)|²`, and the zero-counting density is
`(1/π) log(q|t|/2π) + O(1)`, so the decay rate is set by how fast `|F_k(t)|²` falls — which is
a property of the **basis**, not of the L-function:

| basis | value at `±L` | `F_k(t)` as `t → ∞` | `\|F_k\|²` | **predicted exponent** |
|---|---|---|---|---|
| `EVEN_D` `cos((2k+1)πx/2L)` | `0` (continuous) | `(−1)^k 2a_k cos(tL)/(a_k²−t²) = O(t^{−2})` | `O(t^{−4})` | **−3** (times `log qγ_K`) |
| `ODD` `sin(kπx/L)` | `0` (continuous) | `−(−1)^k 2b_k sin(tL)/(b_k²−t²) = O(t^{−2})` | `O(t^{−4})` | **−3** (times `log qγ_K`) |
| `EVEN` `cos(kπx/L)` | `(−1)^k ≠ 0` (**jump**) | `L[sinc((w−t)L)+sinc((w+t)L)] = O(t^{−1})` | `O(t^{−2})` | **−1** (times `log qγ_K`) |

```
discrepancy(K) ≍ γ_K^{−3} log(q γ_K)      for EVEN_D and ODD
discrepancy(K) ≍ γ_K^{−1} log(q γ_K)      for EVEN   (the basis jumps at ±L)
```

The brief's `γ_K^{−3}` is stated "for modes vanishing at `±L`", which is exactly the first two
rows. A fit exponent near the row's prediction is the pass condition; an exponent near `0`, or
a plateau, is a term error. **A shallow exponent in the `EVEN` basis alone is not a failure**,
and the ζ control is run in all three sectors so that the three predicted rates can be read off
a known-good case in the same run.

**Sign-flip controls.** Each run also assembles, and reports, the matrices obtained by
flipping exactly one term: `LogQ → −LogQ` and `Arch → −Arch`. If Gate A's discrepancy does
not blow up under each flip separately, the test does not discriminate that term and Gate A
is reported as **uninformative for it**, not as passed.
