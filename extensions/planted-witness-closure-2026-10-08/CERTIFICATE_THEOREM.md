# The added-quartet witness certificate

Exploration, 2026-10-08. Read-only outside this directory. No installs, commits, pushes,
uploads, branch changes, resets, background jobs or manuscript edits.

## 1. The exact function

`L = 4/5`, `γ = 14` — **both exactly rational**. Basis: the `EVEN_D` modes of
`src/zeta_window.py`,

```
phi_k(x) = cos(w_k x)/sqrt(L),   w_k = (2k+1)pi/(2L) = (2k+1)(5/8)pi,   k = 0,1,2,3.
```

`w_k L = (2k+1)π/2` exactly, so `phi_k(±L) = 0`: extended by zero each `phi_k` is
continuous on `ℝ`, real and even, supported in `[−L,L]`. Orthonormality
`∫ phi_j phi_k = δ_jk` is the `norm2 = L` convention of the source.

**The frozen coefficient vector**, the previous run's 12-decimal printout **read as exact
rationals** (documented choice; the alternative, re-solving the node condition, cannot give
a rational solution because `w_k` is a rational multiple of `π`):

```
x = ( 99815221908, -354082045264, 928280739889, 54384691304 ) / 10^12 ,
f  = sum_{k=0}^{3} x_k phi_k .
```

`f` is **not** normalised:

```
||f||_2^2 = sum_k x_k^2 = 999999999999981276942897 / 10^24  (exact)
          = 0.999999999999981...  > 0 .
```

Positive scaling cannot change a sign, so the certificate below is proved for the
unnormalised `f`; dividing by the enclosed positive `||f||²` gives the same sign for
`f/||f||`.

**Node residual, recomputed for THIS vector and not inherited:**

```
F(14) = int f(u) cos(14u) du = 4.1293580707791...e-13      (enclosed, width < 1e-70)
```

The previous run's `−4.12e-17` belonged to the un-rounded vector and is **not** used.

## 2. The object certified

The zero configuration is **the nontrivial zeros of ζ together with the four extra zeros**

```
1/2 + delta + 14i ,  1/2 + delta - 14i ,  1/2 - delta + 14i ,  1/2 - delta - 14i ,
```

each counted once (`δ > 0`, so the four are distinct). This is **not** a moved actual zeta
zero. The corresponding Weil functional, derived from the zero side under the manuscript's
real-even convention `F(z) = ∫_{−L}^{L} f(u)e^{izu}du`, is

```
Q_add,delta(f) = Q_zeta(f) + 4 [ A_delta(f)^2 - B_delta(f)^2 ] ,
A_delta(f) = int f(u) cos(14u) cosh(delta u) du ,
B_delta(f) = int f(u) sin(14u) sinh(delta u) du .
```

**Derivation, independent of the proposed formula.** A zero at `½+s` contributes `h` at the
corresponding ordinate; for the quartet the four contributions sum to `4 Re F(14+iδ)²`
(the four members pair into two conjugate pairs, each contributing `2 Re F(14+iδ)²`).
Splitting `e^{i(14+iδ)u}` into even and odd parts and using that `f` is even gives
`F(14+iδ) = A_δ(f) − i B_δ(f)`, whence `4 Re F(14+iδ)² = 4[A_δ² − B_δ²]`. At `δ = 0`:
`A_0 = F(14)`, `B_0 = 0`, so

```
Q_add,0(f) = Q_zeta(f) + 4 F(14)^2 ,
```

which is the correct baseline for *adding* a quartet — the configuration then has a double
on-line zero at 14 superimposed on ζ's zeros. **Multiplicity check:** `4F(14)²` is twice the
`2F(14)²` of a single pair, as it must be.

**Comparison with the count-preserving replacement convention** (`Q_0 = Q_ζ + 2F(14)²`,
quartet replacing a double pair), which is what `paper/weil_window_control.tex` eq. (3)–(4)
and `notes/proofs.md` §B.1 use:

```
Q_rep,delta(f) = Q_zeta(f) - 2F(14)^2 + 4[A_delta^2 - B_delta^2]
Q_add,delta(f) - Q_rep,delta(f) = 2 F(14)^2 = 3.4103196...e-25      (enclosed)
```

**The previous run's script evaluated `Q_rep`, not `Q_add`.** For this vector the two differ
by `3.41e-25`, which is `1.4e-23` of the certified margin — immaterial here, but they are
different functionals and the distinction is recorded, not elided.

## 3. The theorem

> **Theorem (added-quartet witness).** Let `f`, `L = 4/5`, `γ = 14` be exactly as in §1 and
> `Q_add,δ` exactly as in §2. Then
>
> **(a) point.** `Q_add,2/5(f) ∈ [ −0.0146094664641386, −0.0145107657498228 ]`, so
> ```
> Q_add,2/5(f)  <=  -eta_point ,     eta_point = 0.0145107657498227 .
> ```
>
> **(b) interval.** For **every** `δ ∈ [1/4, 49/100]`,
> ```
> Q_add,delta(f)  <=  -eta_interval ,   eta_interval = 0.00245490957806631 .
> ```
>
> Both with all analytic and arithmetic errors accounted for, as itemised in §4–§5.

`[1/4, 49/100] ⊂ (0, 1/2)`, so every member of the quartet lies strictly inside the critical
strip. The interval is the one proposed in the brief; it was **not** narrowed.

Context, enclosed with the same machinery and *not* part of the certificate:
`F'(14) = −0.17626230294801`, `F''(14) = +0.11323208261102521`,
`a_f = 4[F'² + F F''] = 0.12427359776233127`.

## 4. The analytic inputs, named

Four, and only four:

1. **Weil's explicit formula** for `ζ` in the `ψ`-form, i.e. the identity
   `Σ_ρ h(γ_ρ) = h(i/2)+h(−i/2) − 2Σ_n Λ(n)n^{-1/2}g(log n) + (1/2π)∫h(r)[Reψ(1/4+ir/2) − log π]dr`
   for `h = F²`, `g = C` the autocorrelation. **Classical, cited, not reproved here.** It
   applies to this `f`: `f` is continuous with compact support, so `h` is entire and
   `|F(r)|² = O(r^{-4})`, and all three sides converge absolutely.
2. **The archimedean reduction** (application-specific, derived here):
   `ψ(z) = −γ_E + Σ_{m≥0}[1/(m+1) − 1/(m+z)]` gives
   `Re ψ(1/4+ia/2) + γ_E = Σ_{m≥0}[1/(m+1) − 2s_m/(s_m²+a²)]`, `s_m = 2m+½`, hence
   ```
   Arch(f) = -euler C(0) + sum_{m>=0} [ C(0)/(m+1) - 2 int_0^{2L} e^{-s_m u} C(u) du ] .
   ```
   **This uses no digamma and no Lerch transcendent**, which is what makes an enclosure
   possible at all: Limitation 4 of the manuscript records that `mpmath`'s interval context
   has neither.
3. **The archimedean tail bound** (audited; see §5.2).
4. **Elimination of `γ_E`** (audited; see §5.3).

## 5. The complete error budget

Everything below is exact-rational interval arithmetic with outward rounding on a `2^-260`
dyadic grid. **There are no float operations and no libm calls in the certified path.**

### 5.1 Transcendentals, each with its proved remainder

| constant / function | method | proved remainder |
|---|---|---|
| `π` | Machin `π = 16 arctan(1/5) − 4 arctan(1/239)` | alternating series with decreasing terms: error ≤ first omitted term |
| `log 2` | `2 atanh(1/3)` | positive series, remainder ≤ first omitted term `/(1−1/9)` |
| `log 3` | `log 2 + 2 atanh(1/5)` | same |
| `log π` | `e·log2 + 2 atanh((m−1)/(m+1))`, `π = 2^e m` | same, `|(m−1)/(m+1)| ≤ 1/3` |
| `exp(t)`, `t ∈ ℚ` | `Σ_{k<n} t^k/k!` | `|R_n| ≤ |t|^n/n! · 1/(1−|t|/(n+1))`, used only with `n+1 > |t|` |
| `sin t`, `cos t`, `t ∈ ℚ` | Taylor | `|R_n| ≤ |t|^n/n!` (Lagrange, `|f^{(n)}| ≤ 1`) |
| `sin`, `cos` of an **interval** | value at the midpoint widened by the radius | `|sin'| ≤ 1`, `|cos'| ≤ 1` |
| `sqrt` | integer `isqrt` on the grid, then verified | `r² ≤ lo` and `s² ≥ hi` are **checked**, else the script raises |
| `sinh`, `cosh` | `(e^x ± e^{−x})/2` | inherited |

Self-tests against libm agree to `0.05–0.48 ulp` (a **regression check only**; the proof is
the remainder column). One apparent `25 ulp` gap resolved as a difference of *inputs*:
`cos(56/5)` versus `cos(double 11.2)`, the two arguments differing by `7.1e-16`.

### 5.2 The archimedean tail bound — audited and re-derived

```
| sum_{m>=M} t_m |  <=  [ (3/4) C(0) + (1/2) sup_{[0,2L]} |C'| ] / ( M - 3/4 ) ,
t_m = C(0)/(m+1) - 2 int_0^{2L} e^{-s_m u} C(u) du .
```

*Hypotheses, all verified:* `C` is `C¹` on `[0,2L]` (it is a finite combination of
`sin(w_k u)` and `(2L−u)cos(w_k u)`); `C(2L) = 0`; `C(0) = ‖f‖²`.
*Proof.* Integration by parts with `C(2L) = 0` gives `I_m = C(0)/s_m + s_m^{-1}∫e^{-s_m u}C'`,
so `|I_m − C(0)/s_m| ≤ sup|C'|/s_m²`. With `2C(0)/s_m = C(0)/(m+¼)`,
`|t_m| ≤ ¾C(0)/((m+1)(m+¼)) + 2sup|C'|/s_m²`. Then
`Σ_{m≥M}(m+¼)^{-2} ≤ ∫_{M−1}^∞(x+¼)^{-2}dx = 1/(M−¾)` and
`Σ_{m≥M}(2m+½)^{-2} = ¼Σ(m+¼)^{-2} ≤ ¼/(M−¾)`. ∎
*Indexing:* `m` runs from `0`; `M` is the first omitted index. **The bound is correct as
stated in the previous run.**

Numerically: `C(2L) ∈ [−8.1e-78, 8.6e-78]` (the hypothesis, verified, not assumed),
`sup|C'| ≤ 11.436722214306984`, `M = 2^17 − 1 = 131071`, giving

```
archimedean tail allowance  =  +- 4.935033775516166e-05 .
```

This single term is **the whole error budget**: everything else is below `1e-70`. It is
linear in `1/M`, so it can be reduced at linear cost; `M` was chosen as `2^17 − 1` because
`ln(M+1) = 17 log 2` exactly (§5.3).

### 5.3 Elimination of the Euler constant — audited

Rather than enclosing `γ_E`, use `γ_E = Σ_{k≥1}[1/k − ln(1+1/k)]` and the elementary
inequalities `x²/2 − x³/3 ≤ x − ln(1+x) ≤ x²/2` for `x > 0` (both proved by differentiating
`ln(1+x) − x + x²/2` and `x − x²/2 + x³/3 − ln(1+x)`, each increasing from `0`). Summing the
tail with integral bounds gives

```
H_M - euler = ln(M+1) - T ,   T in [ 1/(2(M+1)) - 1/(6 M^2) , 1/(2M) ] ,
```

and with `M + 1 = 2^17`, `ln(M+1) = 17 log 2`. So no harmonic sum and no `γ_E` enclosure is
needed, and `T` contributes a width of `3.8e-6` **inside** the `arch` term, already included
in the figure above.

### 5.4 Support cutoff — proved, not assumed

`C` vanishes outside `|u| < 2L = 8/5`, so the prime sum is over `n` with `log n < 8/5`.
The verifier **proves** `e^{8/5} ≤ 4.953033 < 5`, hence `log 5 > 2L`, hence the prime powers
are exactly `n = 2, 3, 4` with `Λ = log 2, log 3, log 2`. `log 4 = 2 log 2 < 8/5` is likewise
enclosed.

### 5.5 The assembled enclosure

```
pole   = 2 (sum_k x_k p_k)^2 ,  p_k = (-1)^k 2 w_k cosh(L/2)/(w_k^2+1/4)/sqrt(L)
       = 0.45715833231743186                (width 4.8e-77)
arch   in [ 1.2100349293664769 , 1.2101336300807928 ]   (width 9.9e-05)
logpi  = 1.1447298858493788 * ||f||^2       (width 5.9e-78)
prime  = 0.51733913650311314                (width 1.7e-76)
--------------------------------------------------------------------
Q_zeta(f) in [ 0.0051242393314170015 , 0.0052229400457328018 ]   width 9.870e-05
```

`Q_ζ(f) > 0`, consistent with Weil positivity — a negative value here would have signalled
an error, not a discovery.

### 5.6 The quartet, and how the interval is covered

`A_δ`, `B_δ` are evaluated in closed form. With `aL = (2k+1)π/2 ± 56/5` every endpoint
trigonometric value is an **exact** multiple of `sin(56/5)` or `cos(56/5)`, so only those two
transcendentals enter (plus `cosh(δL)`, `sinh(δL)`, and `π` through `a²`).

At `δ = 2/5`: `4A² ≤ 3.327e-4` (**retained, not dropped** — it works against the
certificate), `4B² ≥ 2.006640616e-2`.

For the interval the brief's requirement is met by **subdivision, not sampling**:
`[1/4, 49/100]` is split into `K = 256` closed subintervals; interval arithmetic on each
returns an enclosure valid for **every** `δ` in that subinterval, and the 256 closed
subintervals cover the whole interval with no gaps. All 256 upper bounds are negative; the
worst is on the first, `[1/4, 803/3200]`, giving `η_interval` above.

**A single interval evaluation over the whole box fails**, returning `+0.0279`: `δ` occurs in
several places in the closed forms (`δ`, `δ²`, `cosh(δL)`, `sinh(δL)`) and independent
interval occurrences inflate the result. That attempt is recorded in `certificate.json` as
`interval_single_box_attempt` — it is valid but useless, and subdivision is what closes it.

### 5.7 Margins against the error budget, separated

| statement | margin `η` | error allowance | `η` / error |
|---|---|---|---|
| point, `δ = 2/5` | `1.45108e-2` | `4.93503e-5` | **294** |
| interval-wide, worst at `δ = 1/4` | `2.45491e-3` | `4.93503e-5` | **50** |

These are the honest figures. **No `~6800×` factor is claimed**; that number came from
`δ = 1/2` against an understated allowance (see `RECONCILIATION.md` §2.4).

## 6. What is NOT claimed

- **It does not say `Q_ζ(f) < 0`.** The certificate says `Q_ζ(f) ∈ [+0.00512, +0.00523]`;
  it is the *added quartet* that makes the total negative.
- **It does not move a verified zeta zero.** `γ = 14` is a chosen rational ordinate, not a
  zero of `ζ`; the configuration adds four zeros and removes none.
- **It does not exclude an actual off-line zero.** Negativity of a Weil functional for a
  *fabricated* configuration is the detection direction of the classical criterion. No
  `L`-function is asserted to realise this spectrum.
- **No onset threshold is certified.** `0.204201571681` of the previous run is **not** a
  certified root; no root enclosure was performed.
- **No novelty** is claimed for the digamma series, the projection identity, Machin's
  formula, interval arithmetic, or validation machinery. The application-specific
  derivations actually established here are: the added-quartet functional and its
  multiplicity check (§2), the archimedean elementary-integral reduction with its proved
  tail (§4.2, §5.2), the `γ_E` elimination (§5.3), the exact endpoint reduction that makes
  the `β`'s integer multiples of `π` and the `aL`'s exact (§5.6), and the two enclosures of
  §3.
