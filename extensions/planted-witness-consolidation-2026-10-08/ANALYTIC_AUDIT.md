# Analytic audit

Everything the arithmetic rests on, proved or explicitly cited. Labels: **[P]** proved here,
**[C]** classical, cited; **[V]** verified numerically as a consistency check only.

Throughout: `L = 4/5`, `U = 2L = 8/5`, `γ = 14`, all exactly rational;
`φ_k(x) = cos(w_k x)/√L` with `w_k = (2k+1)π/(2L) = (2k+1)(5/8)π`, `k = 0,1,2,3`;
`x = (99815221908, −354082045264, 928280739889, 54384691304)/10¹²`; `f = Σ x_k φ_k`,
extended by `0` outside `[−L, L]`.

## 1. The basis, its endpoints and the zero extension **[P]**

`w_k L = (2k+1)π/2`, so `cos(w_k L) = 0` **exactly**. Hence `φ_k(±L) = 0` and the
zero-extended `φ_k` is **continuous on ℝ**, real, even, supported in `[−L,L]`.
`φ_k'(x) = −w_k sin(w_k x)/√L` with `sin(w_k L) = (−1)^k ≠ 0`, so `φ_k'` has a jump at `±L`:
`f ∈ C⁰(ℝ)` with `f'` of bounded variation, and `f ∉ C¹(ℝ)`. That regularity class is what
§3 and §4 use; nothing here assumes more.

Orthonormality: `∫_{−L}^{L} cos(w_j x)cos(w_k x)dx = L δ_{jk}` because
`w_j ± w_k = (j±k+…)·(5/4)π` makes both `sin((w_j−w_k)L)` and `sin((w_j+w_k)L)` vanish
(`(w_j ± w_k)L` is an integer multiple of `π`), and the diagonal gives `L + sin(2w_kL)/(2w_k) = L`.
So `‖f‖₂² = Σ_k x_k²` **exactly**, with the Gram matrix equal to the identity:

```
||f||_2^2 = 999999999999981276942897 / 10^24 = 0.999999999999981276942897   (exact)
```

**The witness is not renormalised.** Everything below is stated for this `f`. Under
`f ↦ cf` with `c > 0` the functional scales by `c²`, so the *sign* is invariant but the
*values* are not; a unit-normalised version would have all form values multiplied by
`1/‖f‖₂² = 1 + 1.87×10^{−14}`, which is recorded and not used.

## 2. The autocorrelation, and `C(2L) = 0` **exactly** **[P]**

`C(u) := ∫_ℝ f(x) f(x−u) dx` (the convention of `src/zeta_window.py`, where `C_{jk}` is the
symmetrised `∫ φ_j(x)φ_k(x−u)dx`).

> **Lemma 1 (support).** `C(u) = 0` for `|u| ≥ 2L`; in particular `C(2L) = 0` exactly.
>
> *Proof.* `f(x) ≠ 0` requires `x ∈ [−L,L]` and `f(x−u) ≠ 0` requires `x ∈ [u−L, u+L]`. For
> `u = 2L` the two conditions give `x ∈ [−L,L] ∩ [L,3L] = {L}`, a set of Lebesgue measure
> zero, and the integrand is bounded; so the integral is `0`. For `u > 2L` the intersection
> is empty. ∎

This is an identity, not an estimate. The previously reported `|C(2L)| < 9×10^{−78}` is a
**[V]** consistency check on the closed form, and is *not* what justifies anything.

**Canonical form.** Using `C_terms` of the source with every `β` an integer multiple of `π`
(so `sin β = 0`, `cos β = ±1` exactly — see §2.1), `C` reduces on `[0, 2L]` to

```
C(u) = sum_{k=0}^{3} [ A_k sin(w_k u) + B_k (U - u) cos(w_k u) ] .
```

Two identities follow immediately and are used as structural checks:
`C(0) = U Σ_k B_k = ‖f‖₂²` **[P,V]** and `C(U) = Σ_k A_k sin(w_k U) = 0` **[P]**, since
`w_k U = (2k+1)π`. The second is Lemma 1 again, now visible in the closed form.

### 2.1 Why every `β` is an integer multiple of `π` **[P]**

In `C_terms(j,k)` the phases are `β = ±p_1 L` and `β = ±p_2 L` with `p_1 = w_j − w_k`,
`p_2 = w_j + w_k`. Then `p_1 L = (j−k)π` and `p_2 L = (j+k+1)π`. Hence `sin β = 0` and
`cos β = (−1)^{j−k}`, `(−1)^{j+k+1}` **exactly**. This is special to `L = 4/5` with this
basis and is what makes the reduction exact rather than approximate.

### 2.2 Regularity for the integration by parts **[P]**

On `[0, 2L]` the closed form shows `C` is real-analytic: a finite combination of
`sin(w_k u)` and `(U−u)cos(w_k u)`. So `C ∈ C^∞([0,2L])`, `C(0) = ‖f‖²`, `C(2L) = 0`, and

```
int_0^U e^{-s u} C(u) du = C(0)/s - e^{-sU} C(U)/s + (1/s) int_0^U e^{-su} C'(u) du
                         = C(0)/s + (1/s) int_0^U e^{-su} C'(u) du ,
```

the boundary term at `U` vanishing **by Lemma 1**, not by a bound. Valid for every `s > 0`.

### 2.3 A proved bound for `sup|C'|` — not a sampled maximum **[P]**

Differentiating the canonical form,

```
C'(u) = sum_k [ A_k w_k cos(w_k u) - B_k cos(w_k u) - B_k (U-u) w_k sin(w_k u) ] ,
```

so for every `u ∈ [0,U]`, using `|cos|,|sin| ≤ 1` and `0 ≤ U−u ≤ U`,

```
sup_{[0,U]} |C'|  <=  sum_k ( |A_k| w_k + |B_k| + |B_k| U w_k ) .
```

This is the triangle inequality applied to the closed form — a **proof**, valid at every
point, with no sampling. The verifier evaluates the right-hand side in interval arithmetic
and obtains `sup|C'| ≤ 11.436722214306984`.

## 3. The explicit formula: version, hypotheses, applicability

**Version used [C].** Weil's explicit formula for `ζ` in the `ψ`-form: for `h` even,
holomorphic in `|Im r| ≤ 1/2 + ε` and `O((1+|r|)^{-1-η})` there, with
`g(u) = (1/2π)∫h(r)e^{-iru}dr`,

```
sum_rho h(gamma_rho) = h(i/2) + h(-i/2) - 2 sum_{n>=2} Lambda(n) n^{-1/2} g(log n)
                       + (1/2pi) int h(r) [ Re psi(1/4 + i r/2) - log pi ] dr .
```

Primary references as used by this project: Bombieri, *Remarks on Weil's quadratic
functional in the theory of prime numbers, I* (ref. `[Bombieri]` of the note) for the
quadratic functional, and Zhu (ref. `[Zhu]`) for the geometric-side assembly on a finite
window, which `src/zeta_window.py` implements. **The formula is cited, not reproved**, and
it is the only non-elementary input to the certificate.

**Applicability to this `f` [P].** Take `h(r) = F(r)²` with `F(z) = ∫_{−L}^{L}f(u)e^{izu}du`.
Since `f ∈ L¹` has compact support, `F` is **entire**, and `h` is entire and even
(`F` is even because `f` is). Decay: integrating by parts once, using `f(±L) = 0`,

```
F(r) = -(1/(ir)) int_{-L}^{L} f'(u) e^{iru} du ,
```

and once more, with `f'` of bounded variation and `f'(±L^∓) ≠ 0`,

```
F(r) = (1/(ir)^2) [ f'(u) e^{iru} ]_{-L}^{L} - (1/(ir)^2) int f'' e^{iru}  =  O(r^{-2}) ,
```

so `h(r) = O(r^{-4})` on `ℝ`. In the strip `|Im r| ≤ 1/2 + ε` the same two integrations by
parts give `|F(σ+iτ)| ≤ C e^{|τ|L}(1+|σ|)^{-2}`, so `h = O((1+|r|)^{-4})` uniformly there.
Both hypotheses hold with `η = 3`. Consequently all three sides converge absolutely:
`Σ_ρ|F(γ_ρ)|²` converges because the zero-counting density is `O(log T)` and
`|F|² = O(T^{-4})`; the prime sum is **finite** (§3.1); and
`∫|F|²|Reψ(1/4+ir/2)|dr` converges because `Reψ = O(log(2+|r|))`.
No approximation or continuity argument is needed — `f` is directly admissible — and **no
statement about zero locations is used anywhere**.

### 3.1 Prime-power term and support cutoff **[P]**

`g(log n) = C(log n)` and `C` vanishes for `|u| ≥ 2L = 8/5` (Lemma 1), so the prime sum is
over `n` with `log n < 8/5`. The verifier **proves** `e^{8/5} ≤ 4.953033 < 5`, hence
`log 5 > 8/5`, and `log 4 = 2 log 2 < 8/5` is likewise enclosed. So

```
prime = 2 [ log2 * 2^{-1/2} C(log2) + log3 * 3^{-1/2} C(log3) + log2 * 2^{-1} C(log4) ] ,
```

`Λ(2)=Λ(4)=log 2`, `Λ(3)=log 3`, and no other `n` contributes. **The sum is finite and
complete**, not truncated.

### 3.2 Pole and `log π` terms **[P]**

`h(i/2) + h(−i/2) = 2F(i/2)²` because `F(−i/2) = F(i/2)` for even `f`. With
`F_k(i/2) = (−1)^k 2 w_k \cosh(L/2)/(w_k^2 + 1/4)` (the source's `F_at_half` for `even_d`),
divided by `√L` for the orthonormal mode,

```
pole = 2 ( sum_k x_k p_k )^2 ,   p_k = (-1)^k 2 w_k cosh(L/2) / ((w_k^2 + 1/4) sqrt(L)) .
```

`(1/2π)∫h(r)dr = ‖f‖₂²` by Plancherel, so the `log π` term is `log π · ‖f‖₂²`. Both match
`src/zeta_window.py` term for term.

## 4. The archimedean term as elementary integrals **[P]**

From `ψ(z) = −γ_E + Σ_{m≥0}[1/(m+1) − 1/(m+z)]` **[C]**, taking `z = 1/4 + ia/2` and real
parts, and using `Re(m+1/4+ia/2)^{-1} = (m+1/4)/((m+1/4)²+(a/2)²) = 2s_m/(s_m²+a²)` with
`s_m = 2m+1/2`,

```
Re psi(1/4 + i a/2) + euler = sum_{m>=0} [ 1/(m+1) - 2 s_m/(s_m^2 + a^2) ] .
```

Since `2s_m/(s_m²+a²) = 2∫_0^∞ e^{-s_m t}\cos(at)\,dt` and
`(1/2π)∫h(r)\cos(rt)dr = C(t)`, interchanging the (absolutely convergent) `t`-integral with
the `r`-integral gives

```
Arch(f) = -euler C(0) + sum_{m>=0} [ C(0)/(m+1) - 2 int_0^{2L} e^{-s_m u} C(u) du ] ,
```

the upper limit being `2L` by Lemma 1. **No digamma and no Lerch transcendent appear**, which
is what makes an enclosure possible: Limitation 4 of the note records that `mpmath`'s
interval context has neither.

### 4.1 The closed forms, and large `s_m` **[P]**

With `E = e^{-s U}`, `w_k U = (2k+1)π` so `\cos(w_kU) = −1`, `\sin(w_kU) = 0` exactly, and

```
int_0^U e^{-su} sin(w u) du   = w (1 + E)/(s^2 + w^2) ,
int_0^U e^{-su} (U-u) cos(w u) du = sU/(s^2+w^2) - (1+E)(s^2-w^2)/(s^2+w^2)^2 .
```

Both are exact for every `s > 0`; no cancellation is amplified, because the only subtraction
is `(s²−w²)` with `s², w²` of comparable size only for the single index with `w_k ≈ s_m`, and
there the enclosure is still exact rational. **Large `s_m` does not invalidate either
formula**: `E ∈ (0,1)` and the denominators grow. For `m ≥ 20`, `s_m U ≥ 64.8 > 90\log 2`
fails — so the verifier uses the proved enclosure `E ∈ [0, 2^{-90}]` only for `m ≥ 20`,
where `s_m U ≥ 64.8` and `e^{-64.8} < 2^{-93} < 2^{-90}`; for `m < 20` it encloses `E` by the
`exp` series. **Checked:** `64.8/\log 2 = 93.5 > 90`.

### 4.2 The tail bound, audited **[P]**

```
t_m := C(0)/(m+1) - 2 int_0^{2L} e^{-s_m u} C(u) du ,
| sum_{m >= M} t_m |  <=  [ (3/4) C(0) + (1/2) sup|C'| ] / ( M - 3/4 ) .
```

*Proof.* By §2.2, `|I_m − C(0)/s_m| ≤ sup|C'|/s_m²`. With `2C(0)/s_m = C(0)/(m+1/4)`,
`t_m = C(0)[1/(m+1) − 1/(m+1/4)] − 2(I_m − C(0)/s_m)`, so
`|t_m| ≤ (3/4)C(0)/((m+1)(m+1/4)) + 2\sup|C'|/s_m²`. Then, since `(m+1)(m+1/4) ≥ (m+1/4)²`,
`Σ_{m≥M}(m+1/4)^{-2} ≤ ∫_{M−1}^{∞}(t+1/4)^{-2}dt = 1/(M−3/4)`, and
`Σ_{m≥M}(2m+1/2)^{-2} = (1/4)Σ_{m≥M}(m+1/4)^{-2} ≤ (1/4)/(M−3/4)`. ∎

*Indexing:* `m` starts at `0`; `M` is the first **omitted** index; the sum actually evaluated
is `m = 0 … M−1`. With `M = 2^17 − 1 = 131071`, `C(0) < 1` and `\sup|C'| ≤ 11.436722214…`,

```
| tail |  <=  ( 0.75 + 5.718361... ) / 131070.25  =  4.935033775516166e-05 .
```

Signs: the bound is two-sided, and the verifier adds the symmetric interval
`[−4.9350…e−5, +4.9350…e−5]`.

### 4.3 The harmonic-number / Euler-constant cancellation, audited **[P]**

`γ_E = Σ_{k≥1}[1/k − \ln(1+1/k)]`, and for `x > 0`

```
x^2/2 - x^3/3  <=  x - log(1+x)  <=  x^2/2 ,
```

both by differentiating (`ln(1+x) − x + x²/2` has derivative `x²/(1+x) ≥ 0`;
`x − x²/2 + x³/3 − \ln(1+x)` has derivative `x³/(1+x) ≥ 0`, each vanishing at `0`).
Summing from `k = N` with `∫`-bounds (`Σ_{k≥N}k^{-2} ∈ [1/N, 1/(N−1)]`,
`Σ_{k≥N}k^{-3} ≤ 1/(2(N−1)^2)`) and telescoping the logs
(`Σ_{k=1}^{N−1}[1/k − \ln\frac{k+1}{k}] = H_{N−1} − \ln N`) gives

```
euler = H_{N-1} - log N + T_N ,   T_N in [ 1/(2N) - 1/(6(N-1)^2) , 1/(2(N-1)) ] .
```

Taking `N = M+1` (so `N−1 = M`): `H_M − γ_E = \ln(M+1) − T_{M+1}` with

```
T_{M+1} in [ 1/(2(M+1)) - 1/(6 M^2) , 1/(2M) ] ,
```

and `M + 1 = 2^17` makes `\ln(M+1) = 17\log 2` exactly. So the archimedean sum is evaluated as

```
Arch = C(0) ( 17 log2 - T ) - 2 sum_{m<M} I_m + tail ,
```

**no harmonic number is summed and `γ_E` is never enclosed.** The `T`-width is
`1/(2M) − 1/(2(M+1)) + 1/(6M²) ≈ 3.9×10^{−11}`, inside the figure of §4.2.

## 5. The independent route (used only by `independent_point_check.py`) **[P]**

### 5.1 A closed form for the same series

With `z_k = 1/4 − i w_k/2`, `s_m = 2(m+1/4)` and `|m+z_k|² = (m+1/4)² + (w_k/2)²`:

```
w/(s^2+w^2) = (1/2) Im (m+z)^{-1},   s/(s^2+w^2) = (1/2) Re (m+z)^{-1},
(s^2-w^2)/(s^2+w^2)^2 = (1/4) Re (m+z)^{-2} .
```

Splitting `I_m = I_m^{(0)} + E_m I_m^{(1)}`, using `C(0) = U Σ_k B_k`, and
`Σ_m[1/(m+1) − (m+z)^{-1}] = ψ(z)+γ_E`, `Σ_m (m+z)^{-2} = ψ'(z)`, **the `γ_E` cancels
identically** and

```
Arch(f) = sum_k [ U B_k Re psi(z_k) + A_k Im psi(z_k) + (B_k/2) Re psi'(z_k) ]
          - 2 sum_{m>=0} E_m I_m^{(1)} ,
I_m^{(1)} = sum_k [ A_k w_k/(s_m^2+w_k^2) - B_k (s_m^2-w_k^2)/(s_m^2+w_k^2)^2 ] .
```

The `E`-series converges geometrically (`E_{m+1}/E_m = e^{-16U/5} = e^{-5.12}`), so `24`
terms plus a geometric bound suffice. **Different formula, different error budget, no `1/M`
tail.**

### 5.2 Rigorous digamma

`ψ(z) = ψ(z+J) − Σ_{j<J}(z+j)^{-1}` (exact recurrence **[C]**), and

```
psi(w) = log w - sum_{j>=0} g(w+j) ,   g(x) = 1/x - log(1+1/x) = sum_{n>=2} (-1)^n/(n x^n) ,
```

proved by telescoping: `Σ_{j=0}^{J'-1}g(w+j) = Σ(w+j)^{-1} − \ln\frac{w+J'}{w}` and
`Σ_{j<J'}(w+j)^{-1} = ψ(w+J') − ψ(w) \to \ln(w+J') − ψ(w)`. Tail:
`|Σ_{j≥J'}g(w+j)| ≤ (1/2)(1−1/R_w)^{-1}/(R_w+J'−1)` with `R_w = \mathrm{Re}\,w`. With
`J = 50`, `J' = 20000` the tail is `≤ 2.50×10^{−5}`; each `g` is summed to `n = 8` with the
proved remainder `1/((n+1)m^{n+1})·m/(m−1)`, `m` an integer lower bound for `|x|`.

### 5.3 Rigorous trigamma

`ψ'(z) = Σ_{j<J}(z+j)^{-2} + ψ'(z+J)`, and for `f(t) = (w+t)^{-2}` the midpoint rule on unit
cells gives `Σ_{j≥0}f(j) = ∫_{-1/2}^{∞}f + E` with
`|E| ≤ (1/24)Σ_j \sup_{\text{cell}}|f''| ≤ (1/24)·2/(R_w−3/2)^3`, hence

```
psi'(w) = 1/(w - 1/2) + E ,   |E| <= 1/(12 (Re w - 3/2)^3) ,
```

which at `R_w = 50.25` is `≤ 7.2×10^{−7}`. Complex `\log w` uses
`\tfrac12\log|w|^2 + i\arctan(\mathrm{Im}\,w/\mathrm{Re}\,w)` with `|\mathrm{Im}/\mathrm{Re}| ≤ 0.137`,
inside the proved `\arctan` series range.

## 6. Error contributions — none is literally zero

| contribution | bound | mechanism |
|---|---|---|
| archimedean tail, `M = 2^17−1` | `4.935034e-05` | §4.2 |
| `T` (Euler cancellation) | `≈ 3.9e-11` | §4.3, inside the above |
| grid rounding, `2^-260`, per operation | `≤ 2^{-260}` outward | exact `Fraction` endpoints |
| `π` (Machin, 60 terms) | `< 1.1e-77` | alternating remainder |
| `\log 2`, `\log 3`, `\log π` | `< 5e-78` | `atanh` series remainder |
| `\exp`, `\sin`, `\cos` | `< 1e-77` | Lagrange / ratio-test remainder |
| `\sqrt{}` | brackets **verified** | `isqrt`, checked or raise |
| node residual `F(14)` | `4.1293580707791e-13`, enclosed | carried, never set to `0` |
| `Q_add − Q_rep = 2F(14)²` | `3.4103196e-25`, enclosed | carried |

Everything other than the archimedean tail is below `10^{-70}`, but **it is included, not
discarded**: the verifier propagates all of it.

## 7. What is cited, and what is not claimed

**Cited [C]:** Weil's explicit formula (§3), the digamma series and recurrence (§4, §5.2),
Machin's formula, the classical Weil positivity criterion. **No novelty** is claimed for any
of these, nor for interval arithmetic, nor for the elementary projection identity behind the
witness.

**Established here [P]:** Lemma 1 and the exact `C(2L) = 0`; the integer-multiple-of-`π`
phase reduction special to `L = 4/5`; the proved `\sup|C'|` bound; the elementary-integral
archimedean representation with its audited tail; the `γ_E` elimination; the `ψ/ψ'` closed
form of §5.1 with its `γ_E` cancellation; and the two enclosures of the certificate.

**Not claimed:** no statement about actual zeta-zero locations, no new RH criterion, no
optimal threshold, no all-height theorem, and no Lean-style formal verification — Python's
integer/`Fraction` implementation and the interpreter remain ordinary computational
dependencies.
