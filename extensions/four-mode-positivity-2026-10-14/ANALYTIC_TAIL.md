# ANALYTIC_TAIL.md — the sharper archimedean error budget

This file holds the **new analytic** content of this extension: the archimedean
tail allowance of the published certificate is lowered from `O(1/M)` to
`O(1/M^2)` and then to `O(1/M^3)` **at the same truncation** `M = 2^17 - 1 =
131071`, and the `gamma_E` interval `T` is sharpened from `O(1/M^2)` to
`O(1/M^3)`. Nothing else in the certified path is touched.

Everything here is elementary: repeated integration by parts of an explicit
finite trigonometric integral, plus integral/midpoint comparison of three
convergent series. No new idea about `zeta` is involved, and none is claimed.
What the sharpening buys is stated in §6.

---

## 1. Setting (published, imported unchanged)

Let `f = sum_k x_k phi_k` with `phi_k(u) = cos(w_k u)/sqrt(L)`,
`w_k = (2k+1)pi/(2L)`, `k = 0..3`, `L = 4/5`, supported in `[-L,L]`; let
`C(u) = int f(t) f(t+u) dt` be the autocorrelation, `U = 2L = 8/5`. Put
`s_m = 2m + 1/2`. The published verifier
([`extensions/planted-witness-closure-2026-10-08/verify_certificate.py`](../planted-witness-closure-2026-10-08/verify_certificate.py),
sha256 `c49858a1…`) writes the archimedean term of the explicit formula as

```
   Arch(f) = -gamma_E C(0) + sum_{m>=0} t_m ,
   t_m     = C(0)/(m+1) - 2 I_m ,        I_m = int_0^U e^{-s_m u} C(u) du ,
```

sums `m < M` exactly, eliminates `gamma_E` through `H_M - gamma_E = ln(M+1) - T`
with `T = sum_{k>M}[1/k - log(1+1/k)]` and `ln(M+1) = 17 log 2` exactly, and
**bounds** the remainder `sum_{m>=M} t_m`. The published bound is

```
   | sum_{m>=M} t_m |  <=  [ (3/4) C(0) + (1/2) sup_{[0,U]}|C'| ] / (M - 3/4) ,
```

i.e. `±4.935033775516166e-05` at `M = 2^17-1`, **and that single term is the
entire published error budget**; the published `T` is enclosed to width
`2/(3M^2) = 3.881e-11`.

Two facts about this basis are used throughout and are **exact**, not numerical:
`w_k L = (2k+1)pi/2`, so `phi_k(±L) = 0` and the zero-extension is continuous;
and `w_k U = (2k+1)pi`, so

```
   sin(w_k U) = 0 ,        cos(w_k U) = -1 .                               (*)
```

The canonical form of `C`, also published and reused, is

```
   C(u) = sum_k [ A_k sin(w_k u) + B_k (U-u) cos(w_k u) ] ,     u in [0,U],
```
a finite trigonometric polynomial, hence **real-analytic on `[0,U]`** — so every
`C^n` hypothesis below holds with nothing to check. Its derivatives are
```
   C'(u)   = sum_k [  A_k w_k cos - B_k cos - B_k (U-u) w_k sin ] ,
   C''(u)  = sum_k [ -A_k w_k^2 sin + B_k ( 2 w_k sin - (U-u) w_k^2 cos ) ] ,
   C'''(u) = sum_k [ -A_k w_k^3 cos + B_k ( 3 w_k^2 cos + (U-u) w_k^3 sin ) ] ,
```
(arguments `w_k u` throughout), giving the sup bounds
```
   sup |C''|  <= sum_k [ |A_k| w_k^2 + |B_k| ( 2 w_k   + U w_k^2 ) ] ,
   sup |C'''| <= sum_k [ |A_k| w_k^3 + |B_k| ( 3 w_k^2 + U w_k^3 ) ] .
```

## 2. The endpoint identities

> **Lemma E.** With `(*)` and `C` even:
> ```
>   (E1)  C(U)   = 0
>   (E2)  C''(U) = 0
>   (E3)  C'(U)  = -C'(0)
>   (E4)  C'(0)  = 0            (hence C'(U) = 0 too)
>   (E5)  C''(0) = -||f'||^2 = - sum_k w_k^2 x_k^2        (exact, diagonal)
> ```

**Proof.** (E1) `C(u) = int f(t)f(t+u)dt` and `supp f = [-L,L]`, so the supports
are disjoint at `u = U = 2L`; in the canonical form, `sin(w_kU) = 0` and
`(U-u) = 0`. (E2) In `C''(U)` every term carries either `sin(w_k U) = 0` or the
factor `(U-u) = 0`. (E3) In `C'`, at `u = U` the `sin` terms and the `(U-u)`
factor drop, leaving `C'(U) = sum_k[A_k w_k cos(w_kU) - B_k cos(w_kU)] =
sum_k[-A_k w_k + B_k]` by `cos(w_kU) = -1`; at `u = 0`, `sin = 0`, `cos = 1`,
`(U-u) = U`, leaving `C'(0) = sum_k[A_k w_k - B_k]`. These are negatives.
(E4) `f` is real and even, so `C` is even; `C` is `C^1`, so `C'(0) = 0`.
(E5) At `u = 0` only the `B_k (U-u) w_k^2 cos` term survives, so
`C''(0) = -U sum_k B_k w_k^2`; the canonical reduction gives `B_k = x_k^2/U`
exactly for this basis, whence `C''(0) = - sum_k w_k^2 x_k^2`. Independently,
`C''(0) = -||f'||^2` for any `C^2` autocorrelation, and
`int phi_j' phi_k' = w_k^2 delta_{jk}` here, which is the same value. ∎

**Checked, not assumed.** The implementation evaluates all five by interval
arithmetic rather than substituting zero: `C(U) = 2.70e-79`,
`C''(U) = 5.0e-77`, `C'(U) = -7.3e-78`, `C'(0) = 3.2e-78`, and
`C''(0) = -88.000960048684` from the canonical form against
`-sum_k w_k^2 x_k^2` to the printed digits (`tests/test_tail3.py`). Where an
identity is used to *drop* a term, the interval value is carried instead and
charged.

## 3. Proposition T — the `O(1/M^2)` tail

> **Proposition T.** If `C in C^2([0,U])` and `C(U) = 0`, then
> ```
>   sum_{m>=M} t_m  =  -(3/4) C(0) S_1  -  2 C'(0) S_2  +  R ,
>   S_1 = sum_{m>=M} 1/((m+1)(m+1/4)) ,   S_2 = sum_{m>=M} s_m^{-2} ,
>   |R|  <=  sup|C''| / (4 (M-3/4)^2)  +  2 |C'(U)| sum_{m>=M} e^{-s_m U}/s_m^2 ,
> ```
> with `(4/3) log((M+1)/(M+1/4)) <= S_1 <= (4/3) log(M/(M-3/4))` and
> `(1/4)/(M+1/4) <= S_2 <= (1/4)/(M-3/4)`.

**Proof.** One integration by parts with `C(U) = 0` gives
`I_m = C(0)/s_m + s_m^{-1} J_m`, `J_m = int_0^U e^{-s_m u}C'(u)du`; a second gives
`J_m = C'(0)/s_m - e^{-s_mU}C'(U)/s_m + s_m^{-1}int_0^U e^{-s_m u}C''(u)du`, and
`|int_0^U e^{-su}C''du| <= sup|C''|/s`. So
```
   2 I_m = 2C(0)/s_m + 2C'(0)/s_m^2 - 2e^{-s_mU}C'(U)/s_m^2 + rho_m ,
   |rho_m| <= 2 sup|C''|/s_m^3 .
```
Since `s_m = 2(m+1/4)` we have `2C(0)/s_m = C(0)/(m+1/4)` and
`C(0)/(m+1) - C(0)/(m+1/4) = -(3/4)C(0)/((m+1)(m+1/4))`, giving the display.
Finally `sum_{m>=M}s_m^{-3} = (1/8)sum_{m>=M}(m+1/4)^{-3} <= 1/(16(M-3/4)^2)`,
so `2 sup|C''| sum s_m^{-3} <= sup|C''|/(8(M-3/4)^2)`; the implementation uses
the weaker constant `1/4`. The `S_1`, `S_2` enclosures are the monotone
integral comparison on `[M-1, inf)` and `[M, inf)`. ∎

By **(E4)** the `S_2` term vanishes and by **(E3)+(E4)** so does the `C'(U)`
term, leaving the allowance `sup|C''|/(4(M-3/4)^2)`. At the published frozen
vector `sup|C''| <= 111.883725`, hence `|R| <= 1.6282e-09`.

## 4. Proposition T3 — the `O(1/M^3)` tail

> **Proposition T3.** If `C in C^3([0,U])` and `C(U) = 0`, then for every `m`
> ```
>   t_m = -(3/4)C(0)/((m+1)(m+1/4)) - 2C'(0)/s_m^2 - 2C''(0)/s_m^3 + e_m ,
>   |e_m| <= 2 sup|C'''|/s_m^4 + 2 e^{-s_m U}( |C'(U)|/s_m^2 + |C''(U)|/s_m^3 ) .
> ```
> With **(E2)–(E5)** this collapses to
> ```
>   sum_{m>=M} t_m  =  -(3/4) C(0) S_1  +  2 ( sum_k w_k^2 x_k^2 ) S_3  +  R' ,
>   S_3 = sum_{m>=M} s_m^{-3} ,      |R'| <= 2 sup|C'''| sum_{m>=M} s_m^{-4} .
> ```

**Proof.** A third integration by parts of the inner integral in the proof of
Proposition T gives
`int_0^U e^{-su}C''du = C''(0)/s - e^{-sU}C''(U)/s + s^{-1}int_0^U e^{-su}C'''du`,
so
```
   I_m = C(0)/s + C'(0)/s^2 + C''(0)/s^3 - e^{-sU}[ C'(U)/s^2 + C''(U)/s^3 ]
         + s^{-3} int_0^U e^{-su} C''' du ,
```
and `|s^{-3}int_0^U e^{-su}C'''du| <= sup|C'''|/s^4`. Subtracting `2I_m` from
`C(0)/(m+1)` as before gives the first display. Then `(E2)`, `(E3)`, `(E4)` kill
the two exponential terms and the `C'(0)` term, and `(E5)` replaces `-2C''(0)/s^3`
by `+2(sum_k w_k^2 x_k^2)/s^3`, which is **exact and diagonal** in `x` — this is
why the result is still a *quadratic form* in `x` with interval coefficients, and
hence polarizable. ∎

Series enclosures used, all by the monotone integral comparison
`int_a^inf g <= sum_{m>=M} g(m+1/4) <= g(a) + int_a^inf g`, `a = M+1/4`:

```
   S_3 = (1/8) sum_{m>=M}(m+1/4)^{-3}  in  (1/8)[ 1/(2a^2) ,  a^{-3} + 1/(2a^2) ] ,
   sum_{m>=M} s_m^{-4} = (1/16) sum (m+1/4)^{-4}  <=  (1/16)[ a^{-4} + 1/(3a^3) ] .
```

`S_1` had to be sharpened too, since after the other two fixes it dominates:

> **Lemma S1.** `g(x) = 1/((x+1)(x+1/4)) = (4/3)[1/(x+1/4) - 1/(x+1)]` is
> positive, decreasing and **convex** on `[M-1/2, inf)`, because
> `g'' = (8/3)[(x+1/4)^{-3} - (x+1)^{-3}] > 0`. By Hermite–Hadamard the
> midpoint rule underestimates, `g(m) <= int_{m-1/2}^{m+1/2} g`, so
> ```
>   S_1 <= int_{M-1/2}^inf g = (4/3) log( (M+1/2)/(M-1/4) ) ,
> ```
> and the midpoint defect `int_{m-1/2}^{m+1/2}g - g(m) = g''(xi)/24 <=
> 1/(9(m-1/4)^3)` gives
> ```
>   S_1 >= (4/3) log((M+1/2)/(M-1/4)) - D ,
>   D := (1/9) sum_{m>=M}(m-1/4)^{-3} <= (1/9)[ (M-1/4)^{-3} + 1/(2(M-5/4)^2) ] .
> ```

Width `3.234e-12`, against `5.821e-11` for the monotone bound of Proposition T.
The one rational logarithm is supplied by the published `iv_log` with its proved
remainder.

## 5. The sharper `gamma_E` interval `T`

> **Lemma T.** For `0 < x <= 1`, `x - log(1+x) = sum_{n>=2}(-1)^n x^n/n` is
> alternating with decreasing terms, so
> `x^2/2 - x^3/3 <= x - log(1+x) <= x^2/2 - x^3/3 + x^4/4`. Hence with
> `P_j := sum_{k>M} k^{-j}`,
> ```
>   P_2/2 - P_3/3  <=  T  <=  P_2/2 - P_3/3 + P_4/4 ,
>   P_2 in [ 1/(M+1/2) - 1/(9M^3) ,  1/(M+1/2) ] ,
>   P_3 in [ 1/(2(M+1)^2)         ,  1/(2(M+1/2)^2) ] ,
>   P_4 in [ 1/(3(M+1)^3)         ,  1/(3M^3) ] .
> ```

**Proof.** `1/k^2 <= int_{k-1/2}^{k+1/2}x^{-2}dx = 1/(k^2-1/4)` gives
`P_2 <= int_{M+1/2}^inf x^{-2} = 1/(M+1/2)`, and the defect
`sum (1/4)/(k^2(k^2-1/4)) <= (1/3)P_4 <= 1/(9M^3)` since `k^2-1/4 >= (3/4)k^2`.
`x^{-3}` and `x^{-4}` are convex, so the midpoint rule underestimates the
integral, giving the upper bounds `1/(2(M+1/2)^2)` and `1/(3M^3)`; the left
comparison `1/k^j >= int_k^{k+1}x^{-j}dx` gives the lower bounds. ∎

Width `1.357e-16`, against the published `3.881e-11` — a factor `2.86e+05` —
and the new interval is **contained in** the published one.

## 6. Measured effect, and the only check that matters

At the published frozen vector `x = (99815221908, -354082045264, 928280739889,
54384691304)/10^12`:

| | published `O(1/M)` | Prop. T `O(1/M^2)` | Prop. T3 `O(1/M^3)` |
|---|---|---|---|
| tail allowance (half-width) | `4.935033776e-05` | `1.6500e-09` | `~2e-14` |
| `Q_zeta` enclosure width | `9.870e-05` | `3.339e-09` | `2.477e-12` |
| `T` width | — | `3.881e-11` (published) | `1.357e-16` |

**Each enclosure is contained in the previous one.** That is the decisive test:
three valid enclosures of the same real number, each inside the last. It is
asserted in `verify_route_o.py` (O1) and in `tests/test_tail_sharp.py`,
`tests/test_tail3.py`, and it would fail loudly on a sign error or a dropped
term. The `O(1/M^2)` budget is what makes the detection certificate of
[`CERTIFICATE_THEOREM.md`](CERTIFICATE_THEOREM.md) possible at all; the
`O(1/M^3)` budget is what makes
[`FOUR_MODE_POSITIVITY.md`](FOUR_MODE_POSITIVITY.md) possible, since the
quantity certified there is `1.75e-10`.

## 7. Scope

* This is **not** a precision increase. There is no floating-point arithmetic
  anywhere in this code; `M` is unchanged at `2^17-1` and not one extra term is
  summed. It is a strictly better analytic bound at the same cost.
* It changes **nothing** else in the explicit formula: the pole term, the prime
  sum (support cutoff `log 5 > 2L`, so exactly `n = 2,3,4`), the `log pi` term,
  the `gamma_E` elimination identity and the exact quartet closed forms are the
  published ones, imported.
* `(E5)` and the collapse in Proposition T3 are **specific to this basis**
  (`w_k U in pi Z`). For a basis without that property the `C'(U)`, `C''(U)`
  terms survive and must be carried; the propositions as stated above still
  hold, with those terms non-zero.
* Further sharpening would need a fourth integration by parts or a larger `M`.
  Neither is needed here, and neither is the binding constraint on any result in
  this extension — see [`NON_DETECTION.md`](NON_DETECTION.md) §4.
