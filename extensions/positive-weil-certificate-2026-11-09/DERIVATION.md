# DERIVATION.md

Normalizations, the reference inequality, the gate actually proved here, and
exactly which finite-to-infinite argument is and is not used.

---

## 1. Normalization of the window Weil form

Taken unchanged from `src/zeta_window.py` and `notes/stage3_assembly.md`
(status IMPORTED, W4). Even sector, `f_k(x) = cos(k pi x/L)` on `[-L,L]` extended
by zero, `int f_j f_k = L delta_{jk}` (`2L` for `k=0`). In the **orthonormal**
basis the form is

```
    M  =  Pole + Arch - LogPi - Prime
```

with, writing `C_{jk}(u) = int f_j(x) f_k(x-u) dx` symmetrised and
`P_k = int f_k(x) e^{-x/2} dx = F_k(i/2)`,

```
    Pole_{jk}  = 2 P_j P_k                                  (sign +1 in even sectors)
    Arch_{jk}  = -(gamma + log 4) C(0)
                 + 2 sum_{m>=0}[ C(0)/(2m+1) - int_0^{2L} C(u) e^{-(2m+1/2)u} du ]
    LogPi_{jj} = log(pi) C(0)
    Prime_{jk} = sum_{n : log n < 2L} (2 Lambda(n)/sqrt n) C_{jk}(log n)
```

(The grouping matters: the terms are `O(1/m^2)`. The variant
`sum_{m>=0}[C(0)/(m+1) - 2 int C e^{-(m+1/2)u}]` printed in earlier revisions of
this file **diverges**, its `m`-th term behaving like `-C(0)/m`; the constant
there also dropped the `log 4` of `int_0^inf (1-e^{-s})/sinh s ds = 2 log 2`.
Nothing computed ever used that expression: the pipeline uses the frequency form
below, via `baseline.H`, and the time-domain evaluator uses
`zeta_window.arch_closed`.)

Two facts matter for error accounting and are easy to lose:

* the **prime sum is finite and exact** — `C_{jk}` is supported in `|u| <= 2L`,
  so only `n` with `log n < 2L` contribute, with no truncation whatsoever;
* the **archimedean sum is in closed form** via digamma, Lerch `Phi` and Hurwitz
  `zeta` (`zeta_window.arch_closed`), so it carries no series-truncation error either.

Therefore defect **(i)** of the brief — analytic truncation when assembling an
entry — is *absent* in this normalization once the special functions are
evaluated in ball arithmetic. Defect **(ii)** — the infinite-dimensional
complement — is untouched by that, and is what remains open. **Bounding (i) does
not discharge (ii), and here (i) was never the obstacle.**

In the frequency form used by the reduction (§3), the same object reads

```
    Q(f) = 2 F(i/2)^2 + (1/2pi) int_R |F(t)|^2 Psi_L(t) dt ,
    Psi_L(t) = Re psi(1/4 + it/2) - log pi - P_L(t) ,
    P_L(t)   = sum_{log n < 2L} (2 Lambda(n)/sqrt n) cos(t log n) .
```

For real even `f`, `F` is real and even, so `(1/2pi)int_R = (1/pi)int_0^infty`,
and `(1/pi) int_0^infty |F|^2 = ||f||_2^2` by Plancherel. Those two factors of 2
are the ones that are easiest to drop; they are used explicitly in §3.

## 2. The reference detector (Phase A)

**Polynomial.** With `c = cos theta` and `cos 2theta = 2c^2 - 1`,

```
    P(theta) = 3 + 4 cos theta + cos 2theta = 2c^2 + 4c + 2 = 2(1 + cos theta)^2 .
```

Verified as an identity in `Q[c]` by exact Chebyshev expansion (`wcore.cheb_to_poly`),
not by sampling. Nonnegativity is then structural: twice a square.

**Inequality.** For `sigma > 1`, `log zeta(s) = sum_n c_n n^{-s}` with
`c_{p^k} = 1/k >= 0`, absolutely convergent, so

```
    log| zeta(sigma)^3 zeta(sigma+it)^4 zeta(sigma+2it) |
        = sum_n c_n n^{-sigma} ( 3 + 4 cos(t log n) + cos(2t log n) )
        = sum_n c_n n^{-sigma} P(t log n)  >=  0 ,
```

using `Re log zeta(sigma+it) = sum_n c_n n^{-sigma} cos(t log n)`. Exponentiating,

```
    zeta(sigma)^3 |zeta(sigma+it)|^4 |zeta(sigma+2it)|  >=  1 ,   sigma > 1.   (*)
```

Only K1–K3 enter. No continuation, no functional equation.

**Boundary theorem (no zeros on `Re s = 1`).** Here the wrapper enters, and is
declared. Suppose `zeta(1 + i t_0) = 0` with `t_0 != 0`. As `sigma -> 1+`:

* `zeta(sigma)^3 = O((sigma-1)^{-3})` — **W2**, the pole is simple;
* `|zeta(sigma + i t_0)|^4 = O((sigma-1)^4)` — **W1**, `zeta` is analytic at
  `1+it_0`, so a zero there is of order `>= 1`;
* `|zeta(sigma + 2 i t_0)| = O(1)` — **W3**, `zeta` is regular at `1 + 2it_0`
  because `2t_0 != 0`.

The product is `O(sigma - 1) -> 0`, contradicting `(*)`. Hence `zeta(1+it) != 0`
for `t != 0`. **PROVED**, with W1–W3 named.

**Not claimed.** No quantitative zero-free region is derived. Reproducing a
classical `1 - c/log t` region with stated constants was not attempted, and no
numerical sampling anywhere is presented as a zero-free-region proof.

## 3. The one-stroke reduction and the gate proved here

**Imported (R1, Zhu Thm 1.1).** If `Psi_L(t) >= beta` for all `t >= T#`, then for
real even `f` supported in `[-L,L]`,

```
    Q(f) >= 2F(i/2)^2 + (1/pi) int_0^{T#} (Psi_L(t) - beta) |F(t)|^2 dt + beta ||f||_2^2 .
```

The derivation is the split `(1/pi)int_0^infty Psi_L|F|^2 = (1/pi)int_0^{T#}(Psi_L-beta)|F|^2
+ (1/pi)int_{T#}^infty (Psi_L-beta)|F|^2 + beta ||f||^2`, discarding the middle
term by the hypothesis. The paper then represents the right side in a Legendre
basis as `beta I + 2 p p^T + C` and uses PSD-ness of a leading block plus
super-exponentially small tails.

**The gate.** The paper takes `beta_paper(T#) = log(T#/2pi) - 1/T# - A_L`, from a
lower bound on the digamma, with `A_L = sum_{log n < 2L} 2 Lambda(n)/sqrt n` and
`sup_t P_L = A_L` (R2).

**PROVED HERE (G1–G3).** Use instead the exact series (W5)

```
    Re psi(1/4 + it/2) = -gamma + sum_{n>=0} [ 1/(n+1) - (n+1/4)/((n+1/4)^2 + t^2/4) ] .
```

Each summand is strictly increasing in `t > 0`, because

```
    d/dt [ -(n+1/4)/((n+1/4)^2 + t^2/4) ] = (n+1/4)(t/2) / ((n+1/4)^2 + t^2/4)^2 > 0 .
```

So `t -> Re psi(1/4+it/2)` is strictly increasing on `t > 0`. Combining with
`P_L <= A_L`,

```
    Psi_L(t) >= beta_sharp(T#) := Re psi(1/4 + i T#/2) - log pi - A_L   for all t >= T# .
```

`beta_sharp >= beta_paper` pointwise, since the paper's expression is a lower
bound for the same digamma. Define the **exact gate**

```
    T_gate(L) := the unique t > 0 with Re psi(1/4 + it/2) = log pi + A_L ,
```

the least admissible `T#`; it is located by certified bisection, legitimate
precisely because of the strict monotonicity just proved.

Measured (Arb, 400 bits):

| `L` | `A_L` | `T_1 = 2pi e^{A_L}` | `T_gate` (sharp) | `T_paper` (paper's root) |
|---|---|---|---|---|
| `4/5` | `2.94197352522` | `119.086556188` | `119.086906076` | `120.082403914` |
| `81/100` | `4.38149855633` | `502.389453522` | `502.389536459` | `503.388460911` |

Note `T_1` is **not** itself admissible: `beta_paper(T_1) = -1/T_1 < 0`, because
`T_1` is the root of `log(T/2pi) - A_L`, which drops the `-1/T`. The ordering is
`T_1 < T_gate < T_paper`, so the sharpened gate admits a slightly smaller window
of integration (about `0.83%` at `L = 4/5`, `0.20%` at `L = 81/100`).

**This is a constant, not a threshold.** `T_gate` and `T_paper` agree to leading
order; both are `2pi e^{A_L}(1+o(1))`. The sharpening does not touch the
doubly exponential growth, and does not lower `A_L`, which R2 says is optimal
for any *pointwise* comb envelope.

## 4. What is NOT done

**R4 is not reproduced.** The actual `L = 0.8` lower bound needs the leading
`N x N` Legendre/Bessel block with `N ~ e L T#/2 ~ 217` at 50 digits, plus the
`eps_D, eps_B` tail constants. That assembly was not performed here (SOURCE_AUDIT §3),
so **no lower bound on `lambda_min(L)` is claimed at any `L`, and whole-window
status is UNRESOLVED throughout.**

What *is* done in that direction is the opposite inequality, which needs nothing:

```
    lambda_min(L) = inf_{f} Q(f)/||f||^2  <=  Q(f_x)/||f_x||^2    for any single admissible f_x,
```

so one frozen rational vector gives a rigorous **upper** bound. Equivalently, for
any finite subspace `V_N`, `lambda_min(L) <= lambda_min(M_N)`. A **negative**
`lambda_min(M_N)` would *refute* whole-window positivity; a positive one proves
nothing about the complement. None was found negative.

---

## 5. Phase 2: the finite-to-infinite comparison, derived in this normalization

Not imported. The imported object is only the *idea* of a one-stroke reduction;
every step below is rederived in the repository's own convention and the
constants are certified in `phase2/src/baseline.py`.

**Setting.** `F(t) = int_{-L}^L f(u)e^{itu}du`, `int_R|F|^2 = 2pi||f||_2^2`,

```
    Q(f) = 2 s F(i/2)^2 + (1/2pi) int_R |F(t)|^2 Psi_L(t) dt ,
    Psi_L = H - P_L ,  H(t) = Re psi(1/4+it/2) - log pi ,
    s = +1 even sectors, s = -1 odd sector  (zeta_window.pole_sign).
```

**Step 1, the gate.** `H` is strictly increasing on `t>0` (§3, G1–G3) and even,
and `P_L <= A_L`. Put `beta := H(T) - A_L`. Then for `|t| >= T`,
`Psi_L(t) = H(t) - P_L(t) >= H(T) - A_L = beta`. At `L=4/5, T=200`:
`beta = 0.51846573324385…> 0`, and **the unfavourable exterior set is empty**.

**Step 2, the split.** Since `int_{|t|>=T}(Psi_L-beta)|F|^2 >= 0`,

```
    Q(f) >= R_T(f) := beta ||f||_2^2 + 2 s F(i/2)^2
                      + (1/2pi) int_{-T}^{T} (Psi_L - beta)|F|^2 dt .
```

**Step 3, the basis.** Legendre, `phi_n(u) = sqrt((2n+1)/(2L)) P_n(u/L)`, so
`F_n(t) = i^n sqrt(2L(2n+1)) j_n(tL)` and `sum_{n>=0}(2n+1) j_n(x)^2 = 1`. This
basis is chosen because `j_n(x)` is super-exponentially small for `n > e x/2`,
which is what makes Step 4 effective; the cosine basis of `src/zeta_window.py`
gives only `O(1/n^2)` and cannot reach the required `1e-20`.

**Step 4, the complement.** For `f` in `V_N^perp` with `||f||=1`, Cauchy–Schwarz
in `n` gives `|F(t)|^2 <= 2L S_N(tL)` with `S_N(x) = sum_{n>=N}(2n+1)j_n(x)^2`,
and `S_N` is increasing on `0 <= x < N`, so

```
    g_N := sup_{V_N^perp} (1/2pi) int_{-T}^{T}|F|^2  <=  (2LT/pi) S_N(TL) .
```

This is the Cauchy–Schwarz surrogate for the `(N+1)`-st Slepian eigenvalue:
lossy, rigorous, and super-exponentially small past the Shannon number
`2TL/pi = 101.859…`.

**Step 5, the block bound.** With `M = sup_{|t|<=T}|Psi_L-beta|`,
`a = ||x_N||`, `b = ||x_perp||`, and `s = +1` so the pole term is `>= 0`:

```
    <R_T x_perp, x_perp>  >= (beta - M g_N) b^2 ,
    |cross term|          <= 2 ( 2|p_N||p_perp| + M sqrt(g_N) ) a b ,
```

hence with `lam_D = beta - M g_N`, `kappa = 2|p_N||p_perp| + M sqrt(g_N)`,

```
    lam_min(R_T) >= ( lam_A + lam_D - sqrt((lam_A-lam_D)^2 + 4 kappa^2) ) / 2 ,
```

which is `> 0` as soon as `lam_A > kappa^2 / lam_D`. Here `M <= 8.83262267769`
by monotonicity of `H` alone (no covering), and `|p|^2 = sum_n p_n^2 =
int_{-L}^{L} e^{-u} du = 2 sinh L` exactly.

**Sector care.** `s = -1` in the odd sector, where the pole term is `<= 0` and
must be carried into `lam_D` rather than discarded. The table in REPORT §7.2 is
the **even** sector. The odd sector was not assembled.

**What remains.** Only `lam_A`, the least eigenvalue of the `N x N` compression
of `R_T`. At `N = 220` the chain closes for any `lam_A >= 1.0146e-28`.

## 6. The optional penalized comparison (pilot)

If a cutoff `T' < T` is used, the exterior set `E = {|t| >= T' : Psi_L < beta}`
is nonempty but, by monotonicity of `H`, contained in the finite band
`T' <= |t| <= T`. With `b >= sup_E (beta - Psi_L)` and `q` bounding
`int_E|F|^2 / int_R|F|^2`,

```
    (1/2pi) int_{|t|>=T'} (Psi_L - beta)|F|^2  >=  - b q ||f||_2^2 ,
```

giving

```
    Q(f) >= (beta - b q)||f||_2^2 + 2 s F(i/2)^2
            + (1/2pi) int_{|t|<T'} (Psi_L - beta)|F|^2 dt .
```

**The interior subtraction stays `beta`, not `beta - b q`.** Changing both terms
loses the comparison. `b q` is a loss allowance that the finite reduction must
absorb along with every other certified error — it is not a proof that the
exterior integral is nonnegative.
