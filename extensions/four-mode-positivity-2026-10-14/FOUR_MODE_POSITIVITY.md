# FOUR_MODE_POSITIVITY.md — `G ⪰ m₀ I` on the four-mode subspace

The main new theorem of this extension. It is **unconditional** — no Riemann
hypothesis, no external positivity bound — and it is **finite-dimensional**: it
concerns the restriction of the finite-window Weil form to one explicit
four-dimensional space of test functions. It is **not** positivity of the Weil
form, and it is **not** evidence for RH. §5 says why, at length.

---

## 1. The object

`L = 4/5`, `gamma = 14`, both exactly rational. The four `EVEN_D` modes of
[`src/zeta_window.py`](../../src/zeta_window.py),

```
   phi_k(u) = cos(w_k u)/sqrt(L) ,   w_k = (2k+1)pi/(2L) = (2k+1)(5/8)pi ,  k = 0..3,
```

on `[-L,L]` and zero outside; `w_k L = (2k+1)pi/2`, so `phi_k(±L) = 0`, each
`phi_k` is real, even, continuous on `R` and supported in `[-L,L]`, and
`int phi_j phi_k = delta_{jk}`. Write `E = span{phi_0,…,phi_3}` and
`f_x = sum_k x_k phi_k` for `x in R^4`.

`Q_zeta` is the finite-window Weil functional of the paper in `psi`-form, over
the **actual** nontrivial zeros of `zeta`, evaluated through the explicit
formula: pole term, archimedean term, `-log pi` term, prime-power term (support
cutoff `log 5 > 2L`, so exactly `n = 2,3,4`). On `E` it is a quadratic form,

```
   Q_zeta(f_x) = x^T G x ,        G = G^T in R^{4x4} .
```

`G` is **not** available in closed form: its entries involve `pi`, `log pi`,
`log 2`, `log 3` and the full archimedean sum. What the certificate produces is
a rigorous **interval enclosure** of `G` and a positive lower bound for its
smallest eigenvalue.

## 2. Theorem O1

> **Theorem O1 (certified, unconditional).** With `L`, `gamma` and the basis
> exactly as in §1,
> ```
>    G  >=  m0 I      with   m0 >= 1.659599e-10 > 0 ,
> ```
> that is, `Q_zeta(f_x) >= m0 ||x||_2^2` for every `x in R^4`, where `||x||_2 =
> ||f_x||_{L^2}`. Consequently, for the added-quartet form
> `H(delta) = G + 4 a(delta) a(delta)^T` one has `H(delta) >= m0 I` for **every**
> `delta >= 0`, with no `delta`-dependence, because `4 a a^T >= 0`.

Certified by `verify_route_o.py`, assertions `O1`, `O2`, `O2b`; produced by
`scripts/route_o1.py`; data in [`results/route_o_step1.json`](results/route_o_step1.json).

## 3. Proof structure

**(a) The form is exactly quadratic with interval coefficients.** By
Proposition T3 of [`ANALYTIC_TAIL.md`](ANALYTIC_TAIL.md),

```
   Q_zeta(f_x)  =  Qsharp3(x)  ±  R'(x) ,
   Qsharp3(x)   =  Qcore3(x) - (3/4)||x||^2 S_1 + 2 S_3 sum_k w_k^2 x_k^2 ,
```

where `Qcore3` is the exactly-summed part `m < M` (pole, `log pi`, prime,
archimedean partial sum with the sharp `T`) and the two tail terms are diagonal
in `x`. **Every term is a quadratic monomial in `x` with interval coefficients**
— that is exactly what endpoint identity `(E5)` buys: the `O(1/M^3)` correction
`-2C''(0)/s^3` is `+2(sum_k w_k^2 x_k^2)/s^3`, still quadratic, instead of an
unstructured remainder.

So `Qsharp3` polarizes: with `e_j` the coordinate vectors,

```
   Gt_{jj} = Qsharp3(e_j) ,
   Gt_{jl} = [ Qsharp3(e_j + e_l) - Qsharp3(e_j) - Qsharp3(e_l) ] / 2 ,   j != l,
```

an interval matrix `Gt = Gt_mid + Delta` with `|Delta_{jl}| <= rad_{jl}`.
**Checked, not assumed:** `Qsharp3(0) = 0`, `Qsharp3(2x) = 4 Qsharp3(x)`, and
`x^T Gt x` reproduces `Qsharp3(x)` on an off-axis test vector — gates `O2b`. A
constant or linear part would break the second and third of these.

**(b) The two losses, bounded in the right dimension.** For the interval radius,
`|x^T Delta x| <= sum_{j,l} rad_{jl}|x_j||x_l| <= (max_j sum_l rad_{jl})||x||^2`
by `|x_j x_l| <= (x_j^2+x_l^2)/2` and symmetry of `rad` — a row-sum (Gershgorin-
type) bound, correct in any dimension and here evaluated at `n = 4`:

```
   sigma3 := max_j sum_l rad_{jl} = 8.5390e-12 .
```

For the analytic remainder, `|R'(x)| <= 2 sup|C'''|(x) * S_4bound`, and
`sup|C'''| <= sum_k[|A_k| w_k^3 + |B_k|(3 w_k^2 + U w_k^3)]` is itself bounded by
the same row-sum argument applied to the bilinear coefficient matrices of
`A_k(x), B_k(x)`:

```
   sup|C'''|(x) <= Kappa3 ||x||^2 ,  Kappa3 <= 5464.0163 ,
   |R'(x)|      <= rho3  ||x||^2 ,   rho3   <= 1.0111e-13 ,
   sigma3 + rho3 <= 8.6401e-12 .
```

The uniform bound `rho3` is cross-checked against the per-vector value on three
vectors — the frozen vector (`2.0845e-14`), `e_0` (`3.1843e-16`) and `e_3`
(`5.6786e-14`) — each well inside it.

**(c) Hence** `G >= Gt_mid - (sigma3 + rho3) I`, and it suffices to certify
`lambda_min(Gt_mid)`.

**(d) Exact rational LDL.** `Gt_mid` is an exact rational matrix:

```
 [ +0.00020786693004  -0.00069201435592  +0.00155044614207  -0.01304829451352 ]
 [ -0.00069201435592  +0.00231381009004  -0.00524697174952  +0.04600504015038 ]
 [ +0.00155044614207  -0.00524697174952  +0.01229684538678  -0.11962982214071 ]
 [ -0.01304829451352  +0.04600504015038  -0.11962982214071  +1.51682448730770 ]
```

(printed rounded; the exact rationals are in `results/route_o_step1.json`). With
`t_cert = 873/5000000000000 = 1.746000e-10`, Gaussian `LDL^T` over `Q` on
`Gt_mid - t_cert I` terminates with pivots

```
   2.0787e-04 ,  1.0008e-05 ,  4.5290e-06 ,  1.0376e-08 ,
```

all strictly positive, so `Gt_mid - t_cert I > 0` by Sylvester's criterion —
an **exact rational** computation, no rounding. Therefore

```
   m0 = t_cert - (sigma3 + rho3) = 1.746000e-10 - 8.6401e-12 = 1.659599e-10 > 0 . ∎
```

## 4. The spectrum, and why the margin is so thin

Each eigenvalue of `Gt_mid` bracketed by bisection on the **exact inertia**
(`#{eigenvalues < t} = n_neg(Gt_mid - tI)`, Sylvester's law of inertia; no
deflation, no floating point), brackets of width `3.87e-42`:

```
   1.746061370e-10 ,   2.711087885e-06 ,   3.846071003e-03 ,   1.527794227e+00
   condition number  <=  8.750e+09
```

So the four-mode Weil form is **nearly singular**: `lambda_min = 1.75e-10`
against `lambda_max = 1.53`. This single number explains the whole shape of this
extension.

> ### 4.1 A correction to the previous round, recorded in full
>
> An earlier scan of ours computed the detection diagnostic `Theta(delta)` from
> `G_core`, the matrix of the **tail-free** functional. That is not the Weil
> matrix; the true `G` includes the archimedean tail,
> ```
>    G = G_core - (3/4) S_1 I + 2 S_3 diag(w_k^2)  ±  rho3 .
> ```
> The omitted diagonal term is `(3/4)S_1 = 5.7221e-06`, and
> ```
>    lambda_min(G_core) = 5.722182e-06 ,       lambda_min(G) = 1.746061e-10 .
> ```
> `lambda_min(G_core)` is almost exactly `(3/4)S_1`: restoring the tail is a
> near-total cancellation, down to the `1e-10` level. **Every `Theta` value and
> the onset bracket `(17/200, 9/100) = (0.085, 0.09)` reported in that scan are
> wrong and are withdrawn**; the corrected diagnostic is in
> [`NON_DETECTION.md`](NON_DETECTION.md) §3.
>
> **The detection certificate is unaffected.** It never used a matrix: it
> evaluates `Q_zeta(f_x) + 4(a.x)^2 - 4(b.x)^2` directly for a frozen `x` with
> the full machinery, tail included. The error was in the *diagnostic locating
> the onset*, and it was wrong in the safe direction — the true onset is
> **lower** than was reported, so the certified detection interval stands
> further from it, not closer.
>
> It is also the reason the error budget had to be sharpened twice. A `1e-10`
> eigenvalue cannot be resolved by a `4.9e-05` allowance, nor by `1.65e-09`;
> only the `O(1/M^3)` budget (`~2e-14`) reaches it. **There is no floating
> precision in this code, so no amount of precision could have substituted for
> that analytic improvement.**

## 5. What Theorem O1 does and does not mean

This is the section to read before quoting the theorem.

* **It is not Weil positivity, and it is not a proof, weakening or conditional
  version of RH.** Weil positivity is the statement that `Q_zeta(f) >= 0` for
  **all** admissible `f`; that statement is *equivalent* to RH. A
  four-dimensional subspace at one window length and one basis is a measure-zero
  slice of that condition. Positivity on it is **consistent with, and expected
  under, RH** — and a *negative* value would have been the informative outcome,
  since it would have contradicted RH. Finding the expected sign is not
  evidence.
* **It is unconditional, which is the actual content.** The value `m0` is
  derived from the explicit formula and elementary analysis only. No RH, no
  quasi-RH, and in particular **no external positivity input**: the
  repository's earlier local work obtained a positive baseline only by quoting
  Zhu's certified `lambda_min >= lambda*(0.8) >= 8.9e-18` with its own
  normalisation and range. Nothing external is used here.
  [`src/exact_inertia.py`](../../src/exact_inertia.py) certifies inertia for
  **rational** function-field Toeplitz data and cannot be applied to this
  matrix, whose entries involve `log pi`, `log 2`, `log 3`, `pi` and the
  archimedean sum.
* **It says nothing about any zero of `zeta`.** `gamma = 14` enters only through
  the quartet of §1 of [`CERTIFICATE_THEOREM.md`](CERTIFICATE_THEOREM.md); `G`
  itself does not depend on `gamma`. No ordinate is verified anywhere in this
  extension.
* **It is not uniform in `L`, in the basis, or in the dimension.** `L = 4/5` and
  four `EVEN_D` modes throughout. `m0` is a property of this `4 x 4` matrix.
  Nothing here suggests a lower bound surviving `L -> infinity` or `dim -> infinity`;
  the measured near-singularity (`cond ~ 8.7e9` already at `n = 4`) suggests the
  opposite, and the repository's own measurements (`||P_E r_gamma||^2` falling
  from `3.60e-2` at `gamma = 14` to `1.13e-7` at `gamma = 280`) point the same way.
* **It is a computer-assisted certificate with ordinary runtime dependencies**
  (CPython, its `int`/`Fraction` implementation, the OS). It is **not**
  proof-assistant formal verification — no Lean, no Coq, nothing formalised —
  and it has **not** been externally peer-reviewed.
* **No novelty is claimed** for Sylvester's criterion, Sylvester's law of
  inertia, exact rational `LDL^T`, Gershgorin/row-sum bounds, polarization of a
  quadratic form, interval arithmetic, or integration by parts. The new content
  is the *endpoint collapse* of
  [`ANALYTIC_TAIL.md`](ANALYTIC_TAIL.md) §4 and the fact that the resulting
  budget is small enough to resolve `lambda_min(G)` at all.
