# DERIVATION.md — assumptions, normalisation, and every identity used

Every statement carries a tier: **T1** derived analytically here, **T2-N** numerical at a
stated precision, **T2-C** rigorously enclosed (Arb ball) with a reproducible certificate,
**T3** interpretation. Exact statements over `Q` are marked **T0** where they are stronger
than T1 (checked by exact rational arithmetic, no floating point anywhere).

Nothing in this note is imported from elsewhere in the repository. Where a quantity also
appears in a parent extension, agreement is used as a **check**, never as an input; the
parent values used as check targets are listed in §0.3.

---

## 0. Normalisation, fixed before anything is computed

### 0.1 The two coordinate systems, kept separate

| symbol | meaning | never confused with |
|---|---|---|
| `s = beta + i gamma` | the Dirichlet variable | — |
| `delta = beta - 1/2` | horizontal zero displacement | `Delta`, the strip half-width in `z` |
| `D = 2 delta^2` | the squared-displacement diagnostic | — |
| `z = -2i(s - 1/2)` | the spectral coordinate of the heat family | — |
| `tau` | the deformation parameter | **not** physical time |
| `u = x` | the kernel variable | — |
| `L` | the detector window half-length | plays no role in this note |
| `M_H` | the harmonic mean of the partition | `H_tau`, the entire function |

So for `s = beta + i gamma`,

```
Re z = 2 gamma ,       Im z = -2 (beta - 1/2) = -2 delta ,
delta = -Im(z)/2 ,     D = 2 delta^2 = (Im z)^2 / 2 .
```

**The critical line is the real axis in `z`, and `D` is half the squared height above it.**

### 0.2 The heat family (the repository's normalisation, reproduced and re-verified)

```
H_tau(z) = int_0^infty e^{tau u^2} Phi(u) cos(z u) du ,       d_tau H = - d_z^2 H
zeta :  Phi(u)   = sum_{n>=1} (2 pi^2 n^4 e^{9u} - 3 pi n^2 e^{5u}) exp(-pi n^2 e^{4u})
                   H_0(z) = xi(1/2 + i z/2)/8
DH   :  Phi_F(x) = 4 e^{3x} sum_{n>=1} n a_n exp(-pi n^2 e^{4x}/5),  a = (1,kappa,-kappa,-1,0)
                   H_0(z) = Lambda_F(1/2 + i z/2)/2
```

The sign of the PDE, derived rather than quoted: `d_tau H = int u^2 e^{tau u^2} Phi cos(zu) du`
and `d_z^2 H = -int u^2 e^{tau u^2} Phi cos(zu) du`, so `d_tau H = -d_z^2 H` identically —
the **backward** heat equation in `z`. Real-zero preservation runs toward **increasing** `tau`.
This is the Polymath15 convention, in which `zeta` returns the de Bruijn threshold `1/2`.

### 0.3 Dictionary verification — run before any dynamics (T2-N)

`src/check_dictionary.py`, dps 90, oscillation-resolving quadrature. **ALL PASS.**

| check | worst residual | tolerance |
|---|---|---|
| `H_0^zeta(z) = xi(1/2+iz/2)/8`, five `z` | `9.5e-89` rel | `1e-25` |
| `d_tau H = -d_z^2 H`, both kernels, `z` up to 171.4, `tau = 0, 0.2` | `1.3e-40` rel | `1e-18` |
| `kappa` from the Gauss sum vs its closed form | `3.9e-62` | `1e-45` |
| `kappa` vs the parent's committed 40 digits | `1.2e-41` | `1e-38` |
| `H_0^DH(z) = Lambda_F(1/2+iz/2)/2`, three `z` | `2.0e-90` rel | `1e-20` |
| `Lambda_F(s) = Lambda_F(1-s)`, three `s` | `7.1e-58` rel | `1e-45` |
| `Re z_1`, `Im z_1` vs the parent's committed values | `1.2e-28` | `1e-25` |
| `tau_DH = 2(sigma_0-1/2)^2` vs the parent's committed `1.602537483559822777787819` | `4.7e-25` | `1e-22` |

**A precision trap, recorded.** The PDE check at `z = 171.4` first **FAILED** with a
relative residual of `6.0e-10`. That was neither an error in the identity nor in the
kernel: at `|z| ~ 171` the integrand oscillates about 170 times per unit of `u`, and
`mp.quad` on a coarse node list silently returned about ten correct digits while `H_0`
itself is of size `1e-30` (the gamma factor). The fix is two nodes per oscillation
(`kernels._nodes`) and a finite-difference step matched to the working precision; the
residual then falls to `1.3e-40`. **A check that fails for want of resolution must not be
recorded as a failed identity**, and the first version of this table would have done so.

---

## 1. The mean-gap diagnostic (Part 2) — T0, exact over `Q`

`src/means.py`. For `0 <= beta <= 1` the reflected-coordinate partition `p = beta`,
`q = 1 - beta` has `p + q = 1` and

```
A   = (p+q)/2   = 1/2
G   = sqrt(p q) = sqrt(1/4 - delta^2)
M_H = 2pq/(p+q) = 2pq = 1/2 - 2 delta^2
D   = A - M_H   = 2 delta^2 .
```

### 1.1 The structural collapse (T0)

```
(i)   A == 1/2   identically               -- this IS the functional equation, nothing more
(ii)  G^2 == A * M_H                       -- the classical AM.HM = GM^2
(iii) hence M_H = 2 G^2, and D = A - M_H = 2 delta^2 .
```

> **Lemma 1.1.** The mean triple `(A, G, M_H)` is a function of the single scalar `D`, and
> `D = 2 delta^2`. **No mean gap carries information beyond squared zero displacement.**

Verified exactly at eight rational `delta`, including `delta_1 = 0.3085171824566374` (the
DH control's) and the endpoints `0` and `1/2`.

### 1.2 Overlap with the committed partition-coordinate audit

`extensions/rh-partition-coordinate-audit-2026-10-07/PARTITION_DICTIONARY.md` already
contains these means under other names: its **altitude** `h = sqrt(pq)` *is* `G`, its
**response** `R = 1/(4pq)` *is* `1/(2 M_H)`, and its `chi = 2|delta|`; its identity
`4h^2 + chi^2 = 1` is re-checked here. That audit's verdict is **MISMATCH**, and its
finding that *monotonicity is coordinate-free* applies verbatim to `A - M_H` and `A - G`.
**The arithmetic and harmonic means are therefore not new coordinates.** What that
dictionary does **not** contain is a non-monotone member, which is §1.3.

### 1.3 The three gaps as functions of `D`, and the non-monotone one (T0)

On `0 <= D <= 1/2`:

```
A - M_H = D                                 identity; strictly increasing
A - G   = 1/2 - sqrt(1/4 - D/2)             strictly increasing, convex
G - M_H = sqrt(1/4 - D/2) - 1/2 + D         NOT monotone
```

> **Lemma 1.3.** `G - M_H` vanishes at `D = 0` (equality, `p = q = 1/2`, the RH point) and
> again at `D = 1/2` (complete dominance, `beta in {0,1}`, `G = M_H = 0`); it is strictly
> concave in `D`, with
> ```
> maximum  G - M_H = 1/8   attained exactly at  D = 3/8 ,  where  G = 1/4  exactly
> ```
> (`delta = sqrt3/4 = 0.4330127...`). Hence it is **two-to-one** onto `(0, 1/8)`.

Both the critical point and the critical value are **rational** — `sqrt(1/4 - 3/16) = 1/4` —
so this is exact, not a numerical optimum. The two-to-one structure is exhibited by solving
`G - M_H = 1/10` exactly: `D^2 - (7/10)D + 11/100 = 0`, discriminant `1/20 > 0`, roots
`D = 0.2381966011` and `0.4618033989`, i.e. `delta = 0.3451062` and `0.4805223`.

> **Consequence.** `G - M_H` is a **strictly worse** diagnostic than `D`: it cannot
> distinguish a shallow off-line zero from a deep one. `A - G` and `A - M_H` are strictly
> monotone reparameterisations and so carry exactly `D`'s information, no more.

### 1.4 Along any specified evolution `tau |-> D(tau)` (T1)

Differentiating §1.3 in `tau` (`A' = 0` since `A` is constant, `G' = -D'/(4G)`):

```
(A - M_H)' = D'                    -- always the sign of D'
(A - G)'   = D'/(4G)               -- always the sign of D'   (G > 0)
(G - M_H)' = D' (1 - 1/(4G))       -- the sign of D' iff G > 1/4 iff D < 3/8
```

So for `D > 3/8` the gap `G - M_H` moves **opposite** to the displacement diagnostic.

### 1.5 The arc-coordinate identity is not a restoring law (T1 + negative control)

With `delta = (1/2) sin theta`, `M_H = 1/4 + (1/4)cos 2theta`, hence

```
d^2 M_H / d theta^2 = -cos 2theta = 1 - 4 M_H .
```

True, and **empty**: `theta` is a coordinate on a *fixed* configuration, not a time, so the
second derivative is a statement about the chart. **Negative control:** every function
`c_0 + c_1 cos 2theta` satisfies `f'' = -4(f - c_0)`, checked at three `theta` for three
`(c_0, c_1)` including `(7, -3)`. The oscillator form is a property of "quadratic in
`sin theta`", not of `zeta`. No attractor, no restoring force, and no dynamics.

### 1.6 Outside the strip (T0)

For `|delta| > 1/2` one has `pq < 0`, so `G` is imaginary and `M_H < 0`: the nonnegative
binary-partition reading is unavailable, while `D = 2 delta^2` remains a perfectly good
squared-displacement diagnostic. The range matters:

```
zeta, trivial strip :  sup D = 2 (1/2)^2 = 1/2      -- exactly the edge of the domain
DH,  proved strip   :  sup D = 2 (sigma_0 - 1/2)^2 = 1.602537483559822777787819  > 1/2
```

with `sigma_0 = 1.395136158235109721061359` the parent extension's T1 abscissa.
So on the one object in this repository with certified off-line zeros, the partition
interpretation is **not available on the whole admissible range** (T2-C, `src/certify.py` C4).

### 1.7 One exact identity worth keeping (T1)

In the normalisation of §0, de Bruijn's threshold is literally the supremum of the
diagnostic:

```
tau_c(F) := Delta^2 / 2   where Delta = sup_rho |Im z_rho| = 2 sup_rho |delta_rho| ,
sup_rho D(rho) = sup_rho (Im z_rho)^2 / 2 = Delta^2/2 .

=>   tau_c(F) = sup_rho D(rho) = sup_rho (A - M_H)(rho) .
```

So "all zeros of `H_tau` are real for `tau >= Delta^2/2`" reads: *all zeros become real once
`tau` exceeds the initial maximal arithmetic–harmonic mean gap*. For `zeta` this is
`1/2`; for DH it is `1.6025374835598227777878195` (certified, §1.6), reproducing the
parent extension's `tau_DH` to all 22 committed digits. This is the **only** place found
where the mean language is not merely a relabelling — and it is a dictionary identity, not
a mechanism.

---

## 2. The chain rule along the heat flow (Part 3) — T1

Let `tau |-> z(tau)` be a branch of simple zeros, `H_tau(z(tau)) = 0`.

### 2.1 Zero velocity

Differentiating the defining identity in `tau`:

```
0 = d/d tau [ H_tau(z(tau)) ] = (d_tau H)(z) + H'(z) z'(tau) ,
```

and `d_tau H = -H''` by §0.2, so

```
                    z'(tau) = H''(z) / H'(z) .                                      (2.1)
```

Hence, with `y = Im z` and `D = y^2/2`,

```
    y'(tau) = Im( H''/H' ) ,      D'(tau) = y * Im( H''(z)/H'(z) ) .                (2.2)
```

Both agree with the brief's statement of them. **(2.1)–(2.2) are exact chain-rule
identities and nothing more**: they presuppose a simple zero and are undefined where
`H' = 0`.

Equivalently in `delta`: `delta = -y/2`, `D = 2 delta^2`, `D' = 4 delta delta' = y y'` — the
same expression, so the `delta` and `z` bookkeeping agree.

### 2.2 Verified, two independent ways (T2-N, T2-C)

At the DH orbit `O_1` (`src/dh_flow.py`):

* **R1** — the chain rule in `s`: `z' = (i/2) Lambda''(s)/Lambda'(s)`, `s = 1/2 + iz/2`, with
  `Lambda'`, `Lambda''` assembled from **analytic** Hurwitz-zeta and digamma derivatives
  (no finite differences);
* **R2** — the quadrature route: `H, H', H''` as `int e^{tau u^2} Phi_F(u) u^k trig(zu) du`,
  differentiating `cos(zu)` under the integral sign. R2 calls only `exp` and `cos` on the
  kernel and shares **no special-function code** with R1.

```
R1  z'(0) = -0.933346302058786209031062 + 1.88068072103117090749032 i
R2  z'(0) = -0.933346302058786209031062 + 1.88068072103117090749032 i
relative deviation  1.07e-46
precision doubling dps 60 -> 120 on D'(0):  agreement 2.8e-47
```

---

## 3. The canonical product, and the zero-interaction form — T1

This is the step the brief requires be done properly, so it is done in the variable where
it is legitimate.

### 3.1 Symmetry, growth, order

`Phi` is real and **even** (equivalent to the theta functional equation; verified in the
parent audit and re-checked here), so

```
H_tau(-z) = H_tau(z) ,        H_tau(conj z) = conj H_tau(z) ,
```

and the zero multiset is invariant under `z -> -z` and `z -> conj z`. In particular a
non-real zero `z_j = x + iy` carries the **orbit** `{z_j, -z_j, conj z_j, -conj z_j}`,
which is four points when `x != 0` and two when `x = 0`.

Growth: `|H_tau(z)| <= int_0^infty e^{tau u^2}|Phi(u)| cosh(u |Im z|) du`. Both kernels
decay **double-exponentially** (`Phi = O(e^{c u} e^{-c' e^{4u}})`), so for `t = |Im z|` the
exponent `ut - c' e^{4u}` is maximised at `u ~ (1/4) log(t/4c')` with value `~ (t/4) log t`.
Hence, for `tau` in a bounded set,

```
|H_tau(z)| <= exp( C |z| log|z| + C' )        =>   H_tau has order exactly 1 .
```

### 3.2 Genus 0 in the squared variable — why this is the right variable

`H_tau` even means `H_tau(z) = h_tau(z^2)` with `h_tau` entire, and
`max_{|w|=R}|h_tau(w)| = max_{|z|=sqrt R}|H_tau(z)| <= exp(C sqrt R log sqrt R)`, so

```
h_tau has order 1/2 .
```

Order `< 1` forces **genus 0**. Therefore, writing `w_k = a_k^2` for the zeros of `h_tau`
(one `a_k` from each `+-` pair) and assuming `h_tau(0) = H_tau(0) != 0`,

```
h_tau(w) = h_tau(0) prod_k ( 1 - w/w_k ) ,                                           (3.1)
```

**absolutely convergent, with no exponential factor and no normalisation left free**
(`sum_k 1/|w_k| < infinity`, since order `1/2` gives convergence exponent `<= 1/2`).
If `H_tau(0) = 0` insert `w^m`; nothing below changes. For `zeta` at `tau = 0`,
`H_0(0) = xi(1/2)/8 != 0`.

> **This is the step that would be illegitimate in `z`.** There `sum_k 1/|z_k|` diverges
> (the zero density is `~ log`), the genus is 1, and Hadamard would force an unknown factor
> `e^{cz}`. Working in `w = z^2` removes it. Reporting `H'/H = sum 1/(z - z_k)` without this
> is exactly the normalisation error the brief warns against.

### 3.3 The interaction form (T1)

From (3.1),

```
H'/H (z) = 2z sum_k 1/(z^2 - a_k^2) ,                                                (3.2)
```

locally uniformly convergent off the zeros. Near a simple zero `a_j` the `k = j` term is
exactly `2z/(z^2-a_j^2) = 1/(z-a_j) + 1/(z+a_j) = 1/w + 1/(2a_j) + O(w)` with `w = z - a_j`,
while Taylor gives `H'/H = 1/w + H''(a_j)/(2H'(a_j)) + O(w)`. Comparing,

```
H''(a_j)/H'(a_j) = 2 [ 1/(2 a_j) + sum_{k != j} 2 a_j/(a_j^2 - a_k^2) ]
                 = 2 sum_{k != j} 1/(a_j - z_k) ,                                    (3.3)
```

the second line being the same thing rewritten over the **full symmetric zero list**
`{+-a_k}` in the **paired order** — `1/(2a_j) = 1/(a_j - (-a_j))` and
`2a_j/(a_j^2-a_k^2) = 1/(a_j-a_k) + 1/(a_j+a_k)`. The pairing is **not optional**: the
unpaired sum `sum 1/(a_j - z_k)` is not absolutely convergent.

### 3.4 The sign identity, and why its tail is harmless (T1)

Taking imaginary parts in (3.3), with `Im(1/(z_j - z_k)) = -(y_j - y_k)/|z_j - z_k|^2`:

> **Identity 3.4 (the sign identity).** At a simple zero `z_j` of `H_tau`,
> ```
>     D'(tau) = -2 y_j  sum_{k != j}  (y_j - y_k) / |z_j - z_k|^2 ,
> ```
> the sum over the full zero multiset. If the zeros lie in a bounded strip
> `|Im z_k| <= Delta` — one of de Bruijn's own hypotheses, and true for `zeta` and for DH —
> then `|y_j - y_k| <= 2 Delta` and `sum_k |z_j - z_k|^{-2} < infinity` (order 1), so the
> sum is **absolutely convergent**.

So, unlike `z'` itself, `D'` has an absolutely convergent zero-interaction expansion: the
bounded strip does the work. The partial sums are nevertheless **not** used for any reported
value — every `D'` in this investigation comes from `H''/H'`, where the whole spectrum is
present implicitly. §4 reports the remainder as a measured residual, which is a check on
(3.4), not an ingredient of anything.

### 3.5 Corollary (S1) — the own-orbit term, and when the sign is forced (T1)

Take `y_j > 0` (if `y_j < 0` replace `z_j` by `conj z_j`; `D` is unchanged). Term by term:

| partner | `y_k` | contribution to `D'` |
|---|---|---|
| `conj z_j` | `-y_j` | **exactly `-1`** |
| `-conj z_j` | `+y_j` | **exactly `0`** |
| `-z_j` | `-y_j` | `-y_j^2/\|z_j\|^2` |
| any zero with `y_k <= y_j` | | `<= 0` |
| a zero with `y_k > y_j` | | `> 0` |

> **(S1a)** The own orbit contributes exactly `-1 - y_j^2/|z_j|^2` when `Re z_j != 0`, and
> exactly `-1` when `Re z_j = 0`. In both cases `<= -1`.
>
> **(S1b)** A zero contributes **positively** to `D'` only if it lies **strictly further
> from the critical line on the same side**.
>
> **(S1c)** Hence if `|y_j| = max_k |y_k|` — `z_j` is at maximal distance from the line —
> then every term is `<= 0` and
> ```
>                           D'(tau) <= -1 .
> ```
> Quantitatively, if `sup_k |y_k| = Delta_tau` and `|y_j| >= Delta_tau - eps`, then
> `D' <= -1 + 2 |y_j| eps S_j` with `S_j = sum_{k != j} |z_j - z_k|^{-2} < infinity`.

**(S1c) is the differential form of de Bruijn's strip bound.** Writing
`Dscr(tau) = sup_j D_j(tau) = Delta_tau^2/2`, (S1c) says the upper derivative of `Dscr` is
`<= -1` wherever `Dscr > 0`, which integrates to `Delta_tau^2 <= max(Delta_0^2 - 2 tau, 0)`
— exactly the threshold `tau_c = Delta^2/2` of §1.7. **I do not claim a complete new proof
of de Bruijn's theorem by this route**: that would need the envelope argument carried
through collisions and the attainment of the supremum, neither of which is done here. What
is claimed is (3.4), (S1a–c), and the observation that integrating (S1c) reproduces exactly
the classical threshold and nothing better.

### 3.5b Corollary (S1d) — a necessary ordinate-separation condition (T1)

This is the sharpest consequence of the sign identity beyond (S1c), and the one new
statement in this investigation that is cheap to apply. Take `y_j > 0` and a competing zero
with `y_k > y_j > 0`. Its contribution to `D_j'` is

```
2 y_j (y_k - y_j)/|z_j - z_k|^2  <=  2 y_j (y_k - y_j) / R^2 ,     R := |Re z_j - Re z_k| ,
```

and `max_{0<a<b} 2a(b-a) = b^2/2` at `a = b/2`, so with `y_k <= Delta_tau`,

```
(one competing zero's contribution to D_j')  <=  Delta_tau^2 / (2 R^2) .
```

For a **single** competitor to overturn the conjugate-partner term `-1` it must exceed `1`:

> **(S1d)** A single competing zero can make `D_j'` exceed `-1` only if
> ```
>                 R = |Re z_j - Re z_k|  <  Delta_tau / sqrt 3 ,
> ```
> and `N` competitors all at ordinate distance `>= R` only if
> `R < Delta_tau sqrt(N/2)`. In the Dirichlet variable, `R = 2|gamma_j - gamma_k|` and
> `Delta_tau = 2 delta_max`, so this reads `|gamma_j - gamma_k| < delta_max/sqrt 3`:
> **a competitor must share the ordinate to within a fraction of the strip half-width.**
>
> **(S1d')** And without optimising over `y_j` — for fixed `y_j` the unconstrained maximum
> of `2 y_j h/(R^2+h^2)` is at `h = R`, with value `y_j/R` — one gets the sharper per-zero
> form
> ```
>                 R < y_j ,       i.e.      | gamma_j - gamma_k | < | delta_j | .
> ```
> If the strip cuts the optimum off (`h = Delta - y_j < R`) then
> `2 y_j h > R^2 + h^2 >= 2Rh` by AM–GM, so `y_j > R` again: **(S1d') is necessary in all
> cases.** The audited condition is the conjunction `R < min(y_j, Delta_tau/sqrt3)`.
>
> **AMENDED by `AUDIT.md` §1.1.** The first version of (S1d) stated
> `R < Delta_tau/sqrt2`, which is **valid but 22.5 % too loose**
> (`phi(1/sqrt2) = 0.732050807569 < 1`), and its `N`-competitor form
> `R < Delta_tau sqrt(N)/2`, which is **wrong** — it contradicts its own `N = 1` case.
> The sharp constant is `1/sqrt3`, attained exactly at `(y_j, y_k, h) =
> (2Delta/3, Delta, Delta/3)` where the contribution is exactly `1`; the stationary-point
> algebra reduces to `3u^2 - 4u + 1 = 0`. (S1d') was available and was missed. No earlier
> conclusion changes: the old constant is valid, so everything it excluded stays excluded.

Necessary, not sufficient. Three readings of it (`src/separation_bound.py`):

1. **The minimal counterexample is confluent in the ordinate.** In §4.1's
   `(z^2+1)(z^2+4)` the competitor sits at `R = 0` — the two orbits share their ordinate
   exactly. That is the configuration the repository's
   `extensions/moderate-gap-spectral-bridge-2026-10-21/06_CONFLUENT_BLIND_SPOTS.md`
   identifies as structural (an off-line zero's functional-equation partner is always at
   ordinate distance exactly `0`); here the partner supplies the `-1`, and a **second,
   deeper** orbit at the same ordinate is what overturns it.
2. **No DH census quartet can be a counterexample.** With `Delta_DH = 1.790272316470219`
   the sharp threshold is `R < Delta/sqrt3 = 1.0336142`, and the per-zero form (S1d') is
   stronger still: the deepest competitor for each of `O_2, O_3, O_4` is `O_1`, outside
   `|gamma_j - gamma_k| < |delta_j|` by factors **189, 1086, 406** respectively
   (`AUDIT.md` §1.1). So the four certified `D'(0) < 0` are now *explained*, not merely
   observed.
3. **For `zeta` the condition is not vacuous, and that is as far as it goes.** The trivial
   strip gives `Delta <= 1`, so the (S1d) threshold is `|gamma_j - gamma_k| < 1/(2 sqrt 2)
   = 0.3535533906` as first stated and `1/(2 sqrt 3) = 0.2886751346` once sharpened. The
   mean spacing `2 pi/log(gamma/2 pi)` reaches the looser value at
   `gamma* = 2 pi exp(4 pi sqrt 2) = 328290294.5`, and Lehmer-type pairs are inside it far
   earlier (`6709/6710` at `gamma = 7005.06` has `|d gamma| = 0.0377`). The per-zero form
   (S1d') is the binding one and is **far** more restrictive for a shallow zero: it demands
   `|gamma_j - gamma_k| < |delta_j|`, so a zero at displacement `10^{-3}` needs a deeper
   zero within `10^{-3}` in ordinate. So (S1d)–(S1d') **do not rule out** non-monotonicity
   for `zeta`; they localise where it could occur. **No claim is made that it does occur** —
   that would be a statement about off-line zeros of `zeta`.

(S1d) uses only the sign identity and the strip half-width, so like everything else here it
does **not** distinguish `zeta` from Davenport–Heilbronn, and it gives no information at
`tau = 0` beyond what the strip width already gives.

### 3.6 Corollary (S2) — a monotone functional, with its hypothesis stated (T1)

Symmetrising (3.4) over unordered pairs — the `(j,k)` and `(k,j)` terms sum to
`-2(y_j-y_k)^2/|z_j-z_k|^2` — gives

> **(S2)** Let `Escr(tau) := sum_j D_j(tau) = (1/2) sum_j (Im z_j)^2`, the **total** mean
> gap over the zero multiset. If the number `m(tau)` of non-real zeros is **finite** — see
> the amendment below; `Escr < infinity` is NOT enough — then
> ```
>   Escr'(tau) = -2 sum_{j<k} (y_j - y_k)^2 / |z_j - z_k|^2   <=  0 ,
> ```
> absolutely convergent, with equality iff all `Im z_j` are equal (which, with the `+-` and
> conjugate symmetry, forces all `Im z_j = 0`).

This **is** a monotone functional, it **is** explicitly defined, its monotonicity **is**
proved, and it **is** rigorously the sum of the zero-displacement diagnostic — the four
conditions the brief imposes before such a thing may be mentioned. It is nevertheless not
progress; §3.7 says why.

Rearrangement is legitimate because only pairs with at least one non-real member contribute
(both real gives numerator `0`), and for each such pair the numerator is bounded while
`sum_k |z_j-z_k|^{-2} < infinity`.

> **AMENDED by `AUDIT.md` §1.3 — this was a real error.** The first version required only
> `Escr(tau) < infinity`. That is **not sufficient, and in the infinite case it is
> incompatible**: every conjugate pair contributes **exactly 1** to the dissipation, so
> `Diss >= m/2` and `m = infinity` forces `Diss = +infinity`, while `Escr` can be finite
> (take `D_j = 2^{-j}`) — a finite left side and an infinite right side. The correct
> hypothesis is `m(tau) < infinity`, **finitely many non-real zeros**, under which
> everything above is automatic. The amendment *strengthens* the filter of `REPORT.md` §6(d):
> Davenport–Heilbronn fails `m < infinity` on the strength of the original Davenport–
> Heilbronn 1936 infinitude theorem alone, with no summability guess and no second source.
> (Bibliographic, not verified here.)

### 3.7 What (S2) buys, exactly — and why it is less than de Bruijn (T1)

Each conjugate pair `{z_j, conj z_j}` with `y_j != 0` contributes **exactly `1`** to the
dissipation sum (`(2y_j)^2/(2|y_j|)^2 = 1`). So with `m(tau)` = the number of non-real
zeros, and *provided* `Escr` is absolutely continuous in `tau` on `[0, Lambda]` — which needs
the zero multiset to vary continuously (Hurwitz, fine) **and** `Escr` to stay finite there,
so that the finitely many non-real zeros can be followed as a group through collisions
(§3.8); this is a hypothesis, not a consequence, and it is the step where the infinite case
differs from the polynomial one —

```
Escr'(tau) <= -m(tau) <= -2   whenever Escr(tau) > 0 ,
=>   Escr(tau) = 0  for  tau >= Escr(0)/2 ,
=>   Lambda  <=  Escr(0)/2  =  (1/2) sum_rho D(rho) .                                (3.5)
```

But every non-real zero comes with a conjugate partner carrying the same `D`, so
`sum_j D_j >= 2 max_j D_j`, hence

```
Escr(0)/2  >=  max_rho D(rho)  =  Delta_0^2/2 ,
```

with equality exactly when the non-real zeros are a single purely imaginary pair.

> **de Bruijn's bound is always at least as strong as (3.5), and strictly stronger as soon
> as there are four or more non-real zeros.** Confirmed on the toy model: `Delta_0^2/2 = 2`,
> `Escr(0)/2 = 5/2`, truth `0.699`.

And the obstruction is now sharp. The **only unconditional** lower bound on the dissipation
is the conjugate-partner term, `1` per pair. Anything larger needs a lower bound on
`sum_{j<k}(y_j-y_k)^2/|z_j-z_k|^2` over the remaining pairs, i.e. a guarantee that off-line
zeros are **not clustered at equal imaginary part** near one another. That is a statement
about the very configuration to be excluded. §4.2 measures how badly a naive density
heuristic fails at supplying it.

### 3.8 Collisions — where branch tracking stops (T1)

At a zero of multiplicity `>= 2`, `H' = 0` and (2.1) is undefined; **it is not used there**.
Three statements, all exact:

1. A real simple branch stays real. `H_tau` is real on `R`, so at a real simple zero
   `H''/H'` is real and `z'` is real. A branch can leave the real axis **only through a
   collision**. (Verified numerically in §5.)
2. The locus of collisions is cut out by `H_tau(z) = H'_tau(z) = 0`, two equations in
   `(z, tau)`; the Newton system uses `d_tau H = -H''` and `d_tau H' = -H'''`, so it needs
   `H'''` and no finite differences.
3. **Local factorisation.** For a cluster of `nu` branches colliding at `(z_*, tau_*)`, the
   elementary symmetric functions of the cluster are analytic in `tau` through `tau_*`
   (they are, up to units, coefficients of the Weierstrass factor of `H_tau` on a small
   disc, given by contour integrals of `H'/H` on a fixed circle that no zero crosses). So
   `sum_{cluster} z_j` and `sum_{cluster} z_j^2` extend analytically, while the individual
   `D_j` need not. In the toy model this is explicit: the two colliding branches are the
   roots of the exact quadratic `w^2 - (12 tau - 5) w + (4 - 10 tau + 12 tau^2) = 0` in
   `w = z^2`, whose coefficients are polynomial in `tau`. `Escr` is continuous through the
   collision; `Escr'` **jumps** (from `-12` to `-8`), because two heights have merged and
   their pair term has fallen to zero.

---

## 4. The controls

### 4.1 Control 1 — the finite polynomial heat flow (T0, exact over `Q`)

`src/toy.py`. If `P_0` is a polynomial then `P_tau = exp(-tau d_z^2) P_0 = sum_m
(-tau)^m/m! P_0^{(2m)}` is a **finite** sum, so `P_tau` is polynomial in `(z, tau)` and
`d_tau P = -P''` is an exact identity in `Q[z,tau]` (checked, residual exactly `0`, for
three `P_0`). Every statement of §§2–3 then holds with a finite Hadamard product, no
convergence question and no tail.

> **Theorem 4.1 (counterexample; exact over `Q`).** Let `P_0(z) = (z^2+1)(z^2+4)`, zeros
> `+- i` and `+- 2i`. Under the flow,
> ```
>     at z = 2i  (the extremal zero):   D'(0) = -19/3   <= -1      -- as (S1c) requires
>     at z = i   (the inner zero):      D'(0) = +1/3    >  0       -- D INCREASES
> ```
> Hence **universal monotonicity of the individual diagnostic `D` is refuted**, not merely
> unproved.

Three independent derivations agree exactly: the interaction identity (3.3)/(3.4) over
`Q(i)`; the closed form `P_tau = z^4 + (5-12tau)z^2 + (4-10tau+12tau^2)` with
`D_in(tau) = (5 - 12tau - sqrt(96tau^2-80tau+9))/4`, whose derivative at `0` is
`(-12 + 80/6)/4 = 1/3` exactly (the surd is `sqrt 9 = 3`, so this is rational); and
root-finding at 14 nodes agreeing with the closed form to `1e-30`. `D_in` **increases
across the whole branch**, from `0.5` to `0.827` — not a single-point artefact.

**The exact threshold.** For `P_0 = (z^2+a^2)(z^2+b^2)` with `0 < b < a`, the inner pair has
`D'(0) = 4b^2/(a^2-b^2) - 1`, so

```
     D'(0) > 0   <=>   a^2 < 5 b^2 ,       with  D'(0) = 0  exactly at  a^2 = 5 b^2 .
```

Checked at five exact rational `(a^2,b^2)` including the boundary `(5,1)`.

**Minimality.** A single orbit can never do this: for a lone quartet `{+-z, +-conj z}`,
(S1a) is the whole sum and `D' = -1 - y^2/|z|^2 < 0`. Two orbits — degree 4 — is minimal.

**The collision, and `Escr`.** The discriminant `96tau^2 - 80tau + 9` vanishes at
`tau_* = (10 - sqrt46)/24 = 0.13406958` (the two imaginary pairs collide; branch tracking
stops) and at `tau_** = (10 + sqrt46)/24 = 0.69926375` (the quartet lands on the real axis,
`Escr = 0`). `Escr(tau) = 5 - 12tau` exactly before the collision and
`sqrt(12tau^2-10tau+4) - (12tau-5)/2` after; `Escr = 0` iff `96tau^2-80tau+9 = 0` with
`12tau > 5`, i.e. exactly at `tau_**`. Compare `Delta_0^2/2 = 2` (de Bruijn) and
`Escr(0)/2 = 5/2` (3.5): neither is sharp, and (3.5) is the weaker, as §3.7 proves.

**Negative control on "entropy".** The other obvious quadratic functional,
`sum_j Re(z_j^2) = sum_j (x_j^2 - y_j^2)`, **increases at exactly `2n(n-1)`** for every
polynomial heat flow of degree `n` — because `d_tau P` has no `z^{n-1}` term, so `e_1` is
conserved and `e_2' = -n(n-1)`. It is monotone, universally, and it is **not** a
displacement diagnostic: it increases even when every zero is real and `D == 0` throughout.
Monotone is cheap; *monotone and tied to `D`* is what (S2) supplies, and (S2) is what §3.7
shows to be weak.

### 4.2 Control 2 — Davenport–Heilbronn, which has certified off-line zeros at `tau = 0`

`src/dh_flow.py`, `src/dh_interaction.py`, `src/certify.py`. `D'(0)` is computed from
`H''/H'` — so the **infinite zero sum is never truncated**.

| orbit | `delta` | `D` | certified `D'(0)` | own-orbit term | cross-orbit rest |
|---|---|---|---|---|---|
| `O_1` | `0.3085171825` | `0.1903657037` | `-1.16044463430610822178313161420` | `-1.0000129598` | `-0.1604316745` |
| `O_2` | `0.1508300806` | `0.0454994264` | `-1.04435944844102968030569764300` | `-1.0000017455` | `-0.0443577029` |
| `O_3` | `0.0743560504` | `0.0110576445` | `-1.01207352431659836487600009985` | `-1.0000001995` | `-0.0120733248` |
| `O_4` | `0.2242576946` | `0.1005830272` | `-1.12778399346959232977056785846` | `-1.0000016107` | `-0.1277823828` |

All four are **T2-C**: Arb balls at 400 bits with radii `<= 4.5e-32`, strictly negative,
with one-sided upper bounds printed as plain decimals rounded up. Each is certified *as an
analytic function of the supplied `z_1`* — `H'(z_1) != 0` is part of the certificate — and
that `z_1` encloses an actual zero rests on the parent census, which is **T2**. This is the
parent extension's own convention and is not upgraded here.

**The interaction identity, checked with its remainder.** `src/dh_interaction.py` computes
169 on-line ordinates on `0 < Im s < 260` independently (sign changes of the real
Hardy-type `Lambda_F(1/2+it)`, bisected to `1e-35`) and sums (3.4) over them and the four
census quartets:

| cutoff `\|Im s\| <` | zeros used | partial cross-orbit sum | remainder vs exact | rel. to `\|D'\|` |
|---|---|---|---|---|
| 25 | 16 | `-0.00046489` | `-0.15997` | 13.8 % |
| 100 | 104 | `-0.14965664` | `-0.010775` | 0.93 % |
| 200 | 256 | `-0.15823179` | `-0.0021999` | 0.19 % |
| 260 | 350 | `-0.15877279` | `-0.0016589` | 0.14 % |

The remainder shrinks monotonically, as an absolutely convergent tail must. It is reported
as a **measured residual, not a bound**, and no conclusion rests on it.

*An independent reproduction, worth recording.* The scan finds **122** on-line ordinates
below `T = 200`, and `122 + 2 x 4 = 130` — exactly the parent extension's contour count and
its reconciliation, obtained here from a separate implementation.

**Why the on-line drag is weak — and a heuristic that must not be used.** A mean-spacing
estimate of `sum_k |z_1 - z_k|^{-2}` over the on-line zeros — computed, not guessed: for a
uniform lattice of spacing `h` in `z`, averaged over the offset, the sum is exactly
`pi/(h|y|)` — predicts a cross-orbit contribution of `-1.3026848`; the truth is
`-0.16043167`, an **overestimate by a factor 8.12**. The reason is structural: the nearest
on-line ordinates to `gamma_1 = 85.69934849` are `83.1088410847` and `87.6474763326`, a gap
of `4.5386352` against a mean spacing of `1.488058` — **3.05 times the mean spacing**. The
off-line quartet sits in a large zero gap. Note the size of the error: the predicted
`-1.3026848` is **larger in magnitude than the exact conjugate-partner term `-1`**, so a valid
version of it would roughly double the decay rate. So the density heuristic is **not** usable
as a lower bound on the dissipation of §3.7, which is exactly the estimate that route needs.

### 4.3 Control 3 — `zeta` itself, with no RH assumption and no planted family

`src/zeta_flow.py`. See `REPORT.md` §5 for the measured output. The structure:

* **(Z1)** On a real branch, `Im(z') = 0` to the working precision at every step: the branch
  is pinned to the axis, `D == 0`, `D' == 0`. Tracking a real `zeta` zero forward in `tau`
  can **never** produce an off-line zero, so "the zero stayed real" is evidence of nothing.
* **(Z2)** A branch leaves the axis only through a collision, located by solving
  `H_tau = H'_tau = 0` for `(z,tau)` by Newton (so `H'` is never divided by at the
  collision itself).
* **(Z3)** The emergent conjugate pair is continued below the collision; the local quadratic
  model gives `Im z ~ sqrt(2(tau_c - tau))`, hence `D ~ (tau_c - tau)` and `D' ~ -1`.
* **(Z4)** But that pair is **not** extremal there, so (S1c) does not apply to it. The exact
  two-zero model `P = z^2 - (g/2)^2` collides at `tau = -g^2/8`, so a pair of gap `g` in `z`
  leaves the axis at that scale; the Lehmer pair (zeros 6709/6710, `gap in z = 0.0754`) has
  scale `-7.1e-4` against `-1.04` for the closest pair among the first twenty zeros — a ratio
  of `1459.3753`. Zeros
  therefore leave the axis in order of **increasing** gap, and the low-height pair is one of
  the last. `src/zeta_lehmer.py`.

> So in **neither** control is the sign-forcing hypothesis verifiable: on the DH side `O_1`
> is extremal only within a finite census and `O_2..O_4` provably are not; here the tracked
> pair provably is not. Extremality quantifies over **all** zeros, so no finite computation
> can establish it — which is exactly why the forced-sign statement integrates only to the
> global strip bound.

---

## 5. What is assumed, and where

| assumption | used in | status |
|---|---|---|
| `Phi` real and even (equiv. the functional equation) | §3.1 onward | T1 (and verified) |
| double-exponential decay of `Phi` | §3.1 growth, §3.2 order | T1 |
| `H_tau(0) != 0` | (3.1); otherwise insert `w^m` | T2 at `tau = 0`; harmless |
| zeros in a bounded strip `\|Im z\| <= Delta` | absolute convergence in (3.4) | T1 for `zeta`, T1 for DH (`sigma_0`) |
| the tracked zero is **simple** | (2.1), (2.2), (3.3) | hypothesis; collisions handled in §3.8 |
| `sum_rho D(rho) < infinity` | (S2), (3.5) | **unproved for `zeta`; not known for DH** |
| the DH census zero coordinates | the *location* of `O_1..O_4` only | **T2**, parent extension |
| `sigma_0 = 1.3951361582...` | §1.6, §1.7, (S1c) for DH | T1, parent extension |

Nothing in this note assumes RH, GRH, Weil positivity, simplicity of `zeta`'s zeros, or any
bound on `Lambda` beyond de Bruijn's own `Delta^2/2`.
