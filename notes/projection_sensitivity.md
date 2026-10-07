# Projected-derivative sensitivity in finite Weil windows — proof notes

Companion to `paper/weil_window_control.tex` §"Projected-derivative sensitivity and
near-null overlap". Prepared 2026-10-07 for the projection revision.

Evidence tiers used throughout: **T1** = proved here or in a cited source; **T2** =
numerically corroborated in high-precision arithmetic (not proved); **T3** =
interpretation. Nothing below is offered as a statement about the Riemann Hypothesis.

---

## A. Conventions

Let `C/F_q` be an honest hyperelliptic curve of genus `g` (smooth, geometrically
irreducible, with a degree-2 map to `P^1`). Write

- `p_n = q^n + 1 - #C(F_{q^n})`, the trace of `Frob^n`;
- `alpha_1, alpha_1-bar, ..., alpha_g, alpha_g-bar` for the `2g` inverse Frobenius
  roots, so `|alpha_j| = q^{1/2}` by the Weil bounds (the RH for curves);
- `alpha-hat = alpha q^{-1/2}`, the **normalised** roots, which lie on the unit circle;
- `alpha-hat_j = e^{i theta_j}` with `theta_j in (0, pi)`, assumed **distinct** and
  none equal to `0` or `pi` (the generic case; `pm 1` would be a real root).

The normalised point-count sequence is

```
t(n) := p_n q^{-n/2} = sum over all 2g normalised roots of (alpha-hat)^n
      = sum_{j=1..g} 2 cos(n theta_j),          t(0) = 2g.                      (A.1)
```

`t` is real and even in `n`. Fix a **window degree** `R >= 0`. The coefficient space is
`R^{R+1}` with `c = (c_0, ..., c_R)`, and

```
P_c(z) := sum_{i=0}^{R} c_i z^i,        T_R[i,k] := t(|i-k|),  0 <= i,k <= R.    (A.2)
```

`T_R` is the real symmetric Toeplitz matrix of the point-count data. This is the
**zero-side** object: by (A.1) its entries are built from the Frobenius angles, i.e.
from the zeros of the zeta function of `C`, via the point counts.

---

## B. The Gram identity (T1)

**Proposition B.1.** For every `c in R^{R+1}`,

```
c^T T_R c = sum_{j=1}^{g} 2 |P_c(e^{i theta_j})|^2.                              (B.1)
```

For every `c in C^{R+1}`,

```
c* T_R c = sum_{j=1}^{g} ( |P_c(e^{i theta_j})|^2 + |P_c(e^{-i theta_j})|^2 ).   (B.2)
```

*Proof.* Since `t` is even, `t(|i-k|) = t(i-k)`. For `c in C^{R+1}`, using (A.1) and
`2 cos(m theta) = e^{i m theta} + e^{-i m theta}`,

```
c* T_R c = sum_{i,k} conj(c_i) c_k t(i-k)
         = sum_j sum_{i,k} conj(c_i) c_k ( e^{i(i-k)theta_j} + e^{-i(i-k)theta_j} ).
```

For either sign `s = pm 1`,

```
sum_{i,k} conj(c_i) c_k e^{s i (i-k) theta}
   = ( sum_i conj(c_i) e^{s i i theta} ) ( sum_k c_k e^{-s i k theta} )
   = conj( P_c(e^{-s i theta}) ) * P_c(e^{-s i theta})
   = | P_c(e^{-s i theta}) |^2,
```

which gives (B.2). If `c` is real then `P_c` has real coefficients, so
`P_c(e^{-i theta}) = conj(P_c(e^{i theta}))` and the two moduli coincide, giving the
factor `2` of (B.1). ∎

**Coefficient and indexing.** The `2` in (B.1) is *not* a normalisation choice: it is the
count of normalised roots in each conjugate pair. The unindexed form of the identity is
`c* T_R c = sum over all 2g normalised roots alpha-hat of |P_c(conj(alpha-hat))|^2`;
grouping the `2g` roots into `g` conjugate pairs and using reality of `c` collapses each
pair's two equal terms into one term with coefficient `2`. For complex `c` the two terms
of a pair are **unequal** and must both be retained, which is why (B.2) has no factor `2`.

**Remark B.2 (what this is, and is not).** (B.1) is the classical zero-side Gram identity:
the quadratic form built from point-count data is a nonnegative combination of squared
moduli evaluated at the Frobenius angles. It is the finite function-field analogue of Weil
positivity, and it is where positivity comes from for free: `T_R ⪰ 0` **always**. It is
**not** an input to the Riemann Hypothesis for number fields, and nothing below transfers
positivity to the zeta side.

**Corollary B.3 (the exact kernel).** Put `E := ker T_R`. Then

```
E = { c in R^{R+1} : P_c(e^{i theta_j}) = 0 for j = 1..g },
```

and for `R >= 2g`, `dim E = R + 1 - 2g`.

*Proof.* `T_R ⪰ 0` by (B.1), so `c^T T_R c = 0` forces every summand to vanish. For real
`c`, `P_c(e^{i theta_j}) = 0` is equivalent to divisibility of `P_c` by the real quadratic
`z^2 - 2 cos(theta_j) z + 1`. The `g` quadratics are coprime (the `theta_j` are distinct
and in `(0,pi)`), so `E = { m * Pi : deg m <= R - 2g }` with
`Pi(z) = prod_j (z^2 - 2cos(theta_j) z + 1)` of degree `2g`. Hence `dim E = R+1-2g` when
`R >= 2g`, and `E = {0}` when `R < 2g`. ∎

`E` is an **exact** kernel: it is cut out by `2g` linear conditions with algebraic
coefficients, not by a numerical threshold. This is the one place in this work where a
null space is exact; the zeta-side object of §E is a **near-null band defined by a
cutoff** and is never called a kernel.

---

## C. The deformation family, and what realises it

Fix one angle, `theta := theta_1`, and for `a in R` define the deformed sequence

```
t_a(n) := t(n) + 2 cos(n theta) ( cosh(n a) - 1 ),    T_a := Toeplitz(t_a) on R^{R+1}.  (C.1)
```

`t_a` is real, even in `n`, and even in `a`; `t_0 = t`.

**Proposition C.1 (realisability; T1).** `t_a` is *not* the normalised point-count
sequence of any curve for `a != 0`. Precisely:

1. The two-point reciprocal set `{ e^{a+i theta}, e^{-a-i theta} }` has power sums
   `2 cosh(na + i n theta) = 2 cos(n theta) cosh(na) + 2i sin(n theta) sinh(na)`,
   which is **non-real** for `a != 0`; so no reciprocal pair realises (C.1).
2. The conjugation- and reciprocal-closed **quadruple** `{ e^{pm a pm i theta} }` has
   power sums `4 cos(n theta) cosh(na)`, real, whose increment over the `a=0` value
   `4 cos(n theta)` — a *doubled* angle — is `4 cos(n theta)(cosh(na) - 1)`, i.e.
   **exactly twice** the increment in (C.1).

So (C.1) is the *symmetrised* (real-part) deformation: a legitimate real symmetric
Toeplitz deformation of the point-count data, but **a deformed spectrum need not come
from a curve**, and this one does not. Consequently every statement below is a statement
about the Toeplitz family (C.1), not about a family of curves. The genuine quadruple
family is covered by the same theorem with the coefficient doubled (§D, Remark D.5);
both were checked numerically (`src/stage1_rank.py`).

**Evenness and the `a^2` expansion (T1).** From `cosh(x) = sum_{m>=0} x^{2m}/(2m)!`,

```
t_a(n) - t(n) = 2 cos(n theta) sum_{m>=1} (na)^{2m}/(2m)!
              = a^2 n^2 cos(n theta) + a^4 n^4 cos(n theta)/12 + O(a^6),
```

convergent for every `a`. Hence `a -> T_a` is an **entire, even** family of real symmetric
matrices:

```
T_a = T_0 + a^2 B + a^4 C + O(a^6),
   B[i,k] = (i-k)^2 cos((i-k) theta),      C[i,k] = (i-k)^4 cos((i-k) theta)/12.   (C.2)
```

Only *after* establishing this evenness is one entitled to write
`T_eps = T_0 + eps^2 B + O(eps^4)`: the absence of an `O(a^3)` term is a consequence of
`cosh` being even, not an assumption.

---

## D. The compressed operator and the projection theorem

For `m = 0,1,2` and `c in R^{R+1}` set

```
S_m(c) := sum_{i=0}^{R} i^m c_i e^{i i theta}.                                    (D.1)
```

Note `S_0(c) = P_c(e^{i theta})`.

**Proposition D.1 (the quadratic form `B`; T1).**

```
c^T B c = 2 Re( S_2(c) conj(S_0(c)) ) - 2 |S_1(c)|^2.                             (D.2)
```

*Proof.* Expand `(i-k)^2 = i^2 - 2ik + k^2` and `cos((i-k)theta) = Re e^{i(i-k)theta}`:

```
sum_{i,k} c_i c_k i^2 e^{i(i-k)theta} = S_2(c) conj(S_0(c)),
sum_{i,k} c_i c_k k^2 e^{i(i-k)theta} = S_0(c) conj(S_2(c)),
sum_{i,k} c_i c_k i k e^{i(i-k)theta} = S_1(c) conj(S_1(c)) = |S_1(c)|^2.
```

Taking real parts and summing with coefficients `1, 1, -2` gives (D.2). ∎

**Corollary D.2 (the compression on the exact kernel; T1).** On `E = ker T_0` one has
`S_0(c) = 0` by Corollary B.3, so `c^T B c = -2|S_1(c)|^2`, and therefore

```
(P_E B P_E)|_E = -2 * P_E ( u u^T + v v^T ) P_E |_E,
   u_i = i cos(i theta),    v_i = i sin(i theta),    i = 0..R.                    (D.3)
```

**Coefficient:** exactly `-2`. **Domain:** `E subset R^{R+1}`, `dim E = R+1-2g`.
**Rank:** at most `2`. **Operator norm:**

```
|| (P_E B P_E)|_E || = 2 * lambda_max( P_E (u u^T + v v^T) P_E )
                     <= 2 ( ||u||^2 + ||v||^2 ) = 2 sum_{i=0}^R i^2
                     = R(R+1)(2R+1)/3,                                            (D.4)
```

using `cos^2 + sin^2 = 1`.

**Remark D.3 (the rank is 2, not 1).** `|S_1(c)|^2 = <c,u>^2 + <c,v>^2` is a sum of **two**
real rank-one forms on the real coefficient space. A modulus square is not automatically
one real rank-one operator: the real and imaginary directions `u` and `v` must both be
retained unless a symmetry genuinely removes one. Here the symmetry does so exactly when
`sin(i theta) = 0` for all `i`, i.e. `theta in {0, pi}`. Verified numerically
(`src/stage1_rank.py`): rank exactly `2` at every tested `(curve, R)` with
`theta in (0,pi)` — `E5` and `H2F3` and `G3F5` at `R = 2g+2` and `R = 2g+5`, with
`dim E` from `3` to `6` — and rank exactly `1` at `theta = 0`.

**Remark D.4 (`P'(z)` versus `i z P'(z)`).** Since `P_c'(z) = sum_i i c_i z^{i-1}`,

```
S_1(c) = z P_c'(z) |_{z = e^{i theta}},      d/dtheta [ P_c(e^{i theta}) ] = i S_1(c).
```

So `|S_1(c)| = |P_c'(e^{i theta})|` — the **moduli** agree — but `S_1(c)` is not
`P_c'(e^{i theta})`: they differ by the unimodular factor `e^{i theta}`, and the angular
derivative differs further by `i`. The distinction is immaterial for the modulus and
**material** for the decomposition (D.3), because the phase rotates the real and
imaginary parts into each other and hence changes `u` and `v` individually.

**Theorem D.5 (degenerate projection law; T1).** Assume

- **(H1)** `T_0 ⪰ 0` with `E := ker T_0 != {0}` — by B.1 and B.3 this holds exactly when
  `R >= 2g`, and then `dim E = R+1-2g`;
- **(H2)** *complementary gap:* `gamma_0 := min { lambda in spec(T_0) : lambda > 0 } > 0`
  — automatic here, `T_0` having rank `2g`;
- **(H3)** *analyticity:* `a -> T_a` is a real-analytic family of real symmetric matrices,
  even in `a` — established in (C.2).

Then, as `a -> 0`,

```
lambda_min(T_a) = a^2 * lambda_min( (P_E B P_E)|_E ) + O(a^4)
                = -2 a^2 * max_{c in E, ||c||=1} |S_1(c)|^2 + O(a^4).             (D.5)
```

*Proof sketch.* By (H3) and Rellich's theorem for analytic symmetric families, the
eigenvalue branches of `T_a` are real-analytic in `a` near `0`. The `dim E` branches
emanating from the eigenvalue `0` have expansions whose first-order coefficients are the
eigenvalues of the compression of the first-order perturbation to `E`; here the
first-order perturbation in the variable `a^2` is `B`, so those coefficients are the
eigenvalues of `(P_E B P_E)|_E`. Evenness of the family in `a` forces every odd-order
coefficient to vanish, so the next correction is `O(a^4)`. The remaining `2g` branches
start at the positive eigenvalues of `T_0`, which are bounded below by `gamma_0 > 0`, so
for `a` small they do not attain the minimum. Taking the smallest branch gives (D.5), and
Corollary D.2 converts `lambda_min((P_E B P_E)|_E)` into `-2 max_E |S_1|^2`. ∎

**Remainder dependence.** The `O(a^4)` constant is controlled by two mechanisms: the
diagonal quartic term `P_E C P_E` of (C.2), and the second-order coupling to the
complement, of size `|| P_E B P_{E-perp} ||^2 / gamma_0`. It therefore degrades as the
complementary gap `gamma_0` shrinks and as the window degree `R` grows, since
`||B|| = O(R^2)` by (D.4). The statement (D.5) is an asymptotic law in `a` at fixed
`(R, theta, curve)`; it is **not** uniform in `R`.

**Remark D.6 (the genuine quadruple family).** For the conjugation-closed quadruple of
Proposition C.1(2) the increment is twice (C.1), so the same proof gives
`lambda_min = -4 a^2 max_E |S_1|^2 + O(a^4)`. Verified to the same relative accuracy
(`src/stage1_rank.py`, panel A).

### Degenerate cases, kept separate

**(i) `E = {0}`.** If `R < 2g` then `T_0 ≻ 0` and (D.5) does not apply. Instead
`lambda_min(T_a) = lambda_min(T_0) + O(a^2)` by ordinary non-degenerate perturbation
theory: a first-order shift of a strictly positive eigenvalue, with no projection
interpretation. Covered by `src/stage1_belowd.py`.

**(ii) Vanishing compressed coefficient.** If `max_{c in E} |S_1(c)|^2 = 0`, i.e.
`P_E u = P_E v = 0`, then (D.5) degenerates to `lambda_min(T_a) = O(a^4)`.
**No nonzero fourth- or sixth-order coefficient may be inferred from the vanishing of the
quadratic one.** Determining the `a^4` coefficient requires the next order of degenerate
perturbation theory on `E` — the Schur-complement term
`-P_E B P_{E-perp} (T_0|_{E-perp})^{-1} P_{E-perp} B P_E` together with `P_E C P_E` — and
that combination can itself vanish. No claim is made about this case; it was not observed
in any tested cell.

---

## D'. Regression against the paper's Theorem 6 (T1 + T2)

The paper's Theorem 6 moves a single **real** reciprocal pair `{rho, rho^{-1}}`,
`rho = e^a`, and states

```
lambda_min(T) = -2 a^2 sum_{i=0}^{R} (i - R/2)^2 + O(a^4)
              = -binom(R+2,3) a^2 + O(a^4).                                       (D'.1)
```

That family is **exactly the `theta = 0` member of (C.1)**: with `theta = 0`,

```
t_0(n) = 2 for all n   =>   T_0 = 2 * 1 1^T,  rank 1,  E = { c : sum_i c_i = 0 },
Delta t(n) = 2 (cosh(na) - 1)   = the increment of (C.1) at theta = 0.
```

Here `v_i = i sin(0) = 0`, so the rank drops to `1` (Remark D.3) and `u_i = i`. Since
`P_E` is the projection off the all-ones vector,

```
(P_E u)_i = i - (1/(R+1)) sum_{k=0}^R k = i - R/2,
max_{c in E, ||c||=1} |S_1(c)|^2 = ||P_E u||^2 = sum_{i=0}^{R} (i - R/2)^2
                                 = R(R+1)(R+2)/12,
```

the last step being the paper's `S_2 - S_1^2/(R+1) = R(R+1)(R+2)/12` (notes/proofs.md §A).
Theorem D.5 then gives

```
lambda_min(T_a) = -2 a^2 * R(R+1)(R+2)/12 + O(a^4)
                = -R(R+1)(R+2) a^2 / 6 + O(a^4)
                = -binom(R+2,3) a^2 + O(a^4),                                     (D'.2)
```

**recovering Theorem 6's coefficient exactly.** In the paper's `eps = rho - 1`
normalisation (Remark A6 of `notes/proofs.md`) the same statement reads
`-R(R+1)(R+2) eps^2 / 6 + O(eps^3)`, the regression target of the revision brief.

So the paper's "centred second moment of the window", `sum_i (i - R/2)^2`, **is** the
squared norm of the projected derivative-evaluation functional `||P_E u||^2`. The two
derivations are independent: Theorem 6 computes `lambda_min` of the rank-two matrix
`u w^T + w u^T` directly, whereas Theorem D.5 reaches the same number through the
compression of `B` to `ker T_0`. Checked at `R = 1, 2, 4, 8, 16, 32, 64` in 60-digit
arithmetic, agreement exact to `< 10^-50` (`src/stage1_regression.py`, panel 1).

**Why other members of the family have different coefficients.** For `theta in (0,pi)` the
coefficient is `-2 max_{c in E} |S_1(c)|^2`, which depends on `theta` and on `E` and is
*not* `binom(R+2,3)`: the rank is `2` rather than `1`, and `E` is cut out by `2g`
conditions rather than one. The regression therefore confirms the general law at its one
degenerate member; it does not make the coefficient universal. In particular **no universal
cubic scaling in `R` is claimed** for `theta != 0`.

---

## D''. What "the error scales like `eps^2`" means (T1 for the definition, T2 for the values)

The preregistered kill criterion **K1** (`extensions/projection-formula/PREREGISTRATION.md`,
committed `8ff0545` on 2026-10-07 16:16:42 -0400, **alone and before any computation**;
the code and data commit `fb06e78` is 20 minutes later) reads verbatim:

> **K1.** The relative error of prediction 2 must `-> 0` like `a^2` as `a -> 0`, at every
> `R >= d`.

So the measured quantity is the **relative error of the predicted minimum eigenvalue**,

```
rel(a) := | lambda_min(T_a) - pred(a) | / | pred(a) |,
pred(a) := -2 a^2 max_{c in E, ||c||=1} |S_1(c)|^2.                               (D''.1)
```

Because `pred(a) = Theta(a^2)` while the true residual is `Theta(a^4)` by Theorem D.5,
`rel(a) = Theta(a^2)`. It is therefore:

- **not** the raw residual `|lambda_min - pred|`, which is `Theta(a^4)`;
- **not** the raw residual divided by `a^2`, which is `Theta(a^2)` but carries the units of
  an eigenvalue;
- it **is** the dimensionless relative error of the `a^2` coefficient.

**Exact anchor.** At `R = 1`, `theta = 0`, the paper's Remark A7 gives the closed form
`lambda_min = -4 sinh^2(a/2)` exactly, while `pred = -a^2`. Hence

```
rel(a) = | 4 sinh^2(a/2) - a^2 | / a^2 = a^2/12 + O(a^4),     rel(a)/a^2 -> 1/12.
```

Measured `rel(a)/a^2 = 0.0833...` at `a = 10^-2, 10^-3, 10^-4` — matching `1/12` to the
printed digits. The analogous measured constants at `theta = 0` are

| `R` | 1 | 2 | 4 | 8 | 16 | 32 |
|---|---|---|---|---|---|---|
| `rel(a)/a^2` | 0.0833 | 0.3333 | 1.1333 | 3.9333 | 14.3333 | 54.3333 |

stable across the three decades of `a` (`src/stage1_regression.py`, panel 2). The growth
in `R` is the `||B|| = O(R^2)` remainder dependence noted after Theorem D.5, and is the
reason (D.5) must not be read as uniform in `R`.

**On the 44-case function-field check.** `src/stage1_ff.py` reports `rel` in the sense of
(D''.1) for 4 curves × 11 window degrees × 3 values of `a`. That is the quantity K1
names, and it is the quantity tabulated in `data/stage1_ff.json`. The `dps = 120`
repeat (`data/stage1_ff_dps120.json`) exists to confirm that the reported `rel` is a
property of the mathematics and not of the arithmetic.

---

## E. The zeta-side overlap model (T2 — a numerical model, not a theorem)

Everything in this section is **numerical**. It is **not** a continuum theorem, **not** a
certified positivity statement, and **not** an established prolate projection law. The
model is defined first, in full, and only then compared with measurement.

### E.1 Definition of the model

**Basis family.** `EVEN_D`: `f_k(x) = cos((2k+1) pi x / (2L))`, `k = 0..N-1` — the
quarter-wave even basis, which vanishes at `x = pm L`. This is the basis of the paper's
§4.3 and of `src/stage2_zeta.py`. It is **not** the `EVEN` basis `cos(k pi x / L)` used by
the parity control of §G, and it is not the `ODD` basis. The near-null band is
basis-dependent, so the basis is part of the model.

**Dimensions and precision.** One cell per window, exactly as in §4.3:

| `L` | 0.8 | 1.0 | 1.3 | 1.6 | 2.0 |
|---|---|---|---|---|---|
| `N` | 16 | 18 | 20 | 22 | 26 |
| `mpmath` dps | 70 | 80 | 90 | 100 | 120 |

**Gram matrix and orthonormalisation.** For `EVEN_D`, `int_{-L}^{L} f_j f_k dx = L delta_jk`
exactly, so the basis is already orthogonal and `norm2(k) = L`. Orthonormalisation is
division by `sqrt(norm2(k))` — exact, with no Gram–Schmidt and no ill-conditioning.

**The three forms, kept distinct.** Writing `gamma = gamma_1 = 14.1347...` and
`F(t;c) = int f_c(x) e^{itx} dx` for the coefficient vector `c`:

- `Q_zeta,N` — the `N`-dimensional Galerkin matrix of the **actual** zeta Weil form on the
  window, assembled from the geometric side (`build_matrix`);
- `Q_0,N := Q_zeta,N + 2 F(gamma)^2` — the **augmented counterfactual baseline** of the
  paper's eq. (3), realised as `build_matrix + move_block(delta=0)`; `move_block` at
  `delta = 0` equals `+2 p p^T` with `p` the representer of `c -> F(gamma;c)`;
- `Q_rest := Q_0,N - 4 F(gamma)^2 = Q_zeta,N - 2 F(gamma)^2` — the **remaining-zero form**
  of the paper's Theorem 8.

These are three different matrices and are **nowhere identified**. The sweep perturbs and
measures `Q_0,N` throughout: the baseline eigenvalue used in `C(delta)` is
`lambda_min(Q_0,N)`, never `lambda_min(Q_zeta,N)`.

> **Observed coincidence, not a substitution.** `lambda_min(Q_0,N)` and
> `lambda_min(Q_zeta,N)` agree to all 8 printed digits at every `L`, because the window
> minimiser very nearly annihilates `F(gamma_1)` so the rank-one bump `2F(gamma)^2` barely
> moves the floor. This is an empirical observation about these cells. **The minimum of
> `Q_zeta,N` is not substituted for the minimum of `Q_0,N` anywhere**, and the agreement is
> not used to justify any step.

**Near-null cutoff.** For each `delta`,

```
tau(delta) := 4 K(L,gamma) delta^2,
E_(N,tau) := span { eigenvectors of Q_0,N with eigenvalue <= tau(delta) }.
```

This is a **cutoff-defined near-null band**, not a kernel. `Q_0,N` is strictly positive
definite in every cell (`lambda_min` from `2.5e-17` at `L=0.8` down to `1.1e-67` at
`L=2.0`); the band has no exact null vector, and is never called one. Contrast §B, where
`ker T_R` is exact and algebraic.

**Displacement range.** `delta = 10^-1, 10^-2, ..., 10^-12` (12 decades), plus
`delta = 0.02`, the paper's §4.3 value, used for the preregistered criterion K2.

**Coefficient extraction.** `C(delta) := (lambda_min(Q_0,N) - lambda_min(Q_delta,N))/delta^2`,
a **single finite-difference ratio at each `delta` independently** — not a fitted slope,
not a regression over `delta`, and with no free parameter.

**Overlap ratio.** With `v` the representer of `c -> F'(gamma;c)` (so
`<v,c> = F'(gamma;c)`),

```
r_(N,tau) := || Pi_{E_(N,tau)} v ||^2 / || Pi_basis v ||^2.                        (E.1)
```

> **Deviation from preregistration, retained deliberately.** The preregistration
> (`PREREGISTRATION.md`, line 19) defines
> `cos^2 theta_N = || Pi_N r_gamma ||^2 / || r_gamma ||^2` with the **continuum** norm
> `|| r_gamma ||^2 = K(L,gamma)`. The implementation divides instead by
> `nr2 = sum_a v[a]^2 = || Pi_basis v ||^2`, the squared norm of the **finite-basis
> projection** of the representer. Per the revision brief the implemented definition is
> **retained**, and the omitted tail is reported separately rather than the denominators
> being interchanged:

| `L` | `N` | `||Pi_basis v||^2 / K` | omitted tail | prediction inflation `1/ratio` |
|---|---|---|---|---|
| 0.8 | 16 | 0.93432600 | 6.567e-2 | 1.0703 |
| 1.0 | 18 | 0.92998416 | 7.002e-2 | 1.0753 |
| 1.3 | 20 | 0.98716931 | 1.283e-2 | 1.0130 |
| 1.6 | 22 | 0.97845048 | 2.155e-2 | 1.0220 |
| 2.0 | 26 | 0.99999619 | 3.808e-6 | 1.0000 |

> **Consequence (material).** The prediction as computed is
> `pred = 4K * (proj2/nr2)`, which mixes the *continuum* amplitude `4K` with the
> *finite-basis* denominator and is therefore inflated by `1/ratio` — **up to 7.5 % at
> `L = 1.0`**, i.e. larger than the headline agreement it is being compared against. The
> consistently normalised prediction is `4 * proj2 = 4K * proj2 / K`. Both are tabulated in
> §E.3 and in `data/stage5_audit_norm.json`; neither is presented as the other.

### E.2 The full second-order perturbation

`move_block(delta)` implements the planted displacement of the zero at `gamma`. Expanding,

```
Q_delta - Q_0 = delta^2 P + O(delta^4),
P[a,b] = -4 ( v[a] v[b] + ( p[a] s[b] + p[b] s[a] )/2 ),                           (E.2)
```

where `p, v, s` are the representers of `c -> F(gamma;c)`, `F'(gamma;c)`, `F''(gamma;c)`.
As a quadratic form this is

```
c -> -4 delta^2 [ F'(gamma;c)^2 + F(gamma;c) F''(gamma;c) ],                       (E.3)
```

i.e. the **full** second-order perturbation, `F F''` term included. On the proved side the
manuscript already carries this term: Theorem 7's budget bounds
`4|F'(gamma)^2 + F(gamma) F''(gamma)| <= c(L)` by Cauchy–Schwarz, and Theorem 8 moves the
`F F''` allowance into `O(delta^4)` by completing the square
`4F^2 - 4 delta^2 F F'' = (2F - delta^2 F'')^2 - delta^4 F''^2`.

Two predictions are therefore distinguished throughout, and **the `F F''` term, the
baseline variation, the coupling to the complement and the fourth-order remainder are not
discarded in order to obtain a squared norm**:

- **Squared-norm heuristic (H1):** `pred_H1 = 4 K(L,gamma) * r_(N,tau)`. Keeps only the
  `v v^T` part of (E.2). **Discards** the `F F''` term, the variation of the baseline
  eigenvalues inside the band, the coupling to `E-perp`, and the `delta^4` remainder.
- **Compressed operator:** `C_comp = (lambda_min(Q_0,N) - lambda_min(Lambda_E + delta^2 Pi_E P Pi_E))/delta^2`,
  where `Lambda_E` carries the exact baseline eigenvalues of the band. **Retains** the
  `F F''` term, the baseline variation, and the within-band coupling. **Omits** only the
  coupling to `E-perp` and the `delta^4` remainder.

### E.3 Measured agreement, and what the "≈ 3 %" refers to

Over the 56 `(L, delta)` cells with `dim E_(N,tau) >= 1`, relative error `|C - pred|/|C|`:

| prediction | median | mean | within 3 % |
|---|---|---|---|
| H1 squared-norm, as implemented (`/nr2`) | **4.06 %** | 5.278 | 25 / 56 |
| H1 squared-norm, consistently normalised (`/K`) | **1.94 %** | 4.898 | 35 / 56 |
| compressed operator (retains `F F''`) | **0.23 %** | 0.0148 | 48 / 56 |

**The "≈ 3 %" agreement is a property of the compressed second-order operator, not of the
squared-norm law.** Quoting `≈ 3 %` for `C = 4K cos^2 theta` *overstates* that law's
accuracy; quoting it for the compressed operator *understates* it. The committed
`extensions/projection-formula/REPORT.md` reports `3.02e-2` as a median over the 51 of 58
Stage-5 cells with `C delta^2 >> lambda_floor`; that is a different denominator from the
table above and the two are not interchangeable. Both denominators are now stated
explicitly wherever a percentage appears.

**Exclusions, and relative errors near zero crossings.** The *means* above are dominated by
two cells in which `C(delta)` has already collapsed while the preregistered cutoff still
admits a two-dimensional band:

| `L` | `delta` | `dim E` | `C` | `rel` (H1) | `rel` (compressed) |
|---|---|---|---|---|---|
| 0.8 | `1e-5` | 2 | `1.390e-4` | `1.287e+02` | `3.708e-05` |
| 1.0 | `1e-10` | 2 | `6.456e-5` | `1.627e+02` | `3.651e-06` |

Both are the last `delta` before `dim E` drops: **cutoff lag, not a failure of the
mechanism** — the compressed prediction is accurate to `10^-5`–`10^-6` at those same
cells. A further four cells (`L = 0.8`, `delta <= 1e-9`) have `dim E = 0`, where the
prediction is identically zero and the relative error is `1` by construction; they are
excluded from all three columns. Medians are reported for exactly this reason: a relative
error taken near a zero crossing of the prediction is not a meaningful quantity, and no
percentage in this work is computed across one.

### E.4 The `delta`-staircase claim — WITHDRAWN as uninformative

The committed report records `K3: 11 of 12 C-steps align with an eigenvalue scale within a
factor 3 -> PASS`. On audit this claim does not survive, and it is **excluded from the
upload-ready claims**.

**Event definition, as implemented** (`src/stage2b_k2k3.py`):

- `delta` is swept on a **decade grid** `10^-1 ... 10^-12`, giving **11** consecutive-decade
  brackets per window, 55 in total.
- A *C-step* is a bracket `(d_{i-1}, d_i)` with
  `|C(d_{i-1}) - C(d_i)| / max(|C(d_{i-1})|, 10^-300) > 0.20`.
- The denominator **12 is the observed number of such steps** summed over the five windows
  (`4 + 5 + 2 + 1 + 0`). It is **not** a preregistered list of twelve expected steps.
- A step is scored *aligned* if some `delta_k = sqrt(lambda_k / 4K)` lies in
  `[min(lo,hi)/3, max(lo,hi)*3]`. With `hi = lo/10` that window is `[lo/30, 3 lo]`: a
  factor **90** in `delta`, i.e. a factor **8100** in `delta^2`. The preregistered K3
  criterion was "within a factor 3 in `delta^2`" — **27 times tighter** than what was run.

**Null model.** The test is informative only if a randomly placed step would usually fail
it. Scoring all 55 available decade brackets (`src/stage5_audit_staircase.py`):

| `L` | `N` | #`delta_k` | #C-steps | aligned (implemented) | aligned (preregistered) | null rate (impl.) | null rate (prereg.) |
|---|---|---|---|---|---|---|---|
| 0.8 | 16 | 16 | 4 | 4/4 | 3/4 | 7/11 = 64 % | 6/11 = 55 % |
| 1.0 | 18 | 18 | 5 | 4/5 | 2/5 | 9/11 = 82 % | 5/11 = 45 % |
| 1.3 | 20 | 20 | 2 | 2/2 | 2/2 | 11/11 = 100 % | 9/11 = 82 % |
| 1.6 | 22 | 22 | 1 | 1/1 | 1/1 | 11/11 = 100 % | 11/11 = 100 % |
| 2.0 | 26 | 26 | 0 | 0/0 | 0/0 | 11/11 = 100 % | 10/11 = 91 % |
| **total** | | | **12** | **11/12 = 92 %** | **8/12 = 67 %** | **49/55 = 89 %** | **41/55 = 75 %** |

**Verdict.** Under the implemented criterion, `11/12 = 92 %` sits against a null pass rate
of `89 %`: the test carries essentially no information. Under the preregistered
factor-3-in-`delta^2` criterion the result is `8/12 = 67 %`, **below** its own `75 %` null.
The cause is structural: with 16–26 eigenvalues of `Q_0,N` spread across the sampled
decades, almost any bracket contains some `delta_k`. **H2 is not supported by this test**,
and no staircase claim is carried into the manuscript. This is a correction to the
committed report, recorded here rather than by silently editing it; what is withdrawn is
the *alignment claim*, not the raw `C(delta)` values, which stand.

---

## F. The prolate question (Section 6) — what the sources actually say

**Status of each source.**

| source | status | what it says about a prolate / near-null identification |
|---|---|---|
| CCM arXiv:2310.18423, *Zeta zeros and prolate wave operators* | **CHECKED** (full text) | Introduces a **semilocal** analogue of the prolate wave operator; the archimedean prolate operator is the square of the scaling operator plus the grading of orthogonal polynomials. Positive part realises low-lying zeros, negative part is the Sonin space. No identification of a near-null subspace of a *truncated window Weil form* with a prolate space. |
| CCM arXiv:2511.22755, *Zeta Spectral Triples* | **CHECKED** (full text) | §7, verbatim: *"the observation of [4] that the eigenfunction associated with the lowest eigenvalue of `QW_lambda` is well approximated by prolate spheroidal wave functions."* `[4]` = Connes–Consani, *Spectral triples and ζ-cycles*, Enseign. Math. **69** (2023) 93–148. The prolate wave operator there is `PW_lambda := -d_x (lambda^2 - x^2) d_x + (2 pi lambda x)^2`, a deformation of the harmonic oscillator, whose eigenfunctions are eigenfunctions of the compression of the Fourier transform by the projection `P_lambda` onto `[-lambda, lambda]`. |
| Connes arXiv:2602.04022, survey | **CHECKED** (full text) | §6.4, verbatim: *"In [25], we gave a construction numerically justified of the eigenvectors associated to the first minuscule eigenvalues of `QW_lambda`, using prolate spheroidal wave functions associated to the interval `[-lambda, lambda]`. In particular this gives an educated guess for an approximation of the eigenvector associated to the smallest eigenvalue `eps(lambda)` of `A_lambda`."* The construction is `k_lambda = E(h_lambda)` with `h_lambda` the combination of `h_{0,lambda}, h_{4,lambda}` of vanishing integral, on which *"`QW_lambda` takes non-zero, but extremely small values"*. §6.6 lists as a **remaining step** that the smallest eigenvalue of `QW_lambda` be shown *simple with even eigenvector* — known only *"for the prolate wave operator"*. |
| Zhu arXiv:2608.24827 | **CHECKED** (full text) | Prolate enters only as the **rate constant**: *"`pi^2 / ln` is the Landau–Widom rate governing the plunge of the eigenvalues of time–band limiting (prolate) operators beyond the Shannon number"*, in the empirical law `-ln lambda_*(L) = C N(T*)/ln N(T*) (1+o(1))`, `C = 20.13... ≈ 2 pi^2`, `T* = 2 pi e^{2L}`. No subspace identification. |
| Connes–Consani, Enseign. Math. **69** (2023) 93–148 (`[4]`/`[25]` above — the actual origin of the observation) | **NOT CHECKED** | not freely accessible |
| Connes–Moscovici, PNAS **119** (2022) e2123174119 | **NOT CHECKED** | — |
| Connes–Consani, Selecta Math. **27** (2021) Paper 77 | **NOT CHECKED** | — |
| CCM, *Riemann Zeros via Weil Forms: From Prolate Functions ...* (ref. [31] of the survey) | **NOT CHECKED** | no arXiv identifier recovered |

**Conclusion (T3).** Prior art **does** exist for "prolate spheroidal wave functions
approximate the lowest eigenvector(s) of a truncated Weil form" — Connes–Consani, reported
in Connes's 2026 survey §6.4 and in CCM arXiv:2511.22755 §7. It is explicitly
**numerically justified** and described by its authors as an **educated guess for an
approximation**; it is not a proved identification, and even the simplicity of the lowest
eigenvalue is an open remaining step. It concerns the **multiplicative** window
`[lambda^-1, lambda]` and the operator `QW_lambda` built from the map `E` and the near
intersection of `P_lambda` with `P_lambda-hat` — a different construction from the
**additive** window `[-L,L]` Galerkin form `Q_0,N` in a cosine basis used here. No
identification between the two objects is established in any checked source.

Accordingly:

- The near-null object of this work is described as **"projection onto the computed
  near-null subspace `E_(N,tau)` of `Q_0,N`"**, with its cutoff stated. The phrase "the
  near-null space is the prolate space" is **not** used.
- The prolate quantities `1 - chi_k` were **not computed** (Stage 4, deviation 4). They are
  the right objects — Connes's survey Figure 1 plots `log(1 - chi_2(sqrt x))` — and
  computing them at the `10^-50` scale needs dedicated prolate code. Their absence is a
  scope limitation, stated as such.
- **Absence of a statement from these sources is not proof of novelty**, particularly with
  four relevant sources marked NOT CHECKED.

> **Correction to a committed claim.** `extensions/projection-formula/REPORT.md`, Stage 6
> table, attributes to CCM arXiv:2511.22755 the statement *"the space of eigenvectors of
> the `k` lowest eigenvalues of `QW_lambda` corresponds to the prolate projection
> `Pi(lambda,k)`"* and concludes *"the near-null space **is** prior art"*. The primary
> source does not say this. It says (i) the **single** lowest eigenfunction, not a
> `k`-dimensional eigenspace; (ii) *"is well approximated by"*, not *"corresponds to"*;
> (iii) it is an **observation attributed to [4]**, not a proved identification. The
> notation `Pi(lambda,k)` does not occur in that paper. Recorded here rather than by
> silently editing the committed report; the overstated row is not carried forward.

---

## G. Parity: the negative control (Section 7)

### G.1 What Zhu actually claims, and what he actually leaves open

Zhu arXiv:2608.24827 §6, **verbatim**:

> "Second, the parity asymmetry is systematic. [...] (ii) The bottom of the window spectrum
> alternates parity, `lambda^even_1 < lambda^odd_1 < lambda^even_2 < lambda^odd_2` at every
> scanned `L`, the same alternation as the prolate eigenfunctions of time–band limiting.
> (iii) The sector ratio `lambda^odd_1 / lambda^even_1` grows almost exactly geometrically,
> `≈ 10^(2.07 L + 1.32)` over the scanned range (residuals of `log10` below 0.06). Its
> mechanism is not the pole sign, which in fact favors the odd sector (the pole term adds
> `+2c^2` to the even form and `-2s^2` to the odd one); on the zeros side the asymmetry
> must instead reflect the constraint `F_o(0) = 0`, which denies odd test functions the
> zero-free frequency interval `[0, gamma_1)` in which even minimizers park most of their
> mass. **We leave its quantitative law, and whether the odd sector obeys a shifted version
> of (5), to a separate study.**"

Three things follow, and they must be kept apart.

1. The `F_o(0) = 0` reading is **Zhu's own proposed mechanism**, asserted as *"must
   instead reflect"* — a stated expectation, not a proved claim.
2. The open question Zhu **actually** leaves is *"its quantitative law, and whether the odd
   sector obeys a shifted version of (5)"* — i.e. the **law of the sector ratio** and the
   **Landau–Widom behaviour of the odd floor**. It is **not** "why does parity matter".
   Attributing the latter to Zhu would misstate the source.
3. Zhu's ratio law reproduces his own Table 1 to `max |residual in log10| = 0.0590`,
   consistent with his stated "below 0.06" — verified independently here.

**Zhu Table 1** (his words: *"reference values, not certificates"*; converged exact
assembly in a Legendre basis per sector):

| `L` | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 | 1.0 | 1.1 | 1.2 |
|---|---|---|---|---|---|---|---|---|
| `lambda^even_1` | 9.34e-7 | 1.61e-9 | 4.18e-13 | 1.65e-17 | 4.14e-23 | 5.88e-30 | 2.04e-38 | 7.94e-49 |
| `lambda^odd_1` | 1.94e-4 | 5.97e-7 | 2.56e-10 | 1.57e-14 | 7.22e-20 | 1.49e-26 | 7.68e-35 | 4.76e-45 |
| ratio | 208 | 371 | 612 | 949 | 1746 | 2535 | 3769 | 5997 |

### G.2 The control, revalidated (T2)

The control: in the `EVEN` basis `cos(k pi x / L)` one has `int f_k = 0` for every `k >= 1`
and `2L` for `k = 0`, so imposing `F(0) = int f = 0` is exactly *deleting the constant
mode* (confirmed numerically, `max_{k>=1} |F_k(0)| <= 2.2e-81`). The constrained-even floor
is then compared with the odd floor.

The original run used a single `N` per window. It is revalidated here at three `N` per
window and checked against Zhu's reference values (`src/stage3_revalidate.py`):

| `L` | `N` | even floor | constrained even | odd floor | constr/odd | even / Zhu | odd / Zhu |
|---|---|---|---|---|---|---|---|
| 0.5 | 14 | 1.109e-6 | 1.256e-2 | 2.386e-4 | 52.6 | 1.188 | 1.230 |
| 0.5 | 20 | 1.014e-6 | 1.227e-2 | 2.301e-4 | 53.4 | 1.086 | 1.186 |
| 0.5 | **26** | 9.785e-7 | 1.217e-2 | 2.241e-4 | **54.3** | 1.048 | 1.155 |
| 0.6 | 14 | 2.330e-9 | 8.478e-5 | 7.077e-7 | 119.8 | 1.447 | 1.185 |
| 0.6 | 20 | 2.014e-9 | 7.541e-5 | 7.006e-7 | 107.6 | 1.251 | 1.174 |
| 0.6 | **26** | 1.864e-9 | 7.143e-5 | 6.917e-7 | **103.3** | 1.158 | 1.159 |
| 0.8 | 16 | 3.199e-17 | 1.246e-11 | 2.586e-14 | 481.9 | 1.939 | 1.647 |
| 0.8 | 22 | 2.612e-17 | 7.241e-12 | 2.319e-14 | 312.2 | 1.583 | 1.477 |
| 0.8 | **28** | 2.113e-17 | 6.963e-12 | 2.096e-14 | **332.1** | 1.281 | 1.335 |
| 1.0 | 18 | 5.270e-27 | 1.680e-21 | 1.882e-24 | 892.6 | **896.2** | **126.3** |
| 1.0 | 24 | 3.595e-29 | 5.269e-23 | 5.740e-26 | 917.9 | 6.114 | 3.852 |
| 1.0 | **30** | 1.346e-29 | 2.870e-23 | 3.292e-26 | **871.6** | 2.289 | 2.210 |

**Two findings.**

1. **The `L = 1.0` row of the original run was badly unconverged.** At `N = 18` its even
   floor sits `896x` above Zhu's reference value and its odd floor `126x` above; at
   `N = 30` the factors are `2.3` and `2.2`. Every floor approaches Zhu's values **from
   above**, as finite-`N` variational upper bounds must. The absolute floors from the
   single-`N` run are therefore **not** upload-ready numbers and are not quoted as such.
2. **The ratio is far better conditioned than the floors, and the verdict is robust.**
   `constr/odd` moves only from `892.6` to `871.6` at `L = 1.0` across three orders of
   magnitude of floor convergence. At the largest tested basis the range is
   **54x to 872x**, against a preregistered tolerance of `3x`. The originally reported
   "53x to 893x" came from the smallest `N` per window; the range is `N`-dependent, the
   failure is not.

### G.3 A second, sharper control: the mass test (T2, new)

Zhu's stated mechanism is not the point condition `F_o(0) = 0` alone but **access to the
interval** `[0, gamma_1)`: odd test functions are *"denied the zero-free frequency interval
`[0, gamma_1)` in which even minimizers park most of their mass"*. Our first control
imposes only the **single linear condition** `F(0) = 0`, which is strictly weaker than
oddness. So the interval reading was tested directly, by measuring for each minimiser

```
m := int_0^{gamma_1} |F(t)|^2 dt / int_0^infty |F(t)|^2 dt,     denominator = pi
```

(Plancherel: `F_k(t) = int f_k e^{itx} dx` and `||f||_2 = 1` give
`int_R |F|^2 = 2 pi`; `|F|^2` is even. Verified numerically to `1e-4` in the even sector and
`1e-8` in the odd — the even-sector limit is the slow tail of the quadrature, not the
identity.)

| `L` | `N` | even minimiser | constrained-even minimiser | odd minimiser |
|---|---|---|---|---|
| 0.5 | 26 | 0.99987 | 0.96432 | 0.99668 |
| 0.6 | 26 | 0.99995 | 0.99077 | 0.99894 |
| 0.8 | 28 | 0.99999 | 0.99755 | 0.99970 |
| 1.0 | 30 | 0.99999 | 0.99882 | 0.99986 |

**The odd minimiser parks 99.67 %–99.99 % of its `|F|^2` mass on `[0, gamma_1)` —
essentially as much as the even minimiser.** Within the tested range the odd sector is
therefore **not** denied that interval in the mass sense, so the mass-parking reading is not
what separates the sectors here. The reason is structural and is a scope limitation, not a
refutation of Zhu's intuition at large `L`: for `L <= 1.0` the basis frequencies
`k pi / L <= ~ 10 pi` and `gamma_1 = 14.13`, so `[0, gamma_1)` contains nearly all the
spectral mass of **both** parities, and cannot discriminate between them.

### G.4 Scoped conclusion

> **Within the tested setting, imposing `F(0) = 0` on the even sector does not reproduce the
> odd-sector floor; the single-condition explanation is unsupported.**

Additionally, within the tested range `L <= 1.0`, the mass-parking refinement of that
explanation is not the discriminator either: both parities concentrate essentially all
their `|F|^2` mass on `[0, gamma_1)`.

**What this does not show.**

- It does **not** show that parity is irrelevant. The sectors differ by two to four orders
  of magnitude and follow Zhu's geometric ratio law `≈ 10^(2.07 L + 1.32)` closely. Parity
  carries near-null structure that the single linear constraint does not capture; what that
  structure is remains unidentified here.
- It does **not** solve Zhu's open question. What Zhu leaves open is the **quantitative law**
  of the sector ratio and whether the odd floor obeys a shifted Landau–Widom law. Neither
  is addressed, and neither is claimed.
- It is a **finite, basis-dependent observation** at four windows in one even basis and one
  odd basis, with the floors themselves only variational upper bounds.
