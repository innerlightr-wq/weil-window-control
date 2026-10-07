# Exact Reducing Projections versus Perturbative Spectral Projections

*Repository-only working note: mathematical context, not a new RH theorem or a claim of novelty.*

Companion to [`REPORT.md`](REPORT.md) in this directory. Checks run for this note are in
[`PROJECTION_CONTEXT_CHECKS.md`](PROJECTION_CONTEXT_CHECKS.md) / `check_projection_context.py`.

Throughout, **(a)** marks a statement quoted from a source, **(b)** an elementary deduction checked
here, **(c)** a qualification prompted by this discussion. The historical exploration report is not
rewritten; qualifications to it are recorded below as (c).

## The distinction

Two different things get called "a projection", and they answer different questions.

- The **inheritance theorem** concerns a projection that *commutes with the full operator* and
  therefore preserves an entire spectral sector exactly, at every order.
- The **second-moment calculation** concerns a projection *onto a baseline kernel*, which
  determines a leading perturbative response and nothing beyond it.

These are related uses of projection, **not interchangeable theorems**. One splits an operator; the
other reads off a first coefficient.

## 1. The exact reducing projection

**(a) Source.** De Jesús, *Exact Spectral Inheritance in Prime-Power Resolution Towers: A Splitting
Theorem for Affine-Cocycle Transfer Operators, with Application to the Collatz Carry Sum*,
1 October 2026, [doi:10.5281/zenodo.23088511](https://doi.org/10.5281/zenodo.23088511), 9 pp.
Read from the Zenodo PDF (`prime_power_resolution (1).pdf`, 304,292 B, md5 `a731d033…`).

**Definition 2.** `π_e : K_{e+1} ↠ K_e` is reduction mod `p^e`, a surjective homomorphism whose
fibres all have the same size `m_e := |K_{e+1}|/|K_e| = |ker π_e| ∈ {1,p}`. With

```
J_e : H_e -> H_{e+1},   J_e f := f ∘ π_e          (pullback)
S_e : H_{e+1} -> H_e,   (S_e f)(ū) := Σ_{u ∈ π_e^{-1}(ū)} f(u)      (fibre sum)
```

the normalised Hermitian inner products give `⟨J_e f, J_e g⟩_{e+1} = ⟨f,g⟩_e`, so `J_e` is an
isometry, and **`S_e ∘ J_e = m_e · id`**.

**Theorem 6 (pullback intertwining).** `L_{e+1,pt,x} J_e = J_e L_{e,t,x}` for every `e ≥ 1` and
every `t ∈ Z/p^e Z`. Remark 7 records that the proof never uses `gcd(t,p) = 1`: it **lifts a
lower-level residue `t` by one factor of `p`**. The twist restriction `t → pt` is part of the
statement and is preserved here; nothing below extends it to arbitrary primitive fine-level twists.

**Lemma 8 (pushforward intertwining).** `S_e ∘ L_{e+1,pt,x} = L_{e,t,x} ∘ S_e`, proved
independently of Theorem 6.

**Theorem 9 (splitting).** `Z_{e+1} := ker S_e` satisfies `Z_{e+1} = J_e(H_e)^⊥` exactly, is
invariant under `L_{e+1,pt,x}`, and `H_{e+1} = J_e(H_e) ⊕ Z_{e+1}` is an orthogonal direct sum of
invariant subspaces, so **`L_{e+1,pt,x} ≅ L_{e,t,x} ⊕ N_{e+1,pt,x}`**, with
`dim Z_{e+1} = (m_e − 1)|K_e|`. **Corollary 10:** the spectra are a multiset union with full Jordan
structure inherited, and `ρ(L_{e+1}) = max{ρ(L_e), ρ(N_{e+1})}` exactly. **Corollary 11:** the same
for singular values, the adjoint being block-diagonal in the same basis.

**(b) Elementary consequences, checked here.** With the normalised inner products,
`⟨h, J_e g⟩_{e+1} = (1/|K_{e+1}|) Σ_ū (S_e h)(ū) \bar g(ū)` and `|K_{e+1}| = m_e |K_e|`, so

```
J_e* = S_e / m_e ,        P_old := J_e J_e* = J_e S_e / m_e .
```

`P_old² = J_e (S_e J_e) S_e / m_e² = J_e (m_e id) S_e / m_e² = P_old`, and `P_old* = P_old`: it is
the **orthogonal projection onto observables constant on each reduction fibre**. Using *both*
intertwinings, one on each side,

```
P_old L_{e+1} = (1/m_e) J_e (S_e L_{e+1}) = (1/m_e) J_e L_e S_e      [Lemma 8]
L_{e+1} P_old = (1/m_e) (L_{e+1} J_e) S_e = (1/m_e) J_e L_e S_e      [Theorem 6]
=>  [P_old , L_{e+1,pt,x}] = 0 .
```

**Why the second intertwining is essential.** These transfer operators **need not be normal**, and
the source says so explicitly: Theorem 6 alone gives only that `J_e(H_e)` is invariant, hence
block-*upper-triangular* structure, "for which the complementary block need not be invariant". It
is Lemma 8 that makes the orthogonal complement invariant — equivalently, that makes `J_e(H_e)`
invariant under the adjoint too. The check in `check_projection_context.py` includes a control: an
operator satisfying `L J = J L` but not `S L = L S` has `[P_old, L] ≠ 0`.

## 2. The perturbative spectral projection

Write `epsilon` for the perturbation parameter, to avoid collision with the inheritance paper's
arithmetic parameters `e, p, t, x`.

Fix a **finite-dimensional** weighted space with weights `w_j > 0`, let `X` be the real coordinate
observable and `1` the constant observable, and use the weighted inner product throughout. Put

```
(T_eps f)(X_i) = Σ_j 2 cosh(eps (X_i − X_j)) f(X_j) w_j ,
A_eps = G + T_eps ,   G = G* >= 0 ,   A_0 = G + 2|1><1| ,   E = ker(A_0) .
```

Assume `E ≠ {0}` and that `A_0` has a positive spectral gap `γ > 0` on `E^⊥`. Since
`2cosh(z) = 2 + z² + O(z⁴)`,

```
T_eps = T_0 + eps² B + O(eps⁴),   B with kernel (X−Y)² ,
B = |X²><1| + |1><X²| − 2|X><X| .
```

For `v ∈ E` both `G v = 0` and `⟨1,v⟩ = 0` (each term of `A_0` is `⪰ 0`), so `E` lies in the
mean-zero subspace and the first two terms of `B` annihilate, giving

```
P_E B P_E = −2 |P_E X><P_E X| .
```

With `X_c = X − mean_mu(X)` one has `P_E X = P_E X_c` (because `E ⊥ 1`), so the checked
finite-dimensional consequence is

```
lambda_min(A_eps) = −2 eps² ||P_E X_c||² + O(eps⁴) .
```

The remainder depends on the fixed geometry, on `G` and on the gap `γ`; **no uniformity in a
growing window follows**, and none is claimed. If `E = {0}` the statement is different in kind: a
strictly positive baseline stays positive for all sufficiently small `eps`, so no negative
eigenvalue is created at all. An exact kernel is never replaced here by a tiny positive numerical
eigenvalue: every `E` used in the checks is computed by exact rational linear algebra.

**The decisive contrast.** `P_E` commutes with `A_0` by construction, but **generally does not
commute with `A_eps`**. It therefore controls the leading response without separating the perturbed
operator; higher-order interaction with `E^⊥` remains, and indeed supplies the `O(eps⁴)` term.

## 3. The shared second-moment language, kept apart

**(b)** For a centred observable `v`, orthogonality of `P_old` gives immediately

```
||v||² = ||P_old v||² + ||(I − P_old) v||² .
```

Under the tower's uniform probability measure `P_old` is conditional averaging over fibres, so this
separates the variance **visible at coarse resolution** from the within-fibre remainder. That is a
Pythagorean identity, **not a new conservation principle**.

Contrast `||P_E X_c||²`. This measures the coordinate variance **accessible to the baseline's
zero-eigenvalue directions** — not variance visible at coarse resolution. The two are different
projections of different objects answering different questions.

> A total squared norm measures how much is present. A projected squared norm measures how much is
> present in the particular sector relevant to the question.

**Norm additivity does not imply spectral decoupling.** If `v = v_old + v_new`, then

```
|v><v| = |v_old><v_old| + |v_new><v_new| + |v_old><v_new| + |v_new><v_old| ,
```

and the cross terms are generally nonzero. A new perturbation can couple the sectors even where the
baseline splits exactly. Conversely, a nonzero cross term in one component does not by itself
establish coupling in the total operator: it may cancel against other components.

## 4. Two qualifications (c), verified before recording

**A. Equality with the full moment.** Let `V_mu = ||X_c||²`. Since `P_E` is an orthogonal
projection, `||P_E X_c|| ≤ ||X_c||` with equality **iff `X_c ∈ E`**. It is *not* necessary that `E`
be the entire mean-zero subspace: `G` may constrain other directions while leaving `X_c` untouched.
Checked: 4 nodes, `G = c|v><v|` with `v ⊥ 1` and `v ⊥ X_c`, giving `||P_E X_c||² = V_mu = 10` with
`dim E = 2 < 3 = dim(mean-zero)`. The leading response coefficient keeps its factor 2:
**`2 ||P_E X_c||²`**.

**B. Vanishing quadratic response.** `P_E X_c = 0` does **not** generally force a sixth-order
leading term. Counterexample, confirmed exactly here: counting measure on `X = (−1,0,1)ᵀ`,
`P_1 = 1 1ᵀ/3`, `P_X = X Xᵀ/2`, `G = 6 P_1 + P_X`, `(T_eps)_ij = 2cosh(eps(X_i − X_j))`. Then
`A_0 = 4J + XXᵀ/2`, so `Spec(A_0) = {0, 1, 12}`, `E = span{(1,−2,1)ᵀ}`, `P_E X = 0`, and

```
lambda_min(A_eps) = + eps⁴/6 + O(eps⁶)        (POSITIVE and QUARTIC)
```

verified in 60-digit arithmetic (`|r − 1/6|` falls by a factor 100 per decade, the expected
`O(eps²)` correction). **This refutes a general sixth-order claim.** The `a⁶` behaviour in
[`REPORT.md`](REPORT.md) case (C) remains correct **for that example**; it is not general.

### The quartic coefficient, with its hypotheses *(corrected 2026-10-07 — see the note below)*

**Hypotheses.** Finite-dimensional real space, Euclidean inner product (counting measure);
`G = G* ⪰ 0`; `A_0 = G + 2|1><1|`; `E = ker(A_0) = span(v̂)` **one-dimensional** with `‖v̂‖ = 1`;
a positive spectral gap on `E^⊥`; and `<X,v̂> = 0`, so the quadratic coefficient vanishes. Writing
`A_eps = A_0 + eps²B + eps⁴D + O(eps⁶)` with `B_ij = (X_i−X_j)²` and `D_ij = (X_i−X_j)⁴/12`, and
`A_0⁺` for the Moore–Penrose inverse, second-order Rayleigh–Schrödinger theory for the isolated
simple eigenvalue `0` gives

```
c4 = <v̂, D v̂> − <B v̂, A_0⁺ B v̂> .
```

Since `<1,v̂> = <X,v̂> = 0`, one has `B v̂ = <X²,v̂> 1` and `<v̂, D v̂> = ½|<X²,v̂>|²`, hence

```
c4 = |<X², v̂>|² · eta ,        eta := 1/2 − <1, A_0⁺ 1> .
```

The branch from `0` is the **minimum** eigenvalue for small `eps`: the other eigenvalues stay above
`γ/2` once `‖eps²B + eps⁴D + …‖ < γ/2`, while this branch is `O(eps⁴)`. *(This scalar formula is
not extended to a higher-dimensional kernel; that would require analysing the compressed operator,
and no such generalisation is claimed.)*

**The bracket criterion.** Put `z = A_0⁺1`, `s = <1,z>`, `eta = 1/2 − s`. Since `v̂` is mean-zero,
`1 ⊥ E = ker A_0`, so `1 ∈ range(A_0)` and `A_0 z = 1`. Pairing that with `z`,

```
s = <z, A_0 z> = <z,Gz> + 2s² ,   so   <z,Gz> = s(1 − 2s) ≥ 0 ,   giving   0 < s ≤ 1/2,  eta ≥ 0.
```

> **`eta = 0` if and only if `ker(G)` contains a vector of nonzero mean.**
>
> *(⇒)* If `s = 1/2` then `A_0 z = Gz + 2s·1 = 1` gives `Gz = (1−2s)·1 = 0`, so `z ∈ ker G`, and
> `<1,z> = s = 1/2 ≠ 0`.
> *(⇐)* If `w ∈ ker G` has `<1,w> ≠ 0`, pair `A_0 z = 1` with `w`:
> `<1,w> = <A_0 z, w> = <z, A_0 w> = <z, Gw + 2<1,w>1> = 2<1,w> s`, so `s = 1/2`. ∎

**`G·1 = 0` is sufficient but *not* necessary.** It is sufficient because then `1 ∈ ker G` and
`<1,1> = ‖1‖² ≠ 0`. It is not necessary: Control A in
[`PROJECTION_CONTEXT_CHECKS.md`](PROJECTION_CONTEXT_CHECKS.md) has `G = gg^T`, `g = 1 + X`, with
`G·1 = (0,3,6) ≠ 0` yet `A_0⁺1 = (5/12, 1/6, −1/12)`, `s = 1/2`, `eta = 0`, `c4 = 0` — because
`ker(G) = g^⊥` contains `(1,0,0)`, of mean 1.

**Two independent mechanisms.** `c4` is a product, so

```
c4 = 0   <=>   <X², v̂> = 0   OR   eta = 0 .
```

Control D exhibits the first with the bracket strictly positive: `X = (0,1,2,3)`,
`G = I − vv^T/20`, `v = (1,−3,3,−1)`, where `<1,v> = <X,v> = <X²,v> = 0`, `<1,A_0⁺1> = 4/9`,
`eta = 1/18 > 0`, and `c4 = 0` through the moment factor alone. Control C is the original
sufficient condition (`G = 3XX^T`, `G·1 = 0`, `eta = 0`). Control B is the positive-quartic
example (`s = 1/4`, `|<X²,v̂>|² = 2/3`, `c4 = 1/6`).

**What `c4 = 0` does and does not give.** Under these hypotheses it gives
`λ_min(A_eps) = O(eps⁶)`. It does **not** establish that sixth order is the first nonvanishing
order, nor fix its sign. For Controls A and D the sixth-order ratios were actually computed here
(70-digit, converging to `−1/3` and `−1/10`); for the original `REPORT.md` case (C) the value
`−0.8` is **carried over** from that report's own numerics and was not re-run.

## 5. Comparison, relevance and limits

| | selects | identity actually established | what does **not** follow |
|---|---|---|---|
| **Quotient-induced projection onto the inherited sector** `P_old = J_e S_e/m_e` | observables constant on reduction fibres | `[P_old, L_{e+1,pt,x}] = 0`, hence the exact orthogonal splitting `L_{e+1,pt,x} ≅ L_{e,t,x} ⊕ N_{e+1,pt,x}` with spectra and singular values a multiset union (Thm 9, Cors 10–11), for lifted twists `t → pt` | nothing about non-lifted primitive fine twists; nothing about any *perturbation* of `L`; no transfer to a different operator family |
| **Support projection onto an old spatial window** | test functions supported in a smaller window | monotonicity of the window infimum: `λ*` is non-increasing in `L` because the test space grows | **no** reducing decomposition — the restricted and complementary pieces are not invariant subspaces of the Weil form; the nested-window inheritance experiment failed here, and that is consistent |
| **Spectral projection onto a baseline kernel** `P_E` | the zero-eigenvalue directions of `A_0` | the leading response `λ_min(A_eps) = −2eps²‖P_E X_c‖² + O(eps⁴)` in fixed finite dimension | `P_E` does **not** commute with `A_eps`; no splitting, no control beyond leading order, no uniformity in a growing window |

**The failed nested-window inheritance experiment does not contradict the projected-moment
theorem**: the projections are different objects. One is a quotient-induced projection commuting
with the operator; the other a support projection with no such property; the third a baseline
spectral projection controlling one coefficient.

**Connection to Remark 10.** The note's Remark 10 separates (a) `K(L,γ)` as the exact squared norm
of derivative evaluation, (b) its attainment by `f ∝ u sin(γu)`, and (c) the fact that neither makes
the minimiser of the *entire* Weil form attain that coefficient. The projected-moment statement is
exactly the quantitative form of (c): the response is `4‖P_E(x sin γx)‖² ≤ 4K(L,γ)`, with equality
only if `x sin γx` lies in the baseline kernel. **Maximising derivative evaluation does not
necessarily minimise the whole Weil form.**

**Carried over, not rechecked:** [`REPORT.md`](REPORT.md) §4 records that the toy quartet projection
reduces `K` by only ~0.3%, against a measured `C/4K` of 0.419–0.963 — so **the toy projection does
not quantitatively explain the measured saturation gap**. That finding is attributed to the existing
report and was not re-run for this note.

**Limits, stated explicitly.**

- No quotient-induced reducing projection for the actual Weil operator is constructed here.
- No exact Weil kernel and no uniform spectral gap is established by this comparison.
- Applying a finite-matrix result to the unbounded Weil form needs a domain/form-theoretic argument
  that is not supplied.
- The inheritance theorem supplies **context**, not the missing positivity input nor the
  remaining-zero estimate of Theorem 8.
- Cubic window growth remains measure/normalization-dependent (`coeff ~ L^{s+2}` with mass `~L^s`).
- This comparison does **not** reopen the closed inheritance, scalar-coupling or quasi-RH-transfer
  investigations.

## 6. Provenance

- **Inheritance paper** — read in full from the Zenodo PDF; Definition 2, Theorem 6 (+ Remark 7),
  Lemma 8, Theorem 9, Corollaries 10–11 quoted above by their own numbers.
- **Weil manuscript** — `paper/weil_window_control.tex` at HEAD `12fde24` (the Retraction-18
  revision): Theorem 7 (relative perturbation estimate), Theorem 8 (conditional, with
  `Q_rest = Q_0 − 4F(γ)² = Q_ζ − 2F(γ)²`), Remark 10 (the (a)/(b)/(c) distinction).
- **Exploration** — [`REPORT.md`](REPORT.md) and its scripts in this directory.
- **Background that is established, not ours:** the rank-two identity
  `λ_± = a·b ± ‖a‖‖b‖` is published (Alfakih, arXiv:1609.07055, §2.2 Prop. 2.2, p. 7), and the
  spectral compression uses standard Hermitian (Rellich/Kato) degenerate perturbation theory.
  **Greenbaum–Li–Overton (arXiv:1903.00785) is first-order theory for *simple* eigenvalues of
  general matrices and does not prove the degenerate statement used here.** Kato's monograph:
  **NOT CHECKED** (unavailable here). No new literature search was performed for this note.

---

## Correction note — 2026-10-07

The first version of §4 B of this note stated that the quartic bracket "vanishes exactly when
`G·1 = 0`". **That was too strong under the standing assumption `G ⪰ 0`, and is withdrawn.**
`G·1 = 0` is sufficient; the exact condition is that **`ker(G)` contain a vector of nonzero mean**.
Control A (`G = gg^T` with `g = 1 + X`) refutes the old "only if" direction: `G·1 ≠ 0` while
`eta = 0`.

Also sharpened in the same pass: the hypotheses are now stated (finite dimension, Euclidean inner
product, one-dimensional kernel, positive gap, `<X,v̂> = 0`); bracket cancellation (`eta = 0`) is
distinguished from full-coefficient cancellation (`c4 = 0`, which also occurs when `<X²,v̂> = 0`);
and the claim that `c4 = 0` implies a *nonzero* sixth-order term has been removed — `c4 = 0` gives
only `O(eps⁶)`.

The correction was prompted by discussion and then **derived and verified independently here** in
exact rational arithmetic; the earlier text is not silently erased. Nothing else in this note
changes: the inheritance material of §1, the perturbative projection of §2, the Pythagorean/cross-
term discussion of §3, the comparison table and the scope limits of §5 all stand as written.
Status remains **INTEGRATE-AS-CONTEXT**.
