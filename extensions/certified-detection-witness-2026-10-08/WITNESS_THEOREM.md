# The nodal projected witness: statement, proof, and what it costs

Labels: **[C]** classical / elementary, **[S]** asserted by the manuscript or notes with its
location, **[D]** derived and checked here, **[H]** hypothesis.

## 1. Setting

`L > 0`, `γ ∈ ℝ`. Admissible `f`: real, even, `supp f ⊆ [−L,L]`, `‖f‖₂ = 1`
(**[S]** paper §4.2). `F(z) = ∫_{−L}^{L} f(u)e^{izu}du`, so for real even `f`

```
F(γ)   =  <f, c>,     c(u) = cos(γu)
F'(γ)  = -<f, r>,     r(u) = u sin(γu)
F''(γ) = -<f, d>,     d(u) = u^2 cos(γu)
```

all three of `c, r, d` real and even **[D]**, inner product the real `L²([−L,L])` one. Let
`E ⊆ L²([−L,L])` be a finite-dimensional subspace of real even functions and `P_E` the
`L²`-orthogonal projector. Put `a = P_E c`, `b = P_E r`.

## 2. The construction and the lemma

> **Lemma 1 (nodal projected witness).** **[C — elementary linear algebra, no novelty claimed]**
> Suppose `a ≠ 0` and set
> ```
> w = b - (<b,a>/||a||^2) a .
> ```
> If `w = 0` the construction supplies no witness. Otherwise let `f_node = w/‖w‖₂`. Then
> `f_node ∈ E`, `‖f_node‖₂ = 1`, and
> ```
> (i)   F_node(gamma) = 0 ;
> (ii)  |F_node'(gamma)|^2 = ||w||_2^2 =: K_eff ;
> (iii) K_eff = ||P_E r||^2 - |<P_E r, P_E c>|^2 / ||P_E c||^2 ;
> (iv)  K_eff = max { |F'(gamma)|^2 : f in E, ||f||_2 = 1, F(gamma) = 0 } .
> ```
> If `a = 0`, take `w = b`; then (i) holds trivially and (ii), (iv) hold with
> `K_eff = ‖P_E r‖²`.

**Proof.** **[D]** For `f ∈ E`, `⟨f,c⟩ = ⟨f,P_E c⟩ = ⟨f,a⟩` and `⟨f,r⟩ = ⟨f,b⟩`, since `f ⊥ (1−P_E)`.

(i) `⟨w,a⟩ = ⟨b,a⟩ − (⟨b,a⟩/‖a‖²)‖a‖² = 0`, so `F_node(γ) = ⟨f_node,a⟩ = 0`.

(iii) `‖w‖² = ‖b‖² − 2(⟨b,a⟩²/‖a‖²) + (⟨b,a⟩²/‖a‖⁴)‖a‖² = ‖b‖² − ⟨b,a⟩²/‖a‖²`.

(ii) `⟨w,b⟩ = ‖b‖² − ⟨b,a⟩²/‖a‖² = ‖w‖²`, so `F_node'(γ) = −⟨f_node,b⟩ = −‖w‖²/‖w‖ = −‖w‖`,
whence `|F_node'(γ)|² = ‖w‖² = K_eff`.

(iv) The constraint set is `{f ∈ E : ‖f‖ = 1, ⟨f,a⟩ = 0}`, i.e. the unit sphere of
`E ∩ a^⊥`. On it `|F'(γ)|² = ⟨f,b⟩² = ⟨f, P_{a^⊥}b⟩²`, maximised by Cauchy–Schwarz at
`f ∝ P_{a^⊥}b = w`, with value `‖w‖²`. ∎

This is the ordinary "project out one linear constraint" computation; it is recorded as a lemma
only because every term of the detection budget must refer to the same `f`. **No general theorem
is claimed.**

> **Corollary 2 (the three expected simplifications).** **[D]** On `f_node`:
> ```
> F(gamma) F''(gamma) = 0 ,     Q_0(f_node) = Q_zeta(f_node) ,     a_f = 4 K_eff ,
> ```
> the second because `Q_0 = Q_ζ + 2F(γ)²` (**[S]** paper eq. (3)) and the third because
> `a_f = 4[F'(γ)² + F(γ)F''(γ)]` (**[S]** paper eq. (4)).

All three are **checked numerically** in `certify_witness.py`; `F(γ)` comes out at the rounding
level of the coefficient vector, and the residual is charged, not assumed zero — see §5.

## 3. The witness inequality, in the correct direction

For one fixed admissible `f` put `b_f = Q_0(f)`, `a_f = 4[F'(γ)² + F(γ)F''(γ)]`. The
**certificate must be an upper bound**:

```
Q_delta(f)  <=  b_upper - a_lower delta^2 + E_upper(delta) ,
        b_f <= b_upper ,   a_f >= a_lower > 0 .
```

If the right-hand side is `< 0` for some `δ`, then `Q_δ(f) < 0` and `f` is a negative witness
**for that planted form**. Theorem 7 of the manuscript runs the other way (`Q_δ ≥ Q_0 − B_L(δ)`),
and a *lower* sensitivity bound going negative proves nothing. Nothing below combines an infimum
attained by one vector with a curvature attained by another.

## 4. Two error routes, and which one the verdict uses

**Route T (Taylor).** `E_upper(δ) = (16/3)L⁵δ⁴e^{2Lδ}` **[S]** `proofs.md` §B.3, coefficient
re-derived in `TARGET_AND_SOURCES.md` §3. Valid for every `δ ≥ 0`; informative only while the
`δ⁴` term stays below `a_fδ²`.

**Route X (exact, and what is used).** The finite-displacement quartet is exact, with no Taylor
error **[D]**, verified in the previous extension on three test functions:

```
Q_delta(x) = Q_zeta(x) - 2<x,c>^2 + 4<x,a_delta>^2 - 4<x,b_delta>^2 ,
a_delta(u) = cos(gamma u) cosh(delta u) ,   b_delta(u) = sin(gamma u) sinh(delta u) ,
<x,a_delta> = Re sum_k x_k F_k(gamma + i delta) ,
<x,b_delta> = -Im sum_k x_k F_k(gamma + i delta) .
```

The term `+4⟨x,a_δ⟩²` is **positive** and is retained; dropping it would manufacture negativity.
At `δ = 0`: `a_0 = c`, `b_0 = 0`, so `Q_0(x) = Q_ζ(x) + 2⟨x,c⟩²` ✔, and the response is even in
`δ` ✔ — both checked as controls.

Route X removes the remainder from the budget entirely. The only error left is the enclosure of
`Q_ζ(x)` itself, i.e. the archimedean tail bound of `TARGET_AND_SOURCES.md` §4, plus double
rounding. **So the verdict is not an artefact of a loose Taylor allowance.**

## 5. Residual-node charge, not an assumption

`f_node` is stored as a float coefficient vector, so `F_node(γ) = ⟨x,a⟩` is zero only to rounding.
Rather than treating it as zero, the residual is charged exactly: `Q_0 = Q_ζ + 2⟨x,c⟩²` and
`a_f = 4[F'² + F(−F'')]` are **both evaluated with the actual `⟨x,c⟩`**, and Route X uses
`−2⟨x,c⟩²` as computed. No term is set to zero by hand.

## 6. What the construction costs — the three comparisons

```
K        = ||r||^2        the paper's extremal value, K(L,gamma)        [S] Theorem 8 eq. (5)
||P_E r||^2               what survives projection into E               [D]
K_eff    = ||P_E r||^2 - <P_E r,P_E c>^2/||P_E c||^2                    [D] Lemma 1(iii)
```

`K_eff ≤ ‖P_E r‖² ≤ K` always **[C]**. The two ratios `K_eff/‖P_E r‖²` (price of the node) and
`‖P_E r‖²/K` (price of the finite basis) are reported separately in `results.json`, because they
are different costs and only the first is attributable to the nodal constraint.

**Scope reminders carried forward [S].** `K` is the exact squared norm of the
derivative-evaluation functional, **not** automatically the sharp coefficient of the re-optimised
minimum (paper Remark 10(a)–(c)). `K/L³ → 1/3` only as `|γ|L → ∞`; it is **not** bounded below by
`1/3` for every `γ`. Nothing here claims the measured saturation `C/(4L³/3) → 1` of §4.3, which
is T2.
