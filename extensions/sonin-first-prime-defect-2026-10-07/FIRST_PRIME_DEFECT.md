# FIRST_PRIME_DEFECT.md

**One analytic feasibility test for `S = {∞, 2}`.** Dated 2026-10-07. Labels: **C99-**,
**CC2021-**, **CCM2024-** as fixed in `SOURCES_AND_DEPENDENCIES.md`. Every claim below is tagged
**[SOURCE]** (proved in a checked source), **[STD]** (standard operator identity, verified here),
**[NEW]** (derived here), **[NUM]** (numerical evidence only) or **[OPEN]**.

> **AMENDED 2026-10-07 (later pass, after review).** §3.2's reach claim was **overreached and is
> narrowed below**, and §3.4's parity reading is **withdrawn**. The metric `J^*J` has finite Laurent
> support; the projection `Π_S` contains `K^{-1}`, and **the inverse does not**. Explicitly,
> `1/M(s) = 2 + 2√2 cos(s log2) + 2 cos(2s log2) + ⋯` already carries the `k = 2` ("4") harmonic and
> all higher ones. So **the restriction of the entire construction to `L < log 2` is NOT proved**
> and is retracted as a theorem; what survives is a statement about the bare multiplier. Parity is
> likewise not a barrier: `∂_σ log M_σ(s)|_{σ=1/2} = 2log2 Σ_{k≥1}2^{-kσ}cos(ks log2)` is **even**
> in `s` and carries exactly the weights `log2·2^{-k/2}`. Verified in `amend_checks.py` /
> `amend_output.txt`. **The verdict SOURCE/DOMAIN GAP is unchanged** — the arithmetic comparison was
> not completed either way. §§1–2, 4, 5, 6 stand.

---

## 1. Fixing the spaces, the map, the projection, and the norms

`S = {∞, 2}`; `λ > 0` the Sonin cutoff; `L > 0` the additive test half-width. `Π` always denotes a
projection, `S` always a set of places.

`H_∞ = L²(R)_ev`, `H_S = L²(X_S)^{K_S}`, each with its **own** inner product.
`M_∞ = S_λ(R,e_∞)`, `M_S = S_λ(X_S,α)` **[SOURCE: CCM2024-Def 4.5]**.
`P` = orthogonal projection of `H_∞` onto `M_∞` (CC2021's `S`, renamed to free the letter).
`Π_S` = orthogonal projection of `H_S` onto `M_S`.

`J := θ_S : H_∞ → H_S`, `θ_S(f) =` class of `σ_S ⊗ f`, `σ_2 = ε_0 − ½ε_1`. By
**[SOURCE: CCM2024-Prop 4.3(i)]** `J` is bounded with bounded inverse; by
**[SOURCE: CCM2024-Thm 4.6]** `J(M_∞) = M_S` **with equality** (the proof establishes surjectivity).
By **[SOURCE: CCM2024-Prop 4.7(iii)]** `⟨θ_S f | η_S g⟩ = ⟨f|g⟩`, i.e. `θ_S^*η_S = Id`, so
`η_S = (θ_S^*)^{-1}`. **`η_S` is never used below in place of `J^{-1}`**, and `ϑ_S(g)` — the
integrated scaling representation — is kept distinct from both.

### 1.1 The common unitary picture

`U_∞ := F_μ w_∞` and `U_S := F_μ w_S` are **unitary** onto `L²(R)`
**[SOURCE: CCM2024, proof of Prop 4.7(iii)]**, with `F_μ(g)(s) = ∫ g(λ)λ^{-is}d^*λ`. Transporting
both sides by these unitaries preserves the **original** inner products. Write `P` for
`U_∞ P U_∞^{-1}`, etc. Then by **[SOURCE: CCM2024-Prop 4.6(ii), eq. 57]** with `S∖{∞} = {2}`:

```
J  =  multiplication by  m(s) := 1 - 2^(-1/2 - i s)  =  L_2(1/2 + i s)^(-1)
```

and from the **proof** of that proposition, `w_S(σ_S⊗f)(1×λ) = g(λ) − 2^{-1/2}g(λ/2)` with
`g = w_∞(f)`, together with `F_μ(g(·/p))(s) = p^{-is}F_μ(g)(s)`. So, as an operator identity,

```
J = Id - 2^(-1/2) D^(-1),    D := multiplication by 2^(i s)  =  the scaling  λ ↦ 2λ .   [SOURCE]
```

**The metric.** `H := J^*J =` multiplication by

```
M(s) := |m(s)|^2 = |L_2(1/2+is)|^(-2) = 3/2 - sqrt(2) cos(s log 2)
                 = 3/2 - 2^(-1/2) ( D + D^(-1) )  as an operator.          [SOURCE + STD]
```

**[NUM]** verified to `6.7e-16` over 400 nodes, and the operator form exactly (`checks.py` §3–4).
`min M = 3/2−√2 = 8.578644e-02`, `max M = 3/2+√2 = 2.914214e+00`, `mean M = 3/2` (the cosine
averages out). `M(s) > 0` everywhere, since `m(s) = 0` would need `2^{-1/2} = 1`; the minimum is
attained only on the discrete set `s ∈ (2π/log2)Z`, spacing `9.064720`.

### 1.2 The actual orthogonal projection, and the two objects it is not

**[STD] Lemma 1.** Let `J` be bounded with bounded inverse, `P` an orthogonal projection with range
`M`, `H = J^*J`, and `K := (PHP)|_M`. Then `K ≥ ‖J^{-1}‖^{-2} Id_M > 0`, so `K` is boundedly
invertible on `M` (closed range, coercive), and with `K^{-1}` extended by zero on `M^⊥`,

```
Π := J P K^(-1) P J^*
```

is the **orthogonal** projection onto `J(M)`.

*Proof.* `Π^* = Π` since `K = K^* > 0`. `Π² = JPK^{-1}(PJ^*JP)K^{-1}PJ^* = JPK^{-1}KK^{-1}PJ^* = Π`.
For `v ∈ M`, `ΠJv = JPK^{-1}PJ^*Jv = JPK^{-1}Kv = Jv`. And `ran Π ⊆ J(M)`. ∎

Here `K ≥ (3/2−√2)Id_{M_∞}` because `M(s) ≥ 3/2−√2`, so the coercivity hypothesis **holds for
`S = {∞,2}`** and is not assumed.

**[NUM] Control 2 of the brief**, exact over `Q` (`checks.py` §1), with a `J` that is neither
unitary nor commuting with `P`: `Π` is self-adjoint ✓, idempotent ✓, fixes `J(M)` pointwise ✓,
annihilates `J(M)^⊥` ✓, `tr Π = dim M = 2` ✓. By contrast `JPJ^{-1}` is idempotent but **not**
self-adjoint, and `JPJ^*` is self-adjoint but **not** idempotent (`tr = 15 ≠ 2`); neither equals `Π`.
Both safeguards of the brief hold: `JBJ^* ≥ 0` for `B ≥ 0`, and `tr(JPJ^{-1}) = tr P = 2`
(similarity invariance — non-unitarity does **not** erase it). **Control 5:** for a singular `J`,
`det K = 0` and the formula is inapplicable, so the hypothesis is doing work.

`P` is **not** assumed to commute with `H`, with `D`, or with scaling. `J`, `H`, `D` and
`F := U_S ϑ_S(f)U_S^{-1} =` multiplication by `\hat f(s)` all commute with **each other** (they are
multiplication operators in `s`) but with none of them does `P` commute.

---

## 2. The exact reduction  **[NEW]**

**Proposition 2.** With the above notation, and subject to the trace-class hypothesis (T) below,

```
Tr( vartheta_S(f) Pi_S )  =  Tr_{M_inf} ( K^(-1) · P ( M F ) P ) ,
      K = ( P M P )|_{M_inf} = (3/2) P  -  sqrt(2) · P C P ,
      C := ( D + D^(-1) ) / 2 = ( vartheta(2) + vartheta(2)^(-1) ) / 2 ,   ||C|| = 1 .
```

*Proof.* `Tr(F·JPK^{-1}PJ^*) = Tr(PJ^*FJP·K^{-1})` by cyclicity, and `J^*FJ =` multiplication by
`\bar m \hat f m = M\hat f`, since multiplication operators commute. ∎

**(T) Trace-class hypothesis, stated and not hidden.** The manipulation needs
`ϑ_S(f)Π_S` trace class and the cyclic exchange legitimate. At the archimedean place this is
**[SOURCE: CC2021-Thm 4.6]**, which asserts `Tr(ϑ(f)S)` for all `f ∈ C_c^∞(R_+^*)` with no support
restriction. **No checked source establishes (T) semi-locally.** It is **[OPEN]** and is carried as
a hypothesis, not as a fact.

**The companion positivity.** For any orthogonal projection, `Π_S = Π_S^2 = Π_S^*` gives

```
Tr( vartheta_S(g) Pi_S vartheta_S(g)^* ) = || vartheta_S(g) Pi_S ||_HS^2  >= 0 ,
   and  = Tr( vartheta_S(g^* * g) Pi_S )   by cyclicity and idempotence.          [STD]
```

So `P_S(g) := Tr(ϑ_S(g)Π_Sϑ_S(g)^*) ≥ 0` **free of charge**, needing no `ϑ_S`-invariance of `M_S` —
exactly the archimedean mechanism. Finiteness is again (T). **[D1 of TARGET.md: conditional on (T).]**

### 2.1 Unconditional Neumann expansion of `K^{-1}`  **[NEW]**

`K = (3/2)[P − (2√2/3)·PCP]` and `‖PCP‖ ≤ ‖C‖ = 1`, while `2√2/3 = 0.9428090 < 1`. Hence

```
K^(-1) = (2/3) * sum_{n>=0} ( (2 sqrt2 / 3) P C P )^n           norm-convergent, UNCONDITIONALLY.
```

This is a real structural gain and it is `S`-free in the following sense: convergence needs only the
universal bound `‖PCP‖ ≤ 1`, no spectral information about the Sonin space. The margin is thin —
the worst-case tail after `n` terms is `r^{n+1}/(1−r)` with `r = 0.9428090`, i.e. `16.5` at `n=0`
and `0.380` only by `n = 64` (`reach.py` §9). **Control:** the expansion is sharp, not slack — for a
multiplier `|1 − a e^{-it}|²` the ratio is `2a/(1+a²)`, which is `0.942809` at `a = 2^{-1/2}`,
`0.994475` at `a = 0.9`, and exactly `1` at `a = 1`, where the multiplier acquires a real zero and
the series **diverges**. The convergence is therefore owed specifically to `a = p^{-1/2} < 1`, i.e.
to the critical line sitting strictly inside the half-plane of absolute convergence.

### 2.2 The baseline, and where any gain must come from

Ambient worst case, from `m ≤ M ≤ Mx` alone:

```
(m/Mx) P_inf(g)  <=  P_S(g)  <=  (Mx/m) P_inf(g) ,     Mx/m = 33.97056 ,
||K^(-1)|| <= 1/(3/2 - sqrt2) = 11.65685 .
```

Mean-field value `‖K^{-1}‖ ≈ 2/3 = 0.6666667`. **Gap factor `17.48528`.** By Proposition 2 the entire
gap is controlled by one quantity:

```
K >= ( 3/2 - sqrt(2) ||P C P|| ) P ,       C = ( vartheta(2) + vartheta(2)^(-1) ) / 2 .
```

**[NOT ESTABLISHED IN THIS AUDIT]** Is `‖PCP‖ < 1`, quantitatively? *(Relabelled 2026-10-07:
this is not "genuinely open" until its precise formulation, literature status and possible
approximate invariant vectors have been checked, none of which this pass did. Lack of exact
invariance does not by itself give a uniform norm gap. Note also that **invertibility of `K` needs
only `‖PCP‖ ≤ 1`**, since `2√2/3 < 1`; a strict inequality would only sharpen the constant
`11.65685`, not enable the expansion.)* Equivalently: can the Mellin image of the Sonin space
`S_λ(R,e_∞)` concentrate its `L²` mass near the discrete set `s ∈ (2π/log 2)Z` where `M` attains its
minimum? Note `ϑ(2)S_λ(R,e_∞) ⊄ S_λ(R,e_∞)` and `⊅` — scaling by 2 maps `S_λ` into neither
`S_λ` nor its complement (it sends the vanishing conditions to `|x| < 2λ` and `|x| < λ/2`) — so
`PCP` is a genuine compression; but a compression of a norm-one multiplication operator can have
norm one with no eigenvector. **Not settled by any checked source**, and the finite-dimensional
samples in `checks.py` §2 do **not** certify it: finite-dimensional positivity does not certify
continuum positivity, and no quadrature, basis-truncation, tail or inverse-error budget for `B_λ`
was established. At `‖PCP‖ = 1` the targeted bound degenerates exactly to the ambient one
(`checks.py` §7).

---

## 3. What the arithmetic side actually is  **[NEW — and this is where the route fails]**

Proposition 2's right-hand side is governed by `MF`. Since `M = (3/2) − 2^{-1/2}(D + D^{-1})` and
`D` is the scaling `λ ↦ 2λ`, `MF` is multiplication by the Mellin transform of

```
(3/2) f  -  2^(-1/2) [ f(·/2) + f(2·) ]  —  exactly three dilates, at  2^0, 2^(+1), 2^(-1).
```

### 3.1 The one case where the identity closes

If — and only if — `PCP = 0`, so that `K = (3/2)P`, Proposition 2 becomes
`Tr(ϑ_S(f)Π_S) = (2/3)Tr(P(MF)P)`, and each of the three pieces is an **archimedean** Sonin trace,
to which **[SOURCE: CC2021-Thm 4.6]** applies *with no support restriction*:

```
Tr( vartheta_S(f) Pi_S ) = W_inf(f) - (sqrt2/3)[ W_inf(f(·/2)) + W_inf(f(2·)) ]
                         +   E(f)   - (sqrt2/3)[   E(f(·/2))   +   E(f(2·))   ] ,
```

with `E(h) = ∫ h(ρ^{-1})ε(ρ)d^*ρ` from CC2021 eq. (85) — **independently defined, not by
subtraction**. `PCP = 0` says `Pϑ(2)P = −(Pϑ(2)P)^*`, which is **not** established anywhere and is
false for a generic pair; it is recorded as the mean-field idealisation, not as a hypothesis met.

**The arithmetic side of this identity is purely archimedean.** There is no `W_2` term. That is not
an artefact of `PCP = 0`: restoring `K^{-1}` by §2.1 replaces the three clean traces by interleaved
products `Tr(PD^{k_1}P D^{k_2}P ⋯ FP)`, which **no checked source evaluates at all** — so the only
regime in which the reduction reaches an arithmetic statement is the one that reaches an
archimedean one.

### 3.2 Reach of the **bare metric multiplier** — and what this does *not* prove  **[NARROWED]**

**Proposition 3 (reach of `J^*J`, as proved).** For any finite `S ∋ ∞`, `H = θ_S^*θ_S` is
multiplication by `∏_{p∈S∖{∞}} |L_p(½+is)|^{-2} = ∏_p (1 − p^{−½−is})(1 − p^{−½+is})`, a Laurent
polynomial of degree **at most 1 in each `p^{is}`**. Hence the dilates of `f` that `H F` can
generate are exactly the rationals `d/d′` with `d, d′` **coprime squarefree** divisors of
`∏_{p∈S∖{∞}} p`; the dilate `p^k` with `k ≥ 2` does not occur **in `H F`**.

*Proof.* Each factor contributes only the monomials `p^{0,±is}`; a product of degree-`≤1` Laurent
polynomials in independent variables has exponents in `{−1,0,1}` per variable. ∎

**[NUM]** Verified by exact rational expansion for `S∖{∞} = {2}, {2,3}, {2,3,5}`: 3, 9 and 27
monomials, every exponent in `{−1,0,1}`, dilates `{1/2,1,2}`, `{1/6,…,6}`, `{1/30,…,30}` — no
non-squarefree numerator or denominator (`reach.py` §8a). **This part is correct and stands.**

> **RETRACTED 2026-10-07 — the inference to the full construction.** An earlier version of this
> section concluded a *structural ceiling* `L < log 2 = 0.6931472` for the whole construction, on
> the ground that the first non-squarefree prime power `n = 4` switches on at `2L > log 4`. **That
> inference is invalid.** `Π_S = J P K^{-1} P J^*` contains the **inverse** `K^{-1}`, and finite
> Laurent support of `J^*J` says nothing about the support of `K^{-1}`. The uncompressed scalar case
> makes it explicit: with `r = 2^{-1/2}`, `θ = s log 2`, the Poisson-kernel identity
> `(1 − 2r\cosθ + r²)^{-1} = (1−r²)^{-1}(1 + 2Σ_{k≥1} r^k \cos kθ)` gives
>
>   `1/M(s) = 2 + 2√2 cos(s log 2) + 2 cos(2 s log 2) + √2 cos(3 s log 2) + ⋯`
>
> **The `k = 2` harmonic — the "4" shift — is present with coefficient exactly `2`, and every
> higher one is present too.** **[NUM]** `amend_checks.py` §A: the identity holds to `9.0e-14` over
> 400 nodes and the coefficients are `2, 2.8284271, 2, 1.4142136, 1, 0.7071068` for `k = 0,…,5`.
>
> With the projection inserted the relevant objects are the **compositions** `PCPCP, …`, not free
> dilations, so the finite model is only indicative — but it is enough to kill the support argument.
> **[NUM]** `amend_checks.py` §C: in a finite one-step shift model with compression, `(PCP)` has
> entry `(0,2) = 0` while `(PCP)²` has entry `(0,2) = 0.25`, i.e. **two-step reach appears at
> `n = 2`**; the control `C = Id` generates no shift at all. That model is **not** the Sonin space
> and shows nothing about coefficients.
>
> **What is therefore true and what is not.** `Π_S` is **not** shown to produce the correct Weil
> prime-power coefficients — cancellations in the full assembly are not ruled out in either
> direction. But **no general prime-power obstruction is proved either**, and the claim that the
> construction is confined to `L < log 2` is withdrawn. The correct replacement:
>
> > **The bare metric multiplier contains only squarefree-ratio shifts. Whether the full
> > compressed-inverse projection supplies the required prime-power contributions remains
> > unresolved by that support calculation.**
>
> **Consequence for the project, restated honestly.** Nothing above excludes Zhu's window `L ≤ 0.8`
> on reach grounds. The reason the route does not deliver there is the one in §5 and §4: the
> arithmetic comparison is not completed, not that higher harmonics cannot appear.

### 3.3 On the calibration window the failure is *not* one of support

Honest refinement, because it cuts against the easy version of the argument. On
`(log 2)/2 < L < (log 3)/2` the autocorrelation `h = g^* * g` has `supp h ⊂ [−2L, 2L]` with
`2L < log 3 = 1.0986123 < log 4 = 1.3862944`, so the prime-2 block carries **only** the `k = 1`
atoms, of weight `(log 2)·2^{-1/2} = 0.4901291`:

```
W_2(h) = (log 2) 2^(-1/2) [ h(2) + h(1/2) ]        on this window, exactly.      [NUM/STD]
```

These are precisely the shifts `2^{±1}` that `M` supplies. **So on the calibration window the shift
supports coincide**, and Proposition 3 does **not** exclude the route there. What would then be
needed is an equality of *weights and signs* between the atoms the identity of §3.1 produces at
`h(2^{±1})` — which come from propagating `W_∞`'s own atom at `u = 1` through the three shifts — and
`−(log 2)2^{-1/2}`. **That comparison was not completed**: it requires the exact `h(1)`-coefficient
of `W_∞` in CC2021's principal-value normalisation (C99-Thm VII.4 fixes `∫′` by the condition that
the `α_v`-Fourier transform vanish at 1), which this pass did not pin down. Stating a numerical
mismatch here without that constant would be exactly the invented-constant error this repository has
already made once. **Recorded as [OPEN], not as a result.** The route is therefore *bounded above*
by `L < log 2` and *undecided* on `(log2)/2 < L < (log3)/2`.

### 3.4 Where the prime-2 information sits  **[INTERPRETATION AMENDED]**

The `W_2` data is not absent from the metric. Two expansions, both elementary and both verified:

```
d/ds     log M(s)               =  2 sum_{k>=1} (log 2) 2^(-k/2) sin( k s log 2 )      [ODD in s]
d/dsigma log M_sigma(s)|_(1/2)  =  2 sum_{k>=1} (log 2) 2^(-k/2) cos( k s log 2 )      [EVEN in s]
```

with `M_σ(s) = |1 − 2^{−σ−is}|²`, so that `log M_σ(s) = −2Σ_{k≥1} k^{-1}2^{-kσ}cos(ks log2)`.
**[NUM]** the first to `5.4e-09` (`checks.py` §6), the second to `3.3e-10` with weights
`0.4901291, 0.3465736, 0.2450645, 0.1732868` (`amend_checks.py` §B).

> **WITHDRAWN 2026-10-07.** An earlier version of this section offered "`M` is even in `s` while the
> `W_2` generator is odd" as the obstruction. **That is not a barrier.** It is a property of the
> `s`-derivative alone; the `σ`-derivative of the same local factor is **even** in `s` and carries
> exactly the even prime-power weights. This is nothing more than the expansion of the local Euler
> factor — **not** a new arithmetic theorem, and **not** a construction of a corresponding family of
> geometric operators: CCM2024's `θ_S` is defined in its own critical-line setting, and introducing
> `M_σ` algebraically does not produce geometric operators off `σ = ½`.

**The obstruction, as it should be stated.** The missing ingredient is **not** the presence of
prime-power information in the local factor — it is plainly there, in both parities. It is an
**independently justified operation that extracts that information inside the required trace
formula, with the correct coefficients and a controlled remainder.** That is what `Tr(K^{-1}P·MF·P)`
has not been shown to do, and what no checked source supplies.

## 4. Route A: the only characteristic-zero arithmetic identity, and its defect  **[NEW assembly]**

Because §3 blocks the `θ_S` route from the arithmetic, the arithmetic must come from C99. For
`k = Q`, `S = {∞,2}`, **[SOURCE: C99-Thm VII.4]** gives, with `R_Λ := P̂_Λ P_Λ`:

```
Tr( R_Lambda U(h) ) = 2 h(1) log'(Lambda) + W_inf(h) + W_2(h) + o(1) ,   Lambda -> infinity,
```

`W_v(h) = ∫′_{k_v^*} h(u^{-1})/|1−u| d^*u`. **This is the semi-local Weil sum, in characteristic
zero, proved.** But `R_Λ` is a product of two orthogonal projections: **not self-adjoint, not
positive**, and per the source gate no sign is assigned to it.

**[STD] Lemma 4 (two-projection defect).** For orthogonal projections `P_Λ, P̂_Λ`, Halmos' five-part
decomposition gives `ran∩ran`, `ran∩ker`, `ker∩ran`, `ker∩ker` and a generic part on which the pair
is unitarily equivalent to `[[1,0],[0,0]]` and `[[c²,cs],[cs,s²]]`, `0 < c < 1`. On the four trivial
parts `P̂_ΛP_Λ` equals `Q_Λ` (the orthogonal projection onto `ran P_Λ ∩ ran P̂_Λ`); on the generic
part `P̂_ΛP_Λ` is `[[c²,0],[cs,0]]` and `Q_Λ` is `0`. Hence with `G_Λ := P̂_ΛP_Λ − Q_Λ`:

```
Tr(G_Lambda) = sum_j c_j^2  >= 0 ,      || G_Lambda ||_1 = sum_j c_j
```

(each generic block is rank one with `‖·‖_1 = √(c⁴+c²s²) = c`). **[NUM]** `Tr(P̂P) − Tr Q = Σc_j²`
confirmed on four random projection pairs in dimensions 6 and 8, with a commuting control giving an
empty generic part and `Tr(P̂P) = dim(ran ∩ ran)` exactly (`checks.py` §2; the trace-norm formula
is the elementary deduction from the 2×2 block form, **not** separately verified numerically).

**Consequence [NEW].** `Q_Λ` *is* an orthogonal projection, so `Tr(U(g)Q_ΛU(g)^*) = ‖U(g)Q_Λ‖²_{HS} ≥ 0`
and equals `Tr(Q_ΛU(g^**g))`. Combining with C99-Thm VII.4, for `h = g^* * g`:

```
0 <= || U(g) Q_Lambda ||_HS^2 = 2 h(1) log'(Lambda) + W_inf(h) + W_2(h) - E_Lambda(h) + o(1),
     E_Lambda(h) := Tr( G_Lambda U(h) ) ,    G_Lambda = Phat_Lambda P_Lambda - Q_Lambda .
```

`E_Λ` is **defined by two-projection geometry alone**, with no reference to `W` — it satisfies the
brief's independence requirement, unlike a difference. **But the ACCOUNTING GATE fails on
evaluability:**

1. **`E_Λ` is not independently evaluable.** It needs the principal angles between `P_Λ` and `P̂_Λ`
   on `L²(X_{\{∞,2\}})`. **No checked source computes or bounds them**, and nothing here does.
2. **The divergence is uncancelled.** `h(1) = ∫|g|²d^*x > 0` for `g ≠ 0`, so `2h(1)log′Λ → +∞` and
   the inequality is vacuous unless `E_Λ(h)` is shown to absorb it. That is the whole difficulty,
   not a technicality.
3. **Trace class.** `R_ΛU(h)` is trace class by C99; **`Q_ΛU(h)` is not known to be**.
4. **Wrong cutoff corner.** C99's `B_Λ = {f : f(x) = 0 & \hat f(x) = 0 for |x| > Λ}` is doubly
   **low**-pass — the `ran ∩ ran` corner. The Sonin space is doubly **high**-pass — the `ker ∩ ker`
   corner. These are *different* subspaces of the same two-projection system. Expressing
   `Π_{Sonin} = Id − Q_{ran∩ran} − Q_{ran∩ker} − Q_{ker∩ran} − Q_{gen}` is co-finite relative to the
   identity, so `Tr(Π_{Sonin}U(h))` cannot be read off C99-Thm VII.4 term by term. **No checked
   source supplies the translation.**
5. **And in characteristic zero `Q_Λ` has no proved trace formula at all.** The upgrade from `R_Λ`
   to `Q_Λ` is **[SOURCE: C99-Cor VIII.2]**, stated "in the case of positive characteristic", and its
   engine **[SOURCE: C99-Lem VIII.1]** is a function-field count: `Λ = q^N`, `mod(C_S) = q^Z`,
   `2log′Λ = (2N+1)log q`, `∫′_{R_v^*}χ_v(u)/|1−u|d^*u = −f_v log q_v`, `|d| = q^{2−2g}` with `g`
   the genus, and the exact tally `(2N+1)l − fl + (2−2g)l`; it leans on "since `ξ` is locally
   constant, its Fourier transform has compact support", which fails at `v = ∞`. **`[P̂_Λ, P_Λ] = 0`
   is not proved for `k = Q`** in any checked source. The global formula itself is explicitly
   unproved — "**are not able to prove (16) directly for arbitrary `h`**" — and in positive
   characteristic **[SOURCE: C99-Thm VIII.5]** it is *equivalent* to RH, so it is not an input.

---

## 5. Verdict on the Target Lemma, dependency by dependency

| dep | needed | status |
|---|---|---|
| D1 | `P_S(g)` an ordinary finite positive trace | **CONDITIONAL on (T)**: positivity is free `[STD]`, finiteness is `[OPEN]` semi-locally |
| D2 | `Π_S` the actual orthogonal projection | **ESTABLISHED** `[STD]` + coercivity `K ≥ (3/2−√2)Id` from `[SOURCE]` |
| D3 | `Tr(ϑ_S(f)Π_S) = W_S(f) + E_S(f)`, `E_S` independent | **NOT ESTABLISHED.** Route B closes on archimedean `W_∞` only in the mean-field case `PCP = 0`; with `K^{-1}` restored the terms are interleaved products `Tr(PD^{k_1}PD^{k_2}P⋯FP)` that no checked source evaluates. *(Amended 2026-10-07: these **can** carry higher prime-power shifts — see §3.2 — so this is an unevaluated expression, **not** a proved obstruction.)* Route A has an independently *defined* `E_Λ` that is not independently *evaluable*, in the complementary cutoff corner, and in char 0 only for the non-self-adjoint `R_Λ` |
| D4 | a bound on `E_S` | **NOT ESTABLISHED**; not found in the sources checked |
| D5 | `n = 2` the only active prime power, not vacuous | **ESTABLISHED** `[NUM]`: on `(log2)/2 < L < (log3)/2` the active set is exactly `{2}`; at `L = (log2)/2` it is **empty** and the test is vacuous, so the window must be open; Zhu's `L ≤ 0.8` gives `{2,3,4}` |

**The Target Lemma is not obtained, and `E_S ≤ η P_S` was never reached, because `E_S` could not be
formed.** Per the stop rule, no substitute construction was tried after the failure.

## 6. Comparison against the published estimate

CC2021-Thm 6.11 delivers `W_∞(g*g^*) ≥ Tr(ϑ(g)Sϑ(g)^*) − c|ĝ(0)|²` with `c = 4γ/log2 = 16.98658…`
on `supp g ⊂ [2^{-1/2},2^{1/2}]`, `ĝ(−i/2) = 0`. **Nothing above improves or resolves it.** `Π_S`,
`K^{-1}`, the Schur-complement form and the Neumann series are **tools** by the brief's own
standard; the one quantity that would convert them into a sharper constant, `‖PCP‖`, is `[OPEN]`,
and even settling it would bound `P_S` against `P_∞` — not against `W_∞ + W_2`.

No step normalises by an unproved gap of the Weil form or by a nonexistent bounded `L²` operator
norm of it (Retraction 18 safeguard). No step uses a finite-matrix projected-moment or quartic
result; those remain diagnostic motivation only, kept out of the derivation as the brief requires.
The diverging ambient condition numbers across growing `S` are **not** used to rule anything out:
Proposition 3 is a fixed-`S` statement about the bare multiplier, true for every finite `S`
separately, and no many-prime limit is taken. Nor is a restricted fixed-`S` estimate excluded by it.
