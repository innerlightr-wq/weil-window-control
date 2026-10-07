# TARGET.md

**Written before any numerical work.** Dated 2026-10-07.

## Symbols (fixed once; `S` is a set of places, `Π` is a projection)

`S = {∞, 2}`, a finite set of places of `Q` containing `∞`. `λ > 0` the **Sonin cutoff**.
`L > 0` the **half-width of the additive test support**: `g ∈ C_c^∞(R_+^*)` with
`supp g ⊂ [e^{-L}, e^{L}]`, so in the additive coordinate `u = log x`, `supp g ⊂ [−L, L]` and
`supp(g^* * g) ⊂ [−2L, 2L]`. `λ` and `L` are different quantities and are never identified.

`H_∞ = L²(R)_ev` and `H_S = L²(X_S)^{K_S}`, each with its **own original inner product**.
`M_∞ = S_λ(R, e_∞) ⊂ H_∞` and `M_S = S_λ(X_S, α) ⊂ H_S` the Sonin subspaces (CCM2024-Def 4.5).
`P` = orthogonal projection of `H_∞` onto `M_∞` (this is CC2021's `S`, renamed).
`Π_S` = **orthogonal** projection of `H_S` onto `M_S`.
`J = θ_S : H_∞ → H_S`, bounded with bounded inverse (CCM2024-Prop 4.3(i)), **not** isometric;
`J(M_∞) = M_S` by CCM2024-Thm 4.6. `η_S = (J^*)^{-1}` by CCM2024-Prop 4.7(iii) — never substituted
for `J^{-1}`. `ϑ_S(g)` is the integrated scaling representation on `H_S`, distinct from `J` and from
`η_S`. `h := g^* * g` with CC2021's involution `g^*(x) = \overline{g(x^{-1})}`.

Standard projection formula to be verified, not assumed: with `H = J^*J` and
`K = (P H P)|_{M_∞}`, `Π_S = J P K^{-1} P J^*`, where `K^{-1}` acts on `M_∞` and is extended by
zero on `M_∞^⊥`. Its purpose is to retain the metric correction that `JPJ^{-1}` and `JPJ^*` miss.
Both are *tools*, not the result.

## The statement to be attempted

> **Target Lemma (one-prime Sonin discrepancy).** Let `S = {∞,2}` and `λ > 0`. Put
> `P_S(g) := Tr(ϑ_S(g) Π_S ϑ_S(g)^*)`. Then there is an **independently defined** functional `E_S`
> and a constant `η ≤ 1` such that, for every `g` in an explicit nontrivial admissible class
> `A(L,λ)` with `(log 2)/2 < L < (log 3)/2` — so that `n = 2` is the **only** active positive prime
> power and the prime block is not identically zero —
>
>   `E_S(g) ≤ η · P_S(g)`,
>
> and such that `P_S(g) − E_S(g)` equals the repository's Weil form `Q_ζ` restricted to `A(L,λ)`,
> up to terms displayed explicitly. Where `P_S(g) = 0` the claim is read as `E_S(g) ≤ 0`, not as a
> quotient.

Sign convention fixed once: `Q_S(g) := P_S(g) − E_S(g)`, with pole / contact / finite-prime terms
either inside `Q_S` or displayed separately — never absorbed into notation.

## Dependencies the Lemma needs, and where each must come from

| # | needed | candidate source |
|---|---|---|
| D1 | `P_S(g)` is an **ordinary finite positive** trace | Gram identity `Tr(ϑ_S(g)Π_Sϑ_S(g)^*) = ‖ϑ_S(g)Π_S‖²_{HS}`, valid for any orthogonal projection; finiteness must be checked, not assumed |
| D2 | `Π_S` is the **actual** orthogonal projection onto `M_S` | the formula above, verified (self-adjoint, idempotent, correct range), plus `K` boundedly invertible |
| D3 | an identity `Tr(ϑ_S(f)Π_S) = W_S(f) + E_S(f)` with **`E_S` defined independently of `W_S`** | **this is the open link.** Two candidate routes, below |
| D4 | a bound on `E_S` with explicit support, normalisation, quantifiers and finite-place dependence | nothing in the checked sources |
| D5 | `n = 2` is the only active prime power, and does not vanish identically | window arithmetic on `(log2)/2 < L < (log3)/2`; must avoid the boundary `L = (log2)/2` where the prime block is empty |

## The two candidate routes for D3, named in advance

**Route A — the geometric defect (Connes 1999).** C99-Thm VII.4 is proved for `k = Q`, `S = {∞,2}`,
but for `R_Λ = P̂_Λ P_Λ`, a product of two orthogonal projections. Writing `Q_Λ` for the orthogonal
projection onto `ran P_Λ ∩ ran P̂_Λ` and `G_Λ := P̂_Λ P_Λ − Q_Λ`, the identity
`Tr(Q_Λ U(h)) = 2h(1)log′Λ + Σ_{v∈S}W_v(h) − Tr(G_Λ U(h)) + o(1)` would define the discrepancy
`E_Λ(h) := Tr(G_Λ U(h))` **purely by two-projection geometry**, with no reference to `W`. To be
checked: the structure of `G_Λ`, its trace norm, whether `Q_Λ U(h)` is trace class, whether the
divergent `2h(1)log′Λ` can be cancelled, and whether C99's cutoff convention (`|x| > Λ`, doubly
low-pass) is the Sonin convention (`|x| < λ`, doubly high-pass) or its complement.

**Route B — the metric defect (CCM2024).** Realise everything in the common unitary picture using
`U_∞ = F_μ w_∞` and `U_S = F_μ w_S`, both unitary. Then `J` becomes multiplication by
`m_2(s) = 1 − 2^{−1/2−is}` (CCM2024-Prop 4.6(ii), eq. 57 — to be confirmed as a *modulus-squared*
for `H`, and **not** reused for `η_S`, whose multiplier is the conjugate-inverse). Compute
`Tr(ϑ_S(f)Π_S)` exactly and ask whether CC2021-Thm 4.6 then delivers `W_∞ + W_2`.

## Why this is not the desired positivity renamed

The target asserts an inequality between two **separately computable** quantities: `P_S` is a
Hilbert–Schmidt norm of an explicit operator, and `E_S` is prescribed by projection geometry
(Route A) or by the explicit multiplier (Route B) — in neither case by subtracting `W_S`. If `E_S`
can only be obtained as `Tr(R_ΛU(h)) − Tr(Q_ΛU(h))` with the first term replaced by its arithmetic
value, the Lemma is circular and must be reported as such rather than recorded as a bound.

## Normalisation and constants fixed in advance

Additive coordinate `u = log x`; `d^*x = dx/x`. Mellin pairing `F_μ(g)(s) = ∫ g(λ)λ^{-is}d^*λ`
(CCM2024, proof of Prop. 4.6(ii)), under which **multiplication by `p^{-is}` is the scaling
`λ ↦ λ/p`**. The repository's Weil form is `M = Pole + Arch − LogPi − Prime` with `Prime` assembled
from von Mangoldt terms with `log n < 2L` (`src/zeta_window.py:148,214`); `W_∞` is the
`Arch − LogPi` block alone. Calibration window `(log 2)/2 < L < (log 3)/2`, i.e.
`0.3465735… < L < 0.5493061…`. The boundary `L = (log 2)/2` is **excluded**: there the prime block
vanishes identically and the test is vacuous.

## Comparison the result must beat

The baseline is the ambient worst-case bound obtained from `m ≤ |m_2(s)|² ≤ M`, i.e. the condition
number of the multiplier, together with CC2021-Thm 6.11's published rank-one majorant estimate.
A new expression for `K^{-1}`, a Schur complement, or an orthogonal projection counts as a **tool**
only; to count as a result it must improve or resolve a specific source-level estimate.

## Acceptable alternative outcome

A quantitatively sharper bound on an actual correction term, with the remaining allowance stated
explicitly. Such a bound is useful even if it does not prove positivity.

## Stop rule

If an essential identification fails, record the precise failing term or counterexample and stop
this route. Do not conclude that projection methods in general are impossible, and do not substitute
a different construction after seeing the failure.
