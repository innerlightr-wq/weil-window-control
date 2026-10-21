# 06 — Against the confluent blind-spot obstruction

Answers the brief's §7 (A–E). The governing prior assessment is the 2026-10-15
**OUTCOME D, DECISIVE NO-GO** (`~/scratch/rh-blind-spot-covering-2026-10-15-142629`);
Route D is not reopened here. The question asked is narrower: **does the moderate-gap
theorem supply a new hypothesis that addresses `NG-2` or `NG-3` directly?**

## 0. The exact prior definitions, recovered from source

From `PROBLEM_FORMULATION.md` §2 (line 61):

```
   kappa(Gamma, gamma_*, L)
        = sup { F'(gamma_*)^2 : f in A_L , F(gamma_i) = 0 for all i }
        = dist^2( s , C(Gamma) )                              [INTERP_PROOF Thm 1]
```
with `A_L` the admissible class (real, even, `supp f ⊂ [−L,L]`, `‖f‖₂ = 1`), `s` the
displacement signal and `C(Γ)` the span of the competing spectral signals.
`κ` is **already a supremum over the whole admissible class**, so any bound on it
binds every witness family at once, of any cardinality, adaptive or not.

Confluent cluster bound (`weil-confluent-cluster-theorem-audited`, audited GREEN):
```
   kappa(Gamma)  <=  ( prod_i h_i )^2 * 2 L^{2n+1} / ( (2n+1) (n!)^2 )
   sharp law:    kappa  =  ( prod_i h_i )^2 / (n!)^2 * D_n
```
with `h_i` the gaps from `γ_*` to the `n` competing ordinates. The **exponent
`2(n−1)` is Batenkov-known**; only the constant is ours; never headline the exponent.

Decisive prior no-gos, verbatim content:
* `NG-0` (identity): a witness measures only `t² − β²` — vertical ordinate offset and
  horizontal zero displacement enter at the same order, same coefficient, opposite
  signs.
* `NG-2` (cluster, decisive): one extra ordinate at distance `g` gives
  `δ ~ η√(2S/D₂)/g → ∞`; `n = 2` suffices, cluster size need not grow with `M`.
* `NG-2'` (collision, sharp): at multiplicity `ν ≥ 2` the margin is **exactly 0**,
  because translating a multiple zero gives precisely `s = c^{(1)} ∈ J_{ν−1}`.
* `NG-3` (portability): `δ = η√(1 + S/(2κ)) > η` always.

## 1. (A) Separated critical-line nodes — do guaranteed moderate gaps exclude the
worst confluent degeneracies for the *selected* nodes?

**Yes for the selected nodes, and it is useless.** `[VI → conditional]`

On a run `j ∈ W_ℓ(T,c)` the `ℓ−1` interior ordinates are simple, on the line, and
mutually separated by `≥ 2πc/log T`. In the normalised variable the gaps satisfy
`h_i ≥ c`, so the *upper* bound above no longer collapses:
`κ ≤ c^{2n}·2L^{2n+1}/((2n+1)(n!)²)` is no longer driven to 0 by `∏h_i → 0`.

But (i) an **upper** bound on `κ` is useless — detection needs a **lower** bound; and
(ii) the matching lower bound was never completed: the Hermite–Genocchi remainder
noted in `weil-confluent-cluster-theorem-audited` is **still open** (dependency D3 of
§01). So the best available statement is

> *Conditional on the unfinished Hermite–Genocchi lower bound*, a run of `ℓ ≤ 4`
> `c`-moderate gaps gives, at a positive proportion (`≥ 0.0175` for `ℓ = 4`) of
> indices, `κ ≳ c^{2n}D_n/(n!)²` for the `n ≤ ℓ−1` selected nodes.

**class VI, not IV.** And it is for *on-line* nodes only — which brings us to (B).

## 2. (B) A nearby hypothetical off-line quartet — can its displacement be detected
without an uncontrolled loss of sensitivity?

**No, and the loss is now quantified twice over.** `[IV]`

* **By Proposition 1 of §05**, an off-line pair has internal ordinate gap exactly 0.
  So in `κ(Γ)` the competing-signal set `Γ` always contains a node at distance 0 from
  the target — the worst case of the cluster bound, `∏h_i = 0`, is not an adversarial
  choice but a **structural certainty**. The moderate-gap theorem cannot alter it,
  because the gap it would have to enlarge is identically zero.
* **Independently**, the off-line detection strength was measured (and certified) in
  my own program at `M(δ) = δ⁴·D₂`, `D₂ = 2911.050006543213352847841125028374273138
  ± 1e-37`. The present framework independently gives the same exponent: the budget
  cost of an off-line pair is `(κ_ρ − 1)² = 4V_μ²δ⁴ + O(δ⁶)` (03 §3, `verify.py`
  E5–E6). Two unrelated constructions, same exponent 4.

**Correction to a tempting but wrong transfer.** `NG-2'` says the margin is *exactly*
0 at multiplicity `ν ≥ 2`. It is tempting to apply it verbatim to an off-line pair,
since the pair also shares an ordinate. **That is wrong.** `NG-2'` is about a
*vertical* collision (a repeated zero), where the displacement direction `s` lands
exactly in the jet `J_{ν−1}` spanned by the collision. An off-line pair's collision
is in the ordinate only; the displacement direction is `cosh(δ·)`, which is *not* in
the on-line jet, and detectability is `δ⁴D₂ > 0`, not 0. The correct statement is
"strictly positive but two orders weaker in the small parameter", which is erratum E3
of `weil-offline-detect` and is preserved here. `[III, corrected]`

## 3. (C) Clusters outside the selected window — can remote or unselected nodes
recreate a near-singular direction?

**Yes, and the far field is now measured on both sides.** `[III + V]`

* Prior measurement: the far field `γ₂..γ₅₀` changes `D₂ = 2911.050` to
  `D₂^Q = 2897.142`, i.e. removes only `0.478 %`; the exponent is unchanged. So the
  blind spot is not an artefact of ignoring other zeros.
* New measurement (04 §3.1): in the ARTZ framework the *interference between
  off-line pairs* at separations of order the mean spacing is enormous — the
  normalised `g`-direction correlation reaches `−0.751` at separation 1 and `−0.884`
  at depth `a = 9.42`. Unselected off-line nodes do not merely recreate a
  near-singular direction; they make the budget functional log-additive.
* And §05 §4: adversarial *on-line* placements do **not** help the adversary (they
  raise `trΨ(P)`). The leverage is purely in other off-line pairs.

## 4. (D) Arbitrarily close target ordinates — does the theorem restrict the target
from an almost-collision?

**It imposes no restriction, and the question is moot: for an off-line target the
collision is exact and unavoidable.** `[IV, proved]`

This is the sharpest new statement in this file. In the earlier work, "the target
participates in an almost-collision" was the *adversarial* case (`NG-2`). Here it is
the **only** case. The functional equation supplies the colliding partner for free.
A theorem that guaranteed *every* ordinate were `c`-separated from every other would
simultaneously prove that there are no off-line zeros and no multiple zeros — i.e. it
would be strictly stronger than RH + simplicity. So no gap theorem can ever relax
this, and the deficiency is not a gap in the literature.

## 5. (E) Multiple zeros — does multiplicity produce exact degeneracies not controlled
by moderate gaps elsewhere?

**Yes; `NG-2'` applies verbatim and is untouched.** `[III]`

Multiple zeros sit in `R_2` and contribute to `tr Q` exactly as off-line pairs do
(the "flat charge 4" of Alpöge–Furman §7.2(c)); ARTZ's framework cannot distinguish
a multiple on-line zero from an off-line pair of depth `→ 0` — AF say so explicitly:
> "Replacing the doubles by off-line pairs of depth → 0 (**spectrally the same**)
> gives the extremal for Theorem A."

So the extremal configuration of the whole method is the one on which the detector's
margin is exactly 0. **The method's extremal and the detector's blind spot are the
same configuration.** That is the cleanest statement of why the two programs cannot
help each other. `[IV, proved from II + III]`

## 6. Counterexamples to the over-strong uniform statements

| over-strong statement | counterexample / refutation |
|---|---|
| "moderate gaps give a uniform lower bound on `κ`" | needs the unfinished Hermite–Genocchi remainder (D3); as stated it is class VI |
| "`λ_min(G_d) ≥ λ(c) > 0` uniformly in `d` under separation `c`" | **false**: `λ_min ≤ √(πm)(2sin²(πc/2))^m`, 05 Prop 2 `[certified]` |
| "off-line pairs are a near-collision, hence `NG-2'` gives margin 0" | **false**: margin is `δ⁴D₂ > 0`; see §2 |
| "the slack `Σ(q_j−2)²` is additive over off-line pairs" | **false**: log-additive, 04 §3.1 `[V]` |
| "the gap theorem localises a run near any prescribed ordinate" | **false**: positive proportion only, and off-line ordinates are excluded (§05 Prop 1) |

## 7. Node-span vs derivative-jet: the hypotheses are not interchangeable

Preserved verbatim from the prior work, and relevant here. The spectral node-span
version (`κ = dist²(s, C(Γ))` with `C(Γ)` spanned by *distinct* ordinates) and the
derivative-jet version (`C(Γ)` spanned by a jet `c^{(0)},…,c^{(ν−1)}` at a repeated
ordinate) have different rank and independence hypotheses, and the confluent theorem
moves between them only through the Newton basis. ARTZ's construction is a
*node-span* object (distinct `z ∈ Z_T`, with multiplicities carried as scalar weights
`m_z`, not as jets). **Do not read ARTZ's multiplicity weights as a jet.** The
Sylvester block (44) is `diag(m_x, 2m_z, −2m_z)` — scalar weights on independent
vectors, not derivatives. `[III, preserved]`
