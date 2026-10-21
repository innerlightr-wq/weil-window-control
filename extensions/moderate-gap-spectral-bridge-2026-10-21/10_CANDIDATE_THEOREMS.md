# 10 — The candidate theorems T1–T7: proofs, counterexamples, missing lemmas

Each candidate gets: the statement reached, its provenance class and verification
tier, and either a proof, a counterexample, or the exact missing lemma.

---

## T1 — Exact Gram-framework correspondence

**Reached.** `[class II for the existence, class IV for the quantitative content]`

> **Statement.** ARTZ's `G` is the finite compression, to
> `W = span{f_x, g_z, h_z} ⊂ L²(−1/2,1/2)`, of the zero side of Weil's Hermitian form
> in the coordinate `z_ρ = i(ρ−1/2)log T/(2π)`. The dilation `u ↦ u/log T` identifies
> `W` isometrically (up to the Jacobian) with the real even sector of `L²(−L_eff,L_eff)`,
> `L_eff = (log T)/2`, used by my finite-window program. Under this identification the
> horizontal displacement of a zero appears **only** through the single scalar
> ```
>     kappa_rho = int eta_tau^2(u) cosh( 2 delta_rho u log T ) du  >= 1 ,
> ```
> via `‖g_ρ‖² = (1+κ_ρ)/2`, `‖h_ρ‖² = (κ_ρ−1)/2`, `⟨g_ρ,h_ρ⟩ = 0`,
> `‖g_ρ‖² − ‖h_ρ‖² = 1`; and
> ```
>     kappa_rho - 1  =  2 delta_rho^2 V_mu + O(delta_rho^4) ,
> ```
> where `V_μ` is the **centred second moment of the window measure** — the same
> constant, and the same description, as the coefficient of my published `δ²` law.

**Attribution.** The *existence* of the correspondence is Alpöge–Furman's own stated
mechanism ("a finite compression of Weil's Hermitian form, with Sylvester's law of
inertia handling off-line pairs"), [AF, abstract; §3]. It is **not** new. What is new
is the `κ_ρ` identity, the Gram table of 02 §3, and the coefficient match.

**Proof.** 02 §3, 03 §1. `[certified: verify.py D1, D4–D7, I1–I3, I5, I6]`

**What it does not preserve** (each individually fatal, 03 §3): positivity, moment
order, arithmetic input, support scale. **Verdict: D.**

---

## T2 — Moderate-gap conditioning bound

**Refuted in the uniform form; a certified *upper* bound replaces it.**
`[class IV, certified]`

> **Statement (Proposition 2 of 05).** For equispaced nodes `x_j = jc`,
> `0 ≤ j ≤ d−1`, `m = d−1`, and `G_d = [K_τ(x_i−x_j)]`,
> ```
>     lambda_min(G_d)  <=  ( 2 sin(pi c/2) )^{2m} / C(2m,m)   ~  sqrt(pi m) ( 2 sin^2(pi c/2) )^m ,
> ```
> with `2sin²(πc/2) < 1` for every `c < 1/2`; at `c = 1/4` the rate is
> `0.2928932…`. Hence under the separation ARTZ guarantees, conditioning degrades
> exponentially in the block length.

**Proof.** Binomial test vector; 05 §2. `[certified H1–H4]`
**Counterexample to the hoped-for statement.** The hoped-for `λ_min ≥ λ(c) > 0`
uniformly in `d` is false by the display above; empirically at `c = 1/4`,
`λ_min = 9.97e−2, 2.91e−4, 5.96e−7, 1.11e−9, 1.98e−12` for `d = 2,4,6,8,10`.

**The deeper refutation** is Proposition 1 of 05 (anti-localisation): off-line zeros
have internal ordinate gap exactly 0, so the moderate-gap theorem certifies
separation only at ordinates where there is provably no off-line zero.
`[class IV, proved; verify.py J1–J5]`

**Verdict: E as an RH transfer; C as a standalone finite-configuration estimate.**

---

## T3 — Displacement-sensitive spectral defect

**Reached, new, unconditional — and not summable.** `[class IV]`

> **Statement (the refined budget).** With ARTZ's notation, for every fixed
> `τ ∈ (0,1/200]` and `0 < c ≤ 1/4`,
> ```
>     tr Q  +  K_tau(c)^2 |Rtilde_1|  +  sum_{q_j>0} ( q_j - 2 )^2
>              <=  ||G||_F^2 - N(T)  =  ( C_{eta_tau} - 1 + o(1) ) N(T) ,
> ```
> where `q_j` are the positive eigenvalues of `Q`. For a single simple off-line pair
> of depth `δ`, `Σ_j(q_j−2)² = (κ(δ) − 1)² = 4V_μ²δ⁴ + O(δ⁶)`.

**Proof.** 04 §1; the term is the exact slack in `q² ≥ 4q−4` which ARTZ's Lemma 5
discards. Recovered at zero cost. `[proved; certified C1–C3, E1–E4; checked on 7
synthetic multisets, verify.py F4]`

**What it buys, and what it does not.**
* Rigorous, unconditional, for a hypothetical isolated off-line zero:
  `2m + (m(1+κ)−2)² ≤ (C_{η_τ}−1+o(1))N(T)`. **Vacuous** by a factor
  `(2/π)(C_MT−1)log(T/2π) = 0.2084930…·log(T/2π)`, which exceeds 1 for `T > 761`
  `[certified G1–G5]`.
* **Missing lemma, and it is refuted, not merely missing:**
  `(H-add)`: `Σ_j(q_j−2)² ≥ c_0 Σ_l(κ_l−1)²` for some absolute `c_0 > 0`.
  Measured false — at mean-spacing separation the left side grows like
  `(0.98 + 1.33 ln k)(κ−1)²` while the right side grows like `k(κ−1)²`
  `[V; results/additivity.txt]`. And `(H-add)` with any `c_0>0` would imply the
  **Density Hypothesis** `N(σ,T) ≪ T^{2(1−σ)+ε}` for every `σ ∈ (1/2,1)`
  `[certified G6, G7]`, which is open for `σ ∈ (1/2, 25/32)`. So `(H-add)` is not a
  gap to be filled.

**Verdict: C.**

---

## T4 — Combined spacing/displacement inequality

**Closed.** `[class IV]`

Any inequality `Σ_ρ D(ρ)·Φ(local spacing at ρ) ≤ budget` with `Φ` bounded below on
ordinary spacings fails on configuration C4 of 09 §2: separation exactly twice the
mean, zero local density excess, and the left side collapses to `O(log k)`. The
interference is **not** a near-collision effect; it survives where there is no
collision. A `Φ` vanishing at ordinary spacing makes the inequality vacuous.
`[refuted by explicit configuration]`

**Verdict: E.**

---

## T5 — Improved certified finite-window witness

**Not reached, and disqualified on scaling grounds.** `[class VI]`

My certified witnesses are at **fixed** support: `L = 4, d = 184`
(`Tail ≤ 7.339 < 4D²κ = 14.771`, `δ₀ = 0.29369`), `L = 8, d = 80`
(`δ₀ → 0.003011`, then `4.0e−20` after imposing all signed odd moments),
and the geometric-side certificate `Q_ζ(f) ⪰ 2.237e−36` (width `1.07e−51`).
Every route from ARTZ needs `L_eff ≍ log T → ∞` (08 §4 item 1) and a mode count
`d ≍ (log T)²` (item 2). **Nothing here improves a fixed-`L` certificate, a detector
constant, or a tail estimate.** No attempt was made to rerun the certified
machinery, in accordance with the brief's instruction not to modify it.

**Verdict: E (not applicable).**

---

## T6 — Conditional pointwise theorem

**One honest conditional statement exists; its hypothesis is not independently
meaningful.** `[class VI]`

> **Statement.** Assume `(H-add)` with constant `c_0 > 0`. Then for every fixed
> `σ ∈ (1/2,1)` and `ε > 0`, `N(σ,T) ≪_{ε,τ,c_0} T^{2(1−σ)+ε}`.

This is a *conditional density* theorem, not a pointwise one, and it fails the
brief's own test for a useful conditional theorem: the missing hypothesis `(H-add)`
is (i) **measured false** (04 §3.1) and (ii) **equivalent in strength to the Density
Hypothesis**, i.e. "RH-adjacent in another form". Both disqualify it.

A genuinely pointwise conditional statement would need, in addition, a localisation
hypothesis placing a hypothetical off-line zero inside a certified block — and 05
Proposition 1 shows that **no gap hypothesis can ever supply one**, because the
required gap is identically zero. So the conditional route is not merely unproved;
its missing hypothesis is impossible for structural reasons.

**Verdict: E.**

---

## T7 — No-go theorem

**Reached, with three independent proofs.** `[class IV]`

> **Theorem (no uniform pointwise exclusion).** The spectral-budget transfer cannot
> yield a uniform exclusion of an individual off-critical-line zero. Three
> independent reasons:
>
> 1. **Deficit.** The only pointwise bound available is
>    `2m + (m(1+κ)−2)² ≤ (C_{η_τ}−1+o(1))N(T)`. Since `β ≤ 1` forces `δ ≤ 1/2` and
>    `κ ≤ cosh(δ log T)` for every admissible window, the maximal left side is
>    `≈ T/4` while the budget is `≈ 0.0521·T·log(T/2π)`. The deficit is exactly
>    `(2/π)(C_MT−1)log(T/2π)`, i.e. **one power of `log T`** at `δ = 1/2` and
>    `T^{1−2δ}log T` for smaller `δ`. The budget is extensive (one unit per zero);
>    the single-zero gain is at most `T^{2δ} ≤ T`, sub-extensive by one logarithm.
>    No choice of window changes this: the bound `κ ≤ cosh(δ log T)` uses only
>    `η_τ² ≥ 0`, `∫η_τ² = 1`, `supp η_τ² ⊂ (−1/2,1/2)`. `[certified G1–G5]`
> 2. **Anti-localisation.** An off-line zero `ρ` always has the partner `1−ρ̄` at
>    ordinate distance **0**. Therefore it is never an interior point of a run of
>    positive gaps, and ARTZ's Proposition 1 builds runs precisely by discarding the
>    clusters that contain it. Any gap theorem strong enough to separate `ρ` from its
>    partner would already prove there is no such `ρ`. `[proved; verify.py J1–J5]`
> 3. **Non-additivity.** The only `δ`-sensitive functional in the chain,
>    `Σ_j(q_j−2)²`, is log-additive over off-line pairs at ordinary spacing, with the
>    measured law `0.98 + 1.33 ln g` (max residual `0.088`, versus `1.040` for the
>    best linear fit). The mechanism is the alignment of the `g`-directions as the
>    depth grows (`|corr| = 0.751` at separation 1, depth `a = 2π`; `0.884` at
>    `a = 3π`). `[V, with the mechanism proved]`

**Verdict: E — and this is the useful outcome the brief anticipated.**

---

## Summary table

| candidate | reached? | class | tier | verdict |
|---|---|---|---|---|
| T1 Gram-framework correspondence | yes (existence already in [AF]) | II + IV | certified | **D** |
| T2 moderate-gap conditioning | refuted; replaced by a ceiling | IV | certified | **E** / C standalone |
| T3 displacement-sensitive defect | yes, new inequality | IV | proved + certified | **C** |
| T4 combined spacing/displacement | refuted by explicit configuration | IV | diagnostic + proved | **E** |
| T5 improved certified witness | no (scaling) | VI | — | **E** |
| T6 conditional pointwise theorem | hypothesis refuted and DH-equivalent | VI | — | **E** |
| T7 no-go theorem | yes, three independent proofs | IV | certified + proved + V | **E** (useful) |
