# 08 — Parameter scaling: exact conversion and the failure modes

Answers the brief's §9.

## 1. The conversion, derived

ARTZ normalise by the local mean spacing `2π/log T`. Writing `L := log T/(2π)`:

```
    z_rho    = i (rho - 1/2) L           so   Re z = - gamma L ,  Im z = delta L
    u_paper  in (-1/2, 1/2)              (the conjugate variable)
    u_mine   = u_paper * log T           so   u_mine in ( -(log T)/2 , (log T)/2 )
    L_eff    = (log T)/2                 the half-window my program would call L
    a        = delta * log T             the unique dimensionless displacement
```
Derivation of the middle line: ARTZ's phase is `e^{-2πi z u_paper}` with
`z = i(ρ−1/2)L`; mine is `e^{-iγ u_mine}`. Matching `2π z u_paper` against
`γ u_mine / 1` for an on-line zero (`z = −γL`) gives `u_mine = 2π L u_paper =
u_paper log T`. `[IV, proved]`

## 2. The table

| ARTZ parameter | value / range | my parameter | correspondence | assumption needed | uniform? | failure mode |
|---|---|---|---|---|---|---|
| `τ` (their `δ`) | `0 < τ ≤ 1/200` | — (I have no cutoff; my window is the witness) | none | τ fixed **before** `T → ∞` | yes in `T`, at fixed `τ` | `τ → 0` first is outside Lemma 3.2 |
| `η_τ`, `supp ⊂ (−1/2,1/2)` | fixed | `f`, `supp ⊂ [−L,L]`, `‖f‖₂=1` | `u ↦ u/log T` | — | — | my certificates are at `L ∈ {4/5, 4, 8}`, fixed; here `L_eff ≍ log T` **grows** |
| support length `1` | fixed | `2L` | `2L ↔ log T` | — | — | a transfer needing fixed `L` cannot use this |
| `K_τ(0) = 1` | normalisation | `‖f‖₂ = 1` | mass `1` vs mass `~L` | window mass `L^s`, `s=0` vs `s=1` | — | the `δ²` coefficient power changes `L²` ↔ `L³` (07 §5) |
| mode count — | ARTZ have none | `d` (modes) | `d ↔ ℓ+1 ≤ 5` for a run; `d ≍ L log T/(2π)` for a window | — | no | `λ_min` decays like `(2sin²(πc/2))^{d-1}` (05 Prop 2) |
| `c` (gap parameter) | `0 < c ≤ 1/4` | normalised node separation | identical (both normalised) | — | yes | `c ≤ 1/4` is 4× **sub**-Nyquist |
| `ℓ` (run length) | `2 ≤ ℓ ≤ 4` | `n` competing nodes `= ℓ−1 ≤ 3` | interior nodes only (05 Prop 1) | — | yes | endpoints not certified |
| `N(T)` | `(T/2π)log(T/2π)+O(T)` | — | the budget scale | Riemann–von Mangoldt | yes | budget is extensive; single-zero gain is not |
| `δ_ρ = β−1/2` | not named in ARTZ | `δ` | `Im z_ρ = δL` | — | **no** | only `a = δ log T` is dimensionless |
| `κ_ρ` | not named in ARTZ | — (new, 02 §3) | `κ = ∫η_τ²cosh(2au)du` | — | yes as a function of `a` | quadratic only for `a = O(1)` |
| `C_{η_τ}` | `→ C_MT = 1.32749929632…` | — | budget constant | Lemma 3.2 | yes at fixed `τ` | no explicit rate `C_{η_τ} → C_MT` is given |
| `x` in Montgomery's `F(x,T)` | `1 ≤ x ≤ T` | — | `x = T^{u}`, `|u| ≤ 1` | BGSTB | yes | support `> 1` needs the pair correlation conjecture |

## 3. What happens as `T → ∞`

Three regimes of the single parameter `a = δ log T`, with the depth scalar
`κ(a) = ∫η_τ²(u)cosh(2au)du`:

```
   a = O(1)       (delta = O(1/log T))   kappa = 1 + 2 M_2 a^2 + O(a^4)   bounded
                   <-- BGSTB's own hypothesis |beta-1/2| <= c/log T
   a -> infinity, delta -> 0             kappa ~ c e^{(1-2 tau) a} / a    exponential in a
   delta fixed > 0                       kappa ~ T^{delta} * e^{-C_tau sqrt(log T)}
                                            (saddle point at u = 1/2 - O(1/sqrt(a));
                                             sub-polynomial, not polynomial, correction)
```
`[IV, proved; the `a=O(1)` branch certified E5]`

**The first regime is the only one in which the `δ²` law is the right expansion, and
it is exactly the regime BGSTB had to assume to get kernel positivity.** This is the
single most important scaling fact in the whole investigation: the `δ²–L³`-style
sensitivity analysis and the unconditional pair-correlation framework agree *only*
inside the hypothesis that the latter was built to avoid.

## 4. Failure modes, explicitly

1. **Growing support.** Any improvement here requires `L_eff = (log T)/2 → ∞`. My
   certified witnesses (`L = 4, d = 184`; `L = 8, d = 80`; `L = 4/5, γ = 14`) are at
   *fixed* `L`. A transfer that needs growing support is not an improvement of an
   existing certificate; it is a different theorem. **Disqualifies Candidate T5.**
2. **Growing mode count.** Annihilating the competing ordinates in a window of
   half-length `L_eff` needs `d ≳ 3N` modes where `N ≍ L_eff log T/π` is the ordinate
   count (`weil-projected-coercivity-closed`: minimal span `d = N+1` realises only
   `6e−6` of `κ`). So `d ≍ (log T)²`, and `λ_min` of the node Gram is then below
   `(0.293)^{(log T)^2}`. **Disqualifies the conditioning route independently of
   §05 Prop 1.**
3. **Parameter-dependent error term.** Lemma 3.2's error is `O_τ(1/√log T)N(T)`, with
   an unquantified `τ`-dependence, and no rate is given for `C_{η_τ} → C_MT`. Every
   statement here is therefore made at **fixed `τ = 1/200`** and the constant
   `C_{η_τ}` is carried symbolically, with `C_MT` used only where the direction of
   the inequality makes it safe.
4. **Order of limits.** ARTZ's own order is `T → ∞` then `τ → 0`. Reversing it is not
   licensed. Every bound in 04 and 05 respects the licensed order.
5. **Support `> 1`.** Any attempt to enlarge the Fourier support beyond `[−1,1]`
   (equivalently `x > T`) requires Montgomery's pair correlation conjecture or
   Hardy–Littlewood (AF §7.2(a)). That is a conditional route and is excluded by the
   brief's stop rules.
