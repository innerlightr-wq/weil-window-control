# 07 — The second-order displacement expansion

Answers the brief's §8, item by item. This was flagged as the highest-risk route;
it delivers the most, and still not enough.

## 1. Which quantity should be expanded

For a conjugate off-line configuration `ρ_± = 1/2 ± δ + iγ` the normalised
coordinates are `z_± = (−γ ± iδ)L`, `L = log T/(2π)`. There is exactly one scalar in
ARTZ's construction that depends on `δ`:

```
    kappa(delta)  =  K_tau( 2 i delta L )  =  int eta_tau^2(u) cosh( 2 delta u log T ) du .
```

Everything else — the trace, the inertia, the block structure, `Ψ`, `n_±` — is
`δ`-blind (04 §2, `verify.py` E2). The signed defect contribution of the pair is

```
    Q_off(delta) = 2m ( v_g v_g^T - v_h v_h^T ) ,
    ||v_g||^2 = (1 + kappa)/2 ,  ||v_h||^2 = (kappa - 1)/2 ,  <v_g, v_h> = 0 ,
    spec = { m(1+kappa) , -m(kappa-1) } .
```
`[IV, certified: verify.py D4–D7, E1]`

So the object to expand is **not** `Q_off` as an operator (its trace and inertia are
constant in `δ`; only its *eigenvalue magnitudes* move) but the scalar `κ(δ)`, and
through it the Frobenius-type functionals.

## 2. Is an even expansion justified?

Yes, and exactly. `cosh` is even and entire, `η_τ²` is even with `∫η_τ² = 1` and
compact support, so `κ` is an even entire function of `δ` with

```
    kappa(delta)  =  1 + sum_{p >= 1} (2 delta log T)^{2p} M_{2p} / (2p)! ,
                     M_{2p} := int u^{2p} eta_tau^2(u) du  >  0 ,
```
every coefficient positive, no odd terms, convergence for all `δ` (entire). The
symmetry hypothesis the brief asks about (`Q_off(δ) = Q_off(0) + δ²H + O(δ⁴)`) is
therefore satisfied *in the strong sense*: the expansion is not asymptotic, it is a
convergent even power series. `[IV, proved]`

```
    kappa(delta) - 1  =  2 delta^2 V_mu  +  O(delta^4) ,
    V_mu := (log T)^2 M_2  =  centred second moment of the window measure
            in the unnormalised frequency variable v = u log T  (total mass 1).
```

For the operator: `Q_off(δ) = Q_off(0) + δ²H + O(δ⁴)` with
`H = 2m·V_mu·(v_{g}v_{g}^T + v_{h}v_{h}^T)|_{δ=0}·(1 + …)` — note the **plus** sign:
both the `g` and `h` norms grow with `δ`, in opposite spectral directions, so `H` is
*not* sign-definite as an operator while the two eigenvalue magnitudes both increase.
`[IV, proved]`

## 3. Does the leading term have a definite sign?

**In the trace: no — it cancels identically.**
`tr Q_off(δ) = m(1+κ) − m(κ−1) = 2m` for every `δ`. The `δ²` terms cancel exactly.
This is the algebraic heart of the Sylvester blindness.

**In the Frobenius norm: yes, strictly positive.**
`‖Q_off(δ)‖²_F = 2m²(1 + κ²) = 4m² + 8m²V_μ δ² + O(δ⁴)`, monotone increasing.
`[certified E1]`

**In the budget: yes, strictly positive, and at fourth order.**
The quantity actually available in the inequality chain is the discarded Lemma-5
slack
```
    D(rho)  =  ( m(1+kappa) - 2 )^2 ,    D = (kappa-1)^2 = 4 V_mu^2 delta^4 + O(delta^6)
                                                     (for m = 1).
```
`[IV, certified E3, E5, E6]`

## 4. Does the coefficient relate to my `δ²` sensitivity constant?

**Yes — it is the same constant.** My published law
(`doi:10.5281/zenodo.23250932`, §"The δ²–L³ law") states that its coefficient "is
twice the centred second moment of the window". Here `κ − 1 = 2δ²V_μ` with `V_μ` the
centred second moment. Checked exactly for `η² ≡ 1`: `κ − 1 = (2πy)²/6` and
`2δ²V_μ = δ²(log T)²/6` coincide identically, since `2πy = δ log T`
`[certified I1, I2]`. The function-field companion identity
`2Σ_{i=0}^R(i − R/2)² = \binom{R+2}{3}` is reproduced exactly for `R = 1..200`
`[certified I4]`.

## 5. Is the support-length power `L³`, another power, or normalisation-dependent?

**Normalisation-dependent, and the bookkeeping is now explicit.** With window mass
`~L^s` the coefficient scales as `L^{s+2}`:

| framework | window | mass | coefficient | power |
|---|---|---|---|---|
| ARTZ / Lamzouri | `η_τ²`, probability density | `1` | `2V_μ` | `L²` |
| my published note | unnormalised, `‖f‖₂ = 1` on `[−L,L]` | `~L` | `c(L) = 8L³(1/3 + 1/√5) ≈ 6.244L³` | `L³` |

`[IV, certified I5, I6]`. This confirms, and makes precise, the earlier finding that
the `δ²–L³` law is **not universal**: the power is `s+2`, set by the window's mass
normalisation, and `L³` is the `s = 1` member of the family. Any claim that the
paper's framework "exhibits the `L³` law" would be wrong; it exhibits the `s = 0`
member, `L²`.

## 6. Do multiplicity and close-zero limits commute with the small-`δ` expansion?

**No, in both cases, and the failures are of different kinds.** `[IV]`

* **Multiplicity.** `κ` is independent of `m`, but the eigenvalues are `m(1+κ)` and
  `−m(κ−1)`, so `D(ρ) = (m(1+κ)−2)²` and `lim_{δ→0}D = (2m−2)²`, which is `0` for
  `m = 1` but `4(m−1)² > 0` for `m ≥ 2`. The two limits `δ → 0` and `m → 1` do not
  commute in the sense that `D` does not vanish as `δ → 0` unless `m = 1`. An
  off-line pair of depth `→ 0` is spectrally a multiple on-line zero (AF §7.2(b)),
  which is why the method cannot separate them.
* **Close zeros.** Measured, and the failure is severe: at fixed `δ`, letting a
  second off-line pair approach at separation `D` changes the *total* slack
  non-monotonically and, at mean-spacing separation, reduces the per-pair slack by a
  factor that grows with the number of pairs (04 §3.1). So
  `lim_{D→0}lim_{δ→0} ≠ lim_{δ→0}lim_{D→0}`. **This is the obstruction**, and the
  mechanism is the alignment of the `g`-directions (04 §3.2), not a conditioning
  failure.

## 7. Is the resulting estimate uniform in `T`?

**No, and the non-uniformity is the decisive point.** The only dimensionless
displacement is `a = δ log T`. All of the above is uniform in `T` *as a function of
`a`*, but:
* the expansion `κ − 1 ≈ 2δ²V_μ = a²M_2·2` is valid only for `a = O(1)`, i.e.
  `δ = O(1/log T)` — precisely BGSTB's own positivity hypothesis
  `|β−1/2| ≤ c/log T`;
* for `a → ∞` (any fixed `δ > 0`) the correct behaviour is
  `κ ≈ T^{δ}·(sub-polynomial)`, not the quadratic;
* the order of limits in the papers is: fix `τ`, let `T → ∞`, then `τ → 0`. Any
  statement needing `τ → 0` *before* `T → ∞` is outside the proved regime.

See 08 for the full table.

## 8. Can the trace or Frobenius constraint absorb the displacement term?

* **Trace: no.** `tr Q` is exactly `δ`-independent (§3).
* **Frobenius: yes for one pair, no for many.** `‖G‖²_F ≥ ‖Q_+‖²_F` (elementary:
  `‖P − Q_-‖²_F ≥ 0` and `tr(PQ_+) ≥ 0`, `verify.py` F5), which for one isolated pair
  gives `(m(1+κ))² ≤ (C_{η_τ}+o(1))N(T)` — rigorous, and vacuous by a factor
  `0.2085·log(T/2π)` (04 §5). For many pairs the Frobenius constraint absorbs only
  `O(log k)` pairs' worth, so it cannot carry a count.

## 9. Verdict for Candidate T3

A genuine second-order coefficient **does** exist, is **positive**, is **computed**,
and **coincides with the certified Weil-form coefficient** of my published law. The
budget cost is its square, giving the `δ⁴` exponent independently confirmed by
`weil-offline-detect`'s `M(δ) = δ⁴D₂`.

**But the resulting inequality is not summable over zeros** (04 §3), so it yields no
new arithmetic restriction. **Verdict: C — a restricted, certified, new quantitative
estimate; not B, because it does not improve the RH machinery; not E, because the
coefficient and the inequality are real and nontrivial.**
