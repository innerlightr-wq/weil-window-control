# 09 — Pointwise versus average: what the unconditional estimates can and cannot do

Answers the brief's §10 (mandatory) and Final Question 4.

## 0. The honest one-line answer

**The result remains distributional.** ARTZ contains exactly three zeta-specific
inputs, and none of them localises anything:

| zeta-specific input | what it is | can it localise? |
|---|---|---|
| Lemma 1 / (12) (BGSTB, via Lamzouri) | **one scalar**: `Σ_{z,s}K_τ(z−s)² = (C_{η_τ}+o(1))N(T)` | no — a single global second moment |
| the functional equation | `Z_T = Z̄_T` | no — and it *helps* the adversary (§05 Prop 1) |
| Riemann–von Mangoldt + `S(t) ≪ log t` | over a normalised ordinate window of length `H`, the zero count is `H + O(log T)` | only the local *density*, see §3 |

Everything else in the paper is linear algebra on an arbitrary conjugation-invariant
finite multiset. A single scalar constraint plus a density constraint cannot pin an
individual zero; the four configurations below show exactly how much room is left.

## 1. Ground rules for the configurations below

Each `C_i` is a **finite conjugation-invariant multiset of complex numbers** with
distinct elements, i.e. it satisfies the formal hypotheses of ARTZ §3–§4 and of
Lamzouri Proposition 2.1 (linear independence of `{f_x, g_z, h_z}` follows from
distinctness). The matrix identities `tr G = N`, `‖G‖²_F = Σ_{z,s}K(z−s)²`, the
inertia count, and the chain of 04 §1 all apply verbatim and were checked numerically
(`verify.py` F1–F5).

**No claim is made that any `C_i` is `Z_T` for any `T`.** In particular the *value* of
`‖G‖²_F` for a synthetic configuration is generally not `≈ C_MT·N`; what the
configurations establish is what the *inequalities* permit, i.e. which conclusions
cannot be drawn from `tr`, `‖·‖²_F`, and the inertia alone. Where a configuration is
only illustrative, it is labelled so.

## 2. The four configurations the brief asks for

### C1 — many favourable blocks coexist with rare unfavourable ones
`[formal hypotheses satisfied; verify.py block F case 2]`
24 simple on-line points at integer (= mean) spacing, plus one off-line pair at
depth `a = 2π`. `N = 26`, `n = 24`, `k = 1`, `r = 0`.
* Runs: every interior index of every `ℓ`-run is one of the 24 on-line points;
  the proportion of indices starting a 4-run is `≈ 21/26 = 0.81`, far above ARTZ's
  guaranteed `0.0175`.
* The off-line pair contributes `tr Q = 2` and `Σ_j(q_j−2)² = (κ−1)² = 1300.2`.
* **And its own two indices are excluded from every run** — they are the only pair of
  equal ordinates in the list.

A positive-proportion theorem with a bound of `0.0175` is therefore consistent with
`k` as large as `0.164·N` off-line pairs (the trace bound), each of them sitting in a
locally unfavourable position. **Favourable density and pointwise failure coexist.**

### C2 — a small off-line population does not violate the proportion theorem
`[illustrative]`
Scale C1: `n = 1000` on-line, `k = 10` off-line pairs spread out. The four-gap
proportion is still `≈ 0.98`, while 10 off-line pairs exist. Even `k = 0.1638·N`
(the extremal allowed by `tr Q ≤ (C_MT−1)N`) leaves the four-gap bound positive.
Nothing in Theorem 1 degrades.

### C3 — an individual off-line zero in a poorly conditioned local configuration
`[formal hypotheses satisfied; results/cluster_smooth.txt]`
Two off-line pairs at depth `a = 2π`, normalised ordinate separation `0.7`, embedded
in on-line points. Measured: total slack `0.421 × (κ−1)²` — **less than one isolated
pair** — so the second pair's presence *reduces* the total budget charge. Each pair
individually is "detected" at `0.21 × (κ−1)²`, a factor `4.8` below its isolated
value. Extending to 16 pairs at separation `0.5–0.7` gives total
`0.589–0.657 × (κ−1)²`, i.e. per-pair suppression by a factor `24–27`.
**A hypothetical off-line zero can be placed where the framework charges almost
nothing for it.**

### C4 — a global second moment gives no individual lower bound
`[formal hypotheses satisfied; results/additivity.txt]`
`g` off-line pairs of equal depth `a = 2π` at normalised separation **2** — twice the
mean spacing, i.e. a configuration with **no local density excess at all**. Measured
total slack, in units of one isolated pair:

| `g` | 1 | 2 | 4 | 8 | 16 | 32 |
|---|---|---|---|---|---|---|
| total | 1.000 | 1.864 | 2.871 | 3.801 | 4.583 | 5.623 |
| per pair | 1.000 | 0.932 | 0.718 | 0.475 | 0.286 | 0.176 |

`total = 0.9799 + 1.3315·ln g`, max residual `0.0880`; the best linear fit has max
residual `1.0402`. **The global functional grows logarithmically while the population
grows geometrically**, so no individual lower bound follows from the global bound.
This is the configuration class that decides the whole investigation (04 §3).

## 3. The one local arithmetic restriction, and why it does not save the route

The only local constraint zeta supplies is the counting one. From
`N(t) = (t/2π)log(t/2π) − t/2π + 7/8 + S(t) + O(1/t)` with the unconditional
`S(t) ≪ log t`, a normalised ordinate window of length `H` contains `H + O(log T)`
zeros. Hence:

* a cluster of `g` off-line **pairs** (i.e. `2g` zeros) packed into normalised length
  `s·g` has local density `2/s`, so the excess over the mean is `(2/s − 1)·s·g =
  (2 − s)g` zeros, which forces `g = O(log T)` whenever `s < 2`;
* **at `s = 2` the excess is exactly 0** — and that is precisely configuration C4,
  the log-additive one.

So the counting restriction excludes the *dense* clusters (`s < 2`) and permits
exactly the configuration that kills the method. This is not a coincidence: `s = 2`
is where the density constraint becomes vacuous, and the measured interference is
still strong there because the `g`-direction correlation at separation 1–2 is
`−0.75` to `+0.43` at depth `a = 2π` (04 §3.2).

**Conclusion.** ARTZ supplies no arithmetic restriction specific to the actual zeros
that eliminates the adversarial possibilities. Its content is a positive-proportion
statement and it stays one.

## 4. Scope discipline, applied to every bound in this package

| bound | every individual zero | positive proportion | averaged functional | selected finite configuration | hypothetical isolated off-line zero |
|---|---|---|---|---|---|
| `tr Q ≤ (C_MT−1+o(1))N(T)` | — | **yes** | — | — | — |
| ARTZ (80) with `K_τ(c)²\|R̃_1\|` | — | **yes** | — | — | — |
| Theorem 1, `w_ℓ/N ≥ …` | — | **yes** | — | — | — |
| `Σ_j(q_j−2)²` budget (new, 04 §1) | — | — | **yes** | — | — |
| `2m + (m(1+κ)−2)² ≤ …` (04 §3 B1) | — | — | — | — | **yes**, and vacuous |
| `λ_min(G_d) ≤ √(πm)(2sin²(πc/2))^m` (05) | — | — | — | **yes** | — |
| anti-localisation (05 Prop 1) | **yes** (a structural exclusion) | — | — | — | — |
| `Σ_l(κ_l−1)² ≤ …` | — | would be | — | — | — | **not available** |

The last row is the one the transfer needed, and it is the one that does not exist.

## 5. Why the Density-Hypothesis computation matters anyway

04 §3 shows that the summed form of the new inequality is *equivalent in strength* to
the Density Hypothesis for every `σ ∈ (1/2,1)`. That is worth recording for two
reasons.

1. **It explains the literature.** Alpöge–Furman, Lamzouri and ARTZ all extract the
   **trace** of the signed defect (exactly additive, `verify.py` E2) and never the
   Frobenius slack. The reason is now explicit: the slack is the `δ`-sensitive part,
   and making it additive would prove a famous open problem. Nothing has been
   overlooked.
2. **It is a falsification test that the route passed.** Before measuring anything, I
   had a derivation whose conclusion was too strong to be true. The measurement in
   04 §3.1 located the failing step (additivity), quantified it
   (`0.98 + 1.33 ln g`), and identified its mechanism (alignment of the
   `g`-directions at depth). **The obstruction is not a missing lemma; it is a
   measured property of the construction.**
