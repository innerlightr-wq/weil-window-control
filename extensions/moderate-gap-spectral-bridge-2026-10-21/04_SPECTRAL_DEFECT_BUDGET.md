# 04 — The signed spectral defect budget, and a displacement-sensitive refinement

Answers the brief's §5 (Routes A–D) and Final Question 3.

## 1. Independent reconstruction of the budget inequality

Conventions are ARTZ's: `Z_T = {i(ρ−1/2)log T/(2π)}` with multiplicity,
`G = P + Q`, `P` from the simple critical zeros, `R̃_1` the simple critical zeros
lying within `c` of another one, `n_+(Q) = r + k`, `tr Q = N(T) − n`.

Chain, rederived here line by line (each step re-checked, `verify.py` C1–C5, E1–E4):

```
 (i)   ||P+Q||_F^2 = ||P-Q_-||_F^2 + 2 tr((P-Q_-) Q_+) + ||Q_+||_F^2 ,  tr(Q_- Q_+) = 0
 (ii)  tr(P Q_+) >= 0                                     [P, Q_+ both PSD]
 (iii) ||Q_+||_F^2 = sum_j q_j^2 = 4 tr Q_+ - 4 n_+(Q) + SUM_j (q_j - 2)^2     <-- identity
 (iv)  ||P-Q_-||_F^2 >= sum_j (p_j - qtilde_j)^2           [von Neumann]
 (v)   sum_j (p_j - qtilde_j)^2 >= trPsi(P) - 4 tr Q_-     [ARTZ Lemma 4]
 ---------------------------------------------------------------------------
 (vi)  ||G||_F^2 >= trPsi(P) + 4 tr Q - 4 n_+(Q) + SUM_j (q_j - 2)^2
 (vii) tr Q >= 2 n_+(Q)   =>   4 tr Q - 4 n_+(Q) >= 2 tr Q
 (viii) trPsi(P) >= n + K_tau(c)^2 |Rtilde_1|               [ARTZ Lemma 12]
 ---------------------------------------------------------------------------
 (ix)  tr Q + K_tau(c)^2 |Rtilde_1| + SUM_j (q_j - 2)^2  <=  ||G||_F^2 - N(T)
```

Line (ix) is ARTZ (80) **plus the boxed term**. ARTZ's Lemma 5 drops
`Σ_j(q_j−2)²` at step (iii) by using `q² ≥ 4q−4`; the term is pure slack and is
recovered at **zero cost**. This is the only new inequality in this package.

```
   NEW (class IV, proved):
        tr Q  +  K_tau(c)^2 |Rtilde_1|  +  sum_{q_j > 0} (q_j - 2)^2
                 <=  ||G||_F^2 - N(T)  =  ( C_{eta_tau} - 1 + o(1) ) N(T) .
```

Checked on seven synthetic conjugation-invariant multisets including adversarially
placed on-line zeros and clustered off-line pairs: holds in all seven
(`verify.py` F4) `[diagnostic]`; each step of (i)–(ix) is `[certified]` or `[proved]`.

### Normalisation audit
* `tr G = N(T)` exactly (ARTZ (27)); reproduced to `1e-6` relative on all seven
  configurations (`verify.py` F3).
* `‖G‖²_F = Σ_{z,s∈Z_T}K_τ(z−s)²` with the **bilinear** pairing (03 §2 Q8);
  reproduced to `1e-6` relative (`verify.py` F1).
* The available asymptotic excess is `(C_{η_τ} − 1)N(T)`, with
  `C_MT − 1 = 0.3274992963205883…` `[certified A1]`. In absolute terms
  `(C_MT−1)N(T) = 0.3275·(T/2π)log(T/2π)(1+o(1)) ≈ 0.05212·T·log(T/2π)`.
* ARTZ's `tr(Q) + K_δ(c)²|R_close| ≤ ‖G‖²_F − N(T)` as quoted in the brief is
  exactly (80); the brief's `|R_close|` is ARTZ's `|R̃_1|` — **simple critical zeros
  in close pairs**, not off-line zeros. The off-line zeros sit in `tr Q`.

## 2. Route A — direct transfer of `Q` to a Weil defect

**Verdict: exact identity available, but it is `δ`-blind.** `[IV, certified]`

For one off-line conjugate pair of multiplicity `m`, `Q`'s contribution is
`2m(v_g v_g^T − v_h v_h^T)` with `⟨g,h⟩ = 0`, `‖g‖² = (1+κ)/2`, `‖h‖² = (κ−1)/2`.
Hence its nonzero spectrum is `{m(1+κ), −m(κ−1)}` and

```
    tr          = 2 m                                  (independent of kappa)
    ||.||_F^2   = 2 m^2 (1 + kappa^2)                  (grows like kappa^2)
    inertia     = (1,1)                                (independent of kappa)
```
`[certified E1, E2]`. So **the trace and the inertia — the only two features ARTZ,
Lamzouri and Alpöge–Furman use — are exactly the `δ`-blind ones**, and the
Frobenius norm is the only `δ`-sensitive one. This is the precise sense in which the
brief's warning ("do not infer horizontal sensitivity from a rank or inertia count
alone") is correct: it is not a caution, it is a theorem.

## 3. Route B — a weighted defect

**Verdict: a nonnegative `δ`-sensitive quantity exists and the inequality holds; it
is rigorous for a single isolated pair and FAILS to be additive in general.**

Define, for an off-line pair `ρ` of multiplicity `m_ρ`,
```
    D(rho) := ( m_rho (1 + kappa_rho) - 2 )^2 ,
    kappa_rho = int eta_tau^2(u) cosh( 2 delta_rho u log T ) du  >= 1 ,
    so  D(rho) = (kappa_rho - 1)^2  when m_rho = 1
              = 4 V_mu^2 delta_rho^4 + O(delta_rho^6)  .
```
`D ≥ 0`, `D = 0` iff `δ_ρ = 0` (for `m=1`), and `D` is monotone in `|δ_ρ|`.
`[IV, certified E3, E5, E6]`

**(B1) Single-pair form — rigorous, unconditional.** If exactly one off-line
conjugate pair occurs with `0 < γ ≤ T` and all other zeros there are simple and on
the line, then `Q = 2m(v_gv_g^T − v_hv_h^T)` and (ix) gives
```
    2 m + ( m(1+kappa) - 2 )^2  <=  ( C_{eta_tau} - 1 + o(1) ) N(T) .
```
No isolation hypothesis is needed for this, because `Q` excludes the simple critical
zeros by construction. **Scope: a hypothetical isolated off-line zero.** `[IV, proved]`

**(B2) The summed form — requires additivity, and additivity is FALSE.**
What one would want is
```
    (H-add)      sum_j (q_j - 2)^2  >=  c_0  sum_{off-line pairs} (kappa_l - 1)^2
```
for some absolute `c_0 > 0`. Under `(H-add)` one gets at once
```
    sum_{off-line pairs, gamma<=T} (kappa_l - 1)^2  <=  c_0^{-1}(C_MT - 1 + o(1)) N(T) ,
```
and since `κ ≥ ∫_{|u|≤1/2−τ}η_τ²cosh(2δu log T)du ≫_τ T^{δ(1−2τ)}/(δ log T)`, this
would give, for every fixed `σ ∈ (1/2,1)` and every `ε>0`,
```
    N(sigma, T)  <<_eps  T^{2(1-sigma)+eps}  --- the DENSITY HYPOTHESIS.
```
The exponent identity `1 − 2δ = 2(1−σ)` is exact `[certified G6]`, and the resulting
exponent beats Selberg 1946 (`1 − (σ−1/2)/4`) and Ingham (`3(1−σ)/(2−σ)`) at every
`σ ∈ (1/2, 3/4]` `[certified G7]`:

| σ | slack route | Density Hyp. | Selberg 1946 | Ingham |
|---|---|---|---|---|
| 0.51 | 0.980 | 0.980 | 0.9975 | 0.98658 |
| 0.55 | 0.900 | 0.900 | 0.9875 | 0.93103 |
| 0.60 | 0.800 | 0.800 | 0.9750 | 0.85714 |
| 0.70 | 0.600 | 0.600 | 0.9500 | 0.69231 |
| 0.75 | 0.500 | 0.500 | 0.9375 | 0.60000 |

The Density Hypothesis is open for `σ ∈ (1/2, 25/32)`. **Therefore `(H-add)` must
fail, and the only question is how badly.** I measured it.

### 3.1 Measured non-additivity `[V, diagnostic]`

`Σ_j(q_j−2)²` depends only on the Gram of `{g_l, h_l}` over the off-line pairs
(`P` and the on-line zeros do not enter `Q`), so this is self-contained.
Admissible ARTZ window, `τ = 1/200`; `g` pairs of equal depth at normalised ordinate
separation `sep`; `ideal := (κ−1)²` is one isolated pair's cost.
Full tables: `results/cluster_smooth.txt`, `results/cluster_scan.txt`,
`results/additivity.txt`.

At depth `a = δ log T = 2π` (`κ = 37.0587`), total slack in units of `ideal`:

| sep \ g | 1 | 2 | 4 | 8 | 16 | 32 |
|---|---|---|---|---|---|---|
| 0.5 | 1.000 | 0.431 | 0.669 | 0.689 | 0.657 | — |
| 0.7 | 1.000 | 0.421 | 0.760 | 0.537 | 0.589 | — |
| **2.0** | 1.000 | 1.864 | 2.871 | 3.801 | 4.583 | 5.623 |
| 3.0 | 1.000 | 1.686 | 2.641 | 4.057 | 6.519 | — |

Two separate regimes, both fatal:

* **Saturation (sep ≈ 0.5–0.7).** Sixteen off-line pairs cost **0.59–0.66** of a
  single isolated pair. The cluster total does not grow with `g` at all.
  Per-pair cost at `g=16`: `0.0368 × ideal` — a factor **27** below `(H-add)` with
  `c_0 = 1`.
* **Logarithmic growth at mean-spacing separation (sep = 2).** This is the decisive
  case because `sep = 2` is *twice* the mean spacing: such a configuration has **no
  local excess over the Riemann–von Mangoldt density** and therefore violates no
  counting constraint. The fit is
  ```
        total / ideal  =  0.9799 + 1.3315 ln g ,   max residual 0.0880  (g = 2..32)
  ```
  against a best linear fit with max residual `1.0402` — **twelve times worse**.
  Growth is logarithmic, not linear.

Consequence: with `k` off-line pairs of depth `δ` at mean-spacing separation,
(ix) gives only
```
    0.98 + 1.33 ln k   <=   (C_MT - 1 + o(1)) N(T) / (kappa - 1)^2 ,
    i.e.   k  <=  exp( C N(T) / (kappa-1)^2 ) ,
```
which exceeds `N(T)` for every `δ ≤ 1/2` and every `T`. **The refined budget yields
no bound whatever on the number of off-line zeros.** `[IV, proved from V]`

### 3.2 The mechanism: alignment, not conditioning `[IV]`

`g_z(u) = η_τ(u)e^{-2πiuξ}cosh(2π y u)` with `y = δL`. As the depth `a = 2πy` grows,
`η_τ²(u)cosh(2πyu)` concentrates at `|u| → 1/2`, so the functions `g_z` for different
ordinates `ξ` differ only by the phase `e^{-2πiuξ}` *sampled on a shrinking set* —
and therefore become increasingly **collinear**. The normalised correlation
```
    corr(g_z, g_z') = [ Re K(D + 2iy) + Re K(D) ] / (1 + kappa)
```
stays large in magnitude at separations `D` of order the mean spacing, and grows
with the depth `a = δ log T`:

| `a = δ log T` | `κ` | `D = 1` | `D = 2` | `D = 5` |
|---|---|---|---|---|
| 0.628 | 1.0617 | `+0.0397` | `−0.0144` | `+0.0068` |
| 1.571 | 1.4256 | `−0.0557` | `+0.0071` | `+0.0050` |
| 3.142 | 3.4145 | `−0.3303` | `+0.1087` | `−0.0076` |
| 6.283 | 37.059 | **`−0.7508`** | `+0.4342` | `−0.0858` |
| 9.425 | 548.69 | **`−0.8839`** | `+0.6520` | `−0.2020` |

`[V, diagnostic; results/additivity.txt]` — at `D = 1`, exactly the mean spacing,
`|corr| → 1` as the depth grows. Consequently `Q_+` has **one** large eigenvalue plus
`g−1` small ones: the fraction of `Σ(q_j−2)²` carried by the top eigenvalue falls
from `1.000` (`g=1`) to `0.312` (`g=32`) while the *total* only triples. A sum of
squares of `g` aligned contributions is `O(1)·(largest)²`, not `g·(each)²` — hence
log-additivity.

> This is the same phenomenon my own program isolated and named:
> `weil-kappa-is-alignment-not-conditioning` — the close-pair blind spot is an
> **alignment** effect, not a conditioning effect, and no spectral invariant
> (`tr`, `tr^{-1}`, `det`, `λ_min`, `cond`) can see it. It reappears here, in a
> construction that shares none of that machinery, as the exact reason ARTZ,
> Lamzouri and Alpöge–Furman extract only the **trace** (which *is* exactly
> additive, `verify.py` E2) and never the Frobenius slack. `[III + IV]`

**This also answers, in the negative, the question the Density-Hypothesis
computation raises**: there is no unnoticed major corollary sitting in the discarded
slack. The slack is real, new and `δ`-sensitive, but it is sub-additive by exactly
the amount that makes it useless for counting.

## 4. Route C — a combined spacing/displacement penalty

**Verdict: no such inequality survives the near-collision test.** `[IV]`

The desired object would bound detectability and interference together. Every
candidate of the form `Σ_ρ D(ρ)·Φ(local spacing) ≤ budget` with `Φ` bounded below
fails for the same reason as `(H-add)`: §3.1's `sep = 2` configurations have
*entirely ordinary* local spacing and still collapse the left side to `O(log k)`.
The interference is not a near-collision effect — it survives at twice the mean
spacing, where no collision is present. Any `Φ` able to rescue the inequality would
have to vanish at ordinary spacing, i.e. the inequality would say nothing.

## 5. Route D — the spectral surplus, with scopes kept separate

The brief requires each bound to carry its scope. Here they are.

| bound | statement | scope |
|---|---|---|
| S1 | `tr Q ≤ (C_MT−1+o(1))N(T)`, i.e. `≥ 2 − C_MT = 0.6725` of zeros simple and on the line | **positive proportion** (this is AF / Lamzouri Theorem) |
| S2 | `tr Q + K_τ(c)²\|R̃_1\| ≤ (C_MT−1+o(1))N(T)` | **positive proportion** (ARTZ (80)) |
| S3 | `… + Σ_j(q_j−2)² ≤ …` | **averaged spectral functional** (new, §1) |
| S4 | `2m + (m(1+κ)−2)² ≤ (C_MT−1+o(1))N(T)` | **a hypothetical isolated off-line zero** (§3 B1) |
| S5 | `N(σ,T) ≪ T^{2(1−σ)+ε}` | **would be** a positive-proportion/density statement — **not available**, `(H-add)` false (§3) |
| S6 | `λ_min(G_d) ≤ …` under separation `c` | **a selected finite configuration** (see 05) |

**S4 is the only pointwise statement, and it is vacuous.** With `β ≤ 1`, i.e.
`δ ≤ 1/2`, and `κ ≤ cosh(δ log T) ≤ (T^{1/2}+T^{-1/2})/2` for every admissible
window, the maximal left side is `≈ T/4`, while the budget is
`(C_MT−1)(T/2π)log(T/2π)`. The ratio is

```
    budget / (max gain)  =  (2/pi)(C_MT - 1) log(T/2 pi)  =  0.208509... * log(T/2 pi) ,
```
which exceeds 1 as soon as `log(T/2π) > 4.7963`, i.e. `T > 761`:

| `T` | budget `(C_MT−1)N(T)` | max gain `(1+κ)²` | budget / gain |
|---|---|---|---|
| `10³` | `264.2` | `283.2` | `0.933` |
| `10⁶` | `6.243e5` | `2.510e5` | `2.487` |
| `10¹²` | `1.344e12` | `2.500e11` | `5.378` |
| `10³⁰` | `3.505e30` | `2.500e29` | `14.02` |

`[certified G1–G5]`. The known zero-free region (`δ ≤ 1/2 − c/log T`) multiplies the
gain by `e^{-2c}` and leaves the deficit at `≍ log T`. **The method misses a
pointwise exclusion by exactly one power of `log T`**, at the extreme `δ = 1/2`, and
by `T^{1−2δ}log T` for smaller `δ`. The budget is *extensive* (one unit per zero,
`≍ T log T`), the single-zero gain is at most `T^{2δ} ≤ T`: sub-extensive by one
logarithm. The method is intrinsically a counting method, and the deficit is not a
constant that could be improved by a better window.

## 6. Answers

* **Route A:** exact identity, `δ`-blind. The trace and inertia are provably the
  `δ`-blind invariants.
* **Route B:** the quantity `D(ρ) = (m(1+κ_ρ)−2)²` exists, is rigorous for an
  isolated pair, and is **not additive** — measured log-additive at mean-spacing
  separation, with the saturation constant `≈ 1.33` per `e`-fold.
* **Route C:** closed; the interference is not a collision effect.
* **Route D:** no usable surplus. The only pointwise bound is vacuous by
  `0.2085·log(T/2π)`.
* **Final Question 3:** the `δ²` sensitivity *can* be inserted into the budget — that
  is S3 and it is new — but it cannot be summed, so it yields no new arithmetic
  restriction. **Verdict for Candidate T3: C (restricted), with the summation step
  refuted.**
