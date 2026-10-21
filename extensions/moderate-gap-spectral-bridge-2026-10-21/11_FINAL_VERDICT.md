# 11 — Final verdict

## Overall: **E — the proposed transfer is obstructed**, with one **C**-level new
## inequality and one **D**-level exact correspondence preserved.

Per-route verdicts, as the brief requires:

| route / question | verdict | one-line reason |
|---|---|---|
| **Q1** exact Gram ↔ Weil correspondence | **D** | exists exactly, but it is Alpöge–Furman's own construction; preserves neither positivity, moment order, arithmetic input, nor support scale |
| **Q2** moderate gaps → detector conditioning | **E** | off-line zeros have internal ordinate gap **0**, so they are structurally excluded from every certified run; and `c ≤ 1/4` is 4× sub-Nyquist, so `λ_min` decays exponentially in the block length |
| **Q3** `δ²` sensitivity into the spectral budget | **C** | the free Lemma-5 slack `Σ_j(q_j−2)²` is genuinely `δ`-sensitive, with `(κ−1)² = 4V_μ²δ⁴+O(δ⁶)`; new, unconditional, costless |
| **Q4** any new constraint on an individual off-line zero | **E** | the one pointwise bound is vacuous by exactly `0.2084930…·log(T/2π)`; and the summed form is log-additive, not additive |
| **Q5** a standalone spectral theorem worth preserving | **yes, three** | 05 Prop 1, 05 Prop 2, and the non-additivity law of 04 §3 |
| Route A (direct transfer of `Q`) | **D** | exact identity, provably `δ`-blind in trace and inertia |
| Route B (weighted defect) | **C**, summation step **E** | `D(ρ) = (m(1+κ_ρ)−2)²` exists and is rigorous for one isolated pair; `(H-add)` measured false |
| Route C (detector penalty) | **E** | refuted by an explicit configuration with ordinary local spacing |
| Route D (spectral surplus) | **E** | surplus is extensive, single-zero gain is sub-extensive by one `log T` |
| T1 … T7 | D, E, C, E, E, E, E | see 10 |

**Verdict A is not awarded.** No noncircular new constraint on an individual
hypothetical off-critical-line zero was obtained. **Verdict B is not awarded**: the
new estimate does not improve the existing RH machinery (it cannot, at fixed
support — 08 §4).

## The strongest justified mathematical results

**R1 (new inequality, class IV, proved).** For every fixed `τ ∈ (0,1/200]`,
`0 < c ≤ 1/4`,
```
   tr Q  +  K_tau(c)^2 |Rtilde_1|  +  sum_{q_j > 0} ( q_j - 2 )^2
            <=  ||G||_F^2 - N(T)  =  ( C_{eta_tau} - 1 + o(1) ) N(T) ,
```
a strict strengthening of Angón–Ross–Tang–Zhao (80) obtained at **zero cost** by
retaining the slack their Lemma 5 discards at `q² ≥ 4q − 4`. For a single simple
off-line pair of depth `δ = β − 1/2`,
```
   sum_j (q_j - 2)^2 = ( kappa(delta) - 1 )^2 ,
   kappa(delta) = int eta_tau^2(u) cosh( 2 delta u log T ) du = 1 + 2 delta^2 V_mu + O(delta^4),
```
`V_μ` the centred second moment of the window measure. **This is the only place in
the whole Alpöge–Furman / Lamzouri / ARTZ framework where the horizontal displacement
of a zero appears**; the trace and the inertia — the two invariants those papers use —
are provably `δ`-blind (`tr = 2m` for every `κ`).

**R2 (anti-localisation, class IV, proved).** Every off-critical-line zero has a
partner with the *same ordinate*. Hence no off-line zero is an interior point of a
run of positive gaps; a run of `ℓ` consecutive `c`-moderate gaps certifies `ℓ−1`,
not `ℓ+1`, simple critical-line zeros; and the moderate-gap theorem certifies
separation exactly at the ordinates where there is provably nothing to detect. **No
gap theorem can ever be localised at a hypothetical off-line zero**, because the gap
it would have to enlarge is identically zero.

**R3 (sub-Nyquist conditioning ceiling, class IV, certified).** For nodes with
minimum separation `c` and a kernel band-limited to `(−1/2,1/2)`,
```
   lambda_min(G_d)  <=  ( 2 sin(pi c/2) )^{2(d-1)} / C(2d-2, d-1)
                     ~  sqrt(pi d) ( 2 sin^2(pi c/2) )^{d-1} ,
```
with rate `0.2928932…` at `c = 1/4`. The theorem's entire range `0 < c ≤ 1/4` lies
in the exponentially ill-conditioned regime. Conditioning gets **worse**, not better,
as the separated block lengthens.

**R4 (the quantitative no-go, class IV, certified).** The pointwise deficit is
exactly one power of `log T`:
```
   budget / ( maximal single-pair gain )  =  (2/pi)( C_MT - 1 ) log(T/2 pi)
                                          =  0.2084930... * log(T/2 pi) ,
```
`> 1` for `T > 761`, growing without bound, and unimprovable by any choice of window
(the bound `κ ≤ cosh(δ log T)` uses only `η_τ² ≥ 0`, `∫η_τ² = 1`, support in
`(−1/2,1/2)`).

**R5 (non-additivity, class V with a proved mechanism).** The `δ`-sensitive
functional is not additive over off-line zeros. At depth `a = δ log T = 2π` and
normalised ordinate separation 2 — **twice the mean spacing, hence with no local
density excess** — the total slack for `g` pairs is
```
   total / ( kappa - 1 )^2  =  0.9799 + 1.3315 ln g ,   max residual 0.0880 (g = 2..32),
```
against a best linear fit with max residual `1.0402`. Mechanism: as the depth grows,
`η_τ²(u)cosh(2π y u)` concentrates at `|u| → 1/2`, so the `g`-directions of distinct
off-line pairs **align** (normalised correlation `−0.751` at separation 1, depth
`a = 2π`; `−0.884` at `a = 3π`), and a sum of squares of `g` aligned contributions is
`O(1)·(largest)²`, not `g·(each)²`.

**R6 (why nothing was overlooked, class IV, certified).** If `R5` were false — if the
slack were additive with any positive constant — then `R1` would give
`N(σ,T) ≪_ε T^{2(1−σ)+ε}` for every `σ ∈ (1/2,1)`: the **Density Hypothesis**, open
for `σ ∈ (1/2, 25/32)`, with an exponent strictly better than Selberg 1946 and
Ingham throughout `(1/2, 3/4]`. The non-additivity is therefore *forced*, and it
explains why Alpöge–Furman, Lamzouri and ARTZ all extract the trace (exactly
additive) and never the Frobenius slack.

**R7 (cross-framework confirmation of the `δ⁴` exponent, class IV).** The budget cost
of an off-line pair is `(κ−1)² = 4V_μ²δ⁴ + O(δ⁶)`. My own program independently
measured and certified the off-line discrimination strength as `M(δ) = δ⁴·D₂`,
`D₂ = 2911.050006543213352847841125028374273138 ± 1e−37`. Two constructions sharing
no machinery give the same exponent 4, for the same reason: the displacement enters
through the even function `cosh(δ·)`, so the *response* is `O(δ²)`, and the budgeted
quantity is the square of the response. The `δ²` coefficient is, in both frameworks,
**twice the centred second moment of the window** — and the apparent `L³` of my
published law is the `s = 1` (unnormalised-window) member of a family whose general
power is `L^{s+2}`; ARTZ's probability window is the `s = 0`, i.e. `L²`, member.

## Three highest-value next steps

1. **Obtain BGSTB, Acta Arith. 214 (2024) 357–376, and check Lemma 5's hypotheses
   against the `s = 2` configurations of 09 §2 C4.** This is the single remaining
   dependency (D1 of §01). Everything in this package was derived *from* that lemma
   as restated by Lamzouri and Goldston–Sahay; I did not see the original. The
   specific question: whether the error term `O_f(1/√log T)·N(T)` and the implied
   constant are uniform over configurations whose symmetric-diagonal terms
   `Σ_ρ T^{2δ_ρ}` are large. If they are, `R1` and `R6` stand as written; if the
   lemma carries an additional restriction on `β`, `R6` becomes moot and `R1` applies
   only in the narrow-strip regime `a = O(1)`. **Cheap, decisive, and the only item
   that could change a stated result.**
2. **Settle the additivity exponent of `Σ_j(q_j−2)²` rigorously.** `R5` is a
   measurement on `g ≤ 32` pairs at a single depth with one admissible window. The
   question — *what is the exact growth order of
   `sup{Σ_j(q_j−2)² : g pairs of depth a, normalised separation ≥ s}` in `g`?* — is a
   self-contained problem in the harmonic analysis of the kernel `K_τ`, independent of
   RH, and a clean `Θ(log g)` theorem would convert `R5` from class V to class IV and
   close Route B permanently. The alignment identity of 04 §3.2 is the right starting
   point, and it is the same mechanism as
   `weil-kappa-is-alignment-not-conditioning`'s `κ = β(g*)² + β'(g*)²d/(2κ₂)`.
3. **Finish the Hermite–Genocchi lower bound for `κ(Γ)` and then stop.** It is
   dependency D3 and it is the only thing standing between 06 §1 (class VI) and a
   genuine class-IV statement. But note the ceiling: even completed, it certifies
   `κ` only at on-line blocks, which 05 Proposition 1 shows can never contain the
   zero being tested. So its value is as a *standalone* interpolation theorem — the
   confluent note at `~/scratch/rh-confluent-zenodo-note-2026-10-15-141204` — and
   **not** as RH progress. The 2026-10-15 no-go continues to govern.

## Not done, deliberately

* No publication manuscript was written (the brief forbids it without approval).
* No repository, branch, commit, Zenodo record, or `CITATION.cff` was touched;
  `~/GitHub/weil-window-control` remains at `725f204` with only its pre-existing
  untracked build artefacts.
* No Lean build was run and no formal verification is claimed, for either my work or
  the Lean certificates reported by Alpöge–Furman and by Lamzouri's Appendix A.
* No researcher was contacted.
* No statement here depends on the location of any actual zero of `ζ`; every
  computation is on synthetic spectra or is exact arithmetic on the paper's own
  constants.
