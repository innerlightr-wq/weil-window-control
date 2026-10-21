# RESULT — moderate zero gaps vs finite-window Weil positivity

**Question.** Do the techniques of Angón–Ross–Tang–Zhao, *Unconditional Bounds for Two,
Three and Four Consecutive Moderate Gaps between Zeros of the Riemann Zeta Function*
(arXiv:2610.22927v1, 16 Oct 2026) yield new results for the finite-window Weil-positivity
program of this repository?

**Verdict: NO USABLE TRANSFER.** One exact structural correspondence (already the
antecedents' own construction), one new but unsummable inequality, and a quantified
obstruction. Verification: **68/68 checks**, every sign-critical one in exact-rational
interval arithmetic, Python standard library only (`results/verify.txt`).

## What was obtained

**(R1) A new, costless strengthening of the paper's budget inequality.** `[T1]`
Their Lemma 5 discards the slack in `q² ≥ 4q − 4`. Retaining it gives, for every fixed
`τ ∈ (0,1/200]` and `0 < c ≤ 1/4`,

```
   tr Q  +  K_tau(c)^2 |Rtilde_1|  +  sum_{q_j > 0} ( q_j - 2 )^2
            <=  ||G||_F^2 - N(T)  =  ( C_{eta_tau} - 1 + o(1) ) N(T) ,
```

and for a single simple off-line pair of depth `δ = β − 1/2` the boxed term equals
`(κ(δ) − 1)²` with

```
   kappa(delta) = int eta_tau^2(u) cosh( 2 delta u log T ) du = 1 + 2 delta^2 V_mu + O(delta^4) .
```

**`κ` is the only place in the whole Alpöge–Furman / Lamzouri / ARTZ framework where a
zero's horizontal displacement appears.** The two invariants those papers use are provably
the `δ`-blind ones: the off-line block has spectrum `{m(1+κ), −m(κ−1)}`, so `tr = 2m` for
every `κ`, while `‖·‖²_F = 2m²(1+κ²)` grows.

`V_μ` is the **centred second moment of the window measure** — the same constant, and the
same description, as the coefficient of this repository's `δ²` law (N1). The apparent `L³`
there is the `s = 1` member of a family whose general power is `L^{s+2}` for window mass
`~L^s`; ARTZ's probability window is the `s = 0`, i.e. `L²`, member.

**(R2) Anti-localisation.** `[T1, proved]`
Every off-critical-line zero `ρ` has the partner `1 − ρ̄` at **ordinate distance exactly 0**.
Hence no off-line zero is an interior point of a run of positive gaps; a run of `ℓ`
consecutive `c`-moderate gaps certifies `ℓ − 1`, not `ℓ + 1`, simple critical-line zeros;
and ARTZ's Proposition 1 builds its runs precisely by discarding the clusters containing
the off-line and multiple zeros. **The moderate-gap theorem certifies separation only at
ordinates where there is provably nothing to detect, and no gap theorem can ever change
that**, because the gap it would have to enlarge is identically zero.

**(R3) A sub-Nyquist conditioning ceiling.** `[T1, certified]`
For nodes of minimum separation `c` and a kernel band-limited to `(−1/2,1/2)` (Nyquist
spacing 1), the binomial test vector gives

```
   lambda_min(G_d)  <=  ( 2 sin(pi c/2) )^{2(d-1)} / C(2d-2, d-1)
                     ~  sqrt(pi d) ( 2 sin^2(pi c/2) )^{d-1} ,
```

with rate `0.2928932…` at `c = 1/4`. The theorem's entire range `0 < c ≤ 1/4` is 4×
oversampled, so conditioning gets **worse**, exponentially, as a separated block lengthens.

**(R4) The pointwise deficit is exactly one logarithm.** `[T1, certified]`
Since `β ≤ 1` forces `δ ≤ 1/2` and `κ ≤ cosh(δ log T)` for *every* admissible window,

```
   budget / ( maximal single-pair gain )  =  (2/pi)( C_MT - 1 ) log(T/2 pi)
                                          =  0.2084930... * log(T/2 pi) ,
```

`> 1` for `T > 761` and growing. The budget is extensive (one unit per zero); a single
zero's gain is at most `T^{2δ} ≤ T`, sub-extensive by one logarithm. No choice of window
changes this.

**(R5) The `δ`-sensitive functional is not additive — measured.** `[T2]`
At depth `a = δ log T = 2π` and normalised ordinate separation **2** — twice the mean
spacing, hence no local excess over the Riemann–von Mangoldt density — the total slack for
`g` off-line pairs is

```
   total / (kappa - 1)^2  =  0.9799 + 1.3315 ln g ,   max residual 0.0880  (g = 2..32),
```

against a best linear fit with max residual `1.0402`. Mechanism, proved: as the depth
grows, `η_τ²(u)cosh(2πyu)` concentrates at `|u| → 1/2`, so the `g`-directions of distinct
off-line pairs **align** (normalised correlation `−0.751` at separation 1, `−0.884` at
`a = 3π`), and a sum of squares of `g` aligned contributions is `O(1)·(largest)²`.

**(R6) Why nothing was overlooked.** `[T1, certified]`
Had the slack been additive with any positive constant, (R1) would give
`N(σ,T) ≪_ε T^{2(1−σ)+ε}` for every `σ ∈ (1/2,1)` — the **Density Hypothesis**, open for
`σ ∈ (1/2, 25/32)`, with an exponent strictly better than Selberg (1946) and Ingham
throughout `(1/2, 3/4]`. Non-additivity is therefore forced, and this is exactly why
Alpöge–Furman, Lamzouri and ARTZ all extract the trace (exactly additive) and never the
Frobenius slack.

**(R7) Independent confirmation of the `δ⁴` exponent.** `[T1]`
The budget cost of an off-line pair is `(κ−1)² = 4V_μ²δ⁴ + O(δ⁶)`. This repository's
`weil-offline-detect` study independently certified the off-line discrimination strength as
`M(δ) = δ⁴·D₂`, `D₂ = 2911.050006543213352847841125028374273138 ± 1e−37`. Two
constructions sharing no machinery give the same exponent, for the same reason: the
displacement enters through the even function `cosh(δ·)`, so the response is `O(δ²)` and
the budgeted quantity is its square.

## Non-claims, for this extension specifically

1. **No proof or advance of RH**, and nothing here assumes RH. No statement depends on the
   location of any actual zero of ζ; every computation is exact arithmetic on the paper's
   own constants or is run on *synthetic* conjugation-invariant multisets.
2. **The Density Hypothesis is not claimed.** (R6) is a *falsification test* that the route
   passed: a derivation too strong to be true, whose failing step was then located,
   quantified and explained. The conditional statement in `10_CANDIDATE_THEOREMS.md` (T6)
   has a hypothesis that is both measured false and DH-equivalent.
3. **(R5) is a measurement, not a theorem** — one depth, one separation, one admissible
   window, `g ≤ 32`. The growth order in `g` is open (`11_FINAL_VERDICT.md`, next step 2).
4. **Dependency not closed:** BGSTB, *Acta Arith.* 214 (2024) 357–376 was **not obtained**;
   its Lemma 5 is used as restated by Lamzouri and Goldston–Sahay. (R1) and (R6) rest on
   that restatement. See `ASSUMPTIONS_AND_DEPENDENCIES.md` D1 — it is the only open item
   that could change a stated result.
5. **The existence of the Gram ↔ Weil correspondence is not new**: it is Alpöge–Furman's own
   stated mechanism ("a finite compression of Weil's Hermitian form, with Sylvester's law of
   inertia handling off-line pairs"). What is new is quantitative — `κ`, the Gram table, and
   the shared `δ²` constant.
6. **No Lean.** No build was run, for this work or for the Lean certificates reported by
   Alpöge–Furman and by Lamzouri's Appendix A.
7. **No `mpmath`, `numpy` or `sympy`** on the machine used; every sign-critical conclusion is
   therefore in exact-rational interval arithmetic (`scripts/iv.py`), and everything
   double-precision is labelled `[diagnostic]`.

## Verdicts, per route

| route | verdict |
|---|---|
| exact Gram ↔ Weil correspondence | **D** — exact, but already the antecedents' construction; preserves neither positivity, moment order, arithmetic input, nor support scale |
| moderate gaps → detector conditioning | **E** — anti-localisation (R2) and sub-Nyquist decay (R3), independently |
| `δ²` sensitivity into the spectral budget | **C** — (R1) is new, unconditional and costless |
| any new constraint on an individual off-line zero | **E** — (R4) and (R5) |
| a standalone estimate worth keeping | **yes** — (R2), (R3), (R5) |

Full detail: `README.md` for the layout, `11_FINAL_VERDICT.md` for the verdicts and the
three next steps, `THEOREM_STATUS.md` for every statement with its provenance class and
verification tier.
