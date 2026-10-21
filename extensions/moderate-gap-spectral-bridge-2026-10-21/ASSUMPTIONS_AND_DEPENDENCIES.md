# Assumptions and unresolved dependencies

Listed in decreasing order of how much a stated result depends on them.

## D1 — BGSTB Lemma 5 is taken as restated, not read in the original `[the one that matters]`

The sole analytic input is the unconditional second-moment asymptotic

```
   sum_{rho,rho' : 0<gamma,gamma'<=T} K_tau( i(rho-rho') log T / (2 pi) )^2
        =  ( C_{eta_tau} + O_tau( 1/sqrt(log T) ) ) N(T) ,
```

quoted by ARTZ as their Lemma 1 / (12), by Lamzouri as his Lemma 3.2 (resting on his
Lemma 3.1 = "Lemma 5 of [BGSTB]"), and underlain by the unconditional Montgomery
theorem in the form given by Goldston–Sahay §3:
`F(x,T) = Σ_{ρ,ρ'} x^{ρ−ρ'}W(ρ−ρ')`, `W(u) = 4/(4−u²)`, `F ≥ 0`, `F(x)=F(1/x)`,
`F(x,T) = (T/2πx²)log²T(1+O(1/√log T)) + (T/2π)log x + O(T√log T)` for `1 ≤ x ≤ T`.

I did **not** obtain Acta Arith. 214 (2024) 357–376. I verified that Lemma 3.1 and
the `F(x,T)` asymptotic are mutually consistent and that the `f(0) + 2∫₀¹αf(α)dα`
shape emerges from `∫_{-1}^{1}f(u)F(T^u,T)du` as it should (02 §5), and I recomputed
`C_MT` and `C_{η}` independently (B8, A1). That is as far as I went.

**Why it matters.** `R1` and `R6` of `11_FINAL_VERDICT.md` are derived *from* this
statement. In particular `R6` observes that the symmetric-diagonal terms of that sum
are `Σ_ρ T^{2δ_ρ}·W(2δ_ρ)`, which grow like `T^{2δ}` per off-line pair, so the
asymptotic has density-type content whose extraction is blocked only by the
non-positivity of the `γ ≠ γ'` terms — exactly the obstruction Goldston–Sahay name
("a positivity that fails for zeros far from the line"). If BGSTB's lemma carries an
additional restriction on `β` that the restatements omit, then `R6` becomes moot and
`R1` applies only in the narrow-strip regime `a = δ log T = O(1)`. **`R2`, `R3`, `R4`,
`R5` and `R7` are unaffected either way**, since they are linear algebra, harmonic
analysis, or measurement.

## D2 — Montgomery–Taylor extremality is cited, not re-derived

`C_MT = 1/2 + (1/√2)cot(1/√2)` is recomputed to width `3.3e−40` here, and
`C_{η} = (A∗A)(0) + 2∫₀¹α(A∗A)(α)dα` is checked against it by quadrature to
`1.7e−7`. But the *extremality* of the Montgomery–Taylor window among admissible
windows is taken from Montgomery (1975) and Carneiro–Chandee–Littmann–Milinovich
(2017, Cor. 14) as ARTZ do. Nothing in this package would change if a better window
existed — `R4`'s bound `κ ≤ cosh(δ log T)` is window-independent.

## D3 — the confluent Hermite–Genocchi lower bound is still unfinished

Noted as "remaining" in `weil-confluent-cluster-theorem-audited`. Consequently every
statement of the form "moderate gaps give `κ(Γ) ≥ …`" is **class VI**, not IV
(06 §1). The upper bound
`κ ≤ (∏h_i)²·2L^{2n+1}/((2n+1)(n!)²)` is available and audited GREEN; the exponent
`2(n−1)` is Batenkov-known and must never be headlined.

## D4 — no Lean work

Alpöge–Furman state that their results are formally verified in Lean 4 and Lamzouri
records Lean certificates produced on short notice by the Axiom Math team. I did not
build, run, or audit either, and I claim no Lean verification for anything in this
package.

## D5 — no zeta zero data, and no numerical library

`numpy`, `sympy`, `mpmath` and `pip` are unavailable on this machine. All arithmetic
is Python stdlib: exact `Fraction` interval arithmetic for sign-critical statements
(`scripts/iv.py`, with Taylor series and explicit tail bounds), and double precision
elsewhere, always labelled `[diagnostic]`. Eigenvalues are computed by a hand-rolled
cyclic Jacobi routine (`scripts/linalg.py`) and are used only for `[diagnostic]`
statements; no sign-critical conclusion rests on a floating-point eigenvalue.

## D6 — the non-additivity law is a measurement, not a theorem

`R5` (`total/ideal = 0.9799 + 1.3315 ln g`, max residual `0.0880` for `g = 2..32`)
is established at **one depth** (`a = 2π`), **one separation** (`2`), **one
admissible window** (`τ = 1/200`), and `g ≤ 32`. The mechanism (alignment of the
`g`-directions) is proved and the correlation table is measured across five depths,
but the growth *order* in `g` is not proved. Next step 2 of `11_FINAL_VERDICT.md`
names the exact problem. **Do not quote `1.33` as a constant of nature.**

## D7 — the `η² ≡ 1` (sinc) window is not admissible

`scripts/cluster_scan.py` and blocks D–F of `verify.py` use `η² ≡ 1` on `(−1/2,1/2)`,
for which `K` has closed forms at complex arguments. This window is **not** in ARTZ's
admissible class (`η` must be `C_c^∞` and vanish at `±1/2`). It is used only for
*structural* statements — the Gram identities, the inertia, the slack identity, the
chain of 04 §1 — which hold for any nonnegative even window of unit mass. Every
*quantitative* claim about clusters is repeated with the genuinely admissible window
in `scripts/cluster_smooth.py` and `scripts/additivity.py`, and the two agree
qualitatively (saturation at sub-mean separation, growth at `sep = 1`).

## D8 — literature not read, used only as attribution

Bourin–Lee (2012) trace superadditivity; Bhatia, *Matrix Analysis*; Gonek–Sahay
arXiv:2412.15481; Aryan; Montgomery (1975); Montgomery–Taylor; Titchmarsh;
Carneiro–Chandee–Littmann–Milinovich (2017). The convex trace superadditivity
`trΨ(A+B) ≥ trΨ(A) + trΨ(B)` for `A,B ⪰ 0` and `Ψ` convex nondecreasing with
`Ψ(0)=0` is used as ARTZ use it, and is **not** re-proved here. Selberg's 1946
density theorem (`N(σ,T) ≪ T^{1−(σ−1/2)/4}log T`) and Ingham's
(`N(σ,T) ≪ T^{3(1−σ)/(2−σ)}log⁵T`) are used only for the comparison table of 04 §3
and were taken from a literature search, not from the originals; the Guth–Maynard
range for the Density Hypothesis (`σ ≥ 25/32`) likewise.
