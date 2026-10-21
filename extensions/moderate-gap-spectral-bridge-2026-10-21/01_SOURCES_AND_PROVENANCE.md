# 01 — Sources, provenance, and theorem status

All work in this package is confined to
`~/scratch/rh-moderate-gap-spectral-bridge-2026-10-21-102437/`.
Nothing outside it was modified. No repository, branch, commit, issue, or Zenodo
record was touched; no researcher was contacted.

## 1. The new paper

| field | value |
|---|---|
| title | *Unconditional Bounds for Two, Three and Four Consecutive Moderate Gaps between Zeros of the Riemann Zeta Function* |
| authors | César Angón, Erick Ross (deceased), Ze Sen Tang, Boyang Zhao |
| id | arXiv:2610.22927v1 [math.NT], 16 Oct 2026, 21 pages |
| local file | read read-only from `~/Downloads/`; **not redistributed here** (see `sources/SOURCE_HASHES.txt`) |
| **sha256** | `9ac30e29221ff9832f303e817ecb3ff54b33bb1a41e9d1159b5b29e963e38fb5` |
| text | extracted locally with `pdftotext -layout`; not redistributed |

Referred to below as **ARTZ**. The paper carries an explicit "AI assistance"
statement (ChatGPT GPT-5.6 Sol / GPT-6 Astra used for exploration, intermediate
calculations, and revision; authors take responsibility). Recorded as a fact of
provenance, not as a judgement.

### Antecedents actually read (not merely cited)

| key | source | sha256 |
|---|---|---|
| **L** | Y. Lamzouri, *A new proof that more than 2/3 of the zeros … are simple and on the critical line*, arXiv:2609.02882v2 | `305df7fcbbf96e61ccca2cbae13dcb2b7f9251b054d02dd0c78551982621046b` |
| **AF** | L. Alpöge, R. Furman, *More than two thirds of the zeta zeros are simple and on the critical line*, arXiv:2608.13637v2 | `6de3b156342e7b4a802c34f8ef40432567e9dabe006938da04233f19fc4ef444` |
| **GS** | Goldston–Sahay (with Baluyot, Suriajaya, Turnage-Butterbaugh), *Pair Correlation of Zeros … I*, arXiv:2501.14545 | `0a16b0b0f19b06490111e048a6841f4fe68717c546fae4cb7a44555dc0f96cb3` |

**The PDFs are third-party copyrighted material and are not redistributed in this
repository**; only their fingerprints are shipped, in `sources/SOURCE_HASHES.txt`.

Cited by ARTZ but **not** obtained: BGSTB, *An unconditional Montgomery theorem …*,
Acta Arith. 214 (2024) 357–376 (doi:10.4064/aa230612-20-3) — read only through the
restatements in **L** (Lemma 3.1) and **GS** (§3, "Montgomery Theorem (MT)").
Bourin–Lee 2012, Carneiro–Chandee–Littmann–Milinovich 2017, Gonek–Sahay
arXiv:2412.15481, Aryan, Montgomery 1973, Montgomery–Taylor: **not obtained**;
used only as attributions. This is recorded as a dependency, see §5 below.

## 2. My own RH program — authoritative sources

| # | object | authoritative location | status |
|---|---|---|---|
| 1 | the δ²–L³ sensitivity law (paper) | `~/GitHub/weil-window-control/paper/weil_window_control.tex`, §"The δ²–L³ law" | **published**, doi:10.5281/zenodo.23250932 (concept …21109955; previous version …23212998) |
| 2 | repository | `~/GitHub/weil-window-control`, HEAD `725f204` ("Pin the reproduction checkout to the publication commit"), working tree has only untracked build artefacts | read-only here |
| 3 | four-mode finite-window positivity | repo `extensions/`, GitHub main `0c7d834` | published, no Zenodo |
| 4 | planted-quartet negative witnesses | `extensions/planted-witness-{closure,consolidation}-2026-10-08` | in the published paper, Thm 11 |
| 5 | exact interval enclosures / certified detection witness | `extensions/certified-detection-witness-2026-10-08`; `~/scratch/rh-geometric-side-certificate-2026-10-15-105550` | certified |
| 6 | interpolation / annihilation, detector strength | `~/scratch/rh-blind-spot-covering-2026-10-15-142629/{PROBLEM_FORMULATION,INTERP_PROOF,NO_GO}.md` | proved + diagnostic |
| 7 | near-singular Gram, close-zero obstruction | `~/scratch/rh-local-invariant-gram-study-2026-10-15-120213` | proved |
| 8 | confluent interpolation, divided differences, node span | `~/scratch/rh-confluent-{proof-closure,publication-audit}-2026-10-15-*`; note at `~/scratch/rh-confluent-zenodo-note-2026-10-15-141204` | audited **GREEN**, note Zenodo-ready, **not deposited** |
| 9 | off-line vs on-line detection strength | `~/scratch/weil-offline-detect-2026-10-17-083545` | certified + diagnostic |
| 10 | **Pro's final negative assessment of uniform Route D covering** | `~/scratch/rh-blind-spot-covering-2026-10-15-142629/{VERDICT,NO_GO}.md`, 2026-10-15, 48 files | **OUTCOME D, DECISIVE NO-GO** |

### Chronology, and which assessment governs

```
2026-10-07  four-mode positivity; intrinsic-positivity bridge; OpenAI quasi-RH audit
2026-10-08  planted quartets; certified detection witness
2026-10-14  baseline weighted witness; publication of the delta^2-L^3 note
2026-10-15  quasi-RH path triage -> geometric-side certificate -> local-invariant Gram
            -> confluent proof closure -> publication audit (GREEN)
            -> BLIND-SPOT COVERING: OUTCOME D, DECISIVE NO-GO       <-- governs
2026-10-16  GUE blindness
2026-10-17  off-line detection strength (erratum E3: earlier blanket claim withdrawn)
2026-10-21  this package
```

**The 2026-10-15 no-go governs.** Its decisive content, quoted from `NO_GO.md`:

* `NG-0` (identity): a finite-window witness measures only `t² − β²`; the vertical
  offset of an ordinate and the horizontal displacement of a zero enter at the same
  order with the same coefficient and opposite signs.
* `NG-2` (cluster, decisive): one extra ordinate at distance `g` gives
  `δ ~ η√(2S/D₂)/g → ∞`; cluster size need not grow.
* `NG-3` (portability): `δ = η√(1 + S/(2κ)) > η` always — the exclusion returned is
  strictly worse than the localisation supplied.
* `κ(Γ, γ*, L) = sup{F'(γ*)² : f ∈ A_L, F(γᵢ) = 0 ∀i} = dist²(s, C(Γ))` is **already
  a supremum over the whole admissible class**, so any bound on it binds every witness
  family at once.

Route D is **not reopened** in this package. §06 and §07 below test only whether the
new paper supplies a *new hypothesis* that addresses `NG-2` / `NG-3` directly; the
answer is no, for a reason that is itself a theorem (§06 J1–J5).

## 3. Notation collision — fixed once, here

ARTZ and Lamzouri both use **δ** for the *smoothing cutoff* of the test function
(`η_δ ∈ C_c^∞(−1/2,1/2)`, `η̃_δ = 1` on `|x| ≤ 1/2 − δ`). My RH program uses **δ**
for the *horizontal displacement* `Re ρ − 1/2`. Throughout this package:

| symbol | meaning |
|---|---|
| **τ** | ARTZ's / Lamzouri's cutoff parameter (their δ), `0 < τ ≤ 1/200` |
| **δ** | horizontal displacement `β − 1/2` of a zero (my convention) |
| **a** | the *normalised* displacement `a = δ · log T` — the only dimensionless one |
| **L** | `log T / (2π)`, ARTZ's normalisation constant |
| **κ** | `K_τ(2i·Im z_ρ) = ∫ η_τ²(u) cosh(2δ u log T) du`, the depth scalar (§04) |
| **c** | the moderate-gap parameter, `0 < c ≤ 1/4` |
| **ℓ** | run length, `2 ≤ ℓ ≤ 4` |

Confusing the two δ's would make every statement in §07 and §08 wrong; the scripts
assert the distinction (`iv.py` / `verify.py` use `tau`, `delta`, `a` explicitly).

## 4. Provenance classes used in this package

| class | meaning |
|---|---|
| **I** | classical Fourier / matrix identity |
| **II** | proved by ARTZ or an antecedent (L, AF, GS, BGSTB, Montgomery) |
| **III** | previously established in my RH program |
| **IV** | exact new consequence of combining the frameworks, proved here |
| **V** | numerical observation made here |
| **VI** | conditional or open statement |

and, orthogonally, a verification tier: `[certified]` (exact rational / rigorous
interval arithmetic), `[proved]` (pen-and-paper proof given here), `[diagnostic]`
(double precision), `[cited]`.

## 5. Explicit assumptions and unresolved dependencies

1. **D1 — BGSTB Lemma 5 / Lamzouri Lemma 3.2 are taken as stated.** The only
   analytic input is the unconditional second-moment asymptotic
   `Σ_{ρ,ρ'} K_τ(i(ρ−ρ') log T/2π)² = (C_{η_τ} + O_τ(1/√log T)) N(T)`. I did not
   obtain the Acta Arithmetica paper and did not re-prove it. **§04 §5 and §09 §4
   show that this single statement, combined with elementary matrix algebra, has
   Density-Hypothesis-strength consequences.** That is the main finding of this
   package and it is the one item needing expert resolution.
2. **D2 — Montgomery–Taylor extremality is cited, not re-derived.** `C_MT` is
   recomputed to 40 digits here (`verify.py` A1) and `C_{η} = (A*A)(0) + 2∫₀¹α(A*A)(α)dα`
   is checked against it by quadrature (B8), but the *extremality* of the
   Montgomery–Taylor window among admissible windows is taken from
   [Montgomery 1975], [CCLM 2017, Cor. 14] as ARTZ do.
3. **D3 — the confluent lower bound for `κ` remains unfinished.** The
   Hermite–Genocchi remainder noted in `weil-confluent-cluster-theorem-audited`
   was never completed, so any statement of the form "moderate gaps give
   `κ ≥ …`" is **class VI**, not IV (see §06 §3).
4. **D4 — no Lean work was run and none is claimed.** AF state that their results
   are formally verified in Lean 4 and Lamzouri records Lean certificates produced
   by the Axiom Math team; I did not build, run, or audit either. Recorded as cited.
5. **D5 — no zeta zero data was used.** `mpmath`, `numpy`, `sympy` and `pip` are
   unavailable on this machine; every computation here is Python-stdlib exact
   rational / interval arithmetic or double precision on *synthetic* spectra.
   No statement in this package depends on the location of any actual zero of ζ.
