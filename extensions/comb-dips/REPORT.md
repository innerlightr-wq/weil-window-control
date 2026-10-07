# weil-comb-dips — report

Reconnaissance on the dip set of Zhu's prime-comb symbol
(Zhu, *Weil positivity in compact windows*, [arXiv:2608.24827](https://arxiv.org/abs/2608.24827)).

**Verdict: NULL.** The dips are generic. The useful output is the Q1 measurement:
Zhu's worst-case certification threshold `T1(L) = 2π e^{A_L}` is **essentially sharp** —
the true last dip sits at `t0 = (0.45 … 0.92) · T1` at `β = 0`, never below `T1/10`.
There is no certification shortcut here, and no arithmetic structure in the dips beyond
what the formula forces.

Tiers: **T1** = exact/proved, **T2** = numerical with stated precision, **T3** = interpretation.

---

## 0. What was asked and what was preregistered

`PREREGISTRATION.md` was committed (`9eb4c62`) before any computation. Two questions:

- **Q1** — how far does Weil positivity actually need certifying? True last dip `t0(L,β)`
  versus Zhu's worst-case threshold `T1(L)`.
- **Q2** — do the dips carry arithmetic structure beyond what the formula forces?

Definitions, verbatim from the preregistration:

```
Psi_L(t) = Re psi(1/4 + it/2) - log pi - P_L(t),   P_L(t) = sum_{log n < 2L} (2 Lambda(n)/sqrt n) cos(t log n)
A_L = sum 2 Lambda(n)/sqrt n     T1(L) = 2 pi e^{A_L}     T*(L) = 2 pi e^{2L}
D_L(beta) = { t >= 2 pi : Psi_L(t) < beta }           t0(L,beta) = sup D_L(beta)
```

---

## 1. Stage 1 — implementation and validation (**GATE PASSED**)

**T1.** `H(t) = Re ψ(1/4 + it/2) − log π` is **strictly increasing**: termwise,
`Im ψ'(1/4 + it/2) < 0` (proof in `notes/monotonicity.md`). This is what lets the
branch-and-bound use `min H = H(a)` *exactly* on a cell, with Lipschitz slack needed only
for the comb.

**T2.** `A_L` matches all six of Zhu's tabulated values, max `|diff| = 3.4e-4` (his rounding).
Zhu §7 reproduced at `L = 1.19`, `β = 0.46948`:

| quantity | this work | Zhu §7 |
|---|---|---|
| dip measure | 18.226 | ≈ 18 |
| max depth | 1.9854 at `t = 1549.93` | ≈ 2.0 near 1550 |
| `t0` | 8476 | ≈ 8.48e3 |
| dip intervals | 40 | — |

**Precision.** Plain float64 fails: argument reduction leaks `|t log n| · 2^-53`, giving a
worst error of **1.311e-08** at `L = 1.6` — four orders past spec. Fixed by exact
(Dekker two-product) reduction mod 2π, applied **twice**: the trap is that `q · 2π_hi` is
itself inexact once `q ~ 8e6`. Final accuracy **3.55e-15 at every L**, at 1.7× cost.
Rigorous bound: `src/symbol.py: float64_arg_error_bound`.

---

## 2. Stage 2 — certified scan (**the Q1 answer**)

Three-pass certified branch-and-bound with rigorous two-sided cell bounds
(`src/scan.py`). Every `t0` is reported as a bracket; all brackets are `< 8e-9` wide.
Boundaries re-verified at mpmath 50 dps (`data/stage2_verification.json`).

### `t0(L,β) / T1(L)`

| L | active primes | `A_L` | `T1(L)` | β=0 | β=0.25 | β=0.5 |
|---|---|---|---|---|---|---|
| 0.55 | 2,3 | 2.2488 | 5.95e1 | **0.7747** | 1.0650 | 1.0694 |
| 0.6  | 2,3 | 2.2488 | 5.95e1 | **0.7747** | 1.0650 | 1.0694 |
| 0.7  | 2,3 | 2.9420 | 1.19e2 | **0.9149** | 0.9168 | 1.2982 |
| 0.8  | 2,3 | 2.9420 | 1.19e2 | **0.9149** | 0.9168 | 1.2982 |
| 0.9  | 2,3,5 | 4.3815 | 5.02e2 | **0.6146** | 0.9566 | 1.1736 |
| 1.0  | 2,3,5,7 | 5.8525 | 2.19e3 | **0.7088** | 1.1191 | 1.4173 |
| 1.19 | 2,3,5,7 | 7.0750 | 7.43e3 | **0.7469** | 0.7470 | 1.1412 |
| 1.2  | 2,3,5,7,11 | 8.5210 | 3.15e4 | **0.6295** | 0.7715 | 0.7715 |
| 1.4  | …,13 | 10.2903 | 1.85e5 | **0.7162** | 0.7162 | 1.1213 |
| 1.6  | …,17,19,23 | 14.3232 | 1.04e7 | **0.4494** | 0.8136 | 0.8136 |

**T2 — headline.** At `β = 0` the ratio stays in **[0.449, 0.915]**. It never approaches the
`K2` kill line `T1/10`. **K2 does not trigger at any L.**

**T3.** Zhu writes that any one-stroke certificate "must resolve frequencies up to
`2π e^{A_L}` … a doubly exponential threshold that no pointwise bound on the prime comb
lowers." This measurement is the complementary statement: it is not merely that *pointwise*
bounds cannot lower the threshold — **the dips genuinely run out that far**. The best
achievable saving from knowing `t0` exactly is a factor of about 2.2 in frequency, i.e.
~1 bit off a doubly exponential quantity. No route to certification is opened.

---

## 3. Stage 3 — null models

**N1** (torus): each prime's phase `θ_p` i.i.d. uniform, prime powers sharing `θ_p`; the law
of `P_L` obtained by FFT-convolving the per-prime block laws. **N2** (surrogate): `log p →
log p + U(−σ, σ)`, harmonics `k·(new log p)`, same weights, 200 surrogates per cell.

### P1 — **HOLDS** (T2)

`t0` real/N1 ratios lie in **0.52 – 1.90** at every `L` and `β`; dip-measure ratios
**0.84 – 0.96**. All within the preregistered factor 3.

### P2 — holds on the band test, with a documented low-`t` deviation (T2)

- **4 of 24 cells** have ≥1 statistic outside the N2 95% band; **5 of 72** statistics
  (≈3.6 expected by chance). Unremarkable in count.
- **But the direction is systematic**: real dip measure is above the surrogate band centre
  in **20 of 24 cells**.

**K3 protocol** (3 fresh seeds × 2 spreads × doubled precision) run as preregistered; the
deviation survives all three. Splitting at `tsplit = 2π e^{A/2}` localises it completely:

| | low-`t` p | high-`t` p |
|---|---|---|
| five K3 configurations | 0.015, 0.015, 0.015, 0.095, 0.030 | 0.245, 0.130, 0.265, 0.200, 0.075 |

**In the asymptotic region the real comb is statistically indistinguishable from random
frequencies.** The whole excess is at low `t`.

### Dose-response (added beyond the preregistered design — see deviations)

`L = 1.6`, `β = 0`, 8 spreads, ≥100 surrogates each (`figures/dose_response.png`):

| σ | excess % | low-`t` excess % | p (low `t`) | overlap |
|---|---|---|---|---|
| 0.3 | 7.15 | 12.20 | 0.000 | 0.0385 |
| 0.1 | 6.70 | 11.91 | 0.010 | 0.0410 |
| 0.03 | 6.42 | 11.23 | 0.000 | 0.0545 |
| 0.01 | 5.96 | 9.57 | 0.000 | 0.0838 |
| 0.003 | 4.69 | 7.34 | 0.000 | 0.1389 |
| 0.001 | 3.44 | 5.14 | 0.000 | 0.2149 |
| 0.0003 | 2.78 | 3.01 | 0.010 | 0.3268 |
| 0.0001 | 1.48 | 0.64 | 0.200 | 0.4495 |

Comma detunings among the active primes (`src/commas.py`, height ≤ 200): `9/8` = 0.1178,
`25/24` = 0.0408, `128/125` = 0.0237, `81/80` = 0.0124; effective preservation scales
`s* = d / sqrt(Σe²/3)` are 0.0566, 0.0189, 0.0054, 0.0037 respectively.

**Verdict on the user's rule: EXPLAINED.** The excess fades monotonically as σ drops below
the comma detunings and is no longer significant (p = 0.20) at σ = 1e-4. It is **not** a
candidate anomaly.

**T3 — caveat, and it is a real one.** This test cannot cleanly attribute the fade to
commas. The excess lives at low `t`, and at low `t` the surrogate converges *trivially* to
the real symbol once `σ·t ≪ 1`; at `L = 1.6` the low-`t` region is `[2π, 8099]`, so trivial
convergence sets in around `σ ≈ 1.2e-4` and sweeps upward. The preservation scales `s*`
(0.0037 – 0.057) sit inside that sweep. The `overlap` column makes the confound explicit:
it rises 0.0385 → 0.4495 over exactly the range where the excess falls. So: **EXPLAINED by
known structure, with commas not uniquely implicated.**

---

## 4. Stage 4 — minima catalogue

Every local minimum of `Ψ_L` with `Ψ_L < 1`, located by certified B&B followed by
vectorised golden section (`src/stage4_minima.py`), with depth, curvature, phase vector and
per-prime alignment share `s_p = block_p(t)/A_p`. Note the sign: `Ψ = H − P`, so a dip is
where `P` is **large**; `s_p = 1` means the `p`-block is maximally aligned to deepen the dip.
20 N2 surrogates per `L`; N1 computed on the torus with the matched conditioning.

### Global statistics, real vs N2 95% band

| L | n minima | spacing CV | Fano | frac `Ψ<0` | deepest `Ψ` |
|---|---|---|---|---|---|
| 0.55 | 13 [10, 17] | 0.558 [0.21, 0.83] | — | 0.308 [0.15, 0.50] | −1.065 [−1.94, −0.90] |
| 0.7 | 18 [16, 27] | 0.651 [0.33, 1.95] | — | 0.333 [0.24, 0.46] | −1.527 [−2.79, −1.33] |
| 0.9 | 45 [38, 56] | 1.130 [0.91, 3.15] | 1.0 [0.9, 1.0] | 0.400 [0.21, 0.44] | −1.954 [−4.10, −1.88] |
| 1.0 | 101 [88, 118] | 2.095 [1.51, 3.92] | 2.3 [2.3, 3.2] | 0.416 [0.28, 0.43] | −2.559 [−4.69, −2.45] |
| 1.19 | 176 [145, 199] | 1.981 [1.60, 4.77] | 5.6 [5.7, 7.9]\* | 0.375 [0.34, 0.42] | −2.442 [−5.67, −2.99]\* |
| 1.2 | 362 [288, 386] | 2.821 [2.88, 7.62]\* | 23.7 [20.4, 26.1] | 0.356 [0.33, 0.42] | −3.236 [−7.40, −3.93]\* |
| 1.4 | 699 [655, 795] | 4.604 [3.73, 9.19] | 94.1 [85.2, 108.2] | 0.405 [0.35, 0.39]\* | −3.429 [−7.27, −4.02]\* |
| 1.6 | 3366 [3227, 3599] | 11.68 [9.32, 30.0] | 1394 [1311, 1493] | 0.386 [0.36, 0.38]\* | −5.011 [−9.90, −6.25]\* |

`*` = outside the N2 95% band. (The large Fano values are the density trend, not clustering;
the surrogates carry the same trend, which is why the comparison is to the band and not to 1.)

**T2.** Two mild, opposite deviations, both at large `L`:
- `frac Ψ<0` sits just **above** the band at `L = 1.4, 1.6` — the same low-`t` measure excess
  as Stage 3, seen again.
- the **deepest minimum is consistently *shallower*** than the surrogate band for `L ≥ 1.19`.
  Random frequencies occasionally draw a near-commensurable pair out of `U(−0.3, 0.3)`,
  producing a long near-resonance and a record-deep minimum. The real primes never do.
  **T3:** if anything, the real `log p` are *less* commensurable than random frequencies —
  the opposite direction from "hidden arithmetic structure".

### P3 — **FAILS as stated** (T2)

P3 predicted the deepest dips would be carried by the highest-weight primes (2, 3, then 5, 7)
in N1's proportions. They are not. Spearman correlation between block weight `A_p` and
alignment share in the deepest 1%:

| L | 0.9 | 1.0 | 1.19 | 1.2 | 1.4 | 1.6 |
|---|---|---|---|---|---|---|
| ρ(weight, share) | +0.50 | +0.40 | +0.80 | **−0.80** | **−0.54** | **−0.62** |

At `L ≥ 1.2` the correlation is **negative**: the deepest minima are carried slightly more by
the *light* primes. **T3** — this is a property of the formula, not of the primes. A heavy
prime's block is a sum over several powers (2, 4, 8, 16 for `p=2` at `L=1.6`), so it cannot be
driven to its maximum by a single phase; a light prime contributes one cosine, trivially
alignable. The constraint falls on the heavy blocks, so they end up *less* aligned at the
deepest minima. **N2 surrogates reproduce this**: 4 of 41 real/band comparisons fall outside
the 95% band (≈2 expected), with no systematic direction.

Carrier sets of the deepest 1% (primes aligned past 1/2) are near-complete sets at every `L`
— e.g. at `L = 1.6`: `{3,5,7,11,13,17,23}`, `{2,3,5,7,11,13,19,23}`, `{2,3,5,7,11,17,19,23}`.
A deep dip needs almost everything aligned, which is exactly what N1 says.

*Caveat on power:* "deepest 1%" is 1 minimum for `L ≤ 1.0` and 2–7 for `L = 1.19 … 1.4`.
Only `L = 1.6` (34 minima) carries real statistical weight.

---

## 5. Stage 5 — calibration on the two-prime case

### (a) `{2,3}` exactly (T1)

For `0.549 < L < 0.693` the active set is `{2,3}` and
`Ψ_L(t) = H(t) − a2 cos(t log 2) − a3 cos(t log 3)` with
`a2 = 2 log2 / √2 = 0.980258143`, `a3 = 2 log3 / √3 = 1.268568201`, `A = 2.248826345`.
Envelope period `2π / log(3/2) = 15.496242`. `log3/log2 = 1.584962501`, convergents
`1, 1/2, 2/3, 5/8, 12/19, 41/65, 53/84, 306/485`.

All 5 dips in `[2π, 59.5433]` are accounted for by the beat: the deficit `A − P` at the
alignment times `t = 2πm/log 2` is smallest at `m = 5` (0.1389) and `m = 2` (0.6569), and the
deepest dip (`Ψ = −1.0654`) sits at `m = 2`; the last dip (`t0 = 46.129`) sits at `m = 5`.
`data/stage5_alignments.json`.

### (b)/(c) P4 — **FAILS as stated**, but for a null-model reason (T2)

P4 predicted the `{2,3}` beat explains most dip measure for `L < 0.7` with the explained
fraction *decaying* as primes enter, at a rate matched by N2.

The `{2,3}` **dip** set is useless as a model past `L = 0.7` — it is capped at
`T1(0.55) = 59.54` while the dips of `L = 1.6` run to 4.7e6, so a covered-fraction test is
zero for trivial reasons. The right object is the `{2,3}` **alignment** set
`S_23(c) = { t : P_23(t) ≥ c (a2+a3) }`, defined for all `t`. At `c = 1/2`
(chance density 0.18703):

| L | real frac in `S_23` | N1 prediction | ratio | N2 95% band | in band |
|---|---|---|---|---|---|
| 0.6 | 0.7093 | 0.7731 | 0.918 | [0.602, 0.888] | yes |
| 0.7 | 0.7543 | 0.7563 | 0.997 | [0.667, 0.944] | yes |
| 0.9 | 0.6839 | 0.6266 | 1.091 | [0.621, 0.985] | yes |
| 1.0 | 0.6108 | 0.5808 | 1.052 | [0.609, 0.814] | yes |
| 1.19 | 0.7162 | 0.5840 | 1.226 | [0.667, 0.793] | yes |
| 1.2 | 0.7345 | 0.5838 | 1.258 | [0.686, 0.842] | yes |
| 1.4 | 0.7205 | 0.5841 | 1.234 | [0.702, 0.759] | yes |
| 1.6 | 0.7245 | 0.5770 | 1.256 | [0.711, 0.786] | yes |

**T2.** The real fraction does **not** decay — it is flat at ≈0.72 from `L = 0.7` to `L = 1.6`,
while the N1 pointwise prediction falls 0.77 → 0.58. So P4's decay prediction is wrong.
But the real value is **inside the N2 band at every single `L`**, and the bands are centred on
the same flat ≈0.73–0.75.

**T3.** The failure is in the N1 baseline, not in the primes. N1 treats the blocks as
independent *pointwise*; a dip is an *interval*, and which intervals exist is set by the
heaviest, lowest-frequency blocks. Surrogates inherit that interval structure and therefore
reproduce the real value exactly. Nothing arithmetic is being measured here.

---

## 6. Stage 7 — prior art

| Work | Overlap | Difference |
|---|---|---|
| Zhu, [arXiv:2608.24827](https://arxiv.org/abs/2608.24827), §7 | computes the dip set, measure, depth and `t0` at the single point `L = 1.19`, `β = 0.46948` — reproduced here as the Stage 1 gate; states `T1 = 2π e^{A_L}` and that the comb constant is optimal by Weyl equidistribution | does not scan `t0(L,β)` across `L`, and does not report `t0/T1`. Zhu's claim is that *no pointwise bound* lowers the threshold; the measurement here is that the *actual* last dip is within a factor ≈2.2 of it |
| Groskin, [arXiv:2607.02828](https://arxiv.org/abs/2607.02828) | finite Guinand–Weil dictionary and a two-sided certification rule with budget `B_T ~ (2N+1) ρ log T / (π² T)` for the same CvS/CCM truncation | works on the zero side and the Galerkin matrix; says nothing about the geometric-side comb's dip set |
| Odlyzko–te Riele 1985 (Mertens) | the same lattice technique that Stage 6 would use: reduce "make many `\|γ y − ψ\|` small" to CVP and run LLL/Babai | aligns **zeros** `γ`; Stage 6 would align **`log p`**. Dual lattice, same method. The precedent is for the *method*, not the object |
| Resonance method (Voronin; Soundararajan; Hilberdink) and large values of Dirichlet polynomials (e.g. [arXiv:2509.09771](https://arxiv.org/abs/2509.09771)) | the general problem of making `Σ a_n n^{-1/2} cos(t log n)` large | asymptotic lower bounds for `ζ` extremes, not a certified catalogue of where a *fixed finite* comb exceeds a *fixed* trend |

No prior work found that scans `t0(L,β)/T1(L)`, or that tests the dip set of this comb against
a null model. Both appear to be new; both come out negative.

---

## 7. Deviations from the preregistration

1. **The brief's "one `L` per distinct prime-power set" is wrong.** `L = 0.55` and `0.6` share
   `{2,3}`; `0.7` and `0.8` share `{2,3,4}`. `STAGE2_L` has 10 entries but only **8 distinct
   symbols**. All per-symbol work (Stages 3–5) is done on the 8.
2. **Interval arithmetic.** python-flint / Arb were unavailable. Rigorous interval digamma was
   replaced by **mpmath at 50 dps plus an explicit error budget**; every reported `t0` is a
   two-sided bracket from the certified B&B, re-verified at 50 dps
   (`data/stage2_verification.json`). The preregistration's "flag which ones are" is answered:
   **none** use Arb; all use the 50-dps + explicit-budget substitute.
3. **The dose-response and the low-`t`/high-`t` region split were added beyond the
   preregistered design**, at the user's instruction mid-project, after the K3 protocol flagged
   the P2 direction. They are reported in §3 and are *not* part of the preregistered test
   battery. The `overlap_fraction` confound control in §3 is a further addition of my own.
4. **The P4 statistic was redefined** (§5b): the preregistered "{2,3} explained fraction" via
   the `{2,3}` *dip* set is vacuous past `L = 0.7`. Replaced by the `{2,3}` *alignment* set with
   an N1 conditional baseline. The redefinition was made before seeing the result but after
   seeing that the original is degenerate; it is a deviation either way.
5. **Stage 4's "deepest 1%"** is a single minimum for `L ≤ 1.0`. The preregistered statistic has
   essentially no power below `L = 1.4`; reported anyway, with the caveat.

---

## 8. Kill criteria — disposition

- **K1** (P1–P4 all hold → NULL): **not literally met.** P1 and P2 hold; **P3 and P4 fail as
  stated.** But both fail identically for random-frequency surrogates, so the failures are
  about the predictions and the N1 baseline, not about the primes. The substantive conclusion
  is the one K1 contemplates.
- **K2** (`t0 ≤ T1/10` for some `L ≤ 1.4` → certification-feasibility estimate, then STOP):
  **does not trigger.** Minimum ratio over `L ≤ 1.4` at `β = 0` is 0.6146 (`L = 0.9`); over all
  `L`, 0.4494 (`L = 1.6`). No feasibility estimate written, no certificate attempted.
- **K3** (statistic outside the N2 band → 3 seeds, doubled precision, alternative surrogate
  before calling it structure): **invoked and completed.** The P2 dip-measure direction survived
  all three, was localised entirely to low `t`, and the dose-response then marked it
  **EXPLAINED** (§3), not a candidate anomaly.

---

## 9. Verdict

**NULL.**

- **Q1 — answered.** `t0(L,0)/T1(L) ∈ [0.449, 0.915]` over `L = 0.55 … 1.6`. Zhu's worst-case
  threshold is essentially sharp. The practical certification gain from exact knowledge of the
  last dip is a factor of ≈2.2 in frequency — nothing against a doubly exponential
  `T1 ~ 2π e^{4e^L}`. **No certification shortcut exists on this route.**
- **Q2 — answered, negatively.** No arithmetic structure was found. P1 holds. P2 holds on the
  preregistered band test; its one systematic deviation is confined to low `t`, fades
  monotonically under dose-response, and is confounded with trivial surrogate convergence —
  **EXPLAINED**. P3 and P4 fail as stated, and fail the same way for random frequencies. The
  one clean large-`L` signal points the *other* way: the real comb's deepest minimum is
  consistently **shallower** than random surrogates achieve.

**Not claimed.** That no `L` outside `[0.55, 1.6]` behaves differently; that the low-`t` excess
has no cause (only that this experiment cannot separate commas from trivial convergence); that
Stage 6 LLL extrapolation would agree — it was not run.

---

## Reproducing

```
python3 src/stage1_validate.py     # gate
python3 src/stage2_scan.py         # t0(L,beta), certified
python3 src/stage2_verify.py       # mpmath 50-dps boundary re-verification
python3 src/stage3_nulls.py        # N1, N2
python3 src/stage3_k3.py           # K3 protocol + region split
python3 src/stage3_dose.py ; python3 src/stage3_dose_fine.py
python3 src/stage4_minima.py       # minima catalogue (~2 min)
python3 src/stage5_twoprime.py ; python3 src/stage5c_p4.py
python3 src/fig_dose.py ; python3 src/fig_summary.py
```

Figures: `figures/summary.png` (Q1, P3, P4), `figures/dose_response.png` (dose-response and
its confound).
