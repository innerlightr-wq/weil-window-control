# NULL

The preregistered null outcome, written under kill criterion **K1**.

`PREREGISTRATION.md` §K1: *"if P1–P4 all hold → write `NULL.md` (the dips are generic; report
`t0/T1` per `L` as the practical certification gain) and stop. This is a valid final outcome."*

## K1 was not literally met

P1 and P2 hold. **P3 and P4 fail as stated.** They fail *identically for random-frequency
surrogates*, so they are failures of the predictions and of the N1 pointwise baseline, not
signals from the primes:

- **P3** predicted the deepest dips are carried by the highest-weight primes. The Spearman
  correlation between block weight `A_p` and deep-dip alignment share is **negative** at
  `L ≥ 1.2` (−0.80, −0.54, −0.62). A heavy prime's block sums over several powers and cannot
  be maximised by one phase; a light prime contributes one cosine. 4 of 41 real-vs-band
  comparisons fall outside the N2 95% band (≈2 expected), with no systematic direction.
- **P4** predicted the `{2,3}` explained fraction decays as primes enter. It is **flat at
  ≈0.72** from `L = 0.7` to `L = 1.6` while the N1 prediction falls 0.77 → 0.58 — and the real
  value is inside the N2 band at **every** `L`.

The substantive conclusion is nevertheless the one K1 contemplates, so this file is written.

## The dips are generic

No arithmetic structure was found in the dip set of `Ψ_L` beyond what the formula forces.

- **P1 holds.** `t0` real/N1 ratios 0.52–1.90, dip-measure ratios 0.84–0.96 (spec: factor 3).
- **P2 holds** on the preregistered band test: 5 of 72 statistics outside the N2 95% band,
  ≈3.6 expected by chance.
- Its one systematic feature — real dip measure above the band centre in 20 of 24 cells —
  survived the full **K3** protocol but is confined entirely to **low `t`** (low-`t` p =
  0.015–0.095; high-`t` p = 0.075–0.265). In the asymptotic region the real comb is
  statistically indistinguishable from random frequencies.
- The dose-response test then marked it **EXPLAINED**: the low-`t` excess fades monotonically
  12.20% → 0.64% (p: 0.000 → 0.200) as the surrogate spread falls 0.3 → 1e-4, i.e. below every
  comma detuning among the active primes. It is not a candidate anomaly. The test cannot
  uniquely implicate commas, because at low `t` the surrogate also converges trivially to the
  real symbol once `σ·t ≪ 1`; see `REPORT.md` §3.

## The practical certification gain

`t0(L, 0) / T1(L)`, `T1(L) = 2π e^{A_L}`, certified two-sided, brackets `< 8e-9`, boundaries
re-verified at mpmath 50 dps:

| L | 0.55 | 0.6 | 0.7 | 0.8 | 0.9 | 1.0 | 1.19 | 1.2 | 1.4 | 1.6 |
|---|---|---|---|---|---|---|---|---|---|---|
| `t0/T1` | 0.775 | 0.775 | 0.915 | 0.915 | 0.615 | 0.709 | 0.747 | 0.629 | 0.716 | 0.449 |

Range **[0.449, 0.915]**. **K2 does not trigger** at any `L`: the ratio never approaches the
`T1/10` kill line.

**The gain is a factor of ≈2.2 in frequency.** Against `T1 ~ 2π e^{4e^L}` that is about one bit
off a doubly exponential quantity. Zhu's worst-case threshold is sharp in practice, and
**no certification shortcut exists on this route.**

## Stop

Per K1, work stops here. Stage 6 (LLL extrapolation to `L ≥ 1.8`) was not run; by the
preregistration's own rule its output would have been non-certified estimates in any case.

Full detail, tiering, deviations log and prior art: `REPORT.md`.
