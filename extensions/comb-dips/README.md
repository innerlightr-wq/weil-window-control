# comb-dips — where the prime comb actually dips

An extension to `weil-window-control`. **Not part of the v1.0 technical note**
([doi:10.5281/zenodo.23195177](https://doi.org/10.5281/zenodo.23195177)); nothing here has
been through peer review or a DOI.

Epistemic tiers used throughout: **T1** exact or proved, **T2** numerical with stated
precision, **T3** interpretation.

## Purpose

Two questions about the symbol of Zhu, *Weil positivity in compact windows*
([arXiv:2608.24827](https://arxiv.org/abs/2608.24827)):

```
Psi_L(t) = Re psi(1/4 + it/2) - log pi - P_L(t),  P_L(t) = sum_{log n < 2L} (2 Lambda(n)/sqrt n) cos(t log n)
A_L = sum 2 Lambda(n)/sqrt n     T1(L) = 2 pi e^{A_L}
D_L(beta) = { t >= 2 pi : Psi_L(t) < beta }      t0(L,beta) = sup D_L(beta)
```

- **Q1** — measure the **true last dip** `t0(L)` against Zhu's worst-case certification
  threshold `T1(L)`. How far does Weil positivity actually need certifying?
- **Q2** — do the dips carry **arithmetic structure** beyond what the formula forces?

The design, predictions P1–P4 and kill criteria K1–K3 were fixed in
[`PREREGISTRATION.md`](PREREGISTRATION.md), committed before any computation. A null result
was declared in advance to be a valid outcome.

## Verdict: NULL

See [`NULL.md`](NULL.md) for the preregistered null artifact and [`REPORT.md`](REPORT.md) for
the full tiering, the deviations log and the prior-art table.

### Q1 — `t0(L,0) / T1(L)` (T2; certified brackets `< 8e-9`, boundaries re-verified at 50 digits)

| L | 0.55 | 0.6 | 0.7 | 0.8 | 0.9 | 1.0 | 1.19 | 1.2 | 1.4 | 1.6 |
|---|---|---|---|---|---|---|---|---|---|---|
| `t0/T1` | 0.775 | 0.775 | 0.915 | 0.915 | 0.615 | 0.709 | 0.747 | 0.629 | 0.716 | 0.449 |

Range **[0.449, 0.915]**. The ratio never approaches the `T1/10` kill line, so **K2 does not
trigger** at any `L`.

**Zhu's threshold is sharp in practice.** Knowing the last dip exactly buys a factor of ≈2.2 in
frequency — about one bit off `T1 ~ 2π e^{4e^L}`. **No certification shortcut exists on this
route** (T3). Zhu's own claim is that no *pointwise* bound on the comb lowers the threshold;
this is the complementary measurement — the dips genuinely run out that far.

### Q2 — no arithmetic structure found

- **P1 holds** (T2). `t0` real/N1 ratios 0.52–1.90, dip-measure ratios 0.84–0.96; spec was a
  factor of 3.
- **P2 holds** on the preregistered band test (T2): 5 of 72 statistics outside the N2 95%
  band, ≈3.6 expected by chance. Its one systematic feature — real dip measure above the band
  centre in 20 of 24 cells — survived the full **K3** protocol but is confined entirely to
  **low `t`**; in the asymptotic region the real comb is statistically indistinguishable from
  random frequencies.
  The dose-response test marks it **EXPLAINED**: the low-`t` excess fades monotonically
  **12.20% → 0.64%** (p: 0.000 → 0.200) as the surrogate spread `σ` falls 0.3 → 1e-4, below
  every comma detuning among the active primes. **Caveat (stated, not resolved):** the test
  cannot uniquely implicate comma-type near-relations, because at low `t` the surrogate also
  converges *trivially* to the real symbol once `σ·t ≪ 1`, and the two ranges overlap. The
  overlap control is plotted alongside the excess in `figures/dose_response.png`.
- **P3 and P4 fail as stated** (T2) — and fail **identically for random-frequency surrogates**.
  These are failures of the predictions and of the N1 pointwise baseline, not signals.
  - P3: Spearman correlation between block weight `A_p` and deep-dip alignment share goes
    *negative* at `L ≥ 1.2` (−0.80, −0.54, −0.62). A heavy prime's block sums over several
    powers and cannot be maximised by one phase; a light prime contributes one cosine (T3).
  - P4: the `{2,3}` explained fraction does not decay — see below.

### Deepest-minimum observation (T2)

For `L ≥ 1.19` the real comb's **deepest minimum is consistently shallower** than the N2
surrogate 95% band (e.g. `L = 1.6`: `Ψ_min = −5.011` against a band of `[−9.90, −6.25]`).

**T3, untested:** this is likely a surrogate artifact rather than a property of the primes.
Under a `±0.3` perturbation of each `log p`, surrogate frequencies can *collide* — `log 16`
and `log 17` differ by only 0.061, and harmonics amplify the perturbation (`log 16 = 4 log 2`
moves by `4σ`). A collision merges two frequencies into one of double weight and produces a
long near-resonance with a record-deep minimum, which the fixed real frequencies never
achieve. **This hypothesis was not tested.** Recorded here because, if correct, the one clean
large-`L` signal points *away* from hidden arithmetic structure, not toward it.

### Two-prime beat calibration (T2)

For `0.549 < L < 0.693` the active set is `{2,3}` and the beat picture is exact (T1):
`a2 = 2 log2/√2 = 0.980258143`, `a3 = 2 log3/√3 = 1.268568201`, envelope period
`2π/log(3/2) = 15.496242`; all 5 dips in `[2π, 59.54]` are accounted for by the convergents of
`log3/log2`.

Measured against the `{2,3}` alignment set `S_23(1/2)`, the beat explains **≈72% of dip
structure at every `L` from 0.7 to 1.6** — flat, not decaying as P4 predicted. But the N2
surrogates give the same flat ≈0.73–0.75 and the real value is inside the band at every `L`.
**T3:** this is a weight-concentration effect of the formula — the two heaviest, lowest-
frequency blocks set which intervals are dips at all — and is not arithmetic.

### Stage 6 not run

LLL extrapolation to `L ≥ 1.8` was **not run, by decision**. Given the above the expected
outcome is the same genericity, and by the preregistration's own rule its output would have
been non-certified estimates (LLL cannot prove the absence of a dip).

## Non-claims

1. **Nothing here bears on the Riemann Hypothesis.** It measures a property of a finite
   trigonometric sum against a smooth trend.
2. **All null-model conclusions are T2**, not proofs. "No structure found" is a statement about
   two specific null models (N1 torus, N2 random-frequency surrogate) at specific `L`, not a
   theorem that none exists.
3. Nothing is claimed about `L` outside `[0.55, 1.6]`.
4. The low-`t` excess is marked EXPLAINED, **not** explained-by-commas: this experiment cannot
   separate comma destruction from trivial surrogate convergence.
5. The deepest-minimum collision hypothesis is T3 and untested.

## Prior art

Checked against Zhu §7, Groskin ([arXiv:2607.02828](https://arxiv.org/abs/2607.02828)),
Odlyzko–te Riele 1985, and the resonance / large-values-of-Dirichlet-polynomials literature.
Zhu computes the dip set, measure, depth and `t0` at the **single point** `L = 1.19`,
`β = 0.46948` — reproduced here as the Stage 1 gate. **No prior scan of `t0(L,β)/T1(L)` across
`L` was found**, and no prior test of this dip set against a null model. Full table in
[`REPORT.md`](REPORT.md) §6.

## Layout

```
PREREGISTRATION.md  Stage 0, committed before any computation (incl. the deviations it fixed)
REPORT.md           full results, tiering, deviations log, prior art
NULL.md             the preregistered null artifact (K1)
src/                symbol + certified scanner + null models + catalogue + calibration
data/               all results: CSVs, JSONs, and the certified dip intervals (.npy)
figures/            summary.png (Q1, P3, P4), dose_response.png (dose-response + confound)
notes/              monotonicity.md -- the T1 proof that H' > 0
scripts/            quick_check.sh
```

## Reproduction

Needs `numpy`, `scipy`, `mpmath`, `numba`, `matplotlib`. Runtimes measured on an Apple
Silicon laptop.

```
make comb-dips-quick                   # from the repo root -- gate + 2 Stage 2 rows   ~15 s
```

Full pipeline, from `extensions/comb-dips/`:

```
python3 src/stage1_validate.py     # Stage 1 gate: A_L table, Zhu section 7          ~12 s
python3 src/stage2_scan.py         # Stage 2: certified t0(L,beta), all cells         ~30 s
python3 src/stage2_verify.py       # 50-dps boundary re-verification                  ~2 min
python3 src/stage3_nulls.py        # N1 and N2 tables (200 surrogates x 24 cells)    ~14 min
python3 src/stage3_k3.py           # K3 protocol: 3 seeds x 2 spreads + region split ~20 min
python3 src/stage3_dose.py         # dose-response, 3 cells x 5 spreads              ~12 min
python3 src/stage3_dose_fine.py    # dose-response, L=1.6 beta=0, 8 spreads          ~15 min
python3 src/stage4_minima.py       # minima catalogue + N1/N2 comparison             ~100 s
python3 src/stage5_twoprime.py     # {2,3} beat, exact                                ~20 s
python3 src/stage5c_p4.py          # P4 against the N1 conditional baseline           ~5 min
python3 src/fig_summary.py         # figures/summary.png                               ~5 s
python3 src/fig_dose.py            # figures/dose_response.png                         ~5 s
```

Every script resolves its paths relative to `extensions/comb-dips/`, so they run from that
directory with no configuration.

Sources of the work: originally developed in a standalone repository, imported here through
commit `91f4a5c`.

## License

Code (`src/`, `scripts/`) MIT; prose, data and figures CC BY 4.0 — same split as the parent
repository, see the root [`LICENSE`](../../LICENSE) and
[`LICENSE-CC-BY-4.0`](../../LICENSE-CC-BY-4.0).
