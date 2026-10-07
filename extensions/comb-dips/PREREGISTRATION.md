# Preregistration — weil-comb-dips

**Written and committed before any computation.** Any later change must be logged as a
deviation in `REPORT.md`, with the date and the reason.

**Goal.** (Q1) measure how far Weil positivity actually needs certifying — the true last
dip `t0(L)` of Zhu's symbol versus his worst-case threshold `T1(L)`; (Q2) determine whether
the dips carry any arithmetic structure beyond what the formula forces.

**Expectation going in:** the dips are generic (explained by the null models). **A null
result is a valid, publishable outcome.** Do not chase anomalies without the gates below.

---

## Definitions (Zhu, arXiv:2608.24827, eqs. (3) and Lemma 3.1)

```
Psi_L(t) = Re psi(1/4 + it/2) - log pi - P_L(t),
P_L(t)   = sum_{n : log n < 2L} (2 Lambda(n)/sqrt n) cos(t log n),
A_L      = sum_{log n < 2L} 2 Lambda(n)/sqrt n,
T1(L)    = 2 pi e^{A_L}.
```

Note: `Psi_L` depends on `L` only through the set `{n : log n < 2L}`; it changes only at
`L = ½ log n`. Treat "L-tracking" as before/after each prime-power entry.

**Dip threshold.** For `beta` in `{0, 0.25, 0.5}`, the dip set is
`D_L(beta) = {t >= 2 pi : Psi_L(t) < beta}` and `t0(L,beta) = sup D_L(beta)`. Beyond the
crude envelope (Zhu Lemma 3.1) no dip exists, which bounds the scan range.

---

## Null models

**N1 (torus measure, analytic).** With phases `theta_p` uniform and independent per prime
(prime powers share `theta_p`: `block_p(theta) = sum_k a_{p^k} cos(k theta)`), compute the
tail probability `pi_L(x) = Pr[P_L(theta) > x]` by numerically convolving the per-prime
block distributions (cross-check by Monte Carlo, >= 1e7 samples). Predict the expected dip
measure in windows `[a,b]` as `integral pi_L(H(t) - beta) dt` with
`H(t) = Re psi(1/4 + it/2) - log pi`, and a predicted `t0` as the `t` beyond which the
expected remaining dip measure falls below the median dip-interval length.

**N2 (random-frequency surrogates).** >= 200 surrogates per `L` replacing each `log p` by
an independent uniform real in `[log p - 0.3, log p + 0.3]` (keep prime-power harmonics as
`k * (new log p)`, keep all weights and the digamma term). Run the identical Stage-2
pipeline on each. Report the surrogate distribution of `t0`, dip measure, dip count.

---

## Predictions (preregistered)

- **P1:** `t0(L,beta)` agrees with the N1 prediction within a factor 3 for `L >= 0.9`.
- **P2:** real `t0`, dip measure and dip count lie inside the N2 surrogate 95% bands.
- **P3:** the deepest dips are carried by the highest-weight primes (2, 3, then 5, 7), in
  the proportions N1 predicts.
- **P4:** the `{2,3}` beat explains most dip measure for `L < 0.7` and its explained
  fraction decays as primes enter, at a rate matched by N2.

---

## Kill criteria

- **K1:** if P1–P4 all hold → write `NULL.md` (the dips are generic; report `t0/T1` per `L`
  as the practical certification gain) and stop. This is a valid final outcome.
- **K2:** if `t0(L,beta) <= T1(L)/10` for some `L <= 1.4` → write a certification-feasibility
  estimate for that `L` (Legendre order `N ≈ e^L t0 / 2`, precision from the floor
  `lambda*(L)`, compare with Zhu §7) and **STOP** before building any certificate.
- **K3:** if a statistic falls outside the N2 band → before calling it structure: rerun with
  3 fresh surrogate seeds, double precision on the affected `L`, and the alternative
  surrogate (`log p` perturbed by `±0.1` instead of `±0.3`). Report it only if it survives
  all three, and label it **T2 "candidate anomaly"**, not a result.

---

## Stage gates (also binding)

- **Stage 1 GATE:** do not proceed until Zhu's `A_L` table and the §7 reproduction at
  `L = 1.19`, `beta = 0.46948` both match. Report the match.
- **Stop after Stage 5** and summarize before Stage 6.
- Stage 6 results for `L >= 1.8` are **non-certified estimates** (LLL cannot prove absence
  of dips) and must be labelled so.
- Any `t0` used for certification must be re-verified near the boundary with interval
  arithmetic (python-flint / Arb digamma); flag which ones are.

## Output

`REPORT.md`: tiered findings (**T1** exact, **T2** numerical with precision, **T3**
interpretation), the `t0/T1` table, null-model comparisons with plots, the Stage-5
calibration result, a deviations log, and an explicit verdict:
**NULL / CERTIFICATION-FEASIBLE / CANDIDATE ANOMALY.**
