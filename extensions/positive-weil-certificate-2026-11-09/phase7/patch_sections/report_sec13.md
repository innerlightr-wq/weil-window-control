
---

# 13. PHASE 7b (2026-11-12) — a second window, and what it tests

Same machinery at a second window, to make the comparison with
arXiv:2611.13129v1's `lambda^2 = 5` row as close as this code allows.

## 13.1 The window

`L~ = 80471895621705018/100000000000000000`, strictly below `L* = (log 5)/2` by
`7.3e-18`, so `lambda_-^2 = 4.999999999999999927`. The point of choosing it
below `L*` is that `log n < 2 L~` keeps `n = 5` out of the prime comb, so
`A_L~ = 2.941973525223620455504` is **identical** to `A_{4/5}`,
`T_gate = 119.086906076` is unchanged, and `T = 400` remains admissible — the
run costs one assembly (`3678 s`, peak RSS `393.3 MiB`) rather than the ten
hours an `A_L` jump would force. (`A_L` jumps by `2 log 5/sqrt 5` for
`L > L*`, not at `L*`: at `L*` itself the strict inequality still excludes
`n = 5`, and in the time domain that term carries `C(2L) = 0`. So `L*` is itself
a cheap window — see §13.4.)

## 13.2 The certificates

Constants through `phase5/src/constants.py` at `prec = 2000`, all passing
`require_tight`; blocks `220x220`, max radius `1.36e-48`;
`beta_400 = 1.211613695058064549472`, `kappa = 7.756529961646111179258e-29`,
`d_even = d_odd = 1.211613695058064549472`.

| sector | certified `a` = `mu` (lower endpoint) | Rump `||I-G||_inf` |
|---|---|---|
| even | `8.273438682556152343750000e-18` | `5.05672988450e-11` |
| odd | `7.885009765625000000000000e-15` | `4.01037156749e-16` |

## 13.3 What it tests

**Monotonicity, rigorously.** `eps_lambda` is non-increasing in `lambda`
(arXiv:2611.13129v1 Prop. 8.17, unconditional), so with
`4.95303 < lambda_-^2`:

```
    mu_even(lambda_-) = 8.2734e-18  <=  eps(lambda_-)  <=  eps(e^{4/5})  <=  1.674e-17
```

which holds. This is the only falsification test the pair of runs supports: a
lower bound at the *larger* window exceeding our own certified upper bound at
the *smaller* one would contradict Prop. 8.17 together with
`phase5/certs/upper_highN.json`. It does not. **A lower bound at `lambda_-`
exceeding an external upper bound at exactly `lambda^2 = 5` would not be
contradictory** — the smaller window legitimately has the larger minimum — and
an earlier draft of this section wrongly asserted otherwise.

**Two-window scaling — a DIAGNOSTIC, not a certified scaling result.** The two
certified lower bounds are **certificate floors, not eigenvalues** (each sits at
roughly `90.6 %` of its `eps` in the even sector, and the two happen to sit at
similar fractions, which is an observation), and the rungs below are
**high-precision values, not interval enclosures**. The agreement is evidence
that the two runs and the dictionary are mutually consistent; nothing is
inferred from it:

| | `mu(lambda^2=5.000)/mu(lambda^2=4.953)` | exact rung ratio | rel. diff |
|---|---|---|---|
| even | `0.577351` | `delta_4(5)/delta_4(4.953) = 0.579867` | **`0.43 %`** |
| odd | `0.569027` | `delta_6(5)/delta_6(4.953) = 0.592859` | `4.0 %` |

The even agreement to `0.43 %` is consistent with the exponent structure
`-log delta_4 = 4 pi lambda^2 - 9 log lambda - l_4`, and uses none of the
paper's eigenvalues — but it compares floors with uncertified rungs and is a
diagnostic only.

**Matched-window ratios**, against the paper's own printed
`delta_4(sqrt 5) = 1.2983880087183625604e-18`,
`delta_6(sqrt 5) = 2.02613988982e-15`:

```
    eps_lambda(1)/delta_4      ~>  6.3720849      (paper, N=960: 7.06)
    eps^odd_lambda(1)/delta_6  ~>  3.8916413      (paper, N=960: 4.84)
```

**These are NUMERICAL COMPARISONS, not established inequalities.** The
denominators are evaluated at `lambda^2 = 5` while the numerators are certified
at `lambda_-^2 = 5 - 7.3e-17`, and they are high-precision values rather than
outward enclosures. Consistency only. The paper's values are **rounded printed numbers**, not
outward enclosures, so they bound nothing; and the residual window mismatch
(`7.3e-18` in `L`, about `8.4e-16` relative in `eps` by the derivative
estimate) is an **uncertified** transport.

## 13.4 What is still needed for a genuine head-to-head

1. An `arb`-parameterized window. `block.Symbol`, `baseline.{comb_terms, A_L,
   M_sup}` and `block.assemble` all take a rational `L_num/L_den`, and
   `L* = (log 5)/2` is irrational. Re-parameterizing is a change to
   sign-critical code and must be validated by re-running the full 176-check
   sequence at `L = 4/5` through the new signature and requiring bit-identical
   certificates before use. Once done, the matched run costs one assembly.
2. An outward enclosure for the paper's side, rather than a printed decimal.
3. A quadrature-free **upper** bound at `lambda_-` (the
   `phase5/src/upper_highN.py` route needs no gate). By monotonicity it would
   bound `eps(lambda)` for every `lambda >= lambda_-`, including
   `lambda^2 = 5`, rigorously — the one piece that would let the comparison
   bite from this side.

Three defects in the new driver code for this phase, all precision ownership or
search logic, all failures to certify rather than false certifications, are
recorded as `RETRACTIONS.md` R15 with the rule they imply.
