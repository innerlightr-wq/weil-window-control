# RESULT — verdict: **COMPLETE ADDED-QUARTET WITNESS CERTIFICATE**

Exploration, 2026-10-08. No installs, commits, pushes, uploads, branch changes, resets,
background jobs or manuscript edits.

## Answer to the closing question

**The precise object certified.** For the explicit admissible function `f` of
`CERTIFICATE_THEOREM.md` §1 — four exact rational coefficients in the `EVEN_D` span,
`L = 4/5`, `γ = 14`, both exactly rational — the Weil functional of the zero configuration

> the nontrivial zeros of `ζ`, **together with** the four extra zeros
> `1/2 ± δ ± 14i`, each counted once,

namely `Q_add,δ(f) = Q_ζ(f) + 4[A_δ(f)² − B_δ(f)²]`.

**The point and interval covered, and the rigorous margins.**

```
delta = 2/5                 :  Q_add in [ -0.0146094664641386 , -0.0145107657498228 ]
                               so  Q_add,2/5(f) <= -eta_point,   eta_point = 1.45108e-2
delta in [1/4, 49/100]      :  Q_add,delta(f)  <= -eta_interval, eta_interval = 2.45491e-3
                               for EVERY delta in the interval (256-fold subdivision)
```

| | margin `η` | total error allowance | `η`/error |
|---|---|---|---|
| point `δ = 2/5` | `1.45108e-2` | `4.93503e-5` | **294** |
| interval-wide (worst at `δ = 1/4`) | `2.45491e-3` | `4.93503e-5` | **50** |

**Does anything about actual zeta-zero locations follow? No.** The configuration *adds* four
zeros to `ζ`'s own and removes none; `γ = 14` is a chosen rational ordinate, not a zero of
`ζ`. The certificate does **not** say `Q_ζ(f) < 0` — it says `Q_ζ(f) ∈ [+0.00512, +0.00523]`,
strictly positive, consistent with Weil positivity. Negativity of a Weil functional for a
fabricated configuration is the **detection** direction of the classical criterion; no
`L`-function is asserted to realise this spectrum, and no actual off-line zero is excluded.

## Why this is a *complete* certificate

Every number is an interval with exact `Fraction` endpoints, rounded outward on a `2^-260`
dyadic grid. **No float operation and no libm call occurs in the certified path.** `π`,
`log 2`, `log 3`, `log π`, `exp`, `sin`, `cos`, `sqrt`, `sinh`, `cosh` are produced by
truncated series with the proved remainders tabulated in `CERTIFICATE_THEOREM.md` §5.1;
`sqrt` brackets are *verified* and the script raises if they fail. Input, arithmetic,
analytic-tail and parameter-range errors all propagate into the final strict upper bound.

The whole error budget is one term: the archimedean tail allowance `±4.935e-5` at
`M = 2^17 − 1`. Everything else is below `1e-70`. That term is linear in `1/M` and can be
reduced at linear cost.

Two things that would otherwise have blocked certification were removed analytically rather
than by precision:

- **No digamma and no Lerch transcendent.** The archimedean term is a series of *elementary*
  closed-form integrals. Limitation 4 of the manuscript records that `mpmath`'s interval
  context has neither; this is a route around that, with the `O(1/M)` tail as its price.
- **No Euler-constant enclosure.** `H_M − γ_E = ln(M+1) − T` with
  `T ∈ [1/(2(M+1)) − 1/(6M²), 1/(2M)]`, proved from `x²/2 − x³/3 ≤ x − ln(1+x) ≤ x²/2`, and
  `M + 1 = 2^17` makes `ln(M+1) = 17 log 2`.

**The one analytic input not reproved** is Weil's explicit formula itself (§4.1 of the
theorem file), cited as classical. Its applicability to this `f` is argued, not assumed.

## The interval was closed by subdivision, not by sampling

A **single** interval evaluation over the whole `δ`-box returns `+0.0279` — valid but
useless, because `δ` occurs in several places in the closed forms (`δ`, `δ²`, `cosh(δL)`,
`sinh(δL)`) and independent interval occurrences inflate the result. That attempt is
recorded in `certificate.json` as `interval_single_box_attempt`. The certificate uses `256`
**closed** subintervals covering `[1/4, 49/100]` with no gaps; interval arithmetic on each
is valid for every `δ` in it. All 256 upper bounds are negative; the worst is on
`[1/4, 803/3200]`.

The positive term `4A_δ²` is **retained throughout** — at `δ = 2/5` it is `+3.327e-4`,
working against the certificate — as is the node residual. Nothing approximately vanishing
was treated as exactly vanishing: `F(14) = 4.1293580707791e-13` for the frozen rational
vector, enclosed, and the two candidate functionals differ by exactly `2F(14)² = 3.41e-25`.

## What was corrected from the preceding run

Four substantive items, detailed in `RECONCILIATION.md` §2:

1. **The Taylor-regime claim was wrong.** "`R ≤ 2K_effδ²` fails at every `δ`" is false; it
   **holds for `δ ≲ 0.166`** (ratio `0.330` at `δ = 1/10`, `0.804` at `3/20`, crossing `1`
   between `0.166` and `0.167`). The real obstruction is **non-overlap** with the range
   where the signal beats the baseline (`δ ≳ 0.204`).
2. **The basis-size claim was overstated.** The bases are nested, so the four-mode witness
   survives into `N = 8, 16` by zero padding with the same value — invariance that is
   structural here, since the computation never references `N`. The worse `N = 8, 16`
   numbers describe that *selection rule*, not an impossibility.
3. **The uniform `θ < 1/4` proposal does not reach small `δ`** (it competes with `δ²`, so a
   fixed `θ` needs `δ² > θ`), and for a fixed unmodulated span `K_eff ≤ ‖P_{E_N}r_γ‖² → 0`
   as `|γ| → ∞` — enclosed at `3.60e-2, 1.10e-5, 2.13e-6, 1.13e-7` for `γ = 14, 70, 140,
   280`. So no all-height result follows from a fixed low-frequency basis.
4. **The error allowance was understated and the object was the wrong one.** The previous
   `±3.76e-6` was an entry-wise aggregation at a different `M` (consistent with the present
   formula); the `~6800×` margin was computed at `δ = 1/2` and is replaced by the separated
   `294×` / `50×` above. And the previous script evaluated the **count-preserving
   replacement** functional, not the added-quartet one — they differ by `2F(14)²`, bounded
   here rather than ignored.

**Not claimed:** the onset decimal `0.204201571681` is **not** a certified threshold; no
root enclosure was performed.

## Commands and actual results

```
$ python3 verify_certificate.py            # ~5 min first run; then uses qzeta_cache.json
[1] Q_zeta(f): certified enclosure  (this is the slow step, ~5 min)
[2] node residual and the two functionals
[3] point certificate at delta = 2/5
    Q_add,2/5 = [-0.014609466464138512, -0.014510765749822711]   negative: True
[4] interval certificate on delta in [1/4, 49/100] by subdivision
    worst subinterval upper bound = -0.002454910 on [0.250000, 0.250937]   all negative: True
[5] corrections A and C

  GATE  Q_zeta enclosure is strictly positive                PASS
  GATE  point delta=2/5 strictly negative                    PASS
  GATE  interval [1/4,49/100] strictly negative everywhere   PASS
  GATE  arch tail allowance below the point margin           PASS

wrote certificate.json
```

The verifier **fails rather than relaxing a bound**, and did so twice during development:
once on an `Iv(16, 3)` interval literal (`CertificationFailure: empty interval [16, 3]`),
once on a `sin`/`cos` argument outside its proved-truncation range
(`no usable truncation for t=1120.0`). Both were fixed at the source rather than by widening
a tolerance.

### Assembled enclosure

```
pole   = 0.45715833231743186                            width 4.8e-77
arch   in [1.2100349293664769, 1.2101336300807928]      width 9.9e-05
logpi  = 1.1447298858493788 * ||f||^2                   width 5.9e-78
prime  = 0.51733913650311314                            width 1.7e-76
Q_zeta in [0.0051242393314170015, 0.0052229400457328018] width 9.870e-05
||f||^2 = 999999999999981276942897/10^24  (exact)
sup|C'| <= 11.436722214306984       C(2L) in [-8.1e-78, 8.6e-78]  (hypothesis verified)
F(14) = 4.1293580707791e-13   F'(14) = -0.17626230294801   F''(14) = +0.11323208261102521
a_f = 0.12427359776233127
```

Regression (no weight on the verdict): the previous float value `0.0051682432` and the
independent spectral value `0.0051678677` both lie inside the certified `Q_ζ` enclosure;
constant self-tests agree with libm to `0.05–0.48 ulp`.

## Files

| file | what |
|---|---|
| `CERTIFICATE_THEOREM.md` | the exact function, the object and its multiplicity check, the theorem, the four analytic inputs, the complete error budget |
| `verify_certificate.py` | the self-contained verifier: exact-rational interval arithmetic, proved transcendentals, hard gates |
| `certificate.json` | rational endpoints and status for every quantity, both certificates, and the two corrections |
| `qzeta_cache.json` | the cached `Q_ζ` enclosure (delete to force a recompute) |
| `RECONCILIATION.md` | confirmed / corrected / still unchecked, against the preceding run |
| `RESULT.md` | this file |
| `SOURCE_HASHES.txt` | md5 of every input read |

## Unchanged-file check and final state

```
md5sum -c SOURCE_HASHES.txt   ->   all 11 OK  (see below)
HEAD = cb562f7   branch = main   (unchanged)
untracked: extensions/certified-detection-witness-2026-10-08/   (preceding run, untouched)
           extensions/newtonian-weil-response-2026-10-08/       (earlier run, untouched)
           extensions/planted-witness-closure-2026-10-08/       (this run)
```

Nothing tracked was modified.
