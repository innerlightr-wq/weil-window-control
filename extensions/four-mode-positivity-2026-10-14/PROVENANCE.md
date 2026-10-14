# PROVENANCE.md

## 1. What is imported, and from where

Everything in the certified path except the tail bounds is the repository's
published code, **imported, not copied and not edited**:

| file | sha256 | role |
|---|---|---|
| [`extensions/planted-witness-closure-2026-10-08/verify_certificate.py`](../planted-witness-closure-2026-10-08/verify_certificate.py) | `c49858a1d50afbfefaaa63a4914c0c8bd2eb88269c21b4dd15e24b03deb510f0` | the whole certified path: interval class, every proved transcendental, canonical autocorrelation, pole / archimedean / `log pi` / prime terms, `gamma_E` elimination, support cutoff, exact quartet closed forms |
| [`extensions/planted-witness-consolidation-2026-10-08/independent_point_check.py`](../planted-witness-consolidation-2026-10-08/independent_point_check.py) | `bb72086ecb6e41c7450ac2e685149942fdff9785652de4744e81f5bc164ef4a8` | the separately coded digamma/trigamma implementation, for the cross-check |
| [`extensions/planted-witness-closure-2026-10-08/CERTIFICATE_THEOREM.md`](../planted-witness-closure-2026-10-08/CERTIFICATE_THEOREM.md) | `58ab96b98509a06257b7f60b1836607404ed1decf2486d2cecb2651655a47866` | the published theorem statement and its §5.2 tail bound |
| [`extensions/planted-witness-closure-2026-10-08/certificate.json`](../planted-witness-closure-2026-10-08/certificate.json) | `26639c9110b2c425aacf76f528cf6aae6cbaabfc9f5bb923236534228de7304c` | the frozen published vector and enclosures |
| [`extensions/planted-witness-consolidation-2026-10-08/interval_cells.json`](../planted-witness-consolidation-2026-10-08/interval_cells.json) | `2d5491a9ac534b7c32765a220d762c9aac39cca4496925f5f6147acafec9a84c` | the published 256 cells on `[1/4, 49/100]` |
| [`extensions/planted-witness-consolidation-2026-10-08/ANALYTIC_AUDIT.md`](../planted-witness-consolidation-2026-10-08/ANALYTIC_AUDIT.md) | `35b27b06a3920fafc0b49a61d7668e58650ca716ba765d5f7c2f3388b4711562` | the audit of the published analytic steps |
| [`extensions/planted-witness-consolidation-2026-10-08/REPRODUCE.md`](../planted-witness-consolidation-2026-10-08/REPRODUCE.md) | `a5148cae15b8ff99fca249e8d04b5ea43dbfd221750e9384d33e9424c88a4f5f` | the published reproduction instructions |
| [`notes/proofs.md`](../../notes/proofs.md) | `09b887b605f9d5e47fb084a68965015314b82ee380fded6f5b473b78a2c83d87` | §B.3 Taylor coefficient, §B.5 regime separation |
| [`src/zeta_window.py`](../../src/zeta_window.py) | `5d8676e4efefb545c7bb2d73efaacf54011fbbe4bf4d9be3cf7a5db5863ca272` | the `EVEN_D` basis definition |

`scripts/weil_core.py` resolves the first of these **relative to its own
location** — no absolute paths, so a fresh clone works unchanged — and
**verifies its sha256 at import time**, raising if it differs. The same is done
for the second in `scripts/indep_crosscheck.py`. So the published file cannot be
silently substituted, and this extension cannot run against a modified copy.

What `weil_core.py` adds is only:

1. the coefficient vector `x` becomes a parameter instead of a module constant;
2. the published archimedean loop is rearranged into its `x`-**independent**
   part — `laplace` is linear in `(A,B)`, so
   `sum_{m<M} laplace = sum_k [A_k SA_k + B_k SB_k]` exactly. An algebraic
   rearrangement of the published loop, **not a new bound**;
3. the resulting split `Q_zeta = Qcore + tail` is exposed, so the tail can be
   replaced by the sharper bounds of [`ANALYTIC_TAIL.md`](ANALYTIC_TAIL.md).

Validated against the published run by `tests/test_core_vs_published.py` (16
checks): `pole`, `logpi`, `prime` to all printed digits; `arch` and `Q_zeta`
enclosures **containing** the published ones; the published tail allowance
`4.935033776e-05`; the point bound `Q_add,2/5 <= -0.0145107657498227` with the
recomputed `hi` equal to the published value; the node residual `4.129358e-13`;
`b(0) = 0` exactly; and the multiplicity limit `Q_add,0 = Q_zeta + 4F(14)^2`.
These are reproduced **as inequalities from frozen inputs**, not read as stored
labels.

## 2. Frozen data in this directory

| file | contents |
|---|---|
| `results/certificate_new.json` | Theorem A: 69 closed cells on `[9/100, 1/4]`, the one frozen rational witness, and each cell's exact enclosures and margins |
| `results/route_o_step1.json` | Theorem O1: the exact rational `Gt_mid` and radii, `sigma3`, `Kappa3`, `rho3`, `t_cert`, the LDL pivots, `m0`, and the `S_1`, `S_3` enclosures |
| `results/route_o_step2.json` | Theorem O2: 6072 closed cells on `[0, 3831/64000]` as `(lo, hi, w, tau, tau')` only — `N` and `Loss` are **recomputed** by the verifier, never read |
| `results/theta_table.json` | the `Theta` diagnostic of [`NON_DETECTION.md`](NON_DETECTION.md) §3, for both `G` and the withdrawn `G_core` |
| `results/indep_crosscheck.json` | the independent implementation's enclosures on the new witness, with its truncation parameters, and the comparison against the primary route |
| `results/arch_sums_cache.json` | the one cached object (see §4) |

## 3. Error budget, per statement

| quantity | limiting error source | size |
|---|---|---|
| `Q_zeta`, published route | archimedean tail `[(3/4)C0 + (1/2)sup\|C'\|]/(M-3/4)`, `M = 2^17-1` | `±4.935033776e-05` |
| `Q_zeta`, Theorem A | the `O(1/M^2)` remainder `sup\|C''\|/(4(M-3/4)^2)` | `±1.6500e-09` |
| `Q_zeta`, Theorems O1/O2 | the `O(1/M^3)` remainder `2 sup\|C'''\| S_4bound` | `~2e-14` |
| `gamma_E` elimination `T`, published | `[1/(2(M+1)) - 1/(6M^2), 1/(2M)]` | width `3.881e-11` |
| `gamma_E` elimination `T`, sharpened | Lemma T | width `1.357e-16` |
| polarization radius of `G` | row-sum bound, `sigma3` | `8.5390e-12` |
| transcendentals, quartet | proved remainders, outward rounding | `< 1e-70` |
| independent cross-check | `JG` g-series tail + trigamma midpoint | `~1.1e-04` on `Q_zeta` here |

There is **no floating-point arithmetic anywhere** in the certified path: exact
rational interval arithmetic with outward rounding on a `2^-260` dyadic grid
(the independent route uses its own `2^-200` grid). Raising "precision" is
meaningless here; the budget is analytic, which is why it had to be improved
analytically.

## 4. Caching

One cached object: the `x`-independent archimedean sums `SA_k, SB_k`, in
`results/arch_sums_cache.json`. Its key records `M`, the interval grid depth in
bits, `L`, `gamma`, the dimension, the basis string, the quantity, **and the
sha256 of the published verifier**; on any mismatch it is recomputed (~4 min).
No other cached object enters any proof path, and no stored verdict is read —
each verifier recomputes every enclosure it asserts and exits non-zero on
failure.

## 5. Environment

```
Python 3.14.7 (main, Aug 10 2026, 07:46:56) [GCC 16.1.1 20260728]
third-party packages:  NONE — mpmath, numpy, sympy, scipy all MISSING,
                       and none was installed
```

The repository's `requirements.txt` pins `mpmath`, `numpy`, `matplotlib`; **none
of them is present on the machine where this extension was produced, and nothing
here imports any of them.** Consequently there was no floating-point route
available at any stage, including exploration: every number in this directory is
exact rational interval arithmetic in the CPython standard library.

Resources: one worker, no nested workers, `RLIMIT_AS = 4 GiB` set in every
script, peak RSS `16–20 MiB` throughout (the matrices are `4x4`; the cost is the
`131071`-term rational sum, not memory).

## 6. External sources

**No external paper is a mathematical dependency of any statement in this
directory.** Nothing is imported, quoted as a hypothesis, or relied on from any
preprint; in particular no external positivity bound is used — Theorem O1 is
derived from the explicit formula alone. The classical facts used are listed in
[`METHODS_AND_LIMITATIONS.md`](METHODS_AND_LIMITATIONS.md) §1(a) and are
textbook.

## 7. Repository state when this extension was prepared

Prepared in a separate clean clone, so the author's working checkout was never
written to. Verified before and after: `HEAD =
725f2040c32de7bc5f757dd1f7ef2d212713b8fb` on branch `main`, with its 18
untracked paths unchanged. No reproduction command was run inside that checkout,
so no tracked certificate was regenerated or overwritten there, and no earlier
frozen certificate in this repository is modified by this extension — the
published `[1/4, 49/100]` interval and the `delta = 2/5` point bound are
**recomputed and asserted**, never edited.
