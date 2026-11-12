
---

# 15. PHASE 7d (2026-11-12) — a rigorous comparison of the two windows' minima

With `eps_0 = eps(e^{4/5})` and `eps_- = eps(lambda_-)`, pairing the lower
endpoint at each window with the upper endpoint at the other gives

```
    mu_- / U_0   <=   eps_- / eps_0   <=   U_- / mu_0 .
```

This uses **only certified endpoints** — no prolate denominator, no parameter
transfer, no monotonicity transport, nothing from arXiv:2611.13129v1. Taking all
four inputs as the exact decimals that are packaged (published `mu`, which round
down; published `safe_decimal` upper bounds, which round up):

```
    mu_-/U_0 = 8.27343868254787e-18 / 1.674e-17  = 0.494231701466...
    U_-/mu_0 = 9.475e-18 / 1.43300048827856e-17  = 0.661200053838...
```

hence, outward-rounded,

```
    0.494  <=  eps(lambda_-) / eps(e^{4/5})  <=  0.662 ,
```

i.e. **the true spectral minimum falls by between `33.8 %` and `50.6 %`**
between `lambda^2 = 4.95303242439511480365` and
`lambda^2 = 4.999999999999999927`. Consuming the full-precision endpoints
instead tightens it to `[0.494490, 0.661184]`, same 3-decimal outward form.

This is a statement about the **true minima** of two infinite-dimensional forms,
with every input certified and every division taken in the safe direction. It is
stronger than the floor-vs-rung diagnostic of §13.3, and it is the first
quantitative statement here relating the window Weil form at two windows.

**It does not establish any scaling law.** The rung ratio
`delta_4(5)/delta_4(4.953) = 0.579867` lies inside the bracket, at `51.2 %`
across. That is a **numerical rung-ratio reference**, not a prediction of
Hypothesis 3.4: the hypothesis compares `eps_lambda(j)` with `delta_{4j}` up to
constants `c <= C` and a factor `log lambda`, and says nothing about the ratio of
two Weil minima equalling the ratio of their rungs — that needs the extra
assumption that the proportionality factor is the same at both windows.
Narrowing the bracket below `0.058` (each window's own bracket ratio below about
`1.05`, against `1.1676` at `L = 4/5` and `1.1452` at `lambda_-`, i.e. roughly a
threefold narrowing) would test **that stronger approximation**, not the
hypothesis.

## 15.1 Two qualifications, retained

* **The comparison with `9.1663e-18` is numerical consistency only.** That
  figure is a **reconstructed finite-compression estimate** — the lowest
  eigenvalue of the paper's `N = 960` log-Fourier compression, recovered from a
  printed `-log eps = 39.231` — lying inside our bracket at a nearby window. It
  is **not** an independently known value of the infinite-dimensional minimum,
  and it is never used as a bound.
* **At exactly `lambda^2 = 5` there is only the transported upper bound.**
  `eps_lambda <= 9.475e-18` for every `lambda >= lambda_-` (Prop. 8.17). The
  lower bound runs the other way: `eps_lambda >= 8.27343868254787e-18` for every
  `lambda <= lambda_-`, i.e. on the whole window range
  `0 < L <= 0.80471895621705018` — an interval, not a point — but it does
  **not** extend past that endpoint, and in particular not to `lambda^2 = 5`.

## 15.2 Release acceptance

`~/scratch/rh-connes-continuation-2026-11-12/work/acceptance.py` applies this
patch to a worktree at `7515396` and re-verifies **every** release claim from
the packaged artifact files rather than from values held in memory: the 6 patched
files and all 176 repository checks on the patched tree; the 8 Schur margins
re-certified from their own stored `a`, `d`, `kappa`, `mu` (with `mu = a`
required to fail); the 4 conditioned certificates re-certified by rebuilding the
rank-`K` bump from the stored frozen rational vectors; the 2 `lambda_-` blocks
from their packaged checkpoints; the `lambda_-` upper bound from its stored
exact rationals; the two-window ratio from the four packaged decimals; and the
even-sector count certificate. Stage **B2** separately re-runs the two
**improved** unconditioned `L = 4/5` block bounds
(`A - 1.43300048828e-17 I > 0`, `A - 1.3857e-14 I > 0`) from the frozen blocks,
which stage A does **not** cover: `phase4/verify_t400.py` tests the committed
`1.42e-17` / `2.4e-17` shifts, and `certify_final` assumes `A >= a I` rather than
checking it. It reverts the patch, checks the worktree clean, and exits non-zero
on any failure, so it is usable as a release gate.
**Result: PASS, 42 checks, 0 failed.**
