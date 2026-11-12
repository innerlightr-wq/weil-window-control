
---

# 14. PHASE 7c (2026-11-12) — a two-sided bracket at `lambda_-`, and what transports

The quadrature-free upper route (`phase5/src/upper_highN.py`'s two-stage
structure, with `src/certify.py::Win` and `L` parameterized) run at `L~`.
`N = 160` even cosine modes; exploratory `mpmath` stage locates the vector and
certifies nothing; the frozen exact rationals are enclosed in Arb inside
`owned(prec = 2920)`. **No quadrature anywhere**, so there is no R14-type
remainder to bound.

```
    [CERTIFIED]  Q(f_160)/||f_160||^2   = [9.47476951083049032849709e-18 +/- 8.10e-43]
    [CERTIFIED]  outward upper endpoint =  9.47476951083049032849709e-18
    safe short decimal (rounds UP)      =  9.475e-18
```

## 14.1 The bracket

```
    8.27343868254787e-18  <=  eps(lambda_-)  <=  9.475e-18 ,
        lambda_-^2 = 4.999999999999999927 ,  even sector
```

width `1.2013e-18`, ratio `1.1452`. The second window this repository brackets
from both sides.

| window | certified lower | certified upper | width |
|---|---|---|---|
| `lambda^2 = 4.95303242439511480365` (`L = 4/5`) | `1.43300048827856e-17` | `1.674e-17` | `2.40e-18` |
| `lambda^2 = 4.999999999999999927` (`L~`) | `8.27343868254787e-18` | `9.475e-18` | `1.20e-18` |

## 14.2 The paper's value lies inside it

arXiv:2611.13129v1 reports `eps_lambda(1)` at `lambda^2 = 5`, `N = 960`, as
`-log eps = 39.231`, i.e. `9.1663e-18` (equivalently `7.06 x delta_4`):

```
    8.2734e-18   <   9.1663e-18   <   9.4748e-18          (74.3 % of the bracket)
```

A **reconstructed finite-compression estimate** — the lowest eigenvalue of their
`N = 960` log-Fourier compression, recovered from a printed `-log eps` — falls
inside a two-sided certificate at a window differing by `7.3e-18` in `L`. Their
basis, normalization and codebase are all different from ours, which is what
makes the agreement informative; but this is **numerical consistency only**, not
an inequality between certified objects, and their figure is **not** an
independently known value of the infinite-dimensional minimum.

## 14.3 What transports upward, and what does not

By Prop. 8.17 (`eps_lambda` non-increasing in `lambda`, unconditional), the
**upper** bound transports to every larger window:

```
    eps_lambda  <=  9.475e-18     for every lambda >= lambda_- ,
                                  INCLUDING lambda^2 = 5 exactly.
```

A rigorous upper bound at the paper's own window, with no parameter transfer
and no reliance on their numerics; `3.37 %` weaker than their reconstructed
`9.1663e-18`, which is the right side to be on.

**The lower bound does not transport upward.** `mu(lambda_-) = 8.2734e-18` says
nothing about `eps` at `lambda^2 = 5` or beyond — monotonicity runs the other
way. A certified *lower* bound at exactly `lambda^2 = 5` still needs the
`arb`-window reparameterization of §13.4 item 1.

## 14.4 Ratio comparison, two-sided — NUMERICAL

With `delta_4(sqrt 5) = 1.2983880087183625604e-18`:

```
    eps_lambda(1)/delta_4(lambda)  in  [ 6.3720849 , 7.2973329 ]        (paper: 7.06)
```

two-sided now, and containing `7.06`. **A numerical comparison, not an
established enclosure**: `delta_4` is high-precision rather than
interval-enclosed, and it is evaluated at `lambda^2 = 5` while the bracket is
certified at `lambda_-^2 = 5 - 7.3e-17`. The `L = 4/5` pair
`[6.3998587, 7.4722703]` has the matching window but the same missing
enclosure. Interval-enclosing `delta_{4j}` remains the gating item for every
ratio in §§12-14.

A last diagnostic: the two certified **upper** bounds scale as
`U(lambda_-)/U(e^{4/5}) = 0.56629` against the rung ratio
`delta_4(5)/delta_4(4.953) = 0.57987`, a `2.3 %` agreement — same status as
§13.3's `0.43 %`, i.e. a diagnostic.
