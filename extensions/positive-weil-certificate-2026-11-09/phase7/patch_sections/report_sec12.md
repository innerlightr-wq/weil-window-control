
---

# 12. PHASE 7 (2026-11-12) — constants sharpened from the frozen blocks, conditioned coercivity, and the Connes-paper dictionary

No matrix was reassembled. `phase4/certs/block400_{even,odd}.json` were consumed
as committed; only the Rump shift and a rank-`K` rational bump were applied.
Full record, scripts, logs and artifacts:
`~/scratch/rh-connes-continuation-2026-11-12/notes/CONTINUATION_REPORT_2026-11-12.md`,
whose §10 records the corrections made after same-day external review.

## 12.1 The sharpened unconditioned constants

| quantity | committed | **2026-11-12** | factor | Rump `||I-G||_inf` |
|---|---|---|---|---|
| even-sector `a` at `T=400` | `1.42e-17` (`mu = 1.41e-17`) | **`1.43300048828e-17`** (`mu = 1.43300048827856e-17`) | `1.0163` on `mu` | `1.12012e-12` |
| odd-sector `a` at `T=400` | `2.4e-17` (`mu = 1.47e-17`) | **`1.3857e-14`** (`mu = 1.38569999999861e-14`) | **`577`** | `5.65580e-18` |

`a` is the Rump-certified shift on the finite block; `mu` is the whole-window
constant that `certify_final` then certifies, and it is **strictly below** `a`
by construction — the Schur test requires `(a-mu)(d-mu) > kappa^2`. Every `mu`
below is published as a short exact decimal with its margin; see §12.4.

New whole-window bracket at `L = 4/5`, the **even** sector still binding:

```
    1.43300048827856e-17  <=  lambda_Q(4/5)  <=  1.67312553407068186807726837736e-17
```

short quotable form `1.4330e-17 <= lambda_Q(4/5) <= 1.674e-17` (the upper
endpoint rounds **up**: `1.6731e-17` would not be a bound);
width `2.4013e-18`, down from `2.64e-18`. Parity separation margin
`1.38403e-14` (was `6.2687e-18`), ratio `mu_odd/U = 828`. The odd constant now
**exceeds** Zhu's Thm 6.2 value `8.206e-15` by a factor `1.689`; Zhu's
simplicity and evenness results are untouched and remain his.

**The `T = 400` even route is exhausted.** The certified floor
`1.43300048828e-17` and the independently measured ceiling
`R_400(f_400) = 1.4330724047355850e-17` (`phase4/certs/cross.json`) pin
`lambda_min(A_even, T=400)` to a relative width of `5.0e-5` — not the `~1.6%`
of §10.3, which measured the gap from the *committed* `mu` rather than from
`lambda_min(A)`. Confirmed a third way: inverse iteration on the midpoint block
gives `lambda_1(A_even) = 1.43307240474e-17`, matching `R_400(f_400)` to twelve
digits, and `lambda_1(A_odd) = 1.38574029573e-14`. §11.6's `T = 800` headroom is
correspondingly `+16.4%`, not `+18.3%`; **CONSOLIDATE at `T = 400`** stands,
strengthened.

The tests are not vacuous: the same routine **fails** just above each certified
shift (`1.43308715820e-17` even, `2.4576e-14` odd).

## 12.2 Certified conditioned coercivity, **by parity sector**

For any explicit vectors `v_i` and any `c >= 0`,
`A + c sum_i v_i v_i^T >= a I` implies `<A f, f> >= a ||f||^2` for every `f`
with `<v_i, f> = 0`, since the rank-`K` bump vanishes there. One Rump test on a
rank-`K` **exact-rational** perturbation of a frozen sector block therefore
certifies conditioned coercivity; by Courant-Fischer it also certifies
`lambda_{K+1}(A) >= a`. The `v_i` were chosen by inverse iteration on the
midpoint matrix — a step needing no rigour, since only the interval test is a
certificate.

| sector | conditions | certified `a` | `||I-G||_inf` | sector `mu` | also |
|---|---|---|---|---|---|
| **even** | 4 | `0.61755371093750` | `2.6216530e-32` | **`0.617553710936882`** | `lambda_5(A_even) >= 0.61755371` |
| **even** | 5 | `0.85349273681641` | `2.1380206e-25` | **`0.853492736815556`** | `lambda_6(A_even) >= 0.85349274` |
| **odd** | 3 | `0.10680580139160` | `2.1597066e-37` | **`0.106805801391493`** | `lambda_4(A_odd) >= 0.10680580` |
| **odd** | 4 | `0.78413009643555` | `9.2955279e-29` | **`0.784130096434765`** | `lambda_5(A_odd) >= 0.78413010` |

The complement constant (`d = beta - M g_N` even, `beta - M g_N - 2|p_perp|^2`
odd, the pole term carried) and the coupling
`kappa = 1.3083437712853882590e-29` are **unchanged**; only `a` changed.

**Combined over both sectors** (the form decouples by parity, and complex `f`
splits as `Q(u) + Q(v)` with each part splitting again):

| conditions | certified lower bound on `Q_zeta`, **both sectors** |
|---|---|
| 4 even + 3 odd = 7 | `Q_zeta(f) >= 0.106805801391493 ||f||^2` |
| 4 even + 4 odd = 8 | `Q_zeta(f) >= 0.617553710936882 ||f||^2` |

So the smallness of `lambda_Q(4/5)` is carried by **four even and four odd
explicit directions out of 440**. Conditioning on them raises the certified
constant by sixteen orders of magnitude. **Four conditions on the combined space
cannot do this** — min-max forbids it, since the odd sector carries three
minuscule eigenvalues of its own; the per-sector counts matter and are stated
above.

**Scope.** This proves that *some* explicit subspace of the stated dimension
works. It does **not** identify that subspace with the `E`-images of the
zero-integral prolate combinations (arXiv:2611.13129v1 Fact 5.1), which is
untested here. Removing directions changes the minimum over what remains; it
does not raise `lambda_Q(4/5)`, which stays in the §12.1 bracket. And it does
not escape the gate: `d` still needs `T > T_gate(L) = 2 pi e^{A_L}(1+o(1))`
with `A_L ~ 4 e^L`, so the route remains doubly exponentially capped in `L`.

## 12.3 The dictionary to arXiv:2611.13129v1 (Andrews-Connes-La Gatta)

`lambda = e^L` and `QW_lambda(g,g) = Q(f)` for `g(u) = f(log u)` — **proved**,
the two archimedean normalizations differing by exactly
`int_0^inf (1-e^{-s})/sinh s ds = 2 log 2`, and verified to ~30 digits by an
independent implementation of their eqs. (6)(7)(8) against
`phase5/certs/conventions.json` (`phi_0`, `phi_2`, `phi_0 - phi_2`), pole
moments included. Hence `eps_lambda(1) = lambda_Q^{even}(L)` and
`eps^odd_lambda(1) = lambda_Q^{odd}(L)` with `lambda = e^L`; our window
`L = 4/5` is `lambda^2 = 4.95303242439511480365`, just below their smallest
tabulated `lambda^2 = 5`.

Their eigenvalues are eigenvalues of compressions, so by min-max they bound the
limiting eigenvalues only from above — the authors state this (their §4.1
"Convergence in `N`", Fact 4.1) and treat their values as numerical evidence
with that limitation. What their method cannot do is **certify** the
infinite-dimensional lower bound. Combining §12.1 with the exact Theorem 2.1
rungs at `lambda = e^{4/5}` (Bouwkamp's method, App. A, validated against their
printed 30-digit `chi_0`, `chi_4` at `lambda^2 = 5`, and reproducing their
`-log delta_4 = 41.185`, `-log delta_6 = 33.833`, `-log delta_10 = 21.261`):

```
    delta_4 = 2.2391127080494187067e-18     delta_6  = 3.4175764349345855425e-15
    delta_8 = 2.4465397668381e-12           delta_10 = 9.3660403962553e-10

    6.3998587 <= eps_lambda(1)/delta_4(lambda) <= 7.4722703        (even)
                 eps^odd_lambda(1)/delta_6(lambda) >= 4.0546277    (odd)
```

against their `7.06` and `4.84` at `lambda^2 = 5`, both upper bounds. This is an
explicit certified enclosure at a stated window with a reproducible
certificate — **no priority is claimed**, and `delta_4`, `delta_6` are
high-precision and validated but **not yet interval enclosures**, which is the
one uncertified input and the gating item for every ratio quoted here.

Independent corroborations that fell out of the eigenvector step:

* the counts below `10^{-2}`: **exactly four** even (both bounds certified, see
  §12.5) and **at most three** odd — the paper's Table 3 at `lambda^2 = 5`
  records `4` and `3`;
* the ratios `lambda_j(A)/delta_{4j} = 6.400, 2.954, 2.734, 1.511` (even,
  `j = 1..4`) and `lambda_j(A)/delta_{4j+2} = 4.055, 2.590, 2.634` (odd) all sit
  in or just below the paper's ranges. **These are NUMERICAL**: `lambda_j(A)` are
  eigenvalues of a compression of `R_T`, and since `Q >= R_T` they bound
  `eps_lambda(j)` from neither side for `j >= 2`;
* their sharp test eq. (33) gives
  `log(delta_8/delta_4) - log(lambda_2/lambda_1) = 0.773044` against their
  `0.771` at `lambda^2 = 5` and their reported range `0.53-0.83` over
  `5 <= lambda^2 <= 300`;
* Fact 5.3 item 3 predicts that the odd sector, after `lambda^2 - 2 = 3`
  conditions, retains "one transitional eigenvalue `~ 0.08-0.10`" which the
  pole condition removes, "multiplying the minimum by a factor between `2.1`
  and `9` (`lambda^2 = 5`)". We certify that transitional eigenvalue at
  `>= 0.106805801391493` and the factor at `7.342`.

The cross-sector test shows the one discrepancy and it is accounted for:
`log(lambda_1(A_odd)/lambda_1(A_even)) = 6.8742` against exact
`log(delta_6/delta_4) = 7.3306`, shortfall `0.4564` versus the paper's `0.378`;
the excess `0.079` is exactly `log` of the ratio of the two sectors' `R_T`
deficits (`0.9065` even, `0.8378` odd).

Classification of the paper's contribution: Theorem 2.1 and Appendix A are
**usable** (exact rungs); Proposition 8.17's monotonicity half is
**unconditional and usable**; Theorem 8.14 + Lemma 8.16 + Proposition 8.18 are a
**conditional reduction** resting on (H2') — a hypothesis about smoothed
semilocal pair vectors — plus (R1), (R2) (the latter explicitly left open by the
authors) and an asymptotic-family hypothesis that no finite set of windows can
supply; §6-§7 (ground state, Mellin zeros) give **no transfer**, and their
Theorem 3.1/3.2 upper bound is weaker than what is already certified here.
**Nothing in the paper lower-bounds `eps_lambda` at any `lambda`.**

Theorem 8.14 is nevertheless the paper's real contribution to this project,
because it delivers the complement constant `d = (log sigma)/2` with no `T` and
no `A_L` — i.e. it is a gate-free route to the one estimate that makes this
repository's chain doubly exponentially expensive in `L`. That, and not its
numerics, is what would scale.

## 12.4 The Schur margin, published explicitly

`certify_final` requires the three strict signs `a > mu`, `d > mu`,
`(a-mu)(d-mu) > kappa^2`, so `mu` must sit strictly below `a` by at least the
coupling penalty. Its supremum satisfies `(a-mu*)(d-mu*) = kappa^2`, i.e.
`a - mu* = kappa^2/d` to leading order:

```
    a - mu*  =  1.41280e-58   at L = 4/5        (kappa^2/d = 1.4127964e-58)
             =  4.96559e-57   at lambda_-       (kappa^2/d = 4.9655891e-57)
```

reproducing `kappa^2/d` to six digits in both cases — a check on the Schur
formula itself. A `mu` published at that distance from `a` cannot be serialized:
its margin is `1e-41` of `a` relatively, far below any reasonable print width,
so a round-tripped value would be a ball straddling `a` and would **fail**
replay. Each `mu` above is therefore

```
    mu  =  a * (1 - 10^-12),  truncated DOWN to 15 significant decimal digits,
```

giving up `10^-12` of `a` relatively — about `10^30` times the true margin — so
it certifies with large room while remaining an exact short decimal that
replays bit-for-bit. All eight certificates, their margins, and the fact that
`certify_final` **fails** at `mu = a` in every case, are recorded in
`~/scratch/rh-connes-continuation-2026-11-12/certs/schur_margins.json` and
audited by `work/margin_audit.py`, which exits non-zero unless both directions
hold.

## 12.5 The even-sector count, certified both ways

The conditioned certificates bound the counts from **above** only:
`eps_lambda(5) >= 0.6175537109` gives at most four even eigenvalues below
`10^{-2}`, and `eps^odd_lambda(4) >= 0.1068058014` at most three odd. For the
lower side, `src/certify.py::Win(4,5).matrix(N)` is the Weil form on the first
`N` even cosine modes in an **orthonormal** basis, in closed form through
digamma / Lerch `Phi` / Hurwitz `zeta` as Arb balls — a **compression of `Q`**,
so by min-max `lambda_j(M) >= eps_lambda(j)` and its eigenvalues are *upper*
bounds. For an explicit 4-dimensional `W` of frozen rationals,

```
    eps_lambda(4)  <=  lambda_4(M)  <=  lambda_max(M|_W)  <  c ,
```

certified by the `4x4` test `c G - B > 0`, `G = W^T W`, `B = W^T M W`. At
`N = 160`, `prec = 2920`:

```
    eps_lambda(4)  <=  0.00311839891820  <  10^-2        CERTIFIED
```

Hence **exactly four even eigenvalues of `Q` below `10^{-2}` at `L = 4/5`**, both
bounds certified. `Win` implements the even basis only, so in the odd sector only
**at most three** is certified; the observed `lambda_j(A_odd)` are numerical.

The same compression gives certified upper bounds `eps_lambda(1..5) <=
1.6731255e-17, 8.4227418e-12, 6.4517875e-7, 0.0031183989, 0.62402164`; running
the even conditioning at `K = 1, 2, 3` would supply matching lower bounds and
bracket the first five even eigenvalues individually. Left to a later version.
