# CERTIFICATE_THEOREM.md — the extended detection interval

Extends the certified added-quartet detection interval of the paper (§6.1,
Theorem 11) from `[1/4, 49/100]` down to `[9/100, 49/100]`, in the same
four-mode space, for the same functional, with the same machinery and a sharper
error budget.

---

## 1. The object, frozen

`L = 4/5`, `gamma = 14`, both **exactly rational**. Basis: the four `EVEN_D`
modes of [`src/zeta_window.py`](../../src/zeta_window.py),

```
   phi_k(u) = cos(w_k u)/sqrt(L) ,  w_k = (2k+1)pi/(2L) = (2k+1)(5/8)pi ,  k = 0..3,
```

on `[-L,L]`, extended by zero. `w_k L = (2k+1)pi/2` exactly, so `phi_k(±L) = 0`
and each `phi_k` is continuous on `R`, real, even, supported in `[-L,L]`;
`int phi_j phi_k = delta_{jk}`. Write `f_x = sum_k x_k phi_k`.

The zero configuration is the nontrivial zeros of `zeta` **together with** four
added zeros `1/2 ± delta ± 14i`, each counted once; **no actual zero is removed
or moved**. The corresponding Weil functional is

```
   Q_add,delta(f_x) = Q_zeta(f_x) + 4[ A_delta(f_x)^2 - B_delta(f_x)^2 ] ,
   A_delta(f) = int f(u) cos(14u) cosh(delta u) du ,
   B_delta(f) = int f(u) sin(14u) sinh(delta u) du ,
```

equivalently `x^T[G + 4aa^T - 4bb^T]x`, with `4[A^2-B^2] = 4 Re F(14+i delta)^2`.
At `delta = 0`: `b = 0`, so `Q_add,0 = Q_zeta + 4F(14)^2`, the multiplicity
limit — retained as a control. **No node condition is imposed and no on-line
contribution is subtracted.** The repository also carries the count-preserving
replacement form `Q_rep,delta = Q_add,delta - 2F(14)^2` (paper eqs. (3)–(4));
that is a **different functional** and is not used here.

## 2. The theorem

> **Theorem A (certified).** With `L`, `gamma`, the basis and `Q_add,delta`
> exactly as in §1, there is a family of **69 closed rational cells** with
> disjoint interiors, consecutive and gap-free, whose union is exactly
> `[9/100, 1/4]`, and for each cell a **frozen rational** coefficient vector
> `x^(i) in Q^4` (listed in [`results/certificate_new.json`](results/certificate_new.json);
> here **one** vector suffices for all 69 cells), such that for every `i` and
> **every** `delta` in cell `i`
> ```
>    Q_add,delta( f_{x^(i)} )  <  0 .
> ```
>
> Combined with the inherited certificate —
> `Q_add,delta(f) <= -0.00245490957806631` for every `delta in [1/4, 49/100]`
> with the published frozen vector — the certified detection domain is
> ```
>    delta  in  [ 9/100 , 49/100 ]   =   [ 0.09 , 0.49 ] ,
> ```
> extending the previously certified left endpoint from `1/4` down to `9/100`.

Certified by `verify_new_certificate.py` (five assertion groups, non-zero exit
on any failure); constructed by `scripts/build_cells.py`.

The one witness is
```
 x = ( -403124850533/1000000000000, 1, 327861203727/500000000000,
       14340376641/1000000000000 ) ,
 ||x||^2 = 796343583578364461794543/500000000000000000000000 = 1.5926871671567…
```
(not normalised; every bound below is for this exact rational vector).

**Because one vector suffices this is a single-witness interval statement.** The
69 cells are a subdivision needed only because a single interval box over
`[9/100, 1/4]` inflates — `delta` occurs independently in `delta`, `delta^2`,
`cosh(delta L)` and `sinh(delta L)` — and returns a useless positive bound. The
published run recorded the same phenomenon on `[1/4, 49/100]`, where `256`
subcells were required. No claim is made that one vector works everywhere a
priori; it is certified cell by cell, **as a cell**, by interval arithmetic in
`delta` over the whole closed cell — not by sampling midpoints.

## 3. Margins and cells

Worst unit-`L^2` margin `1.276748e-06` on cell 7; best `1.915756e-04` on cell
68. Per-cell tail allowance `1.7364e-09` throughout. 68 of 69 cells reuse the
previous cell's witness; 0 cells needed halving during construction.

| cell | `delta` range | `Q_add` upper bound | unit-`L^2` margin | tail allowance | margin / allowance |
|---|---|---|---|---|---|
| 0 | `[9/100, 361/4000]` | `-1.087604e-05` | `6.8287e-06` | `1.736e-09` | `6263` |
| 1 | `[361/4000, 453/5000]` | `-1.046039e-05` | `6.5678e-06` | `1.736e-09` | `6024` |
| 2 | `[453/5000, 9109/100000]` | `-9.879781e-06` | `6.2032e-06` | `1.736e-09` | `5690` |
| 3 | `[9109/100000, 1434/15625]` | `-9.069444e-06` | `5.6944e-06` | `1.736e-09` | `5223` |
| 7 | `[5997709/62500000, 6153959/62500000]` | `-2.033461e-06` | `1.2767e-06` | `1.736e-09` | `1171` |
| 66 | `[15216459/62500000, 15372709/62500000]` | `-2.768097e-04` | `1.7380e-04` | `1.736e-09` | `159411` |
| 67 | `[15372709/62500000, 15528959/62500000]` | `-2.835662e-04` | `1.7804e-04` | `1.736e-09` | `163302` |
| 68 | `[15528959/62500000, 1/4]` | `-3.051201e-04` | `1.9158e-04` | `1.736e-09` | `175715` |

(the first cells, the worst-margin cell, and the last cells; all 69 with their
witnesses are in `results/certificate_new.json`. The margin/allowance column
compares the raw `Q_add` upper bound with the raw tail half-width.)

## 4. The join, and the error budget

**The join at `1/4` is certified from both sides** (`verify_new_certificate.py`,
assertion group 3): the last new cell's witness is still negative *at*
`delta = 1/4` (upper bound `-3.335481e-04`) and the inherited published frozen
vector is negative there too (upper bound `-2.521632e-03`). So the two certified
sets genuinely overlap at the shared endpoint; the union is not a splice across a
gap.

Exact-rational interval arithmetic with outward rounding on a `2^-260` dyadic
grid; **no floats and no libm in the certified path**. Inherited unchanged from
the published verifier: the explicit formula in `psi`-form, the archimedean
elementary-integral reduction, the `gamma_E` elimination, the proved prime-power
support cutoff (`log 5 > 2L`, so exactly `n = 2,3,4`), every transcendental with
its proved remainder, and the exact quartet closed forms.

**One new analytic ingredient**, Proposition T of
[`ANALYTIC_TAIL.md`](ANALYTIC_TAIL.md): the archimedean tail is bounded to
`O(1/M^2)` instead of `O(1/M)` at the same truncation `M = 2^17-1`, by two
integrations by parts. Effect at the published frozen vector:

| | published | this certificate |
|---|---|---|
| tail allowance (half-width) | `±4.935033776e-05` | `±1.6500e-09` |
| `Q_zeta` enclosure width | `9.870e-05` | `3.339e-09` |

and the new enclosure is **contained in** the published one — the decisive
consistency check (assertion group 5). The `gamma_E` term `T` contributes width
`3.881e-11`; all transcendentals and the quartet contribute `< 1e-70`.

> **The limiting error source for every cell is the `O(1/M^2)` archimedean tail
> remainder `sup|C''|/(4(M-3/4)^2)`, and nothing else.** There is no floating
> precision anywhere in this code, so no amount of precision could have
> substituted for the analytic improvement. The published `4.94e-05` allowance
> cannot see the `-2.19e-06` minimum eigenvalue at `delta = 9/100`; the new
> `1.65e-09` allowance resolves it with a factor `>1000` to spare.
>
> **How stark this is**, on the new witness: the published `O(1/M)` budget gives
> `Q_zeta ∈ [-3.265e-05, +9.977e-05]`, which does not even determine the
> **sign**; the `O(1/M^3)` budget gives `[2.444734052852e-05,
> 2.444734443683e-05]`.

What now limits extension to smaller `delta` is **not** the error budget but the
vanishing of the margin at the four-mode detection onset — see
[`NON_DETECTION.md`](NON_DETECTION.md) §3.

## 5. Independent cross-check, and exactly how far it reaches

The repository's separately coded implementation
([`independent_point_check.py`](../planted-witness-consolidation-2026-10-08/independent_point_check.py),
sha256 `bb72086e…`) evaluates the archimedean term through a **completely
different closed form** — digamma/trigamma at `z_k = 1/4 - i w_k/2`, with the
Euler constant cancelling identically and a `g`-series/midpoint error budget —
instead of `131071` Laplace integrals with a `1/M^j` tail. Run on the **new**
witness (`scripts/indep_crosscheck.py`):

| `delta` | independent `Q_add` upper bound | primary `Q_add` enclosure | both negative |
|---|---|---|---|
| `1/5` | `-1.067391765e-04` | `[-2.025882246e-04, -2.025882207e-04]` | yes |
| `1/4` | `-2.377002990e-04` | `[-3.335493471e-04, -3.335493432e-04]` | yes |

and on `Q_zeta` of the new witness the primary `O(1/M^3)` enclosure
`[2.444734052852e-05, 2.444734443683e-05]` is **contained in** the independent
`[8.689196092281e-06, 1.202963886872e-04]`. Both enclose the same real number,
the sharper inside the looser, and the independent route reaches a negative
upper bound of its own — not merely matching intermediates.

**What is shared, stated plainly:** the frozen mathematical input (`x, L,
gamma`); the classical explicit formula; the canonical `(A_k, B_k)` reduction of
the autocorrelation (recoded there, not imported); and the idea of exact
rational interval arithmetic (re-implemented there). This is a **second
implementation, not external peer review.**

**What it does NOT reach.** Its budget is about `1.1e-04` on `Q_zeta` here —
roughly `90x` the worst margin of this certificate (`1.28e-06`, at the left
endpoint). So it corroborates the **right** part of `[9/100, 1/4]` only. **The
new left endpoint `delta = 9/100` rests on the primary implementation alone**,
and closing that would need either the sharper tail ported into the independent
route or a third implementation. Neither was done.

## 6. What is NOT claimed

* **Not** `Q_zeta(f_x) < 0`. For the published frozen vector
  `Q_zeta in [+0.005167865935080, +0.005167869273874] > 0`; for the new witness
  `Q_zeta in [+2.444734052852e-05, +2.444734443683e-05] > 0`. It is the added
  quartet that makes the total negative. See
  [`FOUR_MODE_POSITIVITY.md`](FOUR_MODE_POSITIVITY.md) for what positivity on
  this subspace does and does not mean.
* **Nothing about actual zeros of `zeta`.** `gamma = 14` is a chosen rational
  height, **not** a verified ordinate; the quartet is **hypothetical**, added to
  the zero multiset with nothing removed. A negative Weil functional for a
  fabricated configuration is the *detection direction* of the classical
  criterion and locates nothing. **No claim that RH has been solved, weakened,
  or established conditionally.**
* **No all-height or all-window statement.** `L = 4/5` and four modes
  throughout. The repository's own measurement — `||P_E r_gamma||^2` falling from
  `3.60e-2` at `gamma = 14` to `1.13e-7` at `gamma = 280` — is the reason to
  expect no uniformity in height.
* **No threshold, no monotonicity, no uniqueness, no sharpness.** `9/100` is
  where the construction stopped, not a boundary;
  [`NON_DETECTION.md`](NON_DETECTION.md) §4 states the open gap.
* **Not Lean verification.** A computer-assisted certificate run by CPython; no
  proof assistant, nothing formalised, no external peer review.
* **No novelty** for the rank-one/Rayleigh criterion, Sherman–Morrison, the
  Schur complement, interval arithmetic, polarization, or any coordinate change
  — see [`METHODS_AND_LIMITATIONS.md`](METHODS_AND_LIMITATIONS.md) §2.
