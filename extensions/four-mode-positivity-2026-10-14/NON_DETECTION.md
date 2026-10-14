# NON_DETECTION.md — certified absence of a four-mode witness at small `δ`

The complement of the detection certificate. Detection needs **one** vector with
`Q_add,delta(f_x) < 0`; non-detection needs **every** vector to satisfy
`Q_add,delta(f_x) >= 0`, so it is a matrix statement and it needs Theorem O1 of
[`FOUR_MODE_POSITIVITY.md`](FOUR_MODE_POSITIVITY.md) as input.

**Read §4 before quoting anything from this file.** The certified interval
reached here is **shorter than the interval that was targeted**, and the
numerical onset is **not** a certified threshold.

---

## 1. The criterion

With `a(delta)_k = int phi_k(u) cos(14u) cosh(delta u) du`,
`b(delta)_k = int phi_k(u) sin(14u) sinh(delta u) du`,

```
   Q_add,delta(f_x) = x^T M(delta) x ,
   M(delta) = G + 4 a a^T - 4 b b^T = H(delta) - 4 b b^T ,  H = G + 4 a a^T .
```

> **Proposition 1 (classical; Schur complement / Sherman–Morrison /
> Birman–Schwinger).** If `H > 0` then `M = H - 4bb^T` has a negative direction
> **iff** `Theta := 4 b^T H^{-1} b > 1`; the maximiser is `x_* = H^{-1}b`.

No novelty is claimed for Proposition 1 — it is the same rank-one threshold
already recorded in the author's earlier local notes, reused verbatim, not
renamed. Theorem O1 supplies the hypothesis `H > 0` unconditionally.

## 2. Theorem O2 — the certified cells

For a closed cell `[d1, d2]` of width `w = d2 - d1`, write `A_1 = a(d1).x`,
`A_d = a(delta).x` and likewise for `b`. Since
`a(delta) - a(d1) = int_{d1}^{delta} a'(s) ds`, we get `|A_d - A_1| <= w Ahat ||x||`
with `Ahat := sup_cell ||a'||`, hence
`|A_d^2 - A_1^2| <= 2 w Ahat ||x|| |A_1| + w^2 Ahat^2 ||x||^2`, and by AM–GM
(`2uv <= tau u^2 + v^2/tau`, any `tau > 0`)

```
    4 A_d^2  >=  4(1-tau ) A_1^2  -  4 w^2 Ahat^2 (1 + 1/tau ) ||x||^2 ,
   -4 B_d^2  >= -4(1+tau') B_1^2  -  4 w^2 Bhat^2 (1 + 1/tau') ||x||^2 .
```

**The AM–GM step is what makes this feasible.** It charges the cell-induced loss
against the `4A_1^2` term *already present* in `M`, instead of against
`||a|| ||a'||`: `||a|| = 0.88`, while the component of `a` along the near-null
direction of `G` is `~4e-05`. Bounding at the matrix level instead
(`8||a_1|| eta_a`) forces cell widths `~2.6e-09` and `~2.5e+07` cells — four
orders of magnitude worse, and infeasible on this machine.

So for every `delta` in the cell and every `x`,

```
   x^T M(delta) x  >=  x^T N x  -  Loss ||x||^2 ,
   N    = Gt_mid + 4(1-tau) a_1 a_1^T - 4(1+tau') b_1 b_1^T          (exact rational),
   Loss = (sigma3 + rho3) + 4 w^2 [ Ahat^2 (1+1/tau) + Bhat^2 (1+1/tau') ] + tiny,
```

`tiny` covering the `~1e-70` interval widths of `a_1, b_1` at the rational `d1`.
The derivative bounds are elementary and rigorous, using `sinh(t)/t` increasing
(so `|sinh(delta u)| <= |u| sinh(delta L)/L` for `|u| <= L`),
`int|phi_k| u^2 du <= (2/3)L^{5/2}` and `int|phi_k| |u| du <= L^{3/2}`:

```
   Ahat <= 2 (2/3) L^{3/2} sinh(d2 L) ,        Bhat <= 2 L^{3/2} cosh(d2 L) .
```

Certifying `lambda_min(N) >= Loss` by exact rational `LDL^T` on `N - Loss I`
then gives `M(delta) >= 0` on the whole closed cell.

> **Theorem O2 (certified).** There is a family of **6072 closed rational cells**
> with disjoint interiors, consecutive and gap-free from `0`, whose union is
> exactly
> ```
>    [ 0 , 3831/64000 ]  =  [ 0 , 0.05985937… ] ,
> ```
> listed in [`results/route_o_step2.json`](results/route_o_step2.json), such that
> for every `delta` in that interval
> ```
>    Q_add,delta(f_x)  >=  0      for EVERY x in the four-mode subspace,
> ```
> equivalently `M(delta) >= 0`. **No negative witness exists in this subspace
> there** — certified non-detection.

The stored data holds only `(lo, hi, w, tau, tau')` per cell; `N` and `Loss` are
**recomputed** by `verify_route_o.py` (assertions `O3`, `O4`) from those and from
`results/route_o_step1.json`, so nothing in the file is trusted.

## 3. The numerical onset — a diagnostic, not a theorem

`Theta(delta) = 4 b^T H^{-1} b` evaluated at the **midpoint** of the certified
interval matrix for `G` and at the midpoints of `a(delta), b(delta)`. This is
midpoint arithmetic: **a numerical diagnostic, not a certificate.** The
certificates are Theorem O1 and Theorem O2.

| `delta` | `Theta` (true `G`) | inertia of `M` | `Theta` (`G_core`, **withdrawn**) |
|---|---|---|---|
| `13/200` = 0.0650 | `0.906281` | `(0,0,4)` | `0.552010` |
| `133/2000` = 0.0665 | `0.951002` | `(0,0,4)` | `0.578432` |
| `27/400` = 0.0675 | `0.981512` | `(0,0,4)` | `0.596415` |
| `137/2000` = 0.0685 | `1.012584` | `(1,0,3)` | `0.614694` |
| `7/100` = 0.0700 | `1.060257` | `(1,0,3)` | `0.642672` |
| `17/200` = 0.0850 | `1.611474` | `(1,0,3)` | `0.960296` |
| `9/100` = 0.0900 | `1.827636` | `(1,0,3)` | `1.082008` |
| `1/10` | `2.314318` | `(1,0,3)` | `1.350467` |
| `1/5` | `14.604889` | `(1,0,3)` | `6.439811` |
| `1/4` | `34.910858` | `(1,0,3)` | `11.517440` |

`Theta` crosses `1` and the inertia of `M` flips **together**, as Proposition 1
requires; `H > 0` at every row, by exact `LDL`. So, numerically,

```
   onset bracket  (27/400, 137/2000) = (0.0675, 0.0685) .
```

**This is a crossing bracket of a midpoint computation. It is not a certified
sharp threshold.** Monotonicity of `Theta` in `delta` is neither assumed nor
proved, no uniqueness of the crossing is claimed, and the values themselves are
not rigorous enclosures. Data: [`results/theta_table.json`](results/theta_table.json),
produced by `scripts/theta_table.py`.

The last column is kept on purpose: it is the **withdrawn** `G_core` computation
whose crossing at `(0.085, 0.09)` was reported in an earlier round of this work.
See [`FOUR_MODE_POSITIVITY.md`](FOUR_MODE_POSITIVITY.md) §4.1.

## 4. What was targeted, what was reached, and what is left open

Stated plainly, because the shortfall matters more than the result.

* **Targeted:** certified non-detection up to `13/200 = 0.065`, just under the
  corrected numerical onset.
* **Reached:** `3831/64000 = 0.05985937…`, with 6072 cells, when the run hit its
  wall-clock budget (1550 s) with `complete = False`. The cell widths shrink as
  `delta` grows — `Loss` grows like `w^2 cosh^2(delta L)` while the margin
  `lambda_min(N)` falls toward zero at the onset — so the remaining `0.0052` is
  the expensive part, not the cheap part.
* **Therefore the gap `(0.05986, 0.09)` between certified non-detection and
  certified detection is NOT closed**, and the numerical onset bracket
  `(0.0675, 0.0685)` lies **inside** that gap, uncertified from both sides.
* **Not claimed, anywhere:** a sharp threshold; a certified onset; that `0.065`
  was reached; monotonicity; that the gap is small *because* the onset is in it.
* **Open, and recorded as ongoing work rather than as a result:** pushing the
  certified non-detection endpoint up to and past the onset bracket. It is a
  matter of compute at the present cell widths, not of a missing idea, but it
  was not done, so it is not claimed.

## 5. Scope of Theorem O2

* It concerns **only** the four-dimensional subspace spanned by the four
  `EVEN_D` modes at `L = 4/5`, and only the additive quartet at `gamma = 14`.
  **A larger basis could detect at smaller `delta`; nothing here excludes
  that.** Non-detection in a finite basis is not positivity — see non-claim 4 of
  the repository README.
* It is **not** a statement about `zeta`. `Q_add,delta >= 0` on a subspace for a
  **fabricated** zero configuration says nothing about the location of any
  actual zero.
* The quartet is artificial: four zeros at `1/2 ± delta ± 14i` are **added** to
  the true zero multiset and none is removed. `gamma = 14` is a prescribed
  rational height, not an observed ordinate of `zeta`.
