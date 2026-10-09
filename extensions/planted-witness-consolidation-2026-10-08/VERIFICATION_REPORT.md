# Verification report: the added-quartet witness, independently recomputed

**Date** 2026-10-08 · **Verdict** `VERIFIED AND CONSOLIDATED`

Nothing below was taken from a cached value. The certificate was rebuilt cold from the
frozen rational input, and the decisive inequality was then re-derived and re-evaluated by a
**separately implemented** program that uses a different closed form, a different error
budget and its own interval arithmetic.

---

## 1. What was recomputed, and what was independently cross-checked

### 1.1 Cold recomputation (`verify_certificate_consolidated.py`)

Run from an empty directory with **no cache present**:

```
real 5m25.343s      EXIT=0      all four gates PASS
Q_zeta step alone: 293.8 s
```

All four gates passed and the enclosures reproduced the previously cached ones exactly.
Three defects were found and fixed while rebuilding:

| # | defect | status |
|---|---|---|
| 1 | the cache key omitted `gamma`, so changing the detection height would have silently reused a stale `Q_ζ` | **fixed**; key is now `(x, xd, M, grid, L, gamma, ndim, src)` |
| 2 | the δ-interval certificate was emitted without machine-checkable cell data | **fixed**; `interval_cells.json` now records all 256 cells with exact rational endpoints, plus gap and endpoint assertions that raise |
| 3 | my own description said adjacent cells "share an endpoint overlap of 1/3200" | **corrected**: they *abut*, `hi(I_j) = (803+3j)/3200 = lo(I_{j+1})`, overlap exactly 0 |

Cell tiling re-verified independently of the verifier: `lo(I_0) = 1/4`, `hi(I_255) = 49/100`,
zero gaps, uniform width `3/3200`.

### 1.2 Separate implementation (`independent_point_check.py`)

```
real 3m17.482s      EXIT=0      END-TO-END NEGATIVE UPPER BOUND: True
```

Standard library only; it does **not** import the verifier and does **not** read
`certificate.json`, `cold_certificate.json` or any cache. It differs where it matters:

* **archimedean term** — the verifier sums `131071` Laplace integrals and bounds the tail by
  `[¾C(0) + ½ sup|C'|]/(M − ¾)`. The independent program evaluates the *same* series in
  closed form through `ψ` and `ψ'` at four complex points (`ANALYTIC_AUDIT.md` §5.1), where
  **`γ_E` cancels identically** — no harmonic number, no `1/M` tail. Its error budget is a
  `g`-series tail for `ψ` and a midpoint-rule bound for `ψ'`.
* **quartet term** — the verifier uses the `cosh`/`sinh` split
  `A_δ = ∫f\cos(14u)\cosh(δu)`, `B_δ = ∫f\sin(14u)\sinh(δu)`. The independent program uses
  the complex-sinc route `A_δ − iB_δ = Σ_k x_k[\mathrm{sinc}(w_k−ζ) + \mathrm{sinc}(w_k+ζ)]L/\sqrt L`,
  `ζ = γ+iδ` — a different algebraic path.
* **interval arithmetic** — re-implemented, `2^{-200}` grid instead of `2^{-260}`, its own
  `exp`/`sin`/`cos`/`sqrt`/`log`/`atan` with their own proved remainders.

Shared and declared as such: the frozen input; the classical explicit formula; and the
canonical `(A_k, B_k)` reduction of the autocorrelation, recoded here and cross-checked
against `C(0) = Σx_k²` and `C(2L) = 0`. **This is not external peer review.**

### 1.3 The cross-check

| quantity | cold run | independent run | |
|---|---|---|---|
| `pole` | `0.45715833231743186` | `0.45715833231743186` | agree |
| `log π` term | `1.1447298858493788` | `1.1447298858493788` | agree |
| `prime` | `0.51733913650311314` | `0.51733913650311314` | agree |
| `arch` | `[1.2100349293664769, 1.2101336300807928]` | `[1.2100719057373517, 1.2101354876052193]` | **overlap** |
| `Q_ζ(f)` | `[0.0051242393314170015, 0.0052229400457328018]` | `[0.0051612157022915481, 0.0052247975701591636]` | **overlap**, intersection width `6.17e-05` |
| `A_{2/5}` | `−0.0091200379509347295` | `−0.0091200379509347295` | agree (two algebraic routes) |
| `B_{2/5}` | `+0.070827971459836173` | `+0.070827971459836173` | agree (two algebraic routes) |
| `Q_{add,2/5}(f)` | `[−0.014609466464138512, −0.014510765749822711]` | `[−0.014572490093263965, −0.014508908225396348]` | **overlap**, both strictly negative |

Neither `Q_ζ` enclosure contains the other — they *overlap*, which is what two different
error budgets should do. The independent one is narrower (`6.36e-05` vs `9.87e-05`) precisely
because the `ψ`-route carries no `1/M` tail. **Each program on its own yields a strictly
negative upper bound**: `−0.014510765749822711` and `−0.014508908225396348`.

### 1.4 Audits run as separate programs

* `audit_primitives.py` — **30/30 properties pass.** Containment for 400 random rationals
  (including negatives); `add/sub/mul/neg/sq/div` enclose the exact value; `mul` and `sq`
  correct in all four sign quadrants and on straddling intervals; division by an interval
  containing `0` **refused** in all four forms; `sqrt` brackets verified by squaring
  (`lo² ≤ x ≤ hi²`), never against libm; `exp(t)exp(−t) ∋ 1` and `exp(t/2)² = exp(t)` as
  interval statements; `sin²+cos² ∋ 1`, parity, and both double-angle identities;
  `exp(\log 2) ∋ 2`, `exp(\log 3) ∋ 3`, `exp(\log 2+\log 3) ∋ 6`, `exp(\log π) ∋ π`; the
  prime cutoff `e^{8/5} < 5` and `\log 4 < 8/5`.
  Found and documented: `_sincos_rat` has **no range reduction**; with `nmax = 1200` its
  proved truncation reaches `|t| ≈ 371`, and beyond that it **raises** rather than returning
  an unsound enclosure (checked at `t = 1120`). The largest argument the certificate uses is
  `γ_{max}L = 224`, inside the range.
* `audit_cache_invalidation.py` — **11/11 pass.** The unmutated cache is accepted; each of
  ten single-field mutations (any `x_k` by one unit in the last place, `x` denominator, `L`,
  `γ`, `M`, grid, `ndim`, source tag, a deleted field, a missing key block) is **rejected**
  and forces a recompute. A further control: a cache whose key matches but whose *numbers*
  are corrupted is not caught by the key — the key is an input fingerprint, not an output
  checksum — but **is** caught downstream by the gates, which raise
  `one or more gates failed; no certificate is claimed`.

### 1.5 The analytic audit

`ANALYTIC_AUDIT.md` now proves, rather than samples, the four things the arithmetic leans on:

1. **`C(2L) = 0` exactly**, from the support: the two translates overlap in the single point
   `{L}`, a null set, so the integral vanishes. The earlier `|C(2L)| < 9×10^{−78}` is a
   consistency check on the closed form and is **not** what justifies the boundary term.
2. **`\sup_{[0,U]}|C'| ≤ Σ_k(|A_k|w_k + |B_k| + |B_k|Uw_k)`**, by the triangle inequality on
   the closed form — valid at every point, no sampling. Value `≤ 11.436722214306984`.
3. **Explicit-formula applicability**: `f` is continuous with `f(±L) = 0` (because
   `w_kL = (2k+1)π/2` makes `\cos(w_kL) = 0` *exactly*), `F` is entire and even, and two
   integrations by parts give `h = F² = O((1+|r|)^{-4})` uniformly on `|\mathrm{Im}\,r| ≤ 1/2+ε`.
   All three sides converge absolutely. **No statement about zero locations is used.**
4. **The tail and the Euler-constant cancellation**, with indexing and signs: `m` starts at
   `0`, `M` is the first omitted index, `M+1 = 2^{17}` so `\ln(M+1) = 17\log 2` exactly, and
   `γ_E` is eliminated by `H_M − γ_E = \ln(M+1) − T` with `T ∈ [1/(2(M+1)) − 1/(6M²), 1/(2M)]`
   — so **`γ_E` is never enclosed and no harmonic number is summed**.

It also records that the phase reduction is exact only because every `(w_j ± w_k)L` is an
**integer** multiple of `π` at `L = 4/5` (verified for all 16 pairs), and that the prime sum
is **finite and complete** at `n ∈ {2,3,4}`, not truncated.

---

## 2. What exact theorem is now supported

Let `L = 4/5`, `U = 2L`, `γ = 14`, `φ_k(x) = \cos\big((2k+1)\tfrac{π x}{2L}\big)/\sqrt L` for
`k = 0,1,2,3` on `[−L,L]` and `0` outside, and

```
f = sum_{k=0}^{3} x_k phi_k ,
x = (99815221908, -354082045264, 928280739889, 54384691304) / 10^12 ,
||f||_2^2 = 999999999999981276942897/10^24 .
```

Let `Q_ζ` be the finite-window Weil functional `Pole + Arch − LogPi − Prime` evaluated on
`h = F²`, and for `δ > 0` let

```
Q_add,delta(f) = Q_zeta(f) + 4 ( A_delta^2 - B_delta^2 ) ,
A_delta = int f(u) cos(gamma u) cosh(delta u) du ,   B_delta = int f(u) sin(gamma u) sinh(delta u) du ,
```

the value of the same quadratic form after **adding** a quartet of zeros at
`1/2 ± δ ± iγ` to the zero multiset while leaving the true zeros in place.

> **Theorem (certified).** With the data above,
>
> 1. `Q_ζ(f) ∈ [0.0051242393314170015, 0.0052229400457328018]`, in particular `Q_ζ(f) > 0`;
> 2. `Q_{add,2/5}(f) ≤ −0.014510765749822711 < 0`;
> 3. `Q_{add,δ}(f) ≤ −0.0024549095780663143 < 0` for **every** `δ ∈ [1/4, 49/100]`.
>
> Both (2) and (3) are enclosures in exact rational interval arithmetic with outward
> rounding; (3) is a union of 256 abutting closed cells, each enclosure valid for every `δ`
> in its cell — not a grid of sample points. Claim (2) is additionally confirmed by a
> separately implemented program with an independent archimedean closed form, which obtains
> `Q_{add,2/5}(f) ≤ −0.014508908225396348 < 0`.

**What this means.** A four-mode window test certifies positivity for the actual zero set
and strict negativity once a quartet off the critical line at displacement
`δ ∈ [1/4, 1/2)` is *added*. The margin is `η_{point} = 1.45×10^{-2}` at `δ = 2/5` and
`η_{interval} = 2.45×10^{-3}` uniformly on `[1/4, 49/100]`, i.e. `294×` and `50×` the
`4.94×10^{-5}` archimedean tail allowance.

**What this does not mean.** See §4.

### 2.1 Two functionals, kept distinct

`Q_{rep,δ} = Q_ζ − 2F(14)² + 4(A_δ² − B_δ²)` is the **count-preserving** variant, which moves
a quartet off the line instead of adding one. The two differ by exactly `2F(14)²`:

```
F(14) = 4.129358070779164e-13 ,   Q_add - Q_rep = 2F(14)^2 = 3.410319615341804e-25 ,
```

far below every margin, so both are negative here — but they are different objects and the
certificate is stated for `Q_add`. `F(14)` is **enclosed and carried**, never set to zero.

---

## 3. Three claims of mine that were overextended, now corrected

**(A) "`R ≤ 2K_{eff}δ²` fails at every positive `δ`" — false.** It holds for small `δ`. The
certified ratio `R/(2K_{eff}δ²)` is

| `δ` | `1/10` | `3/20` | `166/1000` | `167/1000` | `1/5` | `2/5` |
|---|---|---|---|---|---|---|
| ratio | `0.330055` | `0.804476` | `1.010798` | `1.024652` | `1.549295` | `8.534307` |
| holds | yes | yes | no | no | no | no |

so the Taylor correction is controlled up to `δ ≈ 0.166` and fails above it. The real
obstruction is therefore **not** that the bound never holds; it is that the region where the
Taylor correction is controlled (`δ ≲ 0.166`) and the region where the certificate is
negative (`δ ≥ 1/4`) **do not overlap**. That gap, `(0.166, 0.25)`, is the honest statement.

**(B) "`N = 8` and `N = 16` do not detect" — overstated.** In `src/zeta_window.py`,
`freq(k, L, EVEN_D) = (2k+1)π/(2L)` depends on `k` and `L` only, **not on `N`** (verified in
source). So at a fixed `L = 4/5` the bases are **nested**: the `N = 4` basis is literally the
first four elements of the `N = 8` and `N = 16` bases, and zero-padding `x` to length `8` or
`16` gives the identical function `f`, hence the identical certificate. The larger bases
therefore *do* detect this witness. What is true is only that the *search* at `N = 8, 16` did
not find a witness — a statement about the search, not about the basis.

**(C) The `θ < 1/4` route and the `γ → ∞` behaviour.** The `θ < 1/4` proposal does not reach
small `δ`, so it does not close the `(0.166, 0.25)` gap. And `K_{eff} → 0` as `γ → ∞`: the
certified `‖P_E r_γ‖²` falls off fast,

| `γ` | `14` | `70` | `140` | `280` |
|---|---|---|---|---|
| `‖P_E r_γ‖²` | `3.601682e-02` | `1.100552e-05` | `2.128466e-06` | `1.129069e-07` |

so this four-mode window gives **no uniform-in-`γ`** statement, and nothing here extends to
large height.

**(D) Dependency control, recorded rather than hidden.** A single interval evaluation over
the whole box `δ ∈ [1/4, 49/100]` returns `Q_{add} ≤ +0.02790648537258906` — *not* negative.
`δ` occurs independently in `δ`, `δ²`, `\cosh(δL)` and `\sinh(δL)`, so one box inflates the
result by about an order of magnitude. The 256-fold subdivision is what makes (3) work, and
the failed single-box attempt is kept in `cold_certificate.json` as a documented control.

---

## 4. What remains outside scope

* **No statement about the actual zeros of `ζ`.** The certificate says a *hypothetical added*
  quartet is detected. It says nothing about whether such a quartet exists, and the
  positivity in (1) is a computation on a four-dimensional window, not evidence for RH.
* **No new RH criterion**, and no improvement to the classical Weil criterion.
* **No optimal threshold.** `δ ≥ 1/4` is what this witness reaches; the region
  `δ ∈ (0, 1/4)` is uncovered, and the gap `(0.166, 0.25)` separates the controlled-Taylor
  region from the certified-negative region.
* **No all-height theorem.** `γ = 14` only; `K_{eff} → 0` as `γ → ∞` (§3C).
* **No formal verification.** The arithmetic is exact-rational and interval-sound, and the
  analytic steps are proved in `ANALYTIC_AUDIT.md`, but CPython, its `int`/`Fraction`
  implementation and the operating system remain computational dependencies. There is no
  Lean or Coq artifact.
* **Not external peer review.** §1.2 is a second implementation by the same author in the
  same session. It removes the single-implementation risk; it does not remove the
  single-author risk.
* **The classical inputs are cited, not reproved**: Weil's explicit formula (Bombieri; Zhu
  for the geometric-side assembly) and the digamma series and recurrence. No novelty is
  claimed for them, for interval arithmetic, or for the elementary projection identity behind
  the witness.
