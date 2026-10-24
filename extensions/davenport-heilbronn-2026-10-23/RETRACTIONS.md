# Retraction and correction log

Entries are append-only. Nothing is deleted; a superseded statement stays visible with
the correction attached.

---

## R1 — `dh.N_count` used the two-sided Riemann–von Mangoldt normalisation

**Stated (first draft of `src/dh.py`):** `N(T) = (T/π) log(5T/(2πe))` as the main term
for the zero count, copied in form from the parent `dirichlet_window.N_count`.

**Wrong because:** that is the count over `|Im ρ| < T`. The census counts
`0 < Im ρ < T`, which is half as many. The first contour run showed winding numbers
`8, 21, 54, 130` against "main terms" `15.8, 42.7, 107.5, 259.1` — a factor of almost
exactly 2 — which is what exposed it.

**Corrected to:** `N(T) = (T/2π) log(5T/(2πe))`, the one-sided count. The deviations
then become `+0.08, −0.36, +0.26, +0.71, +0.45`, which is the expected `O(log T)`
fluctuation.

**Blast radius:** diagnostic only. `N_count` is never used as data — the counts come
from the argument principle. No reported number changed.

---

## R2 — the first off-line search contour ran too close to the critical line

**Stated (first census implementation):** off-line zeros located by recursive box
subdivision starting from `0.505 ≤ Re s ≤ 1.4`.

**Wrong because:** at `Re s = 0.505` the contour passes within `0.005` of the on-line
zeros, which is a small fraction of the mean spacing; `arg F` then swings by nearly `π`
between consecutive samples and the adaptive refinement recursed to its depth limit on
every band. The run made no progress for 25 minutes.

**Corrected to:** the located search runs at `Re s ≥ 0.55`, and the uncovered band
`1/2 < Re s < 0.55` is closed by the **reconciliation** instead: the contour count over
the whole strip must equal `#on-line + 2 · #off-line`, so an unlocated zero in the band
would show up as a reconciliation failure. This is a strictly stronger check than a
separate contour over the band, and it is cheap.

**Blast radius:** method only; no reported number changed.

---

## R3 — `dh.G_at` divided by zero at the removable singularity `z = ± i w`

**Stated:** the closed forms
`∫cos(wu)e^{zu}du = 2[z cos(wL) sinh(zL) + w sin(wL) cosh(zL)]/(z²+w²)` were evaluated
without a special branch.

**Wrong because:** at `z = ±iw` both numerator and denominator vanish and mpmath raises
`ZeroDivisionError`. The unit test in Phase B [B3] deliberately evaluates there and
caught it.

**Corrected to:** an explicit limit branch,
`∫cos(wu)e^{±iwu}du = L + sin(2wL)/(2w)` and
`∫sin(wu)e^{±iwu}du = ±i[L − sin(2wL)/(2w)]`, taken when `z² + w² = 0` exactly. The
test now passes at `z = iw` exactly and at `z = iw + 10⁻²⁵`.

**Blast radius:** none — the condition never arose in a reported computation, because
`H` is evaluated at zero ordinates, not at basis frequencies. The branch exists so that
the unit test can cover the case.

---

## R4 — the first independent archimedean route was too weak to be a check

**Stated:** Phase B's route 2 evaluated
`(1/2π)∫|G(it)|² Re ψ(3/4+it/2) dt` by direct quadrature on `[−600, 600]`.

**Wrong because:** the integrand oscillates with period `π/L`, i.e. about 380 times over
the interval at `L = 2`, and `mp.quad` on a coarse node list resolved it only to about
`10⁻⁶`. A `10⁻⁶` agreement between two routes cannot distinguish a correct assembly from
one with a small error, so it was not a usable check.

**Corrected to:** route 2 now uses Gauss's integral
`ψ(z) = ∫₀^∞ [e^{−x}/x − e^{−zx}/(1−e^{−x})] dx`, which after the `t`-integration gives

```
Arch = ∫₀^∞ [ h(0) e^{−x}/x − e^{−3x/4} h(x/2)/(1−e^{−x}) ] dx
     = ∫₀^{4L}[ … ] dx + h(0) E₁(4L),
```

a one-dimensional **non-oscillatory** integral with an exact tail. The two routes now
agree to `6.2 × 10⁻²¹`, limited by the deliberately omitted piece over `[0, x₀]`
(`x₀ = 10⁻²⁰`), whose analytic bound is `1.7 × 10⁻²⁰`.

**Blast radius:** the weak route was never used for a reported value; it was replaced
before Phase B's gate was declared.

---

## R5 — the original Phase B gate threshold was unreachable by construction

**Stated:** Phase B passed only if `|route 1 − route 2| < 10⁻⁴⁰`.

**Wrong because:** route 2 deliberately omits the integrand over `[0, x₀]` with
`x₀ = 10^{−dps/3} = 10⁻²⁰`, so its own declared budget is `≈ 10⁻²⁰` and the threshold
could never be met. A gate that cannot pass is not a gate.

**Corrected to:** the threshold is `10 ×` the route's own declared `[0, x₀]` bound, which
is the quantity the comparison is actually limited by.

**Blast radius:** none; the measured agreement (`6.2 × 10⁻²¹`) is reported unchanged.

---

## R6 — decimal window lengths were parsed at the ambient precision, not at a fixed one

**Stated:** window lengths were built with `mp.mpf(Ls)` wherever they were first needed.

**Wrong because:** `mp.mpf` rounds a decimal string to the *current* `mp.mp.dps` and the
resulting value is frozen. `src/freeze.py` built `L = mp.mpf('1.8')` before raising the
precision, so `L` carried mpmath's **default 15 digits**; the subsequent precision ladder
at 60 / 100 / 160 digits then evaluated the form at that same slightly-wrong `L` three
times, agreed with itself to `10⁻⁸⁵`, and reported a value that differed from the correct
one in the **13th significant digit**:

```
freeze.py (L at 15 digits)   Q/||g||^2 = -1.68123895523673581163409971438e-16
correct   (L at 220 digits)  Q/||g||^2 = -1.68123895523670103326233298941e-16
```

**How it was found:** `src/phaseD.py` and `src/freeze.py` froze the *same* witness and
disagreed in the 13th digit. A ladder that agrees with itself to 85 digits while two
runs of the same computation disagree in the 13th is a parsing bug, not a precision
limit, and that is what the discrepancy turned out to be.

**Corrected to:** `dh.mpf_hi(x)` parses every decimal constant at a fixed 220 digits
regardless of the ambient precision, and every window length, located zero and
perturbation parameter in `src/` now goes through it. `tests/test_precision.py` fails if
any of them is reconstructible at the ambient precision instead.

**Blast radius.** Only non-binary-exact decimals are affected:
`1.2, 1.4, 1.8, 2.2, 2.4, 2.6` are, `1.0, 1.6, 2.0, 0.5` are **not**
(`1.6 = 2^4/10` is not exact either — see below).

* The **certifiable witness is unaffected**: it sits at `L = 2.0`, which is exactly
  representable in binary at every precision, and the Arb certificate parses `2.0`
  independently anyway.
* The **window scan is unaffected**: inside `phaseD`'s loop the ambient precision was
  already `60` from the second window onward, and the first window is `L = 1.0`, exact.
  Every scanned value was computed and reported at a consistent 60 digits.
* The **smallest-`L` witness at `L = 1.8` was affected** and has been recomputed. Its
  sign and its first twelve digits are unchanged.
* `phaseA`, `phaseB`, `census` and `phaseC` set the precision before constructing any
  constant and were already correct; they were converted anyway so that the invariant is
  uniform and testable.

**Correction to this entry's own first draft:** `1.6` was initially listed as
binary-exact. It is not (`1.6 = 8/5`). It is nonetheless unaffected, for the separate
reason given above — the scan ran at a consistent precision — not because the decimal is
exact.
