# Davenport–Heilbronn control for the finite-window Weil detector

Lab record. Extension of `weil-window-control`; **no DOI, not part of the paper, not
peer-reviewed.** Derivations: [`DERIVATION.md`](DERIVATION.md). Sources:
[`SOURCE_AUDIT.md`](SOURCE_AUDIT.md). Errors made and corrected:
[`RETRACTIONS.md`](RETRACTIONS.md). Reproduction: [`REPRODUCE.md`](REPRODUCE.md).

---

## 0. Non-claims

Written before any number was computed.

1. **This experiment does not prove or disprove RH or GRH.** It computes a quadratic
   form attached to a function that is not `ζ` and has no Euler product.
2. The Davenport–Heilbronn function has known off-line zeros and lacks the Euler
   product of a Dirichlet `L`-function. **Nothing here asserts that the absence of an
   Euler product alone causes off-line zeros.** The two facts are recorded side by side
   because they are both true of this function, not because one is claimed to imply the
   other.
3. **Detection calibrates the method; it does not transfer to `ζ`.** A detector that
   sees this function's off-line zeros is a detector that is not blind; it is not
   thereby a detector that would see `ζ`'s, if `ζ` had any.
4. A positive finite matrix proves positivity **only on its stated test space**, never
   on the whole finite-window function space.
5. Non-detection certifies nothing about positivity outside the searched spaces. Where
   a window shows no negative witness, the statement is "no negative witness was found
   in the declared basis", never "the form is positive there".
6. **Numerical root residuals do not establish census completeness.** A small `|F(ρ)|`
   says a root was found, not that no root was missed.
7. Agreement with published decimal zeros is a cross-check, **not an interval
   certificate**.
8. Synthetic zero modifications are **not claimed to define L-functions**. Nothing in
   Phase F or Phase G attaches an analytic function, an Euler product or a functional
   equation to a modified spectrum.
9. No claim is made about **which parity sector detects first globally**; the searches
   here are finite.
10. "Smallest tested window with a detected negative witness" is never upgraded to
    "true first detection threshold".

### Evidence tiers

| tier | meaning |
|---|---|
| **T1** | an explicitly justified analytic or algebraic statement |
| **T2-N** | high-precision numerical evidence |
| **T2-C** | a computational conclusion with rigorous interval/ball bounds and a reproducible certificate |
| **T3** | interpretation, or an unproved hypothesis |

High precision alone is **T2-N**. T2-C here means an Arb ball enclosure; it is **not**
Lean verification and is nowhere described as such.

---

## 1. Repository safety audit, recorded before any edit

```
repository   /home/elias/GitHub/weil-window-control
branch       main
HEAD         35c7d95  "Fill in the citation: merge commit hash, and tag
                       dirichlet-windows-v1 on this commit"
remote       https://github.com/innerlightr-wq/weil-window-control.git
git status --short
             ?? extensions/newtonian-weil-response-2026-10-08/
             ?? paper/weil_window_control_rev_certified_witness.fdb_latexmk
             ?? paper/weil_window_control_rev_certified_witness.fls
```

Those three untracked paths pre-date this work and **must remain untouched**; they are
re-checked at the end (§13). No branch was switched, nothing was reset, committed,
pushed or published, and no licence was changed. No existing source file or regression
fixture was modified: the two pieces of upstream code this extension needs
(`src/zeta_window.py`, `extensions/dirichlet-windows-2026-10-23/src/dirichlet_window.py`)
are **imported**, never edited, and where different behaviour was needed it is supplied
by an extension-local wrapper.

```
Python        3.14.7 (CPython)      mpmath   1.3.0
python-flint  0.9.0 (Arb)           numpy / scipy / matplotlib: not installed, not used
```

`python-flint` was installed into the repository's `.venv` for this work; it is additive
and does not alter any pinned behaviour (the regressions in §13 reproduce the committed
upstream values under it).

**Read before editing:** `extensions/dirichlet-windows-2026-10-23/REPORT.md` and
`README.md` and its `src/`; `notes/proofs.md`; `notes/stage3_assembly.md`;
`src/zeta_window.py`; and the statements of Theorems 7 and 8 and Corollary 9 in
`paper/weil_window_control.tex`. What they actually say is recorded in
`SOURCE_AUDIT.md` §E1 — in particular Theorem 7 is a **relative** estimate (its
`λ*(L) → 0` substitution was withdrawn in Retraction 18) and Theorem 8's hypothesis is
stated there as **sufficient, not necessary**.

---

## 2. Source audit

Full record in [`SOURCE_AUDIT.md`](SOURCE_AUDIT.md). The three things that matter for
reading what follows:

1. **Nothing here depends on reading Davenport–Heilbronn, Spira, Balanzario–Sánchez-Ortiz
   or Bombieri–Ghosh.** The function, `κ`, the functional equation, the explicit formula
   and the zero census are all derived and computed from scratch. The AMS PDF of
   Balanzario–Sánchez-Ortiz returned **HTTP 403** from this machine and was not read.
2. **The brief's numeric `0.8085171824566374 + 85.6993484853776 i` was used only as a
   Newton seed.** Refinement moved it by `7.8 × 10⁻¹⁵`.
3. **The reported `L = 1.8` witness in the literature is UNVERIFIED here.** The document
   was located — Tariq Ziad, *Localized Constrained Weil Forms: Near-Degeneracy and a
   Davenport–Heilbronn First Kernel*, [doi:10.5281/zenodo.23145455](https://doi.org/10.5281/zenodo.23145455),
   7 October 2026 — but only its landing page was read. It claims "an exactly constrained
   negative witness at `a = 9/5`" for **odd, constrained** test functions. A separate
   Zenodo series reports "positive up to a window `L ≈ 3.6` … `λmin = −0.73` at `L = 4`,
   with a minimizer concentrated on its off-line zero `0.8085 + 85.70i`". Those two
   window numbers differ by exactly a factor 2, which is *consistent with* a `g` versus
   `h = g ⋆ g` support convention but does not establish one. **No quantity in this
   extension was fitted to either value**; the scan range was fixed by the brief before
   the documents were located. §12 Q8 says what is and is not novel here.

---

## 3. Phase A — the function (T1 for the derivation, T2-N for the checks)

`DERIVATION.md` §1. Results in `results/phaseA.txt`.

| check | result |
|---|---|
| `\|τ(χ)\| − √5` | `0` to 60 digits |
| `\|ε_χ\| − 1` | `0` to 60 digits |
| `κ = tan(arg ε_χ /2)` **derived** | `0.2840790438404122960282918323931261690911…` |
| `κ` closed form `[√(10−2√5)−2]/(√5−1)` | agrees, difference `2.7 × 10⁻⁵¹` |
| root-number identity `(1−iκ)ε_χ − (1+iκ)` | **exactly `0`** |
| `(a_1..a_5)` | `(1, κ, −κ, −1, 0)` |

`κ` is **derived from the Gauss sum**, not imported. The closed form is a *check*, and
the exact vanishing of `(1−iκ)ε_χ − (1+iκ)` is the statement that the completed equation
`Λ_F(s) = Λ_F(1−s)` holds — distinct from the individual character's equation
`Λ(s,χ) = ε_χ Λ(1−s, χ̄)`, which exchanges `χ` and `χ̄`.

**Three evaluation routes.** A (the `L`-function combination) and B (the
periodic-coefficient Hurwitz form) agree to `≤ 1.2 × 10⁻⁴⁸` relative; they **share
`mp.zeta(s,a)`** and the sharing is not disguised. Route C recomputes the Hurwitz zeta by
**Euler–Maclaurin written here**, calling only `power`, `bernoulli` and `factorial`, and
agrees to `≤ 1.7 × 10⁻⁴⁶`. Route D sums `Σ a_n n^{−s}` directly at `Re s = 4` and agrees
to `4.1 × 10⁻¹⁸` against a tail bound of `4.2 × 10⁻¹⁴`.

**Functional equation**, 20 reproducible points (fixed decimal strings, avoiding the
critical line, the real axis and the gamma poles), with a **scale-aware** residual so
that a tiny `|Λ_F|` cannot manufacture agreement:

```
worst relative residual at  60 digits : 6.95e-58
worst relative residual at 100 digits : 8.09e-99       <- scales with precision
```

Both absolute and relative residuals are tabulated in `results/phaseA.txt`; the point at
`0.81 + 85.7i`, where `|Λ_F| = 2.1 × 10⁻³¹`, has absolute residual `1.4 × 10⁻⁸⁸` and is
the reason the relative column exists.

**Entireness and the strip (T1).** `F(−1) = −5.9×10⁻⁶²`, `F(−3) = 2.5×10⁻⁶¹`, … — `F`
vanishes at every pole of `Γ((s+1)/2)`, so `Λ_F` is entire and there is **no pole term**.
For `Re s ≥ 2`, `|F(s)−1| ≤ ζ(2)−1 = 0.6449 < 1`; sharpening the same elementary estimate
with the actual `|a_n|` gives the exact abscissa

```
σ_0 = 1.395136158235109721061359…,   Σ_{n≥2}|a_n| n^{−σ_0} = 1,
```

so **all zeros of `Λ_F` lie in `1 − σ_0 < Re s < σ_0`**, i.e. `−0.3952 < Re s < 1.3952`.
This is the statement that replaces an arbitrary right cutoff, and it is what makes the
census contour meaningful. Note how *wide* it is compared with `ζ`: this function is
known to have zeros with `Re s > 1`.

**The selected zero** (T2-N; search seed only):

```
β_1 = 0.8085171824566373855533519606068441278507…
γ_1 = 85.69934848537759217192926770894172903799…
δ_1 = 0.3085171824566373855533519606068441278507…
|F(ρ_1)| = 4.0e-101    |F/F'| = 3.2e-101    |seed − ρ_1| = 7.8e-15
|ρ_1(100 dps) − ρ_1(140 dps)| = 3.2e-106
```

The whole orbit `{β_1 ± iγ_1, 1−β_1 ± iγ_1}` was verified, and in the census the point
carries **local winding number exactly 1** on `|s − ρ_1| = 0.02`: a simple zero of the
analytic function `F`. This bears directly on a contested point in the literature
(`SOURCE_AUDIT.md` B3) and is recorded there.

**PHASE A GATE: PASS.**

---

## 4. Phase B — the arithmetic side (T1 derivation; T2-N checks)

`DERIVATION.md` §§2–3. Results in `results/phaseB.txt`, coefficients in
`data/dh_coefficients.csv`.

The form derived — **with every factor accounted for**, and matching the brief's target:

```
Q_F(g) = h(0) log(5/π) + (1/2π)∫_ℝ |G(it)|² Re ψ(¾ + it/2) dt
         − 2 Σ_{2 ≤ n ≤ e^{2L}} b_n h(log n)/√n .
```

The `2` is the doubling in `Λ_F(s)=Λ_F(1−s)`; `log(5/π)` is the conductor and the gamma
factor in one prefactor; `ψ(¾ + ·)` is the odd parity `a = 1`. **No pole term.** The
`L`-function Weil matrices are **not** linearly combined — logarithmic differentiation
does not commute with linear combination — only the archimedean and conductor machinery,
which depends on the gamma factor and conductor alone, is reused at `a = 1`.

| check | result |
|---|---|
| `b_n` recursion vs **Dirichlet-inverse** route, `n ≤ 200` | `max diff 1.2 × 10⁻⁶⁰` (max `\|b_n\| = 6.78`) |
| `b_2 − κ log 2` | **exactly 0** |
| `b_3 + κ log 3` | **exactly 0** |
| `b_6 − (1+κ²) log 6` | **exactly 0** |
| `b_n` for `5 \| n` | **exactly 0** |
| `b_1` | **0** |

**The non-prime-power control.** For `ζ` the set of `n ≤ 50` that are not prime powers and
have `b_n ≠ 0` is *empty*. Here it is

```
6, 12, 14, 18, 21, 24, 26, 28, 34, 36, 38, 39, 42, 46, 48
```

which is the missing Euler product, visible in the coefficients the detector actually
consumes.

**The cutoff** is decided in integer arithmetic (`n ≤ ⌊e^{2L}⌋ + 2`) and then re-tested
as `log n < 2L` at full precision, never by a binary float. At `L = 2.6`, `n` runs to
**181** with margin `2L − log 181 = 1.50 × 10⁻³`; 122 terms are nonzero.

**Transforms.** `G_at` versus quadrature over 30 `(sector, mode, z)` combinations,
including `z = iw` exactly (the removable `0/0`, which is where the first version crashed
— `RETRACTIONS.md` R3): worst relative discrepancy `5.5 × 10⁻³⁶`. Decay:
`sup_{50≤t≤1050} |G(it)| t² = 16.0` (even), `17.7` (odd) — bounded, so `H = O(t⁻⁴)` and
the contour's horizontal edges vanish.

**Independent scalar route.** The archimedean term is recomputed from **Gauss's integral**
`ψ(z) = ∫₀^∞[e^{−x}/x − e^{−zx}/(1−e^{−x})]dx`, which after the `t`-integration becomes

```
Arch = ∫₀^{4L}[ h(0)e^{−x}/x − e^{−3x/4} h(x/2)/(1−e^{−x}) ] dx + h(0) E₁(4L),
```

a **one-dimensional, non-oscillatory** integral with an exact tail — a genuinely different
representation, not the same series rearranged. Over six `(L, sector)` cases the two
routes agree to `6.2 × 10⁻²¹`, against the analytically bounded `[0, x₀]` piece of
`1.7 × 10⁻²⁰` that route 2 deliberately omits. The quadrature on `[x₀, 4L]` is a
**library estimate and is not claimed as a bound**; the only rigorous pieces are the
`[0,x₀]` bound and the exact `E₁` tail. (The first version of this route used frequency-side
quadrature and was too weak to be a check at all — `RETRACTIONS.md` R4.)

**PHASE B GATE: PASS.**

---

## 5. Orbit conventions (T1 derivation; T2-N unit tests)

`DERIVATION.md` §4, tests in `tests/test_orbits.py`, all passing.

For a non-real off-line zero the **complete** orbit is four points and

```
Q_O(g) = 4 Re( A · conj B ) = 4 ( |P|² − |S|² ),
A = G(δ+iγ), B = G(−δ+iγ),
P = ∫ g cosh(δu) e^{iγu} du,   S = ∫ g sinh(δu) e^{iγu} du.
```

A **signed rank-two** form: `+4|P|²` on the `cosh` channel, `−4|S|²` on the `sinh`
channel. It is emphatically **not** `|G(δ+iγ)|²`, and the test suite carries that as an
explicit control (the two differ in both sectors).

| unit test | both sectors |
|---|---|
| on-line conjugate pair vs brute force | agrees to `2.3e-41`, and is `≥ 0` |
| full off-line quartet vs brute force over 4 shifts | agrees to `5.0e-41`; imaginary part vanishes |
| `Q_O = 4(\|P\|²−\|S\|²)` | agrees to `4.8e-41` |
| **collapse** `δ → 0` → **two DOUBLE** on-line zeros (not a simple pair) | residual `∝ δ²`: `2.1e-8, 2.1e-12, 2.1e-16` at `δ = 10^{−4,−6,−8}` |
| real off-line pair at `γ = 0` (**two** zeros, not four) | agrees to `1.4e-42` |
| its sign | `+` in the even sector, `−` in the odd sector |
| mixed even–odd terms vanish | `5.2e-41` |

The last two lines **recover the `dirichlet-windows` real-pair parity split
independently**, as the `γ → 0` limit of the quartet form rather than as a separate
calculation: at `γ = 0` the `sinh` channel survives only in the odd sector. Multiplicity
is preserved in the collapse — four zeros in, four zeros out.

---

## 6. Finite-height zero census (T2-N)

`src/census.py`, `data/dh_zeros.csv`, `results/census.txt`. Strip
`−0.4 < Re s < 1.4`, whose vertical edges are **provably** zero-free by §3.

**Counting by the argument principle**, as the adaptively refined total change of
`arg F` (the step halves whenever an increment exceeds `0.45π`, so no half-turn can be
skipped):

| `T` | winding number | integer | `(T/2π)log(5T/2πe)` | deviation |
|---|---|---|---|---|
| 25 | `8.0` | 8 | 7.9197 | `+0.080` |
| 50 | `21.0` | 21 | 21.3553 | `−0.355` |
| 100 | `54.0` | 54 | 53.7423 | `+0.258` |
| 150 | `91.0` | 91 | 90.2933 | `+0.707` |
| 200 | `130.0` | 130 | 129.5482 | `+0.452` |

Worst deviation of a winding number from an integer: **`1.3 × 10⁻²³`**. The deviations
from the main term are the expected `O(log T)` fluctuation. (The main term was wrong by a
factor 2 in the first draft — `RETRACTIONS.md` R1 — and the factor-2 pattern in this
column is exactly what exposed it.)

**Four off-line quartets below `T = 200`**, every one with local winding number exactly
`1` on a radius-`0.02` circle, i.e. simple:

| | `β` | `γ` | `δ = β − ½` |
|---|---|---|---|
| `O_1` | `0.808517182456637386` | `85.6993484853775922` | `0.3085171824566374` |
| `O_2` | `0.650830080609737082` | `114.163342730756981` | `0.1508300806097371` |
| `O_3` | `0.574356050450805991` | `166.479305913168156` | `0.0743560504508060` |
| `O_4` | `0.724257694626809780` | `176.702461242855825` | `0.2242576946268098` |

**Reconciliation.** 122 critical-line ordinates from sign changes of the Hardy-type `Z`
(`max |Z(γ)| = 6.3 × 10⁻⁵⁹` after refinement), plus `2 × 4` off-line zeros and their
functional-equation mirrors:

```
122 + 2·4 = 130 = the contour count.   RECONCILED.
```

The reconciliation is also what closes the band `½ < Re s < 0.55` that the located search
deliberately avoids (`RETRACTIONS.md` R2): a zero there would break the identity.

**Status, stated as the brief requires.** A *list of found roots*: yes. A *numerically
reconciled census*: yes. A *rigorously complete census*: **no** — the winding numbers are
integers to `1.3 × 10⁻²³` but are computed by adaptive sampling, not interval arithmetic.
Roots on `Re s = ½` are "consistent with the line at the tested precision"; **no
existence/uniqueness certificate is claimed for any of them**.

---

## 7. Phase C — the gate

`results/phaseC.txt`, `data/gate_checks.csv`.

**Route I — logarithmic-derivative contour.** With `σ_c = 1.42`, to the right of every
zero,

```
Q_F(g) = (2/π) Re ∫₀^∞ (Λ_F'/Λ_F)(σ_c+it) · H(σ_c−½+it) dt,
```

all edges accounted for: the left edge at `Re s = 1−σ_c` is the image of the right edge
under `s → 1−s` and is not a separate term, and the horizontal edges vanish by the decay
below. `Λ_F'/Λ_F` is built from `F` and `F'` **directly** (the Hurwitz representation
differentiated term by term) and **never from the `b_n`**, so this route is independent of
the entire arithmetic side.

**Constrained test functions.** Both bases already satisfy `f(±L) = f''(±L) = 0`; the one
extra linear condition `Σ c_k w_k (−1)^k = 0` gives `g'(±L) = 0` as well, hence
`|H(z)| ≤ e^{2L|Re z|}‖g'''‖₁²/|Im z|⁶` — which is what makes both tails small enough to
see ten digits.

**Route II — zero-side sum** over the census, with the tail above the census height
bounded by that majorant times a **measured** zero density (so the tail is T2-N, not a
proved bound, and is labelled as such).

| `L` | sector | `Q` arithmetic | `\|arith − contour\|` | budget | `\|arith − zero-side\|` | budget |
|---|---|---|---|---|---|---|
| 1.0 | even | `0.71945738530768` | **`1.73e-13`** | `1.3e-6` | `1.85e-11` | `7.1e-5` |
| 1.0 | odd | `0.29008262739593` | — | — | `5.65e-11` | |
| 1.5 | even | `0.26941548159636` | — | — | `1.35e-12` | |
| 1.5 | odd | `0.41008734069111` | — | — | `2.54e-12` | |
| 2.0 | even | `0.20367339222052` | — | — | `1.58e-13` | |
| 2.0 | odd | `0.37942651632836` | — | — | `3.97e-13` | |

Route I is the brief's *optional* comparison and each evaluation of its integrand costs
ten Hurwitz zetas, so it was run in full on one representative case rather than all six;
route II runs on all six.

**Gate control A** — negate the selected quartet on the zero side only:

```
Q_{O_1}(g)             = 8.720121523693614e-11   (resolved on this test function)
discrepancy, quartet present = 3.9707e-13
discrepancy, quartet NEGATED = 1.748e-10
predicted change  2 Q_{O_1}  = 1.74402e-10
observed  change             = 1.74402e-10
|predicted − observed|       = 1.83e-62
```

(Negating a contribution changes the zero side by `−2 Q_{O_1}`, hence the discrepancy by
`+2 Q_{O_1}`. The quartet is resolved on this low-frequency test function only at the
`10⁻¹¹` level — it is not a resonant one — which is why the discrepancy *with* the
quartet, `4 × 10⁻¹³`, is smaller than `Q_{O_1}` itself. That is the point of the control:
the gate's residual moves by exactly the predicted amount, to `10⁻⁶²`.)

**Gate control B** — drop `b_6` at a window containing `n = 6`:

```
b_6 = (1+κ²)log 6 = 1.936356076621039      h(log 6) = −0.2398791703848834
predicted change = +2 b_6 h(log 6)/√6 = −0.3792557128259818
observed change                       = −0.3792557128259818
|predicted − observed| = 1.95e-61
```

### GATE CLASSIFICATION: **PASS-N**

Converged numerical agreement, explicitly uncertified. It is **not PASS-C** because the
quadrature in route I and the zero-density in route II's tail are a library estimate and
a measured quantity respectively, not rigorous bounds. Controls A and B are exact algebra
and are the parts that would catch a wrong factor, a wrong conductor, a wrong parity or a
dropped term — they agree to `10⁻⁶¹`.

---

## 8. Phase D — the window scan and the witness

`data/window_scan.csv`, `results/phaseD.txt`, `results/freeze.txt`.

Declared basis schedule, fixed in `src/config.py` before the scan: level 1 = low block
`{0..7}` plus resonant block `k₀ ± 4`; level 2 = `{0..15}` plus `k₀ ± 10`; level 3 =
`{0..23}` plus `k₀ ± 16`, where `w_{k₀}` is the basis frequency closest to `γ_1`.
Enrichment can only lower `λ_min` (a larger subspace), and **levels 2 and 3 are reported
as enrichments, never folded back into the level-1 scan.**

`λ_min` of `Q_F` restricted to the stated space:

| `L` | even L1 | even L2 | even L3 | odd L1 | odd L2 | odd L3 |
|---|---|---|---|---|---|---|
| 1.0 | `4.486e-5` | `4.035e-5` | `3.883e-5` | `1.773e-2` | `1.356e-2` | `1.331e-2` |
| 1.2 | `1.568e-8` | `8.780e-9` | `8.232e-9` | `1.428e-5` | `9.703e-6` | `9.214e-6` |
| 1.4 | `1.637e-11` | `5.087e-14` | `1.445e-14` | `1.124e-8` | `6.539e-11` | `3.181e-11` |
| 1.6 | `1.973e-13` | `2.550e-19` | `2.661e-22` | `2.919e-11` | `2.154e-16` | `1.096e-18` |
| 1.8 | `2.914e-15` | `7.581e-24` | **`−2.662e-21`** | `1.191e-13` | `1.008e-20` | **`−1.681e-16`** |
| 2.0 | `−0.18873` | `−0.48965` | `−0.54739` | `−0.65865` | `−0.67907` | `−0.71018` |
| 2.2 | `−1.1752` | `−1.2199` | `−1.2387` | `−1.0385` | `−1.1784` | `−1.2149` |
| 2.4 | `−1.4106` | `−1.5030` | `−1.5531` | `−1.5521` | `−1.5903` | `−1.6059` |
| 2.6 | `−2.0290` | `−2.1377` | `−2.1532` | `−1.7127` | `−2.0055` | `−2.0523` |

Across all 54 rows, `|Rayleigh quotient of the returned eigenvector − eigenvalue|` is at
most `5.8 × 10⁻⁶⁰`: the eigensolver is not doing anything the direct evaluation does not
confirm. The Gram matrix of the normalised basis is the identity **exactly** (`max|G−I|`
between `4.8` and `7.4 × 10⁻⁶⁰`, condition number `1.0`), so `Bv = λMv` *is* the ordinary
eigenproblem, nothing was regularised and no direction was discarded.

**The shape of the table is the result.** At basis level 3, `λ_min` falls from
`3.9 × 10⁻⁵` at `L = 1.0` to `2.7 × 10⁻²²` at `L = 1.6` — seventeen orders of magnitude,
all of it positive — then turns negative at `L = 1.8` by only `10⁻²¹` (even) and `10⁻¹⁶`
(odd), and at `L = 2.0` by `0.5` to `0.7`. The crossing is not gradual: the form sits on a
knife edge for most of the scanned range and then falls off it.

**Smallest tested window with a detected negative witness: `L = 1.8`**, in **both**
sectors, at basis level 3. This is *not* a claim about a true first detection threshold,
and no claim is made about which sector detects first globally.

Two witnesses are frozen, for two different purposes (`src/freeze.py`):

* **`certificates/witness_smallest_L.json`** — `L = 1.8`, odd, level 3, `N = 57`:
  ```
  Q_F(g)/‖g‖² = −1.68123895523670103326233298940716630655e-16
  dps 60 / 100 / 160 agree to  1.7e-60  and  3.7e-87
  ```
  The sign is reproducible; the value is `10⁻¹⁶` against a form whose largest eigenvalue
  is of order `10`, far too small to build a ball certificate around at an affordable
  precision. **Reported, not certified.** (Two of my own runs first disagreed in the 13th
  digit here; that was a real bug in how decimal window lengths were parsed —
  `RETRACTIONS.md` R6 — and is fixed.)

* **`certificates/witness.json`** — the **certifiable** witness: the smallest window whose
  negativity is macroscopic, `L = 2.0`, odd sector, level 2, `N = 37`, with every
  coefficient stored as an exact 30-digit decimal string.

---

## 9. The certified negative witness (T2-C)

`certificates/ball_certificate.txt`, produced by `src/certify.py`, which **re-implements
the entire form against Arb** — the cross-correlations, `κ`, the `b_n`, and the
archimedean closed form are all written again, sharing no code with `src/dh.py`. Arb
supplies rigorous enclosures for exactly the three special functions needed: `digamma`,
the Hurwitz `ζ(s,a)` and the Lerch `Φ(w,s,z)`.

```
L = 2.0      sector = odd      N = 37      arithmetic cutoff n ≤ 56 (integer-decided)

Q_F(g_*)/‖g_*‖² = −0.679071346441480590986377586093131755967123639793791536…
                   with a rigorous ball radius ≈ 1.4 × 10⁻¹¹⁶

STRICTLY NEGATIVE : True
```

It agrees with the independent mpmath assembly to all 44 digits mpmath carries, which is
simultaneously a cross-check of both implementations.

**Complete error budget.** Everything is inside that single radius:

| contribution | status |
|---|---|
| rounding of every operation | Arb ball arithmetic |
| `κ` | square roots only: `[√(10−2√5)−2]/(√5−1)` |
| `b_n`, `n ≤ 56` | exact recursion over `κ` and `log n` |
| `log(5/π)`, `log n`, `√n` | Arb elementary functions |
| archimedean closed form | Arb `digamma`, Hurwitz `ζ`, Lerch `Φ` |
| cross-correlations `C_jk(u)` | Arb `sin`/`cos` of exact arguments |
| **truncation of the arithmetic sum** | **none** — `h` vanishes beyond `2L`, the sum is finite and complete |
| **spectral tail** | **none** — no zero is used |
| **quadrature** | **none** — no quadrature is used |
| **root uncertainty** | **none for `Q_F`**; enters `Q_{O_1}` only |

The coefficients are exact finite decimals, so the test function carries no error of its
own: the ball encloses `Q_F` for *this* `g`, not for an approximation to some other `g`.

**Tier: T2-C.** An Arb ball enclosure with a reproducible certificate. **It is not a Lean
verification and is nowhere described as one.**

---

## 10. Phase E — the remaining-zero form

`data/remaining_form.csv`, `results/phaseE.txt`. `Q_F` is evaluated through the
**arithmetic** formula throughout, so every remaining zero — including every off-line zero
above the census height and every zero never located — is included implicitly. Nothing is
replaced by a finite zero list. `Q_{O_1}` is the **exact** rank-2 orbit form of
`src/orbitform.py` (polarisation check against the scalar routine: `4.9 × 10⁻⁶³`).

### The decomposition at the frozen witness (T2-C)

```
‖g_*‖²        =  0.999999999999999999999999999999        (unit by construction)
Q_F(g_*)      = −0.679071346441480590986377586093        ±1.4e-116   certified
Q_{O_1}(g_*)  = −0.985198326439305481793695406946        ±1.2e-117   certified*
Q_rest(g_*)   = +0.306126979997824890807317820853        ±4.4e-116   certified*
additive reconciliation  Q_F − (Q_{O_1} + Q_rest)  = 0   (±3.0e-116)
```

\* certified *as a function of* the supplied `(β_1, γ_1)`; that those enclose an actual
zero rests on the census, which is T2. **The sign of `Q_F(g_*)` does not depend on it** —
`Q_F` is computed from the arithmetic side and knows nothing about where any zero is.

**The selected quartet carries more than the whole of the negativity**, and the rest of
the spectrum partially cancels it: `−0.985` from `O_1`, `+0.306` back from everything
else, net `−0.679`.

### Classification, per window and sector

`λ_min` of the finite restriction of `Q_rest` on the level-2 space:

| `L` | even: `λ_min(Q_F)` | even: `λ_min(Q_rest)` | odd: `λ_min(Q_F)` | odd: `λ_min(Q_rest)` | class |
|---|---|---|---|---|---|
| 1.0 | `4.03e-5` | `4.04e-5` | `1.36e-2` | `1.32e-2` | **R4** |
| 1.2 | `8.78e-9` | `8.47e-9` | `9.70e-6` | `9.55e-6` | **R4** |
| 1.4 | `5.09e-14` | `4.64e-14` | `6.54e-11` | `6.63e-11` | **R4** |
| 1.6 | `2.55e-19` | `2.25e-19` | `2.15e-16` | `2.21e-16` | **R4** |
| 1.8 | `7.58e-24` | `7.02e-24` | `1.01e-20` | `1.09e-20` | **R4** |
| 2.0 | `−0.4897` | `2.03e-27` | `−0.6791` | `1.21e-24` | **R4** |
| 2.2 | `−1.2199` | `1.22e-30` | `−1.1784` | `4.88e-28` | **R4** |
| 2.4 | `−1.5030` | `1.81e-33` | `−1.5903` | `4.26e-31` | **R4** |
| 2.6 | `−2.1377` | `8.47e-36` | `−2.0055` | `8.14e-34` | **R4** |

**Phase E classification: R4 at every tested window and sector — a numerically positive
finite restriction.** Not R3: no interval certificate was produced for the *restriction*
(only for the pointwise value at `g_*`, which is a different object). **Not R6, anywhere.**

The second column is the striking one. Once `Q_F` goes macroscopically negative,
`λ_min(Q_rest)` does not merely stay positive — it collapses to `10⁻²⁴`, `10⁻²⁸`,
`10⁻³¹`, `10⁻³⁶`, i.e. `Q_rest` sits **exactly on the boundary of positivity** with a
near-null direction, while `Q_F` is a long way below it. Every bit of the negativity is
`O_1`.

### Audit against the repository's Theorem 8

Theorem 8 (`\label{thm:zetasharp}`) defines `Q_rest(f) = Q_0(f) − 4F(γ)²` where `Q_0` is
the *augmented counterfactual* baseline of `notes/proofs.md` §B.1, convention 3 — a
**double** on-line pair at `γ` that is then moved off the line. So its baseline contains
the collapsed quartet and its "remaining" form is what is left after removing **four**
zeros.

* **Multiplicity convention: the same.** Four zeros removed, four zeros in our orbit.
* **Baseline: different.** Theorem 8 removes the `δ = 0` limit; we remove the four
  **actual** zeros at `δ_1 = 0.3085`. The difference is exactly
  `Q_{O_1}(δ_1) − Q_{O_1}(0)`, which Phase G measures on the frozen witness:
  `Q_{O_1}(0) = 4|F(γ_1)|² = 4 × 0.014405089² = 0.000830026` (Phase G's `t = 0` row is
  `Q_rest + Q_{O_1}(0) = 0.306127 + 0.000830 = 0.306957`), so the difference is
  `−0.985198326 − 0.000830026 = −0.986028352`.
* **Its hypothesis is FALSE here and is not assumed.** Theorem 8 supposes every other
  nontrivial zero is on the critical line; the census finds three further off-line
  quartets below `T = 200` alone. Theorem 8 therefore does not apply to this function,
  and nothing here invokes it.
* **Does its proof transfer under the weaker hypothesis?** Yes, and for a specific
  reason, not by assumption: the proof uses the hypothesis **only** through the single
  inequality `Q̃ ≥ 0` (in `notes/proofs.md`, Theorem B′: "`Q_δ ≥ 0 + 0 − 4δ²F′² − δ⁴F″² −
  |R_4|`"). No other property of the other zeros' location is used anywhere in it.
  So the conclusion
  ```
  λ_min(Q_δ) ≥ −4K(L,γ)δ² − (2/5 + (16/3)e^{2Lδ})L⁵δ⁴
  ```
  holds verbatim with "all other zeros on the line" replaced by "`Q_rest ≥ 0` on the
  admissible sector" — strictly weaker, and the hypothesis this phase measures.
* **What does NOT transfer.** That a measured positive finite restriction establishes that
  hypothesis. It does not, and **R6 is not reported anywhere in this extension**. R6 would
  require control of the infinite-dimensional complement *and* of its coupling to the
  finite block; a positive Galerkin matrix supplies neither.

---

## 11. Phases F and G

### Masking (`results/phaseF.txt`, `data/masking.csv`)

**(A) Fixed-witness removal — additive** (`Q_F(g_*)/‖g‖² = −0.679071346441480591`):

| orbit | `Q_O(g_*)/‖g‖²` | `Q_F − Q_O` | still negative? | effect |
|---|---|---|---|---|
| `O_1` | `−0.985198326439305482` | `+0.306126979997824891` | **NO** | the whole signal |
| `O_2` | `+0.003155790545691441` | `−0.682227136987172032` | yes | **masks** |
| `O_3` | `−0.000006698888947369` | `−0.679064647552533222` | yes | reinforces |
| `O_4` | `+0.000301403207646796` | `−0.679372749649127387` | yes | **masks** |

Sum of all four: `−0.981747831574914613`; additivity residual exactly `0`.
The signs are reported as they come: `O_2` and `O_4` **mask**, `O_3` **reinforces**. Not
every change is called masking.

**(B) Reoptimised — not an additive attribution.** Removing `O_1` and reoptimising in the
same space gives `λ_min = +1.21 × 10⁻²⁴`: **no other negative direction is found**. That
is a statement about this test space and nothing more — failure to find another witness
does not restore whole-window positivity.

**(C) Finite-spectrum hybrid.** `O_1` kept, `O_2…O_4` moved to `Re s = ½` at their
ordinates with multiplicity preserved (quartet → two **double** on-line zeros):

```
fixed witness : −0.679273022625451886   (was −0.679071346441480591), change −2.02e-4
reoptimised   : λ_min = −0.679273139266006236
```

Removing the other enumerated orbits **deepens** the negative, confirming that on balance
they were masking. This is called a **finite-spectrum hybrid** and nothing more: the
unenumerated off-line zeros above the census height are untouched, so it is **not** "all
other zeros moved to the line", and no Euler product or analytic function is attached to
it.

### Sensitivity (`results/phaseG.txt`, `data/sensitivity.csv`)

**(A) Actual-function behaviour.** The reference scale `L_ref = ½ log(5γ_1/2π) =
2.111203034664309`, against a measured macroscopic onset between `L = 1.8` and `L = 2.0`
and a marginal one at `L = 1.8`. This is quoted as a **reference scale, not a predicted
onset**; its provenance is a bandwidth heuristic that `SOURCE_AUDIT.md` E3 records as
**not verified**. Note `δ_1 L = 0.617` at `L = 2`: the actual displacement is nowhere near
a small-displacement regime, and no small-`δ` law is used in Phase D.

**`K` identities**, both by quadrature and in closed form, at `γ = γ_1`:

| `L` | `K_even` (quad) | `K_even` (closed) | `K_odd` (quad) | `K_even+K_odd − 2L³/3` |
|---|---|---|---|---|
| 1.4 | `0.9039859935888327` | `0.9039859935888327` | `0.9253473397445006` | `1.6e-61` |
| 1.8 | `1.932589242028968` | `1.932589242028968` | `1.955410757971032` | `0.0` |
| 2.0 | `2.675097759734244` | `2.675097759734244` | `2.658235573599090` | `6.2e-61` |
| 2.6 | `5.876354553660202` | `5.876354553660202` | `5.840978779673131` | `−1.2e-60` |

`K_even + K_odd = 2L³/3` confirmed identically in `γ`.

**(B) Controlled perturbation** — `Q_t = Q_rest + Q_{O_1}(tδ_1, γ_1)`, `Q_rest` held
fixed, multiplicities preserved, `t = 0` being the collapsed configuration (two double
on-line zeros). This is an explicitly **synthetic deformation**: the Davenport–Heilbronn
function is not deformed and nothing at `t < 1` is an `L`-function.

The **complete** quadratic coefficient, derived in `DERIVATION.md` §4 by expanding (4.2):

```
dQ_O/d(δ²) = 4[Re(m_0 conj m_2) − |m_1|²] = −4[ |F'(γ)|² + Re(F(γ) conj F''(γ)) ],
```

which for real even `g` is exactly the `−4(F'² + FF'')` of `notes/proofs.md` §B.1. At the
frozen witness the **three different numbers** the brief insists on separating are:

```
COMPLETE coefficient   4[Re(m_0 conj m_2) − |m_1|²]  =  −9.64616003587294
ISOLATED term          −4|F'(γ_1)|²                   =  −9.623582778310451
UNIFORM bound          −4 K(L, γ_1)                   = −10.63294229439636
```

with `|F(γ_1)| = 0.01440509`, `|F'(γ_1)| = 1.55109500`, `|F''(γ_1)| = 0.39182780`.

| `t` | `Q_t(g_*)/‖g‖²` | quadratic model | residual | residual`/t⁴` | `λ_min(Q_t)` |
|---|---|---|---|---|---|
| 1.0 | `−0.679071346441` | `−0.611192015467` | `−6.788e-2` | `−0.06788` | `−0.679071` |
| 0.8 | `−0.308148322531` | `−0.280658367613` | `−2.749e-2` | `−0.06711` | `−0.312340` |
| 0.6 | `−0.032198403258` | `−0.023576641503` | `−8.622e-3` | `−0.06653` | `−0.047758` |
| 0.4 | `+0.158360760485` | `+0.160053162860` | `−1.692e-3` | `−0.06611` | `2.05e-24` |
| 0.2 | `+0.270125667806` | `+0.270231045478` | `−1.054e-4` | `−0.06586` | `2.05e-24` |
| 0.1 | `+0.297768936222` | `+0.297775516133` | `−6.580e-6` | `−0.06580` | `2.05e-24` |
| 0.05 | `+0.304661222649` | `+0.304661633797` | `−4.112e-7` | `−0.06578` | `2.05e-24` |
| 0.0 | `+0.306957006351` | `+0.306957006351` | `−3.4e-26` | — | `2.05e-24` |

**The derived law is confirmed with its remainder.** `residual/t⁴` converges to
`−0.0658`: the remainder is exactly `O(δ⁴)`, as the derivation says, and the quadratic
model is accurate to `4 × 10⁻⁷` at `t = 0.05` and to only `7 × 10⁻²` at `t = 1` — because
`δ_1 L = 0.617` is far outside the regime where it applies. **The small-displacement
approximation is reported with its remainder and is used nowhere else.**

The last column makes the brief's warning concrete: the **eigenvalue** shift is not the
fixed-vector shift. `Q_t(g_*)` turns positive between `t = 0.6` and `t = 0.4`, while
`λ_min(Q_t)` is still exactly at its near-null floor `2.05 × 10⁻²⁴` from `t = 0.4`
downward and only departs from it at `t = 0.6`. An eigenvalue shift is not the isolated
negative derivative term and is not `4Kδ²`.

**(C) Overlaps** with the isolated sensitivity optimiser `u cos(γu)` (odd sector):

```
frozen witness g_*        : 0.9513529401
Q_F ground state          : 0.9513529401
Q_rest ground state       : 9.07e-12
```

The detecting witness is, to 95%, the extremiser of the isolated derivative-evaluation
functional `f ↦ F'(γ_1)` — which is why the sensitivity picture and the detection picture
agree here. The `Q_rest` ground state is orthogonal to it, as it must be: it is the
near-null direction of a different form. A small overlap would not have refuted the
sensitivity identity; `u cos(γu)` optimises the isolated functional, not the full form.

---

## 12. Verdict

### The eight questions

**1. Was `F` and its functional equation independently validated?** Yes. `κ` is derived
from the Gauss sum and matches its closed form; the root-number identity holds exactly;
three evaluation routes agree, one of them (Euler–Maclaurin, written here) sharing no
special-function library code with the others; the completed equation holds with a
scale-aware relative residual of `7 × 10⁻⁵⁸` at 60 digits and `8 × 10⁻⁹⁹` at 100 digits.
Entireness and the zero-containing strip `−0.3952 < Re s < 1.3952` are T1.

**2. Did the arithmetic and spectral/contour evaluations agree?** Yes, **PASS-N**.
Arithmetic vs contour: `1.7 × 10⁻¹³`. Arithmetic vs zero-side: `≤ 5.6 × 10⁻¹¹` over six
cases. The two exact-algebra controls agree to `10⁻⁶¹`. Not PASS-C because the quadrature
and the tail density are not rigorous bounds.

**3. Is there an explicit negative witness, and is it T2-N or T2-C?** Yes, and it is
**T2-C**: `Q_F(g_*)/‖g_*‖² = −0.6790713464414805909863775860931317559671…` with an Arb
ball radius of `1.4 × 10⁻¹¹⁶`, at `L = 2.0`, odd sector, `N = 37`, every coefficient an
exact 30-digit decimal, no quadrature and no truncation anywhere in the computation. The
*smallest tested window with a detected negative witness* is `L = 1.8` in both sectors;
that one is **T2-N**, reported and not certified.

**4. What remains after subtracting the selected quartet?** `Q_rest(g_*) = +0.3061269800`,
certified. `Q_{O_1}(g_*) = −0.9851983264` carries the whole of the negativity and more;
the rest of the spectrum gives back `31%` of it.

**5. Is that conclusion pointwise, finite-dimensional, or whole-window?** The certified
`+0.3061` is **pointwise**. The `λ_min(Q_rest) > 0` of §10 is **finite-dimensional**, at
one basis level, and is classified **R4** (numerically positive finite restriction) —
*not* R3, because no interval certificate was produced for the restriction, and **not
R6**, which is not reported anywhere. Nothing whole-window is established.

**6. Which contributions mask or reinforce the selected signal?** On the frozen witness,
`O_2` (`+3.16e-3`) and `O_4` (`+3.01e-4`) **mask**; `O_3` (`−6.70e-6`) **reinforces**. On
balance the other enumerated orbits mask: moving them to the line deepens the witness by
`2.0 × 10⁻⁴`. Removing `O_1` and reoptimising finds no further negative direction in this
space (`λ_min = +1.2 × 10⁻²⁴`), which is a statement about the space, not about
positivity.

**7. Does the controlled perturbation obey the derived sensitivity law?** Yes, with its
remainder. `residual/t⁴ → −0.0658`, confirming the `O(δ⁴)` tail of the complete quadratic
coefficient `−9.6462`, which is distinct from the isolated `−4|F'|² = −9.6236` and from
the uniform `−4K = −10.6329`. At the actual displacement `δ_1 L = 0.617` the quadratic
model is off by 11%, so the law describes the deformation, not the detection.

**8. What was reproduced, and what is new only here?** *Reproduced:* the existence of
off-line Davenport–Heilbronn zeros near the published values (seed matched to 14 places);
the repository's `γ = 0` real-pair parity split, recovered independently as the `γ → 0`
limit of the quartet form; the repository's `−4(F′²+FF″)` quartet expansion, recovered as
the real-even case of the general complex coefficient; the committed upstream `ζ` and
mod-5 regressions, unchanged. *New only to this investigation:* nothing is claimed as a
first. At least two contemporaneous documents work on the same control, one of them with
an interval certificate (`SOURCE_AUDIT.md` §F), and the window-convention ambiguity
between them is unresolved. What this extension owns is an independent implementation
inside this repository's conventions, with its own derivation, census, gates and
certificate.

### What surprised me

1. **How violent the onset is.** `λ_min` falls seventeen orders of magnitude from
   `3.9 × 10⁻⁵` at `L = 1.0` to `2.7 × 10⁻²²` at `L = 1.6` *while staying positive*, is
   negative by only `10⁻²¹`–`10⁻¹⁶` at `L = 1.8`, and is negative by `0.5`–`0.7` at
   `L = 2.0`. I expected a smooth crossing and budgeted precision for resolving a sign near
   zero; what the window actually does is sit on a knife edge for most of the range and
   then fall off it.
2. **`Q_rest` does not merely stay positive — it sits exactly on zero.** Once detection
   happens, `λ_min(Q_rest)` is `10⁻²⁴`, `10⁻²⁸`, `10⁻³¹`, `10⁻³⁶` at `L = 2.0 … 2.6`,
   decreasing while remaining positive. The remaining-zero form is pinned to the boundary
   of positivity, and all of the negativity is one quartet. I had expected the other three
   enumerated off-line quartets to contribute visibly; their total effect on the witness is
   `3.5 × 10⁻³`, against `0.985` from `O_1`.
3. **The witness is essentially the textbook sensitivity extremiser.** Overlap `0.9514`
   with `u cos(γ_1 u)`. The brief warns that a small overlap would not refute the
   sensitivity identity; the surprise is that the overlap is large, so the two pictures
   coincide here even though `δ_1 L = 0.617` is far outside the perturbative regime where
   they are supposed to.
4. **A ladder that agrees with itself to 85 digits can be wrong in the 13th.**
   `mp.mpf('1.8')` built before the precision was raised gave three "independent"
   precisions that all evaluated the same slightly-wrong window length
   (`RETRACTIONS.md` R6). The only reason it was caught is that two scripts froze the same
   witness and disagreed.
5. **The `b_n` make the missing Euler product visible in the input to the detector.** For
   `ζ` the non-prime-power coefficients vanish; here fifteen of them below `n = 50` do not,
   and `b_6 = (1+κ²)log 6` is exact.

### Open questions (five)

1. Is the `L = 1.8` negativity (`−1.7 × 10⁻¹⁶` odd, `−2.7 × 10⁻²¹` even) stable under
   further basis enrichment, and can the onset be bracketed rather than located between two
   scanned windows? Levels 1→3 lower it monotonically; where it converges is unknown.
2. Does `λ_min(Q_rest) > 0` persist as the basis grows, or is the near-null direction the
   leading edge of a genuine negative one? At `L = 2.6` it is `8 × 10⁻³⁶` — small enough
   that a modest enrichment could change its sign, which would move the Phase E
   classification from R4 to R2.
3. What is the correct window convention in the two located contemporaneous documents? Our
   `L = 2.0` level-3 value `−0.710` is converging toward the `−0.73` one of them reports at
   "`L = 4`", which is consistent with a factor-2 convention and with our finite basis being
   a restriction, but nothing here resolves it.
4. Can the census be made rigorously complete? The winding numbers are integers to
   `10⁻²³` by adaptive sampling; Arb's rigorous integration could in principle certify the
   counts, and then the "numerically reconciled" census of §6 would become a complete one.
5. Is there an off-line quartet with small `δ` and low `γ` that the detector sees *earlier*
   than `O_1`? `O_3` has `δ = 0.074` at `γ = 166`, and the onset depends on both; the scan
   resonates only at `γ_1`.

### Primary verdict

> ## A — certified negative Davenport–Heilbronn witness with independent cross-checks.
>
> One explicitly specified compactly supported test function, `L = 2.0`, odd sector,
> `N = 37`, coefficients stored as exact decimals, with
> `Q_F(g_*)/‖g_*‖² = −0.679071346441480590986377586093…` enclosed in an Arb ball of radius
> `1.4 × 10⁻¹¹⁶` and **strictly negative**. Cross-checked by an independent mpmath
> assembly (44 digits), by a logarithmic-derivative contour integral that never touches the
> `b_n` (`1.7 × 10⁻¹³`), by a zero-side sum over a reconciled census (`≤ 5.6 × 10⁻¹¹`), and
> by two exact-algebra gate controls (`10⁻⁶¹`).

**Phase E remaining-form classification, stated separately and not folded into the
verdict: R4 — a numerically positive finite restriction, at every tested window and
sector.** Not R3 (no interval certificate for the restriction), and **not R6** anywhere.
The single pointwise value `Q_rest(g_*) = +0.3061269800` *is* certified, but a certified
positive number at one point is not a positive form.

### What this does and does not say

**This does not advance RH.** Nothing here is about `ζ`. A detector that sees this
function's off-line zeros is a detector that is not blind; it is not thereby a detector
that would see `ζ`'s.

**Which detector assumptions the control validates.** That the finite-window form
assembled in this repository — archimedean closed form, conductor term, von-Mangoldt-type
arithmetic sum, parity sectors, the quartet orbit convention — is correctly normalised on
a function whose spectrum is *known* to be non-critical: three independent evaluations of
the same quantity agree to ten or more digits, and the form does go negative, at a
magnitude and location the spectrum explains. The same machinery, run on `ζ`, is therefore
not trivially positive by construction.

**Which overly broad positivity claims it rules out.** Any claim that finite-window Weil
positivity follows from a functional equation of this shape plus entireness plus a
zero-containing strip: `F` has all three and its window form is negative by `0.68` at
`L = 2`. And any claim that a positive finite Galerkin restriction of a remaining-zero form
establishes a whole-window hypothesis: here `λ_min(Q_rest) > 0` at every tested window
while `Q_F` is macroscopically negative at five of them, and the positive restriction
establishes nothing about the complement — which is exactly why Theorem 8's hypothesis is
stated as sufficient rather than necessary.

**Which arithmetic distinction from `ζ` remains unresolved.** The one that matters.
`F` differs from `ζ` in that its `b_n` are supported off the prime powers — fifteen
non-prime-powers below 50 — and in that its Dirichlet coefficients are periodic but not
multiplicative. **This experiment does not show that either property is what makes the
form go negative.** It shows only that a function with those properties and with off-line
zeros has a negative window form, and that the negativity at the witness is carried by one
off-line quartet. Disentangling "has off-line zeros" from "has no Euler product" would
require a control with one property and not the other, and this extension does not
construct one.

---

## 13. Final repository check and reproducibility

**Clean-process regressions, before and after** (`results/regression_before.txt`,
`results/regression_after.txt`), reading the expected values **from the committed
upstream CSVs**, never edited to make a run pass:

| regression | rows | worst relative deviation |
|---|---|---|
| parent ζ window, `data/stage3_gateB_lambda_min.csv` | 6 | `4.9 × 10⁻²²` |
| dirichlet-windows `χ_{−3}`, `data/phase3_gateB.csv` | 4 | `4.4 × 10⁻¹⁶` |

Identical before and after, under python-flint being installed into the shared `.venv`.
**ALL REGRESSIONS PASS** in both runs.

**Repository diff.**

```
git diff --stat
  README.md | 2 ++
  1 file changed, 2 insertions(+)

git status --short
   M README.md                                             <- the one pointer line
  ?? extensions/davenport-heilbronn-2026-10-23/            <- this extension
  ?? extensions/newtonian-weil-response-2026-10-08/        <- pre-existing, untouched
  ?? paper/weil_window_control_rev_certified_witness.fdb_latexmk   <- pre-existing
  ?? paper/weil_window_control_rev_certified_witness.fls           <- pre-existing
```

The **only** modification to a tracked file anywhere in the repository is one pointer
line (plus its blank separator) appended to the root `README.md` extension list. The
appender is idempotent: running it again detects the marker and does nothing. Every
other byte of this work is inside `extensions/davenport-heilbronn-2026-10-23/`. The three
pre-existing untracked paths recorded in §1 are unchanged, and no existing source file or
regression fixture was touched.
No branch was switched or reset, nothing was committed, pushed or published, and no
licence was changed.

**Integrity.** `MANIFEST.sha256` lists every file in this extension; `sha256sum -c
MANIFEST.sha256` verifies clean. It contains hashes of this extension's own files only —
**no third-party PDF was downloaded** (`SOURCE_AUDIT.md` §H).

**Configuration.** One recorded file, `src/config.py`: conventions, precisions (60 / 100 /
160 digits), basis definitions and levels, windows, the witness-selection thresholds, and
the numerical budgets. No randomness is used anywhere; `SEED` is recorded in the config so
that the claim is auditable rather than merely asserted. Every script is restartable and
communicates with the others only through `data/*.csv` and `certificates/*.json`; all
substantive loops and all waiting live in script files, not in shell one-liners.

**If resources had run out.** They did not, but the ordering was chosen so that they
could: the certified witness (§9) and its decomposition (§10) were produced before the
masking and sensitivity sweeps, and both would have stood alone. The blockers that remain
are listed as open questions 1–5 in §12, and the one that would most change the record is
question 2 — whether `λ_min(Q_rest)`, already at `8 × 10⁻³⁶` at `L = 2.6`, survives a
larger basis.

---

## 14. The standalone note, and what was added for it

After the record above was complete, a standalone note was written:
[`note/dh_note.tex`](note/dh_note.tex) → `note/dh_note.pdf`, 12 pages,
**doi:10.5281/zenodo.23608288**. It carries the
same non-claims, the same tiers, the same retraction log and a reproduction appendix, and
it states the one result: a certified negative witness with the remaining-zero form
positive. Three things were **added to the extension** while writing it, and they are
recorded here rather than folded silently into the sections above.

**(a) A second certificate, at `L = 1.8`.** The witness that Phase D could only report
(`certificates/witness_smallest_L.json`) was put through `src/certify.py` at 600 bits and
**certifies**:

```
L = 1.8, odd, N = 57, cutoff n ≤ 38 (26 nonzero terms)
Q_f(g)/‖g‖² = −1.68123895523670103326233298940716630655e-16   radius 8.15e-177
STRICTLY NEGATIVE : True
Q_{O_1} = −6.45917960933013e-14,  Q_rest = +6.44236721977777e-14
```

So §8's "reported, not certified" is **superseded for the value**: the negativity at
`L = 1.8` is now T2-C in both sectors' deeper one. What is still *not* certified is the
**onset** — that would need a certified lower bound at `L = 1.6` over a whole space,
which is a different and harder object than a certified value at a point. The note's
non-claim 5 says exactly this, and the decomposition at `L = 1.8` adds one fact worth
keeping: there the cancellation between the quartet and the rest is **99.7397%**, against
**31.07%** at `L = 2.0`. The witness is small near the onset because the two large terms
nearly annihilate, not because either is small.

**(b) A term-level sign-flip control.** `src/gate_signflip.py` →
`data/gate_signflip.csv`. Negating each of the archimedean, conductor and arithmetic
terms of `M_DH = Arch − LogPi + LogQ − Prime` in turn must move `Q(g)` by exactly
`−2 Arch(g)`, `−2 LogQ(g)`, `+2 Prime(g)`. Worst residual over all rows and both sectors:
`3.11 × 10⁻⁶¹`. The conductor row is the one that would catch a wrong conductor, and it
reproduces `LogQ(g) − log(5)·h(0) = 0` exactly. This is a third exact-algebra control
alongside §7's A and B; §7's classification is unchanged (**PASS-N**), since what keeps it
off PASS-C is the quadrature and the measured tail density, neither of which this touches.

**(c) A number-level verifier for the note.** `note/verify_note_numbers.py` rebuilds all
217 printed quantities from `data/`, `results/` and `certificates/`, scans the LaTeX for
numeric literals, and fails if a literal is not derivable, disagrees with its source, or
is printed to **more significant digits than its source carries**. It also regenerates the
four machine-generated tables and compares them, cross-checks the Gate A values against
this file, and requires every number in the literature appendix and the bibliography to
appear verbatim in `SOURCE_AUDIT.md`.

It earned its place immediately: the note's first draft had **nine tables and scalars
truncated where they should have been rounded** — `0.9513529400` for `0.9513529401`,
`1.9363560766210387` for `1.93635607662104`, and so on — plus one genuinely wrong value,
the strip half-width printed as `0.8952` when `σ₀ − ½ = 0.895136`. None of these would
have been caught by reading. All numeric tables in the note are now generated from the
data by `note/make_tables.py` and `note/make_witness_table.py` rather than typed.

**No result in §§1–13 changed.** The new certificate strengthens a tier, the new control
adds a row, and the verifier found only transcription errors in an unpublished draft.
`SOURCE_AUDIT.md` gained §I, recording the references the note cites that the main record
did not.
