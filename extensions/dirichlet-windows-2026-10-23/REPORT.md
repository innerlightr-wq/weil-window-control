# Dirichlet windows: transposing the finite-window Weil detector from ζ to L(s, χ)

Extension of `weil-window-control`. Preregistered the same way as the other extensions:
non-claims first, epistemic tiers on every numbered claim, a retraction log, regression
values in `data/`, and a ζ control reproduced in the same run as every χ measurement.

**Carries no DOI. Not in the paper. Not reviewed.**

Started 2026-10-23. Written incrementally, phase by phase; a partial file is a partial run,
not a draft of a finished one.

---

## 0. Non-claims

Written before any computation in this directory, and not revised afterwards.

1. **No claim about GRH or RH.** Nothing here restricts the location of any actual zero of
   any L(s, χ) or of ζ. The detector is calibrated, not applied.
2. **No claim about the Siegel-zero problem.** Phase 5 measures what the *detector* does
   when a real zero is planted. A statement about the detector's response is not a statement
   about L(s, χ). The table in Phase 5 is labelled accordingly and must not be read as
   evidence for or against the existence of Siegel zeros.
3. **The October 2026 OpenAI quasi-Riemann-hypothesis release is not used as input
   anywhere, and is not cited as established.** It is present on this machine
   (`~/Downloads/paper (1).pdf`, *The Quasi-Riemann Hypothesis: A Zero-Free Half-Plane
   ℜs > 7/8*, OpenAI, 30 September 2026), and its abstract covers exactly the family of
   L-functions studied here. It was deliberately not opened past the title page for this
   work and no quantity in this directory depends on it. The existing audit of that release
   is `extensions/openai-quasi-rh-audit-2026-10-07/` (verdict: NO BRIDGE FOUND); that
   verdict is likewise not used as input here.
4. **Planted spectra are fabricated diagnostics**, not L-functions of anything. A planted
   quartet or a planted real zero is a counterfactual added to, or substituted into, a zero
   multiset; it is not a claim that such a zero exists.
5. **Numerical detections are tier T2 unless certified over ℚ**, and **non-detections
   certify nothing**: a finite basis that finds no negative direction is not evidence of
   positivity. Where a non-detection is structural rather than numerical this is stated and
   proved, and is then labelled T1 as a statement about the detector only.
6. **Zeros of L(s, χ) computed here are T2 unless independently cross-checked.** Ordinates
   are produced by root-finding on a Hardy-type Z function assembled from Hurwitz zeta
   values; the assembly is this extension's own code and has no external certificate.
7. **No new spectral architecture.** The finite-window Weil form is Yoshida (1992) and
   Bombieri (2000); the Dirichlet explicit formula is classical (Weil; see
   Iwaniec–Kowalski Ch. 5, Montgomery–Vaughan Ch. 12). What is new here, if anything, is
   only the calibration.
8. **`λ*(0.8) ≥ 8.9e−18` is Zhu's certified lower bound for ζ** and is used only on its
   documented range and normalisation. **No analogue is available for L(s, χ)**: no
   certified positivity lower bound for any Dirichlet window is quoted or assumed. Every
   χ-side baseline in this directory is a *variational upper bound* on λ*_χ(L), i.e. a
   computed λ_min at finite N, and is never used as a lower bound.
9. **Not attempted:** complex characters; non-primitive characters; imprimitive conductors;
   any Lean formalisation; any exact-rational interval certificate on the χ side.

## 0.1 Tier key

Same as the parent repository.

| tier | meaning |
|---|---|
| **T1** | proved here, or cited from a source verified in full text; arithmetic residuals quoted |
| **T2** | measured numerically at stated precision with a precision-doubling or noise-floor check |
| **T3** | observed pattern, fit, or heuristic; no error control claimed |

## 0.2 Retraction log

Corrections to **this directory** and to statements made *about* this work, appended as they
are found. Empty is a real state and is reported as such.

| # | what was claimed | status | where |
|---|---|---|---|
| — | *(none yet)* | | |

## 0.3 Corrections to the prompt that commissioned this extension

Recorded here because they change what the phases can claim, not to score points.

| # | claim in the brief | finding |
|---|---|---|
| P1 | "The paper's Remark 2.4 predicts the window is structurally blind [at γ = 0] because K(L, 0) = 0" | **There is no Remark 2.4 in the paper**, and no remark anywhere in `paper/weil_window_control.tex` about γ = 0 or about real zeros (`grep -n "gamma=0\|real zero\|Siegel\|blind"` returns only a bibliography line). The *mathematical* content is real but is a consequence of equation (K) in Theorem 8, not of a remark: `K(L, 0) = ∫_{−L}^{L} u² sin²(0·u) du = 0`. Phase 5 therefore has no prediction in the paper to confirm or contradict; it establishes the γ = 0 behaviour from scratch. See Phase 5, where the correct statement turns out to be **strictly stronger** than "the δ² coefficient vanishes". |

---

*(Phases appended below as they complete.)*

# Phase 1 — the window form for L(s, χ)

Derivation: [`notes/dirichlet_window.md`](notes/dirichlet_window.md).
Implementation: [`src/dirichlet_window.py`](src/dirichlet_window.py).
Self-test: [`src/selftest.py`](src/selftest.py) → [`results/phase1_selftest.txt`](results/phase1_selftest.txt).

## 1.1 Characters

| χ | disc. d | q | parity a | χ(−1) | values χ(0..q−1) | ε(χ) |
|---|---|---|---|---|---|---|
| χ₋₃ | −3 | 3 | 1 (odd) | −1 | 0, 1, −1 | +1 |
| χ₋₄ | −4 | 4 | 1 (odd) | −1 | 0, 1, 0, −1 | +1 |
| χ₅ | 5 | 5 | 0 (even) | +1 | 0, 1, −1, −1, 1 | +1 |
| χ₋₇ | −7 | 7 | 1 (odd) | −1 | 0, 1, 1, −1, 1, −1, −1 | +1 |
| χ₈ | 8 | 8 | 0 (even) | +1 | 0, 1, 0, −1, 0, −1, 0, 1 | +1 |

The five are the Kronecker symbols of the fundamental discriminants −3, −4, 5, −7, 8 —
i.e. the complete list of real primitive characters of conductor ≤ 8. Parity is
`a = 0 ⟺ d > 0`. **(T1)** Each is verified complete multiplicative, exactly periodic mod q,
supported off `gcd(n,q) > 1`, primitive (no proper divisor of q is a period), and of the
stated parity; the Gauss sum is `τ(χ) = i^a √q` to 10⁻⁴⁰, hence the root number is
**ε(χ) = +1 for all five**, which is what makes `Λ(½+it, χ)` real and a Hardy `Z` function
available.

## 1.2 The result

> **(D1) [T1]** For `f` real, even, supported in `[−L, L]`, and χ real primitive mod q of
> parity a,
> ```
> Q_χ(f) = (1/2π) ∫ |F(t)|² Ψ_{L,χ}(t) dt,
> Ψ_{L,χ}(t) = Re ψ((1/4 + a/2) + it/2) + log(q/π)
>              − Σ_{log n < 2L} 2 Λ(n) χ(n) n^{−1/2} cos(t log n),
> ```
> with **no pole term**, and `Q_χ(f) = Σ_ρ |F(γ_ρ)|²` if every non-trivial zero is on the line.

**Exactly four changes from ζ**, and the two that the brief flagged as sign/factor risks
resolve as follows.

| | ζ | L(s, χ) | risk resolved |
|---|---|---|---|
| pole | `+2F(i/2)²` | **absent** | — |
| conductor | `−log π` | `+log(q/π)` | coefficient is **1, not ½** (the ½ of `½log(q/π)` in `Λ'/Λ` is cancelled by the two contour lines, exactly as the ½ of `½ψ(s/2)` is cancelled by `½[ψ(z)+ψ(z̄)] = Re ψ(z)`); sign is **positive**, so the conductor pushes the form **up** |
| archimedean | `Re ψ(1/4 + it/2)` | `Re ψ((1/4 + a/2) + it/2)`, i.e. `1/4` even / **`3/4` odd** | in u-space this is exactly two substitutions: `z : 1/4 − iα/2 → (1/4 + a/2) − iα/2` and the truncation prefactor `e^{−U/2} → e^{−(1/2+a)U}` (D2) |
| prime | `Λ(n) n^{−1/2}` | `Λ(n) χ(n) n^{−1/2}`, `p ∣ q` dropped | weights now take **both signs** |

Everything is in u-space (Parseval); the cross-correlation machinery `C_{jk}` is reused
verbatim from the parent `src/zeta_window.py`, so **no oscillatory quadrature enters** and no
basis-level formula was re-derived. Both test-function parity sectors are given; note that
**test-function parity (`sector`) and character parity (`a`) are independent**, and that the
ζ code's `pole_sign(sector)` has no analogue here — there being no pole, the odd sector
differs from the even one only through the basis.

## 1.3 Self-test (all six groups pass)

| # | check | external reference | result |
|---|---|---|---|
| 1 | character axioms incl. primitivity | — | pass, all five **(T1)** |
| 2 | `τ(χ) = i^a √q`, hence `ε(χ) = 1` | classical Gauss sum | `\|τ − i^a√q\| < 10⁻³⁰` **(T1)** |
| 3 | `L(1, χ_d)` vs the **class-number formula** | `2πh/(w√\|d\|)` (d<0), `2h log ε_d/√d` (d>0) | rel. err ≤ 1.9e−41 **(T1)** |
| 4a | `L(3, χ₋₄) = π³/32` (Dirichlet β) | exact | agrees to 24 digits **(T1)** |
| 4b | `L(0, χ_d) = 2h(d)/w` (d<0), `= 0` (d>0) | exact **rationals** 1/3, 1/2, 1, 0, 0 | exact **(T1)** |
| 4c | Dirichlet series at `s = 6`, `N = 4000`, proved tail `2qN⁻⁶` | — | within the bound, all five **(T1)** |
| 4d | `Λ(s, χ) = Λ(1−s, χ)` at 3 off-line complex points | functional equation | rel. err ≤ 3.7e−41 **(T1)** — tests `q`, `a` and `ε` simultaneously |
| 5 | `arch_closed_a(·, ·, a=0)` **≡** parent `arch_closed(·, ·)` | parent repository | **max abs diff = 0.0** over 3 sectors × 3 L × 10 pairs; and `a=1` differs by 0.504, so the parity does enter **(T1)** |
| 6 | `\|Im Z_χ(t)\| / max(\|Re Z_χ\|,1)` | functional equation + `ε = 1` | ≤ 2.7e−40 at `t = 1, 5.5, 12.3, 30.7`, all five **(T2)** |

Check 5 is the one that matters structurally: the archimedean generalisation is a **strict
extension** of the parent code — byte-for-byte agreement at `a = 0`, not agreement to a
tolerance — so an archimedean bug on the χ side cannot be a regression on the ζ side, and the
ζ control retained in every later phase will catch it.

Check 3 deserves emphasis as the strongest external anchor available without a network: the
class-number formula pins `L(1, χ_d)` to 41 digits against a closed form containing `h(d)`,
`w`, and the fundamental unit — none of which appears anywhere in this code.

# Phase 2 — Gate A

Scripts: [`src/zeros.py`](src/zeros.py), [`src/zeros_summary.py`](src/zeros_summary.py),
[`src/lmfdb_check.py`](src/lmfdb_check.py), [`src/phase2_gateA.py`](src/phase2_gateA.py),
[`src/phase2_gateA_primeflip.py`](src/phase2_gateA_primeflip.py).
Transcripts: [`results/phase2_zeros.txt`](results/phase2_zeros.txt),
[`results/phase2_lmfdb.txt`](results/phase2_lmfdb.txt),
[`results/phase2_gateA.txt`](results/phase2_gateA.txt).
Data: `data/zeros_<χ>.txt`, [`data/phase2_gateA.csv`](data/phase2_gateA.csv).

## 2.1 Critical-line ordinates **(T2, cross-checked)**

Root-finding on the Hardy-type `Z_χ(t) = e^{iθ_χ(t)} L(½+it, χ)` (real because `ε(χ) = 1`),
adaptive grid at `steps_per_gap = 8`, bisection then `mp.findroot`, `dps_root = 60`,
`T_max = 275`.

| χ | q | a | #ordinates in (0, 275) | γ₁ (40 digits computed, 30 shown) |
|---|---|---|---|---|
| χ₋₃ | 3 | 1 | 170 | 8.03973715568146668171362321417 |
| χ₋₄ | 4 | 1 | 183 | 6.02094890469759665490251152161 |
| χ₅ | 5 | 0 | 193 | 6.64845334472771471612327845998 |
| χ₋₇ | 7 | 1 | 207 | 4.47573828372868313197462848719 |
| χ₈ | 8 | 0 | 213 | 4.89997399700703650103830489920 |

**Completeness.** `found(T) − θ_χ(T)/π` over the ladder `T = 10, 25, …, 275`:
`max |·| = 0.57, 0.57, 0.69, 0.62, 0.61` and total band `1.20, 1.02, 1.18, 1.00, 1.21`
respectively. A missed zero would make this column step down by **1** and stay there; it
oscillates inside a band of width ≈ 1 with no step. **No zero is missed and none is
spurious**, over 966 ordinates.

**LMFDB cross-check — done, the network was available.** For each conductor, every degree-1
primitive L-function was pulled from the LMFDB API and the character was matched **by Euler
factors** (LMFDB's second local coefficient is `χ(p)`; verified on the unambiguous
conductor-3 case) at **19 primes**, with the parity in the label required to equal our `a` —
no Conrey-label guessing. Identification was unique in all five cases.

| χ | LMFDB label | ordinates compared | max abs diff | max rel diff |
|---|---|---|---|---|
| χ₋₃ | `1-3-3.2-r1-0-0` | 114 | 7.00e−28 | 3.51e−30 |
| χ₋₄ | `1-2e2-4.3-r1-0-0` | 122 | 8.00e−28 | 4.02e−30 |
| χ₅ | `1-5-5.4-r0-0-0` | 129 | 4.00e−28 | 2.01e−30 |
| χ₋₇ | `1-7-7.6-r1-0-0` | 140 | 9.00e−28 | 4.52e−30 |
| χ₈ | `1-2e3-8.5-r0-0-0` | 144 | 8.03e−29 | 4.04e−31 |

**649 ordinates** agree with LMFDB's published values to `9.0e−28` worst case — LMFDB stores
about 30 digits, so this is at their truncation level, not ours. The brief asked for the
first five ordinates; all available overlapping ordinates were compared instead. The raw
LMFDB values are kept in [`results/phase2_lmfdb.json`](results/phase2_lmfdb.json) so the
comparison is re-checkable offline. **This is substantially better than the T2 the
non-claims allow for, and the ordinates are used below as cross-checked.**

## 2.2 The predicted decay rate is basis-dependent **(T1)**

Before reading Gate A: `discrepancy(K) = max_{j,k} |M_χ[j,k] − Σ_{n≤K} 2F_j(γ_n)F_k(γ_n)/‖·‖|`
is governed by the tail `Σ_{|γ|>γ_K} |F(γ)|²`, hence by how fast the **basis** transforms
decay — not by the L-function.

| basis | value at ±L | `F_k(t)` | `\|F_k\|²` | predicted exponent |
|---|---|---|---|---|
| `EVEN_D` | 0 (continuous) | `O(t⁻²)` | `O(t⁻⁴)` | **−3** × `log qγ_K` |
| `ODD` | 0 (continuous) | `O(t⁻²)` | `O(t⁻⁴)` | **−3** × `log qγ_K` |
| `EVEN` | `(−1)^k ≠ 0` (**jump**) | `O(t⁻¹)` | `O(t⁻²)` | **−1** × `log qγ_K` |

The brief's `γ_K⁻³` is stated "for modes vanishing at ±L", which is the first two rows. The
`EVEN` basis is expected at **−1** and a shallow exponent there is not a failure. All three
rates are read off the ζ control in the same run.

## 2.3 Gate A: **PASS**, all five characters × three sectors × three L **(T2)**

Fitted exponent of `discrepancy(K)` against `γ_K`, over `K = 20, 40, 80, 160`
(`γ_160 ≈ 190–230` for χ, `334` for ζ), `N = 4`, `dps = 40`:

| | `EVEN_D` (pred. −3) | `ODD` (pred. −3) | `EVEN` (pred. −1) |
|---|---|---|---|
| **ζ control** | −2.74, −2.60, −2.67 | −2.77, −2.83, −2.74 | −0.773, −0.784, −0.775 |
| χ₋₃ | −2.90, −2.92, −2.79 | −2.86, −2.66, −2.76 | −0.823, −0.753, −0.798 |
| χ₋₄ | −2.91, −2.68, −2.71 | −2.94, −2.95, −2.86 | −0.846, −0.878, −0.834 |
| χ₅ | −2.91, −2.80, −2.70 | −2.96, −2.83, −2.87 | −0.854, −0.826, −0.838 |
| χ₋₇ | **−3.02**, −2.98, −2.87 | −2.96, −2.69, −2.75 | −0.860, −0.767, −0.811 |
| χ₈ | **−3.01**, −2.99, −2.91 | −2.99, −2.68, −2.70 | −0.868, −0.765, −0.791 |

(three entries per cell = `L = 0.5, 0.8, 1.0`.)

Every χ exponent is inside the spread the ζ control itself shows, in every sector. Absolute
discrepancies at `K = 160` run `1.3e−5` to `1.6e−4` in `EVEN_D`/`ODD` — the same order as the
parent repository's published ζ Gate A (`1.6e−7` at `L = 0.8`, `N = 4`, `K = 800`; we stop at
`K = 160` because that is where the computed χ ordinates stop).

**Note the ordering.** The discrepancy at fixed `K` grows monotonically with `q`
(`EVEN_D`, `L = 0.8`: `6.1e−6` for ζ, then `1.32e−5`, `1.96e−5`, `2.09e−5`, `2.49e−5`,
`2.67e−5` for `q = 3,4,5,7,8`). That is the `log(qγ_K)` factor in the tail, i.e. the higher
zero density of a larger conductor — a prediction of §2.2 confirmed, not an anomaly.

## 2.4 Sign-flip controls: the test **does** discriminate **(T2)**

Ratio of the `K = 160` discrepancy with exactly one term negated, to the correct one:

| | flip `LogQ` (conductor) | flip `Arch` (archimedean) |
|---|---|---|
| `EVEN_D`/`ODD`, all χ, all L | **×2.6e4 … ×4.4e5** | **×3.1e4 … ×4.4e5** |
| `EVEN`, all χ, all L | ×77 … ×212 | ×118 … ×231 |

Both terms are discriminated everywhere, by four to five orders of magnitude in the bases
that vanish at ±L. The smaller (but still decisive) ratios in the `EVEN` basis are because
its *correct* discrepancy is already `2e−2 … 1.6e−1`, per §2.2. **Gate A is therefore
informative for the conductor term and for the archimedean term separately**, which is
exactly what the brief flagged as the Phase-1 risk; the derivation's two risky constants
(`log q` with coefficient 1 and positive sign; the `1/4 → 1/4 + a/2` digamma shift) are
confirmed numerically and were not adjusted to fit.

## 2.5 Verdict

> **Gate A passes.** The window form (D1) is validated against `Σ_ρ |F(γ_ρ)|²` for every one
> of the five characters, in both test-function parity sectors and in a third (discontinuous)
> basis, at three window lengths, with the predicted decay rate and with the conductor and
> archimedean terms separately discriminated. Phases 3–5 proceed.
