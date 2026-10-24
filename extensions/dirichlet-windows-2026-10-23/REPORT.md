# Dirichlet windows: transposing the finite-window Weil detector from ζ to L(s, χ)

Extension of `weil-window-control`. Preregistered the same way as the other extensions:
non-claims first, epistemic tiers on every numbered claim, a retraction log, regression
values in `data/`, and a ζ control reproduced in the same run as every χ measurement.

**The extension record carries no DOI and is not in the paper. Not reviewed.** The
standalone note drawn from it does carry one: [`note/siegel_note.pdf`](note/siegel_note.pdf),
**doi:10.5281/zenodo.23623473** (§7).

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
   **Added 2026-10-23, and load-bearing:** whether Zhu's infimum is taken over the **even
   sector** or over the **full admissible space** is *not* established here — Zhu's paper was
   not obtained or read for this extension, and the parent paper deliberately leaves the
   question open (`paper/weil_window_control.tex:312–316`). §5.3 shows the two cases have
   different consequences for real zeros, so **no claim is made in either direction**. In
   particular this extension does **not** assert that Zhu's certificate is uninformative
   about the Siegel configuration.
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
| **1** | `src/regression_zeta.py` asserted the ζ-side values against the parent's committed CSVs at relative tolerance `1e−20`, and **reported 14 failures**. | **RETRACTED — the tolerance was wrong, the values are right.** The parent writes eigenvalues as `mp.nstr(lam, 12)`, i.e. **twelve significant digits**, so a full-precision recomputation can agree with the stored *decimal string* only to about `1e−12`. Every "failure" sat at `4e−15 … 1.7e−12`, exactly the parent's printing precision. Tolerance corrected to `5e−12` (half a unit in the parent's last printed digit) with the reason recorded in the file; the suite then passes. **No computed value changed.** The binding precision on any regression against this repository's CSVs is the CSVs' own. | `src/regression_zeta.py`, `results/regression_zeta.txt` |
| **2** | `src/phase2_gateA_primeflip.py` built `L = mp.mpf('0.8')` *before* setting `mp.mp.dps`, so `L` carried only 15-digit precision into a 40-digit run (printing as `0.80000000000000004440…`). | **FIXED.** Precision is now set first. The affected quantities were control *ratios* of order `10⁴`, so no reported conclusion moved — the re-run reproduced the same baseline discrepancy `1.31638e−5` — but the ordering was wrong and is recorded rather than quietly repaired. | `src/phase2_gateA_primeflip.py` |
| **3** | `src/planting.py::quartet_block` initially omitted the factor `h_sign(sector) = −1` for the odd basis, where `zeta_window.F` returns `F̃ = F/i`. | **FIXED before use**, caught by checks (E) and (F) of `selftest_planting.py` (`quartet(γ,0) = 2·pair(γ)` failed by 4.04). No result was computed with the bad version. The parent repository's `stage3e_certify.move_block` has the same latent gap, **harmlessly**, because it is only ever called on the even sectors where `h_sign = +1`; documented in `planting.py::h_sign` and **not** a correction to any parent result. | `src/planting.py`, `src/selftest_planting.py` |
| **4** | The Phase-6 commit message listed a fourth retraction — the `stage3c_even_odd` Connes window — but **the row was never added to this table**. | **FIXED (2026-10-23).** The substance: that file's fourth window is `L = log(13)/2` stored rounded to 10 digits, and `λ_min ≈ 10⁻²⁹` there is hypersensitive because `2L` sits exactly at the `n = 13` prime-power cutoff, so reconstructing `L` from the stored string moved the answer by `1.8e−8` relative. `regression_zeta.py` now substitutes the exact `mp.log(13)/2`. The documentation gap is itself the retraction: a commit message is not the record. | `src/regression_zeta.py` |
| **5** | `verify_siegel_parity.py` (a) first asserted the stored `rank1_eig` column at relative tolerance `1e−25`, which would have reported 12 spurious failures at observed deviations of `4e−11 … 1.3e−10`. | **FIXED before the run that produced the committed transcript.** `data/phase5_siegel.csv` stores that column as `mp.nstr(·,10)`, so the binding tolerance is the **file's** ten significant digits; it is now `1e−9`, with the reason in the file. This is Retraction 1 recurring — the lesson did not transfer the first time. **No computed value changed.** | `src/verify_siegel_parity.py` |

## 0.3 Corrections to the prompt that commissioned this extension

Recorded here because they change what the phases can claim, not to score points.

| # | claim in the brief | finding |
|---|---|---|
| P1 | "The paper's Remark 2.4 predicts the window is structurally blind [at γ = 0] because K(L, 0) = 0" | **There is no Remark 2.4 in the paper**, and no remark anywhere in `paper/weil_window_control.tex` about γ = 0 or about real zeros (`grep -n "gamma=0\|real zero\|Siegel\|blind"` returns only a bibliography line). The *mathematical* content is real but is a consequence of equation (K) in Theorem 8, not of a remark: `K(L, 0) = ∫_{−L}^{L} u² sin²(0·u) du = 0`. Phase 5 therefore has no prediction in the paper to confirm or contradict; it establishes the γ = 0 behaviour from scratch. See Phase 5, where the correct statement turns out to be **strictly stronger** than "the δ² coefficient vanishes". |
| P2 | "The odd-sector pole term … is a real zero at δ = ½ **with opposite sign**." | **Same sign, not opposite** — verified entrywise in §5.4(c): `max\|B − Pole\| ≤ 6.2e−61` in all nine cells while `max\|B + Pole\|` is `0.05 … 11.9`, so the test is decisive and there is **no factor of 2 or ½**. What *is* opposite is the odd sector's pole term relative to the **even sector's** (`pole_sign(ODD) = −1`), and the real-pair perturbation reproduces that flip exactly. On the intended reading the remark is right; read literally it is not, and the sign is exactly what §5.3's consequence rests on. |
| P3 | "the two eigenvalues (R+1) ± ‖u‖‖w‖ are attained on `u+w` and `u−w` respectively, [and] `u+w` is palindromic" | **False for every ρ ≠ 1**, refuted exactly over ℚ in all 20 non-degenerate cells (§5.5). `J u = ρ^R w`, so the parity and eigen-directions are the **ρ^R-weighted** `v± = u ± ρ^R w`; `T_R(u+w)` is parallel to `u+w` only when `ρ^{2R} = 1`. The brief's *identity* `T_R = ½[(u+w)(u+w)ᵀ − (u−w)(u−w)ᵀ]` is nevertheless **exact** (25/25 cells). At ρ = 1 both claims hold and `λ₋ = 0`, so the brief states the degenerate case. |
| P4 | Zhu's certificate listed among the even-sector certificates that "carry no information about the Siegel configuration" | **False — and so was my own replacement for it.** The paper was obtained and read on 2026-10-24: Zhu's certificate is over the **full admissible space**, with no parity or reality restriction (Corollary 6.3). My verdict of "domain not established" was right only as a statement about what had been read; the domain was determined all along. **But the consequence I attached to that branch is also wrong**: full-space positivity does *not* constrain a real zero. See §0.4. |

---

## 0.4 Zhu's certificate: domain resolved, and a consequence withdrawn (2026-10-24)

The paper was obtained (`~/Downloads/2608.24827v2.pdf`, **arXiv:2608.24827v2**, dated
September 3, 2026) and §§1.2, 2, 6 were read. This settles the question §5.3 and the note's
§5 had left open, and it **refutes my own speculation about what settling it would mean.**

**What Zhu proves.** Quoted, not paraphrased:

* **Theorem 1.2** (certified positivity, support 1.6): `Q(f) ≥ 8.9 × 10⁻¹⁸ ‖f‖₂²` for every
  real **even** `f` with `supp f ⊆ [−0.8, 0.8]`.
* **Lemma 6.1** (parity splitting): for real `f = f_e + f_o`, `Q(f) = Q(f_e) + Q(f_o)`, where
  `Q(f_e)` carries the pole term `+2(∫ f_e cosh(x/2) dx)²` and
  `Q(f_o) = −2(∫ f_o sinh(x/2) dx)² + (2π)⁻¹∫|F_o|²Ψ_L`; "the pole term changes sign in the
  odd sector". For complex `f`, `Q(f) = Q(Re f) + Q(Im f)`.
* **(14)**: `Q(f) ≥ 8.2065 × 10⁻¹⁵ ‖f‖₂²` for all real **odd** `f`, `supp f ⊆ [−0.8, 0.8]`,
  from `λ_min(M_odd) = 9.1183 × 10⁻¹⁵` and a verified Cholesky at that shift.
* **Theorem 6.2** (simple even ground state, support 1.6):
  `8.9 × 10⁻¹⁸ ≤ λ^even₁ ≤ 2.523 × 10⁻¹⁶`, `λ^even₂ ≥ 2.085 × 10⁻¹²`,
  `8.206 × 10⁻¹⁵ ≤ λ^odd₁ ≤ 2.347 × 10⁻¹⁴` — clearance a factor 32 against the odd sector
  and 8 × 10³ against the second even value.
* **Corollary 6.3** (positivity for arbitrary test functions): for every **complex**
  `f ∈ L²(ℝ)` with `supp f ⊆ [−0.8, 0.8]`, `Q(f) ≥ 8.9 × 10⁻¹⁸ ‖f‖₂²`, "with **no parity or
  reality restriction on f**". Proof: Lemma 6.1, Theorem 1.2 on the even parts, (14) on the
  odd parts.

**So the domain is the full space.** P4's listing of it as even-sector is wrong, and the
note's "this note does not settle it" is superseded.

**And the consequence I predicted for that branch is WITHDRAWN.** §5.3 and the note's §5 both
said that a full-space infimum *would* constrain a real zero, "since positivity on the full
space includes the odd sector where the response is −2S_δ²". **That does not follow.** By the
explicit formula the certified quantity is the sum over **all** zeros, so if a real pair
`½ ± δ` were present, for odd `f`

```
Q(f) = Σ_{remaining zeros} − 2 S_δ(f)²   ≥   8.9 × 10⁻¹⁸ ‖f‖²
```

merely says `Σ_remaining ≥ 8.9 × 10⁻¹⁸‖f‖² + 2S_δ(f)²`. That is satisfiable for every `δ`:
`Σ_remaining` is not independently bounded, and `Σ_remaining = Q(f) + 2S_δ(f)²` is an
identity, so `Q` cannot supply the bound on itself. Zhu says the same thing in his own terms
(§1.1): a proof of `λ*(L) > 0` at a given `L` is "a **finite fragment of RH**", while only a
*negative* `λ*` would disprove RH. **A finite-window positivity certificate carries no
information about whether a real zero exists, whatever its parity domain.**

**Which ingredient is actually missing, stated correctly.** To exclude a real pair one needs
an unconditional **upper** bound on the remaining-zero form, obtained independently of `Q`:
then `Σ_remaining < 8.9×10⁻¹⁸‖f‖² + 2S_δ(f)²` would contradict the certificate. The note's §5
said "lower bound", which is the ingredient for a different inference — attributing an
*observed* negativity to a planted orbit, as the companion DH note does — and not this one.

**What is unaffected.** §5.3's theorem and both corollaries. Even-sector blindness is a
statement about the **response** of the form to a planted real zero, and it is exact in `δ`;
it does not depend on any certificate's domain. What changes is only the scope paragraph that
discussed the certificates.

**Independent confirmation of the mechanism, and the priority position.** Zhu's **Lemma 6.1**
is the `δ = ½` instance of §5.3's theorem: the same factorisation, the same
`+2(cosh)²` / `−2(sinh)²` split, the same sign reversal, derived independently and in print
first (September 3, 2026), and validated numerically there to relative error `2 × 10⁻¹¹`
(`pwl parity recon.py`). §5.3's content beyond it is the generalisation from `δ = ½` — the
**pole** — to an arbitrary `δ ∈ (0, ½)` — a **real zero** — which is what turns the algebra
into a statement about the Siegel configuration. This is the same identification §5.4(c)
already verified from the other side: a pole at `s = 1` and a real zero at `½ + δ` are the
same object at `δ = ½`. No priority is claimed for the parity splitting itself.

**One further note on versions.** v2 **retracts** a claim of an earlier draft at support 2.38
(`L = 1.19`): Remark 3.3 explains the error and Lemma 3.2 shows no such sharpening can exist.
Citations must therefore carry the version.

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
that vanish at ±L. (One substitution is degenerate at `L = 0.8` — χ₋₄ against χ₈ — and that
observation is recorded in §3.5, "Observations not yet explained", with the part of it that
is proved separated from the part that is not.) The smaller (but still decisive) ratios in the `EVEN` basis are because
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

# Phase 3 — Gate B and the baseline

Script: [`src/phase3_gateB.py`](src/phase3_gateB.py). Transcripts:
[`results/phase3_gateB.txt`](results/phase3_gateB.txt) (main sweep),
[`results/phase3_gateB_N20.txt`](results/phase3_gateB_N20.txt),
[`results/phase3_gateB_N24.txt`](results/phase3_gateB_N24.txt) (ζ anchor).
Data: [`data/phase3_gateB.csv`](data/phase3_gateB.csv) (+ `_N20`, `_N24`).

Both bases are complete in `L²_even[−L, L]`:
`EVEN_D = cos((2k+1)πx/2L)` (vanishes at ±L) and `EVEN = cos(kπx/L)` (contains the
constant, jumps at ±L). Every number below is `λ_min(N)`, a **variational upper bound** for
`λ*_χ(L)`; it is never used as a lower bound.

## 3.1 Gate B: **PASS** **(T2)**

* **Monotonicity: 0 violations in 336 rows** (`N = 4 … 16`, four `L`, two bases, six
  functions). `λ_min(N)` is non-increasing in `N` everywhere.
* **Positivity: `λ_min > 0` in every row.** Note this is *expected* under GRH and is
  **not** evidence for it (non-claim 5).
* **Two-basis agreement:** the two independent bases bracket each `λ*_χ(L)` and agree to
  10–15% at `L = 0.8`, with `EVEN` giving the smaller (better) bound for the χ's and
  `EVEN_D` for ζ.
* **ζ control reproduces the parent repository bit-for-bit.** At `L = 0.8`:

  | N | basis | this run | parent `data/stage3_gateB_lambda_min.csv` |
  |---|---|---|---|
  | 16 | `EVEN_D` | 2.541042905e−17 | 2.54104290505e−17 |
  | 20 | `EVEN_D` | 2.30749258e−17 | 2.30749258029e−17 |
  | 24 | `EVEN_D` | 2.270200618e−17 | 2.2702006184e−17 |
  | 16 | `EVEN` | 3.199462217e−17 | 3.19946221698e−17 |
  | 24 | `EVEN` | 2.335713142e−17 | 2.33571314150e−17 |

  Agreement to **all printed digits**, on a different CPython (3.14.7 vs 3.12.2). This is
  the control that a sector or sign bug has not been introduced into the shared code path,
  and it also settles the environment question in `REPRODUCE.md`.

  On Zhu's interval: the convergent `N = 24` value `2.2702e−17` matches Zhu's certified
  **upper** bound `2.27e−17` to three significant figures. The script's `inZhu` flag reads
  `False` because `2.27020e−17 > 2.27000e−17` by `2e−22`; **the parent repository's own
  committed CSV reports `in_zhu_interval = False` at every `N` for the same reason**, so
  this is agreement with the parent, not a discrepancy with Zhu.

## 3.2 The baseline: `λ*_χ(L)` is 4 to 27 orders of magnitude above `λ*_ζ(L)` **(T2)**

`λ_min` at `N = 16`, `EVEN_D`, `dps = 50`:

| L | ζ | χ₋₃ | χ₋₄ | χ₅ | χ₋₇ | χ₈ |
|---|---|---|---|---|---|---|
| 0.5 | 1.071e−6 | 2.275e−2 | 1.481e−1 | 9.558e−2 | 4.599e−1 | 3.717e−1 |
| 0.8 | 2.541e−17 | 1.133e−5 | 6.175e−4 | 8.970e−4 | 6.336e−2 | 3.153e−2 |
| 1.0 | 7.176e−26 | 8.787e−10 | 7.434e−7 | 2.361e−6 | 1.943e−3 | 1.119e−3 |
| 1.3 | 5.908e−35 | 4.689e−19 | 1.970e−14 | 1.256e−12 | 1.227e−7 | 1.361e−7 |

`λ*_χ/λ*_ζ` at the same points: `2e4 … 4e5` at `L = 0.5`, `4e11 … 2e15` at `L = 0.8`,
`1e16 … 3e22` at `L = 1.0`, `8e15 … 2e27` at `L = 1.3`.

> **This is the single largest structural difference, and it is not the conductor term.** It
> is the **missing pole**. On the ζ side the geometric form is
> `2F(i/2)² + Arch − log π − Prime`, and `λ*_ζ(L)` is the residue of a near-exact
> cancellation between a positive rank-one pole term and the rest — which is why it is
> `10⁻¹⁷` at `L = 0.8` rather than `O(1)`. For χ (non-principal) that term is simply **not
> there**, and `log q · I` stands in its place: a multiple of the identity rather than a
> rank-one projector. A multiple of the identity cannot cancel anything, so there is no
> near-cancellation and `λ*_χ(L)` sits where the zeros put it.

Two consequences recorded for later phases:

1. **The parent repository's small-`δ` convention does not transfer.** It measures
   `C = −λ_min(Q_δ)/δ²`, legitimate for ζ because `λ_min(Q_0) ≈ 10⁻¹⁷`. On the χ side
   `λ_min(Q_0)` is `10⁻⁵ … 10⁻¹`, so the baseline must be subtracted (Phase 4 does, and
   reports both).
2. **Detection needs `C(L)δ² > λ*_χ(L)`**, and with `λ*_χ(0.8)` as large as `6e−2` that is a
   real obstruction at the brief's `δ ≤ 3/8` — borne out in Phase 4, where χ₋₇ and χ₈ are
   *not* detected at `L = 0.8` while χ₋₃, χ₋₄, χ₅ are.

## 3.3 The spectral gap, which turns out to matter **(T2)**

`λ₂/λ₁` at `N = 16`, `EVEN_D`:

| L | ζ | χ₋₃ | χ₋₄ | χ₅ | χ₋₇ | χ₈ |
|---|---|---|---|---|---|---|
| 0.5 | 2.0e4 | 65.6 | 10.9 | 20.6 | **4.26** | 6.20 |
| 0.8 | 5.4e5 | 1.4e4 | 1.9e3 | 866 | **17.8** | 58.6 |
| 1.0 | 6.5e5 | 1.8e5 | 7.0e4 | 3.2e4 | 471 | 1.5e3 |
| 1.3 | 5.0e5 | 9.1e5 | 2.1e6 | 7.1e5 | 2.3e5 | 2.1e5 |

The parent note's §B.5 attributes the measured `C(L)` to the *re-optimised* minimiser being
able to travel far from the ground state "because `Q_0` has a near-null space with
`λ₂/λ₁ ∼ 10⁶`". **That near-null space is a ζ feature at `L = 0.8` and a χ feature only from
`L ≈ 1.3`**: at `L = 0.5` the χ gaps are `4 … 66`, four orders of magnitude smaller than
ζ's. By `L = 1.3` all six have converged to `10⁵ … 10⁶` and the distinction disappears. So
the window length at which the χ forms acquire ζ-like mobility is itself conductor-dependent,
and that is where the open observation of §4.5 lives.

## 3.4 Does the conductor enter like height? **(T2 for the residuals, T3 for any coefficient)**

Script: [`src/phase3b_lpred.py`](src/phase3b_lpred.py). Data:
[`data/phase3b_lpred.csv`](data/phase3b_lpred.csv).

Operationalisation — the parent repository's own (`src/stage3e_certify.py`): the window
"resolves" `γ₁` when the **moved-pair** counterfactual at `γ₁` first makes the form
indefinite. `L_onset` is bisected to `±0.005` at fixed `δ = 3/8`, `N = 16`, `EVEN_D`.
The two candidate predictors are

```
L_pred^(1) = ½ log(γ₁/2π)        (height only)
L_pred^(q) = ½ log(q γ₁/2π)      (ANALYTIC CONDUCTOR q γ₁)
```

| | q | γ₁ | `L_onset` | `L_pred^(1)` | `L_pred^(q)` | `onset − L_pred^(1)` | `onset − L_pred^(q)` |
|---|---|---|---|---|---|---|---|
| **ζ** | 1 | 14.1347 | **0.39668** | 0.40538 | 0.40538 | **−0.0087** | **−0.0087** |
| χ₋₃ | 3 | 8.0397 | 0.63926 | 0.12326 | 0.67257 | +0.5160 | −0.0333 |
| χ₋₄ | 4 | 6.0209 | 0.68145 | −0.02132 | 0.67183 | +0.7028 | +0.0096 |
| χ₅ | 5 | 6.6485 | 0.70606 | 0.02825 | 0.83297 | +0.6778 | −0.1269 |
| χ₋₇ | 7 | 4.4757 | 0.86777 | −0.16960 | 0.80335 | +1.0374 | +0.0644 |
| χ₈ | 8 | 4.9000 | 0.83613 | −0.12432 | 0.91540 | +0.9605 | −0.0793 |

**The parameter-free comparison** (the meaningful one — is the residual a constant?), over
the five characters:

| predictor | mean residual | sd | spread |
|---|---|---|---|
| `L_pred^(1)` | **+0.7789** | 0.1923 | 0.5214 |
| `L_pred^(q)` | **−0.0331** | **0.0667** | **0.1913** |

> **`L_onset ≈ ½ log(q γ₁/2π)`, with no free parameter, to within −0.033 ± 0.067 across the
> five characters — and the ζ row, where `q = 1` and the two predictors coincide, has
> residual −0.0087, consistent with the same offset.** So the conductor does enter the
> resolution scale the way height does, through the analytic conductor `qγ₁`, and the test
> holds across all six L-functions including the `q = 1` case.

**Why the two-parameter least-squares fits must not be read as the test.** They are
reported by the script and they *favour the wrong predictor*:

| predictor | α | β | max\|resid\| | rms |
|---|---|---|---|---|
| `L_pred^(1)` | 0.7196 | **−0.8091** | 0.0554 | 0.0280 |
| `L_pred^(q)` | 0.2018 | +0.6986 | 0.1048 | 0.0603 |

The `L_pred^(1)` fit has a **negative slope**: with a free intercept and a free slope, five
points, and a predictor whose values happen to run *opposite* to `L_onset` across this
family (larger `q` ⟹ smaller `γ₁` here), the regression inverts it and then fits the residual
tightly. A predictor that must be given slope −0.81 is not a predictor. The
`L_pred^(q)` fit has slope +0.70 and intercept +0.20, both near the ideal `(1, 0)`.

**What is not claimed.** Five conductors spanning `q = 3…8` cannot separate `½log q` from any
other slowly growing function of `q`, and `λ*_χ(L)` and `γ₁(χ)` move together across this
family, so they are not independently varied. No coefficient is claimed **(T3)**; what is
claimed is the parameter-free residual comparison above **(T2)**.


## 3.5 Observations not yet explained

As in §4.5: reproducible, recorded, and **not** presented as a result. Here the two halves
have different status and are separated deliberately.

**Proved, and therefore not a hypothesis (T1).** At `L = 0.8` the prime side of the window
cannot distinguish χ₋₄ from χ₈: the character-substitution control of §2.4 returns a ratio of
exactly `1.0`. The reason is complete and elementary — `2L = 1.6`, so the window carries only
the prime powers `n < e^{1.6} = 4.953`, namely `n = 2, 3, 4`, and
`χ₋₄ = (0, 1, 0, −1)` and `χ₈ = (0, 1, 0, −1)` **agree on all three**. They first differ at
`n = 5`, which enters once `2L > log 5`, i.e. `L > 0.8047`; at `L = 1.0` the same control
discriminates by `×6.3e4`. **This is explained, exactly, and is filed here only because the
question it raises is not.**

**Not explained (open).** At a window this short the character enters the form *only* through
`q` and the parity `a` — the arithmetic twist contributes nothing that separates these two
characters. That is the sharpest available statement of the confounding that limits §3.4 and
§4.5: across the five real primitive characters of conductor `≤ 8`, the quantities a short
window can actually see are few, and `q`, `a` and `γ₁(χ)` are not independently varied by
this family. Breaking it needs either a larger window (so more prime powers enter and the
twist becomes visible) or characters chosen so that `q` and `γ₁` move oppositely. Neither is
attempted here, and no mechanism hypothesis is offered.

# Phase 4 — the sensitivity law

Scripts: [`src/phase4_sensitivity.py`](src/phase4_sensitivity.py),
[`src/phase4b_sector_extremal.py`](src/phase4b_sector_extremal.py),
[`src/phase4c_decomposition.py`](src/phase4c_decomposition.py).
Data: [`data/phase4_sensitivity.csv`](data/phase4_sensitivity.csv),
[`data/phase4b_sector_extremal.csv`](data/phase4b_sector_extremal.csv),
[`data/phase4c_decomposition.csv`](data/phase4c_decomposition.csv).

## 4.0 What the brief predicted, and what had to change first

`K(L,γ) = ∫_{−L}^{L} u² sin²(γu) du` is the exact value of an extremal problem over real even
`f` supported in `[−L, L]` with `‖f‖₂ = 1`. **That function class does not mention the
L-function**, so `4K(L, γ)` transfers from ζ to `L(s, χ)` unchanged, as the brief says — this
is a tautology, not a measurement, and it is recorded as such **(T1)**.

One thing had to change before `C_χ` could be measured at all. The parent repository defines
`C = −λ_min(Q_δ)/δ²`, which drops the baseline because `λ_min(Q_0) ≈ 10⁻¹⁷` for ζ. Phase 3
shows `λ_min(Q_0)` is `10⁻⁵ … 10⁻¹` on the χ side, so the quantity measured here is the
**baseline-subtracted**

```
C(L) = [ λ_min(Q_0) − λ_min(Q_δ) ] / δ²,
```

with `Q_0` the `δ = 0` member of the *same* family (`Q + 2·pair` for the added-quartet
counterfactual, `Q + pair` for the moved-pair one). The raw parent-convention value is
reported alongside in the CSV (`C_raw`) so the two cannot be confused. On the ζ rows the two
agree to the printed precision, which is the check that the change is a generalisation.

**ζ regression, in the same run:** at `δ = 0.02`, `C(0.8) = 0.310797595` and
`C(1.3) = 2.498698441`. The parent's committed `data/stage3_C_vs_CS.csv` records
`0.3107975951` and `2.498698441`. Bit-for-bit. (`L = 1.0` is not comparable: the parent used
`N = 18`, this run `N = 16`.)

## 4.1 `C_χ/4K(L, γ₁(χ))`, both counterfactuals, kept distinct **(T2)**

Moved-pair convention, `EVEN_D`, `N = 16/16/20`, `dps = 80/90/105`:

| L | δ | ζ | χ₋₃ | χ₋₄ | χ₅ | χ₋₇ | χ₈ |
|---|---|---|---|---|---|---|---|
| 0.8 | 1/4 | 0.673 | **0.0087** | 0.0630 | 0.0223 | **0.550** | 0.348 |
| 0.8 | 3/8 | 0.704 | **0.0119** | 0.0660 | 0.0241 | **0.555** | 0.351 |
| 1.0 | 1/4 | 0.687 | 0.177 | 0.434 | 0.0183 | 0.0963 | 0.0479 |
| 1.0 | 3/8 | 0.716 | 0.196 | 0.751 | 0.236 | 0.101 | 0.0513 |
| 1.3 | 1/4 | 0.935 | 0.728 | 0.675 | 0.688 | 0.501 | 0.736 |
| 1.3 | 3/8 | 0.972 | 0.765 | 0.710 | 0.814 | 0.636 | 0.845 |

The added-quartet convention agrees with the moved-pair one to within 2% in every cell
(full table in the CSV); the two are kept as separate rows throughout, as in §10 of the
paper. At the small `δ = 0.02` used for the law itself the same pattern holds, more sharply:
`0.419` (ζ) against `0.0072, 0.0608, 0.0211, 0.545, 0.345` at `L = 0.8`.

The ratio differs from the ζ case by far more than the noise floor. **That measurement, and
the mechanism proposed for it, are an open observation and are collected in §4.5 —
"Observations not yet explained" — rather than presented as a result of this extension.**
§4.2 is independent of it and is T1.

## 4.2 The extremal is **sector-dependent**, and the two sectors are exactly complementary **(T1)**

Theorem 8 is stated for the even sector, and `K(L,γ)` is the squared norm of
`f ↦ F'(γ)` on the normalised **even** admissible sector. For real **odd** `f`, writing
`F = i F̃` with `F̃(t) = ∫ f(x) sin(tx) dx`, the same extremal problem is

```
K_odd(L, γ) = sup { |F̃'(γ)|² : f real ODD, supp f ⊂ [−L,L], ‖f‖₂ = 1 }
            = ∫_{−L}^{L} u² cos²(γu) du,     attained by f ∝ u cos(γu),
```

and therefore, writing `K_even` for the `K` of Theorem 8,

> **(D4) [T1]** `K_even(L, γ) + K_odd(L, γ) = ∫_{−L}^{L} u² du = 2L³/3` **for every γ**.

Checked at `L = 0.5, 0.8, 1.0, 1.3, 2.0` and `γ = 0, 10⁻⁶, γ₁(χ₋₄), γ₁(χ₋₃), γ₁(ζ), 100`:
the sum equals `2L³/3` to **0.0** (exact cancellation of the closed forms) in 28 of 30 cells
and to `≤ 1.1e−50` in the two `γ = 10⁻⁶` cells; both extremals are attained to `10⁻³⁰`.

Two consequences, both of which matter later:

* The paper's `γ`-uniform bound `K ≤ 2L³/3` **is not a slack bound**: `2L³/3` is the *total*
  over the two sectors, and `K_even` attains it exactly where `K_odd` vanishes.
* `K_even(L, 0) = 0` **forces** `K_odd(L, 0) = 2L³/3`, the `γ`-uniform **maximum**. So at
  `γ = 0` the even sector has zero derivative-sensitivity and the odd sector has the most it
  can ever have. This is Phase 5's subject.

## 4.4 Baseline, not law: where detection happens **(T2)**

`λ_min(Q_δ) < 0`, moved-pair, `EVEN_D`:

| L | δ | ζ | χ₋₃ | χ₋₄ | χ₅ | χ₋₇ | χ₈ |
|---|---|---|---|---|---|---|---|
| 0.8 | 0.02 | yes | no | no | no | no | no |
| 0.8 | 1/4, 3/8 | yes | yes | yes | yes | **no** | **no** |
| 1.0 | 0.02 | yes | yes | yes | yes | no | no |
| 1.0 | 1/4, 3/8 | yes | yes | yes | yes | yes | yes |
| 1.3 | all | yes | yes | yes | yes | yes | yes |

This is purely the Phase-3 baseline: detection needs `C(L)δ² > λ*_χ(L)`, and
`λ*_{χ₋₇}(0.8) = 6.3e−2` against `C δ² = 0.233 · 0.1406 = 3.3e−2`. The two largest
conductors are the two that fail, but **with five conductors that ordering is not evidence
that `q` is the variable** — `λ*_χ(0.8)` and `γ₁(χ)` move together across this family.


## 4.5 Observations not yet explained

Everything in this section is an **open observation**, kept in the record because it is
reproducible and may matter, and kept *out* of the results of this extension because the
mechanism is not established. The tier labels below are the honest ones: the measurements
are T2, the mechanism is a **hypothesis** at T2, and nothing here is T1.

**The observation.** `C_χ/4K(L, γ₁(χ))`, the saturation of the transferred extremal, is
strongly χ-dependent: at `L = 0.8` the five characters span `0.0087 … 0.555`, a factor of
**64**, with ζ at `0.673` sitting inside rather than above that range; by `L = 1.3` all six
have converged to `0.50 … 0.97` and the spread has largely washed out. The discrepancies are
`O(1)` relative against a numerical noise floor of `10⁻⁴⁰`, so they are not a precision
artefact — but "not an artefact" is not the same as "explained".

**What is *not* claimed.** That `q` is the variable. `λ*_χ(L)`, `λ₂(Q_χ)` and `γ₁(χ)` all
move together across these five characters and are not independently varied, so the family is
too small and too correlated to attribute the effect to the conductor, to the first
eigenvalue, to the second, or to the first ordinate.

### 4.5.1 Where the χ-dependence sits **(T2 measurement)**

The parent note's §B.5 attributes `C(L)` to `4(F'(γ)² + F(γ)F''(γ))` evaluated at the
**re-optimised** minimiser of `Q_δ`. Phase 4c measures the three pieces separately at
`L = 0.8`, `δ = 0.02`:

| | γ₁ | `C` | `4(F'²+FF'')` | rel. | `\|F'\|²/K` | `4FF''/C` | `F(γ₁)²` |
|---|---|---|---|---|---|---|---|
| ζ | 14.135 | 0.31080 | 0.31860 | 2.5e−2 | 0.4294 | 3.2e−6 | 2.0e−11 |
| χ₋₃ | 8.040 | 0.0044231 | 0.0044297 | 1.5e−3 | 0.00722 | 1.3e−3 | 1.5e−9 |
| χ₋₄ | 6.021 | 0.046759 | 0.046776 | 3.5e−4 | 0.0606 | 4.3e−3 | 3.9e−7 |
| χ₅ | 6.648 | 0.018368 | 0.018379 | 6.1e−4 | 0.0206 | 2.0e−2 | 2.1e−6 |
| χ₋₇ | 4.476 | 0.229153 | 0.229170 | 7.6e−5 | 0.4802 | 1.2e−1 | 2.5e−3 |
| χ₈ | 4.900 | 0.147947 | 0.147963 | 1.1e−4 | 0.3142 | 8.9e−2 | 7.4e−4 |

Three readings:

1. **The parent's identification of the mechanism transfers exactly.** `C` equals
   `4(F'² + FF'')` at the perturbed minimiser to between `7.6e−5` and `2.5e−2` relative, for
   every character — and the worst agreement is ζ's own.
2. **The cross term is not the explanation.** `4FF''/C ≤ 0.12` everywhere and `≤ 4.3e−3` for
   ζ, χ₋₃, χ₋₄. So the whole χ-dependence is in `|F'(γ₁)|²/K`, the saturation of the
   derivative extremal itself — which reproduces `C/4K` to the cross term, as it must.
3. **The annihilation residual `F(γ₁)²` orders the characters the same way:** `2e−11` (ζ),
   `1.5e−9` (χ₋₃) — the minimiser essentially annihilates `γ₁` — against `2.5e−3` (χ₋₇),
   where it does not. And `|F'(γ₁)|²/K` is **monotone decreasing in `γ₁`** across all five
   characters (`0.480` at `γ₁ = 4.48` down to `0.0072` at `γ₁ = 8.04`), with ζ at
   `γ₁ = 14.13` the single outlier at `0.429`.

The last point identifies the two regimes. For χ, `λ*_χ(0.8)` is large enough that the
minimiser barely moves, so `C_χ` is essentially a property of the **unperturbed ground state
of `Q_χ` and of `γ₁(χ)`** — and since that ground state's transform decays like `t⁻²`,
`|F'(γ₁)|²` falls with `γ₁`. For ζ the minimiser is free to re-optimise inside a near-null
space (Phase 3.3) and reaches a `|F'|²` that has nothing to do with the ground state. The
χ forms acquire that freedom only around `L ≈ 1.3`, which is exactly where §4.5.1's spread
collapses.

### 4.5.2 Proposed mechanism — a HYPOTHESIS, validated on 18 cells **(T2)**

Script: [`src/phase4d_crossover.py`](src/phase4d_crossover.py) (pure arithmetic on the
committed CSVs). Phase 4c also records `|F'(γ₁)|²/K` at the **unperturbed ground state**, so
"did the minimiser move?" is measured rather than inferred:

| L | | `\|F'\|²/K` at the perturbed minimiser | at the ground state | ratio | moved? |
|---|---|---|---|---|---|
| 0.8 | ζ | 0.429436 | **0.000164616** | **×2609** | **YES** |
| 0.8 | χ₋₃ … χ₈ | 0.0072 … 0.480 | identical to ≤ 0.22% | ×1.000 … ×1.002 | no |
| 1.0 | ζ | 0.637101 | 4.20e−5 | ×15164 | **YES** |
| 1.0 | χ₋₃ | 0.0096227 | 0.0015705 | ×6.13 | **YES** |
| 1.0 | χ₋₄, χ₅, χ₋₇, χ₈ | — | ≤ 2.2% apart | ×1.00 … ×1.02 | no |
| 1.3 | ζ, χ₋₃, χ₋₄, χ₅ | 0.129 … 0.847 | ×106 … ×141141 | — | **YES** |
| 1.3 | χ₋₇, χ₈ | 0.0149, 0.0217 | ≤ 7.7% apart | ×1.06, ×1.08 | no |

So it is not a different *constant* but a different **regime**, and the characters cross over
one at a time as `L` grows: ζ already by `L = 0.8`, χ₋₃ by `1.0`, χ₋₄ and χ₅ by `1.3`,
χ₋₇ and χ₈ not yet at `1.3`.

**The criterion.** The minimiser should rotate out of the ground state exactly when the gain
available to the `δ²` term exceeds the price of the first direction it must pay for:

> **(H1) [T2, HYPOTHESIS]** `R(L, δ) = 4δ² K(L, γ₁) / λ₂(Q)`, and the minimiser
> re-optimises **iff `R > 1`**.  This is a criterion validated on 18 cells, **not a
> theorem**; no perturbation-theoretic proof is offered, and what replaces `λ₂` when the
> near-null space has dimension `> 1` is not known (Open question 3).

Checked on all **18** cells of Phase 4c (`L = 0.8, 1.0, 1.3` × six functions, `δ = 0.02`),
with `λ₂` taken from `data/phase3_gateB.csv` and `K` from `data/phase4c_decomposition.csv`:

```
disagreements : 0 / 18
smallest R among the MOVED cells  : 3.715
largest  R among the frozen cells : 0.04821
the two classes are separated by a factor 77.06, with R = 1 strictly inside the gap
```

`R` ranges over 25 orders of magnitude across the table (`9.3e−5` for χ₈ at `L = 0.8` up to
`4.2e+25` for ζ at `L = 1.3`) and the moved/frozen flag follows it without exception. This is
T2 — a criterion validated on 18 cells, not a theorem — but it identifies the variable: the
χ-dependence of `C/4K` is the **second** eigenvalue of `Q_χ`, not the first.

> **Consequence for the parent's §B.5.** Its statement that `C(L)` is a property of the
> re-optimised minimiser "because `Q_0` has a near-null space with `λ₂/λ₁ ∼ 10⁶`" is correct
> for ζ at `L = 0.8` but is *not* a property of the window — it is a property of `λ₂` in
> absolute terms, and on the χ side at the same `L` the minimiser does not move at all, so
> `C_χ` is read off the unperturbed ground state instead.

# Phase 5 — the Siegel corner (γ = 0)

Scripts: [`src/phase5_siegel.py`](src/phase5_siegel.py) (grid, both sectors),
[`src/phase5_onset.py`](src/phase5_onset.py) (odd-sector onset in `L`),
blocks in [`src/planting.py`](src/planting.py), regressions in
[`src/selftest_planting.py`](src/selftest_planting.py).
Data: [`data/phase5_siegel.csv`](data/phase5_siegel.csv),
[`data/phase5_onset_odd.csv`](data/phase5_onset_odd.csv).

**This is a statement about the detector, not about `L(s, χ)`** (non-claim 2). Nothing below
bears on whether any Siegel zero exists.

## 5.1 What the configuration is **(T1)**

A real zero `β = ½ + δ` forces, by the functional equation, the partner `1 − β = ½ − δ`. In
the variable `z` with `ρ = ½ + iz`, that pair is `z = ∓iδ`: **purely imaginary ordinates**,
which is exactly `γ = 0` in the quartet parametrisation. Its contribution to the window form
is `2h(iδ)` with `h(z) = F(z)F(−z)`, and `F(iδ) = ∫ f(x) e^{−δx} dx`, so

```
f even :  F(−iδ) = +F(iδ)   ⟹   contribution = + 2 A_δ(f)²,   A_δ = ∫ f(u) cosh(δu) du
f odd  :  F(−iδ) = −F(iδ)   ⟹   contribution = − 2 S_δ(f)²,   S_δ = ∫ f(u) sinh(δu) du
```

so the perturbation is **rank one with a definite sign that depends on the parity of the test
function**: positive semidefinite on the even sector, negative semidefinite on the odd one.

This is the **same mechanism as the ζ pole term** (`z = ∓i/2`), and the identity is checked,
not asserted: `real_block(L, ·, 1/2, sector)` reproduces the parent repository's pole-term
matrix `2·pole_sign(sector)·P_j P_k/‖·‖` to `1.1e−50` **in all three sectors**
(`selftest_planting.py`, check C), including `pole_sign(ODD) = −1`. Since the parent's pole
term is validated end-to-end by its own ζ Gate A, the Phase-5 block — sector sign and all —
is anchored to already-checked code. Two further degeneration identities hold exactly:
`quartet(γ, 0) = 2·pair(γ)` and `quartet(0, δ) = 2·real(δ)`, in every sector.

## 5.2 The even sector cannot see a Siegel zero. At all. **(T1)**

> *The underlying theorem and its derivation are stated once, in §5.3; its three independent
> verifications are in §5.4 and the function-field mirror in §5.5. This section states the
> consequence for the even sector and reports the 180-row numerical confirmation.*


> **(D5) [T1]** Let χ be real primitive mod q and suppose every non-trivial zero of
> `L(s, χ)` lies on the critical line except a single real pair `½ ± δ`, `δ ∈ (0, ½)`. Then
> for **every** `L > 0` and **every** real **even** `f` supported in `[−L, L]`,
> ```
> Q_χ(f) = Σ_{on-line ρ} |F(γ_ρ)|²  +  2 A_δ(f)²  ≥  2 A_δ(f)²  ≥  0.
> ```
> Hence **no even-sector finite-window Weil-positivity computation, at any window length and
> in any basis, can exclude a Siegel zero.** Equivalently, as a perturbation: the response is
> `+2aaᵀ` with `a_k = ∫φ_k cosh(δu)du/‖φ_k‖`, a rank-one **positive** semidefinite matrix, so
> `λ_min(Q + 2aaᵀ) ≥ λ_min(Q)` for every basis and every `N`.

Note what this is **not**. It is not "the `δ²` coefficient vanishes". The brief's expectation
was that the window is blind because `K(L, 0) = 0`, which is true but is only the statement
that `F'(0) = 0` for even `f` — a vanishing *leading* coefficient. (D5) is stronger and exact
in `δ`: the response is a perfect square of the **wrong sign** to all orders. The window is
not blind at `γ = 0` in the even sector; it is **anti-detecting**, and a planted real zero
makes the form *more* positive.

**Measured, as a check on the theorem rather than as evidence for it: 180 even-sector rows**
(`L = 0.8, 1.4` × `N = 8, 12, 16` × six functions × `δ = 0.1, 0.25, 0.375, 0.45, 0.49`),
**180 with response `≥ 0`, 0 detections.** Zero violations.

> **Scope note for the parent paper, not a correction to it.** The paper's §6.1 / Theorem 11
> certificate, and its Theorem 8, are explicitly even-sector statements (`tex:268`: "all
> statements here are for the even sector"). (D5) says that sector is precisely the one that
> cannot respond to a real zero. Nothing in the paper claims otherwise; what is new is that
> the exclusion is **structural and total**, not a matter of window length.

### A column that must not be misread

`phase5_siegel.py` reports `Theta = 2aᵀQ⁻¹a` and `theta_gt_1` in **both** sectors.
`Θ > 1` is the Schur-complement criterion for `Q − 2aaᵀ` to be indefinite, i.e. the
**odd-sector** detection criterion. In the **even** sector the perturbation is `+2aaᵀ` and
`Θ` is *not* a detection criterion there — it measures only the alignment of the rank-one
direction with the inverse form, and the even-sector rows show `Θ ≫ 1` together with
`detected = no`, which is consistent, not contradictory.

## 5.3 Sign structure **(T1)**

> **Theorem (sign structure).** Let `f` be real, supported in `[−L, L]`, and let a real zero
> pair `{½ + δ, ½ − δ}`, `δ ∈ (0, ½]`, be present in the zero multiset. Its contribution to
> the finite-window Weil form is
>
> ```
>   + 2 A_δ(f)²     if f is even,     A_δ(f) = ∫ f(u) cosh(δu) du
>   − 2 S_δ(f)²     if f is odd,      S_δ(f) = ∫ f(u) sinh(δu) du
> ```
>
> for **every** `L > 0` and in **every** basis. In particular the contribution is rank one,
> positive semidefinite on the even sector and negative semidefinite on the odd sector.

**Derivation.** The zero side of the explicit formula is `Σ_ρ h(z_ρ)` with
`h(z) = F(z)F(−z)`, `ρ = ½ + iz`; the autocorrelation **factorises**, which is the whole
content. A real zero `β = ½ + δ` forces `1 − β` by the functional equation, so the pair sits
at `z = ∓iδ` and contributes `h(iδ) + h(−iδ) = 2h(iδ) = 2 F(iδ) F(−iδ)`. Now
`F(iδ) = ∫ f(x) e^{−δx} dx`, so splitting `e^{−δx} = cosh(δx) − sinh(δx)` and using the
parity of `f`:

```
f even :  F(iδ) = A_δ,   F(−iδ) = + A_δ    ⟹   2 F(iδ)F(−iδ) = + 2 A_δ²
f odd  :  F(iδ) = − S_δ, F(−iδ) = + S_δ    ⟹   2 F(iδ)F(−iδ) = − 2 S_δ²
```

— even `f` gives `A·A`, odd `f` gives `S·(−S)`. The sign is not a convention: it is the
parity of `F` at a purely imaginary argument. ∎

**Consequence, for the even sector.** Suppose every non-trivial zero of `L(s, χ)` lies on the
critical line except one real pair `½ ± δ`. Then for every real even `f` supported in
`[−L, L]`,

```
Q_χ(f) = Σ_{on-line ρ} |F(γ_ρ)|² + 2 A_δ(f)² ≥ 0.
```

**Even-sector finite-window Weil positivity is compatible with a real zero at every
`δ ∈ (0, ½]`.** An even-sector positivity certificate therefore carries no information about
the Siegel configuration — at any window length, in any basis, to any precision. This
applies to **Theorems 7, 8 and 11 of this repository's paper**, all three of which are
labelled even-sector in their statements. It is not a defect in them — none claims otherwise
— but it does bound what they can be used for.

*The numbering was checked rather than assumed*, because this repository has a recorded
history of it drifting (an interim draft once shifted Theorem 8 → 9). Of the eleven
environments sharing the `theorem` counter in `paper/weil_window_control.tex`: `thm:zeta` is
**#7** ("relative perturbation estimate; **even sector**"), `thm:zetasharp` is **#8**
("sharper, conditional; **even sector**"), `cor:uniform` is #9, and `thm:certwitness` is
**#11** ("certified four-mode witness", whose witness is the four quarter-wave `EVEN_D`
modes `cos((2k+1)πx/2L)/√L`). `tex:268` states globally: "all statements here are for the
even sector".

**What is *not* settled: Zhu's certificate.** The brief that commissioned this step listed
Zhu's `λ*(0.8) ≥ 8.9×10⁻¹⁸` alongside those three. That does **not** follow, and this
extension does not establish it either way, because **the domain of Zhu's infimum is left
open by this repository on purpose**: `paper/weil_window_control.tex:312–316` says only that
"a certificate proved on the full admissible space may be quoted here unchanged", precisely
so that the question need not be settled. The two cases are different in kind:

* if Zhu's `λ*(L)` is an infimum over the **even sector only**, then by the theorem above it
  carries no information about a real zero, exactly as for Theorems 7/8/11;
* if it is an infimum over the **full admissible space**, then it *does* constrain one:
  positivity on the full space includes the odd sector, where the response is `−2S_δ(f)²`, so
  `λ*(L) > 0` would force `2 S_δ(f)² ≤ Σ_{on-line}|F̃(γ_ρ)|² − λ*(L)‖f‖²` for every odd `f`
  — an unconditional constraint on real zeros of ζ in the certified range.

Which case holds is a property of Zhu's paper, which was **not obtained or read** for this
extension. It is recorded as Open question 2 and is the single most consequential loose end
here. **No claim is made in either direction.**

**Consequence, for the odd sector.** By (D4) of §4.2, `K_even(L,γ) + K_odd(L,γ) = 2L³/3`
identically in `γ`, and `K_even(L,0) = 0`. Hence `K_odd(L, 0) = 2L³/3`: the `γ = 0` corner is
the point of **maximum** derivative-evaluation sensitivity available at any ordinate, with
constant `2L³/3`. The blindness at `γ = 0` is a property of one sector, not of the window.


## 5.4 Verification of the sign structure **(T2, three independent checks)**

Script: [`src/verify_siegel_parity.py`](src/verify_siegel_parity.py) →
[`results/verify_siegel_parity.txt`](results/verify_siegel_parity.txt),
[`data/verify_siegel_parity.csv`](data/verify_siegel_parity.csv). `dps = 60`; no float64 on
any path. **All three pass.**

### (a) The planted object was the real **pair**, not a degenerate quartet

Decided **from the committed Phase-5 data**, not from a reading of the code. First the
identity `quartet_block(γ=0, δ) = 2 · real_block(δ)` holds to `5.0e−60` (3 sectors × 2 `L` ×
3 `δ`), so the two differ by exactly the factor 2 and the question is well posed. Then 12
stored Phase-5 rows had their `lam_planted` recomputed both ways:

| | reproduced by PAIR | reproduced by QUARTET |
|---|---|---|
| 12 of 12 probed rows | rel. `1.4e−14 … 1.7e−12` | rel. `1.1e−3 … 1.5e−1` |

**Verdict: PAIR, unanimously.** The discrimination is real, not degenerate: for χ₋₃, χ₈ and
for the odd-sector ζ rows the two conventions differ by 0.1% to 15%. (The two agree to
`1e−20` only in the even-sector ζ rows, where `λ_min` is set by something other than the
planted term.) The stored `rank1_eig` column also matches `±2‖A‖²` — the pair value — to the
`10⁻⁹` that the CSV's 10 stored digits permit, rather than `±4‖A‖²`.

> **No factor-of-2 correction is needed and no retraction follows from (a).** Phase 5 planted
> `{½ + δ, ½ − δ}` throughout.

### (b) Rank one, correct sign, eigenvalue `= +2A²` / `= −2S²` — 150 of 150 cells

5 characters × `L ∈ {0.8, 1.0, 1.4}` × 2 sectors × `δ ∈ {0.1, 0.25, 0.375, 0.45, 0.49}`,
`N = 10`. `A` and `S` are computed by **`mp.quad` on the explicit normalised eigenfunction**
`f = Σ v_k φ_k/‖φ_k‖`, so this is an independent check of the closed forms in
`planting.F_at_imag`, not a restatement of them.

| quantity | worst over the 150 cells |
|---|---|
| relative error, eigenvalue vs `+2A²` (even) / `−2S²` (odd) | **8.5e−61** |
| rank-one residual `max_ij \|B_ij − λ v_i v_j\|` | **2.5e−60** |
| largest non-dominant eigenvalue (rank-one-ness) | **1.5e−60** |
| the complementary integral (`S` on even `f`, `A` on odd `f`) | **exactly 0** |
| even-sector eigenvalues | `+3.14 … +6.39`, all `≥ 0` (75/75) |
| odd-sector eigenvalues | `−0.906 … −0.00644`, all `≤ 0` (75/75) |

So the perturbation is rank one to `10⁻⁶⁰`, its single eigenvalue equals the predicted
`±2A²`/`∓2S²` to `10⁻⁶¹`, and the sign is as stated in every cell.

**The matrices do not depend on χ at all.** Grouping the 150 cells by `(L, sector, δ)` leaves
exactly **one** distinct eigenvalue per group: the response to a real zero is a property of
the basis and `δ` only. That is not a side remark — it is why the conclusion is uniform over
the family rather than five separate measurements.

`K_even(L,γ) + K_odd(L,γ) = 2L³/3` was re-checked with `K_odd` obtained **by quadrature**
rather than by subtraction, at `γ ∈ {0, γ₁(χ), 14.1347}` for 5 characters × 3 `L`: worst
deviation **1.6e−61**.

### (c) Pole consistency — exact, and with the sign the data chose

The `δ = ½` real-pair perturbation was compared entrywise with the parent repository's pole
term `2·pole_sign(sector)·P_j P_k/‖·‖`, reporting **both** differences so the data fixes the
sign:

| sector | `pole_sign` | `max\|B − Pole\|` | `max\|B + Pole\|` | vanishes |
|---|---|---|---|---|
| `EVEN_D` | `+1` | `1.9e−62 … 6.2e−61` | `3.28 … 9.13` | `B − Pole` |
| `EVEN` | `+1` | `1.1e−61 … 3.6e−61` | `4.08 … 11.95` | `B − Pole` |
| `ODD` | `−1` | `1.2e−63 … 1.9e−62` | `0.051 … 0.941` | `B − Pole` |

(3 sectors × `L ∈ {0.5, 0.8, 1.3}`.) **`B = Pole` exactly, in all nine cells**, and the
alternative is wrong by `O(1)`–`O(10)`, so the test is decisive. **There is no normalisation
mismatch — no factor of 2 and no factor of ½.**

**A correction to the brief's wording.** The brief described the pole term as "a real zero at
`δ = ½` with opposite sign". The real-pair perturbation at `δ = ½` equals the pole term with
the **same** sign; what is opposite is the *odd* sector's pole term relative to the *even*
sector's, since `pole_sign(ODD) = −1`, and the real-pair perturbation reproduces that flip
exactly. Read that way the brief's remark is right; read literally it is not, and the
distinction matters because §5.3's consequence rests on this sign.

**The χ-analogue is identically zero, and that is the point.** `M_χ` reconstructs as
`Arch − LogPi + LogQ − Prime` with residual **exactly `0.0`**: for non-principal χ there is
no pole term at all. Meanwhile the real-pair perturbation at `δ = ½` is *not* zero for χ
(max entry `2.67`). On the ζ side a pole-like object at `z = ∓i/2` is already present in the
form; on the χ side nothing of the kind is there to cancel against a real zero. `real_block`
is χ-independent, so the pole-consistency check is intrinsically a ζ-side check — which is
exactly why it anchors the construction: the parent's pole term is validated end-to-end by
the parent's own ζ Gate A.

## 5.5 Function-field cross-check **(T1, exact over ℚ)**

Script: [`src/ff_parity_split.py`](src/ff_parity_split.py); note:
[`notes/ff_parity.md`](notes/ff_parity.md). All arithmetic in `fractions.Fraction`; **no
floating point**.

The same configuration — a real pair in a finite window — exists in the parent repository's
function-field control, where RH is Weil's theorem. There it is the rank-two Toeplitz matrix
`T_R = u wᵀ + w uᵀ`, `u_i = ρ^i`, `w_i = ρ^{−i}`, `0 ≤ i ≤ R`. With `J` the flip
`x_i ↦ x_{R−i}`, and because `J u = ρ^R w` and `J w = ρ^{−R} u`, the parity vectors are the
**`ρ^R`-weighted** ones `v± = u ± ρ^R w`, and

| | verified exactly, 25 cells |
|---|---|
| **(FF1)** `T_R = (1/2ρ^R)[v₊v₊ᵀ − v₋v₋ᵀ]` | PASS |
| **(FF2)** `T_R v± = λ± v±`, `λ± = (R+1) ± ρ^{−R}‖u‖²` | PASS |
| **(FF3)** `‖w‖² = ρ^{−2R}‖u‖²` identically, so `ρ^{−R}‖u‖² = ‖u‖‖w‖ ∈ ℚ` | PASS |
| **(FF4)** `λ₊ > 0`, `λ₋ ≤ 0`, equality iff `ρ = 1` | PASS |
| `v₊` palindromic, `v₋` anti-palindromic under `i ↦ R−i` | PASS |

(`ρ ∈ {2, 3/2, 5/4, 7/3, 1}` × `R ∈ {1, 2, 4, 8, 16}`.) So the **`J`-even direction carries
the positive rank-one piece and the `J`-odd direction the negative one** — the same signed
split, with the same attachment of sign to parity, as `+2A²` / `−2S²`. The dictionary is
tabulated in `notes/ff_parity.md`.

> This shows the sign structure is intrinsic to "a real pair in a finite window" and not an
> artefact of the archimedean bookkeeping of `L(s, χ)`. It does **not** transfer the Phase-5
> consequence: in the function field the analogue of RH is a theorem, so there is no
> configuration to exclude. The function field is a control, exactly as in Part I of the
> parent paper.

**Two corrections to the brief, both checked rather than argued.** The brief's identity
`T_R = ½[(u+w)(u+w)ᵀ − (u−w)(u−w)ᵀ]` is **exact** (PASS, 25/25) — it is polarisation, and the
`ρ^R` weights cancel in it. But the brief's further claims, that `(R+1) ± ‖u‖‖w‖` are
"attained on `u+w` and `u−w` respectively" and that `u+w` is palindromic, are **false for
every `ρ ≠ 1`**, refuted in all 20 non-degenerate cells: `T_R(u+w)` is parallel to `u+w` only
when `‖u‖² = ‖w‖²`, i.e. `ρ^{2R} = 1`. The missing factor is exactly `ρ^R`. At `ρ = 1` both
claims become true and `λ₋ = 0`, consistent with the parent's Corollary A3, where strictness
fails precisely at `ρ = 1` — so the brief's version is the degenerate-case statement.
Neither correction touches any parent result: `notes/proofs.md` §A.2 derives the eigenvalues
basis-independently and never names `u ± w` as eigenvectors.

## 5.6 The odd sector is maximally sensitive at γ = 0 **(T1 + T2)**

For real **odd** `f` the same identity reads

```
Q_χ(f) = Σ_{on-line ρ} |F̃(γ_ρ)|²  −  2 S_δ(f)²,
```

which **can be negative**. And by (D4) of Phase 4, the derivative-evaluation extremal in the
odd sector at `γ = 0` is `K_odd(L, 0) = 2L³/3` — the `γ`-uniform **maximum**. So the `γ = 0`
corner is not a blind spot of the window; it is a blind spot of **one sector**, and the
complementary sector is as sensitive there as it ever gets at any ordinate.

`Q − 2aaᵀ` is indefinite iff `Θ(L, δ) = 2aᵀQ⁻¹a > 1`. Measured in the odd sector,
`N = 16`, with an explicit LDL witness (a single trial function with `Q < 0`) recorded for
every detection:

| L = 0.8, `Θ` | δ=0.1 | δ=0.25 | δ=0.375 | δ=0.45 | δ=0.49 |
|---|---|---|---|---|---|
| ζ | 4.4e10 | 2.8e11 | 6.3e11 | 9.0e11 | 1.1e12 |
| χ₋₃ | **2.13** | 13.4 | 30.4 | 43.9 | 52.3 |
| χ₋₄ | 0.064 | 0.401 | 0.910 | **1.32** | **1.57** |
| χ₅ | 0.080 | 0.501 | **1.14** | **1.65** | **1.96** |
| χ₋₇ | 0.0065 | 0.041 | 0.093 | 0.134 | 0.160 |
| χ₈ | 0.0082 | 0.051 | 0.116 | 0.169 | 0.201 |

| L = 1.4, `Θ` | δ=0.1 | δ=0.49 |
|---|---|---|
| ζ | 2.7e39 | 6.5e40 |
| χ₋₃ | 7.8e31 | 1.9e33 |
| χ₋₄ | 4.4e23 | 1.1e25 |
| χ₅ | 3.5e9 | 8.6e10 |
| χ₋₇ | 1.1e4 | 2.8e5 |
| χ₈ | 1.1e4 | 2.8e5 |

* **`Θ` is monotone non-decreasing in `N`** in every cell (`N = 8, 12, 16`), as nested bases
  require — e.g. χ₋₃ at `δ = 0.1`: `1.952 → 2.082 → 2.135`.
* **`Θ ∝ δ²` to 2%** over the whole `δ` range at fixed `(L, N)`: χ₋₃ at `L = 0.8, N = 16`
  gives `Θ/δ² = 213.5, 214.4, 215.8, 216.9, 217.6`. So the odd-sector sensitivity has the
  same `δ²` law as Theorem 8, with `K_odd(L, 0) = 2L³/3` in place of `K_even`.
* **180 of 180 odd-sector rows have response `< 0`; 135 are detections.** Witness supports
  are small (1–8 basis functions), so each is a low-dimensional explicit certificate.

## 5.7 Onset in `L`, odd sector **(T2)**

Smallest scanned `L` at which `λ_min(Q^odd − 2aaᵀ) < 0`, `N = 12`,
`L ∈ {0.5, 0.8, 1.1, 1.4, 1.7, 2.0, 2.3, 2.6, 3.0}`:

| | δ=0.1 | δ=0.25 | δ=0.375 | δ=0.45 | δ=0.49 |
|---|---|---|---|---|---|
| ζ | **0.5** | 0.5 | 0.5 | 0.5 | 0.5 |
| χ₋₃ | 0.8 | 0.8 | 0.8 | 0.8 | 0.8 |
| χ₋₄ | 1.1 | 1.1 | 1.1 | 0.8 | 0.8 |
| χ₅ | 1.1 | 1.1 | 0.8 | 0.8 | 0.8 |
| χ₋₇ | 1.4 | 1.1 | 1.1 | 1.1 | 1.1 |
| χ₈ | 1.4 | 1.1 | 1.1 | 1.1 | 1.1 |

So the answer to the brief's question — *at what `L` and what basis dimension does a real
zero at each `δ` first produce a negative Rayleigh quotient, if ever* — is:

* **even sector: never**, at any `L`, any `N`, any `δ` (D5, a theorem, confirmed on 180 rows);
* **odd sector: by `L = 1.4` with a 12-dimensional basis, for every character and every
  `δ ≥ 0.1`**, and as early as `L = 0.5` for ζ and `L = 0.8` for χ₋₃. By `L = 1.7` the
  margins are `Θ ≳ 10⁹` and by `L = 3.0`, `Θ ≳ 10²⁴`.

The ordering of the onsets is the ordering of `λ*_odd(χ, L)`, i.e. the Phase-3 baseline
again, and the conductor-vs-height caveat of §3.4 applies unchanged: `q` and `γ₁` are not
independently varied across these five characters.

### The honest limitation, stated plainly

The measured `λ_min(Q^odd_χ − 2aaᵀ) < 0` says the planted response `2S_δ(f)²` exceeds the
**computed** on-line sum for some odd `f`. Turning that into an *exclusion* of a real zero
would need an unconditional lower bound on the remaining-zero form `Σ_{on-line}|F̃(γ_ρ)|²` —
which is the same hypothesis that Theorem B′ of the parent needs, and **no Zhu-style
certified lower bound exists for any Dirichlet window** (non-claim 8). So:

* the even-sector statement (D5) is **unconditional and complete** — it closes that route;
* the odd-sector statement is **a detector-sensitivity measurement**, T2, and is *not* an
  exclusion of anything.

The ζ row deserves a specific caveat: `λ*_odd(ζ, 0.8) = 4.98e−14` (which reproduces the
parent's `data/stage3c_even_odd.csv` value `4.98076645922e-14` exactly), so ζ's odd-sector
form is already within `10⁻¹⁴` of singular and *any* negative rank-one perturbation detects.
The ζ detections are therefore nearly vacuous. The informative ones are χ₋₄ and χ₅ at
`L = 0.8`, which clear a baseline of `9.9e−2` and `6.8e−2` respectively with `Θ` only
`1.3 … 2.0` — a genuine `O(1)` margin, not a near-singularity artefact.

# Phase 6 — report

## 6.0 The ζ regression suite: **18/18 PASS, exit 0**

[`src/regression_zeta.py`](src/regression_zeta.py) →
[`results/regression_zeta.txt`](results/regression_zeta.txt). Every asserted value is the
**parent repository's committed one**, recomputed through the shared code path — not a value
recorded by this extension:

| | checks | worst relative deviation | binding precision |
|---|---|---|---|
| R1 | ζ Gate B, `L = 0.8`, `N = 16/20/24`, **both** even bases | `1.7e−12` | the parent stores `mp.nstr(·,12)` |
| R2 | the even/**odd** split at `L = 0.5, 0.8, 1.0, log(13)/2`, `N = 12` — the odd rows are the sector-bug guard | `2.5e−12` | " |
| R3 | `C(0.8) = 0.3107975951`, `C(1.3) = 2.498698441` | `8.7e−11` | the parent stores `mp.nstr(·,10)` |
| R4 | `arch_closed_a(·,·,0)` **≡** parent `arch_closed` | **exactly 0.0** | — |
| R5 | `real_block(δ=1/2)` **≡** parent pole term, all three sectors | `5.3e−51` | — |

Every deviation in R1–R3 sits at the parent's *printing* precision, so the agreement is as
exact as the reference allows; R4 is bit-identical and R5 is at the arithmetic floor. This is
the brief's working rule — keep known-good ζ values in the suite so a sector or sign bug in
the χ code is caught by the ζ control — and it is the reason the χ-side numbers above can be
read as measurements rather than as output of untested code.

## 6.1 Results, by phase, with tiers

| # | Result | Tier | Script | Data |
|---|---|---|---|---|
| **1** | Window form (D1) for real primitive χ: **no pole**, conductor `+log(q/π)` with coefficient **1**, digamma `1/4 → 1/4 + a/2`, primes twisted by `χ(n)` and so of both signs | **T1** | `src/dirichlet_window.py` | — |
| 1 | (D2) the archimedean generalisation is exactly two substitutions; `arch_closed_a(a=0)` **≡** the parent's `arch_closed`, max abs diff **0.0** | **T1** | `src/selftest.py` | `results/phase1_selftest.txt` |
| 1 | `τ(χ) = i^a√q`, hence `ε(χ) = +1`, all five; `Λ(½+it,χ)` real | **T1** | `src/selftest.py` | " |
| 1 | `L(1, χ_d)` = class-number formula to `1.9e−41`; `L(3,χ₋₄) = π³/32`; `L(0,χ_d) = 2h/w` exactly | **T1** | `src/selftest.py` | " |
| 1 | `real_block(δ=1/2)` **≡** the parent's pole-term matrix, `1.1e−50`, **all three sectors** incl. `pole_sign(ODD) = −1` | **T1** | `src/selftest_planting.py` | `results/phase45_selftest_planting.txt` |
| **2** | 966 critical-line ordinates to `T = 275`, `dps_root = 60`; completeness band `≈ 1` against `θ(T)/π`, **no missed or spurious zero** | **T2** | `src/zeros.py`, `src/zeros_summary.py` | `data/zeros_<χ>.txt` |
| 2 | 649 ordinates agree with **LMFDB** to `9.0e−28` (LMFDB's own truncation); characters identified by Euler factors at 19 primes | **T2, cross-checked** | `src/lmfdb_check.py` | `results/phase2_lmfdb.json` |
| 2 | χ₈'s 213 ordinates reproduced **bit-identically** by an independent re-invocation | **T2** | `src/zeros.py` | `data/zeros_chi_8.txt` |
| **2** | **GATE A PASSES**: 5 χ × 3 sectors × 3 L; exponents −2.68…−3.02 (`EVEN_D`), −2.66…−2.96 (`ODD`), −0.75…−0.88 (`EVEN`, predicted −1) | **T2** | `src/phase2_gateA.py` | `data/phase2_gateA.csv` |
| 2 | Gate A per-basis decay prediction (−3 for bases vanishing at ±L, −1 for the jumping basis) | **T1** | `notes/dirichlet_window.md` §6 | — |
| 2 | Four controls discriminate: conductor `×2.6e4…4.4e5`, archimedean `×3.1e4…4.4e5`, prime `×3.1e4…2.0e5`, parity swap `×2.4e4…1.0e5`, χ-substitution `×2.2e4…1.9e5` | **T2** | `src/phase2_gateA_primeflip.py` | `results/phase2_gateA_primeflip.txt` |
| 2 | At `L = 0.8` the prime side **cannot** separate χ₋₄ from χ₈ (identical on all `n ≤ e^{1.6}`); ratio exactly 1.0. Proved; the question it raises is open (§3.5) | **T1** | " | " |
| **3** | **GATE B PASSES**: 0 monotonicity violations in 336 rows, two complete bases, `λ_min > 0` everywhere | **T2** | `src/phase3_gateB.py` | `data/phase3_gateB.csv` |
| 3 | **ζ regression suite: 18/18 PASS**, every value the parent's committed one — Gate B at `N = 16/20/24` both bases, the even/**odd** split at four `L`, `C(0.8)` and `C(1.3)`, all at the parent's printing precision | **T2** | `src/regression_zeta.py` | `results/regression_zeta.txt` |
| **3** | `λ*_χ(L)` is **4 to 27 orders of magnitude above** `λ*_ζ(L)`; cause is the **missing pole**, not the conductor | **T2 (measurement) + T1 (cause)** | `src/phase3_gateB.py` | `data/phase3_gateB.csv` |
| 3 | `L_onset ≈ ½log(qγ₁/2π)` parameter-free to `−0.033 ± 0.067` over five χ, and `−0.0087` for ζ; the height-only predictor is off by `+0.78` and fits with slope **−0.81** | **T2** (T3 for any coefficient) | `src/phase3b_lpred.py` | `data/phase3b_lpred.csv` |
| **4** | `4K(L,γ)` transfers unchanged (χ-free extremal) | **T1** | — | — |
| 4 | `C_χ/4K` spans `0.0087…0.555` at `L = 0.8` vs ζ's 0.673 — **an open observation, §4.5**, not a result of this extension | **T2, open** | `src/phase4_sensitivity.py` | `data/phase4_sensitivity.csv` |
| **4** | **(D4)** `K_even(L,γ) + K_odd(L,γ) = 2L³/3` for **every** γ, both extremals attained; so `K_even(L,0)=0` **forces** `K_odd(L,0) = 2L³/3`, the maximum | **T1** | `src/phase4b_sector_extremal.py` | `data/phase4b_sector_extremal.csv` |
| 4 | `C = 4(F'² + FF'')` at the perturbed minimiser to `7.6e−5 … 2.5e−2` relative, cross term `≤ 12%` (§4.5.1, open) | **T2, open** | `src/phase4c_decomposition.py` | `data/phase4c_decomposition.csv` |
| 4 | **(H1), a HYPOTHESIS not a result** (§4.5.2): the minimiser re-optimises iff `R = 4δ²K/λ₂ > 1`, 18/18 cells, classes separated `×77` | **T2, open** | `src/phase4d_crossover.py` | `data/phase4d_crossover.csv` |
| **5** | **(D5)** a planted real zero contributes `+2A_δ(f)²`, rank-one **PSD**, on even `f`: **no even-sector window computation at any `L` in any basis can exclude a Siegel zero.** Exact in δ, not a vanishing leading coefficient | **T1** | `src/planting.py` | — |
| 5 | 180/180 even-sector rows: response `≥ 0`, **0 detections** | **T2** | `src/phase5_siegel.py` | `data/phase5_siegel.csv` |
| 5 | Odd sector: response `−2S_δ(f)²`, NSD; 180/180 rows negative, **135 detections** with explicit LDL witnesses on 1–8 basis functions | **T2** | " | " |
| 5 | `Θ = 2aᵀQ⁻¹a` monotone in `N`, `∝ δ²` to 2%; odd-sector onset by `L = 1.4`, `N = 12`, for every χ and every `δ ≥ 0.1` | **T2** | `src/phase5_onset.py` | `data/phase5_onset_odd.csv` |
| **5** | **Sign-structure theorem**: a real pair contributes `+2A_δ²` (even) / `−2S_δ²` (odd), rank one, for every `L` and every basis; derivation = parity of `F` at purely imaginary argument | **T1** | `notes/ff_parity.md`, §5.3 | — |
| **V1** | Step 1(a): the planted object was the real **pair** — decided from the committed data, 12/12 stored rows, unanimous; `quartet(0,δ) = 2·real(δ)` to `5.0e−60` | **T2** | `src/verify_siegel_parity.py` | `results/verify_siegel_parity.txt` |
| **V2** | Step 1(b): rank one to `2.5e−60`, eigenvalue `= +2A²`/`−2S²` to `8.5e−61`, correct sign, **150/150 cells**; `A`,`S` by quadrature on the explicit eigenfunction; response **χ-independent** (one eigenvalue per `(L,sector,δ)` group) | **T2** | " | `data/verify_siegel_parity.csv` |
| **V3** | Step 1(b): `K_even + K_odd = 2L³/3` with `K_odd` by **quadrature**, 45 cells, worst `1.6e−61` | **T2** | " | " |
| **V4** | Step 1(c): `real_block(δ=½)` **= pole term entrywise**, 9 cells, `≤ 6.2e−61`, **same sign**, alternative wrong by `O(1)`–`O(10)`; χ pole analogue identically `0.0` while the χ real-pair block has max entry `2.67` | **T2** | " | `results/verify_siegel_parity.json` |
| **V5** | Step 2: function-field signed parity split (FF1)–(FF4), **exact over ℚ**, 25 cells; `J`-even ↔ `+`, `J`-odd ↔ `−`; brief's `u±w` claims refuted in all 20 non-degenerate cells | **T1** | `src/ff_parity_split.py` | `results/ff_parity_split.json` |

## 6.2 What surprised me

Five things did not behave as the ζ case led me to expect.

1. **The γ = 0 blindness is a property of one *sector*, not of the window.** I expected to
   quantify how blind the detector is and to find it blind to leading order. Instead the even
   sector turned out to be blind *exactly and to all orders in δ* — the response is a perfect
   square of the wrong sign — while the odd sector is **maximally** sensitive there, because
   `K_even + K_odd = 2L³/3` forces `K_odd(L,0)` to the γ-uniform maximum. That identity is
   elementary and I did not expect it to be the thing that explained the asymmetry.

2. **The baseline gap is the missing pole, not the conductor.** I assumed `log q · I` would
   be the big structural change. It is the least consequential one: a multiple of the
   identity cannot cancel anything, and `λ*_ζ(L) ≈ 10⁻¹⁷` is a near-cancellation *against the
   pole term*. Remove the pole and `λ*_χ(L)` simply sits where the zeros put it — four to
   twenty-seven orders of magnitude higher.

3. **`C_χ/4K` varies by a factor 64 across five characters at one `L`, and the reason is
   `λ₂`, not `λ₁`.** The parent note attributes `C` to the minimiser's mobility in a
   near-null space; I expected that to be a property of the window. It is a property of the
   *second* eigenvalue in absolute terms, and the criterion `4δ²K > λ₂` gets the moved/frozen
   flag right in all 18 cells with the two classes separated by a factor 77. At `L = 0.8` the
   ζ minimiser moves by a factor **2609** in `|F'|²/K` and **none** of the five χ minimisers
   moves at all.

4. **At `L = 0.8` the window cannot tell χ₋₄ from χ₈.** They agree on every prime power below
   `e^{1.6} = 4.95`, so the prime-substitution control returns ratio exactly `1.0`. At that
   window length the character enters only through `q` and the parity `a`. I had assumed the
   prime twist was the most strongly constrained part of the form; it is the weakest, until
   `2L > log 5`.

5. **In the function field the parity vectors are `u ± ρ^R w`, not `u ± w`.** The flip
   `J : x_i ↦ x_{R−i}` sends `u ↦ ρ^R w`, so it exchanges `u` and `w` only after the `ρ^R`
   weight is restored — equivalently, `J` exchanges the *normalised* `û` and `ŵ` exactly,
   because `‖u‖/‖w‖ = ρ^R` identically. I expected the symmetric/antisymmetric split to be
   the naive one; it is not, and `u+w` is neither an eigenvector nor palindromic for any
   `ρ ≠ 1`. The same identity `‖w‖² = ρ^{−2R}‖u‖²` is what makes `‖u‖‖w‖`, and hence both
   eigenvalues, **rational** while `‖u‖` and `‖w‖` separately are not — which is why Step 2
   can be exact over ℚ at all.

6. **Two of my own checks failed for reasons that were about the *reference*, not the
   computation** (Retractions 1 and 4): the parent's CSVs store 12 significant digits, and its
   Connes window `L = log(13)/2` is stored rounded to 10 — and `λ_min ≈ 10⁻²⁹` there is
   hypersensitive to `L` because `2L` sits exactly at the `n = 13` cutoff. Both looked like
   drift and were not.

## 6.3 Open questions

Five, no more.

1. **Is there a certified analogue on the χ side?** Nothing in this directory is certified
   over ℚ. The parent's §6.1 four-mode witness is exact-rational interval arithmetic; the
   same machinery should transpose, and the natural first target is a certified
   `Q_χ ⪰ m₀ I` on a small odd-sector subspace — which, unlike the even-sector version,
   would have content for real zeros.

2. **Is Zhu's certificate even-sector or full-space?** This is now the sharpest open
   question here, and it is a question about a paper, not a computation. If Zhu's
   `λ*(0.8) ≥ 8.9×10⁻¹⁸` is an infimum over the **full** admissible space, then by §5.3 it
   already constrains real zeros of ζ in the certified range — positivity on the full space
   includes the odd sector, where the response to a real pair is `−2S_δ(f)²`. If it is
   even-sector only, it carries no such information. The parent paper leaves this open
   deliberately; Zhu's paper was not obtained for this extension. Everything the odd-sector
   measurement would need to become an *exclusion* rather than a sensitivity statement is in
   that one bit.
3. **Is the crossover criterion (H1) a theorem?** `4δ²K(L,γ₁) > λ₂(Q)` is validated on 18
   cells. A perturbation-theoretic proof of the threshold, and of what replaces `λ₂` when the
   near-null space has dimension `> 1`, would turn `C(L)` from a measured constant into a
   computable one.

4. **Complex characters.** `ε(χ) ≠ 1`, `Λ(½+it,χ)` is not real, the prime side carries `χ(n)`
   and `χ̄(n)` separately, and the test function must be allowed complex. Does the
   even/odd decomposition that drives (D4) and (D5) survive, or does `χ ≠ χ̄` mix the sectors?

5. **Which other degenerate configurations are sector-protected?** `γ = 0` is one. A zero at
   the central point (`γ = 0` **and** `δ = 0`), a zero of even multiplicity on the line, and
   the `δ → 1/2` limit where `real_block` becomes the pole term are three more places where a
   rank-one response of fixed sign might hide a configuration from one sector.

## 6.4 The answer to the question in the brief

> *Does the finite-window detector behave identically on `L(s, χ)` and ζ up to the baseline,
> or is there a χ-dependent effect beyond the conductor?*

**Not identically. One effect is established, one is the baseline, and one is a reproducible
observation whose mechanism is not.** The *form* transposes with exactly four changes (no
pole; `+log(q/π)` with coefficient one; digamma argument `1/4 + a/2`; `χ(n)`-twisted,
sign-indefinite prime weights), Gate A confirms all four against `Σ_ρ|F(γ_ρ)|²` in both
parity sectors for all five characters, and the extremal `4K(L,γ)` transfers unchanged
because its function class never mentions the L-function.

**The baseline moves enormously** — `λ*_χ(L)` is `10⁴` to `10²⁷` times `λ*_ζ(L)` — and the
cause is the **missing pole**, not the conductor: `λ*_ζ` is the residue of a near-cancellation
against a positive rank-one pole term, and `log q · I`, a multiple of the identity, cannot
play that role. The resolution scale tracks the **analytic conductor** `qγ₁`, parameter-free
to `−0.033 ± 0.067` across all six L-functions including the `q = 1` case (§3.4).

**The one effect established as a theorem is not about the conductor at all: the `γ = 0`
corner is sector-protected.** A real zero's contribution to the window form is rank one and
**positive** semidefinite on even test functions (`+2A_δ²`) and **negative** semidefinite on
odd ones (`−2S_δ²`), for every `L` and every basis — verified rank one to `10⁻⁶⁰`, with the
predicted eigenvalue to `10⁻⁶¹`, in 150 of 150 cells, and independently in the function field
exactly over ℚ. So even-sector finite-window positivity is compatible with a real zero at
every `δ`, and even-sector certificates carry no information about the Siegel configuration;
the odd sector, which `K_even + K_odd = 2L³/3` makes maximally sensitive exactly there,
detects one by `L = 1.4` with a 12-dimensional basis. **This is the result worth carrying
forward, and it is a statement about the detector, not about any `L(s, χ)`.** Its one
load-bearing loose end is whether Zhu's certificate is even-sector or full-space — unresolved
here, and consequential in opposite directions (§5.3, non-claim 8).

**Beyond those, there is a genuine χ-dependence that this extension does *not* explain.** The
saturation of the transferred extremal, `C_χ/4K`, spans a factor 64 at a single `L` against
ζ's 0.673. A criterion that fits it — the minimiser re-optimises iff `4δ²K(L,γ₁) > λ₂(Q)`,
18/18 cells — is recorded as **hypothesis (H1)** in §4.5 and not as a result: five real
primitive characters of conductor `≤ 8` do not vary `q`, `a`, `λ₂` and `γ₁` independently,
and at `L = 0.8` the window cannot even separate χ₋₄ from χ₈ (§3.5). Calling that effect
"the conductor" would be the overreach this report is trying to avoid.

---

## 7. The note

After the record above was complete, a standalone note was written:
[`note/siegel_note.tex`](note/siegel_note.tex) → `note/siegel_note.pdf`, 10 pages,
**doi:10.5281/zenodo.23623473** — deposited on Zenodo as v1 on 2026-10-24 under
CC BY 4.0. It carries the same non-claims, the same
tiers, the same retraction log and a reproduction appendix, and it states the one result
worth carrying forward: **at a real zero the even sector cannot respond at all, and the odd
sector responds with its largest possible negative value.** Committed as `1f7d3ca`, tag
`siegel-note-v1`, which was moved onto the DOI commit so that the tagged state and the
deposited PDF carry the same DOI rather than a placeholder.

`make -C note verify` re-derives every number printed in the note from this extension's
data files, and fails on any disagreement, any undeclared literal, any overstatement of
precision relative to the registry's own significant digits, a bibliography entry not found
in this report or the parent paper, or — pass 3 — a prerequisite that has not been met: it
reads `results/verify_siegel_parity.txt` and refuses to pass unless the three Siegel-corner
checks of §5.4 are recorded there as passed.

Two things were **corrected rather than carried over** while writing it, and they are logged
as R6 in the note's Appendix B:

1. The function-field eigenvectors are `v± = u ± ρ^R w`, with `λ± = (R+1) ± ρ^{−R}‖u‖²`.
   The form `u ± w` holds in exactly 5 of the 25 cells — the `ρ = 1` ones — and fails for
   every `ρ ≠ 1`. This is P3 of §0.3, restated because a second brief repeated it.
2. Zhu's certificate is still **not** established as an even-sector infimum, so §5 of the
   note states both cases and claims neither. This is P4 of §0.3, likewise repeated.

One wording change was made here and in `README.md`: "no Zhu-style certificate **exists**
for any Dirichlet window" overstates an absence, and now reads "**is known**". The note
says the same, and adds the sharper statement the sentence should have rested on: for
general real `f = f_e + f_o` a real pair contributes `2(A_δ(f_e)² − S_δ(f_o)²)`, so lowering
the form requires a nonzero **odd component** — a condition on test functions rather than a
claim about what this framework can contain.
