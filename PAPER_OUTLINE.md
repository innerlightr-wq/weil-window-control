# Paper outline (structure only — no prose drafted yet)

**Working title.** *Where the RH content sits in the finite-window Weil-form program:
a function-field control with planted spectra*

**Status.** Outline only, per instruction. Every row below points at the script and data
file that produced it. Tiers: **T1** exact/proved · **T2** numerical with stated precision
and a convergence check · **T3** interpretation.

---

## Abstract (draft, 190 words)

The Connes–van Suijlekom and Connes–Consani–Moscovici programs reduce the Riemann
hypothesis to properties of a quadratic form restricted to test functions supported in a
finite window. We ask which step of that pipeline carries the arithmetic content, by
running it where the answer is known — curves over finite fields, where RH is Weil's
theorem — and against fabricated spectra with eigenvalues moved off the critical circle. We find that two of the four steps are free of RH content: the reality of
the approximant's zeros (Carathéodory–Fejér, and already implicit in the CvS hypothesis
class), and the simplicity and evenness of the ground state, neither of which is disturbed
by a planted violation. The content sits in the sign of the lowest eigenvalue, and in
convergence of the normalised minimiser *as a function* — the two being linked, since the
free witness is exactly what obstructs Hurwitz convergence to an off-line target. We
identify the window Toeplitz matrix with the Rosati Gram matrix on H¹ of the curve,
and derive the sensitivity law λ_min ≈ −C·δ², with C asymptotically the Cauchy–Schwarz
value (window length)³/6 on both sides.
**Nothing here bears on RH for ζ.**

---

## Non-claims (to appear immediately after the abstract, not buried)

1. Nothing in this work proves, advances, or provides evidence for RH for ζ.
2. Stage 1 runs where RH is a **theorem** (Weil, curves over finite fields).
3. Stage 2 uses **fabricated** spectral data, built by hand to violate `|β| = 1`.
4. Stage 3's Gate A **assumes** RH (it compares against a sum over zeros written
   `1/2 + iγ`, γ real); it validates the assembly, nothing more.
5. Stage 3e is a **diagnostic** built from zero data with a hand-planted quartet.
6. Stage 3b **reproduces at smaller scale** a published computation (Groskin 2605.20224).
7. `λ*(0.8) ≤ 2.2702e−17` is a **variational upper bound** from a finite basis — no
   rigorous lower bound, and the discretisation error is observed monotone, not bounded.
8. The novelty check (§8) was done at **abstract level only**, not by reading the papers.
9. Stage 3f (Davenport–Heilbronn) was **not attempted**; the reason is stated.
10. No Lean formalisation was attempted.

---

## Section map

### §1 Setup and the four steps
The CvS/CCM pipeline; the four candidate loci (i) real zeros, (ii) sign of λ_min,
(iii) simplicity/evenness, (iv) convergence. Statement of the question.
*No results; framing only.*

### §2 The function-field control
| Result | Tier | Script | Data |
|---|---|---|---|
| Window identity `t(n) = q^{n/2}+q^{−n/2}−N_n q^{−n/2}`, derived | T1 | `notes/stage1_derivation.md` | — |
| Weil positivity automatic: `T_R = Σ w_m w_m^H`, `inertia_neg = 0` in all 40 rows | T1 + T2 | `src/stage1.py` | `results/stage1_spectra.csv` |
| Four curves, exact integer point counts; `N₁=6, N₂=18` for the genus-2 curve | T1 | `src/curves.py`, `src/ffield.py` | `results/stage1_curves.json` |
| H3: `N_1..N_g` determine everything (4–9 exact integer predictions per curve) | T1 | `src/stage1.py` | `results/stage1_curves.json` |
| H1 confirmed and RH-free; survives 400 random Toeplitz | T1 + T2 | `src/adversarial.py` | — |
| H2 **restated** (needs distinct θ_j); counterexample `y²=x⁵+x/F₅`, `P=(1+5T²)²` | T1 | `src/adversarial.py` | — |
| `R=2g` kernel vector in closed form `c_k ∝ A_k q^{−k/2}`, forcing evenness | T1 | `src/stage1.py` | `results/stage1_spectra.csv` |
| H4 **downgraded**: not monotone; `λ_min(R)` is the monotone quantity | T2 | `src/stage1.py` | `results/stage1_spectra.csv` |

### §3 The kernel is an ideal
| Result | Tier | Script | Data |
|---|---|---|---|
| `ker T_R = {P* q : deg q ≤ R−d}`, `dim = R+1−d`, matched in every row | T1 + T2 | `src/function_level.py` | `results/stage2_function_level.csv` |
| CvS hypothesis fails for every `R > d`; even/odd kernel dims `⌈(m+1)/2⌉,⌊(m+1)/2⌋` | T1 + T2 | `src/function_level.py` | `results/stage2_function_level.csv` |
| Corollary: the "evenness failure" past `2g` is solver-arbitrary degeneracy | T3 | — | — |

### §4 The Toeplitz–Rosati bridge
| Result | Tier | Script | Data |
|---|---|---|---|
| `T_R` = Gram matrix of `{1, F/√q, …, F^R/q^{R/2}}` under `Tr(φψ†)`; residuals ≤1.1e−39 | T1 + T2 | `src/gram_trace.py` | `results/stage1e_gram_trace.csv` |
| `F*SF = qS` forces `S_jj = 0` when `|β_j| ≠ √q` ⇒ positive-definite polarization iff RH | T1 (**classical**) | `src/gram_trace.py` | — |
| On planted spectra the Gram identity survives, `S` goes indefinite with inertia **equal** to the Stage 2 Toeplitz signature: (1,0,1), (2,0,2), (2,0,4) | T2 | `src/gram_trace.py` | `results/stage1e_gram_trace.csv` |

### §5 Teeth: planted off-circle spectra
| Result | Tier | Script | Data |
|---|---|---|---|
| First-negative window, **certified exactly over Q** (Sylvester/Jacobi on rational surrogates) | T1 | `src/exact_inertia.py`, `src/stage2.py` | `results/stage2_spectra.csv` |
| H1 holds in every planted row, including `λ_min = −642` — the witness is free | T1 + T2 | `src/stage2.py` | `results/stage2_spectra.csv` |
| Simplicity and definite parity never fail under planting | T2 | `src/stage2.py` | `results/stage2_spectra.csv` |
| Exact closed form `λ_min = (R+1) − |u||w|`, verified to 4.5e−80 | T1 | `src/detection_law.py` | `results/stage2_detection_law.csv` |
| Masking by on-line zeros is weak (`R_det` 2→12 as `k` 0→24) | T1 arithmetic, **T3/toy** transfer | `src/adversarial.py` | `results/stage2_masking.csv` |

### §6 Re-scoring step (iv) at the function level
| Result | Tier | Script | Data |
|---|---|---|---|
| **Theorem**: H1 + Hurwitz forbid convergence to a target with an off-circle zero | T1 | — | — |
| Zero *angles* converge (~1e−3 at R=30) even for off-line spectra — free | T2 | `src/figs_function_level.py` | `results/stage2_zeros.csv` |
| Coefficient distance to target at `R=d` is exactly `√2` (orthogonal); distance to the ideal exactly 1 | T1 + T2 | `src/function_level.py` | `results/stage2_function_level.csv` |
| `|ĉ_R(β)|` normalised **rises** to 0.58–0.71 out to `R=36` instead of falling | T2 | `src/function_level.py` | `results/stage2_function_level.csv` |
| Verdict table: (i) free, (ii) content, (iii) free, (iv) free as ordinates / content as functions | T3 | — | — |

### §7 The ζ side as a control
| Result | Tier | Script | Data |
|---|---|---|---|
| Assembly: closed-form `F`, `C`, archimedean `m`-sum via digamma/Lerch/Hurwitz; also as one 1-D integral | T1 | `notes/stage3_assembly.md`, `src/zeta_window.py`, `src/geom_at_vector.py` | — |
| Odd-sector pole sign `−2P_jP_k` derived (brief gives only the even case) | T1 | `notes/stage3_assembly.md` | — |
| **Gate A**: geometric side vs `Σ_ρ|F(γ_ρ)|²` to 1.6e−7 at K=800, decaying at the rate each basis's smoothness predicts; pole-sign control blows it to 5.35 | T2 | `src/stage3_converge.py` | `results/stage3_gateA_explicit_formula.csv` |
| **Gate B**: `λ*(0.8) ≤ 2.2702e−17`, monotone in two independent complete bases; Zhu's normalisation confirmed identical | T2 | `src/stage3_converge.py` | `results/stage3_gateB_lambda_min.csv` |
| Connes's `[1,13]` window: `λ_min = 8.98e−52` at `N=36`/140 digits, `γ₁` to 7.2e−48, zeros tracked to `n≈20` | T2 | `src/stage3_connes.py` | `results/stage3b_connes_zeros_N36.csv` |
| Ground state simple and even at every `L`; gap `2e4`–`4e5` | T2 | `src/stage3_connes.py` | `results/stage3c_even_odd.csv` |
| Block F on ζ: `A`, `M` indefinite, `S` PSD; `A`'s negative index depends on `L` **only** (not `N`, not basis) | T2 | `src/task34_extra.py` | `results/stage3d_inertia_vs_N.csv` |
| `λ_min ≈ −C(L)δ²`; `C/(4L³/3)` rises 0.46→0.95; Cauchy–Schwarz saturated | T1 bound + T2 | `src/delta2_law.py`, `src/task34_extra.py` | `results/stage3_C_vs_CS.csv` |
| Detection onset tracks `L_pred = ½log(γ/2π)` and moves earlier with basis size; `γ₃₀` detected at `0.85 L_pred` | T2 | `src/stage3e_certify.py` | `results/stage3e_certified_onset.csv` |
| `C(L,γ*)` collapses ~12 orders as `γ*/T*` goes 0.3→0.75 | T2 | `src/task1e_C_vs_height.py` | `results/stage3e_C_vs_height.csv` |

### §8 Literature and novelty
Reproduce the Task 6 table verbatim, **with its abstract-level caveat in the same
paragraph, not a footnote.** Confirmations: H1 = CvS Step 1 (Carathéodory–Fejér 1911) and
Pisarenko; H2 = classical exact recovery; Hurwitz = CvS Step 5; Rosati positivity =
classical; §7's Connes window = Groskin 2605.20224 at smaller scale.
Possibly new (subject to the caveat): the function-field control itself, the
Toeplitz–Rosati bridge, the `δ²`/`L³` sensitivity law on both sides, the height
dependence, and the obstruction half of the Hurwitz argument.

### §9 What would have to be true next
Honest limits: no rigorous lower bound on `λ*`; `N ≳ 2n` blocks high planted heights;
the ζ-side sensitivity statements are conditional on planting height, basis and `K`;
Davenport–Heilbronn untouched.

---

## Retractions (to be included as a numbered appendix, not omitted)

Fourteen are logged in `REPORT.md`. The ones a reader must see:

| # | What was retracted |
|---|---|
| 1 | Frobenius radius computed on `P` rather than its reverse (`0.8` instead of `≤1.3e−51`) |
| 3 | Tolerance-based sign test misclassified exact zeros; replaced by an error-based test plus exact rational certificates |
| 4 | **H2 as originally stated was refuted by my own adversarial pass**; restated with a distinctness hypothesis |
| 5 | H4 as a monotone-error claim not supported |
| 6 | **"(iv) is free" retracted** — scored on zero angles, which is the wrong object; at the function level (iv) carries content |
| 8–9 | A 40-digit Stage 3e run discarded; three sector-rename bugs caught by regression values |
| 11 | **"Low-lying zeros only" retracted as structural** — the Round-1 height sweep used a truncated on-line side (biased toward detection) and an ADDED rather than MOVED quartet (wrong sign away from a zero). Redone: `γ₃₀ ≈ 101` is detected at `0.85 L_pred` |
| 12 | `γ₈₀` "detections" below their own noise floor discarded; reported under-resolved, not undetectable |
| 14 | The Task-4b expectation (basis-dependent `A` inertia) was refuted and the refutation reported |

---

## Figures

`stage1_zeros_circle` · `stage2_zeros_stay_on_circle` (the witness-is-free figure) ·
`stage2_iv_rescored` (angles converge / functions do not) · `stage2_detection_law` ·
`stage3_gates` · `stage3_connes_profile` · `stage3_teeth` · `stage3_teeth_height`.
*To add:* onset vs `L_pred` across heights and basis sizes; `C/(4L³/3)` vs `L`.
