# Paper outline (spec for `paper/weil_window_control.tex`)

**Title.** *A δ²–L³ Sensitivity Law for Finite-Window Weil Positivity, with a
Function-Field Control of the Connes–van Suijlekom Pipeline*

**Author.** Elias De Jesús, Independent Researcher. ORCID 0009-0007-0190-9143.
October 2026 — Technical Note.

**Tiers.** T1 exact/proved · T2 numerical, precision stated, convergence checked ·
T3 interpretation. Every statement in the paper carries one, with its script and data file.

---

## Novelty policy (binding on the whole paper)

Only two things are claimed as new:

1. **The δ²–L³ sensitivity law** — a theorem in the function field, an unconditional bound
   and a sharper conditional one on ζ, measured near-saturation of the sharp constant, and
   the shared `(window length)³/6` form.
2. **The planted-spectrum control of the CvS minimiser pipeline in the function field** —
   running the whole pipeline (ground state, its zeros, simplicity, parity, convergence) on
   spectra deliberately moved off the critical circle.

Everything else is **foundation or confirmation** and is cited as such:
Hallouin–Perret (Prop. 5, Thm 6, Thm 36(ii), Lemmas 33–34) is the function-field
foundation; **Bombieri (2000) §13 is the precedent for planted ζ experiments** and already
exhibits a critical window `t_c`; CvS supply Steps 1 and 5 and state the simplicity
difficulty themselves (Remark 2.3, Introduction); Zhu supplies the certified `λ*(0.8)` and
the bandwidth `T*(L) = 2πe^{2L}`.

---

## Abstract (≤250 words) — drafted in the .tex

## Non-claims — immediately after the abstract, not an appendix

1. No proof or advance of RH for ζ.
2. No new spectral architecture: the window form, its Toeplitz structure and its
   Carathéodory–Fejér theory are Hallouin–Perret's and CvS's.
3. Planted spectra are fabricated diagnostics, not L-functions.
4. Numerical detections are T2 unless interval-certified; non-detections certify nothing.
5. `λ*(0.8) ≤ 2.2702e−17` is a variational upper bound, not a certified enclosure.
6. Stage 3b reproduces Groskin's computation at smaller scale.
7. Block F is claimed **only** for the function field (Retraction 15).
8. Davenport–Heilbronn not attempted; no Lean formalisation.

---

## Section map

| § | Content | Tier | Script | Data |
|---|---|---|---|---|
| 1 | The pipeline, the four steps, the question | — | — | — |
| 2 | **Foundation.** HP Prop. 5: the window Toeplitz matrix *is* the Frobenius Gram matrix (`2g` diagonal, `x_n` off). HP Thm 6: `d` = first singular size, kernel generator. HP Thm 36(ii) = step (i). Stated as prior art. | T1 (cited) | — | — |
| 3 | **Control experiment, genuine spectra.** Four curves, exact point counts, `λ_min(R)`, ground-state zeros, parity, `λ_min` monotonicity. | T1+T2 | `stage1.py`, `curves.py`, `ffield.py` | `stage1_spectra.csv` |
| 4 | **Control experiment, planted spectra (CLAIM 2).** First-negative window certified exactly over Q; H1 holds at every planted row; simplicity and parity never fail; masking is weak. Zeros stay on the circle while `λ_min = −642`. | T1+T2 | `stage2.py`, `exact_inertia.py`, `adversarial.py` | `stage2_*.csv` |
| 5 | **Degeneracy beyond `2g`** — corollary of HP's rank structure (`ker T_R = {P*q}`, `dim = R+1−d`, parity split), aimed at the CvS simple-ground-state hypothesis. CvS Remark 2.3 states it in words. | T1+T2 | `function_level.py` | `stage2_function_level.csv` |
| 6 | **Step (iv) at the function level** — Hurwitz obstruction, presented as making CvS's Step 5 concrete. Angles converge; functions do not. | T1+T2 | `function_level.py`, `figs_function_level.py` | `stage2_function_level.csv` |
| 7 | **The δ²–L³ law (CLAIM 1).** Thm A (function field, T1, with the `log ρ` refinement and `O(a⁴)`); Thm B (ζ, unconditional for `L ≤ 0.8`); Thm B′ (ζ, conditional, sharp `8L³/3`); measured `C/(4L³/3) → 0.945`; the shared `(window length)³/6`. | T1 + T2 | `proofs_check.py`, `delta2_law.py`, `task34_extra.py` | `stage3_C_vs_CS.csv`, `stage3_sharp_bound.csv` |
| 8 | **ζ replication as a control.** Gate A (assembly vs zeros side, with the pole-sign control); Gate B (`λ*(0.8) ≤ 2.2702e−17`, two bases); Connes's `[1,13]` window at `N = 36`. | T2 | `stage3_converge.py`, `stage3_connes.py` | `stage3_gate*.csv`, `stage3b_*.csv` |
| 9 | **Detection onset.** δ-resolved and basis-resolved refinement of **Bombieri's `t_c`** (§13, 2000) and **Zhu's `T*(L)`**. Onsets 1.3×, 1.0×, 1.0×, 0.85× `L_pred`; onset moves earlier with basis; `N ≳ 2n` needed. | T2 | `stage3e_certify.py`, `task1e_C_vs_height.py` | `stage3e_*.csv` |
| 10 | **Recorder split.** Structural Block F in the function field (rank-2 hyperbolic pole part); **withdrawn on ζ** — both `A` and `M` positive definite under the other convention. Convention stated. | T2 | `task2_inertia_reconcile.py` | `stage3d_inertia_reconciled.csv` |
| 11 | Relation to prior work | — | — | `PRIOR_ART.md` |
| 12 | Limitations | — | — | — |
| A | Epistemic tiers table | — | — | — |
| B | Retractions log (16) and the traps they document | — | — | — |
| C | Reproduction table | — | — | — |

**Traps documented in Appendix B** (the methods value of the retractions): zeros-sum
truncation bias; the two counterfactuals (moving vs adding a zero) and why the sharp law
needs the first; float64 spurious negatives; interval-library gaps (no digamma/Lerch in
`mpmath.iv`); the certification asymmetry (non-detection certifies nothing); regression
values as the only thing that caught three sector-rename bugs.

## Figures

`stage1_zeros_circle` · `stage2_zeros_stay_on_circle` · `stage2_iv_rescored` ·
`stage2_detection_law` · `stage3_gates` · `stage3_connes_profile` · `stage3_teeth` ·
`stage3_teeth_height`.

## Citations

Hallouin–Perret 2019; Connes–van Suijlekom 2511.23257; Connes–Consani–Moscovici
2511.22755; Connes 2602.04022; Zhu 2608.24827; Suzuki 2606.09096; Groskin 2605.20224 and
2607.02828; Bombieri 2000; Yoshida; Connes–Consani 2021; and the author's two prior Zenodo
notes (Partition Potential DOI 10.5281/zenodo.21108392; Weil Enforcer — DOI placeholder).
