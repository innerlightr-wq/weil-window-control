# weil-window-control

A quantitative sensitivity analysis of the **finite-window Weil quadratic form** used in the
Connes–van Suijlekom and Connes–Consani–Moscovici programs, calibrated in a setting where
the answer is known. For a curve over a finite field the window form is *exactly* the
Frobenius Gram matrix of Hallouin–Perret (Trans. AMS **372** (2019), Prop. 5) and the
Riemann hypothesis is Weil's theorem; we use that as a control, run the whole minimiser
pipeline on it with both genuine and deliberately planted off-circle spectra, and prove a
sensitivity law — moving a zero a distance `δ` off the critical line costs the lowest
eigenvalue `≈ C·δ²` with `C ≍ (window length)³/6` — which holds in the same form on the
ζ side. **Nothing here proves or advances the Riemann hypothesis for ζ.**

📄 **Paper (current):** [`paper/weil_window_control_rev_certified_witness.pdf`](paper/weil_window_control_rev_certified_witness.pdf) — [doi:10.5281/zenodo.23250932](https://doi.org/10.5281/zenodo.23250932) · all versions: [10.5281/zenodo.21109955](https://doi.org/10.5281/zenodo.21109955) · source [`paper/weil_window_control.tex`](paper/weil_window_control.tex) ·
🔍 **Prior art, verified in full text:** [`PRIOR_ART.md`](PRIOR_ART.md) ·
🧾 **Full lab record incl. 20 retractions:** [`REPORT.md`](REPORT.md)

> The previous version, [`paper/weil_window_control.pdf`](paper/weil_window_control.pdf)
> ([doi:10.5281/zenodo.23212998](https://doi.org/10.5281/zenodo.23212998)), is kept as a
> historical record. Cite the current DOI above.

---

## The certified four-mode witness (§6.1, Theorem 11)

The current version adds **one computer-assisted certified result**. Fix the window
half-length `L = 4/5`, the height `γ = 14` exactly, and the four quarter-wave modes
`φ_k(x) = cos((2k+1)πx/2L)/√L`, `k = 0..3`, on `[−L, L]` and zero outside. For one **frozen
rational** coefficient vector `x`
([exact values and enclosures](extensions/planted-witness-consolidation-2026-10-08/cold_certificate.json)), write

```
Q_add,δ(f) = Q_ζ(f) + 4(A_δ² − B_δ²),
A_δ = ∫ f(u) cos(γu) cosh(δu) du,   B_δ = ∫ f(u) sin(γu) sinh(δu) du,
```

the finite-window Weil form after **adding** a quartet of zeros at `1/2 ± δ ± iγ` to the
zero multiset, the true zeros left in place. Then, as certified enclosures:

```
Q_ζ(f)      > 0
Q_add,2/5(f) < −7/500
Q_add,δ(f)   < −1/500     for every δ in [1/4, 49/100]
```

So this four-mode window is positive on the actual zero set and strictly negative once an
off-line quartet is added at displacement `δ ∈ [1/4, 1/2)`. The bounds above are rounded
outward for legibility; the sharp rational enclosures are in
[`cold_certificate.json`](extensions/planted-witness-consolidation-2026-10-08/cold_certificate.json).

**Evidence, stated precisely.**

* The point `δ = 2/5` was checked by **two separate implementations** with different
  archimedean closed forms, different error budgets and no shared code — the second obtains
  an end-to-end negative upper bound of its own, not merely matching intermediates.
* Uniform negativity on `[1/4, 49/100]` was certified by the **primary** interval
  implementation over **256 abutting covering cells** (not sampled points), and
  **cold-recomputed** from the frozen input with no cache present.
* All arithmetic is exact rational interval arithmetic with outward rounding; `π`, `log`,
  `exp`, `sin`, `cos`, `√` are truncated series with proved remainders. No libm value enters
  the certified path.
* This is a **computer-assisted proof with ordinary runtime dependencies** (CPython, its
  `int`/`Fraction` implementation, the OS). It is **not** proof-assistant formal
  verification, and it has **not** been externally peer-reviewed.

**Scope.** The quartet is **artificial** — a counterfactual added to the zero multiset.
Nothing here restricts the location of any actual zero of `ζ`, and **nothing here proves or
advances the Riemann hypothesis**. It is not an optimal threshold (`δ ∈ (0, 1/4)` is
uncovered) and not uniform in height (`‖P_E r_γ‖²` falls from `3.60e−2` at `γ = 14` to
`1.13e−7` at `γ = 280`).

| | |
|---|---|
| Theorem and analytic audit | [`CERTIFICATE_THEOREM.md`](extensions/planted-witness-closure-2026-10-08/CERTIFICATE_THEOREM.md), [`ANALYTIC_AUDIT.md`](extensions/planted-witness-consolidation-2026-10-08/ANALYTIC_AUDIT.md) |
| Primary verifier | [`verify_certificate_consolidated.py`](extensions/planted-witness-consolidation-2026-10-08/verify_certificate_consolidated.py) |
| Separate implementation | [`independent_point_check.py`](extensions/planted-witness-consolidation-2026-10-08/independent_point_check.py) |
| Results | [`cold_certificate.json`](extensions/planted-witness-consolidation-2026-10-08/cold_certificate.json), [`interval_cells.json`](extensions/planted-witness-consolidation-2026-10-08/interval_cells.json) (all 256 cells), [`independent_point_check.json`](extensions/planted-witness-consolidation-2026-10-08/independent_point_check.json) |
| Verification report | [`VERIFICATION_REPORT.md`](extensions/planted-witness-consolidation-2026-10-08/VERIFICATION_REPORT.md), [`REPRODUCE.md`](extensions/planted-witness-consolidation-2026-10-08/REPRODUCE.md) |

### Quick start (certificate only — standard library, no dependencies)

```bash
git clone https://github.com/innerlightr-wq/weil-window-control.git
cd weil-window-control/extensions/planted-witness-consolidation-2026-10-08

python3 audit_primitives.py            # interval-arithmetic properties   ~4 s
python3 audit_cache_invalidation.py    # cache cannot go stale            ~45 s

rm -f regenerated_qzeta_cache.json cold_certificate.json interval_cells.json
python3 verify_certificate_consolidated.py   # COLD, no cache             ~5.5 min
python3 independent_point_check.py           # separate implementation    ~3.5 min
```

Each exits `0` on success and **non-zero** if any gate fails; the gates assert the
inequalities themselves, they do not print a stored label. `make witness-quick`,
`make witness-cold` and `make witness-independent` run the same things from the repository
root.

---

## Non-claims

Read these before the results.

1. **No proof or advance of RH for ζ.** Part I runs where RH is a *theorem* (Weil, curves
   over finite fields); Part II uses fabricated spectra.
2. **No new spectral architecture.** The window form, its Toeplitz structure and its
   Carathéodory–Fejér theory are due to Hallouin–Perret and Connes–van Suijlekom; the
   finite-window Weil form goes back to Yoshida (1992) and Bombieri (2000).
3. **Planted spectra are fabricated diagnostics**, not *L*-functions of anything.
4. **Numerical detections are tier T2 unless certified, and non-detections certify
   nothing** — a finite basis that finds no negative direction is not evidence of
   positivity.
5. `λ*(0.8) ≤ 2.2702e−17` is a **variational upper bound**, not a certified enclosure. The
   certified **lower** bound `λ*(0.8) ≥ 8.9e−18` is Zhu's, used only on its documented range
   (`L ≤ 0.8`) and normalisation (`Q(f)/‖f‖²`). Upper and lower bounds are not interchanged
   anywhere. Outside that range the ζ bound is **relative** to `λ*(L)`; replacing that
   baseline by zero is not asserted (Retraction 18).
6. The ζ-side replication **reproduces published computations** (Groskin) at smaller scale.
7. **Block F is claimed only for the function field** — on ζ windows the indefiniteness of
   `A` and `M` depends on where `−log π` is placed (Retraction 15).
8. **Not attempted:** Davenport–Heilbronn; any Lean formalisation. The §6.1 certificate is
   exact-rational interval arithmetic run by CPython, **not** a proof-assistant artifact.

## What is claimed as new

Only two things, and both are qualified in [`PRIOR_ART.md`](PRIOR_ART.md):

- **(N1)** the `δ²`–`L³` sensitivity law (function-field theorem, ζ-side bounds, measured
  near-saturation, shared `(window length)³/6` form);
- **(N2)** the planted-spectrum control of the CvS minimiser pipeline in the function field.

Everything else is foundation or confirmation.

---

## Results, by tier

| Result | Tier | Script | Data |
|---|---|---|---|
| Window form = Hallouin–Perret Gram matrix; Weil positivity automatic | T1 (cited) + T2 | `src/stage1.py` | `data/stage1_spectra.csv` |
| Planted spectra: first-negative window, **certified exactly over ℚ** | **T1** | `src/stage2.py`, `src/exact_inertia.py` | `data/stage2_spectra.csv` |
| Step (i) holds at every planted row (incl. `λ_min = −642.8`) | T1 + T2 | `src/stage2.py` | `data/stage2_spectra.csv` |
| `ker T_R = {P·q}`, `dim = R+1−d`; CvS hypothesis fails for `R > d` | T1 + T2 | `src/function_level.py` | `data/stage2_function_level.csv` |
| Hurwitz obstruction: angles converge, functions do not | T1 + T2 | `src/function_level.py` | `data/stage2_function_level.csv` |
| **`λ_min = −C(R+2,3)·a² + O(a⁴)`, `a = log ρ`** | **T1** | `src/proofs_check.py` | — |
| **ζ bound, *relative*: `λ_min ≥ λ*(L) − c(L)δ² − O(δ⁴)`, `c(L) = 8L³(⅓+1/√5)`; unconditional where positivity is certified (`L ≤ 0.8`)** | **T1** | `src/proofs_check.py` | `data/stage3_sharp_bound.csv` |
| Sharper conditional bound, constant `8L³/3` | T1 (conditional) | `src/proofs_check.py` | `data/stage3_sharp_bound.csv` |
| `C/(4K(L,γ)) → 0.963`: the exact constant nearly saturated | T2 | `src/task_final_checks.py` | `data/stage3_exact_constant_K.csv` |
| Gate A: geometric side vs `Σ_ρ|F(γ_ρ)|²` to 1.6e−7 | T2 | `src/stage3_converge.py` | `data/stage3_gateA_explicit_formula.csv` |
| Gate B: `λ*(0.8) ≤ 2.2702e−17`, two independent bases | T2 | `src/stage3_converge.py` | `data/stage3_gateB_lambda_min.csv` |
| Connes `[1,13]` window: `λ_min = 8.977e−52` at `N=36` | T2 | `src/stage3_connes.py` | `data/stage3b_connes_zeros_N36.csv` |
| Detection onset tracks `L_pred = ½log(γ/2π)`, moves earlier with basis | T2 | `src/stage3e_certify.py` | `data/stage3e_certified_onset.csv` |
| Recorder split; Block F structural in the function field only | T2 | `src/task2_inertia_reconcile.py` | `data/stage3d_inertia_reconciled.csv` |
| **Four-mode added-quartet witness: `Q_ζ>0`, `Q_add,δ<0` on `[1/4,49/100]`, certified over ℚ** | **T1** | `extensions/planted-witness-consolidation-2026-10-08/verify_certificate_consolidated.py` | `extensions/planted-witness-consolidation-2026-10-08/cold_certificate.json` |
| Same point `δ=2/5`, separately reimplemented (different archimedean closed form) | **T1** | `extensions/planted-witness-consolidation-2026-10-08/independent_point_check.py` | `extensions/planted-witness-consolidation-2026-10-08/independent_point_check.json` |

---

## Reproduction

```bash
git clone https://github.com/innerlightr-wq/weil-window-control.git
cd weil-window-control
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
make quick          # self-test, ~2 min — start here
```

| Target | What it runs | Runtime* |
|---|---|---|
| `make quick` | exact arithmetic + function-field proof checks | **~2 min** |
| `make stage1` | function-field control, genuine spectra | ~2 min |
| `make stage2` | planted spectra, exact over ℚ, detection law, adversarial pass | ~6 min |
| `make gates` | ζ Gate A and Gate B | ~45 min |
| `make connes` | Connes `[1,13]` window at `N=36`/140 digits, even/odd, split | ~60 min |
| `make law` | `δ²`–`L³` law: proof checks and constants | ~25 min |
| `make onset` | detection onset vs `L_pred` and basis size | **hours** |
| `make paper` | compile the PDF (needs a TeX distribution) | ~10 s |
| `make all` | everything except `onset` | ~2.5 h |

\* measured on one core of an Apple Silicon laptop, CPython 3.12, **without** `gmpy2`.
Installing `gmpy2` speeds the high-precision targets up substantially and changes no result.

**Software versions used for the committed data:** Python 3.12.2, mpmath 1.3.0,
numpy 1.26.4, matplotlib 3.9.2, no gmpy2. `make paper` additionally needs `pdflatex`.

**Numerical policy.** Every eigenvalue sign is produced by `mpmath` at a stated precision
with a precision-doubling or noise-floor check, or certified exactly over ℚ. `float64` is
computed only as a deliberate contrast: it returns spurious negatives of size `1e-16` where
the true value is exactly `0`.

---

## Layout

```
src/        implementation modules
scripts/    one runnable script per stage, plus quick_check.sh
data/       all CSV results and the curve metadata
figures/    generated figures
paper/      LaTeX source and the compiled PDF
notes/      derivations (stage1, stage3 assembly, proofs) and a revision note
extensions/ follow-on work; the added-quartet certificate is in the note (§6.1)
```

## Extensions

Follow-on work kept in this repository. It is preregistered and tiered the same way, and
carries no DOI of its own.

**One extension is now part of the paper.** The added-quartet certificate below is §6.1 /
Theorem 11 of the current note ([doi:10.5281/zenodo.23250932](https://doi.org/10.5281/zenodo.23250932)); the rest are
explorations outside it and have had no review.

- [`extensions/planted-witness-consolidation-2026-10-08/`](extensions/planted-witness-consolidation-2026-10-08/) — **IN THE PAPER (§6.1).** The certified four-mode witness: cold
  recomputation, the separately implemented point check, the interval-primitive and
  cache-invalidation audits, the analytic audit, and all 256 covering cells. Summarised at
  the top of this README.
- [`extensions/planted-witness-closure-2026-10-08/`](extensions/planted-witness-closure-2026-10-08/)
  — **provenance for the above.** The first full-form certificate and its verifier, which the
  consolidated run re-derives and whose files the consolidation's `INPUT_SHA256.txt`
  fingerprints. Kept so that manifest is checkable.
- [`extensions/certified-detection-witness-2026-10-08/`](extensions/certified-detection-witness-2026-10-08/)
  — **provenance for the above.** How the frozen rational coefficient vector `x` was selected.
  Not needed to verify the certificate, kept to document where the witness came from.

- [`extensions/openai-quasi-rh-audit-2026-10-07/`](extensions/openai-quasi-rh-audit-2026-10-07/)
  — audit of whether the quasi-RH material released by OpenAI on 2026-10-06 (pinned commit
  `adc7f1241b42`) supplies an estimate replacing the remaining-zero positivity hypothesis of
  Theorem 8. **Verdict: NO BRIDGE FOUND**, scoped to the audited source version and the
  transfer attempted: a restriction on where zeros may lie does not bound the aggregate
  window-weighted contribution of off-line zeros, and the per-zero bound available is uniform
  in the ordinate and so not summable as it stands — a limitation of that estimate, not a
  general obstruction. The audit also found Retraction 18, and verified the `δ²` expansion,
  `K(L,γ)` and its closed form against the TeX and the implementation.
  *Formalization status, exactly as established:* a comparator benchmark stub was inspected,
  the implementation declaration was located separately via `formalization.yaml`, and the full
  Lean proof and dependency chain were **NOT** independently checked — no build was run.
- [`extensions/comb-dips/`](extensions/comb-dips/) — measures the true last dip `t0(L)` of
  Zhu's prime-comb symbol against his worst-case threshold `T1(L) = 2π e^{A_L}`, and tests the
  dip set for arithmetic structure against preregistered null models.
  **Verdict: NULL.** `t0/T1 ∈ [0.449, 0.915]` over `L = 0.55 … 1.6`, so Zhu's threshold is
  sharp in practice and no certification shortcut exists on that route; no arithmetic
  structure survives the null-model tests.

## License

This repository is **dual-licensed**:

- **Code** — **MIT**, see [`LICENSE`](LICENSE). Covers `src/`, `scripts/`, the `Makefile`
  and the code under `extensions/`.
- **Paper and documentation** — **CC BY 4.0**, see
  [`LICENSE-CC-BY-4.0`](LICENSE-CC-BY-4.0). Covers `paper/`, `README.md`, `REPORT.md`,
  `PRIOR_ART.md`, `PAPER_OUTLINE.md`, `notes/`, `figures/`, `data/` and the prose, data and
  figures under `extensions/`.

CC BY 4.0 requires attribution: if you use the paper, the figures or the data, please cite
the paper, [doi:10.5281/zenodo.23250932](https://doi.org/10.5281/zenodo.23250932).

## How to cite

```bibtex
@techreport{DeJesus2026WeilWindow,
  author = {De Jes\'us, Elias},
  title  = {A $\delta^2$--$L^3$ Sensitivity Law for Finite-Window Weil Positivity,
            with a Function-Field Control of the Connes--van Suijlekom Pipeline},
  year   = {2026},
  type   = {Technical Note},
  doi    = {10.5281/zenodo.23250932},
  note   = {Zenodo. Code: \url{https://github.com/innerlightr-wq/weil-window-control}}
}
```

Cite the paper, DOI [10.5281/zenodo.23250932](https://doi.org/10.5281/zenodo.23250932) — see also [`CITATION.cff`](CITATION.cff), whose `preferred-citation` points there. [10.5281/zenodo.21109955](https://doi.org/10.5281/zenodo.21109955) is the **concept DOI** resolving to the latest version; use the version DOI above to cite this one specifically.

The code is distributed through this GitHub repository and has **no separate software DOI**; please cite the paper DOI for both. To pin the exact code, add the commit:

```
git clone https://github.com/innerlightr-wq/weil-window-control.git
cd weil-window-control && git checkout __COMMIT__
```

## Related work by the author

- *The Partition Potential and the Laguerre Tower*, Zenodo,
  [doi:10.5281/zenodo.21108392](https://doi.org/10.5281/zenodo.21108392)
- *Riemann Hypothesis: The Weil Enforcer in Audit Coordinates*, Zenodo,
  [doi:10.5281/zenodo.21115524](https://doi.org/10.5281/zenodo.21115524) — see
  [`notes/weil_enforcer_v3_note.md`](notes/weil_enforcer_v3_note.md) for a scope
  correction to that record arising from this work.
