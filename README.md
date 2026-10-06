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

📄 **Paper:** [`paper/weil_window_control.pdf`](paper/weil_window_control.pdf) ·
🔍 **Prior art, verified in full text:** [`PRIOR_ART.md`](PRIOR_ART.md) ·
🧾 **Full lab record incl. 16 retractions:** [`REPORT.md`](REPORT.md)

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
5. `λ*(0.8) ≤ 2.2702e−17` is a **variational upper bound**, not a certified enclosure.
6. The ζ-side replication **reproduces published computations** (Groskin) at smaller scale.
7. **Block F is claimed only for the function field** — on ζ windows the indefiniteness of
   `A` and `M` depends on where `−log π` is placed (Retraction 15).
8. **Not attempted:** Davenport–Heilbronn; any Lean formalisation.

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
| **ζ bound `λ_min ≥ λ*(L) − c(L)δ² − O(δ⁴)`, `c(L) = 8L³(⅓+1/√5)`** | **T1** | `src/proofs_check.py` | `data/stage3_sharp_bound.csv` |
| Sharper conditional bound, constant `8L³/3` | T1 (conditional) | `src/proofs_check.py` | `data/stage3_sharp_bound.csv` |
| `C/(4K(L,γ)) → 0.963`: the exact constant nearly saturated | T2 | `src/task_final_checks.py` | `data/stage3_exact_constant_K.csv` |
| Gate A: geometric side vs `Σ_ρ|F(γ_ρ)|²` to 1.6e−7 | T2 | `src/stage3_converge.py` | `data/stage3_gateA_explicit_formula.csv` |
| Gate B: `λ*(0.8) ≤ 2.2702e−17`, two independent bases | T2 | `src/stage3_converge.py` | `data/stage3_gateB_lambda_min.csv` |
| Connes `[1,13]` window: `λ_min = 8.977e−52` at `N=36` | T2 | `src/stage3_connes.py` | `data/stage3b_connes_zeros_N36.csv` |
| Detection onset tracks `L_pred = ½log(γ/2π)`, moves earlier with basis | T2 | `src/stage3e_certify.py` | `data/stage3e_certified_onset.csv` |
| Recorder split; Block F structural in the function field only | T2 | `src/task2_inertia_reconcile.py` | `data/stage3d_inertia_reconciled.csv` |

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
```

## Licensing

This repository is **dual-licensed**:

- **Code** (`src/`, `scripts/`, `Makefile`) — **MIT**, see [`LICENSE`](LICENSE).
- **Paper and documentation** (`paper/`, `REPORT.md`, `PRIOR_ART.md`, `PAPER_OUTLINE.md`,
  `notes/`, `figures/`, `data/`) — **CC BY 4.0**, see
  [`LICENSE-CC-BY-4.0`](LICENSE-CC-BY-4.0).

## How to cite

```bibtex
@techreport{DeJesus2026WeilWindow,
  author = {De Jes\'us, Elias},
  title  = {A $\delta^2$--$L^3$ Sensitivity Law for Finite-Window Weil Positivity,
            with a Function-Field Control of the Connes--van Suijlekom Pipeline},
  year   = {2026},
  type   = {Technical Note},
  note   = {Zenodo. Code: \url{https://github.com/innerlightr-wq/weil-window-control}}
}
```

See [`CITATION.cff`](CITATION.cff). A DOI will be minted on the Zenodo release.

## Related work by the author

- *The Partition Potential and the Laguerre Tower*, Zenodo,
  [doi:10.5281/zenodo.21108392](https://doi.org/10.5281/zenodo.21108392)
- *Riemann Hypothesis: The Weil Enforcer in Audit Coordinates*, Zenodo,
  [doi:10.5281/zenodo.21115524](https://doi.org/10.5281/zenodo.21115524) — see
  [`notes/weil_enforcer_v3_note.md`](notes/weil_enforcer_v3_note.md) for a scope
  correction to that record arising from this work.
