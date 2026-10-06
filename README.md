# weil-window-control

Control experiments for the Connes / Connes–van Suijlekom "finite-prime Weil form"
strategy for RH, run in settings where the answer is known.

**Stages 1 and 2 are complete. Stage 3 (ζ) and Stage 4 (literature) are not started.**

Read `REPORT.md` for tiered findings, the deliverable answer, non-claims, and the
retraction log. `notes/stage1_derivation.md` has the exact derivation of the window
identity and the structural facts used throughout.

## Install

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt     # mpmath, numpy, matplotlib
# optional, much faster high-precision arithmetic:
pip install gmpy2
```

## Reproduce

Each command is self-contained and writes to `results/` or `figures/`.

```bash
# Stage 1 - function-field control (exact point counts -> Toeplitz Weil form)
python3 src/stage1.py --dps 50
#   -> results/stage1_spectra.csv, stage1_zeros.csv,
#      stage1_recorder_split.csv, stage1_curves.json
python3 src/figs_stage1.py
#   -> figures/stage1_lambda_min.png, stage1_zeros_circle.png, stage1_h4_convergence.png

# Stage 2 - planted off-circle configurations
python3 src/stage2.py --dps 50 --Rmax 20
#   -> results/stage2_spectra.csv, stage2_zeros.csv, stage2_rho_sweep.csv

# Stage 2 - exact closed form and detection thresholds
python3 src/detection_law.py
#   -> results/stage2_detection_law.csv

# Adversarial pass (attacks H1, H2, and the Stage 2 diagnosis)
python3 src/adversarial.py
#   -> results/stage2_masking.csv

python3 src/figs_stage2.py
#   -> figures/stage2_lambda_and_masking.png, stage2_zeros_stay_on_circle.png,
#      stage2_detection_law.png
```

Runtime: the whole pipeline is a few minutes on one core at `--dps 50`.
`stage1.py` accepts `--dps` (default 60); `stage2.py` accepts `--dps` and `--Rmax`.

## Layout

| path | contents |
|---|---|
| `src/ffield.py` | exact `F_{q^n}` arithmetic (polynomial mod irreducible) |
| `src/curves.py` | brute-force point counts, L-polynomial via Newton + functional equation, exact power sums |
| `src/weilform.py` | Toeplitz assembly, mpmath symmetric eigensolve, parity, polynomial zeros, inertia |
| `src/analysis.py` | one diagnostic row: precision-doubled sign test, float64 contrast, H1/H4 metrics |
| `src/exact_inertia.py` | **exact** inertia over Q by the Sylvester/Jacobi leading-minor rule |
| `src/planted.py` | functional-equation-closed fake spectra (on-line pair / off-line quartet / off-line real pair) |
| `src/detection_law.py` | closed form for `λ_min`, verification, detection thresholds |
| `src/adversarial.py` | attacks: random-Toeplitz H1, repeated-angle H2, masking, PSD survivors |
| `src/stage1.py`, `src/stage2.py` | drivers |
| `notes/stage1_derivation.md` | T1 derivations |

## Numerical policy

Every eigenvalue sign in this repository is produced by `mpmath` at a stated precision
with a precision-doubling convergence check, or certified exactly over Q. float64 is
computed only as a deliberate contrast: it returns spurious values of size `1e-16` to
`1e-15` where the true `λ_min` is exactly 0 (Zhu's warning, reproduced here — see
Finding F8 in `REPORT.md`).

## Non-claims

Nothing in this repository proves, advances, or provides evidence for the Riemann
hypothesis for ζ. Stage 1 runs where RH is a theorem (Weil, for curves over finite
fields); Stage 2 uses fabricated spectral data constructed by hand to be off-line.
