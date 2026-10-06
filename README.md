# weil-window-control

Control experiments for the Connes / Connes–van Suijlekom "finite-prime Weil form"
strategy for RH, run in settings where the answer is known.

**Stages 1–4 complete, plus a Round-2 pass** that corrected the Stage 3e height sweep,
derived the `delta^2` and `L^3` laws, and settled novelty. Stage 3f (Davenport–Heilbronn)
was deliberately not attempted; the reason is in `REPORT.md`. See `PAPER_OUTLINE.md` for
the section map with every result tagged T1/T2/T3 and pointed at its script and data.

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

# Task 1e - T_R as a Gram matrix in the trace form on H^1 (Cor. 14 at general R, g)
python3 src/gram_trace.py
#   -> results/stage1e_gram_trace.csv

# Items 1-2 - kernel-as-ideal, and function-level re-scoring of step (iv)
python3 src/function_level.py --dps 50 --Rmax 14
#   -> results/stage2_function_level.csv
python3 src/figs_function_level.py
#   -> figures/stage2_iv_rescored.png

# Stage 2 - exact closed form and detection thresholds
python3 src/detection_law.py
#   -> results/stage2_detection_law.csv

# Adversarial pass (attacks H1, H2, and the Stage 2 diagnosis)
python3 src/adversarial.py
#   -> results/stage2_masking.csv

python3 src/figs_stage2.py
#   -> figures/stage2_lambda_and_masking.png, stage2_zeros_stay_on_circle.png,
#      stage2_detection_law.png

# ---- Stage 3: the zeta window Weil form on [-L, L] ----
# Gate A (assembly vs the zeros side) + Gate B (lambda*(0.8) from two bases).
# ~45 min; the zetazero calls dominate the first few minutes.
python3 src/stage3_converge.py --N 4,6,8,10,12,14,16,18,20,24 --dps 50 --zeros 800
#   -> results/stage3_gateA_explicit_formula.csv, stage3_gateB_lambda_min.csv

# 3b Connes's [1,13] experiment at high precision (~40 min)
python3 -c "import sys;sys.path.insert(0,'src');import mpmath as mp;mp.mp.dps=110;\
from stage3_connes import part3b; r=[]; part3b(110,24,50,r)"
# 3c even/odd gap and 3d recorder split (~15 min)
python3 src/stage3_connes.py --N 12 --dps 40 --parts c,d
#   -> results/stage3b_connes_zeros.csv, stage3c_even_odd.csv, stage3d_recorder_split.csv

# 3e teeth on the zeta side (zeros-side diagnostic, NOT a certificate)
python3 src/stage3e_teeth.py --K 400 --N 14 --dps 80
#   -> results/stage3e_teeth.csv
# 3e adversarial extension: does the planted HEIGHT matter? (it does - see REPORT.md)
python3 src/stage3e_height.py
#   -> results/stage3e_height_sweep.csv

python3 src/figs_stage3.py
#   -> figures/stage3_gates.png, stage3_connes_profile.png, stage3_teeth.png

# ---- Round 2 ----
# Task 1: certified detection onset vs basis size (geometric on-line side, zero-moving
# perturbation planted at actual zeta ordinates).  Hours at these precisions.
python3 src/stage3e_certify.py --ns 1,5,10,30,80 --deltas 0.1 \
        --factors 0.70,1.00,1.30,1.70,2.20
python3 src/stage3e_certify.py --ns 30,80 --factors 0.85,1.00,1.15 --dps 220 \
        --out stage3e_onset_highprec.csv
python3 src/task1e_C_vs_height.py        # Task 1(e): C(L, gamma_*) vs height
python3 src/delta2_law.py                # Tasks 2, 3: delta^2 expansion, L^3 bound
python3 src/task34_extra.py              # Task 3b: C vs CS;  Task 4b: inertia vs N
```

Stage 3 runs take tens of minutes each at these precisions. Installing `gmpy2` speeds
mpmath up substantially and is worth it here.

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
| `src/gram_trace.py` | `T_R` = Gram of `{F^k/q^{k/2}}` in `Tr(φψ†)`; where positivity of the polarization enters |
| `src/function_level.py` | kernel-as-ideal, kernel parity split, function-level distances M1/M1′/M2 |
| `src/adversarial.py` | attacks: random-Toeplitz H1, repeated-angle H2, masking, PSD survivors |
| `src/zeta_window.py` | the zeta window form: three bases, closed-form `F`, `C`, and archimedean term |
| `src/stage3_converge.py` | Gate A (explicit-formula validation) and Gate B (lambda* convergence) |
| `src/stage3_connes.py` | 3b Connes [1,13]; 3c even/odd gap; 3d recorder split |
| `src/stage3e_teeth.py` | 3e planted off-line quartet on the zeta side |
| `src/stage3e_height.py` | Round-1 3e height sweep — **superseded**, see retraction 11 |
| `src/stage3e_certify.py` | Round-2 3e: LDL witnesses, geometric on-line side, zero-moving perturbation |
| `src/ldl.py` | LDL^T inertia of every leading block + explicit negative-direction witness |
| `src/iv_certify.py` | mpmath-interval certification (zeros-side; see Task 1c for its limit) |
| `src/geom_at_vector.py` | geometric side at ONE trial vector, no N x N matrix |
| `src/delta2_law.py`, `src/task34_extra.py`, `src/task1e_C_vs_height.py` | Tasks 2–4 |
| `src/stage1.py`, `src/stage2.py` | drivers |
| `notes/stage1_derivation.md` | Stage 1 T1 derivations |
| `notes/stage3_assembly.md` | Stage 3 T1 derivations (bases, Parseval, archimedean closed form) |

## Numerical policy

Every eigenvalue sign in this repository is produced by `mpmath` at a stated precision
with a precision-doubling convergence check, or certified exactly over Q. float64 is
computed only as a deliberate contrast: it returns spurious values of size `1e-16` to
`1e-15` where the true `λ_min` is exactly 0 (Zhu's warning, reproduced here — see
Finding F8 in `REPORT.md`).

## Non-claims

Nothing in this repository proves, advances, or provides evidence for the Riemann
hypothesis for ζ. Stage 1 runs where RH is a theorem (Weil, for curves over finite
fields); Stage 2 uses fabricated spectral data constructed by hand to be off-line;
Stage 3's Gate A *assumes* RH in order to check the assembly, Stage 3e is built from zero
data and is a diagnostic only, and Stage 3b reproduces at smaller scale a computation
already published (Groskin, arXiv:2605.20224). No novelty is claimed anywhere; see the
Stage 4 literature check and the full non-claims list in `REPORT.md`.
