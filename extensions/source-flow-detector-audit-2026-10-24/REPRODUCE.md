# Reproducing this audit

Read-only with respect to everything else. This directory imports from
`extensions/davenport-heilbronn-2026-10-23` and
`extensions/channel-compatibility-exploration-2026-10-24` at commit `e2cf715` and modifies
neither; no manuscript, certificate, data file, licence or tag is touched, and nothing is
committed or published.

## Environment

```bash
.venv/bin/pip install 'mpmath==1.3.0'     # already present in the parent repository
```

Python **3.14.7** (CPython), **mpmath 1.3.0**. `python-flint` is **not** used: nothing here
is offered as a certified enclosure, so the Arb path is deliberately not invoked. The one
certified number quoted (`λ_min(Q_0^{DH}) = −0.679…`, ball radius `1.39×10⁻¹¹⁶`) is read
from the parent extension's existing certificate, not recomputed. No `numpy`, no
`float64` on any path to a reported number.

## Commands

```bash
E=extensions/source-flow-detector-audit-2026-10-24

.venv/bin/python $E/src/source_flow_audit.py --dps 40                      # ~12 min
.venv/bin/python $E/src/lambda_lower_bound.py --arcs 3000 --prec 400       # ~1 min
.venv/bin/python $E/tests/test_identities.py
tests/test_lambda_certificate.py  regression fixture for the 6c' certificate                               # seconds
.venv/bin/python $E/tests/test_lambda_certificate.py                       # seconds
```

`lambda_lower_bound.py` is the only script here that uses **python-flint / Arb**, because
it is the only one producing a certified (T2-C) number. It is self-contained: it imports
nothing from the other extensions. Its partition size is chosen from a measurement —
the Hurwitz-zeta sum amplifies an input ball radius by about `106` at this height, so
`--arcs 60` fails (relative enclosure radius `3.30`) and `--arcs 3000` succeeds (`0.065`).
Lowering `--arcs` below roughly `600` will make the certification fail rather than
silently degrade.

`source_flow_audit.py` writes `data/validation.csv` (one row per checked quantity, with a
note recording what each check is for) and `data/audit_log.txt` (the full transcript,
block by block). `--dps` raises the working precision; 40 is sufficient for every reported
number, and §6 of `REPORT.md` says why.

## What each block does

| block | content |
|---|---|
| **A** | coordinate dictionary; verifies `H_0(z) = (1/8)ξ(1/2+iz/2)` and the evenness of `Φ` |
| **B** | verifies `∂_τH = −∂_z²H` by finite differences on the integral representation, both kernels |
| **C** | derives and checks the Davenport–Heilbronn kernel `Φ_F`: transform, evenness (10), decay (11), positivity on `x ≥ 0` |
| **D** | de Bruijn Theorem 7 applied to each object; reproduces `τ ≥ 1/2` for ζ and gives `τ ≥ 1.6025374836` for DH |
| **E** | the detector bridge `Ψ ↔ Σ_ρ H(ρ−½)`, and the reality ⟹ PSD lemma |
| **F** | the feasibility gate on a complete differentiable family (result: **not justified**) |
| **K** | arithmetic content of the evolution at `τ = 0`: the forced `Λ*Λ` convolution and why it is unusable |
| **J** | pointer to the `η` / eigenvalue identities, which live in the test file |
| **G** | adversarial control 1: the toy family `z² + a² − 2τ` |
| **H** | adversarial control 2: manufactured dynamics from an assembled `B` |
| **I** | the obstruction: the DH counterexample, with a tier on each step |
| **L** | positive kernel ⟹ PSD Fourier-correlation matrix, a different form |
| **M** | the two-sided bracket on `Λ_DH` (analytic lower bound) |

`tests/test_identities.py` checks `η′ = 2(u′)^T x − x^T B′ x`, `Qx = 0` at `η = 1`,
`λ_0′ = −η′/‖x‖²`, and the toy family's transition — all against a **synthetic, labelled**
perturbation, because block F did not justify a source-derived one.

## Precision discipline

Every decimal constant goes through `dh.mpf_hi`, which parses at a fixed 220 digits
independent of the ambient precision (see the parent extension's `RETRACTIONS.md` R6 for
why that matters). Finite-difference checks use a step of `10⁻¹⁰` and therefore verify
their identities to about `10⁻²²`, which is what the report states — not to the working
precision.

## Files

| file | produced by |
|---|---|
| `data/validation.csv`, `data/audit_log.txt` | `src/source_flow_audit.py` |
| `data/lambda_bracket.csv`, `data/lambda_bracket.txt` | `src/lambda_lower_bound.py` (§6c′) |
| `data/run.log` | the committed run of the audit script |

```
REPORT.md            verdict first, then the audit in the brief's order
SOURCE_AUDIT.md      what was read, in which version, and which statement came from where
REPRODUCE.md         this file
src/source_flow_audit.py      the audit proper (mpmath)
src/lambda_lower_bound.py     the certified Lambda_DH lower bound (Arb)
tests/test_identities.py
```
