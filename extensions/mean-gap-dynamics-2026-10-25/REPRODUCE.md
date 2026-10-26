# REPRODUCE.md

```
Python 3.12.2 (CPython)      mpmath 1.3.0
python-flint 0.9.0 (Arb)     -- needed only by src/certify.py
numpy 1.26.4                 -- installed, NOT used on any reported path
```

`python-flint` is the version the parent Davenport–Heilbronn extension declares; installing
it is additive and changes no pinned behaviour. Everything except `src/certify.py` runs on
mpmath alone. No randomness is used anywhere; every script is independent and communicates
only through `data/*.json` and `results/*.txt`.

Run from this directory.

| order | command | what it produces | runtime* |
|---|---|---|---|
| 0 | `python3 src/check_dictionary.py` | `results/stage0_dictionary.txt`, `data/stage0_dictionary.json` — the normalisation audit, **run this first** | ~6 min |
| 1 | `python3 src/means.py` | `results/part2_means.txt`, `data/part2_means.json` — exact over `Q` | < 1 s |
| 2 | `python3 src/toy.py` | `results/part4_toy.txt`, `data/part4_toy.json` — the exact counterexample | ~3 s |
| 3 | `python3 src/dh_flow.py` | `results/part4_dh.txt`, `data/part4_dh.json` — `D'(0)` at the four DH orbits, two routes | ~4 min |
| 4 | `python3 src/certify.py` | `results/certificates.txt`, `data/certificates.json` — Arb at 400 bits | ~40 s |
| 5 | `python3 src/alt_evolution.py` | `results/part5_alt.txt`, `data/part5_alt.json` | ~20 s |
| 6 | `python3 src/dh_interaction.py` | `results/part3_interaction.txt`, `data/part3_interaction.json` — the remainder accounting | ~25 min |
| 7 | `python3 -u src/zeta_flow.py` | `results/part4_zeta.txt`, `data/part4_zeta.json` — `zeta` branches and collisions | ~25 min |
| 8 | `python3 -u src/zeta_lehmer.py` | `results/part4_zeta_lehmer.txt`, `data/part4_zeta_lehmer.json` — why the tracked pair is not extremal | ~3 min |
| 9 | `python3 -u src/separation_bound.py` | `results/part3_separation.txt`, `data/part3_separation.json` — the necessary ordinate-separation condition | < 5 s |
| A1 | `python3 -u src/audit1_claims.py` | `results/audit1.txt`, `data/audit1.json` — AUDIT Task 1: the three outstanding claims | ~10 s |
| A2 | `python3 -u src/audit2_constraints.py` | `results/audit2.txt`, `data/audit2.json` — AUDIT Tasks 2 and 3 | < 5 s |
| A4 | `python3 -u src/audit4_lemma.py` | `results/audit4.txt`, `data/audit4.json` — AUDIT Task 4 (needs step 6's `data/dh_online_ordinates.txt`) | ~10 s |
| A5 | `python3 -u src/dh_fixed_tau.py` | `results/dh_fixed_tau_nav.txt` — navigation: the branch at fixed `tau`, numerical winding numbers | ~8 min |
| A6a | `python3 -u src/check_taylor_vs_integral.py` | `results/taylor_vs_integral.txt`, `data/taylor_vs_integral.json` — the representation the certificate rests on, checked at `tau != 0` against the defining integral, plus a wrong-constant control | ~12 min |
| A6 | `python3 -u src/certify_fixed_tau.py` | `results/certificate_fixed_tau.txt`, `data/certificate_fixed_tau.json` — **the certificate**: `Lambda_DH > 0.16`, plus a negative control | ~4 min |
| — | `python3 tests/regression.py` | the regression suite; exits non-zero on any failure | ~5 min |
| — | `python3 tests/verify_numbers.py` | checks every headline number printed in the Markdown against `data/`; fails on disagreement, absence, or **overstated precision** | < 5 s |

\* one core of an Apple Silicon laptop, without `gmpy2`. Installing `gmpy2` speeds the
high-precision scripts up substantially and changes no result.

`DH_SCAN_HI` (default `260`) sets the critical-line scan height in step 6; lowering it to
`120` cuts that step to about 5 minutes at the cost of a larger measured remainder.

Every script exits `0` on success and **non-zero** if any of its own gates fails; the gates
assert the inequalities and residuals themselves, they do not print a stored label.
`tests/regression.py` carries its expected values **inside the test file**, so editing a file
in `data/` or `results/` makes it fail rather than redefining what passing means, and its
last block is a **negative control** that confirms a perturbed digit is actually rejected.

## Numerical policy

* Quadrature of `H_tau` uses two nodes per oscillation of `cos(zu)`. Without this the
  integrand at `|z| ~ 171` returns about ten correct digits silently; the first version of
  the dictionary audit **failed** for exactly this reason (`DERIVATION.md` §0.3).
* All `z`-derivatives of `H_tau` are taken **analytically** (differentiating `cos(zu)` under
  the integral sign); no finite differences enter any reported value. The one place finite
  differences appear is the `tau`-derivative in the PDE check, with the step matched to the
  working precision and reported as such.
* Collision gates use **scale-free** residuals (`|H|/|H''|`, `|H'|/|H''|`), because `|H|` is
  of order `e^{-pi gamma/4}` and an absolute threshold on it is meaningless. The first
  version of that gate used an absolute threshold and reported three genuine double zeros as
  `UNRESOLVED`.
* One-sided certified bounds are printed as plain decimals from the ball's own endpoint,
  rounded **outward** — never the midpoint rounded to nearest.
* The toy model and the mean-gap algebra are exact over `Q` (`fractions.Fraction` and a
  small exact `Q(i)` class), so they need no interval arithmetic and are stronger as they
  stand.
