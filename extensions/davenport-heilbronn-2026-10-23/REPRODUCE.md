# Reproducing every number in this extension

Everything runs from the **repository root** with the parent repository's pinned
environment plus one addition.

```bash
cd <repository root>
python3 -m venv .venv                       # if not already present
.venv/bin/pip install 'mpmath==1.3.0'
.venv/bin/pip install python-flint           # for src/certify.py only
```

| | version used for the committed numbers |
|---|---|
| Python | **3.14.7** (CPython) |
| mpmath | **1.3.0** |
| python-flint | **0.9.0** (Arb; used only by `src/certify.py`) |
| numpy / scipy / matplotlib | **not used** |

All conventions, precisions, windows, basis sizes and cutoffs are in
`src/config.py`. Nothing numeric is hard-coded twice. No randomness is used; the
`SEED` entry in the config exists only so the claim "no randomness" is auditable.

Every script is restartable: each writes its own `results/*.txt` and `data/*.csv` and
reads nothing from another script except through those files.

## Dependency order

```
src/phaseA.py        # the function, kappa, functional equation, the off-line zero
src/phaseB.py        # b_n, cutoff, transforms, the arithmetic form, scalar cross-route
tests/test_orbits.py # orbit conventions (independent of F; can run any time)
src/census.py        # -> data/dh_zeros.csv     (needed by phaseC, phaseF)
src/phaseC.py        # the gate                 (needs data/dh_zeros.csv)
src/phaseD.py        # the window scan          -> certificates/witness.json
src/phaseE.py        # the remaining-zero form  (needs certificates/witness.json)
src/phaseF.py        # masking                  (needs witness.json + dh_zeros.csv)
src/phaseG.py        # sensitivity              (needs witness.json)
src/phaseH.py        # arithmetic-side split    (needs witness*.json)
src/gate_signflip.py # Gate A term-level sign-flip control
src/certify.py       # Arb ball-arithmetic certificate of the frozen witness
tests/regression.py  # parent-repository regressions, before and after
```

## Commands

```bash
E=extensions/davenport-heilbronn-2026-10-23

.venv/bin/python $E/src/phaseA.py                  #  ~3 min
.venv/bin/python $E/src/phaseB.py                  #  ~6 min
.venv/bin/python $E/tests/test_orbits.py           #  ~2 min
.venv/bin/python $E/src/census.py --Tmax 200       # ~25 min
.venv/bin/python $E/src/phaseC.py                  # ~20 min
.venv/bin/python $E/src/phaseD.py                  # ~45 min
.venv/bin/python $E/src/phaseE.py                  # ~15 min
.venv/bin/python $E/src/phaseF.py                  #  ~5 min
.venv/bin/python $E/src/phaseG.py                  # ~10 min
.venv/bin/python $E/src/phaseH.py                  # ~25 min
.venv/bin/python $E/src/gate_signflip.py           #  ~3 min
.venv/bin/python $E/src/certify.py 400             #  ~5 min
.venv/bin/python $E/tests/regression.py            #  ~5 min
```

or, in order, with the first failure propagated:

```bash
./extensions/davenport-heilbronn-2026-10-23/run_suite.sh
```

Each script exits non-zero when its own gate fails. `run_suite.sh` stops at the first
failure rather than continuing past a failed gate, which is the behaviour the brief
requires of the Phase C classification.

## Which script produces which file

| file | produced by |
|---|---|
| `data/dh_coefficients.csv` | `src/phaseB.py` |
| `data/dh_offline_orbit.csv` | `src/phaseA.py` |
| `data/dh_zeros.csv` | `src/census.py` |
| `data/gate_checks.csv` | `src/phaseC.py` |
| `data/window_scan.csv` | `src/phaseD.py` |
| `data/remaining_form.csv` | `src/phaseE.py` |
| `data/masking.csv` | `src/phaseF.py` |
| `data/sensitivity.csv` | `src/phaseG.py` |
| `data/arith_split.csv` | `src/phaseH.py` |
| `data/gate_signflip.csv` | `src/gate_signflip.py` |
| `certificates/witness.json` | `src/phaseD.py` |
| `certificates/ball_certificate.txt` | `src/certify.py` |
| `results/*.txt`, `results/*.json` | the correspondingly named script |
| `logs/*.log` | the background runs that produced the committed results |

## Re-deriving the frozen witness independently

`certificates/witness.json` holds the window, sector, index set, basis frequencies and
rationalised coefficients as **decimal strings**. To recompute the witness value without
using any of this extension's eigenvalue code:

```bash
.venv/bin/python $E/src/certify.py 400
```

which re-implements the whole form against Arb and prints a ball enclosing
`Q_F(g_*)`, `Q_{O_1}(g_*)` and `Q_rest(g_*)`.
