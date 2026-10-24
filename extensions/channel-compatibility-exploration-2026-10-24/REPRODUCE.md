# Reproducing this exploration

This directory is **read-only with respect to everything else**. It imports from
`extensions/davenport-heilbronn-2026-10-23` at commit `e2cf715` and modifies nothing
there; no manuscript, certificate, data file, root README, tag or licence is touched, and
nothing is committed or published.

## Environment

From the repository root, with the parent repository's pinned environment:

```bash
.venv/bin/pip install 'mpmath==1.3.0'      # already present
```

Python **3.14.7** (CPython), **mpmath 1.3.0**. `python-flint` is *not* used here: no
quantity in this exploration is offered as a certified enclosure, so the Arb path is
deliberately not invoked. No `numpy`, no `float64` anywhere on a path to a reported
number.

## Provenance of everything reused

| imported | from | used for |
|---|---|---|
| `dh.build_matrix`, `dh.build_parts` | `davenport-heilbronn-2026-10-23/src/dh.py` | the complete arithmetic-side matrix `Q` |
| `dh.G_at` | same | the transforms that give the channel vectors |
| `dh.mpf_hi` | same | parsing every decimal constant at 220 digits (see that extension's `RETRACTIONS.md` R6) |
| `phaseD.index_set`, `LEVELS` | `.../src/phaseD.py` | the declared basis levels and index sets |
| `certificates/witness.json`, `witness_smallest_L.json` | same | the frozen witness coefficients, re-parsed, for the benchmark |
| `results/census.json` | same | the four off-line quartet locations (**numerical inputs, T2, not certified**) |

## Commands

```bash
E=extensions/channel-compatibility-exploration-2026-10-24

# everything, in order, with the loop in a script file:
./$E/src/run_all.sh

# or individually:
.venv/bin/python $E/src/explore_channels.py --dps 80 --windows 1.6,1.8,2.0 \
                                            --levels 2,3 --witness
.venv/bin/python $E/src/explore_channels.py --dps 80,140,220 --windows 1.8 --levels 3 \
                                            --out channel_ladder_18.csv
.venv/bin/python $E/src/explore_channels.py --dps 80,140 --windows '' --levels '' \
                                            --decomp 2.0:2
.venv/bin/python $E/src/explore_channels.py --dps 80,140 --windows '' --levels '' \
                                            --enrichment 1.8
.venv/bin/python $E/src/explore_channels.py --dps 80 --windows '' --levels '' \
                                            --twochannel 2.0:2
.venv/bin/python $E/src/explore_channels.py --dps 80 --windows '' --levels '' \
                                            --invariance 2.0:2
```

Approximate runtimes on this machine (single core, no gmpy2): the `N = 37` configurations
take a few minutes each at 80 digits, the `N = 57` configurations substantially longer,
and the 220-digit rung of the ladder dominates the total. `run_all.sh` is restartable:
each invocation is independent and writes its own log and CSV.

## Outputs

| file | produced by |
|---|---|
| `data/channel_results.csv` | the main grid (`--windows --levels`) |
| `data/channel_ladder_18.csv`, `data/channel_ladder_20.csv` | the precision ladders |
| `data/enrichment_results.csv` | `--enrichment` |
| `data/twochannel_results.csv` | `--twochannel` |
| `data/invariance_results.csv` | `--invariance` |
| `data/*.log` | the console transcript of each run |

## What the columns mean

`lmin_Q`, `lmin_B`, `lmin_R` are the least eigenvalues of the complete arithmetic form,
the candidate background `B = Q + uu^T`, and the remaining-zero form `R = Q - Q_O`.
`n_neg_Q` is the number of negative eigenvalues of `Q`; a rank-one channel can account
for at most one, so `n_neg_Q >= 2` would mean `B` cannot be positive definite.
`eta` is `u^T B^{-1} u` computed by a linear solve, never by forming `B^{-1}`.
`eta_from_det` is `1 - det(Q)/det(B)`, and `det_identity_rel_residual` is the relative
disagreement between the two — the matrix-determinant lemma, checked rather than assumed.
`criterion_agrees` and `inertia_agrees` are the validation gates against the direct
eigenvalue computation; they are gates, not findings.

## Precision discipline

Every decimal constant — window lengths, the quartet locations, the frozen witness
coefficients — is parsed by `dh.mpf_hi` at a fixed 220 digits, independent of the ambient
working precision, so that no value previously rounded at a lower precision is reused.
The ladder (80 / 140 / 220 digits) tests **stability**, which is not an interval
certificate and is nowhere described as one.
