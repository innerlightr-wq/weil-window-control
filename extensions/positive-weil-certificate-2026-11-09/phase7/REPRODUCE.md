# Reproducing the two-window package

Everything below runs inside this subtree. No `$HOME` path is read and no
network access is used. The release gate was last run from a clean extraction of
`note/weil-two-window-certificate-2.0.0.tar.gz` with `HOME=/nonexistent`.

## Dependencies

| component | version used |
|---|---|
| Python | 3.14.7 |
| `python-flint` | 0.9.0 (FLINT / Arb) |
| `mpmath` | 1.3.0 |
| TeX | any `pdflatex` with `amsmath`, `booktabs`, `hyperref`, `microtype` |

`numpy`, `gmpy2` and `sympy` are **not** used. `mpmath` only *locates* candidate
vectors, which are then frozen as exact rationals and certified in Arb; it
appears in no pass/fail decision of the gate.

## The release gate — run this first

```sh
cd extensions/positive-weil-certificate-2026-11-09
python phase7/acceptance.py          # 61 checks, ~210 s, nonzero exit on failure
```

It verifies the integrated candidate from the packaged files: it applies and
reverts nothing, reads no value held in memory from the run that produced it,
and can be run twice. Stages are documented in the script's docstring. It
reports **61 checks, 0 failed** on this package, and the repository suite inside
it reports 176 passing checks.

From an extraction of the archive:

```sh
mkdir /tmp/stage && cd /tmp/stage
tar xzf .../weil-two-window-certificate-2.0.0.tar.gz
cd extensions/positive-weil-certificate-2026-11-09
python phase7/acceptance.py
```

## What is replayed, and what is not

**Replayed** (all of it, by the gate): every Schur sign; both improved
unconditioned `T = 400` positive-definiteness checks; the four conditioned
rank-`K` bumps, rebuilt from the stored rational vectors; both second-window
blocks; both closed-form upper quotients, re-certified from their stored exact
rationals; the four-dimensional count certificate including its Gram matrix; the
directed-decimal checks; the two-window ratio; and the complement, coupling and
Bernstein-remainder provenance fields.

**Not replayed** — these are *cold rebuilds* and are not part of the gate:

| assembly | script | measured |
|---|---|---|
| `T = 200`, `L = 4/5` | `phase2/src/run_block.py` | 1343 s |
| `T = 400`, `L = 4/5` | `phase4/src/run_block400.py` | 3532 s |
| `T = 400`, `L~` | `phase7/src/run_block_lam5.py` | 3678 s |

The gate consumes their frozen output (`phase4/certs/block400_{even,odd}.json`,
`phase7/certs/block400_lam5_{even,odd}.json`). Checksums establish that the
files are the ones that were produced; they are **not** evidence of numerical
correctness, which is what the replay is for.

## Regenerating individual pieces

```sh
EXT=extensions/positive-weil-certificate-2026-11-09
PP=$EXT/phase2/src:$EXT/phase4/src:$EXT/phase5/src:$EXT/src:src

# the repository suite (phases 1-6), 176 checks
cd $EXT
PYTHONPATH=$PP python verify.py
PYTHONPATH=$PP python phase2/verify_phase2.py
PYTHONPATH=$PP python phase2/verify_phase3.py
PYTHONPATH=$PP python phase3/verify_audit.py
PYTHONPATH=$PP python phase4/verify_t400.py
PYTHONPATH=$PP python phase5/test_prec.py
PYTHONPATH=$PP python note/consistency_check.py
PYTHONPATH=$PP python phase5/src/upper_highN.py --verify

# phase 7 pieces (each writes into phase7/certs/)
PYTHONPATH=$PP python phase7/src/push_shift.py       # push the frozen T=400 shifts
PYTHONPATH=$PP python phase7/src/finalize.py         # outward endpoints + ratios
PYTHONPATH=$PP python phase7/src/margin_audit.py     # the 8 Schur margins
PYTHONPATH=$PP python phase7/src/condition.py        # even conditioned coercivity
PYTHONPATH=$PP python phase7/src/condition_odd.py    # odd conditioned coercivity
PYTHONPATH=$PP python phase7/src/counts_even.py 160  # the even count certificate
PYTHONPATH=$PP python phase7/src/finalize_lam5.py    # second-window certification
                   python phase7/src/upper_lam5.py 160   # second-window upper bound
                   python phase7/src/window_ratio.py     # the two-window ratio
                   python phase7/src/dictionary_check.py # the normalization dictionary
                   python phase7/src/prolate_rungs.py    # exact prolate rungs
```

`phase7/src/lowvecs.py` and `lowvecs_odd.py` are the *numerical* step that
selects the conditioning subspaces; their output is already frozen inside the
certificates, so the conditioned certificates replay without them.

`phase7/src/make_patch.py` is provenance only and refuses to run: it generated
the documentation diff against the pre-integration baseline, which is now
applied. Its section sources are in `phase7/patch_sections/`.

## Rounding discipline

Every lower bound is consumed and published through the **lower** endpoint of its
enclosure, every upper bound through the **upper** endpoint. No ball midpoint
and no printed decimal is used as a bound. Published `mu` values are short exact
decimals at `a(1 - 10^-12)`, which keeps about `10^30` times the Schur margin
`kappa^2/d` and survives serialization; the gate checks the direction of every
short decimal and explicitly rejects the two that were wrong in earlier drafts.
