# Reproducing this extension

Everything here is `mpmath` at a stated precision. No `float64` is on any path whose sign
matters; where a float appears it is a loop index or a printed label.

## Environment

```bash
cd <repo root>
python3 -m venv .venv && source .venv/bin/activate
pip install 'mpmath==1.3.0'          # the version the parent repository pins
```

`numpy` and `matplotlib` are **not** needed for this extension. `gmpy2` is optional and
changes no result, only the runtime.

Measured on this machine: Linux x86-64, CPython 3.14.7, mpmath 1.3.0, no gmpy2, 8 cores
(each phase is single-threaded; the phases below were run concurrently).

> **Note on the pinned version.** The parent repository's committed data was produced with
> mpmath 1.3.0 under CPython 3.12.2. This extension uses the same mpmath but a newer
> CPython. The ζ control values reproduced inside every phase here are the check that this
> makes no difference; mpmath 1.4.x was **not** used for any committed number.

## Order of operations

Phases 1 and 2 are **gates**: if they fail, nothing downstream is meaningful. Run them
first, in order.

```bash
E=extensions/dirichlet-windows-2026-10-23

# --- Phase 1: the form, and six external-reference self-tests        ~2 min
.venv/bin/python $E/src/selftest.py                 # exits non-zero on any failure
.venv/bin/python $E/src/selftest_planting.py        # planting-block regressions

# --- Phase 2: zeros, then GATE A                                 ~25 + ~20 min
.venv/bin/python $E/src/zeros.py --Tmax 275         # writes data/zeros_<chi>.txt
.venv/bin/python $E/src/lmfdb_check.py              # cross-check (needs network)
.venv/bin/python $E/src/phase2_gateA.py --Ls 0.5,0.8,1.0 --Ks 20,40,80,160 \
    --sectors even_d,even,odd --zeta-zeros 160
#   STOP HERE if Gate A fails for any character.

# --- Phase 3: GATE B and the conductor/resolution test               ~40 + ~30 min
.venv/bin/python $E/src/phase3_gateB.py --Ls 0.5,0.8,1.0,1.3 --Ns 4,6,8,10,12,14,16 --dps 50
.venv/bin/python $E/src/phase3b_lpred.py --delta 0.375 --N 16

# --- Phase 4: the sensitivity law                                    ~40 min
.venv/bin/python $E/src/phase4_sensitivity.py
.venv/bin/python $E/src/phase4b_sector_extremal.py  # exact, seconds

# --- Phase 5: the Siegel corner                                      ~40 min + ~1 h
.venv/bin/python $E/src/phase5_siegel.py --Ls 0.8,1.4 --Ns 8,12,16 --sectors even_d,odd
.venv/bin/python $E/src/phase5_onset.py --sector odd --N 12
```

Console transcripts of the runs that produced the committed numbers are in `results/`;
machine-readable results are in `data/` (CSV) and `results/` (JSON).

## What each gate asserts

| script | asserts | fails how |
|---|---|---|
| `selftest.py` | characters, Gauss sums, `L(1,χ_d)` vs the class-number formula, `L(3,χ₋₄) = π³/32`, `L(0,χ_d)` as exact rationals, a bounded-tail Dirichlet series, `Λ(s) = Λ(1−s)` off the line, `Im Z_χ ≈ 0`, and `arch_closed_a(a=0) ≡` the parent `arch_closed` | non-zero exit, named check |
| `selftest_planting.py` | `F_at_imag` vs the parent `F_at_half` and vs quadrature; `real_block(δ=1/2) ≡` the parent **pole-term matrix** in every sector; rank-one and definite-sign; two degeneration identities | non-zero exit, named check |
| `zeros.py` | found-count vs `θ(T)/π` at a ladder of `T` — a jump of 1 is a missed zero | prints `FLAGS`, and the printed `S-ish` column moves by ~1 |
| `phase2_gateA.py` | assembled geometric side vs `Σ_ρ \|F(γ_ρ)\|²`, fitted decay exponent ≈ −3, plus a sign-flip control on the conductor and on the archimedean term separately | exponent far from −3, or a flip that does not blow the discrepancy up (then Gate A is **uninformative** for that term, not passed) |
| `phase3_gateB.py` | `λ_min(N)` monotone non-increasing in two independent complete bases; ζ control in Zhu's interval at `L = 0.8` | `mono=VIOLATED`, or the ζ control leaving `[8.9e−18, 2.27e−17]` |
| `phase5_siegel.py` | **even-sector response `≥ 0` on every row** (the Phase-5 theorem) | printed violation count non-zero |

## Determinism and caching

`src/dirichlet_window.py` caches archimedean special-function values keyed on
`(dps, L, α, kind, a)`, so a precision change can never reuse a stale value. The cache is
in-process only; there is no on-disk cache to invalidate. Zero files in `data/` are the one
persisted intermediate, and `lmfdb_check.py` keeps the LMFDB ordinates it compared against
inside `results/phase2_lmfdb.json` so the comparison is re-checkable offline.

## Not reproduced here

The ζ-side quantities quoted for comparison that are *not* recomputed in these runs are
Zhu's certified interval for `λ*(0.8)` and the parent repository's `C(L)` table in
`data/stage3_C_vs_CS.csv`. Everything else labelled "ζ control" is computed in the same run
as the χ numbers it is compared against.
