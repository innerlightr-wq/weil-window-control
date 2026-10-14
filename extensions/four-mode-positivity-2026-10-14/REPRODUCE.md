# REPRODUCE.md

Standard library only. No `pip install`, no virtual environment, no
`requirements.txt` — nothing here imports a third-party package.

```bash
git clone https://github.com/innerlightr-wq/weil-window-control.git
cd weil-window-control/extensions/four-mode-positivity-2026-10-14
```

Every script resolves its inputs **relative to its own location** and verifies
the sha256 of the published files it imports, so a fresh checkout works with no
configuration and no paths to edit.

## Fast first (≈2 min each)

```bash
python3 tests/test_core_vs_published.py    # reproduces the PUBLISHED baseline as inequalities
python3 tests/test_matrix_controls.py      # generic matrix controls for Proposition 1 and LDL
python3 tests/test_tail_sharp.py           # validates Proposition T   (O(1/M^2))
python3 tests/test_tail3.py                # validates Proposition T3  (O(1/M^3)) and Lemmas E, S1, T
```

The first of these also builds `results/arch_sums_cache.json` (~4 min on the
first run only). Each script prints `PASS`/`FAIL` per check and exits non-zero
on any failure.

## The two verifiers

```bash
python3 verify_new_certificate.py     # Theorem A        ~17 min
python3 verify_route_o.py             # Theorems O1, O2  ~45 min
```

`verify_new_certificate.py` asserts, regenerating every enclosure from the
mathematical input:

1. mechanical coverage — the 69 cells are closed, consecutive, gap-free, and
   their union is exactly `[9/100, 1/4]`;
2. for each cell, an enclosure of `Q_add,delta(f_x)` valid for **every** `delta`
   in that closed cell, with a strictly negative upper bound;
3. the join with the inherited certificate at `delta = 1/4`, from **both** sides;
4. the inherited point bound at `delta = 2/5`, **recomputed, not read**;
5. the sharper tail enclosure is **contained in** the published one.

`verify_route_o.py` asserts:

* `O1` the sharpened `T` and the `O(1/M^3)` tail are each contained in the
  coarser enclosure they replace (nested enclosures of the same number);
* `O2` `Gt_mid` lies inside a freshly polarized interval matrix, `sigma3` is
  recomputed, and `Gt_mid - t_cert I` is positive definite by exact rational
  `LDL` — hence `G >= m0 I`, `m0 > 0`;
* `O2b` `Qsharp3` really is a quadratic form (`Qsharp3(0) = 0`, 2-homogeneity,
  and `x^T Gt x` reproducing it off-axis), and the polarization defect is within
  the row-sum loss in dimension 4;
* `O3` every one of the 6072 cells satisfies `lambda_min(N) >= Loss` by exact
  rational `LDL`, with `N` and `Loss` **recomputed** from `(lo, hi, w, tau,
  tau')` — hence `M(delta) >= 0` throughout the cell;
* `O4` mechanical coverage of those cells, consecutive from `0` with no gap;
* `O5` `Theta < 1` at the right end of the certified non-detection range, and a
  negative witness exists at `delta = 9/100`, so the two certified regimes are
  genuinely different and the gap between them is explicit.

Both exit `0` on success and **non-zero** if any assertion fails. The gates
assert the inequalities themselves; they do not print a stored label.

## Regenerating the certificates from scratch

```bash
python3 scripts/route_o1.py        # -> results/route_o_step1.json      ~5 min
python3 scripts/build_cells.py     # -> results/certificate_new.json    ~10 min
python3 scripts/theta_table.py     # -> results/theta_table.json        ~3 min
python3 scripts/route_o2.py        # -> results/route_o_step2.json      budget-limited
```

`scripts/route_o2.py` runs to a wall-clock budget and records
`complete = false` with the endpoint it actually reached; the published run
reached `3831/64000` in 1550 s with 6072 cells. Re-running it with a larger
budget extends the interval and is the open item of
[`NON_DETECTION.md`](NON_DETECTION.md) §4.

## The second implementation

```bash
python3 scripts/indep_crosscheck.py    # ~8 min
```

Drives the repository's separately coded digamma/trigamma route on the new
witness and then compares the two enclosures explicitly: the primary must sit
inside the independent one and both `Q_add` upper bounds must be negative. It
exits non-zero if they disagree. Its budget (~`1.1e-04` on `Q_zeta` here) is
about `90x` the worst margin of Theorem A, so it corroborates the right part of
`[9/100, 1/4]` only — that limitation is printed by the script, not hidden.

## Cold run

```bash
rm -f results/arch_sums_cache.json
python3 verify_new_certificate.py
```

The cache key carries the sha256 of the published verifier along with every
parameter, so it cannot go stale; deleting it only costs the ~4 min rebuild.
