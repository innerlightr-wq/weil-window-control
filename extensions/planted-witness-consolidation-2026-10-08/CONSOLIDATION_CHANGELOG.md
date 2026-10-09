# Consolidation changelog — 2026-10-08

Verdict: **`VERIFIED AND CONSOLIDATED`**

## Created in this directory

| file | what it is |
|---|---|
| `INPUT_SHA256.txt` | sha256 of all 12 inputs, taken **before** any work; re-verified at the end, all `OK` |
| `verify_certificate_consolidated.py` | the certificate verifier, with four changes (below) |
| `cold_run.log` | the cold run: `real 5m25.343s`, `EXIT=0`, four gates `PASS` |
| `cold_certificate.json` | the certificate, recomputed with no cache present |
| `interval_cells.json` | all 256 δ-cells with exact rational endpoints and per-cell bounds |
| `regenerated_qzeta_cache.json` | new cache, with a full key and a provenance block |
| `independent_point_check.py` | separate implementation, standard library only, shares no code |
| `indep_run.log` | `real 3m17.482s`, `EXIT=0`, `END-TO-END NEGATIVE UPPER BOUND: True` |
| `independent_point_check.json` | its enclosures as exact rationals |
| `audit_primitives.py` | 30 interval-primitive properties — all pass |
| `audit_cache_invalidation.py` | 1 control + 10 mutations + 1 corruption control — all pass |
| `ANALYTIC_AUDIT.md` | the proofs the arithmetic rests on |
| `VERIFICATION_REPORT.md` | the report, the theorem, and the three corrections |
| `REPRODUCE.md` | environment, commands, exit codes, runtimes, settings |
| `SHA256SUMS` | hashes of everything produced here |

## Changes to the verifier, relative to `planted-witness-closure-2026-10-08/verify_certificate.py`

1. **Cache key completed.** It previously omitted `gamma`, so changing the detection height
   would have silently reused a stale `Q_ζ`. Now `(x, xd, M, grid, L, gamma, ndim, src)`.
   This was a real defect, found by audit, not a cosmetic change.
2. **Cache renamed and given provenance** — `regenerated_qzeta_cache.json`, recording
   producer, `sys.version`, `platform.platform()`, UTC time, and the note that a hash
   authenticates file identity, not mathematical correctness.
3. **δ-subdivision made exact and machine-checkable.** Cells are now
   `[(800+3j)/3200, (803+3j)/3200]`, `j = 0..255`, emitted to `interval_cells.json`, with gap
   and endpoint assertions that raise `CertificationFailure`.
4. **Output renamed** `cold_certificate.json`, so the original `certificate.json` is untouched.

## Manuscript changes (`paper/weil_window_control.tex`)

New PDF `paper/weil_window_control_rev_certified_witness.pdf` (16 pages). Built with
`latexmk -pdf`, **no undefined references**, every numeral verified present in the text
layer. `paper/weil_window_control.pdf` is **unchanged and byte-identical**
(`4ed8d7ee…`).

1. **New §6.1 "A certified four-mode witness for an added quartet"** (`sec:certwitness`),
   tier T1, with **Theorem 11** stating the three certified inequalities; paragraphs on how
   the geometric side was enclosed, the independent reimplementation, and what it does and
   does not say.
2. **Limitation 4 amended** — `mpmath` has no digamma or Lerch, so the geometric side could
   not be interval-evaluated. At the point of §6.1 it now is, via elementary integrals. The
   route is specific to that window and does not lift to the sweeps.
3. **Limitation 5 extended** — a non-detection bounds the search, *and not the basis*: the
   quarter-wave frequencies do not depend on `N`, so the bases are nested.
4. **Limitation 7 extended** — Theorem 11 is one ordinate, one window, one basis.
5. **New Limitation 8** — the controlled-Taylor region (`δ ≲ 0.166`) and the
   certified-negative region (`δ ≥ 1/4`) **do not overlap**; the gap is `(0.166, 0.25)`.
6. **New Limitation 9** — no formal verification; CPython remains a dependency; the second
   implementation is not external peer review.
7. **Retractions 19 and 20 added**, count corrected from eighteen to **twenty**:
   * **19** — "`R ≤ 2K_eff δ²` fails at every positive δ" was **false**; it holds for
     `δ ≲ 0.166`. The obstruction is the non-overlap, not a failure of the bound.
   * **20** — "dimensions 8 and 16 do not detect" retracted as a statement about the basis;
     the bases are nested, so they *do* detect. Only the search failed.
8. **Reproduction appendix** — three rows for the new scripts, and a note that they use the
   standard library only, with runtimes and exit codes.

## Not done, deliberately

* No commit, push, tag, branch change, reset, upload or DOI. Nothing was installed.
* `planted-witness-closure-2026-10-08/` was **not modified**; its `certificate.json` and
  `qzeta_cache.json` still hash as recorded.
* `paper/weil_window_control.pdf`, `src/zeta_window.py` and `notes/proofs.md` unchanged.
  `paper/weil_window_control.tex` **is** changed (that is the authorised manuscript edit);
  its pre-edit hash `bd71ec14…` is in `INPUT_SHA256.txt`.
* No new RH research branch proposed.
