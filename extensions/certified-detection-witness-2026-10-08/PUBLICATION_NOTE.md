# Publication note

This directory records **how the frozen rational witness vector `x` was selected**. It is
published as provenance for §6.1 / Theorem 11 of the technical note
([doi:10.5281/zenodo.23250932](https://doi.org/10.5281/zenodo.23250932)). It is **not**
needed to verify the certificate — for that, see
[`../planted-witness-consolidation-2026-10-08/`](../planted-witness-consolidation-2026-10-08/).

`SOURCE_HASHES.txt` is a historical record of the inputs this selection run read, and is
published unedited. Two of its twelve lines do not resolve or match in the published tree:

| line | status |
|---|---|
| `../../extensions/newtonian-weil-response-2026-10-08/RESULT.md` | **not published.** That exploration (verdict: MISMATCH) is unrelated to this theorem and was deliberately left out of the publication set. The hash is retained so the record of what was read stays complete. |
| `../../paper/weil_window_control.tex` | **does not match.** The manuscript was edited after this run — that edit is §6.1 itself. The previous PDF, `../../paper/weil_window_control.pdf`, still matches. |

The other ten resolve and match:

```bash
grep -v -e newtonian -e 'weil_window_control.tex' SOURCE_HASHES.txt | md5sum -c -
```

Checksum records in this repository are never edited to make changed bytes appear
previously verified; differences are documented instead.
