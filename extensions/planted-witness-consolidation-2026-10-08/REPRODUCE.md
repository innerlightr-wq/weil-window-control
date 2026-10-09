# How to reproduce

## Environment actually used

```
python3        3.13 (CPython), standard library only -- no numpy, no mpmath, no gmpy2
platform       Linux 7.2.5-3-omarchy, x86_64
date           2026-10-08 (UTC)
packages       none installed for this run
```

No third-party package is required and none was installed. `mpmath` is used by
`src/zeta_window.py` (the exploratory float code) but **not** by anything in this directory.

## Commands, in order, with observed exit codes and runtimes

```bash
cd extensions/planted-witness-consolidation-2026-10-08

# 1. cold certificate: delete any cache first, so nothing is reused
rm -f regenerated_qzeta_cache.json cold_certificate.json interval_cells.json
time python3 -u verify_certificate_consolidated.py          # EXIT=0   real 5m25.343s
#    -> cold_certificate.json, interval_cells.json, regenerated_qzeta_cache.json
#    -> four gates printed, all PASS

# 2. separate implementation of the point inequality
time python3 -u independent_point_check.py                  # EXIT=0   real 3m17.482s
#    -> independent_point_check.json
#    -> prints "END-TO-END NEGATIVE UPPER BOUND: True"

# 3. interval-primitive audit            (30 properties)
python3 audit_primitives.py                                 # EXIT=0   ~4 s
#    -> "AUDIT PASSED: every interval-primitive property held."

# 4. cache-invalidation audit            (1 control + 10 mutations + 1 corruption control)
python3 audit_cache_invalidation.py                         # EXIT=0   ~45 s
#    -> "CACHE AUDIT PASSED: every input change invalidates the cache."
```

Step 4 copies the verifier and the cache into a fresh `tempfile.mkdtemp()` directory and
mutates only the copy; **the cache in this directory is never modified by it.**

A second run of step 1 *without* deleting the cache loads `Q_ζ` from
`regenerated_qzeta_cache.json` and finishes in about 30 s. That path is for convenience only;
**the cold path is the validation**, and step 1 as written above forces it.

## Exact inputs (frozen; do not renormalize)

```
L      = 4/5                       window half-length
U      = 2L = 8/5
gamma  = 14                        exact integer, not an approximation to a zero ordinate
NDIM   = 4                         basis dimension
basis  = EVEN_D: phi_k(x) = cos((2k+1) pi x / (2L)) / sqrt(L),  k = 0,1,2,3
XNUM   = [99815221908, -354082045264, 928280739889, 54384691304]
XDEN   = 10^12
||f||^2 = 999999999999981276942897 / 10^24     (exact; Gram = I)
```

`‖f‖² ≠ 1` by `1.87×10^{-14}`. This is deliberate: the witness is a frozen rational vector
and **is not renormalized**. The functional is homogeneous of degree 2, so the sign is
invariant under scaling, but the printed values are not.

## Arithmetic settings and tail parameters

| | `verify_certificate_consolidated.py` | `independent_point_check.py` |
|---|---|---|
| endpoints | `fractions.Fraction`, exact | `fractions.Fraction`, exact |
| outward rounding grid | `2^{-260}` | `2^{-200}` |
| series remainder target | `2^{-300}` | `2^{-240}` |
| archimedean method | `M = 131071 = 2^{17}−1` Laplace integrals | `ψ`, `ψ'` closed form |
| archimedean tail | `[¾C(0)+½sup|C'|]/(M−¾) = 4.935033775516166e-05` | `g`-tail `≤ 2.49e-05`, `ψ'` midpoint `≤ 7.2e-07` |
| `ψ` recurrence depth | — | `JREC = 50` (`Re z → 50.25`) |
| `g`-series terms | — | `JG = 20000`, `NG = 8` |
| `E`-series terms | — | `MEXP = 24`, geometric tail |
| `π` | Machin, 60 terms | Machin |
| `exp` | direct series | range reduction `t = 2^j u`, `|u| ≤ 1`, then `j` squarings |
| `sin`/`cos` | Taylor, Lagrange remainder, `nmax = 1200`, no range reduction (refuses `|t| ≳ 371`) | Taylor, `nmax = 1200` |
| `sqrt` | integer `isqrt`, bracket **verified or raise** | same |
| libm | **not used in the certified path** | **not used** |

`exp` in the independent program was rewritten mid-run: the original direct series raised
`r_exp truncation` at `t = −74.4` (needed `n ≈ 400`) and would have suffered alternating-sign
cancellation. The replacement reduces to `|u| ≤ 1`, sums positive terms, and inverts for
negative arguments. It was checked against a 60-digit `decimal` computation, agreeing to 17
significant digits; `math.exp` appeared to disagree by `5.7×10^{-15}` relative purely because
`float("-74.4")` is a *different rational number* from `−74.4`.

## Cache usage

`regenerated_qzeta_cache.json` carries the key

```
(x, xd, M, grid, L, gamma, ndim, src='even_d cos((2k+1)pi x/2L)/sqrt(L)')
```

and a `provenance` block recording producer, `sys.version`, `platform.platform()` and the UTC
time, with the note that *a hash authenticates file identity, not mathematical correctness;
the cold computation is the validation.* `audit_cache_invalidation.py` shows every one of
those fields invalidates the cache when changed.

## Unchanged historical inputs

`INPUT_SHA256.txt` records the sha256 of every input this directory consumed, taken before
any work began:

* all seven files of `extensions/planted-witness-closure-2026-10-08/` (including
  `verify_certificate.py`, sha256 `c49858a1…`, and the original `certificate.json` and
  `qzeta_cache.json`);
* `src/zeta_window.py`, `notes/proofs.md`, `paper/weil_window_control.tex`,
  `paper/weil_window_control.pdf`;
* `extensions/certified-detection-witness-2026-10-08/RESULT.md`.

**Nothing in those directories was modified.** The original cache was not overwritten: the
consolidated verifier writes to `regenerated_qzeta_cache.json` in *this* directory. Re-check
at any time with `sha256sum -c INPUT_SHA256.txt` from the repository root.

---

## Manifest status at publication (read this before running `sha256sum -c`)

`INPUT_SHA256.txt` is a **historical record**: the sha256 of every input as it stood
*before* the verification run, on 2026-10-08. It is published unedited, and two of its lines
will **not** match the committed tree. That is expected, and it is recorded here rather than
corrected in the manifest, because editing a checksum record to match changed bytes would
destroy the only thing it is good for.

| line | status |
|---|---|
| `../../paper/weil_window_control.tex` → `bd71ec14…` | **will not match.** The manuscript was edited *after* the verification run — that edit is precisely §6.1 / Theorem 11, which reports this certificate. The committed source hashes `4ef9646b…`. |
| all 7 files of `../planted-witness-closure-2026-10-08/` | match |
| `../../src/zeta_window.py`, `../../notes/proofs.md` | match |
| `../../paper/weil_window_control.pdf` → `4ed8d7ee…` | match (the previous, superseded paper, unchanged) |
| `../certified-detection-witness-2026-10-08/RESULT.md` | match |

So the honest check is:

```bash
# from this directory: everything except the deliberately-edited manuscript source
grep -v 'weil_window_control.tex' INPUT_SHA256.txt | sha256sum -c -
```

`SHA256SUMS` is the **publication manifest** for the files produced here plus the current
paper PDF. It excludes itself. Verify from this directory:

```bash
sha256sum -c SHA256SUMS
```

### The published paper PDF

`../../paper/weil_window_control_rev_certified_witness.pdf` is **byte-identical to the file
deposited on Zenodo** as doi:10.5281/zenodo.23250932:

```
md5    b98d95e6807e4ea249d7e4bc648c3bbb
sha256 8f782a40edc893958e3e13e995787a5dd38132f712287719cbcaa17d1df983d9
size   527563 bytes
```

Those are the deposited bytes, committed unchanged. The DOI is advertised in `README.md` and
`CITATION.cff`, **not** inserted into the PDF, precisely so the committed file and the
archived file stay identical. `make paper` rebuilds from `weil_window_control.tex`; a rebuild
is not guaranteed bit-identical (LaTeX embeds timestamps), so compare content, not bytes.
