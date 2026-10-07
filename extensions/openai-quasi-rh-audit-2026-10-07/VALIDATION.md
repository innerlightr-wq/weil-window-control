# VALIDATION.md — what was actually run, and what could not be

## Environment limitation, stated first

**The project's own regression scripts could not be run here.** `requirements.txt` pins
`mpmath==1.3.0`, and `mpmath` is absent from this machine with no way to install it: no `pip`
(`No module named pip`), no `uv`, no `pipx`, and no `python-mpmath` system package. Every
stage of `scripts/quick_check.sh` except the finite-field block imports `mpmath`, as do
`src/proofs_check.py`, `src/task_final_checks.py` and all of `make law`.

So the following were **NOT run**, and nothing below should be read as re-certifying them:
`make quick` (LDL inertia self-test, Stage 1 regression, Stage 2 exact-over-Q, function-field
`δ²` law), `make law`, `make gates`, `make connes`, `make onset`. Their committed outputs in
`data/` are untouched and are carried over from the original runs, not reproduced today.

## What WAS run today

### 1. The mpmath-free part of `scripts/quick_check.sh`

```
PASS: Fermat check for F_5, F_25, F_27   (src/ffield.py)
```

### 2. `validate_revision.py` (new, pure standard library, double precision)

Independent re-derivation of exactly the three things the revision touches, plus the
withdrawn inference. Full transcript in `VALIDATION_raw.txt`. Summary:

| Block | Checks | Result |
|---|---|---|
| Perturbation expansion `Q_δ − Q_0 = −4δ²(F′²+FF″) + R_4`, three (L, f) pairs × three δ | 18 | all PASS |
| Remainder bound `\|R_4\| ≤ (16/3)L⁵δ⁴e^{2Lδ}` | (in above) | all PASS |
| Theorem 7's per-`f` inequality `Q_δ(f) ≥ Q_0(f) − B_L(δ)` | (in above) | all PASS |
| Cauchy–Schwarz constants `\|F\|²≤2L`, `\|F′\|²≤2L³/3`, `\|F″\|²≤2L⁵/5`, `4\|F′²+FF″\|≤c(L)` | 8 | all PASS |
| `c(0.8) = 3.197120 ≈ 6.244·0.8³` | 1 | PASS |
| `K(L,γ)` closed form **as bracketed in the TeX** vs the defining integral, 5 windows | 5 | all PASS |
| `4K(L,γ₁)` vs the manuscript's printed table (0.7419, 1.3427, 3.1158, 5.1130, 10.652) | 5 | all PASS |
| `f ∝ u sin(γu)` attains `\|F′(γ)\|² = K`; `K ≤ 2L³/3` (Cor. 9) | 4 | all PASS |
| Logical regression check (below) | 2 | all PASS |

**`ALL CHECKS PASSED`**, exit 0.

### 3. Logical regression check (the point of Retraction 18)

> `X ≥ b − D` permits replacing `b` by `0` only given a separate `b ≥ 0`.

| case | premise `X ≥ b−D` | `b ≥ 0`? | substituted `X ≥ −D` |
|---|---|---|---|
| `b = −1, D = 0, X = −1` | holds | **no** | **fails** |
| `b = 8.9e−18, D = 1e−20, X = 8.8e−18` (the `L ≤ 0.8` case) | holds | yes | holds |

The first row is a statement about **scalars**. It demonstrates the failure of the discarded
inference; it is **not** an example of a negative actual Weil form.

## Build

```
pdflatex -interaction=nonstopmode -halt-on-error \
         -jobname=weil_window_control_revised_2026-10-07 weil_window_control.tex   (×2)
```

`-jobname` was used deliberately so the published PDF is never written to.

| | original | revised |
|---|---|---|
| file | `paper/weil_window_control.pdf` | `paper/weil_window_control_revised_2026-10-07.pdf` |
| sha256 | `2f8001fc…86b60` (**unchanged**) | — |
| pages | 11 | 12 |
| errors | 0 | 0 |
| LaTeX warnings | 0 | 0 |
| undefined refs/citations | 0 | 0 |
| overfull hbox | 2 | 2 |
| underfull hbox | 1 | 1 |

Build diagnostics are **identical** to the original source built the same way. (An interim
draft had one extra overfull box from the long audit-directory path; a break opportunity was
added and the profile returned to baseline.)

## Numbering parity, checked mechanically against the original PDF

```
equation numbers      orig: 1 2 3 4 5     rev: 1 2 3 4 5     MATCH
theorem/cor/remark    orig: 1 2 3 5 6 7 8 9 10 36            MATCH
                      rev : 1 2 3 5 6 7 8 9 10 36
```

An interim draft shifted Theorem 8 → 9 (an inserted numbered remark) and `K` from eq. (5) →
(6) (an inserted numbered equation). Both were caught by this check and reverted by making the
inserted remark a `\paragraph` and the inserted display unnumbered.

## Rendered-page inspection (not only extracted text)

Pages 6 and 11 were rasterised at 140–150 dpi and inspected visually, because `pdftotext`
renders `\left[ \right]` as stray `"` and `#` — the artifact that caused a false finding during
the audit. Confirmed on the rendered page: `K(L,γ)`'s closed form carries its **large square
brackets with the original signs**; `Q_rest` displays correctly; Theorem 8, Corollary 9 and
Remark 10 keep their numbers; the §4.3 table is unchanged; Retraction 18 appears after 17 with
entries 1–17 intact.
