# SOURCES_AND_STATUS.md

Audit date 2026-10-07. Nothing in this directory modifies the existing paper, code or data.

## 1. Local project

| Item | Finding |
|---|---|
| Local checkout before this audit | **None existed.** No directory matching `*weil*` under `/home/elias` (depth 7), and no git remote on the machine matching `weil-window-control` (all 20+ remotes enumerated). |
| Action taken | Cloned fresh to `/home/elias/GitHub/weil-window-control` (the path was empty; nothing overwritten). |
| HEAD | `f6033d7f0dbb8b420510b90fe410fb96705de8f8` — *"Add comb-dips extension: t0/T1 measurement and null-model test (verdict NULL)"* |
| `git status` | clean, tracking `origin/main`, 161 files. Not reset, not switched, not committed to, not pushed. |
| Repo on GitHub | public, `main`, pushed `2026-10-07T01:24:37Z` |

### The Downloads PDF vs the committed PDF — same document

| | `paper/weil_window_control.pdf` | `~/Downloads/weil_window_control.pdf` |
|---|---|---|
| bytes | 474,032 | 477,046 |
| sha256 | `2f8001fc…86b60` | `851e1d3b…c5c73` |
| pages | 11 | 11 |
| CreationDate | Tue Oct 6 14:24:11 2026 EDT | Tue Oct 6 14:24:11 2026 EDT |
| ModDate | 14:24:11 | 14:29:12 |
| Producer | `pdfTeX-1.40.29` | `macOS … Quartz PDFContext, AppendMode 1.1` |

`pdftotext -layout` output is **byte-identical (0-line diff)**. The Downloads copy is the committed
PDF re-saved through macOS Preview five minutes later, which rewrote the container (+3 KB) without
touching content. **Every claim in this audit therefore attributes to the committed version**, and
no version-splitting is needed.

The manuscript carries `doi:10.5281/zenodo.23195177`; title *"A δ²–L³ Sensitivity Law for
Finite-Window Weil Positivity, with a Function-Field Control of the Connes–van Suijlekom
Pipeline"*. *(Recorded after the audit: that is the version DOI of the manuscript as audited.
The Retraction-18 revision was subsequently deposited as `10.5281/zenodo.23212998`; the
concept DOI for all versions is `10.5281/zenodo.21109955`.)*

### Earlier audits in the local workspace — NOT PRESENT

Searched `~/scratch`, `~/GitHub`, `~/Downloads`. The only directories matching
coupling/inheritance/window are **Collatz/EOC work**, not RH:
`collatz-inner-outer-coupling-2026-09-26`, `eoc-inherited-membership-dynamics-2026-09-26`,
`eoc-windowed-placement-note-2026-09-27`, `eoc-rho-certification-2026-09-28`.
**No resolution-inheritance or coupling audit for this project exists locally.** Nothing in this
report reuses one. The already-tested routes of §7 of the task are known to me only through
`extensions/comb-dips/NULL.md` and `REPORT.md` in the repo itself.

## 2. OpenAI sources

**Pinned commit: `adc7f1241b42` (`openai/math`), authored `2026-10-06T21:58:50Z`, message
"Initial commit".** The repository has exactly **one** commit; the two family-003 preprint
directories and `lean/docs/003.md` were each last touched by it. **There are no subsequent
corrections to pin against.** Repo size ≈ 799 MB.

| Source | Status |
|---|---|
| `preprints/The-Quasi-Riemann-Hypothesis-September-30-2026/build/paper.tex` (7/8) | **CHECKED** — fetched, 766,316 B, 16,677 lines |
| `preprints/The-Quasi-Riemann-Hypothesis-October-5-2026/build/paper2.tex` (11/12 alternate) | **CHECKED** — fetched, 169,005 B |
| `lean/docs/003.md` (scope page) | **CHECKED** — read in full |
| `lean/formalization.yaml` | **CHECKED** — used to map family 003 to its proof files |
| `lean/ComparatorChallenges/{QuasiRiemannHypothesis,DirichletSevenEighths,SiegelZeros,HeckeSevenEighths}.lean` | **CHECKED** — read |
| `lean/OAI/NumberTheory/DirichletL/Nonvanishing.lean` | **CHECKED** — read (33 lines) |
| `lean/OAI/NumberTheory/DirichletL/Detector/FinalAssemblyUnconditional.lean` and its dependency closure | **NOT CHECKED** — not fetched, not built |
| `openai.com/index/sharing-ai-progress-in-mathematics/` (announcement) | **NOT CHECKED** — deliberately; the task asks for the mathematics, and the manuscripts and Lean tree were obtained directly |
| The other 371 families | **NOT CHECKED** — deliberately; no load-bearing dependency required them |
| Zhu, Hallouin–Perret, Connes–van Suijlekom, Bombieri, Groskin (my own references) | **NOT CHECKED here** — relied on as cited by my manuscript and `PRIOR_ART.md`; Zhu's certificate is used only as my paper already uses it |

### The four-way distinction the task demands

1. **Manuscript claim** — `paper.tex` Theorem 1.1 (`\label{thm:main}`, line 106), verbatim:
   *"Every finite-order Hecke $L$-function over $F=\mathbb Q(\sqrt{-3})$ has no zero in
   $\Re s>7/8$. The same holds for every Dirichlet $L$-function, including $\zeta(s)$. A pole at
   $s=1$ for a principal character is allowed."* The paper states explicitly that it does **not**
   establish RH.
2. **Stated formalization scope** — `lean/docs/003.md`: θ = 7/8 for ζ and every Dirichlet
   $L$-function, uniformly over all moduli and characters; also finite-order Hecke
   $L$-functions over $\mathbb Q(\sqrt{-3})$; principal-character poles excluded; later
   applications not included. For the Landau–Siegel corollary it says, verbatim, **"No explicit
   value of $c$ is given"** — so that constant is **not effective**.
3. **Proof personally reconstructed** — the *dependency chain and conversion step* only, from
   `paper.tex` §§1–3 + overview and `paper2.tex` §§1–3, §sec:reduction. Set out in
   `BRIDGE_AUDIT.md`. I did **not** reconstruct the large-sieve, theta-reflection or
   detector estimates; those are cited, not verified.
4. **Formal proof actually checked** — **NONE. NOT CHECKED.** I did not build the Lean project
   (≈799 MB plus a full Mathlib build). What I did verify:
   - The `ComparatorChallenges/*.lean` files are **statement stubs ending in `sorry`** — they are
     the benchmark targets, *not* proofs. `QuasiRiemannHypothesis.lean` is 169 bytes and its body
     is `sorry`. **A comparator file is not evidence of a proof.**
   - The real declaration is `OAI.riemannZeta_ne_zero_of_seven_eighths_lt_re` in
     `OAI/NumberTheory/DirichletL/Nonvanishing.lean`, a one-line delegation to
     `SevenEighths.ProbeFinalAssemblyUnconditional.zeta_nonzero`. That file contains no `sorry`,
     `axiom`, `admit` or `native_decide`.
   - GitHub code search reports `sorry` **0** times under `lean/OAI` and **399** times under
     `lean/ComparatorChallenges`. Code search is index-dependent and **not a substitute for a
     build**; treat it as suggestive only.

No theorem number or quotation in this audit is invented; every one is cited to file and line.
