# SOURCE_AUDIT.md

Exact versions of every source used, located by theorem/section, and an explicit
list of what could not be obtained. Filenames were **not** trusted; every paper
was identified from its own title page and DOI line.

---

## 1. Repository state

| item | value |
|---|---|
| repository | `innerlightr-wq/weil-window-control` (local checkout `~/GitHub/weil-window-control`) |
| branch at start | `main` |
| HEAD at start | `ad78c59a6ca7d2a376c86220e583c5ae94e4962b` |
| working tree at start | **dirty, preserved**: `M README.md`; untracked `extensions/newtonian-weil-response-2026-10-08/`, `extensions/odd-sector-residual-margin-2026-11-02/`, `extensions/spectral-isometry-transfer-2026-11-02/`, `extensions/suzuki-function-field-window-2026-11-02/`, two `paper/*.fls`/`*.fdb_latexmk` build artefacts |
| work branch | `positive-weil-certificate-2026-11-09` |
| `AGENTS.md` / `CLAUDE.md` | **none present** at any level — no repository agent instructions exist to follow |
| environment | `.venv` Python 3.14.7, `mpmath 1.3.0`, `python-flint 0.9.0` (Arb). `numpy`, `gmpy2`, `sympy` absent from the venv. |

Nothing outside `extensions/positive-weil-certificate-2026-11-09/` was modified.

## 2. Papers read, by content not filename

`~/Downloads` holds many same-size duplicates; the version actually used is the
one whose title page carries the stated DOI.

| paper | file used | md5 (first 12) | identifying line | version taken |
|---|---|---|---|---|
| weil_window_control_rev_certified_witness | `weil_window_control_rev_certified_witness.pdf` | `b98d95e6807e` | `doi:10.5281/zenodo.23195177` | all 8 copies are byte-identical |
| confluent_blind_spots_rev_rank_counterexample | `confluent_blind_spots_rev_rank_counterexample.pdf` | `7af76dabf151` | created 2026-10-25, 727818 B | the `rev_rank_counterexample` revision, **not** `confluent_blind_spots.pdf` (`dc8c8ed093f4`, 2026-10-16) |
| dh_note | `dh_note (1).pdf` | `bfb941051770` | `doi:10.5281/zenodo.23608288` | **later** of two; `dh_note.pdf` (`dc1e4db2e9c7`) still has `doi:10.5281/zenodo.XXXXXXXX` |
| siegel_note | `siegel_note.pdf` | `e602a93d469b` | `October 2026 — Technical Note, **version 2**`, `doi:10.5281/zenodo.23624161`, naming v1 `…23623473` | **v2**. `siegel_note (1).pdf` (`a0164bf2e705`) is v1 and is **not** used |
| lambda_note | `lambda_note.pdf` | `1bff692cc4e9` | `doi:10.5281/zenodo.23628224` | single version |

**Version-mixing guard, as instructed.** `lambda_note` is an *earlier* DH
heat-flow bracket. The repository's `extensions/mean-gap-dynamics-2026-10-25`
(merged at HEAD, commit message "certify Lambda_DH > 0.16") is a **stronger,
later** result. No constant from either is used in this build, and none is mixed:
this extension does not touch the de Bruijn–Newman constant at all.

## 3. External source for the finite-to-infinite reduction

> Xuefeng Zhu, *Weil positivity in compact windows: a finite reduction, certified
> two-sided bounds, and a Landau–Widom decay law*, arXiv:**2608.24827**.

| | |
|---|---|
| versions that exist | **v1** 25 Aug 2026 (231 KB), **v2** 2 Sep 2026 (460 KB). **No v3** — checked on the arXiv abstract page, 2026-11-09. |
| version used | **v2** |
| what was obtained | **UPDATED 2026-11-10:** the full PDF, `~/Downloads/2608.24827v2 (1).pdf`, md5 `cff639d6119d`, 34 pp, title page `arXiv:2608.24827v2 [math.NT] 2 Sep 2026`. Previously only the abstract page and the rendered HTML (`arxiv.org/html/2608.24827v2`) were held, and the citations below were checked against those. |
| **what was NOT obtained** | the **supplementary materials**. §17 of the paper describes them — source for both pipelines, software versions, zero data, certified matrices and Cholesky factors, run logs, error budgets, SHA-256 hashes, and an independent verifier `verify_certificate.py` returning `Q(f) >= 9.0e-18` on the archived `M200`. They were **not obtained or executed** here, so **Theorem 1.1's matrix step is still not re-run and not independently verified.** The earlier wording "no ancillary files were located" is corrected: they are described as existing, not offered-and-absent. |

Located in the source and used here:

| item | location | how used |
|---|---|---|
| `Q(f) = 2F(i/2)^2 + (1/2pi)int_R |F|^2 Psi_L` with `Psi_L = Re psi(1/4+it/2) - log pi - comb` | §1 | normalization cross-check against `src/zeta_window.py` |
| **Theorem 1.1**, the one-stroke reduction, with `T_1 = 2pi e^{A_L}`, `beta* = log(T#/2pi) - 1/T# - A_L`, Legendre representation `beta I + 2pp^T + C`, errors `eps_D, eps_B` | §1 | IMPORTED as R1/R4; **R4 not executed** |
| **Lemma 3.2**, `sup_t P_L(t) = A_L`, optimal by Weyl equidistribution on `{log p}` | §3 | IMPORTED as R2 |
| **Theorem 1.4**, the doubly exponential threshold, `A_L = (4+o(1))e^L` | §1 | quoted in REPORT as the reason the route is capped |
| parameters at `L = 0.8`: `T# = 200`, `A_L ≈ 2.9420`, `beta* ≈ 0.5135`, `N = 200` even Legendre modes, 50 digits, 800 panels × 32-pt Gauss–Legendre, `lambda_min(M) >= 9e-18` | §1, §4 | **reproduced**: `A_L` and `beta*` recomputed here in Arb (Phase B) |
| abstract bracket `8.9e-18 <= lambda_min(0.8) <= 2.27e-17` | abstract; also Table 3 at `L=0.8` | used only as an external comparison point, never as an expected answer. **Verified in the PDF 2026-11-10** (abstract line 5; Table 3 row `0.8`) |
| **Lemma 2.5**, exact time-domain archimedean term, with the removable limit `c - (3/2)g(0)` and `c != 0` in general | §2.1 | transcribed into this normalization for the compact evaluation of `Q`; **verified verbatim 2026-11-10** |
| **Lemma 6.1**, parity splitting, pole term `2(c^2 - s^2)` | §6.1 | IMPORTED; cited for the sector signs |
| **Theorem 6.2**, simple even ground state at support 1.6 | §6.3 | `8.9e-18 <= lam_even_1 <= 2.523e-16`, `lam_even_2 >= 2.085e-12`, `8.206e-15 <= lam_odd_1 <= 2.347e-14`; ground state **simple and even**. This is **at `L = 4/5`, the window of this build**, and it is STRONGER than the parity separation obtained here (odd constant `2.3e-17`, about 357x weaker). Recorded 2026-11-10; the note's §1.3/§1.4/§5 are corrected accordingly |
| **§6.2**, exact Rayleigh evaluation of the upper bound in the time domain, polynomials of degree <= 239 by Gauss-160 | §6.2 | the **same route** used in this build's §4; credited there |
| numerical reference `lambda_min(Q|window) = 1.656e-17` at `L = 0.8` | §5, from the exact time-domain assembly | **non-certified reference value.** It lies inside this build's certified bracket; used as a consistency check only, never as an endpoint |

**Reproduced independently here:** `A_L = 2.94197352522…` (vs the paper's 2.9420)
and `beta_paper(200) = 0.5134667749…` (vs 0.5135), both as Arb balls, from the
prime-power list `{2,3,4}` derived from `log n < 2L` by certified comparison.

**Not reproduced:** the `200 x 200` Legendre/Bessel assembly, its `lambda_min >= 9e-18`,
the tail constants `eps_D, eps_B`, and therefore the whole-window lower bound.

## 4. Prior audits read before proposing anything

| audit | what it already settles | consequence here |
|---|---|---|
| `four-mode-positivity-2026-10-14` | certified **negative** `Q_add,delta` on `delta in [9/100,49/100]`, 69 rational cells, one frozen vector | re-run unchanged as the **regression control** (passes, exit 0) |
| `planted-witness-consolidation-2026-10-08` | the certified interval path, published cells on `[1/4,49/100]` | provenance for the above |
| `channel-compatibility-exploration-2026-10-24/EXPLORATION.md` | the background-normalized interaction `eta = 1 - det Q/det B` "only reorganizes the problem"; a new estimate must live in a bound on `E^T S^{-1} E` from **analytic** data | **honoured**: no Schur/determinant reformulation is offered here as new content |
| `intrinsic-positivity-bridge-2026-10-07/POSITIVITY_DEFECT_AUDIT.md` | the positivity defect and signed discrepancy | no bridge re-proposed |
| `sonin-first-prime-defect-2026-10-07/RESULT.md` | SOURCE/DOMAIN GAP; the `Tr(theta_S(f)Pi_S)` identity could not be justified; ceiling `L < log 2` **retracted**; parity reading **withdrawn** | the Sonin route is **not reopened**; no new estimate addressing its gap is claimed |
| `davenport-heilbronn-2026-10-23` | the DH normalization `a_n = (1,k,-k,-1,0)` mod 5, `kappa = [sqrt(10-2sqrt5)-2]/(sqrt5-1)` | used verbatim for control B |
| `dirichlet-windows-2026-10-23` | Dirichlet-window sector behaviour | control C scoped to coefficient signs only |
| `source-flow-detector-audit-2026-10-24` | detector/source-flow audit | no flow claim made |

**Chronology rule applied.** `siegel_note` **v2** supersedes v1's withdrawn claim;
`mean-gap-dynamics-2026-10-25` supersedes `lambda_note`'s bracket. Neither
superseded value is used.

## 5. Oscillator repository (secondary, read-only)

`two-mode-oscillator-record-chronologies` was **not** used. Its quadratic-unit /
Pell classifications concern `Z[sqrt D]` units and odd/odd congruence classes;
nothing there transfers to logarithms of primes or to zeta zeros, and no such
transfer is invented. The only methodological habits carried over are generic and
were already repository practice: exact/interval comparisons for sign-critical
decisions, justified rather than assumed cutoffs, and frozen reproducible
certificates.

## 6. Unavailable / not attempted

* Zhu's reproducibility materials — **not located**; matrix step unverified.
* A quantitative classical zero-free region with stated constants — **not
  attempted**; only the boundary theorem is proved.
* `numpy` in the project venv — absent; not needed.
