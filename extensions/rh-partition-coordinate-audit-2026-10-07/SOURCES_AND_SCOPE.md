# SOURCES_AND_SCOPE.md

Dated **2026-10-07**. New files live only in this directory. Everything else in
`weil-window-control` is read-only; the Two-Face manuscript was **not modified**.

## Repository state at the start of this pass (= at the end)

| | value |
|---|---|
| HEAD | `832bbbb` "Add the one-prime Sonin discrepancy test; verdict SOURCE/DOMAIN GAP" |
| branch | `main`, 0 ahead / 0 behind `origin/main` |
| remote | `https://github.com/innerlightr-wq/weil-window-control.git` |
| `git status --short` | clean before this pass |

Local manuscript is **not** ahead of the remote: `paper/weil_window_control.tex` sha256
`bd71ec14f4044455…`, `paper/weil_window_control.pdf` sha256 `4ed8d7eef6d3c7d1…` (md5
`a29f9670844f341c7986cfc1d4382bbe`, byte-identical to Zenodo `10.5281/zenodo.23212998`).

## Local sources read, with exact loci

| source | what was read |
|---|---|
| **`paper/weil_window_control.tex`** | the baseline convention eq. `\label{eq:baseline}` `Q_0 = Q_true + 2F(γ)²`; the perturbation eq. `\label{eq:pert}`; the budget `B_L(δ) = c(L)δ² + (16/3)L⁵δ⁴e^{2Lδ}`, `c(L) = 8L³(⅓+1/√5) ≈ 6.244L³`; **Theorem 7** (`thm:zeta`, relative perturbation estimate); **Theorem 8** (`thm:zetasharp`, sharper/conditional) with `K(L,γ) = sup|F'(γ)|² = ∫_{-L}^{L}u²sin²(γu)du` and its closed form; **Corollary 9** (`cor:uniform`, `K ≤ 2L³/3`); **Remark 10** (`rem:twoconstants`, claims (a)(b)(c)) |
| **`notes/proofs.md` §§B.1–B.3** | the **exact** quartet contribution `4 Re F(γ+iδ)²` versus the on-line pair's `2F(γ)²`; the three conventions and why the third (double pair → quartet) is adopted; the Cauchy–Schwarz bounds; `G(δ) = F(γ+iδ) = ∫f(u)e^{iγu}e^{-δu}du`, `m_k = ∫fu^ke^{iγu}`, `|m_k| ≤ √(2L)L^k`, `|R_4| ≤ 8LΣ_{n≥4}(2Lδ)^n/n! ≤ (16/3)L⁵δ⁴e^{2Lδ}`. **Transform convention: `F(t) = ∫f(u)e^{itu}du`**, so `F'(t) = −∫uf(u)sin(tu)du` and `F''(t) = −∫u²f(u)cos(tu)du` for real even `f`. |
| **`src/zeta_window.py`** | the form assembly `M = Pole + Arch − LogPi − Prime` (`:214`), `von_mangoldt_terms(twoL)` keeping `log n < 2L` (`:148`) |
| `extensions/second-moment-generalization-2026-10-07/PROJECTION_CONTEXT.md` | its corrected §5 and the interpretive principle (magnitude, alignment, response). Carried as **diagnostic context only**; no step below uses the finite-matrix projected-moment or quartic result |
| `extensions/intrinsic-positivity-bridge-2026-10-07/SEMILOCAL_LITERATURE_UPDATE.md` | semi-local Sonin geometry **is** established (CCM arXiv:2310.18423v2); the positive-trace comparison with discrepancy control is not |
| `extensions/sonin-first-prime-defect-2026-10-07/FIRST_PRIME_DEFECT.md`, `RESULT.md` | the Sonin scaling conventions, `ϑ(λ)ξ(v) = λ^{-1/2}ξ(λ^{-1}v)` on `L²(R)_ev`; the amended reach statement and its retraction; verdict SOURCE/DOMAIN GAP |
| `extensions/openai-quasi-rh-audit-2026-10-07/` | carried for its **safeguards only**: Theorem 7 is a *relative* estimate (Retraction 18), "`λ*→0`" was circular, and no step may normalise by an unproved gap of the Weil form or by a nonexistent bounded `L²` operator norm of it. The quasi-RH zero-free half-plane `Re s > 7/8` gives `\|δ\| ≤ 3/8`, used here **only** as a domain |

Prior audits were **not** re-run or re-derived; their results are carried with provenance,
and only the four dependencies actually used above were rechecked against the files.

## Reference manuscript — verified before use

**The Two-Face Problem: A Practical Guide to Translating Binary-Partition Formulas Across
Disciplines**, Elias De Jesús, ORCID 0009-0007-0190-9143, **April 2026 (Revision 2)**.
Local file `~/Downloads/BridgeGuide3 (1).pdf`, md5 `457335493829f3431ee3419ef511972a`,
324,148 B, pdfTeX-1.40.27. **Title page read and matched** before any use. **READ:** abstract,
§2 (conservation law, Thales variables, AM–GM–HM), §3 (the three descriptions), §4 (Gudermannian
bridge, eqs. 14–20), §5 (translation dictionary), §6 (workflow, incl. §6.4 classification and
§6.6 "When Not to Use This Framework"), §7.5 (methodological caveat), §11 (epistemic tiers),
§12 (conclusion), and the Appendix A cautions. §§8–10 and the empirical sections skimmed —
**NOT CHECKED in detail**; nothing below depends on the SXS/SDSS results.

The items the brief names, as the manuscript states them:
`a+b = 1`, `a,b ∈ (0,1)` (eq. 1); `χ = |2a−1| ∈ [0,1)` (eq. 2); `h = √(ab)`, `L = |a−b|/2 = χ/2`,
`R = 1/(1−χ²) = 1/(4h²)` (eq. 3); `h²+L² = ¼` and `4h²+χ² = 1` (eq. 4); the Gudermannian bridge
`θ = gd(ξ) = arctan(sinh ξ)` with `R = 1+J² = sec²θ` and `1/(2h) = cosh ξ = sec θ` (eqs. 14–20);
the translation dictionary of §5 (Fisher `θ` / Poincaré `ξ` / Compound `y = R/2`).

### Two warnings carried from the manuscript, and used below

1. **§7.5, verbatim:** *"Under purely monotonic fits … all candidate symmetric functions of `χ`
   give identical variance-explained at any value of `χ`, because they are strict monotonic
   transformations of one another."* Applied here in its exact-mathematics form: see
   `COORDINATE_RESPONSE_AUDIT.md` §3.
2. **§6.4 and §6.6:** the classification *exact translation / partial structural correspondence /
   true mismatch*, with the instruction that *"If domain, bounds, or singularities do not align
   across translation, the correct diagnosis is likely a mixed-description object or a true
   mismatch"*, and the explicit **"Broken domains"** failure mode. Applied in
   `GENERATOR_AUDIT.md` §4.

The Appendix cautions are the manuscript's own model for the conclusion reached here, e.g.
*"This is an exact statistical-geometry translation. It does not imply that statistical
estimation and gravitational response are physically coupled."*

### NOTATION COLLISION — flagged once, obeyed throughout

The manuscript writes **`γ = 1/(2h) = cosh ξ`** for the "manifold Lorentz factor" (its eq. 6).
The Weil work writes **`γ`** for the **spectral ordinate** of a zero. These are unrelated.
**In every file in this directory `γ` is the spectral ordinate**; the manuscript's factor is
written `LF` where needed. The manuscript's `L = χ/2` also collides with the Weil window
half-width `L`; **`L` here is always the window half-width**, and the manuscript's quantity is
written `χ/2`.

## Not checked

CC2021 / CCM2024 / Connes 1999 were **not re-read** in this pass; the Sonin scaling convention
is taken from the already-audited record in `extensions/sonin-first-prime-defect-2026-10-07/`.
Zhu's certificate, Kato, Weil [W1]: **NOT CHECKED**, as before. No literature search for prior
art on hyperbolic reparameterisations of zero displacement was performed — **absence of a
prior-art note below means "not searched", not "novel"**. The one substantive observation that
could be prior art (§5 of `GENERATOR_AUDIT.md`) is elementary Fourier analysis and is **not
claimed as new**.

## Added in the post-review amendment pass (2026-10-07)

| item | status |
|---|---|
| `V_η := e^{−½tanh(η)A}` has `dV_η/dη = −½sech²(η)·A·V_η`, hence no constant generator in `η`, hence `V_{η₁}V_{η₂} = V_{η₁+η₂}` iff `tanh η₁ + tanh η₂ = tanh(η₁+η₂)` | **standard**; verified in `amend_checks.py` §A (coefficient to `1.7e-11`; group-law defect `2.66e-02` to `5.40e-01`). This replaces the earlier wording that `η` "breaks" the generator identity |
| `D = d/du` is not skew-adjoint on `L²[−L,L]` without boundary conditions; zero-extended window functions are not a translation-invariant space | **standard**; verified in `amend_checks.py` §B (boundary term equals the adjointness defect exactly, `1.6000000`). Scope correction to `GENERATOR_AUDIT.md` §5 |
| Bonneau, Faraut, Valent, *Self-adjoint extensions of operators and the teaching of quantum mechanics*, arXiv:quant-ph/0103153 — momentum-operator domain distinctions | **NOT CHECKED.** Cited as supplied in review, for the domain qualification only; nothing above depends on it |
| DLMF §1.14 (integral transforms) for the spectral-translation / exponential-weighting relationship | **standard reference**, as supplied in review; the identity itself is verified directly |

Neither item changes the verdict. Both narrow a claim.

## Tooling

Python 3 standard library only — `fractions.Fraction`, `decimal`, `math`, `cmath`. **`mpmath` is
not installable on this box** (no `pip`/`uv`/`pacman`), so `src/zeta_window.py` and the `scripts/`
regressions were **not run**; nothing here depends on them. Quadrature is Simpson with
`N = 2000` nodes on `[−L,L]`, stated at every use. No installs, no background jobs, no sweeps.
