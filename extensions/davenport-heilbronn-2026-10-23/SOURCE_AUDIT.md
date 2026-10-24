# Source audit

Compiled 2026-10-23. Every entry records what was actually inspected, by what means,
and what could not be reached. Bibliographic details were checked against the source
or against a primary index, not copied from the brief.

Hashes of files downloaded into this extension are in `MANIFEST.sha256`; where a source
could not be downloaded, that is stated and the entry is marked **not inspected**.

---

## A. The Davenport–Heilbronn construction

**A1. H. Davenport and H. Heilbronn**, *On the zeros of certain Dirichlet series*,
J. London Math. Soc. **11** (1936), 181–185 (and 307–312, part II).
Status: **not inspected** (no accessible full text from this machine).
Used for: the historical attribution only. Nothing in this extension depends on reading
it: the function, the value of `κ`, and the functional equation are **re-derived from
scratch** in `DERIVATION.md` §1 from the Gauss sum of the mod-5 character, and the
derivation is checked numerically (Phase A, [A1], [A3]). The constant obtained,
`κ = [√(10−2√5) − 2]/(√5 − 1) = 0.2840790438404122960…`, agrees with the value quoted
throughout the literature, which is the cross-check that the convention matches.

**A2. The character convention.** The brief fixes `χ(0..4) = (0,1,i,−i,−1)`. This is
primitive and odd mod 5. Verified independently here: `|τ(χ)| = √5` to 60 digits and
`χ(4) = χ(−1) = −1`. Different sources use `χ(2) = ±i`; swapping `i ↔ −i` replaces
`κ ↔ −κ`, which maps `F` to the function with coefficients `(1,−κ,κ,−1,0)`. Both satisfy
the same completed functional equation; the brief's convention is the one used
everywhere here and the sign of `κ` is never inferred from a source.

---

## B. Published off-line zero computations

**B1. R. Spira**, *Some zeros of the Titchmarsh counterexample*,
Math. Comp. **63** (1994), 747–748.
Status: **not inspected** (full text behind the AMS paywall; the bibliographic record
was confirmed through secondary citations). Reported content, per multiple secondary
sources: zeros of the Davenport–Heilbronn function off `σ = 1/2` computed in
`0 ≤ t ≤ 200`.

**B2. E. P. Balanzario and J. Sánchez-Ortiz**, *Zeros of the Davenport–Heilbronn
counterexample*, Math. Comp. **76** (2007), no. 260, 2045–2049.
Status: **bibliographic record confirmed** (AMS volume/number/pages match the search
index); **PDF not inspected** — `https://www.ams.org/mcom/2007-76-260/
S0025-5718-07-01999-0/S0025-5718-07-01999-0.pdf` returned **HTTP 403** from this
machine. Reported method, per the AMS abstract text reachable in search results: the
off-line zeros are obtained by deforming a Dirichlet series with known zeros into the
Davenport–Heilbronn series.

**Consequence for this experiment.** The numeric
`0.8085171824566374 + 85.6993484853776 i` given in the brief is used **only as a Newton
search seed**. It was *not* taken as certified, and the refinement in Phase A [A6] moved
it by `7.8 × 10⁻¹⁵`, i.e. the published digits agree to about 14 places and no further.
The refined value is this extension's own:

```
β = 0.8085171824566373855533519606068441278507…
γ = 85.69934848537759217192926770894172903799…
```

stable to `3 × 10⁻¹⁰⁶` between 100 and 140 working digits, with `|F(ρ)| = 4 × 10⁻¹⁰¹`
and the Newton step `|F/F′| = 3 × 10⁻¹⁰¹`. The *whole orbit* was verified. In the
census the same point carries a **local winding number 1** on `|s − ρ| = 0.02`, which is
the statement that it is a simple zero of the analytic function `F`, not an "approximate
zero".

**B3. T. Liu and J. Wu**, *Davenport–Heilbronn Function Ratio Properties and Non-Trivial
Zeros Study*, arXiv:2503.24275 [**math.GM**], 31 Mar 2025.
Status: **abstract inspected**; full text not retrieved. This preprint distinguishes
"approximate zeros" from "strict zeros" and suggests Spira's off-line points may be the
former. It is in `math.GM`, carries no journal reference, and this extension's own
evidence (previous paragraph: Newton convergence to the precision floor at three
precisions, agreement of three independent evaluation routes at the root, and a local
winding number of exactly 1) is **inconsistent with** the "approximate zero" reading for
the one zero examined here. This is recorded as a *contested point in the literature that
our computation addresses for a single zero*, not as a refutation of that preprint, whose
argument has not been read.

**B4. Independent census.** This extension does not rely on any published zero list. The
counts in `data/dh_zeros.csv` come from the argument principle applied to `F` on a strip
whose edges are **provably** zero-free (`DERIVATION.md` §1), and the on-line ordinates
from sign changes of a Hardy-type `Z`. The published value B1/B2 is compared against the
independently located zero and the two agree; that is a cross-check, not an interval
certificate (non-claim 7).

---

## C. Bombieri–Ghosh

**C1. E. Bombieri and A. Ghosh**, *Around the Davenport–Heilbronn function*,
Russian Math. Surveys **66** (2011), no. 2, 221–270; Uspekhi Mat. Nauk **66** (2011),
no. 2, 15–66. Received 4 October 2010. Record confirmed at
`mathnet.ru/php/archive.phtml?wshow=paper&jrnid=rm&paperid=9410`.
Status: **bibliographic record confirmed; full text not inspected.** Subject matter per
the published abstract: analytic properties of the Dirichlet series of Davenport–Heilbronn
and Titchmarsh, zeros, reciprocals of Dirichlet series, coefficient estimates.
**Nothing in this extension uses a result of theirs**; it is listed because the brief
asked for it and because it is the standard modern reference for the analytic theory of
`F`. In particular, no zero-density or zero-count constant is imported from it: the
counting function used here is the computed argument-principle count, with
`(T/2π)log(5T/2πe)` quoted only as a diagnostic main term.

---

## D. Weil explicit formulas for this function

No published Weil explicit formula specific to the Davenport–Heilbronn function was
located. The formula used here is **derived in `DERIVATION.md` §3** from the Hadamard
factorisation of `Λ_F` and the functional equation, by the standard contour argument, and
is then checked against (i) an independent scalar route with a different representation
of the archimedean kernel (Phase B [B5]) and (ii) a logarithmic-derivative contour
integral that never touches the `b_n` (Phase C). The archimedean closed form is the
parent repository's, reused at parity `a = 1`, which is independently verified there.

---

## E. The repository's own prior work (reused, not re-derived)

**E1.** `weil-window-control`, root `REPORT.md`, `notes/proofs.md`, `src/zeta_window.py`,
at `main = 35c7d95` (tag `dirichlet-windows-v1`). Paper
[doi:10.5281/zenodo.23250932](https://doi.org/10.5281/zenodo.23250932) (version),
concept `10.5281/zenodo.21109955`.

* **Theorem 7** (`\label{thm:zeta}`, "relative perturbation estimate; even sector"):
  `Q_δ(f) ≥ Q_0(f) − B_L(δ)` and `λ_min(Q_δ) ≥ λ*(L) − B_L(δ)`, with
  `B_L(δ) = c(L)δ² + (16/3)L⁵δ⁴e^{2Lδ}`, `c(L) = 8L³(1/3 + 1/√5)`. It is **relative**:
  Retraction 18 of the parent `REPORT.md` withdrew the step that replaced `λ*(L)` by `0`
  outside `L ≤ 0.8`.
* **Theorem 8** (`\label{thm:zetasharp}`, "sharper, conditional; even sector"): defines
  the **remaining-zero form** `Q_rest(f) = Q_0(f) − 4F(γ)² = Q_ζ(f) − 2F(γ)²` and gives
  `λ_min(Q_δ) ≥ −4K(L,γ)δ² − (2/5 + (16/3)e^{2Lδ})L⁵δ⁴`, with
  `K(L,γ) = ∫_{−L}^{L} u² sin²(γu) du` the exact extremal value, attained by
  `f ∝ u sin(γu)`. Its hypothesis — *every nontrivial zero other than the one being moved
  lies on the critical line* — is stated there as **sufficient** for `Q_rest ≥ 0`, and the
  theorem explicitly records that it is *not necessary* and does *not* follow from
  positivity of `Q_ζ`, of `Q_0`, or from a zero-free strip.
* **Corollary 9** (`\label{cor:uniform}`): the `γ`-uniform `−(8L³/3)δ² − O(δ⁴)`, and
  `K(L,γ) → L³/3` as `γ → ∞`.
* **Remark** after Theorem 8: `4K` is sharp *for the derivative-evaluation step only*;
  it is **not** claimed to be the sharp `δ²` coefficient of the optimised minimum
  eigenvalue.
* **Quartet expansion** (`notes/proofs.md` §B.1): three conventions are tabulated; the
  one adopted there is *double on-line pair → one off-line quartet*, perturbation
  `−4δ²(F′² + FF″)`, the only one preserving the zero count. The present extension uses
  the **actual** zero configuration, so it needs no convention choice for `Q_F` itself;
  the convention enters only in Phase G's synthetic deformation, where it is restated.

**E2.** `extensions/dirichlet-windows-2026-10-23` (tag `dirichlet-windows-v1`), its
`REPORT.md`, `README.md`, `notes/dirichlet_window.md`, `src/dirichlet_window.py`,
`src/zeros.py`. Reused here: `arch_closed_a(terms, L, a)` at `a = 1`, the `LogQ`/`LogPi`
split, and the index-set matrix-builder convention. Its **γ = 0 real-pair parity split**
(`+2A_δ(f)²` even, `−2S_δ(f)²` odd) is *recovered independently* here as the `γ → 0`
limit of the quartet form (4.2) and is re-verified in `tests/test_orbits.py`.

**E3. X. Zhu**, *Weil positivity in compact windows: a finite reduction, certified
two-sided bounds, and a Landau–Widom decay law*, arXiv:2608.24827. Status: **not
inspected for this extension.** It is the parent repository's source for the ζ-side
certificate `λ*(0.8) ≥ 8.9 × 10⁻¹⁸`. **It is not used here at all**: it is a statement
about ζ, and no analogue exists for the Davenport–Heilbronn function. The brief's
"Zhu's bandwidth discussion" is therefore recorded as **not verified**, and the
"reference scale" `L_ref = ½ log(5γ₁/2π)` used in Phase G is treated as a *heuristic
scale with unchecked provenance*, never as a predicted onset.

---

## F. The two search leads named in the brief

**F1. "Localized Constrained Weil Forms: Near-Degeneracy and a Davenport–Heilbronn First
Kernel"** — located: **Tariq Ziad**, Zenodo, published 7 October 2026, v1.0, MIT licence,
[doi:10.5281/zenodo.23145455](https://doi.org/10.5281/zenodo.23145455).
**Landing page inspected; PDF not inspected.** Claims, as stated on the record:

> "the same localization framework is applied to the Davenport–Heilbronn function. A
> uniform positive base range is established for small support, together with an exactly
> constrained negative witness at a=9/5."

> "These results imply the existence of a first loss-of-positivity scale
> ½log 2 < b*_DH < 9/5."

Framework: *odd* functions supported in finite symmetric intervals, **constrained**;
`λ₁` via a finite-dimensional head / analytic-tail decomposition with interval
arithmetic, python-flint/Arb plus Lean 4 auxiliary checks.

**F2. "The Davenport–Heilbronn Function as a Negative Control"** — **not located** as a
standalone document. The closest located material is the *Cross-Energy of Dirichlet
Partial Sums* series on Zenodo, in particular **PART V: Anatomy of a Zero, Euler-Product
Interference and the Davenport–Heilbronn Control**
([zenodo 23070453](https://zenodo.org/records/23070453)) and **PART VII: Weil Positivity
as a Prime Cross-Energy** ([zenodo 23070557](https://zenodo.org/records/23070557)).
**Landing pages only; PDFs not inspected.** A sentence surfaced from that series reads:

> "For the Davenport–Heilbronn function, the same form remains positive up to a window
> L ≈ 3.6 and then becomes macroscopically negative, reaching λmin = −0.73 at L = 4, with
> a minimizer concentrated on its off-line zero 0.8085 + 85.70i."

**F3. Status of the reported `L = 1.8` witness — UNVERIFIED.** What was established:
a document exists (F1) claiming a negative witness at `a = 9/5`. What was **not**
established: the normalisation of `a`, which constraints are imposed, whether the
certificate covers the value, and whether `a` is the half-support of the test function or
of its autocorrelation. Note that F1's `a = 9/5` and F2's `L ≈ 3.6` differ by exactly a
factor 2, which is consistent with — but does **not** establish — the two documents using
the two ends of a `g` / `h = g ⋆ g` support convention. In the convention of the present
extension (`g` supported in `[−L, L]`, `h` in `[−2L, 2L]`), either reading places the
reported onset at `L = 1.8` **or** at `L = 3.6`, and nothing here resolves which.

**This extension's scan range `L = 1.0 … 2.6` was fixed by the brief before these
documents were located, and no quantity was tuned to any reported value.** Where our
measured onset is reported, it is reported as "the smallest tested window with a detected
negative witness" in our own stated convention, with the comparison to F1/F2 flagged as
convention-ambiguous.

---

## G. Novelty

**No claim is made that this is the first Weil-form experiment on the Davenport–Heilbronn
function, or the first test of a remaining-zero hypothesis.** Section F shows at least two
contemporaneous documents working on the same control, one of them with an interval
certificate. What this extension can claim as its own is only: an independent
implementation of the classical control inside this repository's conventions, with its own
derivation, its own zero census, and its own gates. Where a result coincides with F1/F2 it
should be read as agreement, and where it differs the convention ambiguity of §F3 must be
resolved before either is called wrong.

---

## H. Access limitations encountered

| source | attempted | result |
|---|---|---|
| Balanzario–Sánchez-Ortiz, Math. Comp. 76 (2007) | AMS PDF | **HTTP 403** |
| Spira, Math. Comp. 63 (1994) | — | paywalled, not attempted beyond the index |
| Davenport–Heilbronn 1936 | — | no accessible full text |
| Bombieri–Ghosh 2011 | mathnet.ru record | record only; article paywalled |
| Liu–Wu, arXiv:2503.24275 | arXiv abstract | abstract only |
| Ziad, zenodo 23145455 | Zenodo landing page | metadata and abstract only |
| Cross-Energy PART V / VII | Zenodo landing pages | metadata only |
| Zhu, arXiv:2608.24827 | — | **deliberately not used**; see E3 |

No third-party PDF was downloaded, so `MANIFEST.sha256` contains hashes of this
extension's own files only.

---

## I. Additional references cited by the note (`note/dh_note.tex`)

Added when the note was written. Each is a **standard reference cited for context or for
a tool**, not a source of any number used in a computation here.

| reference | status | what it is cited for |
|---|---|---|
| E. C. Titchmarsh, *The Theory of the Riemann Zeta-Function*, 2nd ed., revised by D. R. Heath-Brown, Oxford University Press, **1986**, **§10.25** | **not inspected** | the standard textbook exposition of the Davenport–Heilbronn function. Nothing is taken from it; the section number is the one conventionally cited for this function. |
| H. Yoshida, *On Hermitian forms attached to zeta functions*, Adv. Stud. Pure Math. **21**, Math. Soc. Japan, Kinokuniya, Tokyo, **1992**, 281–325 | **not inspected** | prior art for the finite-window Weil form, as recorded in the parent paper's bibliography, from which this entry is copied. |
| F. Johansson, *Arb: Efficient arbitrary-precision midpoint-radius interval arithmetic*, IEEE Trans. Comput. **66** (**2017**), no. 8, **1281–1292** | **not inspected** | the citation for the ball-arithmetic library used by `src/certify.py` (python-flint 0.9.0 wraps Arb). |
| R. Spira, *Some zeros of the Titchmarsh counterexample*, Math. Comp. **63** (1994), no. **208**, 747–748 | **not inspected** (see §B1) | the issue number is added here to complete the record begun in §B1. |
| E. Bombieri, *Remarks on Weil's quadratic functional…*, Lincei (9) Mat. Appl. **11** (2000), no. 3, 183–233 | **not inspected** | prior art for the finite-window Weil form; entry copied from the parent paper's bibliography. |
| E. De Jesús, the parent paper, [doi:10.5281/zenodo.23250932](https://doi.org/10.5281/zenodo.23250932) | **this repository** | the companion paper, whose §4–5 fix the window conventions and whose remaining-zero theorem is audited in the note's §5. |

No numerical value in the note comes from any of these.
