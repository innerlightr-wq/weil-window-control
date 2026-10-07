# REVISION_CHANGELOG.md — 2026-10-07, Retraction 18 revision

Targeted correction and reconciliation driven by the audit in this directory. Nothing
committed or pushed. Title, authorship, DOI, licence, theorem/corollary/remark numbering,
equation numbering, established results and all numerical tables are preserved — verified
mechanically against the original PDF.

## `paper/weil_window_control.tex`

| # | Location | Change |
|---|---|---|
| 1 | Abstract | ζ bound described as **relative** to `λ*(L)`, "unconditional where that positivity is certified (`L ≤ 0.8`)". |
| 2 | Non-claims 5 | Was "Theorems 7–8 rely on it". Now: Theorem 7 uses Zhu's certificate to convert the relative estimate into an unconditional bound for `L ≤ 0.8`; **Theorem 8 does not use it**, it rests on its own hypothesis. |
| 3 | §4.2, new paragraph before Thm 7 | **Defines what was previously only a macro.** `f` admissible; `Q_ζ` the actual ζ form; `Q_0 = Q_ζ + 2F(γ)²`; `b_0(L,γ) = inf Q_0`; `λ*(L) ≤ inf Q_ζ`; budget `B_L(δ) = c(L)δ² + (16/3)L⁵δ⁴e^{2Lδ}`. States that `b_0` and `λ*` are different objects, never identified; that both are **quadratic-form infima**, with no assumption that the Weil form is bounded on `L²`; and that a lower bound proved on a **larger test space is a fortiori** one on the even sector, so a certificate may be quoted unchanged — with its documented normalisation (`Q(f)/‖f‖²`) and range only. |
| 4 | **Theorem 7** | Retitled *relative perturbation estimate*. Chain now explicit: `Q_δ(f) ≥ Q_0(f) − B_L(δ)` ⟹ `inf Q_δ ≥ b_0(L,γ) − B_L(δ)` and `λ_min(Q_δ) ≥ λ*(L) − B_L(δ)` (the latter because `Q_0 ≥ Q_ζ`). |
| 5 | After Theorem 7 | **The offending sentence is gone.** Replaced by the wording specified in the task, plus the explicit `L ≤ 0.8` instance `λ_min(Q_δ) ≥ 8.9×10⁻¹⁸ − B_L(δ)`. |
| 6 | After Theorem 7 | New *unnumbered* paragraph "Validity versus usefulness of the budget": the remainder is a convergent series tail, so `B_L(δ)` is valid at every `δ`; only its informativeness degrades as `Lδ` grows. Unnumbered **so that no theorem number shifts**. |
| 7 | **Theorem 8** | Remaining-zero form displayed: `Q_rest(f) = Q_0(f) − 4F(γ)² = Q_ζ(f) − 2F(γ)²`, written `Q̃`. The hypothesis is kept exactly as it was and labelled a **sufficient** input where the explicit formula applies, **not** an equivalence; and `Q̃ ≥ 0` is stated not to follow from positivity of `Q_ζ`, of `Q_0`, or from a zero-free strip. |
| 8 | Remark 10 | "sharp" qualified as **sharp for the derivative-evaluation step**, with (a)/(b)/(c) separated: `K` is the exact squared norm of the functional; `f ∝ u sin(γu)` attains it; neither implies the minimiser of the entire form attains the same coefficient. Near-saturation kept as T2. |
| 9 | Limitations 1–2 | Item 1 rewritten (relative estimate; substitution not asserted beyond `L ≤ 0.8`, cross-referencing Retraction 18). Item 2 extended with the `Q̃` non-implication. |
| 10 | Limitations, new item | Scoped audit statement: for the audited source version and the transfer attempted, no usable estimate replaces the remaining-zero positivity hypothesis; the per-zero bound available is uniform in the ordinate and so not summable **as it stands** — "a limitation of that particular estimate, not a general obstruction". |
| 11 | Appendix A | **Retraction 18** appended after 17; count "Seventeen" → "Eighteen". No earlier entry renumbered or removed. |

**Not changed:** Theorem 6 and the function-field section; every numerical table; the
dictionary; §§3, 5, 6, 7; the `K(L,γ)` closed form (brackets and signs verified against the
TeX **and** the rendered page image); the `δ²` expansion; `c(L)`; the remainder bound.

## `notes/proofs.md`

- §B.4 now defines `Q_ζ`, `b_0`, `λ*`, `B_L(δ)` before the theorem, with the same separation.
- **Theorem B restated** as the relative estimate; the chain stops at `λ_min(Q_δ) ≥ λ*(L) − B_L(δ)`. The previous continuation `… ≥ −B_L(δ)` is **withdrawn in place**, with a note saying so.
- The "For `L > 0.8` … uses no positivity at all" paragraph is replaced by "**For `L > 0.8` nothing is asserted**", the `δ = 0` circularity spelled out, the a-fortiori test-space remark, and the certificate's normalisation and range.
- **Logical regression check added**: `X ≥ b − D` permits `b → 0` only given a separate `b ≥ 0`; scalar control `b = −1, D = 0, X = −1`, with an explicit note that this is about scalars and is **not** an example of a negative actual Weil form.
- §B.5: `λ_min(Q_δ) → λ*(L) > 0` softened to `→ λ_min(Q_0)`, positive where certified and of unknown sign beyond.

## `REPORT.md`

- Task 3(b): the withdrawn claim marked withdrawn, certificate normalisation added, relative estimate stated.
- Task 3(c): baseline sign qualified identically.
- **Retraction log entry 18** appended after 17 (with the required sentence verbatim), recording cause, what replaces it, blast radius, and that the audit found it.

## `README.md`

- Results row labelled *relative*, with the certified range.
- Non-claim 5 extended to separate the certified **lower** bound (Zhu, `L ≤ 0.8`, `Q(f)/‖f‖²`) from the variational **upper** bound, and to record that the substitution is not asserted.
- Extensions: entry for this audit directory, with the scoped NO BRIDGE FOUND wording and the **exact** formalization status (comparator stub inspected; implementation declaration located via `formalization.yaml`; full Lean proof and dependency chain **NOT** checked, no build run).

## New files in this directory

`REVISION_CHANGELOG.md`, `validate_revision.py`, `VALIDATION.md`. The four completed audit
documents and the two audit scripts are unchanged.

## Deposit and DOI reconciliation (added after the Zenodo upload)

The revision was deposited on 2026-10-07 as **`10.5281/zenodo.23212998`**, file
`weil_window_control_revised_2026-10-07.pdf`, CC BY 4.0. The deposited file is **byte-identical
to the build produced here** (483,883 B, md5 `a29f9670844f341c7986cfc1d4382bbe`), verified
against the Zenodo API. Concept DOI for all versions: **`10.5281/zenodo.21109955`**.

- `paper/weil_window_control.pdf` now **is** that deposited artifact (same bytes), so the
  committed source and the committed build agree and `make paper` reproduces the content.
  The separate dated copy was removed as a duplicate; the pre-revision PDF remains in git
  history at `f6033d7` and on Zenodo as `10.5281/zenodo.23195177`.
- Version DOI updated to `…23212998` in `README.md`, `CITATION.cff` and
  `extensions/comb-dips/README.md`; the concept DOI `…21109955` added to the README headline
  and as a second `identifiers` entry in `CITATION.cff`; `CITATION.cff` `version` 1.0.0 →
  1.1.0 and `date-released` → 2026-10-07.

### One known cosmetic mismatch, left deliberately

**The PDF's title page still reads `doi:10.5281/zenodo.23195177`** (`paper/weil_window_control.tex:38`)
while the file is published under `…23212998`. Changing it would require a rebuild, which would
break byte-identity with the deposited artifact. The title page was therefore left untouched and
the mismatch is recorded here instead. **Suggested fix at the next revision:** put the *concept*
DOI `10.5281/zenodo.21109955` on the title page, so it never goes stale again.

*(Separately: the v1 Zenodo record `…23195177` holds the 477,046-byte macOS-Preview re-save of
the original PDF, not the 474,032-byte pdfTeX build that was committed. Identical text layer,
different container — noted for the record, no action needed.)*
