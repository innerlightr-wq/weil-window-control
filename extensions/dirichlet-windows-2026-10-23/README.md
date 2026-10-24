# Dirichlet windows — the finite-window Weil detector on `L(s, χ)`

Extension of [`weil-window-control`](../../README.md). **Carries no DOI, is not part of the
paper, and has not been peer-reviewed.** Full lab record, non-claims and retraction log:
[`REPORT.md`](REPORT.md). Derivations: [`notes/`](notes/).

---

## The result

For a real primitive character `χ` mod `q`, a **real zero** `β = ½ + δ` of `L(s, χ)` forces
its partner `1 − β = ½ − δ`, so in the variable `z` with `ρ = ½ + iz` the configuration is
the pair `z = ∓iδ` — a pair of *purely imaginary* ordinates, which is exactly `γ = 0` in the
off-line quartet parametrisation. Its contribution to the finite-window Weil form is
`2F(iδ)F(−iδ)` with `F(iδ) = ∫ f(x) e^{−δx} dx`, and because `F` inherits the parity of `f`
this splits by **test-function parity**:

```
f even :  + 2 A_δ(f)²,   A_δ = ∫ f(u) cosh(δu) du       rank one, POSITIVE semidefinite
f odd  :  − 2 S_δ(f)²,   S_δ = ∫ f(u) sinh(δu) du       rank one, NEGATIVE semidefinite
```

Consequently, if every non-trivial zero of `L(s, χ)` lies on the critical line except one
real pair `½ ± δ`, then `Q_χ(f) = Σ_{on-line} |F(γ_ρ)|² + 2A_δ(f)² ≥ 0` for **every** real
even `f` and **every** window half-length `L`. So **even-sector finite-window Weil positivity
is compatible with a real zero at any `δ ∈ (0, ½)`, and no even-sector positivity computation
— at any `L`, in any basis, to any precision — carries information about the Siegel
configuration.** This is strictly stronger than the observation that the `δ²` coefficient
`4K(L, 0)` vanishes: the response is a perfect square of the wrong sign, exactly in `δ`, not
a leading-order degeneracy. The odd sector is the complementary case and is the point of
**maximum** sensitivity at `γ = 0`, with constant `K_odd(L, 0) = 2L³/3`, because
`K_even(L, γ) + K_odd(L, γ) = 2L³/3` identically in `γ`.

The sign attachment is not a bookkeeping artefact: the same signed rank-one split, with the
same attachment of sign to parity, holds **exactly over ℚ** in the function-field control
where the Riemann hypothesis is Weil's theorem ([`notes/ff_parity.md`](notes/ff_parity.md)),
and at `δ = ½` the even- and odd-sector real-pair perturbations reproduce this repository's
ζ **pole-term** matrix entrywise, `pole_sign(ODD) = −1` included.

**Not claimed:** anything about GRH, RH, or the existence of Siegel zeros; any statement
about `L(s, χ)` as opposed to the detector. The odd-sector sensitivity is a *detector*
measurement and not an exclusion, because converting it would require an unconditional lower
bound on the remaining-zero form and no Zhu-style certificate is known for any Dirichlet
window. The October 2026 OpenAI quasi-Riemann-hypothesis release is **not used as input
anywhere**. See [`REPORT.md`](REPORT.md) §0.

---

## The note

A standalone 12-page note, [`note/siegel_note.pdf`](note/siegel_note.pdf),
[doi:10.5281/zenodo.23624161](https://doi.org/10.5281/zenodo.23624161) (version 2), states the one
result with its non-claims, tiers and retraction log:
**at a real zero the even sector of the finite-window Weil form cannot respond at all,
while the odd sector responds with its largest possible negative value.**
`make -C note verify` re-derives every number printed in it from this extension's data
files, and fails on any disagreement, any undeclared literal, any overstatement of
precision, or an unmet prerequisite.

### Erratum to version 1 — read this before citing it

**Version 1** ([doi:10.5281/zenodo.23623473](https://doi.org/10.5281/zenodo.23623473),
10 pages, deposited 2026-10-24) contains one statement that is **withdrawn**. Its §5 said
that if Zhu's certified infimum were over the full admissible space then it *would* constrain
a real zero. Zhu's paper was obtained the same day: the infimum **is** over the full space
(Corollary 6.3, no parity or reality restriction), and the consequence does **not** follow —
the certified quantity is the sum over *all* zeros, so a real pair's negative odd-sector
contribution is absorbed by a remaining-zero sum that nothing bounds. §5 also named a *lower*
bound on the remaining-zero form as the missing ingredient for an exclusion, where it is an
*upper* bound. Details in [`REPORT.md`](REPORT.md) §0.4.

**Version 2** ([doi:10.5281/zenodo.23624161](https://doi.org/10.5281/zenodo.23624161), 12 pages)
corrects both, cites Zhu accurately, and records that his Lemma 6.1
is the `δ = ½` instance of this note's theorem — derived independently and in print first.
**The theorem and both corollaries are unchanged**: they describe the response of the form to
a planted real zero and never depended on any certificate's domain. Version 2 additionally
attributes the sign-by-parity phenomenon to Bombieri (2000) §9, Theorem 9, which version 1
presented without attribution; see [`REPORT.md`](REPORT.md) §0.5.

---

## Reproducing the verification steps

Environment — the parent repository's pinned `mpmath`:

```bash
cd <repository root>
python3 -m venv .venv
.venv/bin/pip install 'mpmath==1.3.0'
```

| | version used for the committed numbers |
|---|---|
| Python | **3.14.7** (CPython) |
| mpmath | **1.3.0** (the version `requirements.txt` pins) |
| numpy / matplotlib | **not required** by this extension |
| gmpy2 | not installed; it changes runtime only, no result |

The ζ-side values in [`results/regression_zeta.txt`](results/regression_zeta.txt) reproduce
this repository's committed `data/stage3_gateB_lambda_min.csv`,
`data/stage3c_even_odd.csv` and `data/stage3_C_vs_CS.csv` to **all printed digits** under
this newer CPython, which is the check that the interpreter version changes nothing.

### Step 1 — the three Siegel-corner checks (~25 min)

```bash
.venv/bin/python extensions/dirichlet-windows-2026-10-23/src/verify_siegel_parity.py
```

(a) decides **from the committed Phase-5 data**, not from a reading of the code, whether the
planted object was the real pair or the degenerate quartet, by recomputing stored
`lam_planted` values both ways; (b) asserts rank one, the correct sign, and
eigenvalue `= +2A²` / `= −2S²` with `A` and `S` obtained by **quadrature** on the explicit
normalised eigenfunction — an independent check of the closed forms in
`planting.F_at_imag` — plus `K_even + K_odd = 2L³/3` with `K_odd` also by quadrature;
(c) compares the `δ = ½` real-pair perturbation entrywise with the parent pole term in all
three sectors, reporting **both** `B − Pole` and `B + Pole` so the data fixes the sign.

### Step 2 — the function-field cross-check (seconds, exact over ℚ)

```bash
.venv/bin/python extensions/dirichlet-windows-2026-10-23/src/ff_parity_split.py
```

Verifies (FF1)–(FF4) of [`notes/ff_parity.md`](notes/ff_parity.md) in `fractions.Fraction`
for `ρ ∈ {2, 3/2, 5/4, 7/3, 1}` and `R ∈ {1, 2, 4, 8, 16}`. No floating point.

### The whole suite

```bash
./extensions/dirichlet-windows-2026-10-23/run_suite.sh
```

Runs, in order: the Phase-1 self-test, the planting-block regressions, Step 1, Step 2, the
Phase-4d crossover criterion, and the 18-value ζ regression suite. Each script exits
non-zero on failure and the runner propagates the first failure. Phase-by-phase
reproduction of the *measurements* (as opposed to the checks) is in
[`REPRODUCE.md`](REPRODUCE.md).

---

## Layout

```
note/          siegel_note.tex -> siegel_note.pdf, its generated tables, the verifier
REPORT.md      lab record: non-claims first, tiers on every claim, retraction log
REPRODUCE.md   how to regenerate every measurement, with dependency order
run_suite.sh   the test suite
notes/         dirichlet_window.md (the form, Phase 1) · ff_parity.md (Step 2)
src/           18 scripts; see REPRODUCE.md for which produces which data file
data/          CSV results and the computed critical-line ordinates
results/       console transcripts and JSON summaries of the runs that produced them
```

## Cite as

> E. De Jesús, *Dirichlet windows: the finite-window Weil detector on `L(s, χ)`*,
> in `weil-window-control`, tag **`dirichlet-windows-v1`** (2026-10-23), which merges
> branch `dirichlet-windows-2026-10-23` into `main` at merge commit
> **`ed5c298fbe2182bd0f92a410440c159145eb3d3b`**.

The tag is the citable object. It is placed one commit after the merge, on the commit that
fills in the merge hash above — a file cannot contain the hash of the commit that adds it,
and a tag whose citation line still read `PLACEHOLDER` would not be citable. Both objects
resolve: `git show dirichlet-windows-v1` for the state, `git show ed5c298fbe21` for the merge.

The note is a separate citable object:

> E. De Jesús, *Parity blindness of finite-window Weil detectors at a real zero*,
> technical note, Zenodo, CC BY 4.0.
> v1: [doi:10.5281/zenodo.23623473](https://doi.org/10.5281/zenodo.23623473),
> published 2026-10-24, repository state tag **`siegel-note-v1`** — **see the erratum above**.
> v2: [doi:10.5281/zenodo.23624161](https://doi.org/10.5281/zenodo.23624161),
> published 2026-10-24, repository state tag **`siegel-note-v2`**.

Cite **v2**. Both are versions of one Zenodo record, so the concept DOI covers the pair and
v1 routes forward to v2; v1 remains resolvable, as a deposited version must.

The parent paper is cited separately:
[doi:10.5281/zenodo.23250932](https://doi.org/10.5281/zenodo.23250932) (version),
[10.5281/zenodo.21109955](https://doi.org/10.5281/zenodo.21109955) (concept).

## License

As the parent repository: **MIT** for code under `src/` and `run_suite.sh`,
**CC BY 4.0** for `REPORT.md`, `REPRODUCE.md`, `README.md`, `notes/`, `data/` and
`results/`.
