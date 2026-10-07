# SOURCES_AND_SCOPE.md

Audit 2026-10-07. New files live only in this directory. Both repositories are treated as
read-only; nothing in either paper, PDF, README, citation file, old report or memory file was
changed.

> **UPDATED 2026-10-07 (same day, later pass).** One source was added and one conclusion narrowed;
> see [`SEMILOCAL_LITERATURE_UPDATE.md`](SEMILOCAL_LITERATURE_UPDATE.md). Semi-local Sonin spaces
> and the relevant structural identifications are established in later work (CCM arXiv:2310.18423v2,
> 2024, added to the source table below). The audit did not locate the specific positive-trace
> comparison with sufficient discrepancy control needed for our application. **This is a gap in the
> verified transfer, not an absence of semi-local geometry.** A citation detail about the Ramanujan
> deposit is also corrected below.

## Repository state (before = after)

| | HEAD | branch | status |
|---|---|---|---|
| `weil-window-control` | `13ad9dd` | `main` | clean, 0 modified, 0 staged, 0 ahead |
| `thales-ramanujan-quarter-gamma` | `14242bc` | `main` | clean, 0 modified, 0 staged |

Both checkouts already existed at the expected paths; nothing was cloned, pulled, reset or
switched.

## Primary sources — what was actually read

| source | status |
|---|---|
| **Connes–Consani, *Weil positivity and Trace formula, the archimedean place*, author-hosted PDF dated **July 4, 2021** (`alainconnes.org/.../Selecta.pdf`, 55 pp, 935,906 B)** | **READ.** This is the version used for every numbered citation below. Introduction and Theorem 1 (p. 2); Theorem 4.6 (p. 24, eq. 84–85); Lemma 6.9 (p. 44, eq. 139–140); Lemma 6.10 (p. 44, eq. 141); Theorem 6.11 (p. 45, eq. 142); Remark 6.12 (p. 45). |
| arXiv:2006.13771 (57 pp) | **Downloaded, NOT used for numbering.** Numbering and hypotheses were **not** transferred between versions; every citation here is to the July 2021 PDF. |
| **Connes–Consani, *Quasi-inner functions and local factors*, arXiv:2008.10974v1 (25 pp)** | **READ** — abstract, Introduction §1, and the two displayed Theorems. |
| **Connes–Consani–Moscovici, *Zeta zeros and prolate wave operators* / *Semilocal adelic operators*, arXiv:2310.18423v2, 4 May 2024 (30 pp)** — *added 2026-10-07* | **READ** — abstract, §1 Introduction, §4 (incl. §§4.6–4.8), reference list. Cited: Def. 4.5; Prop. 4.6(i)–(ii) with eqs. (57)–(58); Prop. 4.7(i)–(iii) with diagram (59); **CCM2024-Theorem 4.6, p. 23 = introductory Theorem 2, p. 5**; footnote 2, p. 5; §4.8; diagram (62). §§2–3, 5–6 skimmed only — **NOT CHECKED in detail**. |
| **Guillera, arXiv:1911.03968 (15 pp)** | **Downloaded.** Not load-bearing: the Ramanujan comparison is used only as a mechanism table (`MECHANISM_DICTIONARY.md`) and no identity from it is transferred. Treated as **NOT CHECKED in detail**. |
| **De Jesús, *Ramanujan's Formula for 1/π and the Thales Integral Witness*** — version DOI `10.5281/zenodo.20672985` (concept `…20665588`), 8 pp, fetched from Zenodo | **READ.** The repository itself contains no manuscript; it points at this deposit. **Citation correction 2026-10-07:** a local copy (`~/Downloads/QuarterGammaRamanujan (1).pdf`, md5 `517e779cef517dc177f51056117829d5`, three identical copies present) is **byte-identical** to the Zenodo PDF already read, but its **title page** reads *"SHARED QUARTER-GAMMA LATTICE AND THE LIMITS OF THE THALES–RAMANUJAN ANALOGY"*, not the Zenodo metadata title quoted here. The Zenodo metadata title and the PDF title differ; the DOI is unambiguous and the content is unchanged, so **nothing in this audit depends on which title is used** — but citations should give both or cite the DOI. |
| `thales-ramanujan-quarter-gamma` repo (36 files: `docs/`, `scripts/`, `results/`, `thales_ramanujan/`) | **READ** — `README.md` and the documented scope. |
| **Local Weil material at `13ad9dd`** | **READ** — `paper/weil_window_control.tex` (Theorems 7–8 as corrected, Remark 10), `src/zeta_window.py` (the form assembly), and `extensions/second-moment-generalization-2026-10-07/` including `REPORT.md`, `PROJECTION_CONTEXT.md`, `PROJECTION_CONTEXT_CHECKS.md`. |
| `extensions/openai-quasi-rh-audit-2026-10-07/` | **PRESENT and read previously**; not re-run here. Its verdict (NO BRIDGE FOUND, scoped) is carried over, not revisited. |
| Nested-window / coupling and prime-comb investigations | `extensions/comb-dips/` is present (verdict NULL, carried over). The nested-window inheritance result is recorded only inside `PROJECTION_CONTEXT.md`; **no separate report exists**, so it is cited as recorded there and not reconstructed. |
| Connes–Consani semi-local trace formula (their ref. [2]) and Yoshida's original | **NOT CHECKED** — not retrieved. Used only as cited by the sources above. *(2026-10-07: CCM2024 cites the semi-local trace formula as its ref. **[7] = A. Connes, Selecta Math. (N.S.) 5 (1999), 29–106**, and states that in it "the Weil quadratic form `Q_n` becomes the trace of a simple operator theoretic expression". Still **NOT CHECKED** directly; reported as cited.)* |
| "New eigenfunctions for the negative part of the Connes–Moscovici prolate spectrum" (2025 item surfaced in search) | **NOT CHECKED** — not retrieved; concerns eigenfunctions, not the positive-trace comparison. Recorded as a pending lead only. |
| The forthcoming paper promised in the CCM2024 abstract (a "second candidate for the semilocal prolate operator") | **NOT LOCATED** as of 2026-10-07. |
| Kato; Zhu's certificate in the original | **NOT CHECKED** — Zhu is used exactly as the manuscript already uses it. |

## Load-bearing statements, with hypotheses

All from the **July 4, 2021** PDF.

**Theorem 1 (p. 2).** *Let `g ∈ C_c^∞(R*_+)` have support in `[2^{-1/2}, 2^{1/2}]` and Fourier
transform vanishing at `i/2` and `0`. Then* `W_∞(g * g*) ≥ Tr(ϑ(g) S ϑ(g)*)`.
Here `ϑ(λ)ξ(v) = λ^{-1/2}ξ(λ^{-1}v)` on `L²(R)_ev`; `S = Π_{S(1,1)}` projects onto the Sonin space
of even functions vanishing on `[−1,1]` together with their Fourier transform; `W_∞ := −W_R`.
**Two vanishing conditions**: at `i/2` *and* at `0`.

**Theorem 4.6 (p. 24).** *For all* `f ∈ C_c^∞(R*_+)`:
`Tr(ϑ(f)S) = W_∞(f) + ∫ f(ρ^{-1}) ε(ρ) d*ρ` (84), with `ε` given by (85) and `ε(ρ^{-1}) = ε(ρ)`.
**No support restriction.** This is the exact identity; the discrepancy is `E(f) := ∫f(ρ^{-1})ε(ρ)d*ρ`,
defined independently and not by subtraction.

**Lemma 6.9 (p. 44).** For unit vectors `φ, ψ` and `a,b,c ≥ 0`,
`B(ξ) = −b|⟨φ|ξ⟩|² + a|⟨ψ|ξ⟩|² + c‖P_φ ξ‖²` is positive **iff**
`a + c ≥ b` **and** `b(a+c) ≤ a(b+c)|⟨φ|ψ⟩|²` (139); then `B(ξ) ≥ ε‖ξ‖²` with `ε` from (140).

**Lemma 6.10 (p. 44).** With `N_I = −2ℓ₁(1+ε₁)(Id − K_I)` on `H = L²(I)`,
`I = [−½log2, ½log2]`, and `γ ≈ 2.94355`: `⟨ξ|N_I ξ⟩ ≤ γ|⟨ξ₀|ξ⟩|²` (141).

**Theorem 6.11 (p. 45).** *`g` supported in `[2^{-1/2},2^{1/2}]`, Fourier transform vanishing at
`−i/2`* (only one condition here). Then
`W_∞(g*g*) ≥ Tr(ϑ(g)Sϑ(g)*) − c|ĝ(0)|²`, `c = 4γ/log 2` (142).
Theorem 1 is the case `ĝ(0) = 0`.

**Remark 6.12 (p. 45).** `λ_max = 1.05158`, `ε₁ ≈ 0.00122`, quoted; **not independently derived
here**.

**arXiv:2008.10974 — what is and is not established.** Established: the product
`u = ρ_∞ ∏ ρ_p` over a finite set of places containing `∞` is quasi-inner (Theorem, §1); the
semi-local Sonin space `S(u(F))` is infinite-dimensional and the spaces form a filtering system
(second Theorem, §1). **Not established**: the paper states it is "meant to be a **first test**
pertaining to the general **strategy** proposed in [5]", *defines* the semi-local Sonin space as
`ker U₂₂`, and explicitly **postpones to a future paper** the proof that this definition
reproduces the semi-local analogue. ~~**There is no semi-local analogue of Theorem 1 / 6.11 — no
semi-local positive trace with a controlled discrepancy against the semi-local Weil
distribution.** This is decisive for `POSITIVITY_DEFECT_AUDIT.md` §4.~~

> **Narrowed 2026-10-07.** The struck sentence was correct about *this* 2020 paper but was written
> as a statement about the literature, which it is not. The postponed geometry appeared in
> **CCM arXiv:2310.18423v2 (2024)**, summarised next. The surviving claim is narrower and is what
> `POSITIVITY_DEFECT_AUDIT.md` §4 now rests on: **no semi-local analogue of CC2021 Theorem 4.6 or
> Theorem 6.11 — i.e. no independently defined semi-local discrepancy `E_S` and no bound on it —
> was found in the checked sources.**

**CCM arXiv:2310.18423v2 (2024) — what is and is not established.** Established: the semi-local
Sonin space is **defined directly**, `S_λ(X_S,α) := {f ∈ L²(X_S)^{K_S} | f(x) = 0 & F_S f(x) = 0 ∀x,
|x| < λ}` (Def. 4.5); `θ_S(S_λ(R,e_∞)) ⊂ S_λ(X_S,α)` (Prop. 4.6(i)); the Mellin multiplier
`F_μ w_S(θ_S f)(s) = ∏_{p∈S\{∞}} L_p(½+is)^{-1}(F_μ w_∞ f)(s)` (Prop. 4.6(ii), eqs. 57–58);
`F_S ∘ θ_S = θ_S ∘ F̃_R` and `⟨θ_S(f)|η_S(g)⟩ = ⟨f|g⟩` (Prop. 4.7(i),(iii)); and **`θ_S` is a
"hilbertian isomorphism" of Sonin spaces, stable under enlarging `S`** (CCM2024-Thm 4.6, p. 23 =
intro Thm 2, p. 5). **Hypothesis to note:** "hilbertian" is defined in their **footnote 2, p. 5** as
*"the underlying topological vector space structure of a Hilbert space"* — so the map is bounded and
invertible, **not** isometric, and §4.8 states that *"the choice of the finite set `S` plays a key
role in fixing the inner product"* and that *"`B_λ` inherits different inner products from its
embedding in `L²(R, ds/|E_S(s)|²)`."* **Not established**: the semi-local positive-trace comparison.
CCM2024 §1 says only *"We **expect** that the use of such operator-theoretic tools in the semilocal
case **opens a way** to handle Weil's positivity as in [10]"* and that the paper *"provides a more
precise **strategy** for addressing the semilocal Weil positivity"*; the abstract defers a second
candidate prolate operator to *"a forthcoming paper"*. The archimedean identification itself holds
*"up to a finite dimensional possible discrepancy"* (§1).

## Prior art, explicitly

**The positive Sonin trace `Tr(ϑ(g)Sϑ(g)*) ≥ 0` and the alignment criterion (Lemma 6.9) are
Connes–Consani's.** Nothing in this directory discovers either. What is done here is to map them
onto the repository's `Q_ζ` with a derived dictionary, and to say precisely where the arithmetic
step is blocked.

## Extraction hygiene

`pdftotext -layout` was used for navigation, but every displayed formula quoted above was read in
the layout output with its brackets and conjugation bars checked against surrounding context;
where the extraction mangled a display (e.g. eq. (1), whose line order is scrambled), the formula
is **not** quoted verbatim and is described instead. This repository has twice been bitten by
extraction artefacts (a dropped `\left[` and a dropped `Decimal` exponent), both recorded in
`extensions/second-moment-generalization-2026-10-07/`.
