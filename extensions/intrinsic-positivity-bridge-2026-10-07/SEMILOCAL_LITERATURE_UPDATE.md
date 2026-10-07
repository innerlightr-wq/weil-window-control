# SEMILOCAL_LITERATURE_UPDATE.md

**Dated 2026-10-07.** Source-reconciliation update to the audit in this directory. No numerical
work was rerun; no manuscript, script, data file or other repository content was touched.

## Verdict: **GEOMETRY UPDATED / ESTIMATE GAP REMAINS**

The earlier availability claim was **too strong** and is amended below. No new usable bound
follows.

---

## 1. The original overbroad statement

`POSITIVITY_DEFECT_AUDIT.md` §4 and `RESULT.md` §2 said, on the strength of arXiv:2008.10974 alone:

> "**The semi-local `{∞, 2}` setting.** arXiv:2008.10974 … **defines** the semi-local Sonin space
> as `ker U₂₂` and states explicitly that it is 'meant to be a **first test** …', **postponing to a
> future paper** even the proof that this definition reproduces the semi-local analogue. **There is
> no semi-local analogue of Theorem 4.6 or Theorem 6.11…**"

That was accurate *about the 2020 paper* but wrong as a statement about the literature: **the
postponed geometric identification was subsequently carried out.** The promised paper exists.

## 2. What the additional source establishes

**C. A. Connes, C. Consani, H. Moscovici, *Zeta zeros and prolate wave operators* — subtitle
*Semilocal adelic operators*, arXiv:2310.18423v2, 4 May 2024, 30 pp.** READ: abstract, §1
Introduction, §4 (incl. §§4.6–4.8), reference list.

> **Naming convention used throughout, as the brief requires:** **CC2021-Thm4.6** is Theorem 4.6 of
> the July 2021 archimedean paper (the trace identity `Tr(ϑ(f)S) = W_∞(f) + E(f)`).
> **CCM2024-Thm4.6** is Theorem 4.6 of this 2024 paper (printed p. 23), a completely different
> statement. They are **not** related.

Verbatim, the abstract: *"We prove the **stability of the semilocal Sonin space under the increase
of the finite set of places** which govern the semilocal framework and describe their relation with
Hilbert spaces of entire functions."*

- **Definition 4.5** gives the semi-local Sonin space **directly**, not as `ker U₂₂`:
  `S_λ(X_S,α) := {f ∈ L²(X_S)^{K_S} | f(x) = 0 & F_S f(x) = 0 ∀x, |x| < λ}`.
- **Proposition 4.6(i)**: `θ_S(S_λ(R,e_∞)) ⊂ S_λ(X_S,α)`.
- **Proposition 4.6(ii), eq. (57)–(58)**: `F_μ w_S(θ_S f)(s) = ∏_{p∈S\{∞}} L_p(½+is)^{-1}·(F_μ w_∞ f)(s)`.
- **Proposition 4.7**: (i) `F_S ∘ θ_S = θ_S ∘ F̃_R`; (ii) the commutative diagram (59) with the two
  **different** measures `ds/|E_∞(s)|²` and `ds/|E_S(s)|²`; (iii) the pairing
  **`⟨θ_S(f) | η_S(g)⟩ = ⟨f | g⟩`**.
- **CCM2024-Thm4.6 (p. 23)** = introductory **Theorem 2 (p. 5)**: *"the map `θ_S` is a **hilbertian
  isomorphism** of the Sonin spaces `θ_S : S_λ(R,e_∞) → S_λ(X_S,α)`."*
- **§4.8**: *"the choice of the finite set `S` plays a key role in **fixing the inner product** in
  the Hilbertian space"*; and *"`B_λ` inherits **different inner products** from its embedding in
  `L²(R, ds/|E_S(s)|²)`."*

**On the general semi-local trace formula.** CCM2024 §1 states: *"the semilocal trace formula
of [7] gives, for each `n`, a Hilbert space theoretic framework in which **the Weil quadratic form
`Q_n` becomes the trace of a simple operator theoretic expression**"*, where **[7] = A. Connes,
*Trace formula in noncommutative geometry and the zeros of the Riemann zeta function*, Selecta
Math. (N.S.) 5 (1999), 29–106**. That source is **NOT CHECKED** here (not retrieved); but its
existence and content are asserted by CCM2024, which is enough to retire the earlier phrasing.

**On the positivity comparison — still a programme.** CCM2024 §1, verbatim:

> *"We **expect** that the use of such operator-theoretic tools in the semilocal case **opens a way
> to handle Weil's positivity as in [10]**. In fact, the operator theoretic aspect of the present
> paper **provides a more precise strategy for addressing the semilocal Weil positivity** by
> comparing … the trace functional associated to the operator — which is automatically positive for
> a selfadjoint operator — with the Weil functional."*

"Expect", "opens a way", "strategy": **proposed, not proved.** Likewise the abstract closes with
*"with the goal of obtaining (**in a forthcoming paper**) a second candidate for the semilocal
prolate operator."*

## 3. Claim-by-claim table

| # | claim | status | locus |
|---|---|---|---|
| 1 | Semi-local Sonin spaces and their ambient Hilbert spaces | **PROVED IN SOURCE** | CCM2024 Def. 4.5, in `L²(X_S)^{K_S}` |
| 2 | Map from the archimedean space; compatibility with Fourier | **PROVED IN SOURCE** | Prop. 4.6(i); Prop. 4.7(i) `F_S θ_S = θ_S F̃_R` |
| 3 | Stability as the finite set of places changes | **PROVED IN SOURCE** | CCM2024-Thm4.6 / Thm 2 — **hilbertian**, see §4 below |
| 4 | Semi-local prolate operator: formal expression, domain, self-adjointness, identification with Sonin space | **PARTLY PROVED / PARTLY PROPOSED.** The semi-local formulation is given and the Sonin space corresponds to the negative spectrum; but the archimedean identification holds only *"up to a finite dimensional possible discrepancy"* (CCM2024 §1), and a **second candidate** operator is deferred to a forthcoming paper. Domain/self-adjointness qualifications are flagged there as *"still at the formal algebraic level"* for the metaplectic description. | CCM2024 §1, §2, §5 |
| 5 | The existing general semi-local trace formula | **PROVED IN SOURCE (not checked here)** — Connes, Selecta Math. 5 (1999) 29–106, as reported by CCM2024 §1. **NOT CHECKED** directly. | CCM2024 ref. [7] |
| 6 | A specifically **positive** Sonin/Gram trace on convolution squares | **PROVED IN SOURCE at the archimedean place only** (`Tr(ϑ(g)Sϑ(g)*) = ‖ϑ(g)S‖²_HS`). **NOT FOUND IN CHECKED SOURCES semi-locally.** | CC2021 Thm 1 |
| 7 | An independent formula connecting that positive trace to the arithmetic Weil functional, with correction terms | **PROVED IN SOURCE at the archimedean place** (CC2021-Thm4.6: `Tr(ϑ(f)S) = W_∞(f) + E(f)`, `E` defined independently by eq. 85). **PROPOSED IN SOURCE semi-locally** ("a more precise strategy", "we expect"). | CC2021-Thm4.6; CCM2024 §1 |
| 8 | A sign or quantitative estimate controlling those corrections | **PROVED IN SOURCE at the archimedean place, on one window** (CC2021 Lemma 6.10 / Thm 6.11, `g` supported in `[2^{-1/2},2^{1/2}]`). **NOT FOUND IN CHECKED SOURCES semi-locally.** | CC2021 Lem. 6.9, 6.10, Thm 6.11 |

**Stages 1–5 exist. Stages 6–8 exist archimedean-only.** The earlier audit collapsed these; this
table separates them.

## 4. "Hilbertian" is not "unitary" — and the consequence

**Footnote 2, CCM2024 p. 5, verbatim:** *"We use the term 'hilbertian' to denote the underlying
topological vector space structure of a Hilbert space."*

So CCM2024-Thm4.6 asserts a **bounded invertible linear isomorphism**, i.e. a homeomorphism of
topological vector spaces — **not** an isometry and **not** a unitary.

| map | property | intertwines | `S`-dependence |
|---|---|---|---|
| `θ_S : S_λ(R,e_∞) → S_λ(X_S,α)` | **bounded, invertible, NOT isometric** | Fourier: `F_S θ_S = θ_S F̃_R` (Prop. 4.7(i)) | multiplier `∏_{p∈S\{∞}} L_p(½+is)^{-1}` on the critical line (eq. 58) |
| `η_S` | the **dual partner**, `η_S ≠ θ_S` and `η_S ≠ θ_S^{-1}` | — | multiplier `∏ L_p(½−is)` |
| `θ_S, η_S` jointly | **`⟨θ_S(f)|η_S(g)⟩ = ⟨f|g⟩`** (Prop. 4.7(iii)), i.e. `θ_S^* η_S = Id`, so `η_S = (θ_S^*)^{-1}` — **not** `θ_S^* = θ_S^{-1}` | — | — |
| `U_S = F_μ w_S : L²(X_S)^{K_S} → L²(R)` | **unitary** (stated in the proof of Prop. 4.7(iii)) | — | — |
| `υ_S : S_λ(X_S,α) → B_λ` | isomorphism **of hilbertian spaces** | diagram (62) | `B_λ` carries **different inner products** per `S` (§4.8) |

**Elementary deduction (checked here, not in the sources).** On the critical line the multiplier has
modulus `∏_p |1 − p^{−1/2−is}|`, which ranges over `[∏(1−p^{−1/2}), ∏(1+p^{−1/2})]` as `s` varies.
So the distortion is bounded for each fixed `S` but grows without limit:

| `S\{∞}` | lower | upper | condition number |
|---|---|---|---|
| `{2}` | 0.292893 | 1.707107 | 5.83 |
| `{2,3}` | 0.123791 | 2.692705 | 21.75 |
| `{2,3,5}` | 0.068430 | 3.896920 | 56.95 |
| `{2,3,5,7}` | 0.042566 | 5.369817 | 126.15 |
| `{2,3,5,7,11,13}` | 0.021486 | 8.927244 | 415.50 |

Since `Σ_p p^{−1/2}` diverges, `∏(1+p^{−1/2}) → ∞` and `∏(1−p^{−1/2}) → 0`. **Fixed-`S` norm
equivalence is not uniform control as `S` grows.** (Recorded, not pursued.)

**Why this blocks the transport of positivity.** For bounded invertible `J` and an orthogonal
projection `P`:

- `J P J^{-1}` is **idempotent but need not be self-adjoint** in the target inner product — so it
  is a projection, but not an *orthogonal* one, and `Tr(A (JPJ^{-1}) A^*)` has no sign.
- `J P J^*` is **positive but need not be idempotent** — so it is not a projection, and it is not
  the Sonin compression.

Neither object is `S_{semi-local}`, and neither supplies an arithmetic trace identity. One may of
course *redefine* the inner product to make `θ_S` unitary — but then the transported metric is a
different one, and the Weil functional in the original normalisation is no longer the thing being
bounded. §4.8 says this in the authors' own words: the inner product depends on `S`.

## 5. The exact remaining gap

Reconstructed archimedean relation, signs fixed (unchanged from the audit):

```
Tr( theta(f) S ) = W_inf(f) + E(f),        E(f) := int f(rho^{-1}) eps(rho) d*rho     [CC2021-Thm4.6, no support restriction]
W_inf(g*g*)      >= Tr( theta(g) S theta(g)* ) - c |ghat(0)|^2,   c = 4 gamma / log 2 [CC2021-Thm6.11]
```

for `g ∈ C_c^∞` supported in `[2^{-1/2},2^{1/2}]`, `ĝ(−i/2)=0` (and `ĝ(0)=0` for Theorem 1);
trace-class/HS hypotheses as in the source. `Pole` and `Σ_p W_p` stay **outside** `W_∞`.

> **MISSING STATEMENT, stated precisely.** For a finite `S ∋ ∞`, a cutoff `λ > 0`, and a test space
> to be specified: an identity
> `Tr_S( ϑ_S(f) S_λ(X_S,α) ) = W_S(f) + E_S(f)` with `E_S` **defined independently** of `W_S`, a
> proof that the left side is an **ordinary** (not merely regularised) positive trace on convolution
> squares `f = g*g*`, and a bound on `E_S` — with explicit test space, support, normalisation,
> quantifiers, and **finite-place dependence** — strong enough to leave `Σ_{p∈S} W_p` controlled.

Nothing in the checked sources supplies this. CCM2024 supplies the **spaces** (stage 1–3) and the
**operator** (stage 4, with qualifications) and calls the comparison a *strategy*.

**Distinctions kept:** ordinary vs regularised trace; self-adjointness vs nonnegativity (a
self-adjoint operator's trace functional is positive *on squares*, which is not the same as the
operator being positive); formal prolate expression vs self-adjoint realisation; compactness vs
sign or small norm; correspondence of spaces (established) vs control of the Weil form (not).

**`S = {∞,2}` as a source-level diagnostic only:** stages 1–5 are available there; stages 6–8 are
not. No numerical first-prime experiment was run.

**Is the missing inequality an intermediate estimate or the target itself?** As stated it is
**intermediate** — it is restricted to one finite `S` and one cutoff `λ`, and would not by itself
give RH, which needs all `S` and all windows with uniform control. The `θ_S` distortion table above
is a concrete reason to expect the uniformity to be the hard part. **Positive-trace comparison is
also not the only conceivable route**, and nothing here says otherwise.

**Pending leads, not followed:** a 2025 item *"New eigenfunctions for the negative part of the
Connes–Moscovici prolate spectrum"* surfaced in search; it concerns eigenfunctions, not the
positive-trace comparison, and is **NOT CHECKED**. The forthcoming paper promised in the CCM2024
abstract (second candidate semi-local prolate operator) had **not** been located as of this search
date. **This is not an exhaustive literature review.**

## 6. What changes for our project

**A. Does the additional source remove the earlier claim that semi-local geometry is unavailable?**
**Yes.** Semi-local Sonin spaces are defined directly, mapped from the archimedean space
compatibly with Fourier, and proved stable under enlarging `S`. The earlier statement is amended.

**B. Does it supply the positive-trace comparison we need?** **No.** The authors call it a
strategy they *expect* to work.

**C. Does it improve the Weil-window result?** **No.** No new bound, and no change to the window
arithmetic: the archimedean comparison still lives at `L = ½log2 ≈ 0.3466`, still inside Zhu's
certified `L ≤ 0.8`, which still contains the active primes `n = 2,3,4`. **Adding a prime to `S` is
not the same operation as enlarging `L`** — the first changes the ambient adelic space, the second
changes the support of the test function. The audited logarithmic dictionary, the Sonin cutoff `λ`
and our half-width `L` all remain distinct, and `g` is still distinguished from `g*g*`.

**D. Does our projected-second-moment calculation contribute a proved bound?** **No — explanatory
context only**, exactly as recorded before. Nothing here forces our finite-matrix perturbation
theorem into this operator family, and similar terminology supplies no map, kernel, gap or
unbounded-form extension.

## 7. Plain language

Our earlier audit said the semi-local version of this machinery did not exist. That was fair about
the 2020 paper we read, which set up the spaces and openly deferred the geometry — but a 2024 paper
by the same authors with Moscovici did the deferred work. The semi-local Sonin spaces are now
defined outright, connected to the archimedean one by an explicit map, and shown to be stable when
you enlarge the set of primes. **So the geometry is there, and we were wrong to say otherwise.**

What is still not there is the part we would actually need. The archimedean argument works because
a certain projection is orthogonal, which makes a trace a sum of squares and therefore positive for
free. The map carrying the archimedean picture to the semi-local one is only an isomorphism of
*topological* vector spaces — the authors' own footnote says "hilbertian" means exactly that — and
it rescales the Mellin transform by the inverse local L-factors. That rescaling is bounded for any
fixed set of primes but gets worse without limit as primes are added. A bounded invertible map
carries an orthogonal projection to something that is still idempotent but no longer orthogonal, or
to something positive that is no longer a projection. Either way the free positivity does not come
along for the ride.

The authors say as much: the 2024 paper describes a "more precise strategy" for the semi-local
positivity and says they "expect" it to open a way. **A strategy is not an estimate**, and we should
not record one as the other.

**Which parts exist, and what exactly is missing?** Existing: the semi-local spaces, the map and its
Fourier compatibility, stability in `S`, the prolate operator (with a finite-dimensional discrepancy
and a second candidate deferred), and the general semi-local trace formula from 1999. Missing: an
independently defined semi-local discrepancy `E_S` in an identity `Tr_S = W_S + E_S`, proof that the
semi-local trace is an ordinary positive trace on convolution squares, and a bound on `E_S` with
explicit quantifiers and explicit dependence on the finite places.

---

## 8. Record

### Files changed (all inside this directory)

| file | before (sha256, 16) | after (sha256, 16) | change |
|---|---|---|---|
| `SEMILOCAL_LITERATURE_UPDATE.md` | *(did not exist)* | *(this file)* | **created** |
| `RESULT.md` | `f37e80934e6ca97a` | `d28392ab5b6287d5` | correction banner after the verdict; §2 bullet 2 amended; two plain-language passages amended |
| `POSITIVITY_DEFECT_AUDIT.md` | `de392701a37db90d` | `2bf820a771b96305` | correction banner in the preamble; §4 item 2 amended (first clause withdrawn, second clause stated precisely, hilbertian caution added); "no semi-local `P_F`" blockquote annotated as unchanged |
| `SOURCES_AND_SCOPE.md` | `4dd4be4347cf0d0a` | `9549b417c583a3ce` | update banner; CCM2024 added to the source table with what was and was not read; three NOT CHECKED / NOT LOCATED rows added; Selecta 1999 identified as CCM2024 ref. [7]; Ramanujan title/DOI citation correction; the 2020-paper paragraph's final sentence struck and replaced with a narrowed claim plus a CCM2024 established/not-established paragraph |
| `MECHANISM_DICTIONARY.md` | `a105fb767fbaf4ff` | `2d429b7f8b9db07f` | note at the head; table row 4 clarified ("open semi-locally" retained, meaning narrowed); one row added to §2 on bounded-invertible images of orthogonal projections |

**Unchanged, verified by hash:** `checks.py` `ff3bfe147e898dc1`, `checks_output.txt`
`063a9c538b7f44dd`. No script was re-run and no numerical result was recomputed.

**Unchanged elsewhere:** `paper/weil_window_control.tex` md5 `a04434a72d5434d5076c1cabcf4e3a5a`,
`paper/weil_window_control.pdf` md5 `a29f9670844f341c7986cfc1d4382bbe` (byte-identical to the
Zenodo deposit `10.5281/zenodo.23212998`). Root `README.md`, all `src/`, `scripts/`, `notes/` and
`data/`, the earlier `extensions/openai-quasi-rh-audit-2026-10-07/`,
`extensions/second-moment-generalization-2026-10-07/` (including both projection notes) and
`extensions/comb-dips/` were **not touched**. The `thales-ramanujan-quarter-gamma` repository was
**not touched**. No persistent-memory file was written.

### Sources read in this pass

| source | what was read |
|---|---|
| **CCM, arXiv:2310.18423v2, 4 May 2024, 30 pp** | abstract; §1 Introduction; §4 incl. §§4.6–4.8 (Def. 4.5, Props. 4.6–4.7, eqs. 57–59, CCM2024-Thm 4.6 p. 23, diagram 62); footnote 2 p. 5; intro Theorem 2 p. 5; reference list. §§2–3, 5–6 **skimmed only** |
| `~/Downloads/QuarterGammaRamanujan (1).pdf` | title page; md5 `517e779cef517dc177f51056117829d5`, byte-identical to the Zenodo PDF already read; three identical copies present |

### Unverified dependencies

- **A. Connes, *Selecta Math. (N.S.)* 5 (1999), 29–106** (the general semi-local trace formula) — **NOT CHECKED**; reported as cited by CCM2024 §1 (its ref. [7]).
- CCM2024 **§§2–3, 5–6** — **NOT CHECKED in detail**; the prolate-operator claims in the table's row 4 rest on §1's own summary of them.
- *"New eigenfunctions for the negative part of the Connes–Moscovici prolate spectrum"* (2025) — **NOT CHECKED**, not retrieved.
- The forthcoming CCM paper (second candidate semi-local prolate operator) — **NOT LOCATED** as of 2026-10-07.
- **This is not an exhaustive literature review.** Absence below is "not found in the sources checked", never "does not exist" — which is precisely the error being corrected here.

### Git state (unchanged by this pass)

```
HEAD                 13ad9dd  "PROJECTION_CONTEXT: add the interpretive principle (magnitude, alignment, response)"
branch               main
origin/main...main   0 ahead, 0 behind
git status --short   ?? extensions/intrinsic-positivity-bridge-2026-10-07/
```

Nothing was staged, committed, pushed, deposited or memorised. The whole extension directory
remains **untracked**, exactly as before this pass.
