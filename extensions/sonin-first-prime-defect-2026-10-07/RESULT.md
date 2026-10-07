# RESULT.md

**One-prime Sonin discrepancy test, `S = {∞, 2}`.** Dated 2026-10-07.

## Verdict: **SOURCE/DOMAIN GAP** (retained)

> **AMENDED 2026-10-07 after review.** The verdict stands, but two claims are narrowed: the
> "arithmetic reach" obstruction is **scoped to the bare metric multiplier** and the ceiling
> `L < log 2` for the whole construction is **retracted**, because `Π_S` contains `K^{-1}` and the
> inverse has infinite Laurent support; and the even/odd parity reading is **withdrawn**. See the
> next section and `FIRST_PRIME_DEFECT.md` §§3.2, 3.4.

The identity the test requires — `Tr(ϑ_S(f)Π_S) = W_S(f) + E_S(f)` with `E_S` independently
**evaluable** — could not be justified, and the semi-local trace interpretation needed in
characteristic zero is unavailable. **No positivity conclusion about the semi-local Weil form
follows.** Scoped to the construction tested: the `θ_S`-transported Sonin projection at `S = {∞,2}`,
and Connes' `R_Λ = P̂_Λ P_Λ` cutoff. Nothing here says projection methods in general are impossible.

## The single strongest mathematical statement established

> **Metric and projection reconstructed; arithmetic comparison not completed; no general
> prime-power obstruction proved.**

Concretely, for `S = {∞,2}` and any `λ > 0`, in the common unitary picture (`U_∞`, `U_S` both
unitary, `θ_S` = multiplication by `L_2(½+is)^{-1}`):

```
Tr( vartheta_S(f) Pi_S ) = Tr_{M_inf}( K^(-1) · P ( M F ) P ) ,
   Pi_S = J P K^(-1) P J^*  the ACTUAL orthogonal projection onto S_lambda(X_S, alpha),
   M = |L_2(1/2+is)|^(-2) = 3/2 - 2^(-1/2)( vartheta(2) + vartheta(2)^(-1) ),
   K = (3/2) P - sqrt(2) · P C P,   C = ( vartheta(2) + vartheta(2)^(-1) ) / 2,   ||C|| = 1,
```

conditional on a trace-class hypothesis (T) that no checked source supplies semi-locally; and
`K ≥ (3/2−√2)Id_{M_∞} > 0`, so `K^{-1}` exists with an **unconditionally convergent** Neumann
expansion of ratio `2√2/3 = 0.9428090`, needing only `‖PCP‖ ≤ 1`.

**What survives as a proved reach statement, and its exact scope.** `θ_S^*θ_S` is multiplication by
`∏_{p∈S∖{∞}}|L_p(½+is)|^{-2}`, a Laurent polynomial of degree `≤ 1` in each `p^{is}`; so `HF`
generates only dilates `d/d′` with `d,d′` coprime **squarefree**. Verified by exact rational
expansion for `{2}, {2,3}, {2,3,5}`.

> **This is a statement about the bare multiplier only.** An earlier version of this file inferred
> from it a structural ceiling `L < log 2` for the whole construction, and that the route cannot
> reach Zhu's window because `n = 4` is active there. **That inference is retracted.** `Π_S`
> contains `K^{-1}`, and the inverse has infinite Laurent support:
> `1/M(s) = 2 + 2√2 cos(s log2) + 2 cos(2 s log2) + ⋯` — the `k = 2` ("4") harmonic appears with
> coefficient exactly `2`, verified to `9.0e-14` (`amend_checks.py` §A), and compositions `(PCP)^n`
> reach two steps at `n = 2` in a finite control (§C). So **no general prime-power obstruction is
> proved**, and equally `Π_S` is **not** shown to produce the correct Weil coefficients:
> cancellations in the full assembly are unresolved in both directions.

The parity reading is likewise withdrawn. `∂_s log M` is odd, but
`∂_σ log M_σ(s)|_{σ=1/2} = 2Σ_{k≥1}(log2)2^{-k/2}cos(ks log2)` is **even** in `s` with exactly the
`W_2` weights `0.4901291, 0.3465736, 0.2450645, …` (verified to `3.3e-10`, `amend_checks.py` §B).
This is the expansion of the local Euler factor, not a new arithmetic theorem and not a family of
geometric operators off `σ = ½`. **The obstruction, correctly stated:** the missing ingredient is
not the presence of prime-power information — it is there in both parities — but an **independently
justified operation that extracts it inside the required trace formula, with the correct
coefficients and a controlled remainder.**

## FINAL QUESTION

> For `S = {∞,2}`, does using the **actual** orthogonal projection and the relevant directions
> produce a mathematical estimate the arithmetic Weil problem did not already have?

**No.** Using the actual projection changes the problem's shape in three checkable ways and in none
of them produces a new arithmetic estimate:

1. **It fixes the metric correctly.** `Π_S = J P K^{-1} P J^*` with `K = (PJ^*JP)|_{M_∞}` is the
   genuine orthogonal projection onto `S_λ(X_S,α)` — verified exactly over `Q` to be self-adjoint,
   idempotent, to fix `θ_S(M_∞)` pointwise and annihilate its complement. It is distinct from
   `JPJ^{-1}` (idempotent, not self-adjoint) and from `JPJ^*` (self-adjoint, not idempotent). The
   coercivity it needs holds here: `K ≥ (3/2−√2)Id`. **This is a standard identity, not a result.**
2. **It keeps positivity.** `Tr(ϑ_S(g)Π_Sϑ_S(g)^*) = ‖ϑ_S(g)Π_S‖²_{HS} ≥ 0`, free, with no
   `ϑ_S`-invariance — subject only to a trace-class hypothesis (T) that no checked source supplies
   semi-locally. **So non-unitarity of `θ_S` was never the obstruction**, contrary to the reading in
   the earlier bridge audit, which this pass corrects: `JBJ^* ≥ 0` for `B ≥ 0`, and
   `Tr(JPJ^{-1}) = Tr P` holds for any bounded invertible `J`.
3. **But the arithmetic comparison is not completed.** The reduction's only closed form is the
   mean-field case `PCP = 0`, where it lands on archimedean `W_∞` alone; with `K^{-1}` restored the
   terms are interleaved products no checked source evaluates. *(Amended 2026-10-07: that is an
   unevaluated expression, **not** a proved obstruction — see the retraction above.)* And the one
   route that *does* carry the semi-local arithmetic in characteristic zero — Connes 1999 — carries
   it for the wrong operator in the wrong corner.

## The two routes, and exactly where each stops

**Route B (CCM2024 `θ_S`).** New exact reduction, conditional on (T):

```
Tr( vartheta_S(f) Pi_S ) = Tr_{M_inf}( K^(-1) · P ( M F ) P ),
   M = |L_2(1/2+is)|^(-2) = 3/2 - 2^(-1/2)( vartheta(2) + vartheta(2)^(-1) ),
   K = (3/2) P - sqrt(2) · P C P,   C = ( vartheta(2) + vartheta(2)^(-1) ) / 2.
```

Two by-products worth keeping. (a) `K^{-1}` admits an **unconditionally convergent** Neumann
expansion around the mean-field value `(3/2)P`, ratio `2√2/3 = 0.9428090 < 1`, needing only the
universal `‖PCP‖ ≤ 1`; the convergence is owed specifically to `p^{-1/2} < 1` and **fails exactly at
`a = 1`**, where the multiplier acquires a real zero. (b) The whole quantitative question collapses
to one quantity: `K ≥ (3/2 − √2‖PCP‖)P`, with worst case `‖K^{-1}‖ ≤ 11.65685` against mean-field
`2/3` — a **gap factor of `17.48528`**, and a two-sided ambient comparison `Mx/m = 33.97056`.
Whether `‖PCP‖ < 1` — whether the Sonin space's Mellin image can concentrate near
`s ∈ (2π/log2)Z`, spacing `9.064720` — is **not established in this audit**. Note `‖PCP‖ ≤ 1`
already suffices for invertibility, so a strict inequality would only sharpen the constant.
**Route B stops at D3** — the arithmetic comparison — not at a reach obstruction.

**Route A (Connes 1999, the newly read dependency).** C99-Thm VII.4 **is** the semi-local arithmetic
identity and **does** hold for `k = Q`, `S = {∞,2}`:
`Tr(R_ΛU(h)) = 2h(1)log′Λ + W_∞(h) + W_2(h) + o(1)`. Assembling it with Halmos' two-projection
decomposition gives an identity for the **orthogonal** projection `Q_Λ` with a discrepancy
`E_Λ(h) = Tr(G_Λ U(h))`, `G_Λ = P̂_ΛP_Λ − Q_Λ`, defined by projection geometry alone — satisfying the
independence requirement, with `Tr G_Λ = Σ_j c_j² ≥ 0` and `‖G_Λ‖_1 = Σ_j c_j` over the nontrivial
principal angles. It stops for four separate reasons, each sufficient:

- `E_Λ` is **not evaluable**: the angle spectrum of `(P_Λ, P̂_Λ)` on `L²(X_{\{∞,2\}})` is computed
  nowhere in the checked sources.
- The leading term **diverges**: `h(1) = ∫|g|²d^*x > 0`, so `2h(1)log′Λ → +∞` and the inequality is
  vacuous unless `E_Λ` absorbs it. That is the difficulty itself.
- **Wrong corner.** C99's `B_Λ` (`f` and `\hat f` vanish for `|x| > Λ`) is the doubly **low**-pass
  `ran ∩ ran` corner; the Sonin space (`|x| < λ`) is the doubly **high**-pass `ker ∩ ker` corner.
  Different subspaces of the same pair; no checked source translates between them.
- **Characteristic zero has no `Q_Λ` formula.** C99-Cor VIII.2 is stated "in the case of positive
  characteristic", and C99-Lem VIII.1 proves `[P̂_Λ,P_Λ] = 0` by a function-field count
  (`Λ = q^N`, `mod(C_S) = q^Z`, `|d| = q^{2−2g}`, genus `g`, tally `(2N+1)l − fl + (2−2g)l`, and
  "since `ξ` is locally constant, its Fourier transform has compact support"). **No checked source
  proves `[P̂_Λ,P_Λ] = 0` for `k = Q`.** The global formula (16) is explicitly unproved — "are not
  able to prove (16) directly for arbitrary `h`" — and in positive characteristic C99-Thm VIII.5
  makes it *equivalent to RH*, so it is not an input.

## What is honestly open rather than closed

- **The arithmetic comparison itself (D3).** No semi-local `Tr = W_S + E_S` with evaluable `E_S`.
  This is the whole gap, and it is where both routes stop.
- **Whether `Π_S` supplies the Weil prime-power coefficients.** Unresolved in **both** directions:
  the support calculation does not exclude it (the inverse carries all harmonics), and nothing
  establishes it (the interleaved products `Tr(PD^{k_1}PD^{k_2}P⋯FP)` are unevaluated, and
  cancellations in the assembly are not ruled out).
- **`‖PCP‖ < 1` quantitatively** — *not established in this audit*, rather than "genuinely open":
  its precise formulation, literature status and possible approximate invariant vectors were not
  checked, and lack of exact invariance does not by itself give a uniform norm gap. Even settled,
  it would sharpen `P_S` against `P_∞` by at most `17.5×` and give **no** Weil bound, since D3 is
  blocked independently.
- **Hypothesis (T):** trace-class-ness of `ϑ_S(f)Π_S`. Archimedean-only in the sources.
- On the calibration window `(log2)/2 < L < (log3)/2` the prime-2 block has **only** the `k = 1`
  atoms, weight `(log2)2^{-1/2} = 0.4901291` — the shifts `M` itself supplies. A weight-and-sign
  comparison there would need the exact `h(1)`-coefficient of `W_∞` in CC2021's principal-value
  normalisation, which was **not** pinned down; **no numerical mismatch is asserted.**

## Novelty, stated plainly

| item | status |
|---|---|
| positive Sonin trace; CC2021-Thm 4.6; Lemma 6.9/6.10; Thm 6.11 | **prior art** (Connes–Consani) |
| semi-local Sonin spaces, `θ_S`, its multiplier, hilbertian isomorphism, stability in `S` | **prior art** (CCM2024) |
| the semi-local trace formula `Tr(R_ΛU(h)) = 2h(1)log′Λ + Σ_v W_v(h) + o(1)` | **prior art** (C99-Thm VII.4) |
| `Π = JP(PJ^*JP)^{-1}PJ^*`; Halmos five-part decomposition | **standard operator identities**, verified here, not discoveries |
| the exact reduction `Tr(ϑ_S(f)Π_S) = Tr(K^{-1}PMFP)` with `K = (3/2)P − √2 PCP` | **derived here** |
| unconditional Neumann convergence, ratio `2√2/3`, and its failure at `a = 1` | **derived here** |
| Proposition 3: squarefree reach **of the bare multiplier `θ_S^*θ_S`** | **derived here**; the inference to a ceiling `L < log 2` for the whole construction is **RETRACTED** (the inverse `K^{-1}` has infinite Laurent support) |
| `∂_s log M` and `∂_σ log M_σ` both carry the full `W_2` weight sequence, in opposite parities | **derived here**; elementary Euler-factor expansion. The earlier "even vs odd" parity *barrier* is **WITHDRAWN** |
| the char-0 / char-`p` split in C99 and the complementary cutoff corner | **read off the source**; not stated in the prior audits |
| `E_Λ = Tr(G_ΛU(h))` as an independently *defined* discrepancy | **assembled here**; not evaluable |

No prior-art claim is made for anything marked prior art, and no improvement over CC2021-Thm 6.11 is
claimed: none was found. A positive result on an already-certified window would in any case be a new
method or a sharper estimate, **not** a newly certified RH range — and no positive result was
obtained here.

## Effect on the repository

**None.** No manuscript claim changes. The `δ²–L³` sensitivity law, Theorem 7 as corrected by
Retraction 18, Theorem 8's conditional status and `K(L,γ)`'s status as a norm are all untouched.
The one correction this pass makes is to an earlier *audit*, not to the paper: the
intrinsic-positivity-bridge audit's reading that non-isometry of `θ_S` blocks positivity transport
was **too strong** — `JBJ^* ≥ 0` and `Tr(JPJ^{-1}) = Tr P`, and the actual orthogonal projection
keeps positivity for free. The real block is the **uncompleted arithmetic comparison**, not positivity
transport — *and not a proved reach obstruction either.* **That correction is recorded here only;
no existing audit file was edited in this pass.**

## Checks actually run

`checks.py` → `checks_output.txt` (119 lines); `reach.py` → `reach_output.txt` (63 lines);
`amend_checks.py` → `amend_output.txt` (the post-review amendment checks). The first two scripts
and their outputs are left **exactly as run**; the corrections to their *interpretation* are
recorded in the Markdown, not by rewriting the record.
Python standard library only — `mpmath` is not installable on this box, so the repository's `src/`
regression scripts were **not** run and nothing above depends on them. No scans, installs or
background jobs.

| § | check | required control | result |
|---|---|---|---|
| 1 | `Π = JPK^{-1}PJ^*` self-adjoint, idempotent, fixes `J(M)`, kills `J(M)^⊥`, `trΠ = dim M` | **non-unitary, non-commuting `J`**, exact over `Q` | PASS |
| 1 | distinguish `JPJ^{-1}` (idempotent, not s.a.) and `JPJ^*` (s.a., not idempotent, `tr = 15`) | control 2 of the brief | PASS |
| 1 | `Tr(JPJ^{-1}) = Tr P = 2` | similarity-invariance safeguard | PASS |
| 1 | singular `J` ⟹ `det K = 0`, formula inapplicable | control 5, outside the hypotheses | PASS |
| 2 | `Tr(P̂P) − Tr Q = Σc_j²` on 4 random pairs (dim 6, 8) | commuting pair ⟹ generic part empty, `Tr(P̂P) = dim(ran∩ran)` | PASS |
| 2 | Jacobi eigensolver on a rank-one PSD matrix | the trace-preservation bug that once returned `−8.16` | PASS (`min ev = −8.0e−17`) |
| 3 | `|m_2(s)|² = 3/2 − √2cos(s log2)`, dev `6.7e-16`; `m = 0.08578644`, `Mx = 2.914214`, mean `3/2` | omission and wrong-normalisation controls change the numbers | PASS |
| 4 | the exact three-term operator form, dev `0`; threshold `L = (log2)/2 = 0.3465735` | — | PASS |
| 5 | active prime powers vs `L`; `{2}` on the calibration window | **boundary control**: at `L = (log2)/2` the block is EMPTY and the test vacuous | PASS |
| 6 | `∂_s log M = 2Σ(log2)2^{-k/2}sin(ks log2)`, dev `5.4e-9`; `M` even, generator odd | 3-mode even fit to the series fails, `L^∞` error `0.9359592` | PASS |
| 7 | targeted vs ambient on the same objects | at `‖PCP‖ = 1` the targeted bound **degenerates** to the ambient one | PASS |
| 8 | exact Laurent expansion for `{2}, {2,3}, {2,3,5}`: all exponents in `{−1,0,1}`, no non-squarefree dilate | — | PASS (the "`n = 4` unreachable" reading is **retracted**) |
| A | `1/M(s) = (1−r²)^{-1}(1+2Σ r^k cos kθ)`, `r = 2^{-1/2}`; coefficient at `k = 2` is `2 ≠ 0` | dev `9.0e-14` over 400 nodes | PASS |
| B | `∂_σ log M_σ|_{σ=1/2} = 2Σ(log2)2^{-k/2}cos(ks log2)`, **even** in `s` | dev `3.3e-10`; parity verified | PASS |
| C | `(PCP)` has `(0,2) = 0` but `(PCP)²` has `(0,2) = 0.25` in a finite shift model | `C = Id` generates no shift | PASS |
| D | `K` invertible from `‖B‖ ≤ 1` alone, `2√2/3 = 0.9428090 < 1` | — | PASS |
| 9 | Neumann ratio `2√2/3 = 0.9428090 < 1` | **diverges exactly at `a = 1`** (real zero) | PASS |

Finite-dimensional positivity does **not** certify continuum positivity and is not used to. No
quadrature, basis-truncation, tail or inverse-error budget for `B_λ` was established, so §2's samples
are not evidence about `‖PCP‖`; no sign claim rests on a numerical agreement; nothing here is an
enclosure.

## Files created (this directory only)

`SOURCES_AND_DEPENDENCIES.md`, `TARGET.md`, `FIRST_PRIME_DEFECT.md`, `RESULT.md`, `checks.py`,
`checks_output.txt`, `reach.py`, `reach_output.txt`, `amend_checks.py`, `amend_output.txt`.
No manuscript draft.
