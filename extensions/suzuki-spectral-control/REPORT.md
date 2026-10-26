# REPORT.md — Suzuki spectral-function control

Lab record. Extension of `weil-window-control`; **no DOI, not part of the paper, not
peer-reviewed.** Preregistration, committed before any computation:
[`PREREG.md`](PREREG.md) (commit `e8b6a9d`). Derivations:
[`notes/`](notes/). Transcripts: [`results/`](results/).

Source: M. Suzuki, *Spectral convergence associated with Weil's quadratic form*,
[arXiv:2610.32704v1](https://arxiv.org/abs/2610.32704), §§1–3 and §6 read in full before
any code was written.

---

## 0. The result, and the non-claims

**Headline: a certified positivity window is an unconditional spectral result.** The RH
hypothesis in the local part of Suzuki's programme supplies `λ_a > 0` and nothing else;
every step of the chain uses only positivity `(2.4)` and von Neumann self-adjointness.
A checkable certificate therefore substitutes for RH at that window, and since `λ*` is
non-increasing in the window width, one certificate at `L₀` covers every `a ≤ L₀`.
Precise statement, hypotheses and exclusions: [`TRANSFER_RULE.md`](TRANSFER_RULE.md).
**Verified instance: `X = 3`.** It is one instance.

## Non-claims

Fixed in PREREG before any number existed.

1. **Nothing here proves or advances RH.** Theorems 1.1 and 1.2 assume it, irreducibly.
   *Amended by Phase 2, not quietly:* the clause as preregistered also said "and Prop 6.1",
   and that part is now **superseded** — Prop 6.1's conclusions are shown to hold without
   RH at `a = ½log 3`, where `λ_a > 0` is certified rather than inferred (§4). **The
   non-claim itself is unchanged**: an unconditional statement at one fixed window is not a
   step toward RH, which needs `λ_a > 0` for *every* `a` and is *equivalent* to it by Weil's
   criterion. Positivity at finitely many windows was never the obstacle.
2. Roots of `W_N` at a **fixed** window are local spectral values of a finite-dimensional
   Riesz approximation — **not** zeros of `W`, **not** zeros of `ζ`.
3. Suzuki's `a → ∞` limit stays **RH-conditional** whatever is found here.
4. Planted spectra are **fabricated diagnostics**.
5. Liu's certificates are **UNREFEREED**; anything resting on them says so every time.

### Evidence tiers

**T1** exact/certified (exact over `ℚ`, or mpmath with precision-doubling agreement at a
stated tolerance) · **T2** numerical · **T3** interpretation. `float64` appears nowhere.
Non-detections certify nothing. Upper and lower bounds are never interchanged.

---

## 1. Verdicts

| phase | verdict | one line |
|---|---|---|
| 0 (GATE) | **PASS** | `a = L`, `λ_a = λ*(L)`, constant 1 — confirmed to 49 digits against the repository's own code |
| 1 | **PASS** | Suzuki's Tables 1 and 2, the forced-zero lattice, and the `39.3044` artifact all reproduced |
| 2 | **PASS** (scoped to `X = 3`) | the three cited inputs located and each verified **RH-free**; the unconditional list is **Prop 6.1 and Lemma 3.2 only** — positivity `(2.4)` and von Neumann self-adjointness, nothing via de Branges |
| 3 | **PASS** | matched-window comparison: `λ_min(G) = 1.7419e−17` at `a = 0.8047` against repo Gate B `≤ 2.2702e−17` and Zhu `≥ 8.9e−18` at `L = 0.8` |
| 4 | **MIXED / NULL** | prediction holds at `R = 3,4` to within the bisection resolution (`≤ 1.5e−35`, a bound — Retraction 30), fails at `R = 5`; rate is square-root, not `δ²`; interpretation unsupported |

---

## 2. Phase 0 — normalization reconciliation (GATE) · **PASS**

Full derivation: [`notes/normalization.md`](notes/normalization.md). Transcript:
`results/phase0.txt`.

**Q0-a: `a = L`. Q0-b: `λ_a = λ*(L)`, constant `1`.** The Fourier conventions are
identical (`û(z) = ∫u e^{izx}dx` on both sides); the spectral sides are identical (sum over
distinct ordinates with multiplicity weights = sum over zeros with multiplicity); the
supports are identical; and Suzuki's `λ_a := inf σ(A_a)` is literally Zhu's
`inf Q(f)/‖f‖²`.

| test | result | tier |
|---|---|---|
| **P0-T1** algebraic identification | pinned exactly, constant `1` | T1 |
| **P0-T2** cross-basis vs `src/zeta_window.build_parts_idx` | worst relative deviation **`6.8e−49`** at `L = 0.5493`, **`8.5e−49`** at `L = 0.8047`, full `\|j\|,\|k\| ≤ 8` block (threshold was `1e−25`) | T1 |
| **P0-T3** spectral side (assumes RH) | deviation sign `= (−1)^{m+n}` in 5/5 pairs; `\|dev\|` ratios `1.551, 1.589` vs predicted `T`-ratios `1.699, 1.733` — the `1/T` truncation tail, to 8.7% | T2 |
| **P0-T4** structural | `G` real to `0.0`, symmetric to `0.0`, index-reversal symmetric to `0.0`, positive definite (`λ_min = 1.04e−7` at `a = 0.5493`, `2.00e−13` at `a = 0.8047`) | T1 |

**P0-T2 is the decisive one.** Suzuki's exponential-basis Gram, assembled here from an
independent closed form, predicts the repository's cosine-basis window form through
`cos(kπx/L) = (√(2L)/2)(χ_k + χ_{−k})`, `f_0 = √(2L)χ_0`. The repository's code was not
written in this extension and was not modified.

**One derivation worth recording.** Suzuki truncates the archimedean series and bounds the
tail at `10^{−188}`. That is **impossible by truncation** here: the summand bracket is
`O(1/k²)`, so `10^{−188}` would need `k ~ 10^{188}` terms. The closed form derived in
`notes/normalization.md` §4 replaces it with **one analytic, non-oscillatory 1-D integral
per entry** — no series tail to defend. Yoshida (5.15)/(5.16), which Suzuki cites for the
same object, was **not obtained**; the assembly here is independent and reconciled against
the repository instead.

---

## 3. Phase 1 — reproduce Suzuki §6 · **PASS**

Transcript: `results/phase1.txt`. Data: `data/phase1.json`. Gram matrices cached in
`data/gram_X*_N18_dps*.json`, built once at `N = 18` per `(X, dps)` and sliced for
`N = 14, 16` (the index sets are sub-blocks at the same `a`).

<!-- BEGIN:PHASE1_TABLES -->
### P1-T5 — the `J` convention in (6.4)

Of Suzuki's six Table-1 values at `X = 5`, the typeset (6.4) — **overline present** — reproduces **6/6**; the overline-dropped misreading reproduces **2/6**, agreeing on the first two roots and then diverging. Adopted: *typeset (6.4), overline present*. See Retraction 25 and preregistration defect D2.

### P1-T1 — Table 1, first six positive roots of `W_18`

| `j` | `X=3` | `X=5` | `X=7` | `X=11` |
|---|---|---|---|---|
| 1 | `14.1347274566` | `14.1347251417` | `14.1347251417` | `14.1347251417` |
| 2 | `21.0226378897` | `21.0220396388` | `21.0220396388` | `21.0220396388` |
| 3 | `25.0192189912` | `25.0108575802` | `25.0108575801` | `25.0108575801` |
| 4 | `30.6272491287` | `30.4248761361` | `30.4248761259` | `30.4248761259` |
| 5 | `33.6051319793` | `32.9350616846` | `32.9350615877` | `32.9350615877` |
| 6 | `38.6158405572` | `37.5861810892` | `37.5861781588` | `37.5861781588` |

**Agreement with every one of Suzuki's ten printed decimal places, at every `X`** (`P1-T1: PASS`).

### P1-T2 — Table 2

| `X` | `d_6^{(14)}` | `d_6^{(16)}` | `d_6^{(18)}` | `κ₂(G_{a,18})` | Suzuki |
|---|---|---|---|---|---|
| 3 | `1.04` | `1.03` | `1.03` | `4.70302e+7` | `1.04, 1.03, 1.03`; `κ₂ = 4.70e7` |
| 5 | `5.21e-6` | `3.0e-6` | `2.93e-6` | `1.84684e+17` | `5.21e-6, 3.00e-6, 2.93e-6`; `κ₂ = 1.85e17` |
| 7 | `1.11e-9` | `1.52e-11` | `5.43e-13` | `8.24353e+25` | `1.11e-9, 1.52e-11, 5.43e-13`; `κ₂ = 8.24e25` |
| 11 | `1.28e-9` | `3.7e-15` | `1.03e-17` | `2.32604e+34` | `1.28e-9, 3.70e-15, 1.03e-17`; `κ₂ = 2.33e34` |

### P1-T3 — the forced-zero lattice

`φ_n(kπ/a) = √(2a)·sincpi(k+n) = 0` whenever `k+n ≠ 0`, so for `|k| > N` **every** `φ_n` vanishes and `W_N(kπ/a) = 0` identically.

| `X` | `N` | `|W_N|` at `k = N+1, N+2, N+5` | min inside / precision floor | preregistered criterion | old `1e−8` relative-to-peak cut |
|---|---|---|---|---|---|
| 3 | 14 | `0.0, 0.0, 0.0` | `3.035e+110` | **pass** | pass |
| 3 | 16 | `0.0, 0.0, 9.13e-137` | `2.52e+110` | **pass** | pass |
| 3 | 18 | `0.0, 0.0, 0.0` | `1.469e+110` | **pass** | pass |
| 5 | 14 | `0.0, 0.0, 0.0` | `2.563e+106` | **pass** | FAIL |
| 5 | 16 | `8.65e-132, 0.0, 2.64e-132` | `6.609e+105` | **pass** | FAIL |
| 5 | 18 | `0.0, 0.0, 0.0` | `3.444e+105` | **pass** | FAIL |
| 7 | 14 | `0.0, 0.0, 0.0` | `9.0e+104` | **pass** | FAIL |
| 7 | 16 | `0.0, 0.0, 4.42e-125` | `4.339e+103` | **pass** | FAIL |
| 7 | 18 | `0.0, 6.81e-124, 0.0` | `7.056e+102` | **pass** | FAIL |
| 11 | 14 | `0.0, 0.0, 0.0` | `1.493e+104` | **pass** | FAIL |
| 11 | 16 | `0.0, 0.0, 6.96e-117` | `4.909e+102` | **pass** | FAIL |
| 11 | 18 | `0.0, 0.0, 0.0` | `1.935e+101` | **pass** | FAIL |

The dichotomy is **not marginal at any row**: outside the lattice `|W_N|` is at the precision floor (and exactly `0.0` at most sampled points), while inside it exceeds that floor by the margin in column four — of order `10^100` throughout. All 12 rows pass the preregistered criterion: **True**.

> **The `1e−8` column is shown because it is the one that was originally implemented, and it passes only 3/12 rows.** That cut is relative to the *peak* of `|W_N|` on `[0,40]`, a span that grows with `X`, so it rejects values like `1.06e+8` that are nonzero by a hundred orders of magnitude. It tested the dynamic range of `W_N`, not whether `W_N` vanishes. See preregistration defect **D4**.

### P1-T4 — Suzuki's extra zero of `W_14` at `30π/log 11`

`30π/log 11 = 39.304376916222` (Suzuki prints `39.3044`) **is** `15π/a`, i.e. lattice index `k = 15 > N = 14`; `|W_14|` there, relative to the local scale, is `0.0`.

**Preregistered corollary, confirmed:** for `N = 18` the first lattice point is `k = 19`, at `38π/log 11 = 49.7855440939`, which is **outside** Suzuki's `[0,40]` search range (`True`) — so no such artifact can appear in his `N = 18` table, and none does.

<!-- END:PHASE1_TABLES -->
---

## 4. Phase 2 — unconditional fixed-window statement · **PASS** (scoped to `X = 3`)

Dependency analysis: [`notes/prop61_dependency.md`](notes/prop61_dependency.md).
The three cited inputs, run to ground:
[`notes/prop61_dependency_chain.md`](notes/prop61_dependency_chain.md).
Transcript: `results/phase23.txt`.

**The reduction.** Every step of Prop 6.1's proof that can be inspected in
arXiv:2610.32704v1 factors through coercivity `(2.4)` alone. RH's role is located
precisely: §2.2 takes `λ = 0` *"under RH"*, against the RH-free construction with `λ` below
`inf σ(A_a)` — so **RH is used to supply `λ_a > 0`**, exactly what a certificate replaces.
Prop 6.1 does **not** cite Theorem 3.3, the one irreducibly RH-dependent result (it builds
the global `H_W`, needing `Q_W ≥ 0` on all of `C_c^∞(ℝ)`); Thm 3.3 serves **Theorem 1.1**,
the `a → ∞` limit.

**The three previously unverified inputs are now verified RH-free.** Both sources were
obtained and read **from page renders at 200 dpi, not from text extraction** — a rule this
extension adopted the hard way, since Retraction 25 was caused by `pdftotext` and `pypdf`
silently dropping an overline in Suzuki's (6.4). The `[2]`/`[10]` reference mapping was
checked against the printed bibliography of 2610.32704v1 rather than assumed.

| # | input | located at | proof uses | verdict |
|---|---|---|---|---|
| 1 | `[2, Lemma 2.2]` — `span{e^{iπnx/a}}` is a form core | **Connes–Consani**, arXiv:2106.01715 p. 11 (Lemma 2.2), completed p. 13 | closed form for `Û^n = O(1/\|s\|)`; `Q_∞(f,f)=∫\|f̂\|²·2∂_tθ/2π` with `∂_tθ ≍ log\|t\|` a **smooth archimedean weight, not a sum over zeros**; elementary `sin²` bounds; rapid decay of Fourier coefficients | **RH-free** |
| 2 | `[10, §3.2, proof of Thm 1]` — `C_c^∞(−a,a)` is a form core | **Suzuki**, arXiv:2606.09096 pp. 14–15 (Lemma 3.1 + inclusion (3.2)) | an explicit cutoff `η_ε`, Plancherel, one integration by parts, a dyadic split at `\|z\|=1/ε`, giving `‖e_n − v_ε‖²_Q ≪ ε log(1/ε) → 0` | **RH-free** |
| 3 | `[10, §3]` — the form-domain characterization | **Suzuki**, arXiv:2606.09096 pp. 12–14 | (2.8), assembled from the **geometric side only** — `Re ψ(¼+iz/2)`, `−log π`, the **finite** prime sum over `n ≤ e^{2a}`, and `r''_{0,a}` — plus Schmüdgen's standard form theory | **RH-free** |

Two checks worth stating because they could have gone the other way:

- **Inclusion vs. equality.** Suzuki states only `𝔇(Q_W^a) ⊂ H^log(−a,a)`, while 2610.32704
  asserts a *set-theoretic equality*. The equality is nonetheless supported: (2.8) is an
  **identity with a bounded error term**, not a one-sided estimate, so it reads in both
  directions. (Suzuki separately notes that the *Gagliardo*-type space of Bombieri is
  "slightly larger" than `𝔇(Q_W^a)` — a different space, not the one being cited, and the
  two are not conflated here.)
- **The numbering difference is not a discrepancy.** 2610.32704 cites the exponential core
  as `[2, Lemma 2.2]`; 2606.09096 cites it as `[3, Proposition 2.3]`. Both are in
  Connes–Consani §2.1 and Prop 2.3 is proved *from* Lemma 2.2. Same RH-free fact.

Corroborating but **not relied on**: arXiv:2606.09096 §1.2 opens *"we emphasize that none of
them depends on RH"*, and across the whole of its §§2–3 the string `RH` occurs **exactly
once** — flagging RH as *equivalent to* `Γ ⊂ ℝ`, i.e. as a hypothesis not being assumed.

> **PASS (scoped).** At `a = ½log 3 = 0.549306…`, where `λ_a > 0` comes from Zhu's certified
> `λ*(0.8) ≥ 8.9e−18` together with monotonicity of `λ*` in the window width
> (`0.5493 < 0.8`, so `λ*(0.5493) ≥ λ*(0.8) > 0`), the statements below hold **without
> assuming RH**.

**The list is closed, and the closure rule is: positivity `(2.4)` and von Neumann
self-adjointness only.** Every entry was checked against the typeset proof; anything routed
through the de Branges structure is excluded and listed separately underneath.

| unconditional at `a = 0.5493` | what carries it |
|---|---|
| **`ℋ(A_a)` is a Hilbert space**, equal to the form domain `𝔇(A_a^{1/2})` with the `Q_W`-norm | `(2.4)` makes `Q_W` positive definite on `C_c^∞`, so the completion is a Hilbert space; the identification needs the `C_c^∞` core (input 2) and the form-domain characterization (input 3) |
| **`𝒟_a u = iu'` is symmetric with deficiency indices `(1,1)`**, and von Neumann gives all self-adjoint extensions `D̄_{a,θ}` (3.2) | §3.1, verified from the p. 8 render: `Ker(𝒟_a^* − z) = ℂu_{z,a}` (3.1) follows from `(2.5)`, uniqueness of the Riesz representative and **density of `C_c^∞(−a,a)` in `ℋ(A_a)`** (input 2); equal deficiency norms `‖u_{i,a}‖² = J_a(i,i) = J_a(−i,−i) = ‖u_{−i,a}‖²` come from reflection. Then von Neumann [7, Thm X.2] |
| **Lemma 3.1** `𝔇(D̄_{a,θ})` reflection-invariant iff `θ ∈ {0,π}`; `0` an eigenvalue of `D̄_{a,0}` but not `D̄_{a,π}` | reflection `ℛw_{a,θ} = e^{iθ}w_{a,−θ}`; and `𝔇(D̄_a) ⊂ H¹₀` because *"applying `(2.4)` to both the functions and their derivatives gives convergence in `H¹` and zero boundary values"* — `(2.4)` only |
| **`W(a,θ;·)` is entire of exponential type `≤ a`** `(3.4)` | `u_{±i,a} ∈ L²(−a,a) ⊂ L¹(−a,a)`, so their transforms are entire with the stated bound |
| **Prop 6.1 (i)** `u^{(N)}_{z,a} → u_{z,a}` in `ℋ(A_a)`, uniformly on compacts | the two form cores (inputs 1–2) + `(2.4)` make `⋃_N 𝒳_{a,N}` dense; Riesz uniqueness gives `u^{(N)} = P_{a,N}u`, `P_{a,N} → I` strongly; `(2.4)`+`(2.5)` give `‖u_{z,a} − u_{w,a}‖ ≤ λ_a^{−1/2}‖e_z − e_w‖_{L²}`, which makes the convergence uniform on compacts |
| **Prop 6.1 (ii)** `J_{a,N} → J_a` uniformly on compacts of `ℂ²` | `J_a − J_{a,N} = ⟨(I−P)u_{w̄,a}, (I−P)u_{z̄,a}⟩_{A_a}` + Cauchy–Schwarz. `J_a` is used **only** through its definition `(2.6)` as that inner product — never through the de Branges reproducing-kernel formula `(2.7)` |
| **Prop 6.1 (iii)** `W_N(a,π;·) → W(a,π;·)` locally uniformly | `(3.3)`, `(6.4)` and **Lemma 3.2** |
| **Prop 6.1 (iv)** `W_N(a,π;0) ≠ 0` for large `N`, and the normalized limit `(6.5)` | same, with `W(a,π;0) ≠ 0` from Lemma 3.2 |
| **Lemma 3.2** `W(a,π;·)` even and entire; `W(a,π;0) ≠ 0`; **all zeros simple**, coinciding with the eigenvalues of `D̄_{a,π}` | evenness from reflection symmetry of `J_a` + `(3.3)`; `W(a,π;0) ≠ 0` from Lemma 3.1 (`0` is not an eigenvalue of `D̄_{a,π}`); **simplicity from the von Neumann decomposition `(3.2)`** `u_{t,a} = u_0 + c_t w_{a,π}` with `c_t ≠ 0`, giving `∂_z W\|_{z=t} = ‖u_{t,a}‖²_{A_a}/c_t ≠ 0` — and `‖u_{t,a}‖²_{A_a} > 0` **is** positivity |
| **the zeros of `W(a,θ;·)` are real** (§3.2, immediately above Lemma 3.2) | `W(a,θ;z) = 0` ⟹ `u_{z̄,a} ∈ 𝔇(D̄_{a,θ})`; Prop 2.1 gives `u_{z̄,a} ≠ 0` — the **no-common-zeros** clause, proved by choosing a test function whose transform does not vanish at `z`, which is RH-free; so `z̄` is an eigenvalue of a **self-adjoint** operator, hence `z ∈ ℝ` |
| **`G_{a,N}` itself** | Suzuki states it outright on p. 14: *"The computation of `G_{a,N}` does not require RH or any information about the zeros of `ξ(s)`."* Confirmed independently by Phase 0 |

**Excluded from the list — everything routed through the de Branges structure.** Not because
it is false, but because it is not established by positivity and self-adjointness:

| **not** claimed unconditional | why |
|---|---|
| **`Ĥ(A_a)` is a de Branges space** (Prop 2.1 as a whole) | its **third** axiom — norm preservation under the nonreal-zero transformation — needs `Γ ⊂ ℝ`; see qualification 1 below |
| the generating function `E_a`, and formula `(2.7)` | these are where the third axiom earns its keep; Prop 6.1 never invokes them |
| **Theorem 3.3**, and **Theorems 1.1 / 1.2** | the global space `H_W` needs `Q_W ≥ 0` on all of `C_c^∞(ℝ)`, i.e. RH |
| `(5.3)`, `(5.4)`, **Corollaries 5.1 / 5.2** | built on `E_0`, `A/E` and Theorem 1.2 |

Two uses of Prop 2.1 **are** on the unconditional list, and they are the two that do not
touch the third axiom: the **no-common-zeros** clause (above), and the ***description* of
`ℋ(A_a)`** cited in §6 — *"These functions belong to `ℋ(A_a)` by its description in the
proof of Proposition 2.1; no membership in `𝔇(D_a)` is required"* (p. 14). That description
is the `H^log` characterization, i.e. **input 3**, already verified RH-free; the
`χ_{n,a}` qualify because `φ_{n,a}(z) = √(2a)·sin(az+nπ)/(az+nπ) = O(1/\|z\|)` makes the
log-weighted integral converge. Prop 2.1's first two axioms are likewise fine — point
evaluation is continuous by the `(2.4)` bound, and conjugation `F ↦ F^♯` preserves the norm
because `Γ = −Γ` and `Γ = conj Γ` hold **unconditionally** with matching multiplicities. It
is the third axiom alone that fails to reduce.

**What the certificate actually buys.** `λ_a > 0` on *some* window is **already a theorem
without RH** — arXiv:2606.09096 Theorem 1.4 gives
`λ_a = log(1/a) + μ₁ − log(2π) + ψ(2) − 1 + O(a)` for "sufficiently small `a`", as does
Yoshida [17, Lemma 2]. But that is **ineffective in three places at once**: the threshold is
unquantified, `μ₁` is unnamed, and the `O(a)` has no stated constant, so it cannot be
invoked at `a = 0.549`. That the asymptotic is far out of range there is visible
numerically — its explicit part is `log(1/0.5493) − log(2π) + ψ(2) − 1 = −1.8158`, an `O(1)`
quantity, against a measured `λ_min(G_{0.5493,6}) ≈ 1.04e−7` (the `N = 6` entry of the ladder in §5; at `N = 18` the same quantity is `6.93e−8` — see the reconciliation there). **The whole content of the
upgrade is that a certified bound supplies an *effective* instance of `(2.4)` at a *named*
window where the unconditional theorem is silent.** It does not prove positivity for all
`a`, which is equivalent to RH by Weil's criterion.

**Three qualifications that survive, and are not erased by the upgrade.**

1. **Prop 2.1's third de Branges axiom still needs `Γ ⊂ ℝ`** — now confirmed from the
   typeset page 5 of 2610.32704. Membership `F_w ∈ Ĥ(A_a)` is unconditional (it goes
   through the `H^log` characterization, an integral over the **real** axis, where the
   Blaschke factor has modulus 1). But the *norm* equality `‖F_w‖ = ‖F‖` is justified by
   "the definition of `Q_W`", whose exact value is the zero-sum
   `Q_W(v₁,v₂) = Σ_{γ∈Γ} m_γ v̂₁(γ) conj(v̂₂(γ̄))` — and arXiv:2606.09096 p. 13 says in terms
   that `γ ∈ Γ` *may be non-real*. **So this step uses RH itself, not `(2.4)`.** It does not
   bite — Prop 6.1 uses `J_a` only through its definition (2.6) and uses Prop 2.1 only for
   the *no-common-zeros* clause, which is proved RH-free — but its existence shows the
   reduction to `(2.4)` is **not automatic**, which is the failure mode PREREG named in
   advance.
2. **The PASS is at `X = 3` only.** `a(X=5) = 0.8047 > 0.8` is outside Zhu's range and stays
   conditional on Liu's **unrefereed** certificate, as does `X = 7`; `X = 11` stays
   numerical only.
3. **Nothing here is about `a → ∞`.** Theorems 1.1/1.2 rest on Thm 3.3 and the global
   `H_W`, which needs RH; arXiv:2606.09096 itself describes that global construction as made
   *"under RH"*. No fixed-window certificate relaxes them, and the roots of `W_N` at fixed
   `a` are spectral values of a finite-dimensional Riesz approximation — **not** zeros of
   `W`, and **not** zeros of `ζ`.

<!-- BEGIN:PHASE2_TABLES -->
### P2-T5 — prime-free positive control

At `a = 0.3` (`2a < log 2`) the von Mangoldt list is **empty** (`0` terms): `G` is archimedean + pole with **no arithmetic input at all**, and is positive definite, `λ_min = 0.0079864309`, `λ_max = 2.4118891`.

### P2-T2/T3/T4 — roots by window, with their certification status

| `X` | `a` | status | precision-doubling shift | min distance to the lattice |
|---|---|---|---|---|
| 3 | `0.54930614` | **CERTIFIED (Zhu, L <= 0.8)** | `3.99e-121` | `0.7101` |
| 5 | `0.80471896` | **CONDITIONAL (Liu, UNREFEREED)** | `4.12e-121` | `0.8068` |
| 7 | `0.97295507` | **CONDITIONAL (Liu, UNREFEREED)** | `3.31e-121` | `0.6459` |
| 11 | `1.1989476` | **NUMERICAL ONLY** | `4.35e-121` | `0.05971` |

Every root is stable under precision doubling (140 → 180 digits) and none lies on the forced lattice. The label carried with each window is:

* `X = 3` — CERTIFIED roots of the finite approximation W_N
* `X = 5` — roots of W_N, CONDITIONAL on Liu (UNREFEREED)
* `X = 7` — roots of W_N, CONDITIONAL on Liu (UNREFEREED)
* `X = 11` — roots of W_N, NUMERICAL ONLY

<!-- END:PHASE2_TABLES -->
---

## 5. Phase 3 — matched-window comparison · **PASS**

### Reconciliation — the two `λ_min(G)` values quoted at `X = 3`

<!-- BEGIN:RECONCILE -->
Phase 0 reports `λ_min(G) = 1.042312874e−7` at `a = 0.5493061443`; Phase 1, 2 and 3 report `6.9320137e−8` at the **same** `a`. **Neither is a slip and they are not in conflict** — they are the same quantity on two different basis sizes, and the smaller one is *required*, not merely permitted.

`G_{a,N}` is the Gram on the `L²`-orthonormal `{χ_{n,a} : |n| ≤ N}`, a `(2N+1)`-dimensional subspace of the admissible class. Enlarging the subspace can only **lower** a minimized Rayleigh quotient, so `λ_min(G_{a,N})` is non-increasing in `N` and bounded below by `λ*(a)` for every `N`. Phase 0 quotes `N = 6` (a `13×13` matrix); Phase 1/2/3 quote `N = 18` (`37×37`).

| `N` | dim | `λ_min(G_{a,N})` | ratio to previous | |
|---|---|---|---|---|
| 2 | 5 | `1.23566369317e-5` | `—` |  |
| 4 | 9 | `3.05582191509e-7` | `0.0247302` |  |
| 6 | 13 | `1.04231287448e-7` | `0.341091` | **← Phase 0 quotes this** |
| 8 | 17 | `1.02560145378e-7` | `0.983967` |  |
| 10 | 21 | `9.4900447736e-8` | `0.925315` |  |
| 12 | 25 | `8.46432475844e-8` | `0.891916` |  |
| 14 | 29 | `7.74780095048e-8` | `0.915348` |  |
| 16 | 33 | `7.31248572668e-8` | `0.943814` |  |
| 18 | 37 | `6.93201372131e-8` | `0.94797` | **← Phase 1/2/3 quote this** |

Strictly decreasing at **every** step (`monotone = True`), as the variational principle requires; `N = 6 → 18` is a factor `0.6651`. Same `a`, same dps (140), same assembly, one script: `src/reconcile_lmin.py`, transcript `results/reconcile_lmin.txt`.

> This is also the finite-`N` shadow of **Connes–Consani Cor. 2.4** — the sequence decreases *to* `λ*(a)` — which is the RH-free warrant for Phase 3's "variational upper bound" language. **Both numbers are upper bounds on `λ*(0.5493)`; neither is a lower bound and neither certifies positivity.** Positivity at this window comes from Zhu.

<!-- END:RECONCILE -->


Transcript: `results/phase23.txt`.

**The `ln κ₂` agreement is demoted to a sanity check and is not the Phase 3 result.**
See Retraction 28. `κ₂ = λ_max/λ_min` by definition, so
`ln κ₂ / (−ln λ_min) = 1 + ln(λ_max)/(−ln λ_min)`, and with `λ_max ≈ 3` and `−ln λ_min`
between 16 and 78 the ratio is **forced to within 7 % of 1 by arithmetic alone**, for any
matrix with a bounded top eigenvalue and a tiny bottom one. It tests that `λ_max = O(1)`.
It carries no information about Zhu's growth law.

**The Phase 3 result is the matched-window comparison.** Three independently computed
numbers, at two windows 0.59 % apart:

| source | value | window | **direction** |
|---|---|---|---|
| this extension, `λ_min(G_{a,18})`, closed-form Gram, Suzuki's exponential basis, dps 140 | `1.7419e−17` | `a = 0.8047189562` | variational **UPPER** bound |
| repo Gate B, `N = 24`, `even_d` basis (second basis `2.3357e−17`) | `≤ 2.2702e−17` | `L = 0.8` | variational **UPPER** bound |
| Zhu, certified | `≥ 8.9e−18` | `L = 0.8` | certified **LOWER** bound |

**Why `λ_min(G)` is an upper bound.** `G` is the Weil form's Gram on an `L²`-orthonormal
basis, so `λ_min(G_{a,N})` is a Rayleigh quotient minimized over an `N`-dimensional
**subspace** of the admissible class; restricting to a subspace **raises** the minimum.
Hence `λ_min(G_{a,N}) ≥ λ*(a)` for every `N`, decreasing to it. That the limit is exactly
`λ*(a)` is **Connes–Consani Corollary 2.4** (arXiv:2106.01715 p. 13, read from the page
render), which is itself RH-free — a by-product of chasing down Phase 2's input 1.

**The window offset is real and is not smoothed over.** `a = 0.804718956…` and `L = 0.8`
differ by `4.719e−3`. These are **different windows**, so this is a consistency comparison,
not a reproduction. `λ*` is non-increasing in the window width, so
`λ*(0.804719) ≤ λ*(0.8)`, and the observed ordering `1.74e−17 < 2.27e−17` (ratio `0.767`) is
in the direction monotonicity predicts — a weak sign check, not a derivation, since both are
upper bounds at different `N` and their ordering is not forced.

> **What this does NOT establish — the Retraction-18 failure mode, named in PREREG.**
> Zhu's **lower** bound is at `L = 0.8`, and `a = 0.8047 > 0.8`. Monotonicity runs the
> **wrong way** to carry a lower bound outward: `λ*(0.8047) ≤ λ*(0.8)`, so `8.9e−18` is
> **not** a lower bound at `a = 0.8047`. **Nothing here certifies `λ_a > 0` at `X = 5`.**
> The Phase 2 PASS is at `X = 3` (`a = 0.5493 < 0.8`), where the direction does work.

**What it does establish.** A third independent implementation — different basis, different
quadrature, a different codebase from the repository's Gate B — lands within a factor `1.30`
of the repo's upper bound at a window 0.6 % away, and inside Zhu's certified interval
`[8.9e−18, 2.27e−17]`. That is a **cross-implementation consistency check on the assembly**,
at tier **T2**.

<!-- BEGIN:PHASE3_TABLES -->
### P3-T2 — sanity check only (see Retraction 28)

| `X` | `a` | `λ_min(G_{a,18})` | `λ_max` | `ln κ₂` | `−ln λ_min` | ratio |
|---|---|---|---|---|---|---|
| 3 | `0.54930614` | `6.9320137e-8` | `3.26014` | `17.6663` | `16.4845` | `1.07169` |
| 5 | `0.80471896` | `1.7418742e-17` | `3.21696` | `39.7574` | `38.589` | `1.03028` |
| 7 | `0.97295507` | `4.1255062e-26` | `3.40087` | `59.6741` | `58.45` | `1.02094` |
| 11 | `1.1989476` | `1.3221089e-34` | `3.07528` | `79.1321` | `78.0087` | `1.0144` |

The ratio sits within **7.17 %** of 1 — but it was always going to. `κ₂ = λ_max/λ_min` identically, so the ratio is `1 + ln(λ_max)/(−ln λ_min)`, and with `λ_max ≈ 3` and `−ln λ_min` between 16 and 78 it is **forced** into `[1.014, 1.072]` by arithmetic. The preregistered factor-3 threshold could not have failed. What the column actually measures is that `λ_max = O(1)` — retained on that basis, and **not** reported as a test of Zhu's growth law. See Retraction 28.

### P3-T3 — the Phase 3 result: matched-window comparison

| source | value | window | direction |
|---|---|---|---|
| this extension, `λ_min(G_{a,18})`, closed-form Gram, dps 140 | `1.7418742e-17` | `a = 0.8047189562` | variational **UPPER** |
| repo Gate B, `N = 24`, `even_d` basis (2nd basis `2.3357e-17`) | `≤ 2.2702e-17` | `L = 0.8` | variational **UPPER** |
| Zhu, certified | `≥ 8.9e-18` | `L = 0.8` | certified **LOWER** |

The windows differ by `0.004719` (0.59 %), so this is a consistency comparison, **not** a reproduction. `λ*` is non-increasing in the window width, so `λ*(0.8047) ≤ λ*(0.8)`, and the observed ordering (ratio `0.7673`) runs the way monotonicity predicts — a weak sign check, since both are upper bounds at different `N`.

> **Zhu's lower bound does not transfer.** It is certified at `L = 0.8`, and `a = 0.8047 > 0.8`; monotonicity runs the **wrong way** to carry a lower bound outward (`λ*(0.8047) ≤ λ*(0.8)`). **Nothing here certifies `λ_a > 0` at `X = 5`.** The Phase 2 PASS is at `X = 3`, where `a = 0.5493 < 0.8` and the direction does work. This is the Retraction-18 failure mode, named in PREREG.

That `λ_min(G_{a,N})` decreases to exactly `λ*(a)` is **Connes–Consani Cor. 2.4** (arXiv:2106.01715 p. 13), itself RH-free — found while chasing Phase 2's input 1.

*gateB_repo, gateB_repo_basis2 and zhu_lower are CITED literature values kept at their published number of significant figures; every other number in this file is computed here and serialized at full working precision.*

<!-- END:PHASE3_TABLES -->
---

## 6. Phase 4 — function-field and planted-spectrum control · **MIXED / NULL**

Construction fixed before computing: [`notes/ff_construction.md`](notes/ff_construction.md).
Prior art recorded before Phase 4: [`notes/prior_2607.24830.md`](notes/prior_2607.24830.md).
Transcript: `results/phase4.txt`.

### 4a — the analogue exists; the preregistered KILL does not fire

`φ_i(w) = w^i` on the unit circle generates the Hallouin–Perret Toeplitz Gram from the
Frobenius angles exactly as `φ_{n,a}(γ)` generates `G_{a,N}` from the ordinates. The
singularity of `T_R` for `R > d` is a **range restriction** (the construction runs on
`R < d`), not an obstruction.

Where the analogue genuinely stops: `W_N` comes from the deficiency structure of `i d/dx`
— indices `(1,1)`, deficiency points `±i`, a `θ`-circle of self-adjoint extensions — and a
finite Toeplitz matrix has none of that. What replaces it is classical and already contains
the Gram: `T_R` is the moment matrix of the Frobenius measure on the circle, `J` is the
Christoffel–Darboux kernel **verbatim**, and the `θ`-family is the **para-orthogonal
polynomials** `B_R(z;β) = Φ_R(z) − βΦ*_R(z)`. So `W_N ↔ B_R` and *zeros real* ↔ *zeros on
the unit circle*.

### 4b — genuine spectra
<!-- BEGIN:PHASE4B -->
Blocks: three **online** pairs at angles `0.7, 1.9, 2.8` — six atoms, `d = 6`.

| `R` | `λ_min(T_R)` | `max‖z|−1|` | max distance to a Frobenius angle |
|---|---|---|---|
| 2 | `4.7300472` | `0.0` | `0.1615` |
| 3 | `4.5509431` | `0.0` | `0.4255` |
| 4 | `4.0279414` | `7.779e-62` | `0.3533` |
| 5 | `1.4283655` | `0.0` | `0.3416` |
| 6 | — (`T_6` singular) | `0.0` | `7.779e-62` |

**All zeros lie on the unit circle to the precision floor at every `R`**, as the classical theorem requires for a positive definite moment matrix. At `R < 2g` they do *not* reproduce the Frobenius angles — they are the **Gauss-quadrature nodes** of the Frobenius measure, which only interlace its support. At **`R = 2g = 6`** — the number of atoms, and the largest `R` for which `B_R`'s matrix `T_5` is nonsingular — the Gauss rule becomes **exact** and the zeros land on the six angles **to `10^{−61}`** — `max distance = 7.779e−62`, against `0.34` at `R < 2g`, sixty orders of magnitude. **Not "exactly", and not `0.000e+00`**: that figure was a `float64` cast rounding `7.78e−62` down (Retraction 29), and `7.78e−62` is the resolution of the computation, not a property of the spectrum.

> So the analogue of Suzuki's `W_N` **reproduces the Frobenius spectrum to `10^{−61}`**, at the one `R` where it can, on the side where RH is Weil's theorem.

<!-- END:PHASE4B -->
### 4c — planted off-circle spectra
<!-- BEGIN:PHASE4C -->
Blocks: two **online** pairs (`0.7, 1.9`) + one **quartet** at angle `2.8` with `ρ = e^δ` — eight atoms, `d = 8`.

| `R` | `δ*` from `λ_min(T_R)` | `δ*` from `λ_min(T_{R−1})` | `δ*` where zeros leave the circle | rel. gap to `T_{R−1}` | departure exponent `p` |
|---|---|---|---|---|---|
| 3 | `0.6822246993101394` | `0.8776671466177126` | `0.8776671466177126` | `≤ 1.14e-35` | `0.7146` |
| 4 | `0.6616174734011838` | `0.6822246993101394` | `0.6822246993101394` | `≤ 1.47e-35` | `0.5688` |
| 5 | `0.5618039833415618` | `0.6616174734011838` | `1.610837505968302` | `≤ 1.43` | `0.706` |

**The preregistered prediction holds at `R = 3, 4` to within the bisection resolution** (`≤ 1.14e-35` at `R = 3`, `≤ 1.47e-35` at `R = 4`) **and fails at `R = 5`.** The two bisections return bit-identical values there, but they share a bracket and a halving sequence, so that is a **resolution bound, not an identity** (Retraction 30). At `R = 5` the zeros stay on the circle to `δ = 1.6108` while positivity is lost at `0.6616` — a factor `2.4`. So the converse of the classical implication (positive definite ⟹ zeros on the circle) is **false**, which corrects the equivalence asserted in `notes/ff_construction.md` §4 (Retraction 26a).

**The departure rate is square-root, not `δ²`.** Fitted `max‖z|−1| ~ (δ−δ*)^p` gives `p = 0.7146, 0.5688, 0.706` across `R` — consistent with `1/2` (a generic square-root bifurcation when two circle zeros collide and leave), with the residual bias explained by the finite offsets used. There is no `δ²` content.

> **Verdict on 4c.** Per PREREG's outcomes this is neither "coincide within 10%" (it fails at `R = 5`) nor a clean offset law (it is *exact* at `R = 3, 4`). And the preregistered *reading* — "the `δ²`–`L³` law governs a second pipeline" — is **not supported either way**: where the thresholds coincide, the coincidence is the classical Szegő implication, not the `δ²`–`L³` law.

<!-- END:PHASE4C -->
---

## 7. Retraction log

Continues the **main** series of the parent `REPORT.md`, which reaches 20. (The
Davenport–Heilbronn extension keeps a separate `R1–R6` series.)

**Retraction 21 — the pole term's conjugation.** `src/suzuki.py` first wrote the pole term
as `φ_m(i/2)conj(φ_n(i/2)) + c.c.` The correct term is `2 Re[φ_m(i/2)φ_n(i/2)]`: the
continuation of `û(t)conj(f̂(t))` is `û(z)conj(f̂(z̄))`, so
`h(i/2) = φ_m(i/2)conj(φ_n(−i/2)) = φ_m(i/2)φ_n(i/2)`. Both versions are real and they
**coincide for a real basis**, so the repository's real-basis code cannot distinguish them
— only the derivation does. Cross-check: the correct form reproduces `zeta_window`'s
`pole_sign = −1` for the odd sector automatically; the wrong one does not. *Found before any
Phase 1 result; no reported number was affected.*

**Retraction 22 — a factor 2 in the prime term.** First written as
`2Λ(n)n^{−1/2}S(log n)`; the correct coefficient is `Λ(n)n^{−1/2}S(log n)` with
`S = C(v)+C(−v)`. The repository's `2Λ/√n · Csym` equals this because `Csym = S/2` for a
real basis. *Found before any Phase 1 result.*

**Retraction 23 — a factor `2L` in the gate's own test code.** P0-T2 used `1/√(2L)` for the
`j = 0` basis coefficient instead of `√(2L)`: `χ_0` is the **constant** `1/√(2L)` and
`f_0 = 1`, so `f_0 = √(2L)χ_0`. This produced a **1.2–3.0 % cross-basis disagreement and an
initial Phase 0 FAIL**. The gate caught a bug in the gate. After the fix the agreement is
`6.8e−49`.

**Retraction 24 — catastrophic cancellation in the `d = 0` archimedean integrand.**
Computed as `S(v) − 2`, a difference of two numbers near `2`, which loses `~|log v|`
digits; tanh–sinh quadrature samples `v ~ 10^{−2·dps}`, and at `v = 10^{−55}` with dps 60
the naive form is wrong by a **relative `2.65e−7`**. Replaced by the half-angle form
`−4sin²(θ/2) − (v/a)cos θ`. Verified **not** to change the dps 40/60/90 values (which agree
to 28 digits), but necessary at the 140/180 digits Phase 1 uses.

**Retraction 25 — a discrepancy I reported in the source that does not exist.** PREREG
asserted that the brief's restatement `J(w,z) = φ(z)ᵀ conj(G^{-1}φ(w))` and "the paper's
(6.4)" `J(w,z) = φ(z)ᵀ G^{-1}φ(w)` were different functions. **They are — but the paper's
is the first one.** The typeset (6.4) carries an **overline over the whole group**
`G^{-1}φ(w)`, which **both `pdftotext` and `pypdf` silently drop**. Rendering page 15 as an
image shows it plainly. The brief transcribed the paper correctly; my text extraction did
not. Measured consequence: the conjugated form reproduces **6/6** of Suzuki's Table-1
values at `X = 5`; the overline-dropped misreading reproduces **2/6**, agreeing on the first
two roots and then diverging (`25.0108575815` against `25.0108575802`, and
`37.5862119696` against `37.5861810892`). *Lesson, recorded: do not infer a discrepancy in a
source from a text-extraction artifact — render the page.*

**Retraction 26 — an overstated equivalence of my own, and an off-by-one.**
(a) `notes/ff_construction.md` §4 asserted that positive definiteness of the moment matrix
is **equivalent** to all zeros of `B_R` lying on the circle, and used that to pre-emptively
weaken Phase 4c's interpretation. The classical theorem is **one-directional**
(positive definite ⟹ zeros on the circle); the converse is false, and `R = 5` exhibits it —
the zeros stay on the circle to `δ = 1.611` while positivity is lost at `0.662`.
(b) The first 4c run compared the zero-departure threshold against `λ_min(T_R)`, but `B_R`
has degree `R` and is built from `T_{R−1}`; against the correct matrix the agreement at
`R = 3, 4` is **exact** rather than `28.7 %` and `187 %` off. Both the wrong and the right
comparison are reported in `results/phase4.txt`.

**Retraction 27 — two serialization-precision floors, same signature.**
Both were found by the same tell: *a value repeating identically across rows that should
differ.*
(a) `γ₆` was transcribed at 10 decimal places, putting a `2.567e−11` floor under the
`d₆ = t_{a,N,6} − γ₆` column. (b) After that was fixed, the roots themselves were being
stored in the JSON at **14 significant figures** and then differenced against `γ₆` at full
precision, putting a `3.287e−13` floor under the same column — which is exactly where the
`X=7, N=18` and `X=11, N=16,18` entries had saturated. Neither was caught by reading the
numbers; both were caught by the repetition.
**Fix:** `γ₆` is now taken from `mp.zetazero(6)` at dps 60, and **every** mpmath quantity
written to a results JSON goes through a single `fp()` helper that serializes at the full
working precision (`mp.nstr(x, mp.mp.dps)`), so `make_report.py` cannot reintroduce the
loss. With that, Table 2 reproduces **all twelve** of Suzuki's `d₆` values and **all four**
`κ₂` values.
**Counterpart, in the opposite direction:** *cited* constants — Zhu's `8.9e−18`, the repo's
`2.2702e−17` and `2.3357e−17` — are deliberately **not** passed through `fp()`. Padding a
two-significant-figure published bound to 140 digits printed
`8.900000000000000405…e−18`, binary noise from the precision the `mpf` happened to be built
at, displayed as though Zhu had certified 140 digits. Full precision is for what we
computed; cited numbers keep the digits they came with.

**Retraction 28 — a Phase 3 "test" that was arithmetic, not evidence.**
PREREG framed the agreement of `ln κ₂` with `−ln λ_min` as a check of Zhu's growth law
`−ln λ*(L) ∼ C·N(T*)/ln N(T*)`, with "agreement within a factor 3 → PASS". **It is not a
test of that law, or of anything else.** Since `κ₂ = λ_max/λ_min` identically,
`ln κ₂ /(−ln λ_min) = 1 + ln(λ_max)/(−ln λ_min)`; here `λ_max ∈ [3.08, 3.40]` and
`−ln λ_min ∈ [16.5, 78.0]`, so the ratio lies in `[1.014, 1.072]` **by construction**. A
factor-3 threshold could not have failed. The quantity measures that `λ_max = O(1)`, which
is worth knowing and is retained — **relabelled as a sanity check**, with `λ_max` now
printed in the table so the near-automaticity is visible rather than implied. The Phase 3
result is the matched-window comparison instead, which can and does carry a direction-of-
bound statement that could have come out wrong.

**Retraction 29 — "exactly" was a `float64` cast.**
Phase 4b's headline was that at `R = 2g = 6` the POPUC zeros *"recover the Frobenius
spectrum **exactly**"*, reported as `max distance = 0.000e+00` and tested with
`worst_angle_distance == 0.0`. The angle distances were being computed through
`float(·)` casts — the one place `float64` had survived in this extension, in a quantity
that carried a headline. Recomputed in mpmath the distance is **`7.779e−62`**, not zero;
`float64` had rounded it down and the exact-zero comparison then passed.
**The claim is now stated as a bound, with a number: the zeros agree with the Frobenius
angles to `10^{−61}`** (`max distance = 7.779e−62`). The words *exactly* and *exact
recovery*, and the figure `0.000e+00`, are **withdrawn** and do not appear anywhere in this
extension's output. `0.000e+00` is not a measurement; `< 10^{−61}` is.
**Fix:** the casts are gone and the tolerance matches the on-circle test (`< 1e−40`).
*The substance survives and the test still passes* — `7.78e−62` against `0.34` at `R < 2g`
is a gap of sixty orders of magnitude, and `7.78e−62` is the same floor the `R = 4`
on-circle residual sits at, so it is the resolution of the computation and not a property
of the spectrum.
**Why it is logged rather than quietly corrected:** this is the third defect in this
extension found by a *displayed* number rather than by reasoning (Retraction 27(a), 27(b),
now 29), and the second caused by a formatting layer rather than by the mathematics. The
pattern is the finding.

**Retraction 30 — the same defect as 29, found by looking for it: Phase 4c's thresholds
do not "coincide exactly".**
Checking Retraction 29's pattern elsewhere turned up a second instance immediately.
Phase 4c reported that at `R = 3, 4` the zero-departure threshold and the
`λ_min(T_{R−1})` positivity threshold *"coincide **exactly** (to the bisection floor)"*,
with a printed relative difference of `0.0` — and the two values **are** bit-identical at
120 digits. That is not evidence of equality. Both come from the same `thresh()` helper,
which bisects on the same bracket `[1e−6, 3]` with the same halving sequence and the same
stopping rule `hi − lo < 1e−35`, then returns `(lo+hi)/2`. **Two different predicates that
flip inside the same final subinterval therefore return the identical number by
construction.** All the agreement establishes is that the two thresholds lie within one
bisection interval of each other.
**Fix:** the claim is now a bound with its provenance stated — the thresholds agree to a
**relative gap `≤ 1.14e−35` at `R = 3` and `≤ 1.47e−35` at `R = 4`**, computed as
`bisection resolution / δ*` — and the transcript says in terms that bit-identical returns
are a shared-bisection-path artifact, not a proof of equality. The word *exactly* is
withdrawn here too.
**This does not change any verdict.** `R = 5` still fails by a factor `2.4` (relative gap
`1.43`, which is nowhere near the resolution), the converse of the classical Szegő
implication is still refuted, and Phase 4c is still **MIXED / NULL** with its preregistered
reading unsupported. What changes is that an exact-equality claim becomes a bound.
**The pattern, now four deep:** Retraction 27(a), 27(b), 29 and 30 are all the same error —
**a number that was really a floor, a cast or a resolution, reported as if it were a
measurement**. Three were caught by a value repeating where it should vary. The fourth was
caught only because I went looking for the third one's signature elsewhere, which is the
only reason it is in this list rather than still in the report.

## 8. Preregistration defects

Recorded rather than silently repaired, because a preregistration that is quietly edited is
not a preregistration.

**D1 — P0-T3's threshold was unreachable.** PREREG fixed `1e−8` absolute on the
arithmetic-vs-zero-side comparison. The zero sum truncated at `NZ` zeros has a tail
`~ 2(2a)³log(T/2π)/(2πT)`, so `1e−8` needs `T ~ 10^8`, i.e. `~10^7` zeros. Replaced —
**before** looking at whether the replacement passes — by the test that actually
discriminates a normalization error from a truncation tail: the deviation must carry the
sign `(−1)^{m+n}` and must fall like `1/T`. A fixed tolerance would not have distinguished
a constant offset from a shrinking tail; this does.

**D2 — P1-T5 rested on a false premise.** See Retraction 25. The test was rewritten, before
the data was examined, as a measurement of how far the overline-dropped misreading
diverges.

**D3 — Phase 4c's interpretation was weakened before the data existed, on a claim that
itself turned out to be too strong.** See Retraction 26(a). The amendment in
`notes/ff_construction.md` §4 stands as written and is corrected here.

**D4 — P1-T3's "bounded away from zero" was implemented as a relative-to-peak cut, and
that cut, not the mathematics, failed 9 of 12 rows.**

PREREG's wording, verbatim and unedited:

> **PASS: `|W_N(kπ/a)|` at the precision floor for all `|k| > N` tested, and bounded away
> from zero at sampled `|k| ≤ N`, at every `(a,N)`.**

**Verdict against that wording: PASS at all 12 `(a,N)` rows.** Taking the two clauses in
turn, as written:

| PREREG clause | measured | rows |
|---|---|---|
| *"at the precision floor for all `\|k\| > N`"* | `0.0` exactly at most sampled points, and `≤ 9.13e−117` relative at every other — against a floor of `10^{−115}` | **12/12** |
| *"bounded away from zero at sampled `\|k\| ≤ N`"* | smallest inside value exceeds that same floor by `1.9e+101` to `3.0e+110` | **12/12** |

Both clauses are resolution-relative, and PREREG names no numeric cut for either. **The
`1e−8` that failed nine rows appears nowhere in PREREG.** It was an implementation choice:
`v/scale > 1e−8`, relative to the **peak** of `|W_N|` on `[0,40]` — a quantity that has
nothing to do with being nonzero, and whose dynamic range *grows with `X`*. It therefore
rejected values like `1.06e+8`, nonzero by roughly a hundred orders of magnitude, and did so
for `X = 5, 7, 11` while passing `X = 3`. It was measuring the dynamic range of `W_N`.

**So this is an implementation that did not match its own preregistration, not a
preregistered threshold retuned after seeing data** — the distinction is the whole point,
and it is why **both cuts are now computed and reported per row** in the P1-T3 table above,
rather than the failing one being deleted. A reader can see exactly what changed and judge
it. The `10^{110}` margin is what makes the PASS safe; the per-row disclosure is what makes
it checkable.

**Why it was not caught earlier:** the nine `FAIL` lines were visible in every run of
Phase 1 and I read past them, because the verdict table already said Phase 1 **PASS** on
the strength of Tables 1 and 2. The generated report block was simultaneously printing
*"All 12 rows pass: **False**"* directly beneath that **PASS**. A contradiction printed in
my own output is the cheapest possible signal, and I still needed it put in front of me.

## 9. Out of scope this round

Rigorous `W − W_N` error bounds via Prop 6.1's residual identity. Certified coercivity
supplies only the **lower-bound** ingredient and the constants are numerically vacuous at
these windows: Liu's `2^{−151}` gives `λ^{−1/2} ≈ 2^{75.5}`, so the Cauchy–Schwarz step in
Prop 6.1's proof carries no usable information here. Recorded as future work.

## 10. Repository effect

**None outside this directory.** No file in `paper/`, no core module in `src/`, and no part
of the v1.0 note was modified. The repository's `src/zeta_window.py`, `src/planted.py` and
`src/weilform.py` are **imported, never edited**. One line is added to the root README's
extension list, in a separate commit so it can be dropped independently.
