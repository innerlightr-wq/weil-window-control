# POSITIVITY_DEFECT_AUDIT.md

All citations are to the **July 4, 2021** PDF (see `SOURCES_AND_SCOPE.md`). Signs are fixed
explicitly throughout. **Everything in §§1–3 is Connes–Consani's; §4 is the arithmetic step and it
is blocked; §5 is what, if anything, this repository adds.**

> **CORRECTION 2026-10-07 — [`SEMILOCAL_LITERATURE_UPDATE.md`](SEMILOCAL_LITERATURE_UPDATE.md).**
> §4 item 2 below asserted that there is **no** semi-local geometry. That was **too strong**:
> semi-local Sonin spaces and the relevant structural identifications are established in later work
> (Connes–Consani–Moscovici, arXiv:2310.18423v2, 2024, Def. 4.5, Props. 4.6–4.7, Thm 4.6 = Thm 2).
> This audit did not locate the specific positive-trace comparison with sufficient discrepancy
> control needed for our application. **This is a gap in the verified transfer, not an absence of
> semi-local geometry.** The amended text is inline at §4 item 2. §§1–3 are unaffected.

## 1. The positive construction, and why it needs no invariance

`ϑ` is the unitary scaling representation of `R*_+` on `L²(R)_ev`, `ϑ(λ)ξ(v) = λ^{-1/2}ξ(λ^{-1}v)`;
`S = Π_{S(1,1)}` is the orthogonal projection onto the Sonin space of even functions vanishing on
`[−1,1]` together with their Fourier transform. For `f ∈ C_c^∞(R*_+)`, `ϑ(f) = ∫f(λ)ϑ(λ)d*λ`.

For `f = g*g*` one has `ϑ(f) = ϑ(g)ϑ(g)*`, hence

```
      Tr( theta(f) S )  =  Tr( theta(g) S theta(g)* )  =  || theta(g) S ||_HS^2   >=  0 ,
```

the middle step by `S = S² = S*` and cyclicity **under the trace-class hypothesis**, the last by
Hilbert–Schmidt. **This is a Gram square, not the compression of a positive operator to an
invariant subspace** — which is exactly why `ran S` need not be `ϑ`-invariant, and indeed it is
not: the paper says so immediately after Theorem 1 ("the scaling action `ϑ` does not restrict to
the Sonin's subspace, it can be *compressed*"). `checks.py` §2 is the finite-dimensional control:
`A S A*` is PSD with all leading minors `≥ 0` for an explicit `A` that does **not** commute with
`S`, and it is PSD precisely because it equals `(AS)(AS)*`.

**This is a third kind of projection**, distinct from both already documented in
`extensions/second-moment-generalization-2026-10-07/PROJECTION_CONTEXT.md`: it is neither an exact
reducing projection (which commutes with the operator) nor a baseline-kernel spectral projection
(which fixes a leading perturbative coefficient). It is a **Gram/compression** projection, positive
by construction and carrying no invariance.

*Caveat recorded, not waved through:* cyclicity needs its hypotheses, and a regularised trace is
not automatically a positive trace. The paper's admissibility and trace-class conditions were
**read but not re-derived here**; they are taken as given.

## 2. The exact identity and the signed discrepancy

**Theorem 4.6 (eq. 84), no support restriction:**

```
      Tr( theta(f) S )  =  W_inf(f)  +  E(f),        E(f) := ∫ f(rho^{-1}) eps(rho) d*rho ,
```

with `ε` given independently by eq. (85) and `ε(ρ^{-1}) = ε(ρ)`. Solving for the Weil side, with
the sign fixed:

```
      W_inf(f)  =  Tr( theta(f) S )  -  E(f) .
```

**`E` is defined independently of `W_∞`, not by subtraction** — this is what makes the identity
useful rather than vacuous. Its sign is *not* assumed favourable: Theorem 6.11 exists precisely
because `E` need not vanish.

**The discrepancy is then estimated, not merely named.** In the proof of Theorem 6.11: put
`k = Y*g` and `ξ(x) = k(e^x)`, so `ξ ∈ H = L²(I)`, `I = [−½log2, ½log2]`; then
`E(f) = ⟨ξ | N_I ξ⟩` with `N_I = −2ℓ₁(1+ε₁)(Id − K_I)` (Proposition 5.5). **Lemma 6.10** gives the
operator inequality

```
      < xi | N_I xi >   <=   gamma | < xi_0 | xi > |^2 ,      gamma ~ 2.94355 ,
```

i.e. `N_I` is dominated *above* by a rank-one positive operator along the normalised constant `ξ₀`.
With `⟨ξ₀|ξ⟩ = (log2)^{-1/2} k̂(0)` and `k̂(0) = −2ĝ(0)` this yields **Theorem 6.11**:

```
      W_inf(g*g*)  >=  Tr( theta(g) S theta(g)* )  -  c |ghat(0)|^2 ,    c = 4 gamma / log 2 ,
```

for `g` supported in `[2^{-1/2},2^{1/2}]` with `ĝ(−i/2) = 0`. Adding `ĝ(0) = 0` gives **Theorem 1**.
(`c = 16.98658…` is derived from the quoted `γ`; the source prints no decimal for `c`, and
`checks.py` §4 verifies only the arithmetic, not an invented value.)

**What compactness does and does not do.** The compactness results (including those of
arXiv:2008.10974) establish *structure* — triangularity modulo compacts, infinite-dimensionality of
the Sonin spaces. **Compact does not mean small in norm, sign-definite, or uniformly controlled as
the window grows**, and no step here uses it that way.

## 3. The alignment theorem, reconstructed and tested

**Lemma 6.9.** For unit `φ, ψ` and `a,b,c ≥ 0`, the form
`B(ξ) = −b|⟨φ|ξ⟩|² + a|⟨ψ|ξ⟩|² + c‖P_φ ξ‖²` is positive **iff**

```
      a + c >= b        and        b(a+c)  <=  a(b+c) |<phi|psi>|^2 .
```

**Two-dimensional reduction, verified.** Since `c ≥ 0`, positivity is decided on
`E = span{φ, ψ}`; in the orthonormal basis `(φ, ψ₂/‖ψ₂‖)` the matrix is
`[[aα²−b, aαβ],[aαβ, aβ²+c]]` with `α = ⟨φ|ψ⟩`, `β = ‖ψ₂‖`, `α²+β² = 1`; trace `a+c−b`,
determinant `aα²(b+c) − b(a+c)`. `checks.py` §3 confirms in exact rational arithmetic that the
stated criterion agrees with "trace `≥ 0` **and** determinant `≥ 0`" on five cases, **including
deliberately failing ones**: with `(a,b,c) = (3,1,1)` fixed, alignment `|⟨φ|ψ⟩|² = 1` passes and
`|⟨φ|ψ⟩|² = 1/10` fails. **Magnitude alone does not settle it; alignment does.**

**How the authors use it.** They apply Lemma 6.9 with `φ = ζ`, `ψ = ξ₀` to bound `N_I` by a
rank-one along `ξ₀`, combining an eigenvalue bound on `K_I` (`λ_max = 1.05158`, `ε₁ ≈ 0.00122`,
Remark 6.12 — quoted, not re-derived here) with the vector overlap `⟨ξ₀|ξ⟩`. The negative
direction is controlled because the available positive rank-one is **aligned** with it and the
complement carries a gap.

**Does our perspective supply a sharper bound for the actual discrepancy operator? No.** Their
defect `N_I` is already reduced to a rank-one majorant, which is the best shape the rank-one lemma
can deliver; sharpening would require a better eigenvalue bound on `K_I` or a better overlap
estimate, both of which live entirely inside their Toeplitz analysis. Our machinery
(`Theorem A` of the second-moment exploration) is a *perturbative expansion in a parameter* for a
different operator family (`2cosh` kernel) answering a different question. **There is no map, and
none is asserted.** Re-deriving a Schur complement or renaming their quantities would not be a new
estimate, and is not offered as one.

## 4. The arithmetic step — attempted, and precisely blocked

**At the Connes–Consani window the decomposition closes, and it is entirely theirs.** Using the
dictionary of `MECHANISM_DICTIONARY.md` §3, for `g ∈ C_c^∞` supported in the open
`(2^{-1/2}, 2^{1/2})` with `ĝ(±i/2) = 0` and `ĝ(0) = 0`:

```
      Q_zeta(g)  =  P_F(g)  +  R_F(g),
      P_F(g)  :=  Pole(g)  +  Tr( theta(g) S theta(g)* )        >= 0   (both terms, independently)
      R_F(g)  :=  - E(g*g*)  -  sum_p W_p(g*g*) .
```

At that window `Σ_p W_p(g*g*) = 0` exactly (no active prime power — `checks.py` §1), and
Lemma 6.10 gives `E(g*g*) ≤ γ|⟨ξ₀|ξ⟩|² = c|ĝ(0)|² = 0`. Hence **`R_F ≥ 0` on that test space**,
and `Q_ζ(g) ≥ 0` with a structural reason. *Every ingredient of this is Connes–Consani's; the only
contribution here is the exact placement into the repository's `Q_ζ` and the verification that the
prime block is empty.*

**The smallest meaningful extension is not available from these sources.** Two candidates were
considered, as instructed:

1. **First active prime.** Requires `L > ½log2`. But `E(f) = ⟨ξ|N_I ξ⟩` and Lemma 6.10 are proved
   **on `H = L²(I)` with `I = [−½log2, ½log2]` specifically**, and Theorem 6.11's support
   hypothesis is `[2^{-1/2},2^{1/2}]`. Enlarging the window invalidates the archimedean comparison
   *and* switches on `W_2` — **two separate conditions, both failing at once**. Nothing in the
   sources extends Lemma 6.10 to a larger interval.
2. **The semi-local `{∞, 2}` setting.** arXiv:2008.10974 establishes that `u = ρ_∞ρ_2` is
   quasi-inner and that `S(u(F))` is infinite-dimensional and filtering. It **defines** the
   semi-local Sonin space as `ker U₂₂` and states explicitly that it is "meant to be a **first
   test** pertaining to the general **strategy**", **postponing to a future paper** even the proof
   that this definition reproduces the semi-local analogue.

   > **Amended 2026-10-07.** The sentence that stood here — *"There is no semi-local analogue of
   > Theorem 4.6 or Theorem 6.11: no semi-local positive trace, and no semi-local identified
   > discrepancy"* — conflated the **geometry** with the **estimate**, and its first clause is
   > withdrawn as to the geometry. The postponed identification has appeared:
   > **CCM arXiv:2310.18423v2 (2024)** defines `S_λ(X_S,α)` directly (Def. 4.5), proves
   > `θ_S(S_λ(R,e_∞)) ⊂ S_λ(X_S,α)` and `F_S θ_S = θ_S F̃_R` (Props. 4.6–4.7), and proves `θ_S` is
   > a **hilbertian isomorphism** of Sonin spaces, stable under enlarging `S` (its Thm 4.6, p. 23 =
   > intro Thm 2, p. 5). CCM2024 §1 further reports a **general semi-local trace formula** in
   > A. Connes, *Selecta Math.* **5** (1999), 29–106 (**NOT CHECKED** here).
   >
   > **What remains missing is the second clause, now stated precisely:** no semi-local analogue of
   > **CC2021 Theorem 4.6** (an identity `Tr_S = W_S + E_S` with `E_S` defined *independently* of
   > `W_S`) and no semi-local analogue of **Theorem 6.11** (a bound on `E_S` with explicit test
   > space, support, normalisation and finite-place dependence) was found in the checked sources.
   > CCM2024 §1 calls this a *"more precise strategy for addressing the semilocal Weil positivity"*
   > which the authors *"expect"* will *"open a way"* — **proposed, not proved**.
   >
   > **And a structural caution.** "Hilbertian" there is **not** unitary: CCM2024 footnote 2 reads
   > *"We use the term 'hilbertian' to denote the underlying topological vector space structure of
   > a Hilbert space."* On the critical line `θ_S` multiplies the Mellin transform by
   > `∏_{p∈S\{∞}} L_p(½+is)^{-1}` (eq. 58), and `⟨θ_S f | η_S g⟩ = ⟨f|g⟩` (Prop. 4.7(iii)) makes
   > `η_S = (θ_S^*)^{-1}`, not `θ_S^{-1}`. A bounded invertible `J` sends an orthogonal projection
   > `P` either to `JPJ^{-1}` (idempotent, not self-adjoint — no sign) or to `JPJ^*` (positive, not
   > idempotent — not the Sonin compression). **The archimedean "positive for free" does not
   > transport along a non-isometric map.** The distortion is bounded for fixed `S` but not
   > uniformly: condition numbers `5.83, 21.75, 56.95, 126.15, …` for `S\{∞\} = {2}, {2,3},
   > {2,3,5}, {2,3,5,7}, …`, diverging because `Σ_p p^{-1/2}` does.

> **Therefore no semi-local `P_F` is written down here.** Inventing one would be exactly the error
> the brief forbids. The arithmetic step is blocked at the source, not by an obstruction we found.
> *(Unchanged 2026-10-07: the 2024 geometry does not license writing one down, because the
> positivity and the discrepancy bound are the parts still missing.)*

**Comparison with existing certificates, before any claim of improvement.** The archimedean-only
positivity established at `L = ½log2 ≈ 0.3466` sits **strictly inside** the window `L ≤ 0.8` where
the manuscript already quotes Zhu's certified `λ*(0.8) ≥ 8.9×10^{-18}` — and Zhu's window **does**
contain active primes (`n = 2,3,4`). The two are different kinds of statement and neither subsumes
the other: Connes–Consani supply a *structural* reason plus an operator inequality valid for a
class of `g`; Zhu supplies a *numerical certificate* for the window form's minimum eigenvalue.
**But as a positivity range, the intrinsic construction does not currently reach as far as the
certificate already in use, and it is archimedean-only where the certificate is not.** No
improvement is claimed.

## 5. Does the repository's projection work contribute? Diagnostically only.

Read against the corrected Theorems 7–8 and Remark 10:

- Our **Theorem 7** is a *relative* perturbation budget, not an independent positivity baseline.
  Connes–Consani's `Tr(ϑ(g)Sϑ(g)*)` is exactly the kind of **independent** baseline Theorem 7
  lacks — but it is theirs, and only at their window.
- Our **Theorem 8** remains conditional on all remaining zeros being on the line; nothing here
  discharges it.
- `K(L,γ)` is a derivative-evaluation **norm**, not the attained coefficient of the minimum
  eigenvalue (Remark 10 (a)/(b)/(c)). Nothing here changes that.
- Our `E = ker A_0` is an **exact** kernel in finite dimension; the Sonin `S` is an infinite-rank
  projection on `L²(R)_ev`. **No transfer of the quartic toy-model formula is made**, and the
  hypotheses do not match.
- Our results are finite-dimensional perturbation theorems; `W_∞` is an unbounded quadratic form.

**What the projection work does supply** is the vocabulary that makes Lemma 6.9 legible — magnitude
(`a, b, c`), alignment (`|⟨φ|ψ⟩|²`), response — which is precisely the interpretive principle
recorded in `PROJECTION_CONTEXT.md` §4. That is an **interpretation of their theorem, not an
estimate**, and it is credited accordingly. **No positivity estimate here comes from our
sensitivity law.**

Closed routes (7/8 transfer, prime-comb, scalar-coupling, nested-window inheritance) were **not**
reopened; no new lemma here changes their premises. Their failure is not treated as a no-go theorem
for positive compression.
