# notes/prop61_dependency_chain.md — the three cited inputs of Prop 6.1, run to ground

Phase 2 upgrade. The companion note `prop61_dependency.md` left the verdict **CONDITIONAL**
because three load-bearing inputs of Prop 6.1 lived in papers not then obtained. Both have
now been read.

**Sources.** They were *not* found in `~/Downloads` (only `2610.32704v1.pdf` was there; a
search of `~`, `~/Desktop`, `~/Documents` found nothing). Both were downloaded from arXiv
and identity-checked against title and author line:

| short | arXiv | title / author | pages |
|---|---|---|---|
| **[S]** | 2606.09096 | *Weil's quadratic form via the screw function*, M. Suzuki | 33 |
| **[CC]** | 2106.01715 | *Spectral triples and ζ-cycles*, A. Connes, C. Consani (= Enseign. Math. **69** (2023), 93–148) | 56 |

So that the exact bytes read here can be checked against anyone else's copy:

```
sha256  75b8512f7081bd69152781c7a20f3b5dad4374fe0dfc80f3014f9223f044c81a   2606.09096.pdf   491037 bytes
sha256  784000e34f3efa4e8fe455f849d9b7e195f2c7b3ee0a0a073836cbd457543280   2106.01715.pdf  2743824 bytes
```

(arXiv may re-issue a version at the same identifier; a hash mismatch means a different
PDF, not necessarily different mathematics — re-check the page numbers below if it occurs.)

**Reference mapping, verified against the printed bibliography of 2610.32704v1** (not
assumed): its `[2]` is Connes–Consani, Enseign. Math. **69** (2023) = **[CC]**; its `[10]`
is Suzuki, *Weil's quadratic form via the screw function*, arXiv:2606.09096 = **[S]**.

**Method.** Per the standing instruction, every statement below was read from a **220/200 dpi
page render**, not from text extraction. This is not pedantry: Retraction 25 in this same
extension was caused by `pdftotext` and `pypdf` silently dropping an overline in Suzuki's
(6.4). Pages rendered and read: **[S]** 3, 5, 12, 13, 14, 15; **[CC]** 11, 13; and
2610.32704v1 page 5.

---

## The blanket statement

**[S]**, §1.2 *Main Results*, page 3, first line, read from the render:

> "Before presenting the main results, we emphasize that **none of them depends on RH**."

This is recorded but **not relied on**. A blanket claim about "the main results" does not by
itself settle whether a *particular lemma inside a proof* is RH-free, so each of the three
inputs was traced to its own proof. The blanket statement is corroboration, not evidence.

Worth noting for contrast, from the same page: **[S]** describes `[14]` (Suzuki, Canad. J.
Math. **75**) as having built the global Hilbert space — the completion of `C_c^∞(ℝ)` under
`Q_W` — **"under RH"**. So the author distinguishes the two cases explicitly, and the
localized `(−a,a)` results are on the unconditional side of his own line.

**Corroborating scan.** Across the whole of §§2–3 of **[S]** (pages 8–16), which contain
all three inputs, the string `RH` occurs **exactly once**: page 13, in the sentence
identifying RH as *equivalent to* `Γ ⊂ ℝ` — i.e. flagging it as a hypothesis **not** being
assumed. No step in §§2–3 invokes it.

---

## Input 1 — `[2, Lemma 2.2]`: the exponential span is a form core

**Where Prop 6.1 uses it.** 2610.32704v1 p. 20: *"By [2, Lemma 2.2], the linear span of
`{e^{iπnx/a} : n ∈ ℤ}` is a core for the quadratic form `Q_{W_a}`."*

**Found at [CC] page 11, Lemma 2.2**, read from the render:

> **Lemma 2.2.** *Let `λ > 1`, and `U ∈ L²([λ^{-1},λ], d*u)` be the function
> `U(u) := u^{iπ/log λ}`. Then the space of Laurent polynomials `ℂ[U,U^{-1}]` is a core for
> the quadratic form `QW_λ`.*

**Coordinate match.** **[CC]** work multiplicatively on `[λ^{-1},λ]`; **[S]** and
2610.32704 work additively on `(−a,a)`. **[CC]** perform the change of variable themselves
on the same page: *"We switch from `ℝ*₊` to the additive group `ℝ` (using the logarithm)
and let `L = 2 log λ` and `ℋ = L²([−L/2,L/2])`. Under this change of variable the function
`U` becomes `U(x) = exp(2πix/L)`."* With `a = log λ`, `L = 2a`, so
`U(x) = exp(2πix/2a) = exp(iπx/a) = e₁(x)`, and `ℂ[U,U^{-1}] = span{e^{iπnx/a} : n ∈ ℤ}`
on `(−a,a)`. **The two statements are the same statement.**

**RH-dependence of the proof.** The proof is a Fourier-analytic density argument and
nothing else. Its ingredients, in order: the Fourier transform of `U^n` evaluated in closed
form,
`Û^n(s) = 2 log λ (−1)^n sin(s log λ)(πn − s log λ)^{-1} = O(1/|s|)`  (2.13);
the representation `Q_∞(f,f) = ∫|f̂(t)|² · 2∂_tθ(t)/2π dt`  (2.15), with `θ` the
Riemann–von Mangoldt angle, so the weight `∂_tθ(t) ≍ log|t|` is a **smooth archimedean
weight**, not a sum over zeros; the elementary bounds `sin²(s)s^{-2} ≤ 2(1+s²)^{-1}` and
`1+(s−a)² ≤ 2(1+a²)(1+s²)`, giving `‖Û^k‖²₁ = O(log|k|)`; and rapid decay of Fourier
coefficients of a `C^∞` function on the circle (**[CC]** p. 13).

> **Verdict: RH-FREE.** No zero of `ζ` is located, counted, or assumed real anywhere in the
> argument. The only arithmetic input is the density `∂_tθ`, which is unconditional.

**Bonus found on the same page — [CC] Corollary 2.4**, which Phase 3 now cites:

> **Corollary 2.4.** *The lower bound of `QW_λ` is the limit, when `N → ∞`, of the smallest
> eigenvalue of the restriction of `QW_λ` to the linear span `E_N` of the functions `U^k`
> for `|k| ≤ N`.*

This is the unconditional justification that `λ_min(G_{a,N})` is a **variational upper
bound on `λ_a` decreasing to it** — exactly the caveat Phase 3 attaches to its matched-window
comparison. It is a corollary of Lemma 2.2 + Prop 2.3 and is likewise RH-free.

---

## Input 2 — `[10, §3.2, proof of Thm 1]`: `C_c^∞(−a,a)` is a form core

**Where Prop 6.1 uses it.** 2610.32704v1 p. 20: *"By [10, Section 3.2, proof of Theorem 1],
`C_c^∞(−a,a)` is also a core for this form"*; and p. 5, in establishing
`ℋ(A_a) = 𝔇(A_a^{1/2})`.

**Found at [S] §3.2 "Completion of the proof", pages 14–15.** §3 of **[S]** is titled
*Proof of Theorem 1.1*, so the citation's "Theorem 1" is **[S]**'s Theorem 1.1. The
structure, read from the renders:

1. `B_a` is the operator on `H¹₀(−a,a)`; `A_a` the form operator. **Lemma 3.1** (p. 14)
   proves `B_a ⊂ A_a`. Its proof uses Input 1 — *"The space `E` is a core for the quadratic
   form `Q_W^a` [3, Proposition 2.3]"* — plus `ê_n(z) = 2sin(az+nπ)/(z+nπ/a)`, the decay
   `v̂(z) = O(|z|^{-1})` for `v ∈ H¹₀`, formula (3.1), and Schmüdgen's
   [12, Prop. 10.5(v)]. **Note:** here the Connes–Consani core is cited as
   **[3, Proposition 2.3]**, whereas 2610.32704 cites **[2, Lemma 2.2]**. Both are in
   **[CC]** §2.1, and Prop 2.3 is proved *from* Lemma 2.2 (**[CC]** p. 13: *"The density of
   `ℂ[U,U^{-1}]` in the domain of `QW_λ` for the graph-norm shown in Lemma 2.2"*). The two
   citations point at the same RH-free fact; the numbering difference is not a discrepancy.
2. §3.2 then needs only the inclusion
   `(3.2)  E ⊂ closure of C_c^∞(−a,a) in ‖·‖_Q`,
   from which `𝔇(Q̄_B) = 𝔇(Q)` and the Friedrichs extension of `B` equals `A`.
3. **(3.2) is proved explicitly on p. 15** by a cutoff construction: `η_ε = 1` on
   `[−a+ε, a−ε]`, `η_ε(±a) = 0`, `|η'_ε| ≪ ε^{-1}`; set `v_ε = e_n η_ε ∈ C_c^∞(−a,a)` and
   `w_ε = e_n(1−η_ε)`. Then `‖w_ε‖²_{L²} ≤ 2ε`; `|ŵ_ε(z)| ≤ 2ε` (3.3); integration by parts
   gives `|ŵ_ε(z)| ≤ M/|z|` (3.4) uniformly in `ε`; splitting the `I_ε` integral at
   `|z| = 1/ε` gives `≪ ε log(1/ε)` on both pieces, hence
   `‖e_n − v_ε‖²_Q ≪ ε log(1/ε) → 0`.

> **Verdict: RH-FREE.** The argument is a cutoff, Plancherel, one integration by parts, and
> a dyadic split. It locates no zero and assumes nothing about `Γ`. Its only non-elementary
> input is Input 1, already cleared.

---

## Input 3 — `[10, §3]`: the form-domain characterization

**Where Prop 6.1's chain uses it.** 2610.32704v1 p. 5: *"By the characterization of the form
domain in [10, Section 3] and the identification `ℋ(A_a) = 𝔇(A_a^{1/2})` above, we have the
set-theoretic equality `ℋ(A_a) = {u ∈ L²(−a,a) : ∫_ℝ (1+log⁺|t|)|û(t)|² dt < ∞}`."*

**Found at [S] pages 13–14.** Read from the renders, the pieces are:

- **Definition** (p. 13): `𝔇(Q_W^a) := {v ∈ L²(−a,a) : |Q_W^a(v)| < ∞}`, a Hilbert space
  under `‖v‖²_{Q_W^a} := Q_W^a(v) + (1−λ_a)‖v‖²_{L²}`. Definitional; the shift `(1−λ_a)`
  makes the norm positive and is exactly where `λ_a` enters.
- **The space** (p. 13):
  `H^log(−a,a) := {v ∈ L²(−a,a) : ∫_ℝ (1+log⁺|z|)|v̂(z)|² dz < ∞}`.
- **The inclusion** (p. 14, first line): `𝔇(Q_W^a) ⊂ H^log(−a,a)`, "by (2.6) and (2.8)".
- **Form-operator characterization** (p. 14): `v ∈ 𝔇(A_a) ⊂ 𝔇(Q_W^a)` iff `w ↦ Q_W^a(v,w)`
  is `L²`-continuous on `𝔇(Q_W^a)`; equivalent to Schmüdgen [12, Def. 10.4], with
  `A_a v` characterized by `Q_W^a(v,w) = ⟨A_a v, w⟩_{L²}`. Standard form theory, and
  [12, Cor. 10.8] gives the bijection between lower-bounded self-adjoint operators and
  lower-bounded closed forms with dense domain.
- **Lower-boundedness and lower semicontinuity of `Q_W^a`** (p. 14) is attributed to
  **[3, Proposition 2.1]** = **[CC]**, "based on the formula (2.8)", lower semicontinuity
  by Fatou.

**Inclusion vs. equality — checked, because it matters.** **[S]** states only `⊂`, while
2610.32704 asserts a *set-theoretic equality*. The equality is nevertheless supported, and
the reason is **[S]** (2.8) on page 12:

> `Q_W^a(v) = (1/2π) ∫_{−∞}^{∞} (log|z| − log(2π)) |v̂(z)|² dz + O(a)·‖v‖²_{L²}`

This is an **identity with a bounded error term**, not a one-sided estimate. For `v ∈ L²` the
`O(a)‖v‖²` term and the bounded-`|z|` region are finite regardless, so
`|Q_W^a(v)| < ∞ ⟺ ∫ log⁺|z| |v̂|² dz < ∞ ⟺ v ∈ H^log`. The equality follows from (2.8)
read in both directions; **[S]** only needed, and only stated, the direction he used.

(**[S]** p. 14 separately notes that the *Gagliardo-type* space of
[2, Thm 5.1] — finiteness of `∫∫|v(x)−v(y)|²/|x−y| + ∫|v|²` — is "slightly larger than
`𝔇(Q_W^a)`". That is a different space from `H^log` and is **not** the characterization being
cited; it is recorded here so the two are not conflated.)

**RH-dependence.** (2.8) is assembled on p. 12 from the *geometric* side of the explicit
formula only: `Re ψ(1/4 + iz/2)`, `−log π`, the **finite** prime sum over `n ≤ e^{2a}`, and
the remainder `r''_{0,a}` with `r̂''₁(z) = O(|z|^{-2})`. **No zero of `ζ` appears in it.**

> **Verdict: RH-FREE.**

---

## Result

| # | input | source, located | RH-free? |
|---|---|---|---|
| 1 | `[2, Lemma 2.2]` — exponential span is a form core | **[CC]** p. 11 (Lemma 2.2), proof completed p. 13 | **yes** |
| 2 | `[10, §3.2, proof of Thm 1]` — `C_c^∞(−a,a)` is a form core | **[S]** pp. 14–15 (Lemma 3.1 + (3.2) cutoff proof) | **yes** |
| 3 | `[10, §3]` — form-domain characterization | **[S]** pp. 12–14 ((2.8), `H^log`, Schmüdgen form theory) | **yes** |

**X = 3 → the Phase 2 condition is met.**

---

## What the PASS is, stated exactly

> **PASS (scoped).** At `a = ½log 3 = 0.549306…`, where `λ_a > 0` is supplied by Zhu's
> certified `λ*(0.8) ≥ 8.9e−18` together with monotonicity of `λ*` in the window
> (`0.5493 < 0.8`, and a larger window minimizes over a larger space, so
> `λ*(0.5493) ≥ λ*(0.8) > 0`), **the conclusions of Prop 6.1 hold without assuming RH**:
> `u^{(N)}_{z,a} → u_{z,a}`, `J_{a,N} → J_a`, `W_N(a,π;·) → W(a,π;·)` locally uniformly,
> and the normalized limit (6.5).

Every step of Prop 6.1's proof factors through `(2.4)` alone (the step-by-step table is in
`prop61_dependency.md` §1), and the three inputs that could not previously be inspected are
now each verified RH-free at their source. The conditional clause of the earlier verdict is
discharged.

### Three qualifications that survive the upgrade, and are not erased by it

1. **Prop 2.1's third de Branges axiom still needs `Γ ⊂ ℝ`.** This was recorded in
   `prop61_dependency.md` §3 and is now **confirmed from the 2610.32704 p. 5 render**. The
   membership `F_w ∈ Ĥ(A_a)` is fine unconditionally — it goes through the `H^log`
   characterization, whose integral is over the **real** axis, where `|F_w(t)| = |F(t)|`
   holds because the Blaschke factor has modulus 1 there. But the *norm* equality
   `‖F_w‖ = ‖F‖` is justified by "the definition of `Q_W`", and the exact `Q_W` is the
   zero-sum **[S]** (3.1) `Q_W(v₁,v₂) = Σ_{γ∈Γ} m_γ v̂₁(γ) conj(v̂₂(γ̄))`, where **[S]** p. 13
   says in terms that `γ ∈ Γ` *may be non-real*. At a non-real `γ` the factor does not have
   modulus 1 and the norms differ. **So this step uses RH itself, not `(2.4)`.** It does not
   bite, because Prop 6.1's chain uses `J_a` only through its definition (2.6) as
   `⟨u_{w̄,a}, u_{z̄,a}⟩` and never invokes the third axiom, the generating function `E_a`,
   or (2.7); Prop 6.1 uses Prop 2.1 only for the *no-common-zeros* clause, which is proved
   RH-free (choose a test function whose transform does not vanish at `z`). It is kept on
   the record because it shows the reduction to `(2.4)` is **not automatic** — the paper
   does contain a step that needs `Γ ⊂ ℝ`.
2. **The window `a = 0.8047` (`X = 5`) is outside Zhu's certified range** `L ≤ 0.8`, and
   `X = 7, 11` further still. The PASS is for `X = 3` **only**. `X = 5, 7` remain labelled
   conditional on Liu's unrefereed certificate; `X = 11` remains numerical only.
3. **Nothing here is about `a → ∞`.** 2610.32704's Theorems 1.1/1.2 rest on its Thm 3.3 and
   the global space `H_W`, which needs `Q_W ≥ 0` on all of `C_c^∞(ℝ)` — i.e. RH — and
   **[S]** itself describes that global construction as made "under RH". No fixed-window
   certificate relaxes them.

### What the certificate actually buys, given [S] Theorem 1.4

**[S]** Theorem 1.4 (p. 5, read from the render) already gives unconditionally:

> *For sufficiently small `a > 0`, the lowest eigenvalue `λ_a` is positive, simple, and
> satisfies `λ_a = log(1/a) + μ₁ − log(2π) + ψ(2) − 1 + O(a)` as `a → 0+`, for some constant
> `μ₁ > 0`.*

So `λ_a > 0` on *some* window is not new (cf. also Yoshida [17, Lemma 2]). But that
statement is **ineffective in three places at once**: the threshold "sufficiently small" is
unquantified, `μ₁` is unnamed, and the `O(a)` has no stated constant. It therefore cannot be
invoked at `a = 0.549`. That the asymptotic is far out of range there is visible
numerically: its explicit part is `log(1/0.5493) − log(2π) + ψ(2) − 1 = −1.8158`, an `O(1)`
quantity, against a measured `λ_min(G_{0.5493,6}) ≈ 1.04e−7` — eighteen orders from the
`8.9e−18` floor and nowhere near the regime where the expansion describes `λ_a`.

> **This is the whole content of the upgrade:** the certified bound supplies an *effective*
> instance of `(2.4)` at a *named* window where the unconditional asymptotic theorem is
> silent. It does not prove positivity of `Q_W` for all `a`, which is equivalent to RH by
> Weil's criterion — and **[S]** p. 5 restates that equivalence explicitly ("the failure of
> RH is equivalent to the existence of some `a > 0` for which `λ_a < 0`").
