# notes/prop61_dependency.md — where RH enters Prop 6.1, step by step

Phase 2 deliverable. **The claim under test, not assumed:** *on windows where `λ_a > 0` is
certified rather than inferred from RH, Prop 6.1 holds unconditionally.* PREREG fixed an
honesty condition: any use of RH not reducible to `(2.4)` withdraws the claim for that part.

Reading of arXiv:2610.32704v1. Three cited inputs live in papers that were **not obtained
when this note was first written** and are marked **[UNVERIFIED]** below.

> **SUPERSEDED IN PART.** Both sources have since been obtained and read from page renders.
> All three inputs are **RH-free**, and the verdict in §4 is upgraded from CONDITIONAL to a
> scoped **PASS at `X = 3`**. See **`prop61_dependency_chain.md`** for the located
> statements, the proofs, the verified reference mapping, and the qualifications that
> survive. The `[UNVERIFIED]` marks below are left standing as the record of what was
> unknown at the time; they are discharged, not deleted.

## 0. What is unconditional before any of this

`A(z) = ξ(1/2 − iz)` is **real-entire and even with no RH input**: `ξ(1−s) = ξ(s)` gives
`A(−z) = A(z)`, and `ξ(s̄) = conj ξ(s)` gives `A(z̄) = conj A(z)`. Hence `Γ = −Γ`,
`Γ = conj Γ`, with matching multiplicities, unconditionally. What RH adds is `Γ ⊂ ℝ`.

The **arithmetic** form is also unconditional: `A_a` on `L²(−a,a)` is defined by the
geometric side of the explicit formula (archimedean + pole + finitely many prime powers),
with no reference to any zero. So `λ_a := inf σ(A_a)` is a well-defined real number with no
RH input, and `G_{a,N}` is computable — which is Suzuki's own point in §6.

**The role of RH in the construction is therefore located precisely:** §2.2 says *"Under RH,
we may take `λ = 0`, obtaining the Hilbert space `H(A_a)`"*, against the earlier
construction of `H(T_a)` with `T_a = A_a − λI` for `λ` below `inf σ(A_a)`, which the
introduction states **"hold[s] without assuming RH"**. So RH is used *to supply
`λ_a > 0`* — exactly the hypothesis a certificate can replace.

## 1. Prop 6.1's own proof, line by line

| step | what it uses | verdict |
|---|---|---|
| "By [2, Lemma 2.2], the span of `{e^{iπnx/a}}` is a core for `Q_{W_a}`" | Connes–Consani, Enseign. Math. **69** = arXiv:2106.01715 | ~~[UNVERIFIED]~~ → **RH-free**, verified at p. 11 |
| "By [10, §3.2, proof of Thm 1], `C_c^∞(−a,a)` is also a core" | Suzuki arXiv:2606.09096 | ~~[UNVERIFIED]~~ → **RH-free**, verified at pp. 14–15 |
| "Hence, **under RH**, (2.4) and the definition of `H(A_a)` imply `∪_N X_{a,N} = H(A_a)`" | the two cores + `(2.4)` + the completion being well defined | **via (2.4) only** — given `λ_a > 0`, `Q_W` is positive definite on `C_c^∞`, so the completion exists and the cores are dense in it |
| "Uniqueness in the Riesz representation gives `u^{(N)}_{z,a} = P_{a,N}u_{z,a}`, `P_{a,N} → I` strongly" | Riesz in a Hilbert space | **RH-free** |
| "`(2.4)` and `(2.5)` imply `‖u_{z,a} − u_{w,a}‖ ≤ λ_a^{−1/2}‖e_z − e_w‖_{L²}`" | `(2.4)` explicitly | **via (2.4) only** |
| "the strong convergence of `P_{a,N}` is uniform on `{u_{z,a} : z ∈ K}`" | the previous line + compactness | **via (2.4) only** |
| "`J_a − J_{a,N} = ⟨(I−P)u_{w̄}, (I−P)u_{z̄}⟩`, then Cauchy–Schwarz" | orthogonality in `H(A_a)` | **RH-free** |
| "the assertions for `W_N` and its normalization follow from (3.3), (6.4), Lemma 3.2" | §3.2 and Lemma 3.2 | see §2 |

**Prop 6.1's proof does not cite Theorem 3.3.** That is the one result in the paper whose
RH-dependence is irreducible (it builds the *global* space `H_W`, which needs `Q_W ≥ 0` on
all of `C_c^∞(ℝ)`, i.e. RH, and the sampling isomorphism `S : H_W → ℓ²(Γ,m)` of `[9]`).
Theorem 3.3 is used for **Theorem 1.1** — the `a → ∞` limit — not for Prop 6.1.

## 2. The supporting results Prop 6.1 does invoke

**`(2.5)`** `⟨u, u_{z,a}⟩_{A_a} = û(z̄)`: the Riesz representative exists because point
evaluation is bounded, and the bound *is* `|û(z)| ≤ (∫_{−a}^{a}e^{−2x Im z}dx)^{1/2}
λ_a^{−1/2}‖u‖_{A_a}` — **(2.4) only**.

**Lemma 3.1** (`D(D̄_{a,θ})` reflection-invariant iff `θ ∈ {0,π}`; `0` an eigenvalue of
`D̄_{a,0}` but not `D̄_{a,π}`): uses `Rw_{a,θ} = e^{iθ}w_{a,−θ}`, the evenness of `u_{0,a}`
and `w_{a,0}`, and *"for a sequence of test functions converging in the graph norm,
applying `(2.4)` to both the functions and their derivatives gives convergence in `H¹` and
zero boundary values"* — **(2.4) only**.

**Lemma 3.2** (`W(a,π;·)` even, `W(a,π;0) ≠ 0`, zeros simple = eigenvalues of `D̄_{a,π}`):
uses reflection symmetry of `J_a`, Lemma 3.1, the von Neumann decomposition (which needs
`‖u_{i,a}‖ = ‖u_{−i,a}‖`, i.e. `J_a(i,i) = J_a(−i,−i)`, from reflection — RH-free), the
boundary form, and `(2.5)`. **via (2.4) only.**

**§3.2's identification** of zeros of `W(a,θ;·)` with eigenvalues: uses *"Proposition 2.1
gives `u_{z̄,a} ≠ 0`"* — the **no-common-zeros** clause, which Prop 2.1 proves by choosing a
test function whose transform does not vanish at `z`. **RH-free.**

## 3. The one place the reduction genuinely FAILS — and why it does not bite

**Proposition 2.1's third de Branges axiom** (invariance under the transformation of
nonreal zeros, *with preservation of the norm*) argues: for `F ∈ Ĥ(A_a)` with `F(w) = 0`,
`w ∉ ℝ`, set `F_w(z) = ((z−w̄)/(z−w))F(z)`; then *"since `|F_w(t)| = |F(t)|` on the real
axis … the **definition of `Q_W`** yields `‖F_w‖ = ‖F‖`."*

That step reads the norm as `Σ_{γ∈Γ} m_γ|F(γ)|²` and needs the sampling points to lie **on
the real axis**, where the Blaschke factor has modulus 1. **If `Γ ⊄ ℝ` the factor has
modulus ≠ 1 at the nonreal sampling points and the norms differ.** So this step uses
`Γ ⊂ ℝ` *essentially* — **not** through `(2.4)`.

> **This is the failure mode PREREG named in advance, and it is real.** It is recorded even
> though it does not bite, because it shows the reduction to `(2.4)` is **not automatic**:
> the paper does contain a step that needs `Γ ⊂ ℝ` itself.

Why it does not bite for Prop 6.1: the chain in §§1–2 never uses the third axiom. `J_a` is
used only through its definition `(2.6)` as `⟨u_{w̄,a}, u_{z̄,a}⟩_{A_a}`, which needs the
Riesz vectors and nothing more. The de Branges *generating function* `E_a` and the formula
`(2.7)` — which is where the axiom earns its keep — appear nowhere in Prop 6.1's proof or
in Lemmas 3.1/3.2.

## 4. Verdict

> **PASS (scoped to `X = 3`), upgraded from CONDITIONAL.** Every step of Prop 6.1's proof
> that can be inspected in arXiv:2610.32704v1 factors through `(2.4)` alone. The three
> cited inputs that could not previously be inspected — `[2, Lemma 2.2]` and `[10, §3.2]`
> (the two form cores) and `[10, §3]` (the form-domain characterization) — have now been
> located in their sources and are **each RH-free**. So on a window where `λ_a > 0` is
> **certified**, the conclusions of Prop 6.1 — `u^{(N)}_{z,a} → u_{z,a}`, `J_{a,N} → J_a`,
> `W_N(a,π;·) → W(a,π;·)` locally uniformly, and the normalized limit `(6.5)` — follow
> **without assuming RH**.
>
> The certified window is `a = ½log 3 = 0.5493`, via Zhu's `λ*(0.8) ≥ 8.9e−18` and the
> monotonicity of `λ*` in the window width. **`X = 5, 7` (a = 0.8047, 0.9730) lie outside
> Zhu's range and stay conditional on Liu's unrefereed certificate; `X = 11` stays numerical
> only.**

Full record, including the verified `[2]`/`[10]` reference mapping, the page renders each
statement was read from, and the inclusion-vs-equality check on the form domain:
**`prop61_dependency_chain.md`**.

**§3 above is NOT superseded.** Prop 2.1's third de Branges axiom genuinely needs `Γ ⊂ ℝ`,
now confirmed from the typeset page 5 of 2610.32704; it simply is not invoked by Prop 6.1's
chain. The PASS is a statement about that chain, not about the whole paper.

**What this does NOT give, in any case.** Theorem 1.1 and Theorem 1.2 — the `a → ∞`
statements — rest on Theorem 3.3 and the global space `H_W`, whose construction needs
`Q_W ≥ 0` on `C_c^∞(ℝ)`, i.e. RH itself. **No fixed-window certificate can relax them**, and
the roots of `W_N` at fixed `a` are local spectral values of a finite-dimensional Riesz
approximation — not zeros of `W`, and not zeros of `ζ`.

**Nor does it give positivity itself.** `λ_a > 0` for *some* window is already a theorem
without RH (`[10]` Thm 1.4, Yoshida [17, Lemma 2]); what those give is ineffective, and what
a certificate adds is an *effective* instance at a *named* `a`. Positivity for **all** `a`
is equivalent to RH by Weil's criterion and is untouched here.
