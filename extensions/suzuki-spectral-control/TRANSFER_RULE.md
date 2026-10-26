# The transfer rule

**A certified positivity window is an unconditional spectral result.**

Two literatures have been running next to each other without being joined in print: the
**positivity-certificate** line (Zhu; Liu; others), which proves `λ*(L) > 0` on explicit
windows by checkable computation, and **Suzuki's local spectral programme**
([arXiv:2610.32704](https://arxiv.org/abs/2610.32704)), which builds a Hilbert space, a
self-adjoint operator and a spectral function `W(a,θ;·)` on the interval `(−a,a)` and states
its results *"Assume RH"*. The observation here is that the RH hypothesis in the local part
of Suzuki's programme is used **only** to supply `λ_a > 0`. Anything that supplies `λ_a > 0`
some other way does the same job. A certificate is such a thing.

---

## Statement

> **Transfer rule.** Let `a > 0`. Suppose
>
> **(H1)** a lower bound `λ*(L₀) ≥ c > 0` is established by a **checkable certificate**, for
> some `L₀ ≥ a`; and
>
> **(H2)** the certificate's normalization is **reconciled** with Suzuki's — same window
> parameter, same quadratic form, same basis scaling.
>
> Then at that `a`, **without assuming RH**:
>
> 1. `ℋ(A_a)`, the `Q_W`-completion of `C_c^∞(−a,a)`, is a Hilbert space, and equals the
>    form domain `𝔇(A_a^{1/2})` with the `Q_W`-norm.
> 2. `𝒟_a u = iu'` on `C_c^∞(−a,a)` is densely defined and symmetric with **deficiency
>    indices `(1,1)`**, and von Neumann's theory gives its self-adjoint extensions
>    `D̄_{a,θ}`, `0 ≤ θ < 2π`.
> 3. The spectral function `W(a,θ;·)` is **entire of exponential type ≤ `a`**, and its zeros
>    are **real** and are exactly the eigenvalues of `D̄_{a,θ}`.
> 4. For `θ = π`: `W(a,π;·)` is **even**, `W(a,π;0) ≠ 0`, and **all its zeros are simple**
>    (Lemma 3.2).
> 5. The finite-dimensional approximations **converge**: `u^{(N)}_{z,a} → u_{z,a}`,
>    `J_{a,N} → J_a`, `W_N(a,π;·) → W(a,π;·)` locally uniformly, and the normalized limit
>    (6.5) holds (Prop 6.1).
>
> `λ*` is **non-increasing in the window width**, so a single certificate at `L₀` discharges
> **(H1)** for *every* `a ≤ L₀` at once.

**Why it holds.** Each step above was traced to its proof and found to use only
**positivity `(2.4)`** and **von Neumann self-adjointness** — never the de Branges structure.
Three cited inputs carry the chain and were previously unverifiable; all three are now
located in their sources and each is **RH-free**:

| input | source | verdict |
|---|---|---|
| `span{e^{iπnx/a}}` is a form core | Connes–Consani, arXiv:2106.01715 p. 11 | RH-free |
| `C_c^∞(−a,a)` is a form core | Suzuki, arXiv:2606.09096 §3.2 | RH-free |
| the form-domain characterization | Suzuki, arXiv:2606.09096 §3 | RH-free |

Full chain, with the page renders each statement was read from:
[`notes/prop61_dependency_chain.md`](notes/prop61_dependency_chain.md).

---

## What is *not* on the list

**Excluded, because it is not established by positivity and self-adjointness:**

- **"`Ĥ(A_a)` is a de Branges space"** — Prop 2.1 as a whole. Its *third* axiom, norm
  preservation under the nonreal-zero transformation, is justified from the definition of
  `Q_W` as the sum over `Γ`, and at a **non-real** `γ` the Blaschke factor does not have
  modulus 1. That step uses `Γ ⊂ ℝ`, i.e. RH itself. It is not invoked by 1–5 above.
- The generating function `E_a` and formula (2.7); Theorem 3.3; Theorems 1.1 and 1.2;
  (5.3), (5.4), Corollaries 5.1 and 5.2.

**It is local, and that is where the limitation bites.** `a → ∞` is where RH lives:
Theorems 1.1/1.2 run through the global space `H_W`, which needs `Q_W ≥ 0` on all of
`C_c^∞(ℝ)` — equivalent to RH by Weil's criterion. **No finite collection of certified
windows approaches that.** `λ_a > 0` for *some* window is, in fact, already an unconditional
theorem (arXiv:2606.09096 Thm 1.4); what it lacks is effectivity, and that is the whole of
what a certificate adds — an explicit instance at a *named* `a`. The rule converts
certificates into spectral statements at **finitely many, individually named** windows. It
does not accumulate toward the limit.

**The novelty is the observation and its verification, not new machinery.** Both ingredients
existed. Nothing here is a new theorem about `ζ`.

---

## Two hypotheses that are not formalities

**(H1) must be a certified *lower* bound.** A variational **upper** bound does not discharge
it — including this extension's own `λ_min(G_{a,N})` values, which are upper bounds at every
`N` (Connes–Consani Cor. 2.4). A non-detection certifies nothing. And **the direction of
monotonicity only runs inward**: a lower bound at `L₀` says nothing at any `a > L₀`. That
error has a number in this repository — Retraction 18 — and Phase 3 contains a worked
instance of *not* making it, at `a = 0.8047` against a certificate at `L = 0.8`.

**(H2) is a gate with teeth, not a box to tick.** Suzuki's `a`, the repository's `L`, and a
certificate's own window parameter are three conventions, and the factors relating the bases
are easy to get wrong in a way that leaves every number finite and plausible. This extension
ran that reconciliation as a blocking **Phase 0 gate** and it **failed on the first attempt**
— a factor `2L` in the gate's own test code (Retraction 23). It passes now, pinned to 49
digits against the repository's own assembly, with the conclusion `a = L`, `λ_a = λ*(L)`,
conversion constant exactly `1`.

---

## Instances

| window | `a` | certificate | status |
|---|---|---|---|
| **`X = 3`** | `½log 3 = 0.5493` | Zhu, `λ*(0.8) ≥ 8.9e−18`, and `0.5493 < 0.8` | **VERIFIED — unconditional** |
| `X = 5` | `½log 5 = 0.8047` | outside Zhu's range | conditional on Liu (**unrefereed**) |
| `X = 7` | `½log 7 = 0.9730` | outside Zhu's range | conditional on Liu (**unrefereed**) |
| `X = 11` | `½log 11 = 1.1990` | none | numerical only |

`X = 3` is the verified instance, and it is one instance. **A rule with one worked example is
an observation, not a method** — a second independent instance is what would make the
transfer rule load-bearing rather than anecdotal.

### Candidate second instance — not claimed, and what the audit must settle

The `jose-paul-b` repository is reported to certify positivity with autocorrelation support
up to `log 7`, i.e. half-width `a ≤ ½log 7 = 0.9730` — which by the monotonicity clause of
**(H1)** would cover `X = 7` **and** `X = 5`, taking three of Suzuki's four windows
unconditional. **This is recorded as a candidate and is asserted nowhere else in this
extension.** It is AI-assisted and checked only by automated audit, so it has not met the
standard the rest of this work is held to. Before it could be listed above, an audit would
have to settle, at minimum:

1. **Support convention.** Is the condition `log n < 2a` or `≤ 2a`, and is the endpoint term
   at `p^r = e^{2a}` included or excluded? Suzuki's §6 notes the contribution at `p^r = X`
   **vanishes**; a certificate that includes it, or that uses a half-open convention the
   other way, is certifying a *different* form.
2. **(H2) in full** — the Phase 0 reconciliation, run against the certificate's own
   normalization, to the same 49-digit standard, with the conversion constant exhibited.
3. **Direction and strictness.** Is the output a certified **lower** bound with an explicit
   positive constant, or a variational upper bound, or a non-detection?
4. **Whether the certificate is checkable independently of its generator** — re-derived from
   its stated inputs rather than re-run by its own code.

Until those are done, `X = 5` and `X = 7` stay exactly as the table says.
