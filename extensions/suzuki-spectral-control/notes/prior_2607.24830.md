# notes/prior_2607.24830.md — what the prior numerical realization already covers

Read before Phase 4, as PREREG requires. **Source inspected: the arXiv abstract page only**
(metadata + abstract). The full PDF was **not** obtained, so everything below is at the
level of the abstract's own claims and must be treated as a bibliographic record, not as a
verified reading. Where this extension's scope depends on what that paper does *not* do,
that is flagged as an inference from the abstract.

## The paper

**T. Kim, Y. Hong, M. Kim, S. Choi, J. Jang, M. Kim**, *A Numerical Realization of
Suzuki's Weil-Quadratic-Form Operator: The Archimedean Spectral Law, its Universality, and
an Operator Form of Weil's Positivity Criterion*, [arXiv:2607.24830](https://arxiv.org/abs/2607.24830).
Submitted 23 Jul 2026 (v1), last revised 28 Sep 2026 (v3). 18 pages.

It realizes Suzuki's **earlier** operator — the one from
[arXiv:2606.09096](https://arxiv.org/abs/2606.09096), reference [10] of the paper this
extension studies — via **P1 finite-element discretization with Richardson
extrapolation**, to **30-digit precision**.

## Its nine stated results

| | claim (as the abstract states it) |
|---|---|
| R1 | the spectrum follows an archimedean law with closed form, `A_k(a) = log(1/a) + log(k−2) + B_0 + O(a)`, `B_0 = log q − 2log 2`, to 30 digits |
| R2 | head coefficient depends on the conductor only; the degree of the `L`-function appears as a logarithmic slope |
| R3 | total spectral intensity follows the prime number theorem |
| R4 | **nontrivial zeros occur in explicit-formula error terms, not as eigenvalues** |
| R5 | the best-match line descends toward the critical line |
| R6 | **"Weil's positivity criterion is realized in operator form: bounded residual growth corresponds to all zeros on the line"** — an injected off-line zero causes exponential blow-up |
| R7 | the lowest eigenvalue is strictly positive with superexponential decay |
| R8 | **the characteristic function `W(a,0;z)` is computed for the first time, with all zeros confirmed real** |
| R9 | indirect traces of GUE statistics in the moment structure |

The authors state the work does not prove RH and that the results are "archimedean and
universal".

## What it already covers, and therefore what this extension must not re-claim

1. **An injected off-line zero is already studied** (R6), and the measured response is
   **exponential blow-up of a residual**. So "inject an off-line zero and see the
   construction react" is **not new**.
2. **A characteristic function with all zeros real is already computed** (R8) — so
   "compute `W` numerically and find real zeros" is **not new** either.
3. **R7** already reports strict positivity of the lowest eigenvalue with superexponential
   decay, which overlaps this repository's own `λ*(L)` material and Zhu's law.

## What this extension adds that the abstract does not claim

Each item is a difference in the *object*, not a difference in presentation.

1. **The phase.** R8 computes `W(a,**0**;z)`. Suzuki's Theorems 1.1 and 1.2 and his §6
   are all at the **distinguished phase `θ = π`**, which is a *different self-adjoint
   extension* — by Lemma 3.1, `θ = 0` has `0` as an eigenvalue and `θ = π` does not, and
   only `θ = 0, π` have reflection-invariant domains. This extension works at `θ = π`.
2. **The discretization.** That paper uses **P1 finite elements + Richardson
   extrapolation** on the operator. Suzuki's §6 uses a **Riesz-vector approximation in the
   exponential basis `χ_{n,a}` against the arithmetic Gram matrix `G_{a,N}`** — a different
   finite-dimensional scheme, and the one his Prop 6.1 is about. This extension implements
   §6.
3. **Suzuki's own tables.** Reproducing Tables 1 and 2 of arXiv:2610.32704v1 (`d_6^{(N)}`,
   `κ₂(G_{a,18})`) is specific to the newer paper and cannot be in a July 2026 preprint.
4. **The forced-zero lattice.** `W_N(kπ/a) ≡ 0` for `|k| > N`, and the identification of
   Suzuki's reported extra zero of `W_14` at `30π/log 11 ≈ 39.3044` as lattice point
   `k = 15`. Specific to the `χ_{n,a}` basis, so it has no analogue in a finite-element
   scheme.
5. **A function-field / Frobenius analogue.** The abstract mentions **no** function-field,
   Frobenius, or Hallouin–Perret construction. Phase 4's attempt to build the `W_N`
   analogue on the Hallouin–Perret Gram matrix — where RH *is* Weil's theorem — appears
   nowhere in it.
6. **The quantitative coincidence test.** R6 reports that an injected off-line zero makes
   a residual blow up exponentially. It does **not**, per the abstract, locate a
   *threshold* in the displacement `δ`, nor compare that threshold with the one at which
   positivity is lost, nor tie either to this repository's `δ²`–`L³` law with its constant
   `C ≍ (window length)³/6`. Phase 4c's preregistered prediction — that loss of root
   reality coincides with loss of positivity at the `δ²`–`L³` threshold, to 10% — is a
   statement about **two thresholds agreeing**, which is a different measurement from
   "blow-up happens".
7. **Certified-window unconditionality.** Whether Prop 6.1 survives with *certified*
   coercivity (Zhu `L ≤ 0.8`; Liu, unrefereed, at half-width 1 and 17/16) in place of
   assumed RH is a question about Suzuki's §6 proof, not about a numerical realization.

## Caution carried into Phase 4

R4 — *"nontrivial zeros occur in explicit-formula error terms, not as eigenvalues"* — is,
if correct for that operator, a warning directly relevant to this extension's
non-claim 2: **a root of a finite-window spectral function is a local spectral value, not
a zeta zero.** It is recorded here as independent support for that non-claim, from a
source this extension did not write.
