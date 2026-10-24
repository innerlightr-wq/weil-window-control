# Source-defined deformation → fixed-window detector: a feasibility audit

A bounded analytical feasibility study. No manuscript, certificate, licence or existing
data file was touched. The audit itself was committed afterwards (`5a23735`), its lower
bound certified and guarded (`c97ae0a`, tag `lambda-lower-certify`), and the note of §10
drawn from it; nothing else in the repository was modified.

**Repository state at the start and end.** `weil-window-control` on `main` at
`e2cf7153ce6903be7c8560ef26c20ea449babbd9`, working tree carrying the three pre-existing
untracked paths (`extensions/channel-compatibility-exploration-2026-10-24/`,
`extensions/newtonian-weil-response-2026-10-08/`, two `paper/*.fdb_latexmk`/`.fls` build
products), all left untouched.

---

## Verdict, in plain language

**D — the bridge is blocked for a specified class of arguments, by a counterexample
inside this repository.**

The source-to-spectrum transfer itself **succeeded**: this deformation supplies genuine
function-level dynamics, and de Bruijn's theorem converts them into an unconditional
statement about the detector. The conclusion is not that source-defined dynamics cannot
help. It is narrower and more useful:

> **This deformation provides genuine function-level dynamics, but the general properties
> controlling that dynamics are too broad to distinguish ζ from Davenport–Heilbronn.**

The theorem that does the work never touches the arithmetic that separates the two
functions, and the Davenport–Heilbronn function satisfies every hypothesis a reversal
would use while **certifiably violating its conclusion**.

Concretely, three findings:

1. **A genuine forward constraint exists, and I can state it exactly.** De Bruijn's
   Theorem 7 needs only Hermitian symmetry, super-exponential decay, and a bounded zero
   strip — **not** kernel positivity (that is a hypothesis of the surrounding framework,
   not of the theorem; §2 of `SOURCE_AUDIT.md`). I derived the Davenport–Heilbronn kernel,
   verified it satisfies all three, and obtained: **every zero of `H^{DH}_τ` is real for
   `τ ≥ 2(σ_0 − ½)² = 1.602537483559822777787819`**, hence `Q^{DH}_τ ⪰ 0` there. For `ζ`
   the same computation returns `τ ≥ 1/2`, reproducing de Bruijn 1950 — which is how I
   know the normalization bookkeeping is right.

2. **That same statement blocks the sought inequality.** `Q^{DH}_0` on the `L = 2.0` odd
   level-2 space has `λ_min = −0.679071346441480590986377586093…`, enclosed in an Arb ball
   of radius `1.39×10⁻¹¹⁶`. So no **universally applicable** backward-preservation
   argument resting only on the shared hypotheses can work. The precise statement is not
   "the inequality cannot exist" but: along the required path from `τ_DH` down to `0`,
   **one of its ingredients must fail** — the inequality, the control on `a(θ)`, the
   differentiability, or the standing assumption `B_τ ≻ 0`. The channel representation
   needs `B_τ ≻ 0` at every `τ`, and zero-tracking can meet collisions; neither can be
   carried along the whole path for free.

   **Positive kernels are not enough.** `Φ_F` is not merely symmetric and rapidly
   decaying — it is positive. So the package
   `positive kernel + reflection symmetry + rapid decay + bounded zero strip` still does
   not deliver positivity at `τ = 0`. That strengthens the counterexample considerably.

3. **The arithmetic is not erased by the deformation — it is simply not used.** The
   deformed kernel is
   `Φ_{F,τ}(x) = e^{τx²} · 4e^{3x} Σ_{n≥1} n a_n e^{−πn²e^{4x}/5}`: the coefficients
   `a_n` are still there at `τ ≠ 0`. What has **not** been established is that the
   deformed family retains a usable Euler product, gamma-factor decomposition or
   prime-side explicit formula. And the only arithmetic input the whole construction
   actually consumes is the **initial strip width** — which is also all the toy family
   `z² + a² − 2τ` consumes, and that family has no arithmetic content at all.

The one place where the source equation does something genuinely different from the
previous explorations: its quadratic term **forces** a `Λ*Λ` divisor convolution rather
than deleting mixed terms by hand. But that series converges absolutely only for
`Re s > 1`, i.e. outside the strip where it would have to be used, and it enters `∂_τQ`
through a sign-indefinite factor. It is another exact, sign-indefinite expression.

---

## 1. The exact source object and parameter domain

```
Φ(u)   = Σ_{n≥1} (2π²n⁴e^{9u} − 3πn²e^{5u}) exp(−πn²e^{4u})          [S1 eq. (3)]
H_τ(z) = ∫_0^∞ e^{τu²} Φ(u) cos(zu) du                                [S1 eq. (4)]
H_0(z) = (1/8) ξ(1/2 + iz/2)                                          [S1 eq. (1)]
∂_τ H  = −∂_z² H                                                      [S1 eq. (31)]
```

Parameter domain: `τ ∈ ℝ`, with the real-zero property preserved **toward increasing `τ`**
and the de Bruijn–Newman constant defined by "all zeros real iff `τ ≥ Λ`". Provenance and
equation numbers in `SOURCE_AUDIT.md`.

**Coordinates, kept strictly separate throughout:** `s` the Dirichlet variable;
`z = −2i(s − ½)` the spectral coordinate, so for `s = β + iγ`, `z = 2γ − 2i(β − ½)` and
`Im z = −2δ`; `τ` the deformation parameter (**not** physical time, and no hidden
mechanism is claimed); `u = x` the kernel variable; `L` the detector half-window in the
test-function variable.

**Verified, not assumed** (`data/audit_log.txt` block A, B):

| check | result |
|---|---|
| `H_0(z) = (1/8)ξ(1/2+iz/2)` at `z = 0, 1, 5, 14.134725` | worst relative residual `7.17×10⁻⁴³` |
| `∂_τH = −∂_z²H`, finite differences, ζ kernel, `τ = 0` and `0.2` | `1.0×10⁻²³`, `2.6×10⁻²⁵` |
| the same for the DH kernel | `2.8×10⁻²²`, `5.2×10⁻²²` |
| `Φ` even (equivalent to the theta functional equation) | `3.9×10⁻⁷¹`, `1.6×10⁻⁶⁴`, `2.3×10⁻⁶⁰` at `x = 0.2, 0.5, 0.9` |

The analytic form of the PDE check is one line and worth recording because it fixes the
sign: `∂_τH = ∫x²e^{τx²}Φcos(zx)dx` and `∂_z²H = −∫x²e^{τx²}Φcos(zx)dx`, so
`∂_τH = −∂_z²H` identically — the backward heat equation in `z`.

*A precision trap, and how it was handled.* At 40 working digits the evenness check
appeared to **fail** at `x = 0.9`, returning a relative residual of `2.07×10⁺⁵`. That is
cancellation, not asymmetry: `Φ(0.9) ≈ 7.5×10⁻⁴⁶` while the series at `x = −0.9` has terms
of opposite sign and size `10⁻¹`, so the check needs working precision well above the
size of the cancellation. The check now sets the precision from each point's own
magnitude (`dps = 70, 70, 105`) and all three resolve, as tabulated. A point that still
did not resolve would be reported as `UNRESOLVED` rather than as a large residual — the
first version of this row recorded `2.07×10⁺⁵` as a "worst residual", which would have
been actively misleading in the CSV. Evenness is in any case *equivalent* to the
functional equation, so this is a bookkeeping test of the substitution, not the source of
the fact.

---

## 2. The Davenport–Heilbronn kernel, derived here

**Not inferred from the `ζ` case.** From the Mellin representation of the completed
function,

```
Λ_F(s) = ∫_0^∞ ω(t) t^{(s+1)/2} dt/t,      ω(t) = Σ_{n≥1} n a_n e^{−πn²t/5},
```

the substitutions `t = e^{4x}`, `s = ½ + iz/2` give `Λ_F(½+iz/2) = ∫_ℝ Φ_F(x)e^{izx}dx`
with

```
┌──────────────────────────────────────────────────────────────┐
│  Φ_F(x) = 4 e^{3x} Σ_{n≥1} n a_n exp(−π n² e^{4x}/5)         │
│  H^{DH}_τ(z) = ∫_0^∞ e^{τx²} Φ_F(x) cos(zx) dx               │
│  H^{DH}_0(z) = (1/2) Λ_F(½ + iz/2)                           │
└──────────────────────────────────────────────────────────────┘
```

Checks (block C):

| | result |
|---|---|
| `2∫_0^∞Φ_F cos(zx)dx = Λ_F(½+iz/2)` at `z = 0, 2, 10, 2γ₁` | `0, 0, 1.7×10⁻⁴¹, 5.1×10⁻²⁰` |
| **(10)** `Φ_F` even, at `x = 0.1, 0.3, 0.6, 1.0` | `2.1×10⁻⁴¹` … `3.4×10⁻²⁹` |
| **(11)** decay | `Φ_F(x) = O(e^{3\|x\|}e^{−(π/5)e^{4\|x\|}})`: double-exponential, so (11) holds for **every** `α > 0` |
| positivity on `x ≥ 0` | margin `0.3608086252 > 0` at `x = 0`; sampled minimum on `[0,3]` positive |

The positivity argument is elementary and worth stating because it is the same shape as
the one Newman–Wu give for `ζ`. **With the normalization corrected** (see the review
response below), with `T = πe^{4x}/5 ≥ π/5`, `|a_n| ≤ 1` and `a_1 = 1`, factor out the
leading term:

```
Σ_{n≥1} n a_n e^{−n²T} ≥ e^{−T} · [ 1 − Σ_{n≥2} n e^{−(n²−1)T} ].
```

The bracket is **increasing** in `T`, so one evaluation bounds it below for every
`x ≥ 0`: at `T = π/5` it is `0.6763199240515301 > 0`, rising to `0.8994`, `0.9950`,
`0.99999939` at `T = 1, 2, 5`. With evenness this gives `Φ_F > 0` on all of `ℝ`, with a
decaying prefactor `e^{−T}`.

The **unfactored** difference `e^{−T} − Σ_{n≥2} n e^{−n²T}` must not be used here: it
equals `0.3608086252` at `T = π/5` but tends to `0` as `T` grows (`2.06×10⁻⁹` at
`T = 20`), so it is not a uniform lower bound. The first version of this report quoted
that number as though it were one.

**Positivity is recorded but is not needed.** De Bruijn's Theorem 7 does not use it
(`SOURCE_AUDIT.md` §S3). It is reported because the surrounding Newman–Wu framework
assumes a positive measure, so its availability means the DH case sits inside that
framework too, not merely inside the single theorem.

### The resulting forward constraint

Every zero of `Λ_F` has `|β − ½| ≤ σ_0 − ½` with `σ_0 = 1.395136158235109721061359`
(**T1**, parent `REPORT.md` §3), and `Im z = −2(β − ½)`, so

```
Δ_DH = 2(σ_0 − ½) = 1.790272316470219442122718
τ_DH = Δ²/2       = 1.602537483559822777787819
```

and **all zeros of `H^{DH}_τ` are real for `τ ≥ τ_DH`.** For `ζ`, `Δ = 1` and
`τ ≥ 1/2` — de Bruijn's value, reproduced.

---

## 3. The bridge to the detector, and what it is worth

`Ψ_{ij}(z) = ½[G_i(iz/2)G_j(−iz/2) + G_j(iz/2)G_i(−iz/2)]`, and with `z = 2γ` this is
exactly the repository's symmetrized `H(ρ−½)`: verified against `dh.H_at` at
`γ = 3.5` and `γ = γ₁` to `1.1×10⁻⁴⁰` and `3.7×10⁻⁴¹` (block E). So
`(Q_τ)_{ij} = Σ_k m_k Ψ_{ij}(z_k)` is the right object, and at `τ = 0` it is the window
form the parent extension already assembled and gated.

**Lemma (reality ⟹ PSD).** For real `g` and **real** `z`, `G_j(−iz/2) = conj G_j(iz/2)`,
so with `v_j = G_j(iz/2)`,

```
Ψ(z) = Re(v v*) = (Re v)(Re v)^T + (Im v)(Im v)^T ⪰ 0.
```

Hence if every zero of `H_τ` is real, `Q_τ` is a sum of positive semidefinite matrices
and `Q_τ ⪰ 0`. Checked at `z = 1.0, 7.0, 171.4`: least eigenvalue at the precision floor
(`−2.3×10⁻⁴¹` or smaller in magnitude), `max|Im Ψ| = 0` exactly.

*Correction.* An earlier draft of this block asserted rank 2. The observed rank is **1**:
in a pure parity sector `G_j(iz/2)` is purely imaginary (odd) or purely real (even) for
real `z`, so one of the two rank-one terms vanishes identically — the same phase removal
as the channel decomposition of the preceding exploration. The conclusion `Ψ ⪰ 0` is
unaffected. The numerics caught this, not the prose.

**Combining §2 and the lemma:** `Q^{DH}_τ ⪰ 0` for `τ ≥ 1.6025374836`, and `Q^ζ_τ ⪰ 0`
for `τ ≥ 1/2`, both unconditionally, on **every** finite window and test space. That is a
real source-derived statement about the detector. It is also useless for the problem at
hand, because `τ = 0` is where the arithmetic lives and the implication runs the other way.

### Feasibility gate

| requirement | status |
|---|---|
| absolute convergence of `Σ_k m_k Ψ(z_k)` | **justified** — `Ψ = O(\|Re z\|^{−4})` for both repository bases (they vanish at `±L`), zero density `O(log)`, and `Ψ` grows only like `e^{L\|Im z\|}` inside a bounded strip |
| complete-orbit pairing and reality | **justified** — orbits `{z, −z, z̄, −z̄}`, `Ψ` even with `Ψ(z̄) = conj Ψ(z)`, so each orbit contributes a real symmetric matrix |
| agreement at `τ = 0` with the existing Weil form | **justified by reuse** — parent Phase C, `5.65×10⁻¹¹` vs the zero side and `1.73×10⁻¹³` vs an independent contour (**T2-N**) |
| differentiability in `τ` of the **complete** family | **NOT justified** |

The last line is the gate. A fixed zero-free contour gives an analytic *finite-contour*
family, and simple zeros move analytically; but differentiating the complete sum term by
term needs a uniform bound on `Σ_k |Ψ′(z_k) z_k′|` with `z_k′ = H″/H′(z_k)`, which is not
controlled near near-collisions, and a fixed contour can be crossed by a zero as `τ`
varies. **The missing lemma** is a uniform-in-`τ` tail bound for `Σ_{|z_k|>R}|Ψ′(z_k)z_k′|`
on a parameter interval, together with a contour that is zero-free throughout it.

Per the brief, the numerical `τ`-sweep branch is therefore **not run**, and no substitute
sweep is run in its place. The two decisive conclusions (§2 and §5) do not depend on the
gate: neither uses differentiability.

---

## 4. The first genuinely additional estimate attempted, and its outcome

### Attempt: a backward comparison inequality

Sought: `η′(θ) ≤ a(θ)(1 − η(θ))` with `θ` oriented to propagate `η ≤ 1` from `τ ≥ τ_DH`
(where §3 gives it unconditionally) down to `τ = 0`, with `a` controlled by the source
hypotheses.

**Refuted.** Chain, with tiers (block I):

| step | statement | tier |
|---|---|---|
| (i) | `Φ_F` satisfies (10), (11), and `Φ_F > 0` | T1 + numerical confirmation |
| (ii) | `Δ_DH = 2(σ_0 − ½)`, from the proved zero-free abscissa | **T1** |
| (iii) | de Bruijn Thm 7 ⟹ all zeros of `H^{DH}_τ` real for `τ ≥ 1.6025374836` | **T1** |
| (iv) | hence `Q^{DH}_τ ⪰ 0` there, by the lemma of §3 | **T1** |
| (v) | but `λ_min(Q^{DH}_0) = −0.679071346441480590986377586093…`, ball radius `1.39×10⁻¹¹⁶` | **T2-C** |
| (vi) | and `Q^{DH}_0` *is* that arithmetic matrix, by the explicit formula | **T1** (the derivation); the `5.65×10⁻¹¹` agreement is a cross-check of two implementations, and the Arb ball certifies the evaluated expression — neither is the source of the identity |

Any such inequality, applied to DH, would give `Q^{DH}_0 ⪰ 0` on the stated space,
contradicting (v)–(vi). **So no universally applicable backward-preservation argument
resting only on (10), (11), kernel positivity, a strip bound and the heat equation can
work** — somewhere along the path from `τ_DH` to `0`, one of its ingredients must fail.

The three roles are kept separate, as they must be: the **identity** `Q_0 = Σ_ρ Ψ(z_ρ)`
rests on its analytic derivation (T1); the `5.65×10⁻¹¹` agreement is a **cross-check** of
two independent implementations of the two sides; and the Arb ball **certifies the
evaluated expression**, not the identity.

**What the refutation leaves open.** It says nothing about a backward inequality using a
hypothesis DH *fails*. The only such hypothesis the parent work supplies is an **Euler
product**: DH has none, `ζ` does. So any surviving route must use the Euler product
essentially — and the parent Phase H already showed that deleting the non-Euler part of
the arithmetic side is not a usable way to do it (`λ_min` of the truncation is negative at
17 of 18 configurations, including at windows where the true form is positive).

### Why neither term of `η′` is controlled

`η_τ = u_τ^T B_τ^{-1} u_τ` gives `η′ = 2(u′)^T x − x^T B′ x` with `x = B^{-1}u`; at a
crossing `η = 1` one has `Qx = (B − uu^T)x = u(1 − η) = 0`, so `x` spans `ker Q`, and for
the eigenvalue through zero `λ_0′ = −η′/‖x‖²`. All three are verified in
`tests/test_identities.py` against a **synthetic, labelled** perturbation
(`3.97×10⁻³³`, `3.59×10⁻⁵¹`, `1.32×10⁻³⁰`) — because §3's gate did not justify a
source-derived derivative, no source-derived one is used.

The research question was whether the source law controls either term. **It does not:**
`u′` needs the selected quartet's velocity `H″/H′(z_1)`, and `B′` needs the velocity of
every other zero. Both require exactly the spectrum the detector was supposed to
constrain. A "prediction" that first solves the full spectral problem is not an
independent explanation, which is the brief's own disqualifier.

---

## 5. The adversarial controls

**Generic heat flow.** `F_τ(z) = z² + a² − 2τ` satisfies `∂_τF = −2 = −∂_z²F`, the same
equation, with zeros `±√(2τ − a²)`: purely imaginary for `τ < a²/2`, real for
`τ ≥ a²/2`. Its initial strip half-width is `Δ = a`, so it becomes real exactly at
`Δ²/2` — it **saturates** Theorem 7. Checked at `a = 1` (`τ_c = 0.5`) and
`a = Δ_DH = 1.790272316470219442122718` (`τ_c = 1.502537…` vs `1.702537…` straddling
`1.6025374836`).

*Hypothesis check.* `F_τ` is **not** the Fourier transform of an integrable, rapidly
decaying kernel — such a transform is **bounded on the real axis**, and a degree-2
polynomial is not. So `F_τ` is a **polynomial heat-flow control reproducing the
strip-shrinking rate**, not an example satisfying every hypothesis of Theorem 7, and not
a counterexample to it. The counterexample inside the intended class is the DH
construction of §2.

*What it establishes:* the entire transition from indefinite to positive semidefinite is
produced by the heat equation and the initial strip width alone. Any argument deriving a
zero restriction from those ingredients must also apply to `F_τ`, and therefore cannot be
arithmetic in content. The only arithmetic input anywhere in §2–§3 is the *value* of
`σ_0`, which is a crude elementary bound, not a deep fact.

**Manufactured dynamics.** For `B ≻ 0`, the flow `x′(t) = −½Bx(t)`, `x(0) = u` gives
`∫_0^∞‖x(t)‖²dt = u^T B^{-1} u = η`. Verified on the actual `L = 2.0` background:
`5.6022332216114890407` both ways, relative difference `6.85×10⁻²⁴`. So `η` *always*
admits a "dynamical" representation, for any positive definite matrix from any source.
The brief's warning is confirmed: this carries no information about where `B` came from,
and no such construction was used in the argument.

**Complete-object separation.** The three kinds of change were never combined: §2–§3
deform the *complete analytic function* at fixed detector space; the preceding
exploration's enrichment changed the *detector space* at fixed function; the parent
Phase F moved a *single orbit* at fixed everything else. Each is reported under its own
heading, and the `η′` identity is the only place where they would meet — which is exactly
where the gate stopped.

---

## 6. Numerical checks and their limitations

All numerics are at **40 working digits** (`data/audit_log.txt`, `data/validation.csv`),
with the identity tests at 50. That is deliberate: no quantity reported here needs more,
because none of them is close to a decision boundary — the two decisive numbers
(`τ_DH = 1.6025…` and `λ_min(Q_0^{DH}) = −0.679…`) are separated from zero by order one,
and the second is quoted from an existing Arb certificate rather than recomputed.

Limitations, stated separately as the brief requires:

* **Quadrature.** `H_τ` is evaluated by `mp.quad` on `[0, ½, 1, 3]`. The omitted tail
  beyond `x = 3` is bounded by the double-exponential decay (`Φ_F(3) ≈ 5.7×10⁻⁴⁴⁴⁰⁸`) and
  is irrelevant; but the quadrature error inside is a **library estimate**, not a bound.
* **Truncation.** The kernels are truncated at `N = 80` (ζ) and `N = 240` (DH) terms; the
  omitted terms carry `exp(−πN²e^{4x}/5)` and are far below the working precision.
* **Finite differences.** The PDE checks use a step `10⁻¹⁰` at 40 digits, so they verify
  the equation to about `10⁻²²`–`10⁻²⁵`, which is what is reported — not to the full
  working precision.
* **Zero localization.** No zero of `H_τ` is localized anywhere in this audit for `τ ≠ 0`,
  and none is needed. The `τ = 0` quartet locations are reused from the parent census and
  remain **T2, numerical inputs, not independently certified** — every statement here that
  mentions `γ₁` inherits that.
* **Linear solves.** Only in the manufactured-dynamics control and the identity tests;
  residuals are reported, and a small solve residual is not offered as a bound on anything
  else.
* **Precision stability is not interval certification.** Nothing in this directory is an
  interval enclosure. The one certified number quoted (`−0.679…`, radius `1.39×10⁻¹¹⁶`) is
  the parent extension's, reused.

No positivity is claimed from a failed search, no Cholesky failure is read as positivity,
no small eigenvalue is regularized away, and no complete spectral result is claimed from a
truncated contour.

---

## 6b. Two layers of positivity, and why one does not reach the other

Added in response to the external review, which supplied the observation.

`Φ_F ≥ 0` forces a positivity that has nothing to do with the zeros. For any nonnegative
integrable `Φ` with transform `Φ̂`,

```
Σ_{j,k} c_j conj(c_k) Φ̂(y_j − y_k) = ∫ Φ(x) |Σ_j c_j e^{i y_j x}|² dx ≥ 0,
```

so the matrix of Fourier correlations is positive semidefinite. For the DH kernel
`Φ̂_F(y) = Λ_F(½ + iy/2)`, so **the matrix `M_{jk} = Λ_F(½ + i(y_j − y_k)/2)` must be
PSD**. Checked at `y = 0, 1, 2.5, 4, 7` (block L): eigenvalues
`2.41×10⁻⁵, 2.98×10⁻³, 0.1417, 1.3638, 5.6575` — PSD, as required.

**And yet the same function's zero-based window form is certifiably negative.** There is
no contradiction: these are different quadratic forms built from different aspects of one
function, one from the kernel's correlations and one from the zeros. The lesson is sharp
and it applies to any layered programme:

> **Positivity in the source's correlation geometry does not transfer to positivity in
> its zero geometry.**

So the missing bridge is not "find something positive". It is "prove that *this*
positive structure controls *that* spectral form" — a strictly harder requirement, and
the one this audit did not meet.

---

## 6c. A two-sided bracket on the DH de Bruijn–Newman constant

Also added in response to the review, which derived the lower bound.

**Upper bound**, from §2: `Λ_DH ≤ 2(σ_0 − ½)² = 1.602537483559822777787819`.

**Lower bound**, by analytic persistence rather than computation. The inputs:

1. `H^{DH}_τ(z)` is jointly analytic in `(τ, z)` — the double-exponential decay of `Φ_F`
   makes the integral locally uniformly convergent;
2. at `τ = 0` there is a zero at `z_1 = 2γ_1 − 2iδ_1` with
   `Im z_1 = −0.617034364913274771106703921214` (a **numerical input**, T2, from the
   parent census), so a disc of radius `r < 0.617034…` about `z_1` does not meet the real
   axis;
3. the argument principle / Hurwitz on that disc gives a zero inside it for all
   sufficiently small `|τ|`, hence a **nonreal** zero, hence `τ < Λ_DH`;
4. forward preservation makes `Λ_DH` well defined.

Therefore

```
┌──────────────────────────────────────────────────┐
│   0 < Λ_DH ≤ 1.602537483559822777787819          │
└──────────────────────────────────────────────────┘
```

This is a **deduction** from the established kernel and the off-real-zero input, **not a
numerical localization**. Quantifying "sufficiently small" needs exactly the locally
uniform estimate that §3's gate records as missing, so **no numerical lower bound is
claimed and no flow numerics were run for it**. It also does not require differentiating
the complete detector family — which is the point the review makes about the endpoint
argument generally, and it is right: the counterexample of §4 needs only the two
endpoints, not the differential law.

The upper bound is a **sufficient** condition from the strip theorem, not a computed
location of the actual transition.

---

## 6c′. The lower bound, certified

`src/lambda_lower_bound.py` → `data/lambda_bracket.csv`, `data/lambda_bracket.txt`.
Arb at 400 bits. This replaces §6c's non-quantitative "`Λ_DH > 0`" with a number, and it
follows the review's §3 argument exactly; it touches neither the global differentiability
material of §3 nor the flow numerics.

### The disc

In the audit's variable the critical line **is** the real axis, and the first off-line
zero of `H_0` sits at `z_1 = 2γ_1 − 2iδ_1`:

| | |
|---|---|
| `Re z_1` | `171.398696970755184343858535418` |
| `Im z_1` | `−0.617034364913274771106703921214` |
| radius `r` | `0.3` |
| margin `\|Im z_1\| − r` | `0.31703436491327477111 > 0` |

so the closure of `D` is disjoint from the real axis.

### `m = min_{∂D} |H_0|`

`∂D` is partitioned into **3000** arcs of width `2π/3000`. Each arc
`{c + re^{iθ} : |θ − θ_m| ≤ Δθ/2}` is contained in the ball of centre `c + re^{iθ_m}` and
radius `2r sin(Δθ/4) = 3.14159251004×10⁻⁴`, because
`|re^{iθ} − re^{iθ_m}| = 2r|sin((θ−θ_m)/2)|`. `H_0 = Λ_F(½+iz/2)/2` is enclosed on each
such ball with Arb's `Γ` and Hurwitz `ζ`. Worst relative enclosure radius over the arcs:
`0.0649`.

```
m enclosure = [6.5935981835929959400×10⁻³⁰ ± 4.52×10⁻⁵⁰]
m ≥ 6.5935981835929959399×10⁻³⁰  > 0
```

so `∂D` passes through no zero. The second line is a **plain decimal floor**: the ball's
lower endpoint rounded *down*, so the inequality is exact. The first version of this
section quoted `…9400×10⁻³⁰`, which is the midpoint rounded to nearest and therefore
asserts `4.5×10⁻⁵⁰` more than the enclosure supports. See the corrections list. The partition was sized from a measurement, not guessed:
the Hurwitz-zeta sum amplifies an input ball radius by a factor `≈ 106` at this height
(that is the genuine `|ζ′|`, not an Arb artefact), so a 60-arc partition gave a relative
enclosure radius of `3.30` and failed; 3000 arcs give `0.065`.

### The zero inside `D` is certified, not assumed

Summing the principal argument increments of the arc enclosures:

```
(1/2π) Σ Δarg = 1.00000000000000 ± 3×10⁻²⁰
worst (|increment| + both arc uncertainties) = 0.13136  <  π
```

Every increment is unambiguous, so the winding number is **exactly 1**: `H_0` has exactly
one zero in `D`, counted with multiplicity. **The census is therefore not an input to this
bound** — the disc centre came from it, but nothing depends on that centre being a zero.

### `M` and `C`

`M = max_{∂D}|Im z| = |Im z_1| + r = 0.91703436491327477111`, `τ_0 = 0.1`.

The perturbation bound rests on one inequality, verified in one line:
`|e^{τx²} − 1| = |∫_0^{τx²} e^t dt| ≤ |τ|x² e^{|τ|x²} ≤ |τ|x² e^{τ_0x²}` for `|τ| ≤ τ_0`.
With `|cos(zx)| ≤ cosh(x·Im z) ≤ e^{Mx}` for `x ≥ 0`,

```
|H_τ(z) − H_0(z)| ≤ |τ| ∫_0^∞ x²|Φ_F(x)| e^{τ_0x² + Mx} dx =: |τ| C   on ∂D.
```

`C` is split at `X = 2`. On `[0, X]`, rigorous Arb integration of the explicit integrand,
with the kernel truncated at `n ≤ 40` and the omitted terms carried as an error ball
(`Σ_{n>N} n e^{−πn²/5} ≤ 2(N+1)e^{−π(N+1)²/5}`, using `e^{4x} ≥ 1` on `x ≥ 0`):

```
∫_0^2 = 0.048291904880526492445
```

For the tail, `|Φ_F(x)| ≤ 8e^{3x}e^{−πe^{4x}/5}` for `x ≥ 1`, and with
`f(x) = x²·8e^{3x}e^{−πe^{4x}/5}e^{τ_0x²+Mx}`,

```
(log f)′(x) = 2/x + 3 − (4π/5)e^{4x} + 2τ_0x + M  ≤  −1   for x ≥ 2,
```

since `(4π/5)e^{8} = 7491.9646` dominates `2/x + 3 + 2τ_0x + M + 1` there and grows far
faster. Hence `f(x) ≤ f(2)e^{−(x−2)}` and `∫_2^∞ f ≤ f(2) = 4.4817×10⁻⁸⁰⁹`. So

```
C enclosure = [0.048291904880526492445 ± 3.83×10⁻²²]
C ≤ 0.048291904880526492446
```

again a plain decimal, here rounded *up*.

### The result

```
┌────────────────────────────────────────────────────────────────────────────┐
│  ε = min(τ_0, m/(2C)) ≥ 6.8268151773112564447×10⁻²⁹                        │
│                                                                            │
│  6.8268151773112564447×10⁻²⁹  ≤  Λ_DH  ≤  1.602537483559822777787819       │
└────────────────────────────────────────────────────────────────────────────┘
```

`ε` is unchanged by the rounding correction: it is computed from the ball's lower
endpoint for `m` in the first place (`arb(mlow.mid() - mlow.rad())`), so only the
*displayed* digits of `m` and `C` moved, each by one unit in the last place.

**Tiers.** `ε` is **T2-C**: `m`, `C` and the winding number are Arb enclosures, and the
identity `H_0 = Λ_F(½+iz/2)/2` is the T1 Mellin derivation of §2. The upper bound is
**T1 given** the strip bound `σ_0 = 1.395136158235109721061359`, which is itself T1 in the
parent extension.

### Why `ε` is so small, and what that says

`ε ≈ 7×10⁻²⁹` is a true certified bound and a quantitatively weak one. The reason is
structural, not a tuning artefact:

* `m ≈ 6.6×10⁻³⁰` is small because `Λ_F` carries the gamma factor, whose modulus at
  `Im((s+1)/2) ≈ 42.85` is of order `e^{−πγ₁/4} ≈ 6×10⁻³⁰`. That is a **prefactor**, not
  cancellation, and it is present for any disc about this zero.
* `C ≈ 0.048` is `O(1)` because the perturbation bound replaces `cos(zx)` by `e^{Mx}` and
  so **discards every bit of oscillatory cancellation** — the same cancellation that makes
  `H_0` itself small.

So the ratio `m/C` inherits the gamma decay on one side and none of it on the other, and
`ε ∼ 10⁻²⁹` for any reasonable disc radius; no choice of `r` changes the order of
magnitude. Improving it would need a perturbation estimate that keeps the oscillation —
a steepest-descent or saddle-point treatment of `∫(e^{τx²}−1)Φ_F(x)cos(zx)dx` — which is
a different and much more delicate undertaking. It was **not** attempted, and the disc and
the bound were **not** loosened to produce a larger number.

**Reproduced.** The whole certification was re-run from scratch on 2026-10-24
(`--arcs 3000 --prec 400`) and its output is **identical** to the committed
`data/lambda_bracket.txt` in every digit, including the Arb radii. A cheap regression
fixture, `tests/test_lambda_certificate.py`, now guards it: it compares every quantity
against a value written into the fixture itself — so editing the data file fails the test
rather than redefining what passing means — and **recomputes** `ε = min(τ_0, m/(2C))` in
exact rational arithmetic instead of reading it, together with `|Im z_1| − r` and
`|Im z_1| + r`. Checked against a negative control: changing one digit of `m` fails it.

What the bound does settle is the qualitative question, with a certificate rather than a
continuity appeal: `Λ_DH` is **bounded away from zero**, so the Davenport–Heilbronn flow
genuinely has a positive transition parameter, and the two-sided bracket of §6c is now
numerical at both ends.

---

## 7. The principal obstruction

**Direction, and the absence of arithmetic away from `τ = 0`.**

The source deformation supplies exactly one unconditional implication —
*bounded initial strip + heat flow ⟹ real zeros, hence `Q_τ ⪰ 0`, for `τ ≥ Δ²/2`* — and
it points away from `τ = 0`. Reversing it is not merely unproved: for the one object in
this repository where the answer at `τ = 0` is known with certainty, the reversal is
**false**. And for `τ > 0`, **no appropriate prime-side representation of the deformed
family has been established** — not "there is none", which would be an impossibility
theorem this audit does not prove. So there is no *identified* arithmetic constraint there
to carry back, and at `τ = 0` the flow gives a derivative whose two terms both require the
full spectrum.

The obstruction in one sentence — and this replaces a sharper-sounding slogan that an
external review correctly judged too strong:

> **The general deformation theorem controls zero locations without using the arithmetic
> properties that distinguish ζ from Davenport–Heilbronn. No estimate transferring those
> distinguishing properties into backward positivity has been established.**

The earlier formulation ("the deformation is informative exactly where the arithmetic is
absent, and the arithmetic is present exactly where the deformation is uninformative") is
**withdrawn**: `Φ_{F,τ}` still contains the `a_n` at every `τ`, so the arithmetic does not
disappear. The problem is that the theorem being applied does not use enough of it.

---

## 8. One next lemma — and whether the evidence warrants it

**The evidence does not warrant expanding the flow numerics**, and I am not proposing a
next experiment in this direction. The forward statement is classical, the backward
statement is blocked for the stated class, and no usable arithmetic structure has been
identified at `τ ≠ 0`.

What this audit does **not** establish is that intrinsic arithmetic dynamics cannot help.
It establishes that *this* deformation, controlled by *these* general properties, operates
in a class too broad to enforce the desired initial condition. Three things have to be
kept apart:

```
a valid evolution  ≠  a useful spectral restriction  ≠  an arithmetic explanation
                                                        of the initial spectrum
```

This audit obtained the first, and a one-directional version of the second. It did not
obtain the third.

If anything from this audit is worth keeping, it is not a lemma about the flow but the
**obstruction instrument** it produced:

> **The filter.** The Davenport–Heilbronn function satisfies reflection symmetry, rapid
> decay, **kernel positivity**, and a bounded zero strip, with an explicit
> `Δ_DH = 1.790272316470219442122718`, while its window form at `τ = 0` is certifiably
> indefinite. So **any proposed enforcement theorem must contain an independently
> established hypothesis that excludes this control.** Invoking kernel positivity,
> symmetry, fast decay, a bounded strip, or generic heat evolution will not do it.

That is the one thing here that will save work later, and it is cheap to apply.

**The methodological direction it points to** is the one already used in this author's
`hypercuboid-character-gateways` project, where the three finite-field character-sum
symmetry classes are first **reduced** — Classes I and II collapse to a classical
two-variable character sum by an elementary projective argument — and only the genuine
residue, Class III, is then handed to the deeper geometric identification. The analogue
here is: **reduce the admissible class of arithmetic objects first, using genuine
arithmetic constraints, and apply the deformation theorem only to what survives.** This
audit did the opposite — it applied a general theorem to an unreduced class — and the DH
control is exactly what an unreduced class lets in.

Sharpening the bracket of §6c (deciding whether `Λ_DH < τ_DH` strictly, or locating it)
would need the off-line quartets tracked under the flow. **It is not attempted**, and per
the review's assessment it should not be until a new, independently justified arithmetic
constraint has been identified first.

---

## 9. Scope statements

* **Nothing here bears on RH or GRH**, and no finite calculation here establishes an
  infinite-dimensional limit. The forward statements `Q^ζ_τ ⪰ 0` for `τ ≥ 1/2` and
  `Q^{DH}_τ ⪰ 0` for `τ ≥ 1.6025…` hold on every finite window, but say nothing at `τ = 0`.
* **"Dynamics" means a specified parameter evolution.** No physical time, and no hidden
  physical mechanism generating zeta, is claimed or implied anywhere.
* **The direction is recorded everywhere it matters.** Real-zero preservation is toward
  increasing `τ`; no statement here infers anything backward toward `τ = 0`.
* **No best-bound claim.** The Polymath numerical upper bound on `Λ` is not used; `Λ ≥ 0`
  is cited only as a remark.
* **A positive finite restriction is not whole-window positivity**, and removing one
  orbit does not leave only on-line zeros — the other off-line quartets remain implicitly
  in the complete form.
* **No novelty is claimed.** De Bruijn's Theorem 7, the heat equation, the zero-velocity
  identity, the derivative of `η` and the rank-one algebra are all textbook or cited. The
  only assembly I have not seen elsewhere is the application of Theorem 7 to the
  Davenport–Heilbronn kernel, and I make no priority claim about it.
* **No manuscript concern arose** during this audit; had one, it would be recorded here
  and nowhere else.

### Corrections made during this audit

1. **"rank 2" → rank 1** in the PSD lemma block (§3). In a pure parity sector one of the
   two rank-one terms vanishes identically. The numerics caught the prose error; the
   lemma's conclusion is unaffected.
2. **`mp.eigsy` on a numerically complex matrix.** `Ψ(z)` carries an imaginary residue at
   the precision floor for real `z`; the eigenvalue check now takes the real part and
   reports `max|Im Ψ|` separately, which is itself a check of the lemma.
3. **A cancellation artefact recorded as a residual.** The `Φ_zeta` evenness row first
   reported `2.07×10⁺⁵` as its "worst relative residual", which reads in the CSV as a
   failed check. It was a 40-digit evaluation of a quantity needing about 105 digits. The
   check now chooses its precision per point and reports `UNRESOLVED` where it cannot
   decide; all three points resolve, at `3.9×10⁻⁷¹`, `1.6×10⁻⁶⁴`, `2.3×10⁻⁶⁰`.
4. **A one-sided bound printed with round-to-nearest (§6c′).** `arb.str(n)` rounds the
   midpoint, so `m ≥ 6.5935981835929959400×10⁻³⁰` asserted `4.5×10⁻⁵⁰` more than the
   enclosure supports — numerically harmless, wrong in principle for a *certified* floor.
   `src/lambda_lower_bound.py` now prints one-sided bounds with `floor_str` / `ceil_str`,
   which round the ball's own endpoint down (for a floor) or up (for a ceiling) and emit
   a plain decimal with no `±` part. `m ≥ …9399×10⁻³⁰`, `C ≤ …446`. **`ε` is unaffected**,
   because it was already computed from `mid − rad`. Found while preparing the note, by
   asking whether each published digit of a one-sided bound was actually implied.

5. **A patch that did not land.** The rank-1 correction above was first applied by a
   backgrounded command that completed *after* the run it was meant to fix, so the
   regenerated transcript still carried the wrong prose. Caught by grepping the new log
   for the stale string instead of assuming the edit had taken effect. The same thing
   happened a second time with this very list.

---

## Review response

An external review checked the cited theorem's hypotheses and parameter convention and
confirmed both, while declining to rerun the scripts or re-audit the kernel construction
step by step. It raised six points. **All six are accepted**, and four of them changed the
mathematics or its presentation rather than only the wording:

| # | the review's point | what changed |
|---|---|---|
| 1 | The uniform positivity bound was presented wrongly: `e^{−T} − Σ_{n≥2} n e^{−n²T}` tends to `0`, so `0.3608` is not a uniform lower bound. The factored form `e^{−T}[1 − Σ_{n≥2} n e^{−(n²−1)T}]` has an **increasing** bracket, `≈ 0.67632` at `T = π/5`. | **§2 and block C4 rewritten.** Verified: bracket `= 0.6763199240515301` at `T = π/5`, rising to `0.99999939` at `T = 5`; the old expression falls to `2.06×10⁻⁹` at `T = 20`. The review's correction is right and my presentation was wrong. |
| 2 | "No such inequality exists" overstates it; the correct statement is that no *universally applicable* backward argument on those shared hypotheses can work, and that somewhere on the path the inequality, the coefficient control, the differentiability, or `B_τ ≻ 0` must fail. | **§4 and block I rewritten** to that statement, including the two reasons it matters (`B_τ ≻ 0` is a standing assumption; zero-tracking can meet collisions). |
| 3 | The slogan "the deformation is informative exactly where the arithmetic is absent…" is too strong: `Φ_{F,τ}` still contains the `a_n`. | **Slogan withdrawn** and replaced; see §7. The replacement names the actual problem — the theorem does not use enough of the arithmetic. |
| 4 | "There is no prime side at all" should be "no appropriate prime-side representation has been established". | **Changed** in §1, §7 and block K. The audit proves no impossibility theorem about every conceivable representation. |
| 5 | The toy polynomial is not a Fourier transform of an integrable rapidly decaying kernel, because such a transform is **bounded on the real axis** and a quadratic is not. Describe it as a polynomial heat-flow control reproducing the strip-shrinking rate. | **§5 rewritten** with the boundedness argument, which is a better reason than the one I gave. |
| 6 | A two-sided bracket is available: nonreal zeros persist for small `τ` by analytic continuity, giving `0 < Λ_DH ≤ 1.6025374836`, without differentiating the complete detector. Also: positive kernels force a PSD Fourier-correlation matrix, which coexists with a negative zero-based witness. | **Two new sections, §6b and §6c**, with the PSD matrix verified (least eigenvalue `2.41×10⁻⁵`). Both are the review's observations, recorded as such. |

The review also noted that the identity, the numerical cross-check and the Arb enclosure
play three different roles. That is now stated explicitly in §4 rather than compressed
into a single tier label.

---

## Closing answer

> **Have we transferred an independently established source constraint into the detector,
> or have we only described how an unconstrained detector evolves?**

**We transferred one, and it is the wrong one — but the transfer itself worked.** De
Bruijn's theorem, applied through a kernel derived and checked here, gives `Q_τ ⪰ 0` for
`τ ≥ 1/2` (ζ) and `τ ≥ 1.602537483559822777787819` (Davenport–Heilbronn), on every finite
window, with no hypothesis beyond reflection symmetry, rapid decay and a bounded initial
strip. That is an independently established source constraint landing on the detector, and
the bridge did not fail at its first step.

What it does not do is distinguish the two functions. The theorem that produces it never
uses the arithmetic that separates `ζ` from Davenport–Heilbronn, and the DH control —
which has a **positive** kernel, reflection symmetry, rapid decay and a bounded strip —
is certifiably negative at `τ = 0`. So the general properties controlling this dynamics
are too broad to enforce the condition we want, and the deformation can manufacture a
property the starting function did not have; finding it at a later parameter is therefore
no evidence it was there initially.

Beyond that one forward implication, what remains is a description of how an
unconstrained detector evolves. The useful residue is the filter of §8 and the
methodological lesson that the class must be reduced by genuine arithmetic constraints
*before* a general deformation theorem is applied, not after.

**Stopping here.** The sharpening of §6c is not pursued, and no flow numerics are
expanded.

---

## 10. The note

A standalone 12-page note, [`note/lambda_note.pdf`](note/lambda_note.pdf),
[doi:10.5281/zenodo.23628224](https://doi.org/10.5281/zenodo.23628224), states the two results
worth carrying forward and nothing else:

> **A two-sided certified bracket**
> `6.8268151773112564447×10⁻²⁹ ≤ Λ_DH ≤ 1.602537483559822777787819`,
> **and the obstruction**: de Bruijn's theorem needs only conjugate reflection, rapid decay
> and a bounded initial strip, all of which Davenport–Heilbronn has — together with a
> positive kernel — while its window form at `τ = 0` is certifiably indefinite. So no
> argument resting on that package alone can carry positivity backward to `τ = 0`.

`make -C note verify` re-derives every number printed in it from this directory's data
files and fails on any disagreement, any undeclared literal, any overstatement of precision
relative to the source's own significant digits, a bibliography number not found in the
record, or an unmet prerequisite — pass 2 refuses to pass unless `data/lambda_bracket.txt`
records the disc, `m > 0` and the unambiguous winding increments, **and**
`tests/test_lambda_certificate.py` passes. 55 literals matched a derived quantity, 196 are
structural or quoted from this record.

**The normalisation audit (§S4 of `SOURCE_AUDIT.md`) is part of this note's claim, not a
footnote to it.** The threshold `Δ²/2` is a statement about the variable the transform is
taken in, and the literature carries at least three such variables for `ζ`, so one theorem
appears in print as `1/8` (Newman 1976, `z = −i(s−½)`), as `1/2` (Polymath15 and
Rodgers–Tao, `z = −2i(s−½)`), and in de Bruijn's own `Ξ(2z)` scaling. **This audit is in the
Polymath15 convention**, pinned by its own check that `ζ` returns `1/2` rather than by
inheritance, and the note states the convention explicitly. The DH number is `2Δ_s²`; in the
Titchmarsh-type scaling it would read `Δ_s²/2 = 0.40063437088995569`, exactly a quarter.

What the note deliberately leaves out: the feasibility gate of §3, the `η′` computation of
§4, the manufactured-dynamics control of §5, and the Phase H experiment, which appears only
as a single citation inside the note's §6 where it is needed to show why the Euler product
does not immediately supply the missing hypothesis.

Two things the note does **not** say, both of which an earlier draft of this record did:

* it makes no claim about the literature. The audit surveyed none in this direction, so the
  absence of a distinguishing estimate is recorded as an absence *in the note*;
* it does not say the heat flow cannot prove RH. It says that one named package of
  hypotheses cannot, because Davenport–Heilbronn satisfies all of them.
