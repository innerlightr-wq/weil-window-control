# GENERATOR_AUDIT.md

The actual generators, and whether an exact connection to off-line displacement exists.
Tags: **[SOURCE]** (in a checked source), **[STD]** (standard, verified here), **[NEW]** (derived
here), **[NUM]** (numerical), **[NO-GO]** (proved negative). `γ` = spectral ordinate throughout.

> **AMENDED 2026-10-07 (later pass, after review).** Three corrections, all adopted; the verdict
> **MISMATCH** for the specific identification is unchanged. **(i)** The wording "`η` breaks the
> generator identity" is **withdrawn** — the identity is a fact about `δ` and is untouched by
> relabelling; what fails is that the relabelled family `V_η = e^{−½tanh(η)A}` has **no constant
> generator in `η`** (§4). **(ii)** A **domain qualification** is added: `[A,D] = −I` holds on a
> *common core*, and `D = d/du` is **not** skew-adjoint on `L²[−L,L]` without boundary conditions —
> zero-extended window functions are not a translation-invariant space (§5.1). **(iii)** A
> relationship between two operators is **not** an identification of them. Verified in
> `amend_checks.py` / `amend_output.txt`. A positive reading is added as §3.1.

---

## 1. Transform convention, fixed from the source

**[SOURCE: `notes/proofs.md` §§B.1–B.3]** `f` real, even, `supp f ⊂ [−L,L]`, `‖f‖₂ = 1`, and

```
F(t) = int f(u) e^{i t u} du ,
F'(t)  = - int u   f(u) sin(t u) du ,
F''(t) = - int u^2 f(u) cos(t u) du        (real even f) .
```

So **`F'(γ)` corresponds to multiplication by `iu`** under this convention, and `F''(γ)` to
multiplication by `−u²`. The brief's question "`iu`, `−iu`, or another" is answered: with
`F(t) = ∫f e^{itu}du`, differentiation in `t` is multiplication by `iu`. (With the opposite sign
convention `F(t) = ∫fe^{-itu}du` it is `−iu`; the repository uses the first, and the real-even
reduction above is the check that fixes it.)

## 2. The two generators, named

| | **intrinsic scaling** (Sonin/CC2021) | **off-line displacement** |
|---|---|---|
| action | `ϑ(λ)ξ(v) = λ^{-1/2}ξ(λ^{-1}v)` on `L²(R)_ev` **[SOURCE]** | `β = ½+δ` ⟹ ordinate `γ ↦ γ + iδ` |
| in `u = log v` | **translation** `u ↦ u + t`, `λ = e^t` | **multiplication** by `e^{−δu}` |
| generator | `D = d/du` | `A = ` multiplication by `u` |
| adjointness | **skew**-adjoint; `e^{tD}` unitary | **self**-adjoint; `e^{−δA}` positive, **not** unitary |
| boundedness | `‖e^{tD}‖ = 1` **on `L²(R)`** | `A` unbounded on `L²(R)`; **bounded on `L²[−L,L]` with `‖A‖ = L`** |
| what `t`/`δ` acts on | the physical/support coordinate | the spectral ordinate, imaginarily |
| domain | `L²(R)`; **not** the fixed window — see §5.1 | `A` bounded ⟸ the window; this is where `L` enters |

**The two operators do not live on the same space**, and §5.1 records exactly what that costs. The
table compares them; it does not place them on a common domain.

`A` is the generator the brief asks for: **already present**, as the multiplication operator dual
to the scaling translation. It was not invented to make a formula work — §3 shows it is forced.

## 3. The exact generator identity  **[NEW statement, [STD] content]**

**Proposition 1.** For `f` supported in `[−L,L]` and `A` = multiplication by `u` on `L²[−L,L]`,

```
F(gamma + i delta)  =  \widehat{ e^{-delta A} f } (gamma)          for every real delta,
```

i.e. the imaginary shift of the spectral ordinate is **exactly** the one-parameter semigroup
`e^{−δA}` applied to the test function.

*Proof.* `F(γ+iδ) = ∫f(u)e^{i(γ+iδ)u}du = ∫(e^{−δu}f(u))e^{iγu}du`. ∎

That is literally `notes/proofs.md` §B.3's own `G(δ) = ∫f(u)e^{iγu}e^{−δu}du`, read as an operator
statement. **[NUM]** verified for three test functions at `δ = 0.05, 0.2, 0.375`, `L = 0.8`,
`γ = 14.1347`: deviations `0` to `3.7e-17` (`checks.py` §3).

**Corollary 2 (the exact quartet functional).** With `Φ(δ) := Q_δ(f) − Q_0(f)`,
**[SOURCE: `notes/proofs.md` §B.1]** gives, with no remainder,

```
Phi(delta) = 4 Re F(gamma + i delta)^2 - 4 F(gamma)^2
           = 2 [ g(delta)^2 + g(-delta)^2 ] - 4 g(0)^2 ,     g(delta) := \widehat{e^{-delta A} f}(gamma) ,
```

the second form because `conj G(δ) = G(−δ)` for real even `f`. `Φ` is therefore **even** in `δ`
(**[NUM]** `|Re G(+δ)² − Re G(−δ)²| ≤ 2.6e-17`), so `Φ(δ) = Ψ(δ²)` with `Ψ` **entire** — `\hat f`
is entire of exponential type `L` because `f` is compactly supported.

**Note the shape.** `A` acts on each factor **before** squaring. So `Φ` is a sum of squares of
`e^{∓δA}`-displaced evaluations, **not** `⟨(e^{δA}+e^{−δA})f, ·⟩`. The brief's probe
`e^{ηA}+e^{−ηA}` is the right *kind* of object but the wrong *placement* and, as §4 shows, the
wrong parameter.

## 3.1 What this actually says, in plain terms  **[the useful outcome]**

Proposition 1 has a reading worth stating, because it is the substantive content of this audit:

> **The window does not move. The contribution from different positions inside it is reweighted.**

Moving the evaluation point off the critical line multiplies the test function by `e^{−δu}`: for
`δ > 0` one side of the window is amplified and the other attenuated, while every position stays
where it is. This is the classical transform relationship between spectral translation and
exponential weighting — **not** a new arithmetic mechanism, and not a dynamics.

It also explains the second moment without any reference to rapidity. Expanding
`e^{−δu} = 1 − δu + δ²u²/2 − ⋯`, the count-preserving quartet comparison cancels the odd powers of
`δ` (Corollary 2: `Φ` is even), and what survives at second order is built from the **first and
second coordinate moments** of the window — which is exactly what `F'(γ)` and `F''(γ)` are
(§1: `F'` ↔ multiplication by `iu`, `F''` ↔ by `−u²`). So:

> **Coordinate moments enter the sensitivity law because exponential reweighting differentiates
> into powers of the coordinate.**

That is a simpler and better-founded explanation than the proposed bridge, and it is the thing to
keep. It does **not** make the response a nonnegative norm: `K(L,γ)` controls the
derivative-evaluation piece, while `F'² + FF''` carries either sign (§7).

## 4. THE FALSIFICATION GATE: `η` is additive for a law the problem does not have  **[NO-GO]**

**The native composition law is plain addition in `δ`:**

```
e^{-delta_1 A} e^{-delta_2 A} = e^{-(delta_1 + delta_2) A} .
```

**[NUM]** `checks.py` §4: composing two displacements and comparing against the single
displacement `δ₁+δ₂` gives deviations `8e-18` to `3.8e-17`; comparing against the Möbius
composite `(δ₁+δ₂)/(1+4δ₁δ₂)` gives `2.1e-04`, `1.5e-03`, `4.5e-03` — **four to fourteen orders of
magnitude larger**, at `(δ₁,δ₂) = (0.1,0.15), (0.2,0.3), (0.375,0.375)`.

Now `2δ = tanh η`, so `η` is additive **exactly when** `δ` composes by Möbius addition
`δ₁ ⊕ δ₂ = (δ₁+δ₂)/(1+4δ₁δ₂)`. It does not. Therefore:

> **No-go 1.** The only native composition on the off-line displacement is `e^{−δA}`'s semigroup
> law, which is **additive in `δ`**. The rapidity coordinate `η = arctanh(2δ)` is additive for
> Möbius addition, which is **not** the native law. So `η` is additive for the wrong operation,
> and the additivity that makes rapidity attractive is **spurious here**.

**Stated precisely, and corrected 2026-10-07.** `η` does **not** "break" Proposition 1 — that
identity is about `δ` and is untouched by relabelling. What fails is weaker and sharper. Put
`V_η := e^{−½tanh(η)A}`, the same operator family written in rapidity. Then

```
d V_eta / d eta  =  - (1/2) sech^2(eta) * A * V_eta ,
```

so the generator **depends on `η`**: the family is not generated by any fixed operator with respect
to additive `η`. **[NUM]** `amend_checks.py` §A: the coefficient matches `−½sech²η` to `1.7e-11` at
`η = 0, 0.3, 0.6584789, 0.9729551`, taking the values `−0.5, −0.4575685, −0.3333334, −0.2187500`.
Consequently `V_{η₁}V_{η₂} = V_{η₁+η₂}` **iff** `tanh η₁ + tanh η₂ = tanh(η₁+η₂)`, which fails:
differences `2.657e-02` at `(0.2,0.3)`, `1.626e-01` at `(0.5,0.5)`, `5.400e-01` at
`(0.9729551, 0.9729551)` — the last being the quasi-RH endpoint `|δ| = 3/8`. It holds only in the
trivial case where one argument is `0`.

> **So: you can relabel the displacement using rapidity, but that does not make the operation obey
> rapidity addition.** The coordinate is exact; the group law is not inherited.

**Searched, per the brief's §8, and not found:** multiplicative scaling composes as
`ϑ(λ₁)ϑ(λ₂) = ϑ(λ₁λ₂)`, i.e. additively in `t = log λ`; Mellin translations compose additively;
the `θ_S` local-factor transfer composes **multiplicatively in the factors** and does not act on
`δ` at all (`extensions/sonin-first-prime-defect-2026-10-07/`); no Möbius action on zero
displacement appears in any checked formula. **This is "not found in the sources checked", not a
proof of nonexistence** — but the positive fact above (the native law *is* additive in `δ`) is a
proof that Möbius addition is not it.

**Domain mismatch, in the manuscript's own vocabulary.** **[NUM]** `checks.py` §4: additive
composition **leaves** the partition domain (`0.3+0.3 = 0.6 > ½`; `0.375+0.375 = 0.75`;
`0.45+0.4 = 0.85`) while Möbius composition **cannot** (`0.4412`, `0.4800`, `0.4942`, all `< ½`).
By the manuscript's §6.4 criterion — *"If domain, bounds, or singularities do not align across
translation, the correct diagnosis is likely … a true mismatch"* — and its §6.6 **"Broken
domains"** failure mode, this is a **true mismatch**. The bounded coordinate is bounded for a
reason that has no counterpart on the operator side.

## 5. The two generators are canonically conjugate, not equal  **[NO-GO]**

**[STD]**, exact over `Q` on polynomials (`checks.py` §5): `(AD − DA)p = up' − (up)' = −p`, so

```
[ A , D ] = - I .
```

Verified for `p = 1`, `p = u`, and an arbitrary quartic; the sign-convention control
`[D,A] = +I` also checked. `A` and `D` are a **Heisenberg pair**: non-commuting, and exactly dual
under Fourier/Mellin — multiplication by `e^{−δu}` is the Fourier-conjugate of translation
`u ↦ u+t`.

> **No-go 2.** The Sonin scaling parameter `t` and the displacement parameter `δ` are **not the
> same kind of variable**. `t` generates a *real* translation of the physical coordinate; `δ`
> generates an *imaginary* translation of the spectral ordinate. They act along conjugate,
> non-commuting directions. **`t = η` is unjustified, and so is `t = δ`** — they are dual, not
> equal. The brief's fallback question ("is one dual to the other under Mellin/Fourier?") is
> answered **yes, exactly**, and duality is not identity.

**A relationship is not an identification.** `[A,D] = −I` is a genuine and exact relationship, and
it is all it is. It supplies **neither** the missing semi-local Sonin–Weil trace comparison (see
`extensions/sonin-first-prime-defect-2026-10-07/`) **nor** any new positivity bound.

### 5.1 Domain qualification — added 2026-10-07

The commutator was verified **on polynomials** (`checks.py` §5). Polynomials restricted to
`(−L,L)`, and `C_c^∞(−L,L)`, are a **common core** for `A` and `D`; that is the correct scope of
the computation, and the following two things are **not** true on the fixed window:

- **`e^{tD}` is not a unitary group on `L²[−L,L]`.** Translation moves support out of the window:
  **[NUM]** for `supp f = [−0.8, 0.8]`, the shifted support `[−0.8+t, 0.8+t]` overlaps the window
  in length `1.5, 1.2, 0.8, 0.0` at `t = 0.1, 0.4, 0.8, 1.6`, so the retained mass fraction is at
  most `0.938, 0.750, 0.500, 0.000`. **Zero-extended window functions are not a
  translation-invariant space**, so there is no translation group *on that space* to be generated.
- **`D` with no boundary condition is not skew-adjoint on `L²[−L,L]`.** Integration by parts leaves
  `[fg]_{−L}^{L}`: **[NUM]** for `f, g` both vanishing at `±L`, `⟨Df,g⟩ + ⟨f,Dg⟩ = −4.9e-13` and
  the boundary term is `0`; for `f = 1+u`, `g = 2−u`, which vanish at neither endpoint, the sum is
  `1.6000000` and **equals the boundary term exactly**. (Momentum-operator domain distinctions:
  Bonneau–Faraut–Valent, arXiv:quant-ph/0103153 — **NOT CHECKED**, cited as supplied in review.)

So `[A,D] = −I` is a statement on a common core, **not** a claim about two globally defined
generators on one fixed-window space. `A` is bounded self-adjoint on `L²[−L,L]` with `‖A‖ = L`;
`D` generates translations on `L²(R)`, which is where the Sonin/scaling action actually lives. The
two operations differ in kind:

| operation | what it does |
|---|---|
| translation in the logarithmic coordinate (`e^{tD}`) | **moves** the function along that coordinate |
| multiplication by `e^{−δu}` (`e^{−δA}`) | **changes the weight** of its values without moving their positions |

So of the brief's §4 outcomes: **outcome 1 (exact generator correspondence) holds for `δ`** with
the generator `A`, and **outcome 4 (incompatibility) holds for `η`** — no fixed `A` satisfies
`e^{ηA}+e^{−ηA}` matching, because `Φ(½tanh η)` would require `A² = u² tanh²η/(4η²)`, which is
`η`-dependent and hence not an operator. Outcome 3 (formal Taylor similarity) is explicitly **not**
claimed as a bridge; see the control in §7.

## 6. `F'(γ)² + F(γ)F''(γ)`: what it is, and what it is not

**[STD]** The calculus identity, exact over `Q` for an arbitrary quartic `F` (`checks.py` §2;
the wrong-sign variant `F'² − FF''` fails, as a control):

```
(1/2) d^2/dgamma^2 F(gamma)^2 = F'(gamma)^2 + F(gamma) F''(gamma) =: A_f(gamma) ,
```

so `ΔQ = −2δ²(d²/dγ²)F(γ)² + O(δ⁴)`. **This is elementary and is not presented as new.** It gives
`A_f` a curvature reading: *half the second derivative of the squared spectral profile*. Whether
that reading is **native** to the scaling/Sonin framework: `d²/dγ²` is multiplication by `(iu)² =
−u² = −A²`, so

```
A_f(gamma) = - (1/2) \widehat{ A^2 (f * f) } (gamma)   up to the convolution convention,
```

i.e. `A_f` is a **second moment of the generator against the autocorrelation**. That is a
legitimate operator reading and it is the honest limit of what the curvature form buys.

**Of the five candidate forms the brief lists, the outcome is:**

| candidate | verdict |
|---|---|
| a derivative of a square | **YES**, exactly: `½(F²)''` |
| a second derivative of `F²` | **YES**, same statement |
| a generator variance | **NO** — a variance is `≥ 0`; `A_f` is not (§7) |
| a quadratic form involving `A²` | **YES in form** (`−½·A²` moment of the autocorrelation), but the form is **indefinite** |
| a projected generator norm **plus an explicit correction** | **PARTLY**: `F'(γ)² = \|⟨f, AP_Lσ_γ⟩\|²` is a genuine projected-generator object, but the correction `F(γ)F''(γ)` is not of that type and carries either sign |

## 7. The projected-generator reading, with its exact scope  **[STD, relabelling]**

**[SOURCE: Theorem 8]** `K(L,γ) = sup_{f admissible}|F'(γ)|² = ∫_{−L}^{L}u²sin²(γu)du`, attained
at `f ∝ u sin(γu)`. Reading `σ_γ(u) := sin(γu)` and `P_L :=` multiplication by `1_{[−L,L]}`:

```
K(L, gamma) = || A P_L sigma_gamma ||^2_{L^2}  —  a PROJECTED generator norm.
```

**[NUM]** quadrature versus the paper's closed form agrees to `2.7e-11` at
`(L,γ) ∈ {0.5,0.8,1.0} × {3, 14.1347}` (`checks.py` §6).

**The projection is not optional.** `‖Aσ_γ‖²` over all of `R` **diverges**: partial values
`3.39e+02`, `3.33e+05`, `3.32e+08` on `[−10,10]`, `[−100,100]`, `[−1000,1000]`. So the window
projection is what makes the coefficient *exist*, which is the sharpest form of
`PROJECTION_CONTEXT.md`'s lesson available here: not "`‖PAψ‖² ≤ ‖Aψ‖²`" but "`‖Aψ‖² = ∞` and only
the projected quantity is defined."

**The brief's §6 identification list, answered concretely.** Hilbert space `L²([−L,L])` (equally
`P_L L²(R)`); generator `A =` multiplication by `u`, **bounded self-adjoint**, `‖A‖ = L`, domain
all of `L²[−L,L]`; projection `P_L` — multiplication by an indicator, **independently justified**
because the support constraint is part of the problem statement, not an added hypothesis; vector
`σ_γ = sin(γ·)`, which is **not** in `L²(R)` — so `AP_Lσ_γ` is defined only after projection, and
the order matters.

**But this is a relabelling, not a result.** `K = ‖AP_Lσ_γ‖²` is the paper's own extremal value
rewritten; the only added content is the identification of `s_γ(u) = u sin(γu)` as `Aσ_γ`. And the
crude bound it suggests is **weaker** than the paper's: **[NUM]** `L²‖P_Lσ_γ‖²` versus `2L³/3`
versus the true `K` at `γ = 3`: `(0.119, 0.0833, 0.0646)` at `L = 0.5`; `(0.618, 0.341, 0.264)` at
`L = 0.8`; `(1.047, 0.667, 0.324)` at `L = 1.0`. The paper's constant wins everywhere.

**Why the full coefficient is not a projected norm — control 4.** **[NUM]** `A_f` takes **both
signs** (`checks.py` §6): positive for `cos(πu/2L)` (`1.56e-04`), the tent (`1.63e-04`) and
`u sin(γu)` (`3.41e-02`); **negative** for `cos(γu)` (`−1.21e-01`), `cos²(γu)` (`−5.24e-03`) and
`f ≡ 1` (`−8.80e-03`), at `L = 0.8`, `γ = 14.1347`. A projected generator norm is `≥ 0`, so
`A_f` **cannot** be one. This is consistent with the manuscript under audit: Theorem 8 completes
the square precisely because `FF''` may be adverse, and §B.2 bounds `|F'²+FF''|` in absolute value.

> **Conclusion of §§6–7.** The projected-moment theorem of `PROJECTION_CONTEXT.md` **remains
> contextual only.** There is a real projected-generator object here — `K(L,γ) = ‖AP_Lσ_γ‖²` — but
> it is the paper's own sharp constant for the *derivative-evaluation step alone*, exactly as
> Remark 10(a)–(c) already warns, and the full sensitivity coefficient is sign-indefinite and
> therefore not of that type. No step above imports the finite-matrix result.

## 8. Fake-generator control  **[NUM, brief control 2]**

Any even analytic `g` is a power series in the squared parameter, so one can always *write*
`g(x) = 2I + x²A² + (x⁴/12)A⁴ + ⋯` by **defining** `A² := a₁I`, `A⁴ := 12a₂I`, …. For a single
operator this is only consistent if `a₂ = a₁²/12`. Witness `g(x) = cos x + e^{−x²}` (even,
analytic, nothing to do with scaling): `a₁ = −2`, `a₁²/12 = 0.333333`, `a₂ = 0.541667` —
**inconsistent**, so no single generator reproduces it (`checks2.py` §9).

**Therefore an even Taylor expansion is not evidence of a scaling generator.** The quartet case is
different in kind: Proposition 1 is an **identity at all orders**, checked as an identity
(`checks.py` §3), not a coefficient match. That is why `δ` qualifies and `η` does not.

## 9. Composition control  **[NUM, brief control 3]**

A system where `η = arctanh(2δ)` is valid algebra and has no native additive dynamics: a single
Bernoulli trial with `p = ½+δ`. For `(p₁,p₂) = (0.6,0.7)`, the natural operations give
`δ(p₁p₂) = −0.0800` (AND) and `δ((p₁+p₂)/2) = +0.1500` (equal mixture, `= (δ₁+δ₂)/2`), while
`δ₁ ⊕ δ₂ = +0.2778` matches **neither**; likewise for `(0.75,0.55)` (`checks2.py` §10). Exact
translation, empty dynamics — the manuscript's own category.
