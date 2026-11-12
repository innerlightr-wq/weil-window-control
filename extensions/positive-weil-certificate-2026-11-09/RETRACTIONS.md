# RETRACTIONS.md

Corrections to the first build, after external review on 2026-11-09. Each entry
states what was claimed, what is true, and how it was checked here.

---

## R1 — The proposed sparsity target is FALSE. **Retracted in full.**

**Claimed** (REPORT.md v1, §6, "single best next step"):

> Find `c > 0`, `eta > 0`, `theta < 1` with
> `meas{ t in [T,2T] : P_L(t) > theta A_L } <= c T^{1-eta}` uniformly in `L`.

**This cannot hold, at any fixed `L`.** The reviewer's argument is correct and I
reproduce it:

Let `p_1 < ... < p_r` be the primes below `e^{2L}` and write the comb as a
function on the torus,

```
    P_L(t) = Pi(t log p_1, ..., t log p_r),
    Pi(x)  = sum_{p_i^k < e^{2L}} (2 log p_i / p_i^{k/2}) cos(k x_i) .
```

`Pi` is continuous with `Pi(0) = A_L`, so for each `theta < 1` the set
`U_theta = {Pi > theta A_L}` is **open and contains 0**, hence `mu(U_theta) > 0`
for Haar measure `mu`. The `log p_i` are linearly independent over `Q` (unique
factorization), so by the **continuous Kronecker–Weyl theorem** the line
`t -> (t log p_1, ..., t log p_r)` equidistributes in `T^r`. Therefore

```
    lim_{T->inf} meas{ t in [0,T] : t in E_{L,theta} } / T  =  mu(U_theta) > 0 ,
```

so `meas E_{L,theta}(T)` has **positive lower density** and no `c T^{1-eta}`
bound with `eta > 0` is possible.

**Checked numerically** (`src/refutation.py`, status NUMERICAL, corroborating a
proof rather than replacing it). At `L = 4/5`, primes `{2,3}`, `A_L = 2.941973525`:

| `theta` | Haar `mu(U_theta)` | time density on `[T,2T]`, `T = 10^2 … 10^6` |
|---|---|---|
| 0.50 | 0.129458 | 0.132470, 0.129095, 0.129440, 0.129425, **0.129480** |
| 0.80 | 0.045786 | 0.047595, 0.046125, 0.045795, 0.045870, **0.045740** |
| 0.90 | 0.022147 | 0.026110, 0.022750, 0.022090, 0.022140, **0.022145** |

The time density converges to the Haar measure to four significant figures over
four decades of `T` and is **flat, not decaying**. Exactly as the proof predicts.

**What was wrong in my reasoning.** I inferred "sparse" from the Weyl argument,
but Weyl equidistribution says the near-resonant set is *visited with a definite
positive frequency*. It is sparse in proportion, never in density, and that
proportion does not shrink with `T`. Lemma 3.2 bounds the supremum; I wanted the
distribution to decay, and it does not — it converges to a positive constant.

**A second, independent obstruction** to any purely density-based route, which
the reviewer also raises and which I confirm: `g(u) = f(u) cos(omega u)` has the
**same support and same parity** as `f`, with
`G(t) = [F(t-omega) + F(t+omega)]/2`. The Fourier profile can therefore be
deliberately translated onto a chosen resonance while the support constraint is
untouched. So a bound on `meas(E)` alone can never control `int_E |F|^2`.

**Replacement direction** (not a claim, a target):

1. effective equidistribution with the **correct positive-density main term**
   `mu(U_theta)` plus an error term, not a power saving;
2. a bound on how much Fourier energy an admissible *window-supported* `f` can
   place on `E` — the real question, and one the modulation remark shows is not
   answered by density;
3. an explicit inequality in which the favourable region absorbs that energy.

The open question is therefore: **can the support restriction on `f` limit
concentration on the unfavourable set enough to improve the weighted lower
bound?** R1 does **not** settle that, and does **not** cap every extension of
the detector.

## R2 — The exact multiplier belongs to the `T_1` scale only. **Scope corrected.**

**Claimed** (REPORT.md v1 §3 table): the `T_gate` row carried ratio `4.2186…`,
which reads as the exact factor.

**True.** The exact factor `exp(2 log 5 / sqrt 5) = 4.218691593787…` applies to
`T_1 = 2 pi e^{A_L}` and to nothing else. A gate that involves the digamma, or
the `-1/T` correction, is not a pure exponential in `A_L` and its ratio is not
exact. Recomputed in Arb:

| gate | `L = 4/5` | `L = 81/100` | ratio | exact? |
|---|---|---|---|---|
| `T_1 = 2 pi e^{A_L}` | `119.0865561877` | `502.3894535220` | `4.218691593787` | **yes**, to `1e-20` |
| `T_gate` (digamma) | `119.0869060760` | `502.3895364590` | `4.218679895322` | **no** |
| `T_paper` (`-1/T`) | `120.0824039145` | `503.3884609107` | `4.192025180219` | **no** |

The reviewer's `502.39/120.08 ≈ 4.18379` is precisely
`T_gate(81/100)/T_paper(4/5) = 4.183706522205`, a cross-definition ratio that
means nothing. My table invited that reading by listing three gates with a
shared "exactly" caption. Only the `T_1` row is exact.

## R3 — The 69 cells cover `[9/100, 1/4]`, not `[9/100, 49/100]`. **Corrected.**

**Claimed** (REPORT.md v1 §2): "69 cells, `Q_add,delta < 0` on `[9/100, 49/100]`".

**True**, from the verifier's own output:

* **recomputed here**: the 69 rational cells covering `[9/100, 1/4]`, each with
  `Q_add,delta(f_x) < 0` throughout the closed cell (worst unit-`L2` margin
  `1.276748e-06` on cell 7);
* **imported**: the published certificate on `[1/4, 49/100]`. Of it, two things
  were recomputed — the join at `delta = 1/4` (inherited frozen vector gives
  `<= -2.521632e-03`) and one interior point `Q_add,2/5 <= -0.0145107657498227`.

The union `[9/100, 49/100]` is certified only by **combining** the recomputed
and the imported pieces. The repository theorem separates them and so should I.

## R4 — A `923 x 923` assembly is an attempt, not a guarantee. **Softened.**

**Claimed** (REPORT.md v1 §3): executing the assembly is "a pure compute cost …
with no new mathematics".

**True.** Satisfying the envelope gate is necessary, not sufficient. Whether the
resulting finite reduction is positive after the tail errors `eps_D, eps_B` and
the complement bound are included is **not** settled by the gate. The right
description is an *attempted* certificate with an unknown outcome.

## R5 — Sequencing. **Accepted.**

The reviewer's recommendation is adopted: **complete the `L = 0.8` lower-bound
reproduction before scaling to `0.81`.** Until the pipeline has produced one
whole-window lower bound, scaling it 75x demonstrates nothing. This supersedes
the "next step" of REPORT.md v1 §6, which R1 has retracted anyway.

## R6 — "no bound on meas(E) can ever control int_E |F|^2" is FALSE. **Corrected.**

**Claimed** (RETRACTIONS R1 and REPORT v2 §6): a second obstruction, namely that
*no* bound on `meas(E)` can control `int_E |F|^2`.

**That is too strong and it is wrong.** Finite measure alone already gives
control, and a uniform *local* occupation bound gives much better control. Both
are elementary and both are proved and tested in `phase2/src/lemmas.py`:

> **Lemma A (finite measure).** For measurable `E` of finite measure and any
> `f in L^2[-L,L]`, with `F(t) = int f e^{itu} du`,
> ```
>     int_E |F|^2  <=  min(1, L meas(E)/pi) * int_R |F|^2 .
> ```
> *Proof.* `|F(t)|^2 <= 2L ||f||_2^2` by Cauchy–Schwarz, so
> `int_E |F|^2 <= meas(E) 2L ||f||_2^2`, and `||f||_2^2 = (1/2pi) int_R |F|^2`
> by Parseval. □

> **Lemma B (uniform local occupation).** If `meas(E ∩ I) <= eps h` for every
> cell `I` of one fixed partition of `R` into intervals of length `h`, then
> ```
>     int_E |F|^2  <=  q int_R |F|^2 ,   q = min(1, eps (1 + 2 L h)) .
> ```
> *Proof.* `G = |F|^2` is `C^1` because `f` and `u f(u)` are both in
> `L^2[-L,L]`, so `F, F'` are in `L^2` and `F` is entire of exponential type `L`.
> On a cell, `h sup_I G <= int_I G + h int_I |G'|`; summing,
> `int_E G <= eps (int_R G + h int_R |G'|)`; and
> `int_R |G'| <= 2||F'||_2||F||_2 <= 2L int_R G` since
> `||F'||_2^2 = 2pi||u f||_2^2 <= L^2 ||F||_2^2`. □

Both are **standard estimates / derived tools**, not new arithmetic results.
Numerically checked on explicit `f` and explicit `E` with `q_B` informative
(`q_B = 0.09…0.26` against actual ratios `0.05…0.10`).

**What survives of the original point, correctly stated:** *global asymptotic
density* supplies no uniform control of weighted energy, because density is a
limiting statement and says nothing about where the mass sits at any finite
scale. Positive limiting density does **not** prevent a useful local occupation
estimate. R1's refutation of the power-saving target is unaffected — that
remains correct.

**The modulation statement, made exact.** Multiplication by `e^{i omega u}`
translates `F` exactly and preserves support, but does **not** preserve real
parity. Multiplication by `cos(omega u)` preserves support **and** real parity
and gives

```
    G(t) = [ F(t + omega) + F(t - omega) ] / 2 ,
```

a **superposition** — with possible interference between the two translates and
a changed `L^2` norm — not an unrestricted concentration theorem. The first
build described it as if it placed the whole profile on a resonance. It does
not.

Tests must not encode the false impossibility as an invariant; verifier group
V7 is restricted to the *density* statement, which is the part that is true.

## R7 — "the remaining blocker is not mathematics". **Withdrawn.**

**Claimed** (REPORT v3 §7.3): "Blocker class: resource / implementation, not
mathematics."

**Withdrawn.** The correct statement is: **the finite block remains uncomputed,
and the success of this reduction is unproved.** Clearing the exterior gate and
bounding the complement does not establish that the compression of `R_T` is
positive definite. `R_T` is obtained from `Q` by *discarding* the exterior
integral — a favourable term — so even positivity of `Q` itself would not
automatically give positivity of `R_T`. The comparison can fail while `Q` is
positive, and nothing proved before this run excluded that.

## R8 — the Schur threshold is not a lower bound, and needs strict inequalities

**Claimed** (REPORT v3 §7.2): the chain "closes as soon as `lam_A >= 1.0146e-28`".

**Corrected on two counts.**

1. **Strictness.** With `A >= a I`, `D >= d I`, `||B|| <= kappa`, positive
   definiteness needs `a > 0`, `d > 0` **and** `a d > kappa^2`. Equality gives
   only a non-negative comparison, not a positive one.
2. **The threshold is not the bound.** `kappa^2/d` is the break-even point, not
   the resulting constant. Using the rounded reported numbers,
   `kappa^2/d ≈ 1.01456382257e-28`; if `a` were exactly `1.0146e-28` the stable
   estimate `(a d - kappa^2)/(a + d)` would be only about `3.62e-33`. (Those two
   figures are explanatory arithmetic on rounded inputs, **not** certified.)

The certificate therefore does **not** use the `2x2` eigenvalue formula, which
cancels catastrophically near break-even. For a target `mu` it certifies the
three strict interval signs

```
    a > mu ,     d > mu ,     (a - mu)(d - mu) > kappa^2 ,
```

which is exactly the condition for
`<(R_T - mu)x, x> >= (a-mu)|x_N|^2 - 2 kappa |x_N||x_perp| + (d-mu)|x_perp|^2 >= 0`.
`a`, `d`, `kappa`, the penalty `kappa^2/d` and `mu` are reported **separately
for each parity sector**, and the whole-window constant is the smaller of the
two.

**Naming.** `N = 220` means Legendre degrees `0..219`, i.e. **110 modes in each
parity sector**, not 220 per sector. Even-sector constants are not reused for
the odd sector: there the pole term enters with `s = -1`, so it is *subtracted*
in `d` rather than dropped as favourable.

**What each constant bounds.** `beta` and `M` concern the scalar symbol
`Psi_L - beta` on `|t| >= T` and `|t| <= T` respectively. `g_N` bounds the
Fourier energy fraction inside `[-T,T]` of functions in `V_N^perp`, i.e. a
property of the *basis complement*, not of `Q`. `a` is a lower bound for the
least eigenvalue of the **exact compression of `R_T`** to `V_N` — not of `Q`,
and not of the infinite comparison operator. Quadrature and rounding errors are
charged **once**, as an entrywise widening of the assembled matrix before the
`LDL^T` factorization, so the certified `a` already contains them.

## R9 — the cosine-basis claim was overstated

**Claimed** (commit message, phase 2): the cosine basis "provably cannot reach
the `1e-20` the coupling needs".

**Corrected.** What was shown is that the Cauchy–Schwarz surrogate for the
complement in the cosine basis decays only like `O(1/N^3)`, so reaching `1e-20`
by that route would need `N` of order `10^7` — a **demonstrated cost**, under a
fixed budget, for one particular bound. That is not a proof that no cosine-basis
argument, or no sharper complement bound in that basis, can reach an arbitrary
positive tolerance. No nonzero obstruction was proved. Legendre was chosen
because its complement bound is cheap, not because the alternative is impossible.

## R10 — the occupation no-go was overstated

**Claimed** (REPORT v3 §7.4): "no occupation bound was consumed in any
certificate, **and none can be** at this margin."

**Corrected.** The first half stands: none was consumed. The second half does
not. Taking `b = M = 8.8326` is a *conservative sufficient* condition on `q`,
not a necessary one: `b` is defined as `sup_E (beta - Psi_L)` over the actual
unfavourable set, which can be far smaller than the global `sup|Psi_L - beta|`,
and the actual deficit can be correspondingly smaller.

Moreover the set that matters, `E = {|t| >= T : Psi_L(t) < beta}`, is **not**
the fixed-threshold prime-comb resonance set of R1, and once the exterior gate
holds it is *bounded* — contained in a finite band by monotonicity of `H`. The
positive limiting density established in R1 for the comb set therefore does
**not** transfer to `E` and is not a no-go theorem for it.

What is true and was measured: the particular pairing of `b = M` with Lemmas A/B
is about 18 orders of magnitude short at this margin. That is a statement about
one sufficient route, not about every occupation estimate.

## R11 — the "lower beta to buy a larger `a`" next step. **Refuted, withdrawn.**

**Proposed** (end of the phase-3 report): lower `beta` below `H(T) - A_L`,
because the paper's smaller `beta*` "makes the interior easier", and recover a
larger certified `a`.

**That has the wrong sign.** Writing the complete comparison form and applying
Parseval,

```
    R_{T,beta}(f) = pole(f) + (1/2pi) int_{|t|<=T} Psi_L |F|^2
                            + (beta/2pi) int_{|t|>T} |F|^2 ,
```

so for `beta_2 >= beta_1`

```
    R_{T,beta_2}(f) - R_{T,beta_1}(f) = ((beta_2-beta_1)/2pi) int_{|t|>T}|F|^2 >= 0 .
```

`R` is monotone **non-decreasing** in `beta`, and so are its spectral infimum
and the least eigenvalue of its compression to any fixed subspace. In matrix
form, on one sector,

```
    A(beta_2) - A(beta_1) = (beta_2 - beta_1)(I - K_T) ,
    (K_T)_{nm} = (1/2pi) int_{-T}^{T} v_n v_m dt ,     0 <= K_T <= I ,
```

the bounds on `K_T` being immediate from
`<K_T y,y> = (1/2pi) int_{-T}^{T}|F|^2 in [0, (1/2pi)int_R |F|^2] = [0, ||y||^2]`.

Verified numerically (`phase3/src/audit.py`): the identity holds entrywise to
`5.2e-143`, `A(0.55) - A(0.45)` is certified positive definite, and the reverse
difference is certified **not** positive definite.

**My error:** I compared only `Psi_L - beta` inside the interior integral and
concluded a smaller `beta` makes it less negative. That ignores the `beta||f||^2`
term, which moves the other way and dominates by exactly the identity above.
Any apparent gain from lowering `beta` would have to come from enclosure quality
or another changed parameter — never from an operator improvement.

Separately: `beta = H(T) - A_L` is a **sufficient** exterior bound, not
automatically the exact infimum of `Psi_L` on `|t| >= T`; no optimality is
claimed for it.

## R12 — the "silently swallowed exception" characterisation. **Corrected.**

**Claimed** (phase-3 report §8.3 and commit message): FLINT's adaptive
integrator "silently swallowed an exception ... and returned a wrong value
rather than failing."

**I could not reproduce that, and the characterisation is withdrawn.** A minimal
reproducer (`phase3/src/callback_repro.py`, python-flint 0.9.0, Python 3.14.7)
gives, on every path tested:

| case | behaviour | safe? |
|---|---|---|
| callback raises a Python exception | propagates as `ValueError` | yes, loud |
| callback performs an invalid op (`1/sqrt` at a ball containing 0) | returns `nan` after exhausting the evaluation limit | yes, detectable |
| callback returns a non-finite value (documented protocol) | returns `nan` | yes |
| the exact unhardened phase-2 callback, on `[1,6]` and on `[1e-5,5]` | raises `SystemError` | yes, loud |

**What is real:** when a Python exception escapes the callback, the failure
surfaces as a `SystemError` ("returned a result with an exception set") only
once control returns to the interpreter, rather than as the original exception
at the point of failure. That is a **binding-level** issue in the callback
adapter — FLINT's C core has no knowledge of Python exceptions and continues on
the adapter's return value. It is not a FLINT core defect and not a defect in
this repository's wrapper.

**What was wrong in my report:** the wrong values I attributed to silent
swallowing (`-0.000799204100841` where the truth is `-0.7110991221349`) came from
a run in which the matrix side still carried the **pole-phase bug** and the two
sides were computed at different `T`; I conflated interleaved output from two
runs. No instance of a silently-wrong finite return has been demonstrated.

**Unaffected:** the certificate itself. `acb.integral` is used **only** in
`rigorous_entry`, which is called only by the verifier's spot-check; the
assembly in `block.assemble` never calls it. Grep-verified.

The hardening is kept regardless: the callback records any failure in a flag and
the integration result is rejected on that flag, never on the returned value
alone.

---

## R13 — the upper endpoint was mis-quoted, twice. **Corrected.**

Found 2026-11-10 while assembling the technical note, by a checker written
specifically to test the **rounding direction** of every quoted bound
(`note/consistency_check.py`, checks U1–U6). Two separate errors, both in the
safe direction numerically, both wrong as stated.

### (a) The ball midpoint was quoted as the bracket's upper endpoint

`REPORT.md` §10, §10.5, §10.6, §11 and `README.md` reported

```
    1.41e-17 <= lambda_Q <= 1.90138291114546440e-17 .
```

`1.90138291114546440e-17` is the **midpoint** of the enclosure
`[1.90138291114546440e-17 +/- 6.03e-35]` recorded in
`phase4/certs/cross.json`. A midpoint is not an upper bound for its own ball:
the ball's upper endpoint is `1.901382911145464406030000e-17`, larger by about
`6e-34`. The overstatement is numerically negligible and changes no conclusion,
but the inequality as printed was false.

**Corrected form.** Quote the bracket as

```
    1.41e-17 <= lambda_Q <= 1.9014e-17 ,
```

which dominates every endpoint in play, and reserve the full-precision value
for the artifact, where it is consumed through `.upper()`.

### (b) "The exact upper endpoint" was computed at 53 bits

`phase5/certs/parity.json` recorded

```
    U_upper_endpoint = 1.9013829111454646892e-17
```

and `REPORT.md` §11.2 and `HYPOTHESES.md` B24 repeated it as *the exact upper
endpoint of the `q(f_400)` enclosure*. It is not. `phase5/src/parity.py`
computed `U = arb(arb(cross["Q_f400"]).upper())` **outside any owned precision
scope**, so the serialized ball was re-parsed at the ambient 53 bits: the
midpoint rounded down, and `.upper()` then inflated the endpoint by about
`2.8e-33`. The true endpoint is

```
    U = 1.901382911145464406030000e-17 .
```

This is the precision-ownership defect of `REPORT.md` §11.1 recurring **in the
very routine that documents it** — and specifically the half that was flagged
there as the harder one: *an input constructed outside the scope carries the
ambient precision's radius.* Re-parsing a serialized ball is exactly such an
input.

**Fixed:** the separation block in `phase5/src/parity.py` now runs inside
`owned(PREC)`, and `phase5/certs/parity.json` is regenerated with the true
endpoint, a corrected `separation = 3.986170888545355939700000e-18`, an explicit
`U_safe_decimal = "1.9014e-17"`, and an `erratum` field naming this entry.

### What is unaffected

Both errors were conservative, so nothing certified is withdrawn:

* the quoted midpoint understated the bound by `6e-34`, and the 53-bit endpoint
  overstated it by `2.8e-33` — the first direction is wrong but the claim it
  supported (`lambda_Q <= 1.9014e-17`) holds with room to spare;
* the parity separation used the **larger** `U`, so `U < mu_odd` was proved with
  a smaller margin than the data actually supports; the corrected separation is
  larger, not smaller;
* `phase4/verify_t400.py` and `phase5/src/pilot800.py` consume
  `cross.json`'s `Q_f400` as a **ball**, and a low-precision re-parse of a ball
  is still a valid (merely wider) enclosure, so every quantity derived from it
  remains valid.

### The general lesson, recorded

A serialized enclosure is a ball that *contains* the quantity. It must be
consumed through the correct endpoint — lower for a lower bound, upper for an
upper bound — and never through its midpoint or its printed decimal. Printing a
ball's midpoint and then using that decimal as a bound is a type error that no
amount of precision detects, because the arithmetic is all correct; only a check
of the **direction** catches it. That check now exists.

---

## R14 — the time-domain evaluator omits its outer quadrature remainder. **Upper bound withdrawn and replaced.**

Found 2026-11-10 while tracing whether a mis-stated polynomial degree had
propagated into an error bound. It had not; something worse had.

### The defect

`phase4/src/compact.py::arch_compact` evaluates

```
   int_0^{2L} 2[e^{-2x} g(0) - e^{-x/2} g(x)] / [1 - e^{-2x}] dx
```

by panelled Gauss–Legendre — geometric panels from `eta = 1e-40` to `1e-3`, then
32 uniform panels of order 32 — and returns

```
   const + tot + head
```

where `tot` is the bare quadrature sum and `head` bounds only the `[0, eta]`
piece. **No term bounds the truncation error of the panelled rule.** The module
docstring asserted that "a Bernstein-ellipse remainder is used"; none is
computed anywhere in the file.

Ball arithmetic does not supply it. Arb encloses the nodes, the weights and
every rounding, so the returned radius is a true bound on *those* errors — but
the difference between a Gauss–Legendre sum and the integral it approximates is
a property of the rule, not of the arithmetic, and no interval operation
produces it. The published radius `+/- 6.03e-35` on `q(f_400)` is therefore
incomplete, not merely optimistic.

### What is withdrawn

Everything the time-domain route produced, as a **certificate**:

| quantity | status |
|---|---|
| `q(f_400) = 1.90138291114546440e-17` | **WITHDRAWN** as a certificate; retained as a measurement |
| `q(f_200) = 3.64859175290407702e-17` | **WITHDRAWN** as a certificate; retained as a measurement |
| upper endpoint `U = 1.901382911145464406030000e-17` | **WITHDRAWN** |
| the bracket `1.41e-17 <= lambda_Q <= 1.9014e-17` | **WITHDRAWN** |
| `e_T = q - r_T` for all `T`, and their ratios | **WITHDRAWN** as certificates; retained as measurements |
| the parity separation argument *as run on `f_400`* | **RE-ESTABLISHED** on another witness, below |

The artifacts are **not** deleted. `phase4/certs/cross.json` and
`phase4/certs/q_f200.json` keep their recorded values; what changes is their
status, recorded here and in `note/CLAIM_INDEX.md`.

### What is unaffected, and why

* **The lower bound `1.41e-17`.** Its interior quadrature carries a certified
  Bernstein-ellipse remainder, folded into the matrix entries exactly once
  (`phase2/src/block.py`, `quad_error_factor` in the certificate meta), and
  `phase4/verify_t400.py` re-derives it. Different routine, complete budget.
* **The odd-sector bound `2.3e-17`**, from the same matrix route.
* **The scalar `R_T(f)` values.** `phase5/src/scalar.py` *does* add a certified
  remainder and publishes it as `*_quad_err`. In particular the ceiling
  `lambda_R(800) <= r_800(f_400)` never involves `Q`.
* **`beta`, `M`, `kappa`, `g_N`, `S_N`,** and every constant in the chain.

### The replacement

The phase-1 route is independent of `compact.py` and performs **no quadrature at
all**: `src/certify.py` evaluates the form in closed form through digamma,
Lerch `Phi` and Hurwitz `zeta`, each as an Arb ball. Replayed under owned
precision from its exact rationals (`phase5/src/replay_upper.py`, 10 checks),
the frozen 28-mode even cosine witness at `L = 4/5` gives, identically at 1200
and 1600 bits,

```
   Q(f_28)/||f_28||^2 = 2.112992611721490779453e-17 +/- 1.9e-354 ,
```

outward upper endpoint `2.112992611721490779453e-17`. The published bracket is
therefore

```
   1.41e-17 <= lambda_Q(4/5) <= 2.1130e-17 ,
```

width `7.03e-18` (was `4.914e-18`). It remains strictly nested inside Zhu's
certified `[8.9e-18, 2.27e-17]`, now `49.1%` narrower rather than `64.4%`.

**The parity separation survives.** `2.1130e-17 < 2.3e-17`, so with the
closed-form witness — which is even —

```
   inf_even <= 2.112992611721490779453e-17 < 2.3e-17 <= inf_odd ,
```

margin `1.870073882785092e-18` (was `3.986e-18`). The conclusion that the even
infimum determines `lambda_Q` is unchanged.

### The lesson

A quadrature routine that returns an interval is not thereby certified. The
interval bounds the arithmetic; the *rule* needs its own remainder, derived and
added explicitly. Two routines in this build compute panelled Gauss–Legendre
integrals: `scalar.py` adds a remainder and publishes it, `compact.py` did not,
and nothing in the test suite compared them on that point. The regression added
with the repair requires the allowance to be present in the certificate and
rejects a run that removes or understates it.

### R14 addendum, same day: the endpoint improves past the withdrawn value

Repairing the quadrature remainder turned out not to be the efficient way to
recover the bracket. The closed-form route has **no quadrature to bound at all**,
so a better trial function in that basis buys a better certified endpoint with
nothing new to prove. Pushing the even cosine basis from 28 modes to 160:

| modes `N` | certified `Q(f_N)/‖f_N‖²` ≤ |
|---|---|
| 28 | `2.1129926117214908e-17` |
| 48 | `1.9204310132e-17` |
| 64 | `1.8541624482e-17` |
| 96 | `1.7476544722e-17` |
| 128 | `1.6972723502e-17` |
| **160** | **`1.6731255341e-17`** |

The published bracket is therefore

```
   1.41e-17 <= lambda_Q(4/5) <= 1.674e-17 ,     width 2.64e-18
```

which is **narrower than the bracket the defect cost us**: the withdrawn
time-domain endpoint was `1.9014e-17`, and `1.6731e-17` beats it by `12.0%`.
Against Zhu's certified `[8.9e-18, 2.27e-17]` the width is now `80.9%` smaller,
a factor `5.23`, and Zhu's own numerical reference estimate `1.656e-17` for this
infimum lies inside it — the certified endpoint sits about `1%` above that
reference, which is the right side to be on.

Parity separation margin grows from `1.87e-18` to `6.268744659293181e-18`.

**The outer quadrature remainder remains open**, and is now off the critical
path: closing it would restore a second, independent route to the bracket rather
than improve it. `phase5/src/outer_convergence.py` records the diagnostic
(refinement difference `+/- 2.71e-35`, not a bound) and the reason the obvious
Bernstein route fails — `max|g| * R^{deg g}` is hopeless at `deg g = 877` — along
with the coefficient-sensitive bound that should replace it.

---

## R15 — precision ownership: a rule for NEW drivers, after three repeats on 2026-11-12

Not a retraction: nothing published here is withdrawn. This records a rule,
because the defect of `REPORT.md` §11.1 and of `R13(a)/(b)` recurred **four
times in one session** in new driver code written outside the repaired path, and each
time it produced a silently vacuous or silently unprovable certificate rather
than an error.

| # | what happened | consequence |
|---|---|---|
| A | a fresh `sector_dk` recomputed `S_tail`, `g`, `pole_vec_tail`, `kappa`, `d` inline, called before any `verify_pd_preconditioned` had set `ctx.prec` | `d = [+/- 3.09e+230]`, `kappa = nan` in the **even** sector only (the odd sector, called after a Rump test, was clean) — so `certify_final` could not close and the sector reported a **vacuous** `mu >= 0`. Exactly the corruption pattern of the disclosed defect in `phase4/certs/finite_block_certificate_T400.json`. |
| B | a shift search that doubles upward from a start guess and breaks on first failure | a start guess `1.4%` above the true `lambda_min` certified **nothing**, with no error |
| C | serialized matrix entries re-parsed with `arb(v)` after `constants()` had exited its `owned(prec)` scope and restored `ctx.prec = 53` | every entry radius inflated from `5.12e-49` to `7.18e-19`, a factor `1.4e30`; against a `lambda_min` of `8.3e-18` **no shift certified at all**, down to `7.5e-30` |
| D | `mu` bisected to the supremum (`a - mu* = kappa^2/d ~ 1e-58`) and then serialized with `.str(25)`, whose truncation radius is `~3e-46` | the published `mu` is a ball **straddling** `a`, so on replay `bool(a > mu)` is False and `certify_final` refuses. The Schur step had been performed correctly at full precision; the **artifact** did not support the claim |

### The rule

1. **Every new driver computes its certificate constants through
   `phase5/src/constants.py::constants(L_num, L_den, T, N, prec=...)`**, never
   inline. That function runs inside `owned(prec)` and applies `require_tight`
   to `beta`, `M`, `kappa`, `d_even`, `d_odd`. In case A it would have **raised**
   instead of writing a vacuous record.
2. **Any re-parse of a serialized ball names its precision**, by setting
   `ctx.prec` immediately before the parse or by parsing inside `owned(prec)`.
   `constants()` restores the ambient precision on exit, so code after it is at
   `53` bits unless it says otherwise. This is `R13(b)` applied to matrix
   deserialization rather than to a scalar endpoint.
3. **A shift search is two-sided.** A one-sided ladder cannot distinguish "the
   matrix is not positive definite at this shift" from "my guess was too high",
   and the two have opposite meanings.
4. **A published constant must be serializable with room to spare.** A `mu`
   bisected to the Schur supremum has margin `kappa^2/d`, which is `~1e-41` of
   `a` relatively and cannot survive any print width. Publish
   `mu = a*(1 - 10^-12)` truncated DOWN to 15 significant digits: it keeps
   `10^30` times the needed margin, loses `10^-12` of the constant, and replays
   exactly. A certificate that cannot be re-verified from its own artifact is
   not a certificate.

### What is unaffected

Nothing certified in phases 1-6 changes. Defects A-C occurred only in new
investigation drivers (`~/scratch/rh-connes-continuation-2026-11-12/work/`), and
all three were failures to *certify*, never false certifications: in each case
`verify_pd_preconditioned` sets its own precision internally, so the shifts it
reported as PASS were valid. Defect A's `a` was a valid pass whose chain was
lost; B and C certified nothing; D's underlying certification was sound and only
its serialization was unusable, now repaired in `certs/schur_margins.json`. The repaired results are in
`certs/certificate_T400_lam5_repaired.json` of that directory.
