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
