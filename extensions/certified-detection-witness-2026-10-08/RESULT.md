# RESULT — verdict: **CERTIFIED FULL-FORM PLANTED DETECTION**

with one arithmetic gap labelled explicitly in §5. Exploration only, 2026-10-08. No manuscript
edit, commit, push, branch change, reset, upload, install or background job.

## The one-line answer

> **We established a negative detection margin for one actual admissible function, with its own
> baseline and all errors accounted for** — not merely a sharpened interpretation of the existing
> numerics. The function is explicit (four coefficients), its *complete* Weil form value is
> computed on the geometric side without assuming any zero location, and the negative interval
> `δ ∈ [0.2042, 0.5]` lies inside the critical strip with a margin `2.55×10^{-2}`, about
> `6.8×10^3` times the proved error allowance. The nodal constraint `F(γ) = 0` is what makes it
> work, and the Taylor budget never would have.

## 1. The witness

`L = 0.8`, `γ = 14` **exactly** (an exact chosen ordinate; no zeta-zero value is used), basis
`even_d` `φ_k(x) = cos((2k+1)πx/2L)/√L` (orthonormal, `Gram = I`), `E = span{φ_0,…,φ_3}` — an
explicitly frozen coefficient span, **no spectral projector is certified anywhere**.

```
f_node = sum_k x_k phi_k ,   x = ( +0.099815221908,  -0.354082045264,
                                   +0.928280739889,  +0.054384691304 ) ,   sum x_k^2 = 1
```

`f_node` is real, even, supported in `[−L,L]`, `‖f_node‖₂ = 1` — admissible in the paper's sense
(§4.2) — and, being in the `even_d` span, it is continuous and vanishes at `±L`.

| quantity | value |
|---|---|
| `F(γ)` | `−4.12×10^{−17}` (machine zero; **charged, not assumed**) |
| `F'(γ)` | `−0.176262303` |
| `F''(γ)` | `+0.113232083` |
| `a_f = 4[F'² + F F'']` | `0.124273598` |
| `4 K_eff` | `0.124273598` (equal, as Corollary 2 predicts) |
| `K = K(L,γ)` | `0.181590042` |
| `‖P_E r‖²` | `0.036016816` |
| `K_eff` | `0.031068399` |
| `b_f = Q_0(f_node) = Q_ζ(f_node)` | `0.005168243` |
| archimedean tail bound | `3.76×10^{−6}` |
| `b_upper` | `0.005172000` |

**What the constraint costs, separated.** `K_eff/‖P_E r‖² = 0.8626` — the node costs **13.7%** of
the available derivative sensitivity. `‖P_E r‖²/K = 0.1983` — the 4-dimensional basis costs
**80%**. So *the nodal constraint is cheap; the small basis is expensive*, and the two costs must
not be conflated. Against that, the node drops the baseline from `Q_ζ = 0.182154` (witness A, same
`E`) to `0.005168` — a factor **35** — and removes the `2F(γ)²` augmentation (`0.2134`) outright.
**That trade is why the witness detects.**

## 2. The certificate (Route X — exact, no Taylor remainder)

The finite-displacement quartet is exact:

```
Q_delta(x) = Q_zeta(x) - 2<x,c>^2 + 4<x,a_delta>^2 - 4<x,b_delta>^2 ,
a_delta(u) = cos(gamma u)cosh(delta u) ,  b_delta(u) = sin(gamma u)sinh(delta u) .
```

The term `+4⟨x,a_δ⟩²` is **positive and retained**, never dropped. Results:

| `δ` | negative part `−4⟨x,b_δ⟩²` | positive part retained | `Q_δ(f_node) ≤` |
|---|---|---|---|
| `0.204201571681` | `−0.005194375` | `+2.237×10^{−5}` | `0.000000000` |
| `0.25` | `−0.007794921` | `+5.035×10^{−5}` | **`−0.002572569`** |
| `0.35` | `−0.015330503` | `+1.944×10^{−4}` | **`−0.009964094`** |
| `0.50` | `−0.031514947` | `+8.185×10^{−4}` | **`−0.025524492`** |

> **Witness certificate (planted, full form).** With `f_node` and `L = 0.8` as above, the Weil
> functional of the zero configuration
> ```
> { nontrivial zeros of zeta }  ∪  { 1/2 ± delta ± 14 i }
> ```
> evaluated at `f_node` is **strictly negative for every `δ ∈ [0.20421, 0.5]`**, with value
> `≤ −0.02552` at `δ = 1/2`.

The interval is **nonempty** and lies inside the critical strip (`|δ| ≤ 1/2`, i.e. `Re s ∈ [0,1]`).
Because `⟨f_node, c⟩ = 0` to `4×10^{−17}`, the `−2⟨x,c⟩²` term contributes `< 10^{−32}`, so the
value equals the Weil functional of that configuration — ζ's own zeros, **plus** four extra
off-line zeros — and not of a signed pseudo-configuration. `Q_ζ` is the **complete** form,
assembled on the geometric side (pole, archimedean, `log π`, prime powers `n = 2,3,4`), so **no
zero sum is truncated** and no tail is dropped.

At the documented planted ordinate `γ₁ = 14.134725141734693790…` the same computation gives
`δ* = 0.188950152679` and `Q_δ ≤ −0.022364899` at `δ = 1/2`. **That row is discovery, not
certification:** this run does not rigorously enclose `γ₁`, and the two ordinates are never
interchanged. At `γ₁` the configuration is the physically meaningful one — ζ's zero at `γ₁`
*moved* off the line, since there `2F(γ₁)²` removes an actual pair — which is why it is worth
reporting even uncertified.

## 3. Why the Taylor route would not have worked

| `δ` | signal `a_fδ²` | `R_upper = (16/3)L⁵δ⁴e^{2Lδ}` | `b_upper − a_fδ² + R_upper` |
|---|---|---|---|
| `0.100` | `0.001243` | `0.000205` | `+0.004134` |
| `0.161` | — | — | `+0.003470` (the minimum) |
| `0.204` | `0.005172` | `0.004195` | `+0.004195` |
| `0.300` | `0.011185` | `0.022877` | `+0.016864` |
| `0.500` | `0.031068` | `0.243088` | `+0.217192` |

**The Taylor budget is positive at every `δ ∈ (0,0.5]`** — it never certifies anything. The
brief's suggested sufficient condition `R_upper ≤ 2K_eff δ²`, i.e. `δ²e^{2Lδ} ≤ 0.035555`, **fails
at every `δ` in the strip** (`0.0489` at `δ = 0.19`, `0.0579` at `δ = 0.2043`). So the
`(16/3)L⁵δ⁴e^{2Lδ}` remainder — whose coefficient was re-derived and verified — is simply too
coarse here, and the result rests entirely on the exact finite-displacement representation. The
remainder bound is *valid* at every `δ`; it is *uninformative* in the range that matters, exactly
the validity-versus-usefulness distinction Theorem 7's discussion already draws.

## 4. Pipeline validation, before any witness claim

`mpmath` and `numpy` are not installed and installs are forbidden, so `src/zeta_window.py` could
not be run. Its archimedean block uses digamma / Lerch Φ / Hurwitz ζ — precisely what
**Limitation 4** of the manuscript records as not interval-evaluable. The archimedean term was
therefore re-expressed as a series of **elementary** integrals with a **proved** tail bound
(`TARGET_AND_SOURCES.md` §4), using no digamma and no Lerch transcendent.

1. `K(L,γ)` closed form (from the **TeX**, see §6) vs quadrature: rel. dev. `2.3×10^{−14}`.
2. **Archimedean term, two independent representations**, `φ_0`, `L = 0.8`: exact m-sum
   `−0.98589817` (tail bound `1.18×10^{−5}`) vs spectral
   `(1/2π)∫|F|²Reψ(1/4+ir/2)dr = −0.98590193` with an independently implemented digamma
   (validated to `10^{−12}` on `ψ(1), ψ(1/4), ψ(1/2)`). **Agree to `3.8×10^{−6}`, inside the bound.**
3. **Against the repository's recorded spectrum**, `N = 4`: `λ₂` runs
   `6.46→3.09→2.786 ×10^{−6}` as `M` runs `2×10^5→2×10^6→10^7`, converging on the recorded
   `2.711086672×10^{−6}`; `λ_min` decays like the truncation bound, `3.75×10^{−6}→3.75×10^{−7}→
   7.52×10^{−8}`, toward the recorded `1.7339346963×10^{−10}`. **`λ_min` is below this pipeline's
   error bar and is not resolved** — which does not affect the target, since `Q_ζ(f_node)` is
   `5.2×10^{−3}`.
4. **The witness's own `Q_ζ`, two independent routes**: exact m-sum `0.0051682432 ± 3.76×10^{−6}`
   versus spectral `0.0051678677`. **Agree to `3.8×10^{−7}`.**
5. Normalisation: `(1/2π)∫|F_x|²dr = 0.9999999999`.
6. **Controls.** `δ = 0`: the exact response reduces to `+2⟨x,c⟩²` with residual **exactly 0**.
   `±δ`: response even in `δ`, deviation **exactly 0**. `Q_ζ(f_node) > 0`, consistent with Weil
   positivity (a negative value would have been a red flag against the pipeline).

## 5. The arithmetic gap, stated plainly

The certificate rests on (i) exact closed-form integrals, (ii) a **proved** archimedean tail bound
`(¾C(0) + ½sup|C'|)/(M−¾) = 3.76×10^{−6}` at `M = 2×10^6`, and (iii) IEEE **double** arithmetic.

**(iii) is not directed-rounding interval arithmetic, so this is not a formal interval
certificate.** Accumulated roundoff over the `2×10^6`-term m-sum is *estimated* at `~10^{−13}`
(not proved), and is independently corroborated by the two-route agreement at `3.8×10^{−7}`. The
detection margin `2.55×10^{−2}` exceeds the proved tail bound by `6.8×10^3` and the roundoff
estimate by `~10^{11}`. **The single step missing for a formal certificate is directed rounding in
the m-sum**; nothing else. Increased precision alone is not certification, and is not claimed as
such.

## 6. One extraction artifact caught, and not turned into a correction

The PDF text layer renders the closed form of `K(L,γ)` with the outer bracket dropped, which
flips two signs and disagrees with quadrature by 2%. The **TeX** (`paper/weil_window_control.tex`
lines 372–375) is `K = L³/3 − [ L²sin2γL/(2γ) + Lcos2γL/(2γ²) − sin2γL/(4γ³) ]`, which reproduces
quadrature to `2.3×10^{−14}`. **The manuscript is correct; the extraction was not.** No correction
is proposed. (Named sources not found: `extensions/density-stiffness-audit/` does not exist, and
the "97.09%/2.91% split" is nowhere in the repository — see `TARGET_AND_SOURCES.md` §2.)

## 7. The arithmetic boundary — what would be needed to reach an actual ζ quartet

Write, with the correct multiplicities, for a hypothetical off-line quartet
`ρ_± = ½ ± δ ± iγ` replacing a simple on-line zero at `γ`:

```
Q_actual(f) = Q_target_quartet(f) + Q_rest(f) ,
Q_target_quartet(f) = 4 Re F(gamma + i delta)^2 = 4[<f,a_delta>^2 - <f,b_delta>^2] ,
Q_rest(f)           = sum over every nontrivial zero OTHER than the four of the quartet .
```

**What is needed is an UPPER bound on `Q_rest(f_node)`.** The target's contribution at
`δ = 1/2` is `−0.0307` (its negative part `−0.03151` plus the retained `+0.00082`), so the
quartet's negative contribution dominates iff

```
Q_rest(f_node)  <  0.0307 .
```

That is an **upper** bound on a sum of squares, i.e. a statement that the *other* zeros place
little mass on `f_node` — and it is **not** the remaining-zero nonnegativity `Q̃ ≥ 0` of
**Limitation 2** and Theorem 8, which is a *lower* bound used for the sharper *lower* estimate.
The two requirements are different in direction and neither implies the other; conflating them is
the trap this project already documents.

**What the geometric-side evaluation does and does not certify.** `Q_ζ(f_node) = 0.005168` is
computed from the pole, the archimedean integral, `log π` and the prime powers — it accounts for
**all** of ζ's zeros without assuming any location. What it certifies is therefore a
counterfactual: *the value of the Weil functional at `f_node` for ζ's own zero configuration,
whatever that is*. It does **not** license substituting that baseline into a different
hypothetical configuration in which one of those zeros has been moved: moving a zero changes
`Q_ζ` itself by `−2F(γ)² + 4ReF(γ+iδ)²`, which is exactly the correction applied here, and the
`−2F(γ)²` piece is legitimate only when there is an actual zero at `γ` to remove. At `γ = 14`
there is none, which is why the certified statement is phrased as *adding* a quartet to ζ's zeros,
not as *moving* one.

**This is detection, one direction of the classical criterion, and nothing more.** Exhibiting an
`f` whose Weil functional is negative for a configuration containing an off-line quartet shows
that such a configuration violates Weil positivity. To turn that into an exclusion one would need
an **independently established positivity statement for the same actual form, domain and window**
that the negative witness contradicts — and `Q_ζ(f_node) = +0.005168 > 0` is precisely *not*
contradicted here, because the configuration tested is not ζ's. Nothing in this run asserts that
ζ has no off-line zero.

**The ONE missing estimate, with its quantifiers.**

> Find `θ > 0`, `L₀ > 0` and `N₀` such that: for every `L ≤ L₀`, every `γ` with `|γ| ≥ 2`, and the
> nodal witness `f_node = f_node(L, γ, N)` of `E = span{φ_0,…,φ_{N−1}}` with `N ≤ N₀`,
> ```
> Q_rest(f_node)  <=  theta * 4 K_eff(L, gamma, N)        with  theta < 1/4 ,
> ```
> where `Q_rest` is the sum over all nontrivial zeros except the quartet at `½ ± δ ± iγ`,
> **uniformly in `δ ∈ [0, 1/2]`**.
>
> Parameter dependence that matters: `θ` must not grow with `γ` (the zero density grows like
> `log γ`, while `K_eff = O(L³)` does not grow with `γ`), and must not grow as `N` increases (the
> measurements here show `Q_ζ` growing from `0.0052` at `N = 4` to `0.4486` at `N = 16` while
> `4K_eff` grows only from `0.124` to `0.670` — i.e. the ratio worsens with `N`, so any uniform
> statement must be at **fixed small** `N`).

With `θ < 1/4` the quartet's `−4K_eff δ²` at `δ = 1/2` would beat `Q_rest`, giving detection for
every such `γ`. This is a single, bounded estimate; no attempt is made here to solve the
remaining-zero problem in general.

## 8. Scope distinctions carried forward

- `K` is the derivative-evaluation norm, **not** the sharp coefficient of the re-optimised
  minimum (paper Remark 10(a)–(c)). Nothing here claims the measured saturation of §4.3, which is
  T2.
- `K/L³ → 1/3` only as `|γ|L → ∞`; it is **not** bounded below by `1/3` for every `γ`. At
  `L = 0.8, γ = 14`, `K/L³ = 0.3547`.
- No Landau–Widom fit, no `97.09/2.91` split, no fitted coefficient of any kind appears.
- Fixed operator rank does not rule out cubic growth of a perturbative coefficient; no rank
  argument is used.
- **Finite-model proof ≠ full-form witness certificate ≠ uniform continuum theorem ≠ arithmetic
  zero-exclusion.** This run delivers the **second**, at one window, one ordinate, one basis size,
  modulo §5's arithmetic gap. It delivers none of the other three, and in particular no negative
  eigenvalue of a truncated model is offered as evidence of anything.

## 9. Files, commands, final state

```
python3 certify_witness.py        # library: closed forms, exact Laplace integrals, tail bound
python3 run_certify.py            # driver -> results.json  (~12 min)
python3 crosscheck_Qzeta.py       # independent spectral cross-check of Q_zeta(f_node)
```

| file | what |
|---|---|
| `TARGET_AND_SOURCES.md` | conventions verified, named sources located or reported missing, the frozen plan, the route around Limitation 4 |
| `WITNESS_THEOREM.md` | Lemma 1 and Corollary 2 with proofs, the inequality direction, the two error routes, the cost comparisons |
| `certify_witness.py` | the library |
| `run_certify.py` | the driver |
| `crosscheck_Qzeta.py` | the independent cross-check |
| `results.json` | 12 witness records, both ordinates, `N ∈ {4,8}` frozen + `N = 16` trend, with validation block |
| `RESULT.md` | this file |
| `SOURCE_HASHES.txt` | md5 of every input relied on |

**Unchanged-file check and final git status are appended below by the run that produced this
file.**

### Unchanged-file check (`md5sum -c SOURCE_HASHES.txt`), run 2026-10-08

```
src/zeta_window.py: OK                                    notes/proofs.md: OK
src/stage3e_certify.py: OK                                paper/weil_window_control.tex: OK
src/delta2_law.py: OK                                     paper/weil_window_control.pdf: OK
src/weilform.py: OK                                       data/stage3_gateB_lambda_min.csv: OK
data/stage3_delta2_law.csv: OK                            data/stage3_gateA_explicit_formula.csv: OK
extensions/second-moment-generalization-2026-10-07/PROJECTION_CONTEXT.md: OK
extensions/newtonian-weil-response-2026-10-08/RESULT.md: OK
```

### Final `git status`

```
HEAD = cb562f7   branch = main   (unchanged)
?? extensions/certified-detection-witness-2026-10-08/     (this run)
?? extensions/newtonian-weil-response-2026-10-08/         (the previous run, untouched)
```

Nothing tracked was modified; no commit, push, tag, branch change, reset, upload, install or
background job.
