# Source audit

What was actually read, in what version, and which statement was taken from where.
Nothing below was inferred from an abstract where the body was needed.

## Primary sources

### S1. D. H. J. Polymath, *Effective approximation of heat flow evolution of the Riemann ξ function, and a new upper bound for the de Bruijn–Newman constant*
arXiv:**1904.12438**. **Body read** via the arXiv HTML rendering. Equations used, by number:

| taken | from |
|---|---|
| `H_0(z) := (1/8) ξ(1/2 + iz/2)` | **eq. (1)** |
| `Φ(u) := Σ_{n≥1} (2π²n⁴e^{9u} − 3πn²e^{5u}) exp(−πn²e^{4u})` | **eq. (3)** |
| `H_t(z) := ∫_0^∞ e^{tu²} Φ(u) cos(zu) du` | **eq. (4)** |
| `∂H_t/∂t = −H_t''` | **eq. (31)**, §3 "Dynamics of zeroes" |
| "no zeroes `H_t(x+iy)=0` with `x∈ℝ` and `y > max(y_0² − 2(t−t_0), 0)^{1/2}`" when `t > t_0` and no zeros with `y > y_0` at `t_0` | **Theorem 3.2**, §3 |
| "there exists a finite constant `Λ` … such that the zeroes of `H_t` have purely real zeroes if and only if `t ≥ Λ`" | §1, before Theorem 1.1 |

The deformation is described there as "the evolution of `H_0` under the **backwards** heat
equation". **No best-bound claim is used anywhere in this audit**, and the project's
numerical upper bound on `Λ` is not relied upon.

### S2. B. Rodgers and T. Tao, *The De Bruijn–Newman constant is non-negative*
arXiv:**1801.05914v5** (submitted 18 Jan 2018, last revised 3 Jul 2021). **Abstract and
the `H_t`/`Φ` definitions read**; they agree verbatim with S1's eqs. (3)–(4). Taken from
it: the existence statement "the zeroes of `H_t` are all real precisely when `t ≥ Λ`",
and the result `Λ ≥ 0` — used **only** for the remark that the forward direction leaves no
margin at `t = 0`, not as an input to any estimate here.

### S3. C. M. Newman and W. Wu, *Constants of de Bruijn–Newman type in analytic number theory and statistical physics*
arXiv:**1901.06596v2**, 13 Mar 2019. **Body read** (PDF extracted to text locally). This
is the source for the *hypotheses*, which is why it was read rather than cited from the
abstract:

* **(10)** `F(−t) = conj F(t)` for all `t`;
* **(11)** `|F(t)| ≤ A exp(−|t|^{2+α})` for some `A, α > 0`;
* **Theorem 7** (attributed to de Bruijn [DB50]): *"Suppose that the function `F`
  satisfies (10), (11) and the zeros of the entire function (9) lie in the strip
  `|ℑz| ≤ Δ`. Then all the roots of the entire function `∫ F(t) e^{λt²/2} e^{izt} dt`
  (16) lie in the strip `|ℑz| ≤ max(Δ² − λ, 0)^{1/2}` (17)."*
* the framework-level setting is "`ρ` **a positive Borel measure** on the real line"
  (§1, eq. (1)), and for `ζ` the paper notes separately that *"It is not hard to see that
  `Φ(u) ≥ 0` for all `u`; this positivity is valid term by term … for `u ≥ 0` and then the
  even-ness of `Φ` implies positivity"*;
* `Λ_DN = Λ_DN(Φ) ∈ (−∞, 1/2]`, with `1/2 ∈ P_Φ` **due to de Bruijn [DB50] — Theorem 7**,
  and `P_Φ` bounded below due to Newman [Ne76] — Theorem 9.

**The point that matters for this audit, and the reason the body had to be read:**
**Theorem 7 does not require kernel positivity.** It needs (10), (11) and a bounded zero
strip. Positivity is a hypothesis of the Newman–Wu *framework* (a positive measure) and is
noted separately for the `ζ` kernel; it is not used by the theorem that supplies the
`Λ ≤ 1/2` bound and the strip-shrinking conclusion. That is what makes the theorem
applicable to the Davenport–Heilbronn function below.

### Normalization reconciliation, checked not assumed
S1 deforms with `e^{tu²}`; S3's Theorem 7 with `e^{λt²/2}`. Hence `λ = 2t`, and (17)
becomes `|ℑz| ≤ max(Δ² − 2t, 0)^{1/2}`: all zeros real once `t ≥ Δ²/2`. For `ζ` the
classical strip gives `Δ = 1`, so `t ≥ 1/2` — **de Bruijn's 1950 value, reproduced**,
which is the check that the `λ = 2t` bookkeeping is right. It also reproduces S1's
Theorem 3.2 with `y_0 = Δ`, `t_0 = 0`.

## Literature values restated (not used as inputs)

The note drawn from this audit restates the two known bounds for `ζ` in order to contrast
them with the DH bracket. Neither enters any estimate here. Both are recorded with the
place they were read, so that no value is carried from memory:

| value | source, as read |
|---|---|
| `Λ ≥ 0` | **S2** abstract, arXiv:1801.05914v5: the de Bruijn–Newman constant is non-negative |
| `Λ ≤ 0.22` | **S1** abstract, arXiv:1904.12438v2 (v2 of 4 Aug 2019, listed as to appear in *Research in the Mathematical Sciences*): *"we are able to obtain a new upper bound `Λ ≤ 0.22` unconditionally"* |
| `Λ ≤ 1/2` | de Bruijn's classical value, as **S3** states it: `1/2 ∈ P_Φ`, "due to De Bruijn [DB50] — see Theorem 7" |

**References resolved from S3's own reference list** (p. 21–22 of the local text
extraction), rather than from memory, since the note cites them:

* `[DB50]` N. G. de Bruijn, *The roots of trigonometric integrals*, Duke Math. J. **17**(3)
  (1950), 197–226. **Not read directly**; used only as the attribution S3 gives for
  Theorem 7, which was read.
* `[Ne76]` C. M. Newman, *Fourier transforms with only real zeros*, Proc. Amer. Math. Soc.
  **61**(2) (1976), 245–251. **Not read directly**; S3's attribution for `P_Φ` bounded
  below.
* `[CSV94]` G. Csordas, W. Smith and R. S. Varga, *Lehmer pairs of zeros, the
  de Bruijn–Newman constant `Λ`, and the Riemann Hypothesis*, Constructive Approximation
  **10**(1) (1994), 107–129. **Not read directly.** S3 describes its method as a zero
  *repulsion* argument on Lehmer pairs, giving `Λ_DN ≥ −4.379×10⁻⁶`. It is therefore
  **not** the background for the Rouché persistence bound of §6c′, which is a local
  argument on one disc and uses no repulsion phenomenon; the note cites it only for what
  the `ζ` lower-bound literature actually did, and says so.

## S4. The normalisation audit (2026-10-24), and two more primary sources read

A post-draft check of the note's upper bound asked a question the earlier audit had answered
only for itself: **in which `z`-variable is the strip half-width measured?** The threshold
`Δ²/2` is a statement about the variable the transform is taken in, and the literature uses
at least three such variables for `ζ`, so one theorem appears in print with three numbers.
Two further sources were obtained and read to settle it.

### S4a. C. M. Newman, *Fourier transforms with only real zeros*, Proc. Amer. Math. Soc. **61** (1976), no. 2, 245–251
**Read** (PDF, and p. 246 read from the page image because the OCR of that page is poor).

* **Theorem 1**, quoted: *"Suppose `f` is a real, even, integrable function on ℝ with
  `|f(t)| = O(exp(−|t|^γ))` as `|t| → +∞` for some `γ > 2`, such that all the zeros of the
  Fourier transform of `f` lie in the strip `|Im z| ≤ Δ`; then for `δ ≥ 0`, all the zeros of
  the Fourier transform of `exp(δt²)f(t)` lie in the strip
  `|Im z| ≤ [max(Δ² − 2δ, 0)]^{1/2}`."*
  **This is the theorem in the multiplier convention this audit uses** — `exp(δt²)`, threshold
  `δ = Δ²/2` — so it, and not the survey's restatement, is what the note now cites.
* **The attribution**, p. 246: *"due to de Bruijn [B, Theorem 13]"*, with `[B]` = the Duke
  1950 paper. **The original theorem number is 13**; `7` is Newman–Wu's own numbering.
* **The convention**, p. 246, quoted: *"all the zeros of `Ξ₀` lie in the critical strip
  `|Im z| ≤ ½`, also implies that (at least) for `b ≤ −⅛`, the zeros of `Ξ_b` are all real."*
  So **Newman's own `z` gives `1/8`**, because his `z = −i(s−½)` and `Δ = ½`.

### S4b. MathSciNet review MR0037351 of de Bruijn (1950)
**Read.** Confirms the strip result in outline and records de Bruijn's own normalisation,
`Ξ(2z) = ∫ φ(t)e^{izt} dt` with `φ(t) = Σ(2n⁴π²e^{9t/4} − 3n²πe^{5t/4})exp(−n²πe^t)` — a
**third** scaling, different again from both of the above.

### The four conventions, and which one this audit is in

| source | `z` in terms of `s` | multiplier | `Δ` for `ζ` | threshold for `ζ` |
|---|---|---|---|---|
| de Bruijn 1950 (per MR0037351) | `Ξ(2z)` | — | — | — |
| Newman 1976, Thm 1 + p. 246 | `−i(s−½)` | `e^{δt²}` | `1/2` | **`1/8`** |
| Newman–Wu, Thm 7 | — | `e^{λt²/2}` | — | threshold `λ = Δ²` |
| **Polymath15 / Rodgers–Tao — and this audit** | **`−2i(s−½)`** | **`e^{τu²}`** | **`1`** | **`1/2`** |

**This audit is in the last row**, and that is pinned by its own check, not asserted: the
same computation applied to `ζ` returns `τ ≥ 1/2`, the value Polymath15 and Rodgers–Tao
quote. **The DH number `1.602537483559822777787819` is `2Δ_s²` with
`Δ_s = 0.895136158235109721061359`**; in the Titchmarsh-type scaling of row 2 the same bound
would read `Δ_s²/2 = 0.40063437088995569`, exactly a quarter. **Polymath15's `Λ ≤ 0.22` is
in this audit's convention**, so it is comparable with `1/2` and not with `1/8`.

None of these sources disagree with each other. They are one theorem in several variables,
and the hazard is **inheriting a number across them**, which is what the note now avoids by
stating its convention and checking it.

### S4c. Reference detail verified, not typed
B. Rodgers and T. Tao, *The de Bruijn–Newman constant is non-negative*, **Forum of
Mathematics, Pi, vol. 8 (2020), article e6, doi:10.1017/fmp.2020.6** — confirmed by resolving
the DOI, having previously been cited here only by its arXiv identifier.

### S4d. Csordas–Smith–Varga: still not consulted
It is cited in the note **only** for the route the `ζ` lower-bound literature took, and that
description is taken from **S3 §2.6**, which was read. No result of CSV94 is used.

## Not used, deliberately

* No other recent RH announcement or preprint was surveyed; none is cited.
* The Polymath numerical upper bound on `Λ` is not used.
* `Λ ≥ 0` (S2) is used only as a remark, never as an inequality input.
* Nothing from the confluent-interpolation work in this repository was needed: no
  projection or directional estimate from it appears in the final argument.

## Reused from this repository (provenance recorded, nothing modified)

| reused | from | tier there |
|---|---|---|
| the window form `Q`, `dh.build_matrix`, `dh.G_at`, `dh.H_at`, `dh.mpf_hi` | `davenport-heilbronn-2026-10-23/src/dh.py` at `e2cf715` | — |
| `σ_0 = 1.395136158235109721061359`, the proved zero-free abscissa of `Λ_F` | that extension's `REPORT.md` §3 / `DERIVATION.md` §1 | **T1** |
| the explicit formula identifying the arithmetic side with `Σ_ρ H(ρ−½)` | its `DERIVATION.md` §3 | T1 derivation |
| its numerical agreement (arithmetic vs zero-side `5.65×10⁻¹¹`, vs contour `1.73×10⁻¹³`) | its Phase C, classification **PASS-N** | **T2-N** |
| `λ_min(Q_0^{DH}) = −0.679071346441480590986377586093…`, Arb ball radius `1.39×10⁻¹¹⁶` | its `certificates/ball_certificate.txt` | **T2-C** |
| Phase H: `λ_min` of the Euler-only truncation negative at 17 of 18 configurations | its `REPORT.md` §15 | T2-N |
| the rank-one channel reduction, `η = u^T B^{-1} u`, `η = 1 − det Q/det B` | `channel-compatibility-exploration-2026-10-24` | T2-N |

## Access notes

The Newman–Wu PDF had to be fetched and converted locally; an automated reader returned
only binary. The arXiv HTML of S1 was readable directly. S2's body was not needed beyond
the definitions and the existence statement, both of which agree with S1.
