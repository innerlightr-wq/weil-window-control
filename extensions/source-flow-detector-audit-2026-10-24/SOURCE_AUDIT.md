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
