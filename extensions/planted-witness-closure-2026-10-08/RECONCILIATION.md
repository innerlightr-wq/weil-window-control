# Reconciliation: confirmed, corrected, still unchecked

What the preceding run (`extensions/certified-detection-witness-2026-10-08/`) asserted,
against what this run can prove. Source hashes of everything read: `SOURCE_HASHES.txt`.

## 1. Confirmed

| assertion | status |
|---|---|
| The nodal projection lemma (`F(γ) = 0`, `\|F'(γ)\|² = K_eff`, constrained extremality) | **confirmed**, elementary linear algebra as it was labelled |
| The exact finite-displacement quartet `4 Re F(γ+iδ)² = 4[A_δ² − B_δ²]`, with `F(γ+iδ) = A_δ − iB_δ` for even `f` | **confirmed**, re-derived independently in `CERTIFICATE_THEOREM.md` §2 |
| The archimedean reduction `Re ψ(1/4+ia/2) + γ_E = Σ_m[1/(m+1) − 2s_m/(s_m²+a²)]` | **confirmed** |
| The archimedean tail bound `[¾C(0) + ½sup\|C'\|]/(M−¾)`, with its hypotheses and indexing | **audited and confirmed correct**; re-derived in `CERTIFICATE_THEOREM.md` §5.2. Its hypothesis `C(2L) = 0` is now *verified* (`\|C(2L)\| < 9e-78`) rather than assumed |
| `Q_ζ(f_node) ≈ 0.00517` | **confirmed and now enclosed**: `[0.0051242393, 0.0052229400]`. The previous float value `0.0051682432` and the independent spectral value `0.0051678677` both lie inside |
| `Q_ζ(f) > 0`, consistent with Weil positivity | **confirmed** |
| `F'(14)`, `F''(14)`, `a_f` | **confirmed and enclosed**: `−0.17626230294801`, `+0.11323208261102521`, `0.12427359776233127` |
| The `K(L,γ)` PDF-extraction artifact (bracket dropped, two signs flipped) — manuscript correct, extraction wrong | **confirmed**, no correction proposed |
| `extensions/density-stiffness-audit/` does not exist; the `97.09%/2.91%` figures appear nowhere | **confirmed**; neither is used |
| A negative value of the planted functional is attainable at this window with this witness | **confirmed, and now proved** — see `CERTIFICATE_THEOREM.md` §3 |

## 2. Corrected

### 2.1 The Taylor-regime claim was wrong *(the brief's correction A)*

The previous run stated that `R_bound ≤ 2K_eff δ²` "fails at every `δ` in the strip". **That
is false**, and the brief's reasoning is right: `R/(2K_effδ²) = (8/3)L⁵δ²e^{2Lδ}/K_eff → 0`
as `δ → 0`. The error was testing only `δ ≥ 0.19`. Enclosed values:

| `δ` | `R/(2K_eff δ²) ≤` | condition holds? |
|---|---|---|
| `1/10` | `0.330055` | **yes** |
| `3/20` | `0.804476` | **yes** |
| `166/1000` | `1.010798` | no |
| `167/1000` | `1.024652` | no |
| `1/5` | `1.549295` | no |
| `2/5` | `8.534307` | no |

So the condition **holds for `δ ≲ 0.166`** and fails above. **The actual obstruction is
non-overlap**: the range where the Taylor remainder is dominated by the signal stops at
`≈ 0.166`, while the range where the signal beats the baseline begins near `≈ 0.204`. The
useful distinction survives in corrected form — *the coarse Taylor test fails to certify
this witness even though its remainder is small at sufficiently small displacement*.

### 2.2 The basis-size claim was overstated *(the brief's correction B)*

The previous run wrote that `N = 8` and `N = 16` "do NOT detect". The bases are **nested**
and describe the same form, so the four-mode witness is available in the eight- and
sixteen-mode spaces **by zero padding**, with the *same* value of every functional; an
optimised minimum over a larger space cannot increase. **Zero-padding invariance is
structural here, not merely numerical**: the computation in `verify_certificate.py` is
expressed in terms of the four nonzero coefficients and never references `N`, so padding
changes nothing identically.

Corrected statement: the worse `N = 8` and `N = 16` numbers are outcomes of **that
particular projected-witness selection rule** (`P_E r` with the node projected out, `E`
growing), not evidence that detection is impossible in larger spaces. Statements about
witness A are restricted to the tested cases; no all-`N` claim is made.

### 2.3 The proposed uniform theorem does not do what was suggested *(correction C)*

`Q_rest ≤ θ·4K_eff` with `θ < 1/4` does **not** yield detection at arbitrarily small `δ`:
the quartet contributes `≈ −4K_eff δ²` at leading order, so beating a *fixed* positive `θ`
requires `δ² > θ` (plus the error allowance). The previous phrasing invited the opposite
reading and is withdrawn.

Also, for fixed `L` and a fixed **unmodulated** span `E_N`, `⟨φ_j, u sin(γu)⟩ → 0` as
`|γ| → ∞` by Fourier decay of the finitely many integrable functions `u φ_j(u)`, hence
`K_eff ≤ ‖P_{E_N} r_γ‖² → 0`. Enclosed:

| `γ` | `‖P_E r_γ‖²` |
|---|---|
| `14` | `3.601682e-02` |
| `70` | `1.100552e-05` |
| `140` | `2.128466e-06` |
| `280` | `1.129069e-07` |

(The `γ = 14` entry is large because `w_3 = 13.744…` nearly resonates with `14`.) This is
distinct from the continuum norm `K ∼ L³/3` and it **obstructs inferring an all-height
result from a fixed low-frequency basis**. A frequency-adapted span changes the premise and
would need a new baseline analysis; that is not begun here. None of this rules out other
detection methods, and none of it settles any actual-zero question.

### 2.4 The error allowance was understated, and the margin factor with it

The previous run quoted `±3.76e-6` for `Q_ζ(f_node)` and a `~6800×` margin-to-error factor.
Two corrections:

- that allowance was an **entry-wise aggregation** `Σ|x_i| E_ij |x_j|` of per-matrix-entry
  tail bounds at `M = 2×10^6`; the direct bound for the *assembled* function at
  `M = 2^17 − 1` is `4.935e-5`. The two are mutually consistent (the bound is linear in
  `1/M`; at `M = 2×10^6` the present formula gives `3.2e-6`), so the difference is the
  truncation `M`, chosen here as `2^17 − 1` so that `ln(M+1) = 17 log 2` exactly;
- the `6800×` figure was computed at `δ = 1/2` and is **not** a uniform margin. The honest,
  separated figures are **294×** at the certified point `δ = 2/5` and **50×** interval-wide
  (worst at `δ = 1/4`).

### 2.5 The object was not the one intended

The previous script evaluated the **count-preserving replacement** functional
`Q_rep,δ = Q_ζ − 2F(14)² + 4[A_δ² − B_δ²]`, whose `δ → 0` baseline is `Q_ζ + 2F(14)²`, while
the intended object is the **added-quartet** functional `Q_add,δ = Q_ζ + 4[A_δ² − B_δ²]`,
baseline `Q_ζ + 4F(14)²`. They differ by `2F(14)²`. For the frozen rational vector

```
F(14) = 4.1293580707791e-13 ,   Q_add - Q_rep = 2F(14)^2 = 3.4103196e-25 ,
```

i.e. `1.4e-23` of the certified margin. The previous conclusion therefore survives the
change of object, but the functionals are different and the previous report conflated them.
The certificate here is stated for `Q_add` and the difference is bounded, not ignored.

### 2.6 The onset decimal is not a certified threshold

`0.204201571681` was a float bisection of a float model. **No root enclosure has been
performed, here or previously**, so it is not reported as a certified threshold. What is
certified is negativity on `[1/4, 49/100]`, which does not pin the onset.

## 3. Still unchecked

1. **The explicit formula itself** is a classical cited input, not reproved
   (`CERTIFICATE_THEOREM.md` §4.1). Its applicability to this `f` is argued (continuous,
   compactly supported, `|F|² = O(r^{-4})`) but the identity is Weil's.
2. **No root enclosure** for the first-negative displacement.
3. **`γ₁` is still not enclosed.** This run deliberately uses `γ = 14` exactly and makes no
   statement at `γ₁`.
4. **Nothing at other windows or ordinates.** `L = 4/5`, `γ = 14`, one four-dimensional span.
5. **The interval was not optimised.** `[1/4, 49/100]` is the brief's proposed target,
   certified as given; the true negative set is larger and is not determined.
6. **Normalisation.** The certificate is for the unnormalised `f`; `‖f‖²` is enclosed and
   positive, so the normalised quotient has the same sign, but no normalised margin is
   quoted.
