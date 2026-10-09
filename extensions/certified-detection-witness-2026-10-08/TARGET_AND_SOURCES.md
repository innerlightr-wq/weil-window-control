# Target and sources

Exploration, 2026-10-08. Read-only outside this directory. No manuscript edits, commits,
pushes, branch changes, resets, uploads, installs or background jobs.

## 0. Honesty about the order of work

Not preregistered, and not labelled as such. A **float prototype was run before this file was
written**, in order to decide whether building certified machinery was worthwhile. What it
showed, and what therefore was already known when the plan below was fixed:

- the nodal constraint costs very little derivative sensitivity (`K_eff/K ≈ 0.995`);
- the baseline cost of both witnesses is `O(0.5)`, so the naive threshold sits at `δ ≈ 0.86`,
  where the Taylor remainder is several times the signal.

So the *sign* of the outcome was known before the certified run. What was **not** known, and is
what the certified work decides, is (i) whether the geometric side can be enclosed at all on this
machine, (ii) the exact `K_eff` and `b_f` with error bars, and (iii) whether the **exact**
finite-displacement quartet (which has no Taylor remainder) changes the conclusion. Item (iii) is
the only route that could have produced a detection interval, and it is run below.

## 1. Workspace state as found

```
repo    weil-window-control (paths below are relative to the repository root)
remote  https://github.com/innerlightr-wq/weil-window-control
HEAD    cb562f7  "Add the RH partition-coordinate audit; verdict MISMATCH; branch closed"
branch  main (= origin/main), the only branch
status  one untracked directory, extensions/newtonian-weil-response-2026-10-08/
        (the completed Newtonian–Weil response diagnostic, from the previous run)
```

## 2. Named sources: located, or reported missing

| the brief names | what actually exists |
|---|---|
| `extensions/density-stiffness-audit/` | **DOES NOT EXIST.** No directory, no file, no string `density-stiffness` anywhere in the repository. The nearest completed work is `extensions/newtonian-weil-response-2026-10-08/` (the stiffness-weighted susceptibility comparison, verdict MISMATCH) and `extensions/second-moment-generalization-2026-10-07/`. Those are what "the completed density/stiffness audit" must refer to, and they are read and reused below. |
| "the corrected Theorem 7" | `paper/weil_window_control.tex` `\label{thm:zeta}`, printed as **Theorem 7**, "relative perturbation estimate; even sector". Corrected by Retraction 18. ✔ |
| "Section 4.3" | printed **§4.3 "How close the measured constant comes"**, `\label{sec:sat}`. ✔ |
| "Limitation 2" | `\label{sec:limits}` item 2: `thm:zetasharp` is conditional on all other zeros being on the line; `Q̃ ≥ 0` follows neither from `Q_ζ ≥ 0` nor from `Q_0 ≥ 0` nor from a zero-free strip. ✔ |
| "the 97.09%/2.91% split" | **NOT FOUND.** No such split is recorded anywhere. The only `2.91` in the repository is `3/2+√2 = 2.914214`, the maximum of the `θ_S` multiplier in `extensions/sonin-first-prime-defect-2026-10-07/`, which is a different quantity. Treated as not existing; nothing is built on it. |
| "a fitted Landau–Widom floor law" | `paper/weil_window_control.tex` line 577 and ref. [11] (Zhu, arXiv:2608.24827). Used **nowhere** below. |
| "the existing projected-derivative / near-null-overlap experiments" | `extensions/second-moment-generalization-2026-10-07/PROJECTION_CONTEXT.md` §4 (the unweighted projection response) and `extensions/newtonian-weil-response-2026-10-08/` (the susceptibility enclosure `S_v`, and the proof that the stiffness weighting is null at `L = 0.8`). ✔ |
| "whether the nodal construction has already been tested" | **It has not.** Searched for `F(gamma)=0`, nodal, node-at-ordinate constructions: the only node experiment on record is a parity/`F(0)=0` one, which the brief itself distinguishes. The `F(γ)=0` constraint is new here. |

Source hashes: `SOURCE_HASHES.txt`.

## 3. Conventions, verified against the source before use

All four verified in `src/zeta_window.py`, `src/stage3e_certify.py`, `notes/proofs.md` §B,
and `paper/weil_window_control.tex` §4.2:

```
F(z) = int_{-L}^{L} f(u) e^{izu} du ;  admissible f: real, even, supp f ⊂ [-L,L], ||f||_2 = 1
Q_zeta  = Pole + Arch - LogPi - Prime                      (build_matrix)
Q_0     = Q_zeta + 2 F(gamma)^2                            (paper eq. (3))
Q_delta - Q_0 = -4 delta^2 [F'(gamma)^2 + F(gamma)F''(gamma)] + R_4      (paper eq. (4))
|R_4| <= (16/3) L^5 delta^4 e^{2 L delta}                  (proofs.md B.3)
basis even_d: phi_k = cos(w_k x)/sqrt(L), w_k = (2k+1)pi/(2L);  Gram = I exactly
```

`Q_zeta` (actual arithmetic form), `Q_0` (augmented baseline), `Q_delta` (planted form) and
`Q̃ = Q_0 − 4F(γ)² = Q_zeta − 2F(γ)²` (remaining-zero form) are kept distinct throughout.

**Remainder coefficient verified.** `proofs.md` §B.3 gives `|c_n| ≤ 2L(2L)^n/n!`, so
`|R_4| ≤ 4·2L·Σ_{n≥4}(2Lδ)^n/n! ≤ 8L·(2Lδ)^4 e^{2Lδ}/24 = (16/3)L^5δ^4e^{2Lδ}`. ✔ Normalisation:
`‖f‖₂ = 1`; domain: all `δ ≥ 0`. Theorem 8 carries an extra `(2/5)L^5δ^4` from the completed
square — that is the *lower*-bound budget and is not used here.

**One extraction artifact caught, and NOT turned into a correction.** The PDF text layer renders
the closed form of `K(L,γ)` as `L³/3 − L²sin2γL/(2γ) + Lcos2γL/(2γ²) − sin2γL/(4γ³)`, which
disagrees with direct quadrature by 2%. The **TeX** (lines 372–375) reads
`K = L³/3 − [ L²sin2γL/(2γ) + Lcos2γL/(2γ²) − sin2γL/(4γ³) ]` — the outer bracket is dropped by
`pdftotext`. With the bracket the formula reproduces quadrature to `2.3×10^{-14}`. **The
manuscript is correct; the extraction was not.** No correction is proposed.

## 4. The environment constraint, and the route around Limitation 4

`mpmath` and `numpy` are **not installed** and installs are forbidden, so
`src/zeta_window.py` cannot be run. Its archimedean block goes through digamma / Lerch Φ /
Hurwitz ζ — which is exactly what **Limitation 4** of the manuscript records as *not*
interval-certifiable ("mpmath's interval context has no digamma and no Lerch transcendent, so the
geometric side could not be interval-evaluated").

So the archimedean term is re-expressed here as a series of **elementary** integrals:

```
Re psi(1/4 + i a/2) + euler = sum_{m>=0} [ 1/(m+1) - 2 s_m/(s_m^2+a^2) ],   s_m = 2m+1/2
Arch(f) = -euler C(0) + sum_{m>=0} [ C(0)/(m+1) - 2 int_0^{2L} e^{-s_m u} C(u) du ]
```

with `C` the symmetrised autocorrelation. Each integral is evaluated in **closed form** (no
quadrature error), and the tail is bounded rigorously by
`|Σ_{m≥M} t_m| ≤ (¾C(0) + ½ sup|C'|)/(M−¾)`, from `|2I_m − 2C(0)/s_m| ≤ 2sup|C'|/s_m²` (integration
by parts; `C(2L) = 0`) plus integral bounds on the two residual series. **No digamma and no Lerch
transcendent appear.** This is a route around Limitation 4, not a circumvention of it: the price
is the `O(1/M)` tail, which is reported with every number.

## 5. Validation of this pipeline, before any witness claim

1. `K(L,γ)` closed form vs direct quadrature: rel. dev. `2.3×10^{-14}` at three `(L,γ)`.
2. `F'_k`, `F''_k` closed forms vs finite differences: `5×10^{-11}` and `2×10^{-4}` (the latter
   finite-difference-limited).
3. **The archimedean term, two independent representations.** For `φ_0` at `L = 0.8`:
   exact m-sum `−0.98589817` (tail bound `1.18×10^{-5}`) versus the spectral route
   `(1/2π)∫|F|²Reψ(1/4+ir/2)dr = −0.98590193`, computed with an independently implemented
   digamma validated to `10^{-12}` on `ψ(1)`, `ψ(1/4)`, `ψ(1/2)`. **Agreement `3.8×10^{-6}`,
   inside the bound.**
4. **Against the repository's recorded spectrum.** `L = 0.8`, `even_d`, `N = 4`:
   `λ_2` runs `6.46×10^{-6} → 3.09×10^{-6} → 2.786×10^{-6}` as `M` runs `2×10^5 → 2×10^6 → 10^7`,
   converging on the recorded `2.711086672×10^{-6}`; and `λ_min` decays like the truncation bound,
   `3.75×10^{-6} → 3.75×10^{-7} → 7.52×10^{-8}`, toward the recorded `1.7339346963×10^{-10}`.
   **`λ_min` itself is below this pipeline's error bar and is NOT resolved here** — which does not
   matter, because the target quantity `Q_ζ(f_witness)` is `O(0.5)`.

## 6. The frozen plan

**Target.** An UPPER bound `Q_δ(f) ≤ b_upper − a_lower δ² + E_upper(δ) < 0` for **one** explicit
admissible `f`, with `b_f ≤ b_upper`, `a_f ≥ a_lower > 0`, every term evaluated at that same `f`.

**Witnesses, both on `E = span{φ_0,…,φ_{N−1}}` (even_d), an explicitly frozen coefficient span —
no spectral projector is certified anywhere.**

- **A** — `P_E r` normalised, `r(u) = u sin(γu)`: the projection into `E` of the paper's own
  extremiser for `K(L,γ)`, i.e. the strongest existing explicit witness, with its **full** `FF''`
  term retained.
- **B** — the nodal witness: `w = P_E r − (⟨P_E r, P_E c⟩/‖P_E c‖²)P_E c`, `f_node = w/‖w‖`, which
  satisfies `F(γ) = 0` exactly in exact arithmetic.

**Ordinates.** `γ = 14` **exactly** (an exact chosen ordinate — the certified case) and
`γ = γ₁ = 14.134725141734693790…`, the documented planted ordinate, **not rigorously enclosed
here** and so labelled discovery-only. The two are never interchanged.

**Window and sizes.** `L = 0.8` (where Zhu's certificate reaches); `N ∈ {4, 8}` frozen, with
`N = 16` reported as a clearly-labelled trend row only.

**Decision rule, fixed now.** The certificate is the **exact** finite-displacement quartet, which
has no Taylor remainder:
`Q_δ(x) = Q_ζ(x) − 2⟨x,c⟩² + 4⟨x,a_δ⟩² − 4⟨x,b_δ⟩²`, with the positive term `+4⟨x,a_δ⟩²` retained,
not dropped. `f` is a certified negative witness iff `Q_ζ(x) + tail + Δ_exact(δ) < 0` for some
`δ`, scanned on `δ ∈ (0, 2]`. The Taylor budget `B_L(δ)` is reported alongside for comparison but
is not what the verdict rests on.

**Not repeated:** the density audit, the Newtonian analogy, the stiffness-weighting comparison,
the quasi-RH transfer. **Not used:** the Landau–Widom fit, the `97.09/2.91` split, any fitted
coefficient.
