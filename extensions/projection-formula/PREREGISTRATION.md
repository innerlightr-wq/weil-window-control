# Preregistration — projection-formula

Committed before any computation. Any later change is logged as a deviation in `REPORT.md`.

## Expectation going in

This is a **response law** — how the window form reacts to a moved zero — **not a positivity
mechanism**. It cannot bear on RH. The planted data are fabricated. Stated here so that a
confirmation is not later read as more than it is.

## The candidate formula

Even sector, planted moved zero at ordinate `γ`, displacement `δ`:

> `λ_min(Q_δ) ≈ λ_floor − 4·K(L,γ)·cos²θ_N · δ²`

where `K(L,γ) = ∫_{−L}^{L} u² sin²(γu) du` (paper Thm 8), `r_γ(u) ∝ u·sin(γu)` is the
representer of `f ↦ F′(γ)`, `N` is the near-null eigenspace of the unperturbed baseline form
`Q_0 = Q_ζ + 2F(γ)²` (paper eq. (3)), and `cos²θ_N = ‖Π_N r_γ‖² / ‖r_γ‖²`.

## Hypotheses

- **H1 (projection).** The measured sensitivity `C(δ) := (λ_min(Q_0) − λ_min(Q_δ))/δ²` equals
  `4K(L,γ)·cos²θ_{N_δ}` to within 10%, where `N_δ` is the span of the eigenvectors of `Q_0`
  with eigenvalue `≤ 4K·δ²` (fixed rule, fixed here).
- **H2 (staircase).** `C(δ)` is a staircase in `δ`: it steps up when `δ²` crosses an eigenvalue
  scale of `Q_0`, and is flat between consecutive scales.
- **H3 (parity).** The odd-sector floor sits above the even one *because* odd `f` has `F(0)=0`.
  Test: the floor of the even sector restricted to `{f even : ∫f = 0}` should be comparable to
  the odd-sector floor.
- **H4 (time–band floor, confirmation only).** `λ_floor` is the prolate/Landau–Widom floor:
  `log λ*(L)` tracks `1 − χ_k` for the relevant `k`, and Zhu's `2π²N(T*)/ln N(T*)`. Prior art;
  nothing new claimed.

## Predictions

1. **Function field (exact).** Moving one conjugate pair off the circle,
   `e^{±iθ_j} → e^{±a}e^{±iθ_j}` with functional-equation partners, changes the trace data by
   `Δt(n) = 2cos(nθ_j)(cosh(na) − 1) = a²n²cos(nθ_j) + O(a⁴)`.
   On `c ∈ ker T_R` the `a²` form reduces to a negative multiple of `|S_1(c)|²`,
   `S_k(c) = Σ_i i^k c_i e^{i i θ_j}`, i.e. minus the squared derivative of the polynomial at
   `e^{iθ_j}` — the discrete analogue of `−4F′(γ)²`.
2. For `R ≥ d = 2g`, degenerate perturbation theory gives
   `λ_min(T_R^{(a)}) = −a²·max_{c∈N,‖c‖=1}(const·|S_1(c)|²) + O(a⁴)`, the max expressible as a
   squared projection norm.
3. On the ζ side, `C(δ) → 4K(L,γ)` as `δ` grows past the near-null scales (`cos²θ_N → 1`), and
   `C(δ)` falls below `4K` when `δ` is small (`cos²θ_N < 1`).

## Kill criteria

- **K1.** The relative error of prediction 2 must `→ 0` like `a²` as `a → 0`, at every `R ≥ d`.
  If not, **STOP** and report: the formula is wrong even where everything is exact.
- **K2.** H1 must match the measured `C` within 10% at `≥ 4` of the 5 windows, at the `δ`
  values used in the paper's §4.3. Otherwise report H1 **failed**.
- **K3.** H2 plateau edges must align with eigenvalue scales of `Q_0` within a factor 3 in
  `δ²`. Otherwise report H2 **failed**.
- **K4.** H3 is supported if the constrained-even floor is within a factor 3 of the odd floor
  at every `L ∈ {0.5, 0.6, 0.8, 1.0}`. Otherwise report H3 **failed**.

## Stage gates

- Stage 1 is a **gate**: if K1 fails, stop before Stage 2.
- Stage 4 is **confirmation only** (prior art: Connes "Letter" Fig. 1; Zhu's Landau–Widom law).
- Tiers: **T1** exact/proved, **T2** numerical with stated precision and a precision-doubling
  check, **T3** interpretation. Regression values kept; retractions logged.

## Non-claims (fixed now)

Response law only. No bearing on RH or on Weil positivity itself. All planted zero data are
fabricated: nothing here is evidence about the actual zeros of ζ. Stage 4 claims nothing new.
