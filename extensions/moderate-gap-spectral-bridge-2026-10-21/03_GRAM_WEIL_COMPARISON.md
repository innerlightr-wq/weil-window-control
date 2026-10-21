# 03 — The two Gram frameworks: exact dictionary and exact obstruction

Answers the brief's §4, questions 1–8, and Final Question 1.

## 0. Summary

There **is** an exact dictionary, and it is not mine: **it is Alpöge–Furman's own
stated construction.** Their abstract says the Riemann hypothesis
> "…classically needed to read the zero side as a positive sum over real ordinates,
> is replaced by a rank–trace inequality applied to **a finite compression of Weil's
> Hermitian form**, with Sylvester's law of inertia handling off-line pairs."

So ARTZ's `G` *is* a finite compression of the zero side of Weil's Hermitian form.
That settles the existence of the correspondence as **class II**, not class IV.

What is new here is **quantitative**: exactly which part of their matrix carries the
horizontal displacement, and that its leading coefficient is the same constant as in
my published `δ²–L³` law. And what is new on the negative side is a precise list of
the four structures the correspondence does **not** preserve — each of which is
individually fatal to a transfer.

## 1. The dictionary

| | ARTZ / Lamzouri / AF | my finite-window program |
|---|---|---|
| ambient space | `L²((−1/2,1/2), η_τ² du)`, real-symmetric sector | `L²([−L,L])`, real even sector |
| conjugate variable | `z = i(ρ−1/2)·log T/(2π)` | `γ` (unnormalised ordinate) |
| scale map | `u_mine = u_paper · log T` | — |
| effective half-window | `L_eff = (log T)/2` | `L` |
| test object | the *window* `η_τ`, fixed | the *witness* `f`, optimised |
| form | `G = Σ_{z∈Z_T} v_{f_z} v_{f_z}^T`, real symmetric, **indefinite** | `Q_ζ(f) = (zero side) − (prime side) − (arch)`, **positivity ⟺ RH** |
| moment order | **second** (`‖G‖²_F = Σ_{z,s}K(z−s)²` = pair correlation) | **first** (one explicit formula) |
| arithmetic input | BGSTB unconditional pair correlation (no prime sums appear) | an explicit prime sum (595 877 prime powers to `e^16` in the certified witness) |
| how an off-line zero shows | signature `(1,1)` block, inertia count | sign change of a quadratic form |
| displacement enters through | `κ_ρ = ∫η_τ²(u)cosh(2δ_ρ u log T)du` only | `F(γ ± iδ)`, i.e. `cosh(δ·)` against the window |

**The one genuinely new identification.** `Im z_ρ = δ_ρ · log T/(2π)`, and the
displacement appears *only* through

```
    kappa_rho = K_tau(2 i Im z_rho) = int eta_tau^2(u) cosh(2 delta_rho u log T) du
              = 1 + 2 delta_rho^2 V_mu + O(delta_rho^4) ,
    V_mu = centred second moment of the window measure in the unnormalised
           frequency variable v = u log T  (mass 1):  V_mu = (log T)^2 int u^2 eta_tau^2.
```
`[IV, certified: verify.py D4–D7, E5, I1–I3]`

My published law (`doi:10.5281/zenodo.23250932`, §"The δ²–L³ law", remark after
Thm `ff`) states that its `δ²` coefficient "is **twice the centred second moment of
the window**; the centring falls out of the optimisation, not from a convention."
That is verbatim the coefficient above. With a window of mass `~L^s` the coefficient
scales as `L^{s+2}`: `s = 0` (ARTZ's probability window) gives `L²`; `s = 1` (my
unnormalised window) gives the published `L³`. Checked exactly for the `η² ≡ 1`
window: `κ − 1 = (2πy)²/6` and `2δ²V_μ = δ²(log T)²/6` agree identically.
`[IV, certified I1, I2, I5, I6]`

**So the two frameworks share one constant, and it is the δ² sensitivity constant.**

## 2. The eight questions of the brief, §4

1. **Are the spaces naturally isomorphic after rescaling?** Yes, by the dilation
   `u ↦ u/log T`, which maps `(−1/2,1/2)` to `(−(log T)/2, (log T)/2)`. The
   isomorphism is isometric up to the Jacobian and maps real-even functions to
   real-even functions. `[I]`
2. **Are the kernels related by an explicit transformation?** Yes:
   `K_τ(ξ) = \hat{η_τ²}(ξ)` and my window's autocorrelation kernel are the same
   object in different variables. The *paper's* kernel is a probability density's
   transform (`K(0)=1`); mine is unnormalised. `[I]`
3. **Does the transformation preserve positivity?** **No, and the failure is
   structural, not technical.** ARTZ's `G` is deliberately indefinite:
   `n_-(G) ≤ n_-(Q) = k` counts the off-line pairs, and that indefiniteness is the
   *mechanism*. My `Q_ζ` is a form whose positivity on all admissible `f` is
   *equivalent to RH*. One matrix is a counter, the other a certificate. `[IV, proved]`
4. **Does it preserve the signed contributions of off-line zeros?** Only in sign,
   not in magnitude. Under the dictionary the off-line pair's `(1,1)` block has
   eigenvalues `m(1+κ), −m(κ−1)`, whose **trace `2m` is independent of `δ`**
   (`verify.py` E2). My form's response to displacement is `O(δ²)` with a
   non-vanishing coefficient. So the dictionary maps a `δ`-sensitive quantity onto a
   `δ`-blind one, unless one keeps the Frobenius norm — §04. `[IV, certified]`
5. **Does it preserve eigenvalue lower bounds?** No. ARTZ never use a lower bound on
   `λ_min`; they use `tr`, `‖·‖_F` and the inertia. My program's detector strength
   `κ(Γ)` is a *supremum over the admissible class*, i.e. an upper bound on what any
   witness can do. The two are not dual. `[III + IV]`
6. **Does it preserve the horizontal displacement parameter?** Yes, but only after
   the rescaling `a = δ log T`, and `a` is the only dimensionless version. BGSTB's
   own positivity hypothesis is `|β−1/2| ≤ c/log T`, i.e. **`a = O(1)`** — see §08.
   `[II + IV]`
7. **Does it preserve the compact-support constraints?** Only with the scale factor.
   ARTZ's support is `(−1/2,1/2)` in `u`, i.e. `(−(log T)/2,(log T)/2)` in my
   variable — the support **grows with `T`**. My certified witnesses live at
   `L ∈ {4/5, 4, 8}`, fixed. A transfer that needs `L ≍ log T` is outside every
   certificate I hold. `[IV, proved; see 08 §3]`
8. **Does it preserve the real-vs-complex distinction?** This is the subtlest point
   and it is where a naive transfer breaks. Two different pairings coexist:
   ```
       Hermitian:  <f_z, f_s>        = K(z - sbar)      (PSD Gram)
       bilinear:   v_{f_z} . v_{f_s} = K(z - s)         (what G is built from)
   ```
   They coincide on the critical line and differ off it. ARTZ's `‖G‖²_F = Σ K(z−s)²`
   uses the **bilinear** one; the identity behind it is
   `Σ_i⟨f_z,ψ_i⟩⟨f_s,ψ_i⟩ = ⟨f_z, f_{s̄}⟩`, valid because every basis function
   satisfies `conj(φ(u)) = φ(−u)`. Conflating the two would put the symmetric-diagonal
   terms `K(z−z̄)² = κ²` in the wrong place. The two *total* sums agree only because
   `Z_T = Z̄_T`. `[IV, proved; verified numerically F1–F3]`

## 3. The four non-preserved structures, as a table

| structure | ARTZ side | my side | preserved? |
|---|---|---|---|
| positivity | indefinite by design, `n_- = k` | positivity ⟺ RH | **no** (Q3) |
| moment order | second (pair correlation) | first (explicit formula) | **no** |
| arithmetic input | BGSTB; no prime sums | explicit prime sum | **no** |
| support scale | `L ≍ log T` | `L` fixed | **no** (Q7) |

The second row is the deepest. **Squaring destroys the sign information that makes a
Weil form a detector.** In ARTZ's framework an off-line zero of depth `δ` costs
`(κ−1)² = 4V_μ²δ⁴ + O(δ⁶)` — a *fourth*-order quantity — because the budgeted object
is the square of a `δ²` response. In my framework the Weil form's *response* to
displacement is `δ²·L^{s+2}`, and the *discrimination* strength against the best
on-line mimic was measured independently at `M(δ) = δ⁴·D₂`
(`~/scratch/weil-offline-detect-2026-10-17-083545`, certified
`D₂ = 2911.050006543213352847841125028374273138 ± 1e-37`).

> **Cross-framework agreement (class IV).** Both frameworks give exponent **4** for
> the budget/discrimination cost of an off-line pair, and in both the exponent is
> `2 × 2`: the displacement enters through `cosh(δ·)`, which is even, so the
> *response* is `O(δ²)`; the budgeted quantity is the square of that response. The
> constants differ (they are the respective windows' second moments) but the
> exponent is framework-independent. This is an independent confirmation of the
> `δ⁴` law from a construction that shares none of its machinery.

## 4. Required outcome

**Dictionary: delivered** (§1–§2), with the existence of the correspondence correctly
attributed to Alpöge–Furman [AF, abstract and §3] and the quantitative content
(§1, the `κ` identity and the shared `δ²` constant) new here.

**Obstruction: delivered** (§3). The correspondence is an isomorphism of Hilbert
spaces and of kernels, but **not** of positivity, moment order, arithmetic input, or
support scale. Verdict for Final Question 1 / Candidate T1: **D — exact structural
correspondence, no independent spectral inequality from the dictionary alone.**
