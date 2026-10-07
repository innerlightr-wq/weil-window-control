# Projection-formula test — report

Does a projection law describe the window form's **response** to a moved zero?

> **Verdict: PARTIAL.** The projection law is confirmed as a response law — K1, K2 and K3 all
> pass, and it reconciles the paper's §4.3 saturation ratios. The parity hypothesis H3 is
> **false**: K4 fails by two to three orders of magnitude.

Tiers: **T1** proved/exact, **T2** computed at stated precision, **T3** interpretation.
Preregistration in [`PREREGISTRATION.md`](PREREGISTRATION.md), committed before any computation.

---

## Stage 1 — function field, exact (GATE): **K1 PASSES**

### (a) The `a²` form (T1, derived)

With `t(n) = p_n q^{−n/2} = Σ_{j=1..g} 2cos(nθ_j)` and `T_R[i,k] = t(|i−k|)`:

> `cᵀ T_R c = Σ_j 2 |P_c(e^{iθ_j})|²`, `P_c(z) = Σ_i c_i z^i`.  **(T1)**

Verified numerically at dps 60: relative error `2.3e−61`, `0`, `2.7e−61`, `4.8e−61` on the four
honest curves. Hence `ker T_R = {c : P_c(e^{iθ_j}) = 0 ∀j}`, of dimension `R+1−2g` for `R ≥ 2g`.

Moving one conjugate pair off the circle with its functional-equation partners gives
`Δt(n) = 2cos(nθ_j)(cosh(na) − 1) = a²n²cos(nθ_j) + O(a⁴)`. Writing `S_k(c) = Σ_i i^k c_i e^{iiθ_j}`
and `m = i−k`, expanding `(i−k)² = i² − 2ik + k²`:

> `Q₂(c) := Σ_{i,k} c_i c_k (i−k)² cos((i−k)θ_j) = 2 Re( S₂(c) conj(S₀(c)) ) − 2|S₁(c)|²`  **(T1)**

On `c ∈ ker T_R` we have `S₀(c) = P_c(e^{iθ_j}) = 0`, so the first term vanishes and

> `Q₂(c) = −2|S₁(c)|² = −2 |P_c′(e^{iθ_j})|²`  **(T1)**

since `S₁(c) = e^{iθ_j} P_c′(e^{iθ_j})`. This is minus the squared derivative of the polynomial at
`e^{iθ_j}` — the discrete analogue of `−4F′(γ)²`, as predicted. The constant is exactly **2**.

### (b) Degenerate perturbation theory (T1)

`T_R` vanishes on `N = ker T_R` and the perturbation is `O(a²)` against an `O(1)` gap, so to
leading order `λ_min` is the least eigenvalue of the perturbation compressed to `N`:

> `λ_min(T_R^{(a)}) = −2a² · μ_max + O(a⁴)`, `μ_max = max_{c∈N, ‖c‖=1} |S₁(c)|²`.  **(T1)**

`|S₁(c)|² = ⟨c,u⟩² + ⟨c,v⟩²` with `u_i = i cos(iθ_j)`, `v_i = i sin(iθ_j)`, so

> `μ_max = λ_max( Π_N (uuᵀ + vvᵀ) Π_N )`, i.e. the squared operator norm of `Π_N` on the
> two-dimensional representer space `span{u,v}` — the `cos²θ` form, with `‖u‖²+‖v‖² = Σ_{i≤R} i²`.

The representer space is **two**-dimensional here, against one (`r_γ ∝ u sin γu`) in the
continuous even sector, because real `c` sees both the cosine and sine parts. **(T3)**

### (c) Verification against exact `λ_min` (T2, dps 60 and 120)

| | |
|---|---|
| rows (`R = d … d+10`, four curves, `a = 10⁻²,10⁻³,10⁻⁴`) | 132 |
| `rel.err(a=10⁻²)/rel.err(10⁻³)` | **100.00 – 100.08** at all 44 `(curve,R)` pairs |
| `rel.err(10⁻³)/rel.err(10⁻⁴)` | **100.00** at all 44 |
| precision doubling, dps 60 → 120 | max relative change **0.000e+00** |

> **K1 PASSES.** The relative error goes to zero exactly like `a²` at every `R ≥ d`, on every
> curve. Regression values in [`data/stage1_ff.json`](data/stage1_ff.json).

### `R < d` (no kernel)

`T_R ≻ 0` there, the perturbation is **non-degenerate**, and first-order PT with the rank-one
projection onto the ground state gives relative error `O(a⁴)` — measured
`1.04e−10 → 1.04e−14 → 1.04e−18`, i.e. a factor `10⁴` per decade, better than the kernel case.
A threshold near-null space therefore does work, degenerating to `N = span(v_min)`. As `R → d`,
`|S₀(v_min)| → 0` (`1.237 → 0.069` for G3F3), a continuous transition into the kernel regime. **(T2)**

---

## Stage 2 — ζ side: **K2 PASSES, K3 PASSES**

Setup verified: `move_block` at `δ=0` equals `+2pp^T`, so `Q_0 = Q_ζ + 2F(γ)²` exactly, and
`Q_δ = build_matrix + move_block(δ)` is continuous at `δ=0`. **(T1)**

`‖r_γ‖²/K(L,γ₁)` in the finite basis = `0.934, 0.930, 0.987, 0.978, 1.000` at `L = 0.8…2.0` —
the basis does not represent `r_γ` exactly, which is part of the residual H1 error. **(T2)**

### K2, at the paper's §4.3 value `δ = 0.02`

| L | C (ours) | C (paper §4.3) | `4K cos²θ_N` | `cos²θ_N` | rel H1 | rel compressed |
|---|---|---|---|---|---|---|
| 0.8 | 0.3107976 | 0.3107976 | 0.3287235 | 0.44308 | **5.77 %** | 1.86 % |
| 1.0 | 0.8530343 | 0.8530343 | 0.893829 | 0.66572 | **4.78 %** | 2.57 % |
| 1.3 | 2.498698 | 2.498698 | 2.73361 | 0.87734 | **9.40 %** | 0.19 % |
| 1.6 | 4.921791 | 4.921791 | 5.044572 | 0.98662 | **2.50 %** | 0.02 % |
| 2.0 | 10.08283 | 10.08283 | 10.07125 | 0.94546 | **0.11 %** | 0.73 % |

Reconciliation with §4.3: `max |C_ours − C_paper|/C_paper = 1.68e−10`.

> **K2 PASSES: within 10 % at 5 of 5 windows.**

### Reconciling the §4.3 ratios

| L | 0.8 | 1.0 | 1.3 | 1.6 | 2.0 |
|---|---|---|---|---|---|
| `cos²θ_N` | 0.443 | 0.666 | 0.877 | 0.987 | 0.945 |
| §4.3 `C/4K` | 0.419 | 0.635 | 0.802 | 0.963 | 0.947 |

> **The saturation ratio of §4.3 *is* the projection fraction** (`figures/cos2_vs_ratio.png`).
> The paper's "the re-optimised minimiser reaches 96.3 % at `L=1.6`" is the statement that at
> `L=1.6` the representer `r_γ` lies almost entirely inside the near-null space. **(T3)**

### K3, the staircase

`C(δ)` plateaus and steps; grey lines in `figures/staircase.png` mark `δ_k = √(λ_k/4K)`.

> **K3 PASSES: 11 of 12 `C`-steps align with an eigenvalue scale within a factor 3**
> (L=0.8: 4/4, L=1.0: 4/5, L=1.3: 2/2, L=1.6: 1/1, L=2.0: 0 steps in range).

### Where the scalar formula fails, and why

At isolated `δ` the fixed threshold rule `λ ≤ 4Kδ²` misassigns `N_δ` and the scalar prediction
is badly wrong — `L=0.8, δ=10⁻⁵`: `C = 1.39e−4` against a prediction of `0.018` (rel 128). At the
*same* points the **compressed** prediction `λ_min(Λ_N + δ²Π_N P Π_N)` is accurate to `3.7e−5`.

> **(T3)** The mechanism — degenerate perturbation theory on the near-null space — is right. The
> scalar `cos²θ` is a coarse proxy for it, and its failures are failures of the preregistered
> threshold rule, not of the mechanism. The compressed form is the one to quote.

---

## Stage 3 — parity (H3): **K4 FAILS**

In the EVEN basis `∫cos(kπx/L)dx = 0` for every `k ≥ 1` and `2L` for `k = 0`, so the constraint
`F(0) = ∫f = 0` is exactly "delete the constant mode" **(T1)**; confirmed numerically,
`max_{k≥1}|F_k(0)| ≤ 2.2e−81`.

| L | even floor | even with `∫f=0` | odd floor | constrained/odd |
|---|---|---|---|---|
| 0.5 | 1.109e−6 | 1.256e−2 | 2.386e−4 | **52.6** |
| 0.6 | 2.330e−9 | 8.478e−5 | 7.077e−7 | **119.8** |
| 0.8 | 3.199e−17 | 1.246e−11 | 2.586e−14 | **481.9** |
| 1.0 | 5.270e−27 | 1.680e−21 | 1.882e−24 | **892.6** |

> **K4 FAILS at every `L`** — 53× to 893×, against a tolerance of 3×. **H3 is false.**

**(T3)** The failure direction is informative: imposing `F(0)=0` on even functions raises the
floor *far above* the odd floor, growing with `L`. So `F(0)=0` is not what sets the odd floor —
odd functions reach much smaller Rayleigh quotients than `F(0)=0`-constrained even ones. Parity
carries near-null structure that the single linear constraint does not capture.

---

## Stage 4 — time–band floor (H4, confirmation only)

Window dictionary: Connes's multiplicative window `[λ⁻¹, λ]` has additive log-length `2 log λ`;
ours is `[−L,L]`, additive length `2L`; hence **`log λ = L`**, no factor-of-2 discrepancy. **(T1)**

| L | 0.8 | 1.0 | 1.3 | 1.6 | 2.0 |
|---|---|---|---|---|---|
| `−ln λ*(L)` | 38.21 | 61.04 | 89.67 | 115.11 | 154.17 |
| Zhu `2π²N(T*)/ln N(T*)` | 56.36 | 77.24 | 142.29 | 270.35 | 636.85 |
| ratio | 0.678 | 0.790 | 0.630 | 0.426 | 0.242 |

**T2 confirmation only.** Both sides grow super-exponentially in `L`, but the ratio **decreases
steadily** (0.68 → 0.24), so at these `L` and `N` the measured floor is not yet tracking the
asymptotic law closely. Two reasons, neither in tension with the law: `λ*(L)` here is a finite-`N`
variational *upper* bound, so `−ln λ*` is an *under*estimate that worsens as `L` grows at fixed
`N`; and the law is asymptotic in `N(T*)`, which is only `≈ 7…117` across this range.
**Nothing new claimed, and this is not offered as a check of the law.**

Incidentally `λ_min(Q_0)` and `λ_min(Q_ζ)` agree to all 8 printed digits at every `L`: the rank-one
bump `2F(γ₁)²` does not move the floor, because the window minimiser already nearly annihilates
`F(γ₁)`. **(T2)**

**Scope reduction (deviation 4):** the prolate quantities `1 − χ_k` were **not** computed. Doing
so at the `10⁻⁵⁰` scale needs dedicated prolate code; the Zhu-law comparison is reported instead,
and the identification "near-null eigenspace = prolate projection `Π(λ,k)`" is cited, not verified
here.

---

## Stage 5 — combined formula

`λ_min(Q_δ) ≈ λ_floor − 4K cos²θ_N δ²`, 60 `(L,δ)` cells:

| regime | rows | median relative error |
|---|---|---|
| `Cδ² ≫ λ*` | 51 | **3.02e−2** |
| `Cδ² ≪ λ*` | 7 | 4.80e−6 |

> **It holds to ≈3 % wherever the response term dominates the floor**, which is the whole regime
> the paper's detection criterion `C(L)δ² > λ*(L)` lives in. In the opposite regime it is accurate
> only trivially: there `λ_min ≈ λ_floor` and the `δ²` term is negligible, so the formula is
> testing nothing. **(T2/T3)**

---

## Stage 6 — prior art (brief)

| source | the specific statement `C = 4K cos²θ_N` | the parity explanation |
|---|---|---|
| Bombieri 2000 §13 | not found | not found |
| CCM arXiv:2511.22755 | not found; **but** "the space of eigenvectors of the `k` lowest eigenvalues of `QW_λ` corresponds to the prolate projection `Π(λ,k)`" is there — the near-null space *is* prior art | not found |
| Zhu arXiv:2608.24827 | not found | §6 treats the parity of the window ground state and resolves the odd sector at support 1.6, but as a certification, not as an `F(0)=0` explanation |
| Groskin 2605.20224 / 2607.02828 | not found | not found |
| Connes 2602.04022 | not found | not found |

Degenerate perturbation theory is textbook. **(T3)** What is not textbook is the identification of
the compression space with the prolate/near-null space *and* of the representer with `u sin(γu)`;
the first half is Connes–Consani–Moscovici, the second is the paper's own Theorem 8. The
combination is what this note tests, and it is a response law throughout.

---

## Deviations from the preregistration

1. **Sector.** The brief said "the EVEN basis that produced §4.3". §4.3 was in fact produced with
   **EVEN_D** (the default of `delta2_law.perturbed_C`). EVEN_D is used here, which is what
   reproduces the §4.3 numbers to `1.7e−10`.
2. **`δ` grid.** The sweep uses `δ = 10⁻¹…10⁻¹²`; K2 requires the §4.3 value, so `δ = 0.02` was
   evaluated separately (`src/stage2b_k2k3.py`).
3. **Stage 1 constant.** The preregistration left the constant open ("up to the constant from `a`");
   it is exactly 2: `λ_min = −2a²·max|S₁(c)|²`.
4. **Stage 4 scope reduction**, as recorded above: `1 − χ_k` not computed.
5. **`C(δ)` convention.** The preregistration defines `C = (λ_min(Q_0) − λ_min(Q_δ))/δ²`; the paper
   uses `−λ_min(Q_δ)/δ²`. They agree to `1.7e−10` at `δ=0.02` because `λ_min(Q_0) ≤ 2.5e−17`.

No retractions: nothing asserted earlier in this note was later found wrong.

---

## Verdict

> **PARTIAL.**
>
> - **The projection law describes the response.** K1, K2, K3 pass. The compressed form
>   `λ_min(Λ_N + δ²Π_N P Π_N)` is accurate to `10⁻²–10⁻⁵`; the scalar `4K cos²θ_N` to ~5 % at the
>   paper's `δ`, and it **explains the §4.3 saturation ratios**.
> - **The parity hypothesis is false.** K4 fails by 53×–893×. `F(0)=0` does not explain the odd
>   floor.
> - Stage 4 is confirmation of prior art, with one quantity not computed.

## Non-claims

This is a **response law**: how the window form reacts to a zero that has been moved by hand. It
is **not** a positivity mechanism and has **no bearing on RH**. All planted zero data are
**fabricated** — nothing here is evidence about the actual zeros of ζ. The floor `λ*(L)` is taken
from the existing certification and from Zhu; no new lower bound is claimed. Stage 4 claims nothing
new. H3 is reported as **failed**, not as a weaker version of itself.

## Reproduction

```
python3 src/stage1_ff.py --dps 60      # gate; also --dps 120 for the doubling check
python3 src/stage1_belowd.py           # R < d
python3 src/stage2_zeta.py             # delta sweep, H1/H2          ~5 min
python3 src/stage2b_k2k3.py            # K2 at delta=0.02, K3        ~4 min
python3 src/stage3_parity.py           # H3 / K4
python3 src/stage45.py                 # floor and combined formula
python3 src/figures.py
```
