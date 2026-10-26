# notes/ff_construction.md — the function-field analogue, decided before computing

Phase 4a deliverable. PREREG: *"Document the construction in `notes/` BEFORE computing; if
no faithful analogue of `φ_{n,a}` exists for that basis, STOP 4a and record why."* The
named risk was that the function-field Gram is singular for `R > d`, so `G^{-1}` in (6.4)
would not exist.

## 1. The objects on each side

| | Suzuki, `ζ` | function field (`src/stage1.py`, Hallouin–Perret) |
|---|---|---|
| spectral set | `Γ ⊂ ℝ` (under RH), ordinates of `ξ` | `{±θ_j}_{j=1..g}`, Frobenius angles — **RH here is Weil's theorem**, so `θ_j ∈ ℝ` is a fact |
| test space | `X_{a,N} = span{χ_{n,a}}`, `\|n\| ≤ N` | `span{z^i}`, `i = 0..R` |
| transform | `φ_{n,a}(z) = √(2a) sin(az+nπ)/(az+nπ)` | `φ_i(w) = w^i`, evaluated on the unit circle |
| Gram | `G_{a,N} = Q_W(χ_m,χ_n)` | `T_R[i,k] = t(\|i−k\|)`, `t(n) = p_n q^{−n/2} = Σ_j 2cos(nθ_j)` |
| spectral identity | `G_{mn} = Σ_γ m_γ φ_m(γ)φ_n(γ)` | `cᵀT_R c = Σ_j 2\|P_c(e^{iθ_j})\|²`, `P_c(z) = Σ c_i z^i` (repo, verified to `2.3e−61`) |

> **The `φ` analogue EXISTS and is exact.** `φ_i(w) = w^i` on the unit circle generates
> `T_R` from the Frobenius angles in exactly the way `φ_{n,a}(γ)` generates `G_{a,N}` from
> the ordinates. The preregistered KILL for "no faithful `φ_n` analogue" does **not** fire.

## 2. The singularity risk, resolved

`ker T_R = {P·q}` with `dim = R+1−d` for `R ≥ d` (repo, `src/function_level.py`), so `T_R`
is **singular for `R > d`** and `T_R^{-1}` does not exist. For `R < d` the repo records
`T_R ≻ 0`. So the construction is available **only on `R < d`**, and that is where Phase 4
runs. The KILL as written ("no faithful analogue of `φ_{n,a}`") does not fire; the
singularity is a *range restriction*, not an obstruction, and is recorded as such.

## 3. Where the analogue genuinely stops — and what replaces it

`W_N(a,π;z) = (z−i)J(−i,z) − (z+i)J(i,z)` is **not** a formula about a Gram matrix. It
comes from the deficiency structure of `D_a = i d/dx` on `H(A_a)`: deficiency indices
`(1,1)`, deficiency vectors `u_{±i,a}` at the two points `±i` off the real axis, and a
circle of self-adjoint extensions `D̄_{a,θ}` whose eigenvalues are the zeros of
`W(a,θ;·)` (§3.2). **The function-field side has no differential operator, no deficiency
indices, and no `θ`-family**, so `±i` and the combination above have no literal analogue.

Transporting `W_N` verbatim would therefore be a notational exercise, not a construction.
But the *structure* does transport, through the classical theory the function-field Gram
already belongs to:

* `T_R` is a **Toeplitz** matrix — the moment matrix of the measure
  `μ = Σ_j 2(δ_{θ_j} + δ_{−θ_j})` on the unit circle.
* `J_{a,N}(w,z) = φ(z)ᵀ conj(G^{-1}φ(w))` is a **reproducing kernel**; its Toeplitz
  counterpart is the **Christoffel–Darboux kernel** `K_R(w,z) = φ(z)ᵀ conj(T_R^{-1}φ(w))`
  — the same formula, the same object.
* The role of the `θ`-circle of self-adjoint extensions is played by the
  **para-orthogonal polynomials on the unit circle (POPUC)**,
  `B_R(z;β) := Φ_R(z) − β Φ*_R(z)`, `|β| = 1`, where `Φ_R` is the monic Szegő polynomial
  and `Φ*_R(z) = z^R conj(Φ_R(1/conj z))`. These are the characteristic polynomials of the
  rank-one unitary (Gauss-quadrature) perturbations, exactly as `W(a,θ;·)` is the spectral
  function of the rank-one self-adjoint perturbations `D̄_{a,θ}`.
* **Classical fact (Jones–Njåstad–Thron; Simon, *OPUC*):** if the moment matrix is
  **positive definite**, all zeros of `B_R(·;β)` lie **on the unit circle** and are simple,
  and they interlace the support of `μ`.

So the faithful statement of the Phase 4 question is:

> `W_N` ↔ `B_R(·;β)`; "zeros real" ↔ "zeros on the unit circle"; `θ = π` ↔ a fixed `β`.

## 4. What this predicts for Phase 4c — recorded BEFORE computing

PREREG 4c predicts that the `δ` at which `W_N` first acquires non-real roots coincides
(to 10%) with the `δ` at which `λ_min` turns negative, and reads a coincidence as
*"the `δ²`–`L³` law governs a second pipeline"*.

**On the function-field side that coincidence is forced by the classical fact above, not
discovered.** Zeros of `B_R` on the circle is *equivalent* to positive definiteness of the
moment matrix (one direction is the theorem quoted; the converse fails only in degenerate
cases). So if the two thresholds agree, it is **one phenomenon counted twice**, and the
agreement is evidence about Szegő theory, not about the `δ²`–`L³` law governing a second
pipeline.

> **Preregistered amendment, made before running 4c.** The outcome "coincide within 10%"
> must therefore **not** be reported as *"the `δ²`–`L³` law governs a second pipeline"* on
> the function-field side. The informative measurements are instead:
> 1. **how exactly** the two thresholds agree (if the mechanism is structural the agreement
>    should be far tighter than 10% — at the level of the eigenvalue solve, not 10%);
> 2. whether the **rate** at which the zeros leave the circle follows the `δ²` law, which
>    is a quantitative statement the structural equivalence does **not** supply.
>
> A tight agreement is then a *consistency check on the implementation*; (2) is the only
> part that could carry new information.

This amendment is recorded here rather than applied silently, because it weakens a
preregistered claim before the data exists.
