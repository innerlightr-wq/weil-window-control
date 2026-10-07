# PROJECTION_CONTEXT_CHECKS.md

Checks for [`PROJECTION_CONTEXT.md`](PROJECTION_CONTEXT.md), run by
`check_projection_context.py` (standard library; exact `Fraction` arithmetic for all matrices,
kernels, projections and reduced inverses, plus the 60/70-digit `Decimal` Jacobi solver imported
**read-only** from the existing `kernel2.py`). No existing script, result or data file was altered.

**Provenance of the 2026-10-07 correction.** The correction to the quartic criterion was *prompted
by discussion*, then **derived and verified independently here** before any documentation changed.
It is not taken on authority. Numerical fitting alone was not used: every structural claim below is
decided in exact rational arithmetic.

## Checks run NOW

### §1 — `[P_old, L] = 0` needs both intertwinings
Synthetic tower step, `|K_e| = 2`, `|K_{e+1}| = 4`, fibre size `m = 2`, with a deliberately
**non-normal** coarse operator. Confirms `S J = m·id`; `P_old = JS/m` idempotent and self-adjoint;
range = fibre-constant functions; both intertwining forms; the vanishing commutator;
non-normality; and a **control** showing that an operator with `L J = J L` but `S L ≠ L S` has
`[P_old, L] ≠ 0`.

### §2 — qualification A
`‖P_E X_c‖² = V_mu` iff `X_c ∈ E`, with an explicit `G` giving equality at `dim E = 2 < 3 =`
`dim(mean-zero)`: equality does **not** require `E` to be the whole mean-zero subspace.

### §3 — the positive-quartic example
`Spec(A_0) = {0,1,12}`, `E = span(1,−2,1)`, `P_E X = 0`, `c4 = 1/6`, and
`λ_min(A_eps) → +eps⁴/6` confirmed in 60-digit `Decimal`.

### §4 — the corrected criterion and four exact controls

Structural identities verified on **every** control: `s = <z,Gz> + 2s²`; `0 < s ≤ 1/2` so
`eta ≥ 0`; the criterion **`eta = 0 ⟺ ker(G)` contains a vector of nonzero mean**; and
sufficiency of `G·1 = 0`.

All nodes carry counting measure and the Euclidean inner product; `1` is the all-ones vector.

| control | `X` | `G` | `ker(A_0)` | `A_0⁺1` | `s` | `eta` | `|<X²,v̂>|²` | `c4` | `G·1 = 0`? | `ker G` has nonzero mean? |
|---|---|---|---|---|---|---|---|---|---|---|
| **A** | `(−1,0,1)` | `g gᵀ`, `g = 1+X = (0,1,2)` | `(1,−2,1)` | `(5/12, 1/6, −1/12)` | `1/2` | `0` | `2/3` | `0` | **no** (`(0,3,6)`) | **yes** (`(1,0,0)`) |
| **B** | `(−1,0,1)` | `6P_1 + P_X` | `(1,−2,1)` | `(1/12,1/12,1/12)` | `1/4` | `1/4` | `2/3` | `1/6` | no | no |
| **C** | `(−1,0,1)` | `3 X Xᵀ` | `(1,−2,1)` | `(1/6,1/6,1/6)` | `1/2` | `0` | `2/3` | `0` | **yes** | yes |
| **D** | `(0,1,2,3)` | `I − vvᵀ/20`, `v = (1,−3,3,−1)` | `(1,−3,3,−1)` | `(1/9,1/9,1/9,1/9)` | `4/9` | `1/18` | `0` | `0` | no | no |

- **A is the refutation**: `G·1 ≠ 0` yet `eta = 0`, so the earlier "bracket vanishes iff `G·1 = 0`"
  is false. It does **not** establish a nonzero sixth-order term.
- **C** is the original sufficient condition, with the **ordinary outer product** `G = 3XXᵀ` on
  `X = (−1,0,1)`; stated explicitly because `REPORT.md` case (C) used a different normalisation
  (4 nodes `x = (0,1,3,4)`, `G = 3|x_c><x_c|` with `x_c` the centred coordinate). Both satisfy
  `G·1 = 0`.
- **D** shows the second, independent mechanism: `c4 = 0` through the **moment factor**
  `<X²,v̂> = 0` while the bracket `eta = 1/18` is strictly positive.

Sixth-order ratios for A and D were **computed here** at 70 digits (`λ/eps⁶ → −1/3` and `−1/10`).
This is reported as a computation on those two examples; `c4 = 0` by itself does not establish a
nonzero sixth-order term nor fix its sign.

## Results carried over, NOT re-run

- `REPORT.md` case (C)'s sixth-order value `−0.8` — from that report's own numerics.
- `REPORT.md` §4's finding that the toy quartet projection reduces `K` by only ~0.3% against a
  measured `C/4K` of 0.419–0.963.
- The exploration's `L^{s+2}` scaling law and the rank-two identity's published status.

## Note on precision

The quartic and sixth-order confirmations cannot be done in `float64`: at `eps = 1e-3` the
eigenvalue is `~1e-13` to `~1e-19` while the matrix entries are `O(1)–O(10)`, so it is formed by
cancellation at the round-off level (a double-precision solve returns `0.173` where the answer is
`1/6`). The checks therefore use the high-precision `Decimal` solver, and test **convergence**
rather than a fixed tolerance. No conclusion is drawn in `float64` about the sign or scale of a
tiny eigenvalue. **Stable digits are not a rigorous enclosure** and none is claimed.

## Transcript

```
==============================================================================
1.  [P_old, L] = 0 needs BOTH intertwinings  (elementary deduction, checked)
==============================================================================
  [PASS] S J = m . id   Definition 2 of the source
  [PASS] P_old = J S / m is idempotent
  [PASS] P_old is self-adjoint
  [PASS] range(P_old) = fibre-constant functions
  [PASS] Theorem 6 form:  L_{e+1} J = J L_e
  [PASS] Lemma 8 form:    S L_{e+1} = L_e S
  [PASS] => [P_old, L_{e+1}] = 0
  [PASS] L_{e+1} is NOT normal (so invariance of the complement is extra information)
  [PASS] control: an operator with L J = J L but S L != L S has [P_old,L] != 0   pullback invariance alone is not enough

==============================================================================
2.  Qualification A:  ||P_E X_c||^2 = V_mu  <=>  X_c in E   (not E = mean-zero)
==============================================================================
  [PASS] chosen v has <v,1> = 0 and <v,X_c> = 0
  [PASS] X_c lies in E = {1,v}^perp
  [PASS] hence ||P_E X_c||^2 = V_mu with dim E = 2 < 3 = dim(mean-zero)   V_mu = 10, dim E = 2, dim mean-zero = 3  => E need NOT be the whole mean-zero subspace
  [PASS] coefficient keeps its factor 2:  2||P_E X_c||^2   = 20

==============================================================================
3.  Qualification B: the proposed counterexample to a GENERAL sixth-order claim
==============================================================================
  [PASS] <1,X> = 0
  [PASS] A_0 = 4J + XX^T/2
  [PASS] Spec(A_0) = {0,1,12}   A_0 1 = 12.1,  A_0 X = 1.X,  A_0 w = 0
  [PASS] E = ker(A_0) = span{(1,-2,1)}   dim 1, since 1 and X carry eigenvalues 12 and 1
  [PASS] P_E X = 0   <X,w> = 0
  [PASS] predicted quartic coefficient = 1/6   |<X^2,vhat>|^2 = 2/3, <1,A_0^+1> = 1/4, coeff = 1/6
    (60-digit Decimal; float64 fails below eps ~ 1e-2 by cancellation)
    eps        lambda_min          lambda_min/eps^4   (predict 1/6 = 0.16666667)
    1e-2        1.666638888e-09       0.1666638888   |r-1/6| = 2.78e-06
    1e-3        1.666666389e-13       0.1666666389   |r-1/6| = 2.78e-08
    1e-4        1.666666664e-17       0.1666666664   |r-1/6| = 2.78e-10
    deviation shrinks by factors ['100', '100'] per decade (expect ~100 = O(eps^2))
  [PASS] lambda_min(A_eps) -> + eps^4/6 : POSITIVE, QUARTIC, not sixth order   final |r-1/6| = 2.8e-10

==============================================================================
4.  The quartic coefficient and the CORRECT bracket criterion   [corrected 2026-10-07]
==============================================================================
    Hypotheses: finite dimension, Euclidean inner product (counting measure),
    G = G* >= 0, A_0 = G + 2|1><1|, E = ker(A_0) = span(vhat) ONE-DIMENSIONAL with
    ||vhat|| = 1, positive gap on E-perp, and <X,vhat> = 0 (quadratic term absent).
    Then  c4 = <vhat,D vhat> - <B vhat, A_0^+ B vhat> = |<X^2,vhat>|^2 * eta,
    with B_ij = (X_i-X_j)^2,  D_ij = (X_i-X_j)^4/12,  eta = 1/2 - <1, A_0^+ 1>.

    SUPERSEDED: an earlier draft of PROJECTION_CONTEXT.md said the bracket vanishes
    'exactly when G.1 = 0'. That is too strong. G.1 = 0 is SUFFICIENT, not necessary;
    the exact condition is that ker(G) contain a vector of NONZERO MEAN. Control A
    below refutes the old statement. (Retained here so the history is not erased.)

    4a.  structural identities, on every control below
  [PASS] [A] s = <z,Gz> + 2s^2
  [PASS] [A] 0 < s <= 1/2 and eta >= 0   s=1/2, eta=0
  [PASS] [A] CRITERION  eta == 0  <=>  ker(G) has a nonzero-mean vector   eta=0, kerG nonzero-mean=True
  [PASS] [A] G.1 = 0  =>  eta = 0  (sufficiency)
  [PASS] [B] s = <z,Gz> + 2s^2
  [PASS] [B] 0 < s <= 1/2 and eta >= 0   s=1/4, eta=1/4
  [PASS] [B] CRITERION  eta == 0  <=>  ker(G) has a nonzero-mean vector   eta=1/4, kerG nonzero-mean=False
  [PASS] [B] G.1 = 0  =>  eta = 0  (sufficiency)
  [PASS] [C] s = <z,Gz> + 2s^2
  [PASS] [C] 0 < s <= 1/2 and eta >= 0   s=1/2, eta=0
  [PASS] [C] CRITERION  eta == 0  <=>  ker(G) has a nonzero-mean vector   eta=0, kerG nonzero-mean=True
  [PASS] [C] G.1 = 0  =>  eta = 0  (sufficiency)
  [PASS] [D] s = <z,Gz> + 2s^2
  [PASS] [D] 0 < s <= 1/2 and eta >= 0   s=4/9, eta=1/18
  [PASS] [D] CRITERION  eta == 0  <=>  ker(G) has a nonzero-mean vector   eta=1/18, kerG nonzero-mean=False
  [PASS] [D] G.1 = 0  =>  eta = 0  (sufficiency)

    4b.  the four controls, exact expected values
  [PASS] A: A_0 = [[2,2,2],[2,3,4],[2,4,6]]
  [PASS] A: ker(A_0) = span(1,-2,1)
  [PASS] A: P_E X = 0
  [PASS] A: G.1 = (0,3,6) != 0
  [PASS] A: A_0^+ 1 = (5/12,1/6,-1/12)
  [PASS] A: <1,A_0^+1> = 1/2, eta = 0, c4 = 0
  [PASS] A: REFUTES 'bracket vanishes iff G.1 = 0'
  [PASS] B: ker(A_0) = span(1,-2,1)
  [PASS] B: <1,A_0^+1> = 1/4
  [PASS] B: |<X^2,vhat>|^2 = 2/3
  [PASS] B: c4 = 1/6  (positive quartic)
  [PASS] C: G.1 = 0 and eta = 0 and c4 = 0
  [PASS] D: ker(A_0) = span(v), v = (1,-3,3,-1) up to sign
  [PASS] D: <1,v> = <X,v> = <X^2,v> = 0
  [PASS] D: <1,A_0^+1> = 4/9, eta = 1/18 > 0
  [PASS] D: c4 = 0 via the MOMENT factor while the bracket is positive

    4c.  c4 = 0 is O(eps^6); the sixth-order coefficient is COMPUTED here, 70 digits
      control A: lam/eps^4 -> -3.333e-09 (c4 = 0 confirmed), lam/eps^6 -> -0.33333334
  [PASS] A: lam_min/eps^4 -> 0, consistent with c4 = 0
      control D: lam/eps^4 -> -1.000e-09 (c4 = 0 confirmed), lam/eps^6 -> -0.10000001
  [PASS] D: lam_min/eps^4 -> 0, consistent with c4 = 0
      (these sixth-order values are computed; c4 = 0 by itself does NOT establish
       a nonzero sixth-order term, nor fix its sign.)

==============================================================================
ALL CHECKS PASSED
==============================================================================
```
