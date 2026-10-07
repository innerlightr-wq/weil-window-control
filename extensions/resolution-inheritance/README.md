# No resolution inheritance in the finite-window Weil form

A short **negative result**. Not part of the v1.0 technical note
([doi:10.5281/zenodo.23195177](https://doi.org/10.5281/zenodo.23195177)); no review, no DOI.

Tiers: **T1** exact/proved, **T2** numerical with stated precision, **T3** interpretation.

## Question

In the author's prime-power transfer-operator work, refining the resolution gives an *exact
reducing* decomposition `L_{e+1} ≅ L_e ⊕ N_{e+1}`. Does the Weil quadratic form show any
analogue when the window is enlarged, `H_L ⊂ H_{L'}`?

This was a reconnaissance run with an explicit falsification gate, and it was asked
adversarially: assume nothing, try to break it.

## Verdict: FAIL — no inheritance, exactly and for a structural reason

## Method

The project's basis `cos(kπx/L)` on `[−L,L]` is **not nested in `L`**, so the question
cannot be asked with `zeta_window.build_matrix` alone. `src/wform.py` evaluates `Q(f,g)` for
arbitrary compactly-supported piecewise-cosine `f,g` (same functional as
[`notes/stage3_assembly.md`](../../notes/stage3_assembly.md), `m`-sum resummed in closed
form). That allows a genuinely nested orthonormal basis:

```
H_{L'} = H_L (+) H_new,   H_L  = cosines on [-L, L], extended by zero
                          H_new = even functions on the annulus L <= |x| <= L'
```

**Validation (T2).** The evaluator reproduces `build_matrix` to full working precision —
rel. `4.2e-31` at dps 30, `4.3e-41` at dps 40 (`data/validate.log`). The closed-form
correlation was separately checked against brute-force quadrature across every piece-type
pair, worst `9.9e-32` at dps 30.

## Results

`E_old = ‖Q^old_{L'} − Q_L‖_F/‖Q_L‖_F`; `E_cross = ‖P_old Q P_new‖_F/‖Q‖_F`;
`σ_x` = top singular value of the cross block over `‖Q‖_op` (invariant under orthonormal
changes of basis *inside* each block).

| L | L' | N_in | N_new | dps | basis | E_old | E_cross | σ_x/‖Q‖ | worst eig mismatch |
|---|---|---|---|---|---|---|---|---|---|
| 0.5 | 0.8 | 4 | 4 | 30 | even | 1.6e−42 | 0.2971 | 0.3819 | 0.491 |
| 0.5 | 0.8 | 6 | 6 | 30 | even | 1.3e−42 | 0.2481 | 0.3239 | 0.071 |
| 0.5 | 0.8 | 8 | 8 | 30 | even | 1.0e−42 | 0.2243 | 0.3245 | 0.136 |
| 0.5 | 0.8 | 12 | 12 | 30 | even | — | 0.1847 | 0.2891 | — |
| 0.5 | 0.8 | 16 | 16 | 30 | even | — | 0.1695 | 0.2859 | — |
| 0.5 | 0.8 | 8 | 8 | **20** | even | — | 0.2243 | 0.3245 | — |
| 0.5 | 0.8 | 8 | 8 | **50** | even | 4.7e−63 | 0.2243 | 0.3245 | — |
| 0.5 | 0.8 | 6 | 6 | 30 | **even_d** | 1.2e−42 | 0.2204 | — | 0.144 |
| 0.8 | 1.0 | 6 | 6 | 30 | even | 1.2e−42 | 0.2233 | — | 0.221 |
| 0.8 | 1.6 | 6 | 6 | 30 | even | 2.4e−43 | **0.4685** | — | **0.447** |
| 1.0 | 1.2 | 6 | 6 | 30 | even | 8.0e−45 | 0.2155 | — | 0.249 |

Coupling by component of the form (L=0.5→0.8, N=6, `data/parts.log`):

| part | ‖cross‖/‖part‖ |
|---|---|
| Pole (rank one) | 0.489 |
| Prime | 0.582 |
| Arch | 0.155 |
| **LogPi** | **exactly 0** |

### 1. `E_old = 0` is a tautology, not inheritance (T1)

It tracks precision exactly (`1e−42` at dps 30, `1e−63` at dps 50). But `Q` is **one fixed,
`L`-independent form**: for `f,g` supported in `[−L,L]` the correlation vanishes past `2L`,
so the extra prime powers and the longer archimedean integral contribute nothing to the old
block. *Any* isometric nesting gives `E_old = 0`.

**Control:** comparing same-*index* blocks of the project's own non-nested bases gives
`E_old` = 0.49–1.25 (`data/naive.log`) — a pure basis artifact, pointing the other way.
Neither number is evidence about the form.

### 2. The coupling is order one and does not decrease (T2)

The Frobenius ratio drifts down only because adding high annulus modes dilutes the
normalisation. The basis-free `σ_x/‖Q‖` **plateaus at ≈ 0.29** (0.38, 0.32, 0.32, 0.29, 0.29
for N = 4…16). At dps **20, 30 and 50 the values are identical to 4 significant figures**.
Nothing here is numerical.

### 3. The spectrum is not inherited (T2)

Under a reducing splitting, `Spec(Q_{L'}) = Spec(Q_L) ⊔ Spec(N)`. It does not hold:
`λ_min(Q_L) = 1.589e−6` has no counterpart in `Spec(Q_{L'})` (nearest off by 7%), while
`λ_min(Q_{L'}) = 2.31e−8`.

### 4. The ground state leaks and stays leaked (T2)

`‖v_new‖` = 0.131, 0.0768, 0.0773, 0.0843, 0.0895 for N = 4…16 — plateaus near 0.09, does
not tend to 0. Overlap of the inner part with `gs(Q_L)` plateaus at 0.989, not 1
(`data/overlap.log`).

### 5. The δ² coefficient is not additive (T2)

`C = 4F'(γ₁)²` split as `4S_I² + 4S_N² + 8 S_I S_N`:

| L → L' | C(L)+4S_N² | C(L') | rel err | cross term `8 S_I S_N` |
|---|---|---|---|---|
| 0.5 → 0.8 | 4.94e−3 | 5.46e−5 | **89.5** | −7.65e−3 (larger than C(L') itself) |
| 0.8 → 1.0 | 3.618e−5 | 3.879e−5 | 0.067 | −1.60e−7 |
| 0.5 → 0.8 (even_d) | 6.93e−4 | 3.67e−3 | 0.81 | +9.73e−4 |

The `0.8 → 1.0` row looks additive only because `‖v_new‖ = 4e−5` there: the new sector is
unoccupied, so it is the trivial case. Even then the increment is `2.6e−6` against the
`4L²ΔL ≈ 0.512` that `L³` scaling would predict.

**On the `L²ΔL` form (T3):** any smooth `C(L) ~ L³` satisfies `C(L') − C(L) ≈ 3L²ΔL` by
Taylor. That is calculus, not inheritance, and it should not be read as structure.

**Caveat (T2 scope):** the `C` here is `4F'(γ₁)²` at the N-truncated ground state, *not* the
project's calibrated `C(L)` (EVEN_D, N = 16–26, perturbed-`λ_min` route in
[`src/delta2_law.py`](../../src/delta2_law.py)). Absolute values are not comparable to that
table; the conclusion is about the decomposition, not the value.

## Why (T1 obstruction, T3 framing)

The pole term `2 P_f P_g`, with `P_f = ∫ f(x) e^{−x/2} dx`, is **rank one**. A rank-one form
is block-diagonal only if one block lies in its kernel, and `e^{−x/2}` has nonzero projection
onto both `H_L` and the annulus. So `E_cross` is bounded below *a priori* — matching the
measured `Pole` cross ratio of 0.489.

The only part that does inherit is `−log π · δ(u)`, whose cross block is **exactly zero**
because it sees only `C(0) = ⟨f,g⟩`. It is a contact term. Every non-local part couples at
O(1).

**T3.** This is structurally unsurprising. The prime-power splitting needs *two*
intertwinings — `J_e` is a pullback along a group quotient `π_e`, and it is the fiber sum
`S_e` that upgrades invariance to *reducing*; invariance alone gives only block-triangularity.
The Weil nesting `H_L ⊂ H_{L'}` is an inclusion of supports with no symmetry behind it and no
`S` map, so there is nothing to make the complement invariant.

## Strongest counterexample

**L = 0.8 → L' = 1.6.** `E_cross = 0.469`; `λ_min` falls `4.138e−12 → 4.143e−17` (factor
10⁵); the old `λ_min` is absent from the new spectrum, nearest off by 45%. Under an exact
splitting it would have to survive unchanged.

## Non-claims

1. Nothing here bears on RH, or on the prime-power tower result — that theorem is unaffected,
   and the analogy was explicitly offered as motivation only.
2. The failure is shown for the **support nesting** `H_L ⊂ H_{L'}`. It does not rule out some
   other decomposition of the Weil form being reducing; none was looked for.
3. The δ² section is the weakest part (see its caveat) and is not converged in N.
4. Even sector only; `L ≤ 1.6`; `N ≤ 16`.

## Reproduction

Needs `mpmath`. From this directory:

```
python3 src/validate.py       # evaluator vs the project's build_matrix        ~6 s
python3 src/run_recon.py      # E_old, E_cross, spectrum, per pair            ~30 s
python3 src/conv.py           # convergence in N and dps, basis-free metric    ~60 s
python3 src/parts_delta2.py   # coupling by component; delta^2 additivity      ~15 s
python3 src/gs_overlap.py     # ground-state leak vs N                         ~45 s
python3 src/naive.py          # the same-index basis artifact control          ~10 s
```

Logs as run are in `data/`.

## License

Code MIT, prose and data CC BY 4.0 — the parent repository's split
([`LICENSE`](../../LICENSE), [`LICENSE-CC-BY-4.0`](../../LICENSE-CC-BY-4.0)).
