# Revision note for *Riemann Hypothesis: The Weil Enforcer in Audit Coordinates*

**Record revised:** De Jesús, E. (2026). *Riemann Hypothesis: The Weil Enforcer in
Audit Coordinates.* Zenodo. [doi:10.5281/zenodo.21115524](https://doi.org/10.5281/zenodo.21115524)

**Status:** erratum / scope correction, arising from the function-field control experiment
in *A δ²–L³ Sensitivity Law for Finite-Window Weil Positivity* (companion note). Nothing
below affects the function-field results of the Enforcer note; the correction is to how far
Block F was claimed to transfer to ζ.

---

## What stands

**Block F is structural in the function field.** In the recorder split `T_R = A − 2M` of the
window Toeplitz form for a curve `X/F_q`, the pole part of `A` is

```
q^{|n|/2} + q^{−|n|/2} = u_i w_j + w_i u_j ,    u_i = q^{i/2},  w_i = q^{−i/2},
```

an **exact rank-2 hyperbolic plane** (signature `(1,·,1)`), whose scale grows like `q^{R/2}`
against a fixed `2g·I`. Consequently `A` is *eventually* indefinite for every window
convention tested — the brief convention from `R = 2`–3, the uniform convention from
`R = 1` — and `M` is indefinite throughout, while `S = T_R` is PSD by Weil's theorem.
This is verified for four curves (`g = 1, 2, 3`) in `data/stage1_recorder_split.csv`.

**Theorem 8 of the Enforcer note is unaffected.** It is a function-field statement and the
hyperbolic-plane mechanism above is exactly what underwrites it.

## What must be narrowed

**On ζ windows, the indefiniteness of `A` and `M` is a bookkeeping artefact.** Writing the
ζ window form as `Pole + Arch − log(π)·I − Prime`, two equally natural splits give the same
`S`:

```
I :  A = Pole + Arch − log(π)·I ,   2M = Prime
II:  A = Pole + Arch ,              2M = Prime + log(π)·I
```

| L | conv I: `A` / `M` | conv II: `A` / `M` |
|---|---|---|
| 0.5 | (1,0,11) / (6,0,6) | **(0,0,12) / (0,0,12)** |
| 0.8 | (1,0,11) / (6,0,6) | **(0,0,12) / (0,0,12)** |
| log 13/2 | (2,0,10) / (5,0,7) | (0,0,12) / (2,0,10) |

(N = 12, 40 digits, identical in all three bases tested — two even, one odd;
`data/stage3d_inertia_reconciled.csv`.) **Under Convention II both `A` and `M` are
positive definite at `L ≤ 0.8`.** So "A and M are both indefinite" is not a property of the
ζ window form; it depends on where the `−log π` term is placed.

**The structural reason the transfer fails.** The function-field `A` inherits a hyperbolic
plane because the poles of `ζ_X` sit at `T = 1` and `T = 1/q`, i.e. **off** the critical
circle, contributing `q^{n/2} + q^{−n/2}`. For `ζ(s)` the poles at `s = 0, 1` contribute, for
**even** test functions,

```
h(i/2) + h(−i/2) = 2 F(i/2) F(−i/2) = +2 F(i/2)² ≥ 0 ,
```

which is **rank one and definite**, not hyperbolic (and rank-one *negative* in the odd
sector, where `F(−i/2) = −F(i/2)`). There is no hyperbolic plane on the ζ side to find.

## Recommended amendment

Wherever the Enforcer note asserts Block F for ζ or for "audit coordinates" generally,
restrict the claim to the function field and add: *on ζ windows the split is not inertia-
invariant; the indefiniteness of A and M depends on the placement of the `−log π` term, and
under one natural convention both are positive definite for `L ≤ 0.8`.* Logged as Retraction 15 in the companion note's retraction log
(doi:10.5281/zenodo.21115524 is the record this amends).

**One further caveat, carried over.** The `(A, M, S)` split was already reported as
basis-dependent in the 389a1 appendix and in Stage 1 of the companion work. This note
sharpens that from "the inertia values move" to "whether `A` and `M` are indefinite at all
moves". Only `S` is invariant.
