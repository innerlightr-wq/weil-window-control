# The parity sign split, in the function field

Step 2 of the 2026-10-23 verification pass. Everything here is **exact over ℚ** in
`fractions.Fraction`; no floating point enters. Script:
[`src/ff_parity_split.py`](../src/ff_parity_split.py). Tier **T1**.

## The correspondence

The function-field control of the parent repository (Theorem 6, `notes/proofs.md` §A) puts a
single off-circle **real** pair of normalised Frobenius eigenvalues `{ρ, 1/ρ}` into the
window and gets the rank-two Toeplitz matrix

```
T_R = u wᵀ + w uᵀ,     u_i = ρ^i,  w_i = ρ^{−i},   0 ≤ i ≤ R.
```

This is the same object as the ζ/Dirichlet-side Siegel configuration: a real pair, not a
quartet. The dictionary is

| function field | finite window on `L(s, χ)` |
|---|---|
| flip `J : x_i ↦ x_{R−i}` | test-function parity `f(u) ↦ f(−u)` |
| `J`-palindromic sector | **even** sector |
| `J`-anti-palindromic sector | **odd** sector |
| `+` rank-one piece, on the palindromic direction | `+2A²`, `A = ∫ f cosh(δu) du` |
| `−` rank-one piece, on the anti-palindromic direction | `−2S²`, `S = ∫ f sinh(δu) du` |
| `λ₊ = (R+1) + ‖u‖‖w‖ > 0` | even-sector response, `≥ 0` |
| `λ₋ = (R+1) − ‖u‖‖w‖ ≤ 0` | odd-sector response, `≤ 0` |
| equality `λ₋ = 0` iff `ρ = 1` (Corollary A3) | zero response iff `δ = 0` |

So the Phase-5 result is not an artefact of the ζ-side bookkeeping. **The same signed
rank-one split, with the same attachment of sign to parity, is visible in the setting where
the Riemann hypothesis is Weil's theorem** — and there it is a finite exact computation over
ℚ rather than a high-precision one.

## The split, and the weight the brief omits

Let `J` be the flip. Since `u_{R−i} = ρ^R w_i` and `w_{R−i} = ρ^{−R} u_i`,

```
J u = ρ^R w,        J w = ρ^{−R} u,
```

so **`J` does not fix `u + w`.** It exchanges `u` and `w` only after the `ρ^R` weight is
restored. Put

```
v₊ = u + ρ^R w          J v₊ = + v₊      (palindromic)
v₋ = u − ρ^R w          J v₋ = − v₋      (anti-palindromic)
```

Then, exactly over ℚ for rational `ρ`:

> **(FF1)** `T_R = (1 / 2ρ^R) [ v₊ v₊ᵀ − v₋ v₋ᵀ ]`
>
> **(FF2)** `T_R v₊ = λ₊ v₊` and `T_R v₋ = λ₋ v₋` with
> `λ± = (R+1) ± ρ^{−R}‖u‖²`
>
> **(FF3)** `‖w‖² = ρ^{−2R}‖u‖²` **identically**, hence `ρ^{−R}‖u‖² = ‖u‖‖w‖`, so `λ±` are
> the parent repository's `(R+1) ± ‖u‖‖w‖` and are **rational** even though `‖u‖` and `‖w‖`
> individually are not.
>
> **(FF4)** `λ₊ > 0` and `λ₋ ≤ 0`, with `λ₋ = 0` iff `ρ = 1`.

(FF3) is what makes the whole thing exact over ℚ, and it is also the reason `J` works: the
normalised vectors satisfy `J û = ŵ` and `J ŵ = û` exactly, because `‖u‖/‖w‖ = ρ^R` is
precisely the weight `J` introduces. `v±` are `‖u‖·(û ± ŵ)`.

**Verified in all 25 cells** (`ρ ∈ {2, 3/2, 5/4, 7/3, 1}` × `R ∈ {1, 2, 4, 8, 16}`),
exactly, zero residual.

## Two corrections to the commissioning brief

Both are checked rather than argued, in `src/ff_parity_split.py`.

1. **The brief's identity is correct.** `T_R = ½[(u+w)(u+w)ᵀ − (u−w)(u−w)ᵀ]` holds exactly
   — it is just the polarisation identity, and the `ρ^R` weights cancel in it. **PASS in all
   25 cells.**

2. **But `u+w` and `u−w` are neither the eigenvectors nor the parity vectors.** The brief
   says the eigenvalues `(R+1) ± ‖u‖‖w‖` are "attained on `u+w` and `u−w` respectively" and
   that `u+w` is palindromic. Both are **false for every `ρ ≠ 1`**, refuted in all 20
   non-degenerate cells. Directly:
   `T_R(u+w) = u[(R+1) + ‖w‖²] + w[‖u‖² + (R+1)]`, which is parallel to `u+w` only when
   `‖u‖² = ‖w‖²`, i.e. `ρ^{2R} = 1`. The failure is exactly the missing `ρ^R`.

   At `ρ = 1` both brief claims become true and `λ₋ = 0` — consistent with the parent's
   Corollary A3, where strictness fails precisely at `ρ = 1`. So the brief's version is the
   degenerate-case statement.

Neither correction touches any parent-repository result: `notes/proofs.md` §A.2 derives the
eigenvalues basis-independently, from the `2×2` matrix `[[p, ‖w‖²], [‖u‖², p]]`, and never
claims `u ± w` as eigenvectors.

## What this does and does not add

It **does** show the signed-parity structure is intrinsic to "a real pair in a finite
window", not a feature of the archimedean bookkeeping of `L(s, χ)`.

It does **not** transfer the Phase-5 consequence to the function field as a statement about
Siegel zeros: in the function field the analogue of RH is a theorem, so there is no
configuration to exclude. The function field is a control here, exactly as it is in Part I
of the parent paper.
