# Derivation: the finite-window Weil form for the Davenport–Heilbronn function

Every displayed identity is numerically checked by `src/phaseB.py` (coefficients and the
scalar form) and `src/phaseC.py` (the gate); the check that closes each step is named in
square brackets.

---

## 1. The function, the character, and `κ`

Let `χ` be the character mod 5 with

```
χ(0), χ(1), χ(2), χ(3), χ(4) = 0, 1, i, −i, −1.
```

`χ` is primitive mod 5 and **odd**: `χ(−1) = χ(4) = −1`. Its Gauss sum is
`τ(χ) = Σ_{n=1}^{4} χ(n) e^{2πin/5}`, with `|τ(χ)| = √5`, and the root number of the odd
primitive character is

```
ε_χ = τ(χ) / (i√5),      |ε_χ| = 1,
```

so that the completed `Λ(s,χ) = (5/π)^{(s+1)/2} Γ((s+1)/2) L(s,χ)` obeys
`Λ(s,χ) = ε_χ Λ(1−s, χ̄)`. **This equation exchanges `χ` and `χ̄`**; it is not
self-reciprocal, and that is the whole reason the Davenport–Heilbronn combination is
formed.

Put

```
F(s) = [ (1 − iκ) L(s,χ) + (1 + iκ) L(s,χ̄) ] / 2                              (1.1)
```

with `κ` real, and `Λ_F(s) = (5/π)^{(s+1)/2} Γ((s+1)/2) F(s)`, so that
`Λ_F = [(1−iκ)Λ(·,χ) + (1+iκ)Λ(·,χ̄)]/2`. Then

```
Λ_F(1−s) = [ (1−iκ) Λ(1−s,χ) + (1+iκ) Λ(1−s,χ̄) ] / 2
         = [ (1−iκ) Λ(s,χ̄)/ε_χ̄ + (1+iκ) Λ(s,χ)/ε_χ ] / 2,
```

and matching the coefficient of `Λ(s,χ)` against (1.1) gives the single real condition

```
(1 − iκ) ε_χ = 1 + iκ        ⟺        ε_χ = (1+iκ)/(1−iκ).                     (1.2)
```

Because `|ε_χ| = 1`, (1.2) is solvable, and writing `θ = arg ε_χ` the solution is

```
κ = tan(θ/2).                                                                  (1.3)
```

The coefficient of `Λ(s,χ̄)` then matches automatically: conjugating (1.2) and using
`ε_χ̄ = conj(ε_χ)` (valid because `χ̄ = conj χ` has `τ(χ̄) = χ(−1) conj τ(χ) = −conj τ(χ)`,
hence `ε_χ̄ = −conj τ(χ)/(i√5) = conj( τ(χ)/(i√5) ) = conj ε_χ`) gives
`(1+iκ) ε_χ̄ = 1 − iκ`, which is the second condition. **[A1]**

Evaluating: `τ(χ) = (ζ₅ − ζ₅⁴) + i(ζ₅² − ζ₅³) = 2i sin 72° − 2 sin 144°`, so
`ε_χ = (2 sin 72° + 2i sin 144°)/√5` and `tan θ = sin 144°/sin 72° = sin 36°/sin 72°
= 1/(2 cos 36°) = 1/φ`, `φ` the golden ratio. The half-angle gives the closed form

```
κ = [ √(10 − 2√5) − 2 ] / (√5 − 1) = 0.28407904384041229602829183239312616909…  (1.4)
```

which is the classical Davenport–Heilbronn constant. It is **derived here, not imported**;
`src/dh.py:kappa()` computes `tan(arg(ε_χ)/2)` from the Gauss sum and `phaseA.py [A1]`
checks it against (1.4). **[A1]**

### Dirichlet coefficients

`(1−iκ)χ(n) + (1+iκ)χ̄(n) = 2 Re χ(n) + 2κ Im χ(n)`, so `a_n = Re χ(n) + κ Im χ(n)` and

```
(a_1, a_2, a_3, a_4, a_5) = (1, κ, −κ, −1, 0),      a_{n+5} = a_n.             (1.5)
```

These are **real** and periodic but **not multiplicative**: `a_2 a_3 = −κ² ≈ −0.0807`
while `a_6 = a_1 = 1`. Hence no Euler product. (Non-claim 2: the absence of an Euler
product is a property of this function, not an asserted *cause* of its off-line zeros.)

Equivalently `F(s) = 5^{−s} Σ_{r=1}^{5} a_r ζ(s, r/5)`; since `Σ_r a_r = 0` the poles of
the five Hurwitz zetas at `s = 1` cancel and `F` is entire, with
`F(1) = −(1/5) Σ_r a_r ψ(r/5)`. **[A2]**

### The strip

For `Re s ≥ 2`, `|F(s) − 1| ≤ Σ_{n≥2} |a_n| n^{−2} ≤ Σ_{n≥2} n^{−2} = ζ(2) − 1 =
0.6449… < 1`, so `F(s) ≠ 0` there. The completed `Λ_F` adds no zeros (the `Γ` factor and
the exponential are zero-free), so `Λ_F(s) ≠ 0` for `Re s ≥ 2`; by
`Λ_F(s) = Λ_F(1−s)` also `Λ_F(s) ≠ 0` for `Re s ≤ −1`. **All zeros of `Λ_F` lie in the
open strip `−1 < Re s < 2`.** This is a T1 statement; it replaces any arbitrary right
cutoff. **[A5]**

Note that this strip is *wide*: the Davenport–Heilbronn function is known to have zeros
with `Re s > 1`, so the strip cannot be narrowed to a neighbourhood of the critical line
the way it can for `ζ`.

### Trivial zeros and entireness of `Λ_F`

`Γ((s+1)/2)` has simple poles at `s = −1, −3, −5, …`. Both `L(s,χ)` and `L(s,χ̄)` vanish
there (the standard trivial zeros of an odd Dirichlet `L`-function), hence so does `F`,
and `Λ_F` is entire of order 1. There is **no pole term** in the explicit formula below,
and the cancellation of the trivial zeros of `F` against the poles of `Γ` is exactly what
makes the spectral sum run over `Λ_F`-zeros only. **[A4]**

---

## 2. Logarithmic derivative

In `Re s > 2` write `−F'(s)/F(s) = Σ_{n≥2} b_n n^{−s}`. Multiplying by
`F(s) = Σ a_n n^{−s}` and comparing coefficients of `n^{−s}` in
`(−F'/F)·F = −F' = Σ a_n log n · n^{−s}` gives the Dirichlet convolution

```
Σ_{d|n} b_d a_{n/d} = a_n log n,      b_1 = 0,
```

whence, since `a_1 = 1`,

```
b_n = a_n log n − Σ_{d|n, 2 ≤ d < n} b_d a_{n/d}.                              (2.1)
```

Convergence: `|F(s) − 1| < 1` on `Re s ≥ 2` makes `log F` analytic there with an
absolutely convergent Dirichlet series in `Re s > 2`, so (2.1) is legitimate in
`Re s > 2`, which is the half-plane the contour in §3 uses.

**Second, independent route.** Build the Dirichlet inverse `c = a^{−1}`
(`c_1 = 1`, `c_n = −Σ_{d|n, d<n} c_d a_{n/d}`) and set `b_n = Σ_{d|n, d≥2} a_d log d ·
c_{n/d}`. `src/dh.py` implements both; they agree to the working precision. **[B1]**

**Exact checks.**

```
b_2 = a_2 log 2 = κ log 2,
b_3 = a_3 log 3 = −κ log 3,
b_6 = a_6 log 6 − (b_2 a_3 + b_3 a_2) = log 6 + κ² log 2 + κ² log 3 = (1 + κ²) log 6.
```

`b_6 ≠ 0` **at a non-prime-power** is the signature of the missing Euler product and is
the control used in gate B. Also `b_n = 0` whenever `5 | n` (every divisor chain hits a
factor with `a = 0`), and in particular `b_5 = b_{10} = 0`. **[B1]**

---

## 3. Transform conventions and the explicit formula

`g` real, supported in `[−L, L]`, continuous on `ℝ` (both bases below vanish at `±L`).

```
G(z) = ∫ g(u) e^{zu} du,   h(v) = ∫ g(u+v) g(u) du,   H(z) = ∫ h(v) e^{zv} dv = G(z)G(−z).
```

`h` is supported in `[−2L, 2L]`, `h(0) = ‖g‖₂²`, and `H(it) = |G(it)|²` for real `g`.
`H` is even. The spectral form is over the zeros of the **entire completed** function:

```
Q_F(g) = Σ_ρ mult(ρ) · H(ρ − ½),      ρ running over zeros of Λ_F.             (3.1)
```

### The contour

`Λ_F'/Λ_F(s) = ½ log(5/π) + ½ ψ((s+1)/2) + F'/F(s)` is holomorphic on and to the right of
`Re s = 2` and, by the functional equation, on and to the left of `Re s = −1` (§1). Let
`σ₀ = 2`. For a sequence of heights `T_j → ∞` avoiding the ordinates of zeros,

```
Σ_{|Im ρ| < T_j} H(ρ−½) = (1/2πi) [ ∫_{σ₀−iT_j}^{σ₀+iT_j} − ∫_{1−σ₀−iT_j}^{1−σ₀+iT_j} ]
                           (Λ_F'/Λ_F)(s) H(s−½) ds  +  (horizontal edges).
```

The horizontal edges vanish along a suitable sequence `T_j`: on `|σ| ≤ 2`,
`|H(σ+iT−½)| = O(T^{−4})` because `g` is continuous with a jump only in `g'` at `±L`, so
`G(z) = O(|Im z|^{−2}) e^{L|Re z|}`; and `Λ_F'/Λ_F = O(log T)` on the standard sequence of
heights. **No limiting term survives.** [Numerically monitored as the zero-side tail in
`phaseC.py`; the statement about the edges is T1 modulo the standard choice of `T_j`.]

Substituting `s → 1−s` in the second integral and using `Λ_F(s) = Λ_F(1−s)` (so
`(Λ_F'/Λ_F)(1−s) = −(Λ_F'/Λ_F)(s)`) together with `H` even gives

```
Q_F(g) = 2 · (1/2πi) ∫_{(2)} (Λ_F'/Λ_F)(s) H(s−½) ds.                          (3.2)
```

### The arithmetic term

With `H(s−½) = ∫ h(v) e^{(s−½)v} dv` and
`(1/2πi)∫_{(σ₀)} e^{s x} ds = δ(x)` against the compactly supported continuous `h`,

```
(1/2πi) ∫_{(2)} n^{−s} H(s−½) ds = ∫ h(v) e^{−v/2} δ(v − log n) dv = h(log n)/√n,
```

so, since `F'/F(s) = −Σ_{n≥2} b_n n^{−s}`,

```
2 · (1/2πi) ∫_{(2)} (F'/F)(s) H(s−½) ds = −2 Σ_{n≥2} b_n h(log n)/√n.          (3.3)
```

`h` vanishes outside `[−2L, 2L]`, so the sum is **finite**: only `log n < 2L`, i.e.
`n ≤ e^{2L}`, contributes. The cutoff is applied as an integer bound `n ≤ ⌊e^{2L}⌋ + 2`
followed by an arbitrary-precision re-test `log n < 2L`, never by a binary-float
comparison. At `L = 2.6`, `2L = 5.2` and `e^{5.2} = 181.27…`, so `n` runs to **181**.
**[B2]**

### The archimedean and conductor terms

`½ log(5/π) + ½ ψ((s+1)/2)` is holomorphic for `Re s > −1`, so the contour may be moved
from `Re s = 2` to `Re s = ½`. With `s = ½ + it`, `ds = i dt`:

```
2 · (1/2πi) ∫_{(½)} [ ½ log(5/π) + ½ ψ(¾ + it/2) ] H(it) dt·i
   = (1/2π) ∫_ℝ [ log(5/π) + ψ(¾ + it/2) ] |G(it)|² dt.
```

`(1/2π)∫ |G(it)|² dt = h(0)` by Parseval, and `Im ψ(¾+it/2)` is odd in `t` while
`|G(it)|²` is even, so the imaginary part integrates to zero. Hence

```
A_{5,odd}(g) = h(0) log(5/π) + (1/2π) ∫_ℝ |G(it)|² Re ψ(¾ + it/2) dt.          (3.4)
```

### The form

```
┌──────────────────────────────────────────────────────────────────────────────┐
│  Q_F(g) = A_{5,odd}(g) − 2 Σ_{2 ≤ n ≤ e^{2L}} b_n h(log n)/√n                │
└──────────────────────────────────────────────────────────────────────────────┘
```

with `A_{5,odd}` as in (3.4). This is exactly the target of the brief, with every factor
accounted for: the `2` in front of the arithmetic sum comes from the doubling in (3.2),
and the conductor enters as `log(5/π)` because the `π` of the gamma factor and the `5` of
the conductor sit in the same prefactor `(5/π)^{(s+1)/2}`.

**There is no pole term** (`Λ_F` is entire, §1) — the one structural difference from the
parent repository's `ζ` form, where `2F(i/2)²` appears.

### Relation to the parent repository

Setting `a = 1` (odd gamma factor) and `q = 5`, this is precisely the matrix assembled by
`extensions/dirichlet-windows-2026-10-23/src/dirichlet_window.build_parts_chi`, with
`Λ(n)χ(n)` replaced by `b_n`:

```
M_DH = Arch_a(a=1) − LogPi + LogQ(5) − Prime_b,      no Pole.
```

`Arch_a(·, L, 1)` is the already-verified closed form for
`(1/2π)∫|G(it)|² Re ψ(¾+it/2) dt`. **The `L`-function Weil matrices are NOT linearly
combined**: `−F'/F` is not a linear combination of `−L'/L(·,χ)` and `−L'/L(·,χ̄)`, because
logarithmic differentiation does not commute with linear combination. Only the
*archimedean and conductor* machinery is reused, and those depend on the gamma factor and
the conductor alone, both of which `F` shares with `L(s,χ)`.

---

## 4. Orbits

For `g` real, `G(z̄) = conj G(z)`, hence `H(z̄) = conj H(z)` and `H(−z) = H(z)`.

**Non-real off-line zero.** The complete functional-equation + conjugation orbit of
`ρ = β + iγ` with `β ≠ ½`, `γ ≠ 0` is

```
O = { β + iγ, 1 − β + iγ, β − iγ, 1 − β − iγ },
```

four distinct points. With `δ = β − ½` the shifts `ρ − ½` are `{ ±δ ± iγ }` and

```
Q_O(g) = H(δ+iγ) + H(−δ+iγ) + H(δ−iγ) + H(−δ−iγ)
       = 2 [ H(δ+iγ) + H(−δ+iγ) ]                                 (H even)
       = 4 Re( A · conj B ),     A = G(δ+iγ),  B = G(−δ+iγ).                   (4.1)
```

Writing `P = ∫ g cosh(δu) e^{iγu} du` and `S = ∫ g sinh(δu) e^{iγu} du`, so `A = P + S`
and `B = P − S`, and `Re(A conj B) = |P|² − |S|²`:

```
┌──────────────────────────────────────────────┐
│   Q_O(g) = 4 ( |P|² − |S|² )                 │                               (4.2)
└──────────────────────────────────────────────┘
```

This is a **signed rank-two** form: `+4|P|²` on the `cosh` channel and `−4|S|²` on the
`sinh` channel. It is **not** `|G(δ+iγ)|²`, which would be `|P+S|²` and is nonnegative.

*Parity split.* For real `g`:

* `g` even → `g cosh(δu)` even, `g sinh(δu)` odd, so `P = ∫ g cosh(δu) cos(γu) du` is real
  and `S = i ∫ g sinh(δu) sin(γu) du` is imaginary;
* `g` odd → `P = i ∫ g cosh(δu) sin(γu) du`, `S = ∫ g sinh(δu) cos(γu) du`.

In both sectors `Q_O(g) = 4[(cosh-channel)² − (sinh-channel)²]`, and in both sectors the
*negative* channel is nonzero as soon as `γ ≠ 0`. At `γ = 0` the orbit degenerates (see
below) and the negative channel survives only in the odd sector — which is the
`dirichlet-windows` real-pair result, recovered here as the `γ → 0` limit of (4.2).

**Collapse.** As `δ → 0`, `A, B → G(iγ)` and `Q_O → 4|G(iγ)|² = 4 H(iγ)` — the
contribution of **two double** on-line zeros at `½ ± iγ`, not of a simple pair
(`2 H(iγ)`). Multiplicity is preserved: four zeros in, four zeros out. **[C-unit]**

**Real off-line pair** (`γ = 0`, `β` real, `β ≠ ½`): the orbit is `{β, 1−β}`, **two**
points, and `Q_O = 2 H(δ)` with `δ = β − ½` real. (Using (4.1) here would double-count.)

**On-line pair:** `{½ ± iγ}`, `Q = 2 H(iγ) = 2|G(iγ)|² ≥ 0`.

---

## 5. Bases

Both sectors use bases that vanish at `±L`, so `g` is continuous on `ℝ`:

```
even:  φ_k(x) = cos((2k+1)πx/2L),    k = 0, 1, 2, …     (parent EVEN_D)
odd:   ψ_k(x) = sin((k+1)πx/L),      k = 0, 1, 2, …     (parent ODD)
```

Both satisfy `∫_{−L}^{L} f_j f_k = L δ_{jk}`, so the Gram matrix in the normalised basis
`f_k/√L` is the identity — **exactly**, not numerically — and the generalised problem
`Bv = λMv` reduces to the ordinary one. (The Gram matrix is still assembled and its
conditioning reported, because the identity is a property of the *full* index set and the
code permits arbitrary index subsets.)

Closed-form transforms, used by `dh.G_at` and unit-tested against quadrature:

```
∫_{−L}^{L} cos(wu) e^{zu} du = 2 [ z cos(wL) sinh(zL) + w sin(wL) cosh(zL) ] / (z² + w²)
∫_{−L}^{L} sin(wu) e^{zu} du = 2 [ z sin(wL) cosh(zL) − w cos(wL) sinh(zL) ] / (z² + w²)
```

Both are entire in `z` (the apparent poles at `z = ±iw` are removable); the code evaluates
them directly and the unit test covers `z` near `±iw`.
