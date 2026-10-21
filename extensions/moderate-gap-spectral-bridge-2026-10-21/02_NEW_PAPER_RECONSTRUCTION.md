# 02 — Reconstruction of Angón–Ross–Tang–Zhao

Everything in this file is **class II** (proved by ARTZ or an antecedent) unless
marked otherwise. Equation numbers are ARTZ's. Recomputed constants carry their
verification tier and the `verify.py` check id.

## 1. The normalised zero multiset `Z_T`

For a nontrivial zero `ρ = β + iγ` with `0 < γ ≤ T` put

```
    L  :=  log T / (2 pi)                                              (ARTZ (2))
    z_rho := i (rho - 1/2) L
```

Writing `β = 1/2 + δ` this is

```
    z_rho = ( -gamma + i delta ) L ,     Re z = -gamma L ,  Im z = delta L .   [I]
```

**This is the first fact that matters for the transfer.** In ARTZ's coordinate the
*ordinate* of a zero is the real part (rescaled by `L`) and its *horizontal
displacement* is the imaginary part (rescaled by the same `L`). A zero is on the
critical line iff `z_rho ∈ R`. The functional equation `ρ ↦ 1 − ρ̄` acts as
`z ↦ z̄`, so `Z_T = Z̄_T` as multisets — which is the only structural hypothesis
Lamzouri's Proposition 2.1 needs.

`Z_T` splits (ARTZ (3)) into
`R_1` (simple critical, `n` elements), `R_2` (multiple critical, `r` distinct),
`R_3` (`k` off-line conjugate pairs), with
`N(T) = Σ_{l≤n+r} m_{x_l} + Σ_{l≤k} 2 m_{z_l}` (ARTZ (5)).

## 2. The window and the kernel

`η_τ ∈ C_c^∞((−1/2,1/2))`, real, even, `∫η_τ² = 1`; `K_τ(ξ) := \hat{η_τ²}(ξ)`,
extended to all of `C` (entire, since `η_τ²` is compactly supported).
ARTZ build `η_τ` from the Montgomery–Taylor window

```
    A(x) = cos( sqrt2 · x ) / ( sqrt2 · sin(1/sqrt2) )   on (-1/2,1/2),   int A = 1
    eta_tau = psi_tau sqrt(A) / sqrt(A_tau) ,  A_tau = int psi_tau^2 A
```

with `ψ_τ ∈ C_c^∞` equal to 1 on `|x| ≤ 1/2 − τ`. As `τ → 0`, `η_τ² → A`.

| quantity | value | check |
|---|---|---|
| `∫_{-1/2}^{1/2} cos(√2 x)dx = √2 sin(1/√2)` | `0.918725369866…` | B1 `[certified]` |
| `Â(0) = 1` | exact | B2 `[certified]` |
| `Â(1/4)` | `0.907203020832…` | B3 `[certified]` |
| `A_τ ≥ sin((1−2τ)/√2)/sin(1/√2) ≥ 1−2τ ≥ 99/100` at `τ=1/200` | `0.9917000761` | B4, B5 `[certified]` |
| `K_τ(1/4) ≥ Â(1/4) − 2/99` | `0.8870010006 > √3/2` | B6 `[certified]` |
| `K_τ(c)² > 3/4` for `0<c≤1/4` | `≥ 0.7867707751` | B7 `[certified]` |

**Correction to a typesetting artefact.** `pdftotext` loses the `1/` inside
`sin(1/√2)` in (9) and mangles the denominators of `Â(c)` in (70). The correct
transform, fixed by `Â(0)=1`, is

```
    Ahat(c) = [ sin(1/sqrt2 + pi c)/(sqrt2 + 2 pi c)
              + sin(1/sqrt2 - pi c)/(sqrt2 - 2 pi c) ] / ( sqrt2 sin(1/sqrt2) ).
```
Lamzouri's `f_0` (his proof of Lemma 3.2) confirms `A = f_0` exactly. `[certified]`

## 3. The real Hilbert space, and where the displacement lives

`f_z(u) := η_τ(u) e^{-2πiuz}`. For `z = ξ + iy`,

```
    f_z(u)  = eta(u) e^{-2 pi i u xi} e^{ 2 pi y u}
    g_z(u)  = (f_z + f_zbar)/2    = eta(u) e^{-2 pi i u xi} cosh(2 pi y u)
    h_z(u)  = (f_z - f_zbar)/(2i) = eta(u) e^{-2 pi i u xi} sinh(2 pi y u) / i
```

`W := span{f_{x_1},…,f_{x_{n+r}}, g_{z_1},…,g_{z_k}, h_{z_1},…,h_{z_k}}`,
`dim W = n + r + 2k`. Every member satisfies `conj(φ(u)) = φ(−u)`, which makes all
inner products real (ARTZ (22) = Lamzouri (2.15)).

**The depth scalar.** Define

```
    kappa_rho  :=  K_tau( 2 i Im z_rho )  =  int eta_tau^2(u) cosh( 2 delta u log T ) du  >= 1 .
```

Then, for a single pair,

```
    ||g_z||^2 = (1 + kappa)/2 ,   ||h_z||^2 = (kappa - 1)/2 ,   <g_z, h_z> = 0 ,
    ||g_z||^2 - ||h_z||^2 = 1 .                                       [IV, certified]
```
(`verify.py` D4–D7; the last identity is ARTZ (19) / Lamzouri (2.11).)
The proof is one line: `cosh² = (1+cosh 2)/2`, `sinh² = (cosh 2 − 1)/2`, and the
cross term is odd. **`κ_ρ` is the only place in the whole construction where the
horizontal displacement of a zero appears.** Everything else — the trace, the
inertia, the block structure — is blind to it.

For the complete Gram, with `Δ = ξ − ξ'`:

```
   <f_x , f_x'> = K(Delta)
   <f_x , g_z'> = Re K(Delta + i y')        <f_x , h_z'> = - Im K(Delta + i y')
   <g_z , g_z'> = ( Re K(Delta + i(y+y')) + Re K(Delta + i(y-y')) ) / 2
   <h_z , h_z'> = ( Re K(Delta + i(y+y')) - Re K(Delta + i(y-y')) ) / 2
   <g_z , h_z'> = -( Im K(Delta + i(y+y')) - Im K(Delta + i(y-y')) ) / 2
```
all real. `[IV, certified D1; corroborated against Simpson quadrature to 5e-15]`

## 4. `G = P + Q`, trace, Frobenius norm, inertia

`α_{ij} := Σ_{z∈Z_T} ⟨f_z,ψ_i⟩⟨f_z,ψ_j⟩`, `G = [α_{ij}]` real symmetric (ARTZ (25)).

```
   G    = sum_{l<=n+r} m_x  v_f v_f^T  +  2 sum_{l<=k} m_z ( v_g v_g^T - v_h v_h^T )   (26)
   tr G = N(T)                                                                          (27)
   ||G||_F^2 = sum_{z,s in Z_T} K_tau(z-s)^2                                            (28)
```

The identity behind (28) is `Σ_i⟨f_z,ψ_i⟩⟨f_s,ψ_i⟩ = ⟨f_z, f_{s̄}⟩ = K(z−s)`
(note: the *bilinear* pairing, not the Hermitian one — these differ off the line,
and `⟨f_z,f_s⟩ = K(z − s̄)`; the two double sums agree because `Z_T = Z̄_T`).
Both (27) and (28) were reproduced numerically on seven synthetic
conjugation-invariant multisets, to 1e-6 relative (`verify.py` F1–F3) `[diagnostic]`.

`P := Σ_{x∈R_1} v_f v_f^T ⪰ 0`, `Q := G − P`, so `tr P = n`, `tr Q = N(T) − n`
(ARTZ (40),(41)) and, by Sylvester's law applied to (44),

```
   n_+(Q) = r + k ,    n_-(Q) = k .                                                     (45)
```

**The eigenvalues behind that inertia count.** For one off-line pair of
multiplicity `m`, isolated in `Q`, the nonzero spectrum of `Q` is

```
   { m (1 + kappa) ,  - m (kappa - 1) } ,   trace 2m ,  ||.||_F^2 = 2 m^2 (1 + kappa^2).
                                                                 [IV, certified E1,E2]
```
The **trace** is `2m` for every `κ`: Sylvester's law and the trace are blind to the
depth. The **Frobenius norm** is not. This gap is the whole content of §04.

## 5. The two analytic inputs

**Lemma 1 / (12) (= Lamzouri Lemma 3.2, resting on BGSTB Lemma 5).** For fixed
`τ ∈ (0,1/4)`,

```
   sum_{rho,rho' : 0<gamma,gamma'<=T}  K_tau( i(rho-rho') log T/(2 pi) )^2
        =  ( C_{eta_tau} + O_tau(1/sqrt(log T)) ) N(T) ,      C_{eta_tau} -> C_MT.
```

`C_η = (η²∗η²)(0) + 2∫₀¹ α(η²∗η²)(α)dα`. With `η² = A`, quadrature gives
`1.327499121` against `C_MT = 1.3274992963205883…` `[diagnostic B8]`.

```
   C_MT = 1/2 + (1/sqrt2) cot(1/sqrt2) = 1.3274992963205883...       [certified, A1]
   C_MT < 4/3                                                        [certified, A2]
```

The underlying object is Montgomery's `F(x,T) = Σ_{ρ,ρ'} x^{ρ−ρ'} W(ρ−ρ')`,
`W(u)=4/(4−u²)`, **with complex differences** — the unconditional form of
Montgomery's theorem of BGSTB, as restated in GS §3:
`F(x,T) ≥ 0`, `F(x,T)=F(1/x,T)`, and
`F(x,T) = (T/2πx²)log²T(1+O(1/√log T)) + (T/2π)log x + O(T√log T)` uniformly for
`1 ≤ x ≤ T`. Lemma 1 is the average `∫_{-1}^{1} f(u) F(T^u,T) du` with
`f = η²∗η²`; the `f(0)` term comes from `x` near 1 and the `2∫₀¹αf` term from
`log x = u log T`. **This identification is reproduced here and matters for §04 §5.**
`[IV, proved]`

**Lemma 9.** `K_τ(c)² > 3/4` for `0<c≤1/4`, `τ ≤ 1/200` (reproduced above).

## 6. The inequality chain

```
   Lemma 4   (p-q)^2 >= Psi(p) - 4q ,   Psi(x) = x^2 (x<=2), 4x-4 (x>2)   [certified C1-C3]
   Lemma 5   ||P+Q||_F^2 >= trPsi(P) + 4 tr(Q) - 4 n_+(Q)                  [class II]
   Lemma 6   trPsi(P) + 2 tr(Q) <= ||G||_F^2
   Lemma 10  trPsi( v_x v_x^T + v_y v_y^T ) = 2 + 2 K(x-y)^2              [certified C4]
   Lemma 11  three close points: trPsi >= 3 + 3K(c)^2                      [certified C5]
   Lemma 12  trPsi(P) >= n + K(c)^2 |Rtilde_1|  and
             tr(Q) + K_tau(c)^2 |Rtilde_1| <= ||G||_F^2 - N(T)             (80)
   Prop 1    w_l(T,c) >= N(T) - ((l+2)/2) tr Q - (l/2)|Rtilde_1| - l
   Thm 1     liminf w_l/N(T) >= 1 - ((l+2)/2)(C_MT - 1)
```

Lemma 11's last step is `4K − 1 ≥ 3K²` on `[1/3,1]`, i.e. `−(3K−1)(K−1) ≥ 0`
`[certified C5]`; Lemma 9 delivers `K ≥ 0.887 ≫ 1/3`.

| ℓ | `(ℓ+4−(ℓ+2)C_MT)/2` | ARTZ |
|---|---|---|
| 2 | `0.3450014074…` | 0.3450014 |
| 3 | `0.1812517592…` | 0.1812517 |
| 4 | `0.0175021110…` | 0.0175021 |
| 5 | `−0.1462484…` | — |

`[certified A3, A4]`. **The ceiling `ℓ ≤ 4` is set by `C_MT < (ℓ+4)/(ℓ+2)`, i.e.
`ℓ < (4−2C_MT)/(C_MT−1) = 4.1069…`, not by the kernel bound** — the condition
`ℓ/(ℓ+2) < K_τ(c)²` is satisfied up to `ℓ = 6`. `[certified A5, A6]`

## 7. What the theorem does **not** say

**A run of `ℓ` consecutive `c`-moderate gaps certifies `ℓ−1`, not `ℓ+1`, simple
critical-line zeros.** `[IV, proved; verify.py J1–J5]`

> *Proposition (anti-localisation).* Let `j ∈ W_ℓ(T,c)`, `ℓ ≥ 2`. Then for every
> interior index `1 ≤ a ≤ ℓ−1` the ordinate `γ_{j+a}` is attained by exactly one
> zero, which is therefore simple and on the critical line. The endpoints
> `γ_j, γ_{j+ℓ}` may be off the line or multiple.
>
> *Proof.* An off-line zero `ρ` has the partner `1 − ρ̄` with the **same ordinate**,
> so an off-line zero always contributes two equal entries to the ordinate list; a
> multiple zero likewise. Both gaps adjacent to an interior index are
> `≥ 2πc/log T > 0`, so no other zero shares `γ_{j+a}`, so the zero there is simple
> and unpaired, hence on the line. For the endpoints only one adjacent gap is
> constrained. ∎

Verified combinatorially on 300 random configurations: 0 interior violations, and
off-line endpoints do occur. The user's caution in the brief is therefore exactly
right, and for exactly this reason — this becomes the decisive obstruction in §05.
