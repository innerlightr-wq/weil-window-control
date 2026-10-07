# COORDINATE_RESPONSE_AUDIT.md

The Weil sensitivity law in every partition coordinate; remainder comparison; monotonicity and
convexity tests. `γ` = spectral ordinate; `L` = window half-width.

---

## 1. The starting point, taken from the current local theorem (not from memory)

**[SOURCE: `paper/weil_window_control.tex` eq. `eq:baseline`, eq. `eq:pert`; Theorem 7;
`notes/proofs.md` §§B.1–B.3]**

```
Q_0        = Q_true + 2 F(gamma)^2                                     (count-preserving baseline)
Q_delta(f) - Q_0(f) = -4 delta^2 ( F'(gamma)^2 + F(gamma) F''(gamma) ) + R_4
A_f(gamma) := F'(gamma)^2 + F(gamma) F''(gamma)
B_L(delta) = c(L) delta^2 + (16/3) L^5 delta^4 e^{2 L delta},   c(L) = 8 L^3 (1/3 + 1/sqrt5)
4 |A_f| <= c(L)   and   |R_4| <= (16/3) L^5 delta^4 e^{2 L delta}
```

and the **exact** functional behind the expansion, with no remainder:

```
Phi(delta) := Q_delta(f) - Q_0(f) = 4 Re F(gamma + i delta)^2 - 4 F(gamma)^2 .
```

`Φ` is even in `δ`, so `Φ = Ψ(δ²)` with `Ψ` entire (`GENERATOR_AUDIT.md` §3). At `L = 0.8`,
`c(L) = 3.1971202`.

**[NUM]** `Φ` against its leading term, `L = 0.8`, `γ = 14.1347`, Simpson `N = 2000`
(`checks.py` §3): ratio `Φ / (−4δ²A_f)` is `1.0000076, 1.0007582, 1.0105784` at
`δ = 0.01, 0.1, 0.375` for `cos(πu/2L)`; `1.0000065, 1.0006454, 1.0090551` for the tent;
`1.0000136, 1.0013591, 1.0192587` for `u sin(γu)`. The quadratic law is the `δ → 0` limit; `Φ` is
the object.

## 2. The leading term in all four coordinates

Substituting `4δ² = χ² = tanh²η = (R−1)/R` (`PARTITION_DICTIONARY.md` §3):

| coordinate | leading term | expansion near the RH point |
|---|---|---|
| `δ` | `ΔQ = −4δ²A_f + R_4` | already quadratic |
| `χ` | `ΔQ = −χ²A_f + O(χ⁴)` | quadratic in `χ` |
| `η` (signed) | `ΔQ = −tanh²(η)·A_f + R̃_4` | `tanh²η = η² − (2/3)η⁴ + O(η⁶)`, so `−η²A_f + O(η⁴)` |
| `R` | `ΔQ = −((R−1)/R)·A_f + R̂_4 = −(1−1/R)A_f + ⋯` | `(R−1)/R = (R−1) − (R−1)² + ⋯`, so **first order** in `R−1` |

The last line is the brief's observation: a quadratic law in `δ` becomes **linear** in `R−1`. It is
also the whole of what the change of coordinates does. `4δ²`, `χ²`, `tanh²η` and `(R−1)/R` are
**the same number**; which power it is written as is a property of the label, not of the form.

**Does any of these simplify the remainder as well, or only the leading term?** Only the leading
term. §3 shows the remainder is unchanged when rewritten exactly and strictly worse when truncated.

## 3. Remainder comparison at the SAME physical displacement  **[NUM, brief control 5]**

`checks2.py` §7, `L = 0.8`:

| `δ` | `B_L(δ)` | via `η`, **exact** | via `R`, **exact** | via `η`, `tanh²η → η²` | via `R`, `(R−1)/R → R−1` |
|---|---|---|---|---|---|
| 0.0500 | 8.0046329e-03 | 8.0046329e-03 | 8.0046329e-03 | 8.0584929e-03 | 8.0856134e-03 |
| 0.1000 | 3.2176288e-02 | 3.2176288e-02 | 3.2176288e-02 | 3.3067800e-02 | 3.3526603e-02 |
| 0.2500 | 2.1000420e-01 | 2.1000420e-01 | 2.1000420e-01 | 2.5660501e-01 | 2.8568766e-01 |
| 0.3750 | 5.1256746e-01 | 5.1256746e-01 | 5.1256746e-01 | 9.6980809e-01 | 1.4749186e+00 |

Maximum deviation between `B_L(δ)` and its exact `η`- and `R`-rewritings across the table: **`0`**.
Truncation ratios against the exact bound: `η`-truncation `1.22190` at `δ = 0.25` and `1.89206` at
`δ = 0.375`; `R`-truncation `1.36039` and `2.87751`. The mechanism is elementary and one-sided:
`tanh²η ≤ η²` and `(R−1)/R ≤ R−1` on the domain.

> **Remainder verdict.** Rewritten exactly, the bound is **the same bound** — no gain, no loss.
> Truncated in `η` or `R`, it is **strictly worse**, by up to `2.9×` at `|δ| = 3/8`. A prettier
> leading term with a worse remainder is not a gain, and that is precisely what truncation in `η`
> or `R` buys.
>
> **Scope, stated 2026-10-07.** This is evidence against **these coordinates for this budget**:
> `B_L(δ) = c(L)δ² + (16/3)L⁵δ⁴e^{2Lδ}` at `L = 0.8` on `|δ| ≤ 3/8`. It is **not** a theorem that a
> change of variables can never improve an estimate or a numerical calculation. No such general
> claim is made or needed.

## 4. Monotonicity transfers; convexity does not  **[NUM + [STD]]**

On `0 ≤ δ < ½` the quantities `δ²`, `χ² = 4δ²`, `η = arctanh(2δ)` and `R = 1/(1−4δ²)` are
**strictly increasing functions of one another** — verified in lockstep over ten nodes
(`checks2.py` §8). Therefore:

> **Monotonicity is coordinate-free.** Any monotonicity statement in one of these coordinates is
> **logically equivalent** to the same statement in the others. **No monotonicity theorem can be
> obtained by changing coordinates.** This is the Two-Face manuscript's own §7.5 point — "all
> candidate symmetric functions of `χ` give identical variance-explained … because they are strict
> monotonic transformations of one another" — in its exact-mathematics form.

**Convexity is the one property a nonlinear monotone change does not preserve**, so it is the only
place a genuine difference could live. It does differ, and that is exactly why it is useless here:

- `x = δ²` is **linear** in `δ²` (so convex and concave);
- the same quantity as a function of `R` is `(R−1)/(4R)`, with `d²/dR² = −1/(2R³) < 0`:
  **concave** for all `R > 0` (`−0.500` at `R = 1`, `−0.0625` at `R = 2`);
- as a function of `η` it is `tanh²(η)/4`, **convex then concave**, with the inflection at
  `tanh²η = 1/3`, i.e. `η = arctanh(1/√3) = 0.6584789`, i.e. **`δ = 0.2886751`** — *inside* the
  physical domain and inside the quasi-RH range `|δ| ≤ 3/8`.

> **Convexity verdict.** A convexity claim about the sensitivity law is **coordinate-dependent**,
> with the `η`-inflection falling inside the physical range. So it is not a statement about the
> Weil form, and **no convexity payoff is available by this route** (brief option C fails).

**And the exact `Φ` has no universal sign or monotonicity anyway.** **[NUM]** `checks2.py` §8 at
`δ = 0, 0.05, 0.1, 0.2, 0.3, 0.375`: `Φ` is negative and decreasing for `cos(πu/2L)`
(`0 → −8.888e-05`) and for `u sin(γu)` (`0 → −1.956e-02`), but **positive and increasing** for
`cos(γu)` (`0 → +6.881e-02`). This matches `A_f`'s sign-indefiniteness (`GENERATOR_AUDIT.md` §7):
positive for three test functions, negative for three others.

## 5. Does the FULL quartet contribution depend on `R` more simply?  **[NO]**

The brief asks for an exact or controlled `Q_R − Q_1 = Φ(R;γ,f)` structurally simpler in `R` than
in `δ`. It does not exist here, and the reason is sharp:

```
Phi = Psi(delta^2) ,   Psi entire ;      delta^2 = (R-1)/(4R)  is RATIONAL in R ;
                                         delta^2 = tanh^2(eta)/4  is MEROMORPHIC in eta .
```

| coordinate | substitution for `δ²` | singularities of the composition |
|---|---|---|
| `δ` | `δ²` | **none** — entire |
| `χ` | `χ²/4` | **none** — entire |
| `η` | `tanh²(η)/4` | **poles at `η = ±iπ/2 + ikπ`** → radius `π/2 = 1.5707963` |
| `R` | `(R−1)/(4R)` | pole at `R = 0` only; **analytic at `R = 1`** |

**An error worth naming, because it is the obvious one and it is wrong.** `δ = ½√(1−1/R)` has a
branch point at `R = 1` — the RH point — and it is tempting to report that as an obstruction. It
is not: because `Φ` is **even** in `δ`, it depends on `R` only through the *rational* `(R−1)/(4R)`,
so `Φ` is single-valued and analytic at `R = 1`. **[NUM]** `checks2.py` §12 exhibits the two
branches of `δ(R)` near `R = 1` and confirms that `Φ` does not see them. No branch point is
claimed.

So the ranking is: **`δ` and `χ` are strictly the best coordinates** (entire); `R` is as good near
the RH point but compactifies the wrong end (`|δ| → ½` becomes `R → ∞`); `η` is **strictly worst**,
with a finite radius of analyticity `π/2` where `δ` has none. All that happens in `R` is the
substitution `δ = ½√(1−1/R)` the brief anticipated — **coordinate translation only.** Again this
ranks *these* coordinates for *this* functional; it is not a general statement about coordinate
changes.

## 6. Coordinate-invariance control  **[NUM, brief control 1]**

Nine combinations (three test functions × `δ = 0.05, 0.2, 0.375`): `Φ` computed from `δ` directly,
from `½ tanh(arctanh 2δ)`, and from `½√(1−1/R)` agree to the printed precision and **all sign
statements are invariant** (`checks2.py` §11). No claimed effect below is an artifact of
reparameterisation — and equally, no claimed effect survives it, because there is none.

## 7. Summary table against the brief's five possible payoffs

| | payoff sought | outcome |
|---|---|---|
| **A** | exact generator bridge: off-line displacement generated by the same operator as intrinsic scaling | **NO.** An exact generator exists — `A =` multiplication by `u`, at parameter `δ` — but it is **canonically conjugate** to the scaling generator, `[A,D] = −I` **on a common core** (see `GENERATOR_AUDIT.md` §5.1 for the domain qualification). Dual, not the same; and a relationship between two operators is not an identification of them. |
| **B** | improved remainder in `η` or `R` | **NO.** Identical if exact, up to `2.9×` worse if truncated (§3). |
| **C** | new monotonicity or convexity in `R` or `\|η\|` | **NO.** Monotonicity is coordinate-free (so nothing is gained); convexity is coordinate-dependent with an inflection inside the physical range (so nothing is a theorem). And `Φ` has no universal sign across `f`. |
| **D** | projected-generator identity for the sensitivity coefficient | **PARTIAL and already in the paper.** `K(L,γ) = ‖AP_Lσ_γ‖²` exactly, and the projection is indispensable since `‖Aσ_γ‖² = ∞` — but this is Theorem 8's own extremal value relabelled, it applies to the derivative-evaluation step only (Remark 10(a)–(c)), and the full coefficient `A_f` is sign-indefinite and so cannot be any such norm. |
| **E** | clean no-go | **YES — this is the result.** Two proved no-gos: the native composition is **additive in `δ`**, not Möbius, so `η`'s additivity is spurious here (and additive composition leaves the partition domain while Möbius cannot); and the two generators are **canonically conjugate**, so the scaling parameter and the displacement parameter are not the same kind of variable. *Amended 2026-10-07: `η` does not break the `δ`-identity — the relabelled family `V_η = e^{−½tanh(η)A}` simply has no constant generator in `η`, so it inherits no addition law.* |
| — | the **positive** outcome | **Off-line displacement is exponential reweighting of the window, not the Sonin scaling action.** The window does not move; the contribution from different positions inside it is reweighted. Coordinate moments enter the sensitivity law because exponential reweighting differentiates into powers of the coordinate. See `GENERATOR_AUDIT.md` §3.1. |
