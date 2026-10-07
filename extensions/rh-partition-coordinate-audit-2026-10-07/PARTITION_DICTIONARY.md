# PARTITION_DICTIONARY.md

The exact `δ ↔ χ ↔ h ↔ η ↔ R` dictionary and its domains. Everything here is **exact algebra**,
verified over `Q` in `checks.py` §1. **Nothing here is a dynamical statement.**

Reminder (see `SOURCES_AND_SCOPE.md`): **`γ` is the spectral ordinate**, never the manuscript's
Lorentz factor; **`L` is the window half-width**, never the manuscript's `χ/2`.

## 1. The partition

A hypothetical zero has `β = ½ + δ`; its functional-equation partner has `1 − β = ½ − δ`. Put

```
a = 1/2 + delta ,      b = 1/2 - delta ,      a + b = 1 .
```

This is the conserved binary partition of the manuscript's eq. (1), and it is the **only**
statement about zeta that this dictionary needs. `a, b ∈ (0,1)` ⟺ `|δ| < ½`.

## 2. The five coordinates

| name | symbol | formula in `δ` | range on `0 ≤ δ < ½` |
|---|---|---|---|
| displacement | `δ` | — | `[0, ½)` |
| asymmetry | `χ = \|2a−1\|` | `2\|δ\|` | `[0, 1)` |
| altitude | `h = √(ab)` | `½√(1−4δ²)` | `(0, ½]`, **decreasing** |
| signed rapidity | `η` | `arctanh(2δ)` | `(−∞, ∞)`; `ξ := \|η\|` |
| response | `R = 1/(4ab)` | `1/(1−4δ²)` | `[1, ∞)` |

## 3. The identities, verified

Exact over `Q` for `δ ∈ {0, 1/10, 1/4, 3/8, −1/5, 49/100}` (`checks.py` §1):

```
a + b = 1 ,        h^2 = a b ,        4 h^2 + chi^2 = 1 ,        R (1 - 4 delta^2) = 1 .
```

To machine precision for `δ ∈ {0, 0.1, 0.25, 0.375, −0.2, 0.49}`:

```
2 delta = tanh(eta) ,     R = cosh^2(eta) ,     R - 1 = sinh^2(eta) ,
delta^2 = (1/4) tanh^2(eta) = (R - 1) / (4 R) .
```

Largest deviation across the table: `1.6e-17` on `δ² = (R−1)/(4R)`, `< 1e-12` on the rest.

Two derived readings, both exact: `h = ½ sech η` and the manuscript's Lorentz factor
`LF = 1/(2h) = cosh η = √R`, so **`R = LF²`**. The Gudermannian leg `θ = gd(η) = arctan(sinh η)`
gives `R = sec²θ` and `χ = sin θ`, which is the manuscript's eqs. (14)–(20) specialised here.
It is recorded for completeness and **is not used** below: no quantity in the Weil problem was
found to be natural in the circular coordinate.

## 4. The RH point

```
RH  <=>  delta = 0  <=>  chi = 0  <=>  eta = 0  <=>  R = 1  <=>  h = 1/2  <=>  a = b = 1/2 .
```

Verified. All five coordinates mark the same point, and each is a strictly monotone function of
every other on the physical domain (`checks.py` §1, `checks2.py` §8) — which is exactly why the
equivalence is automatic and carries no information. The manuscript says this itself (§7.5).

## 5. Domains, kept separate (brief control 6)

| constraint | origin | `δ` | `R` | `\|η\|` |
|---|---|---|---|---|
| `β ∈ (0,1)` | functional equation / trivial strip | `\|δ\| < ½` | `[1, ∞)` | `[0, ∞)` |
| quasi-RH `Re s > 7/8` | the audited OpenAI release, carried as a **domain only** | `\|δ\| ≤ 3/8` | `[1, 16/7]`, i.e. `R ≤ 2.2857143` | `≤ arctanh(3/4) = 0.9729551` |

**`L` is a different variable.** The certified Weil window `L ≤ 0.8`, Zhu's certificate, and the
constants `c(L)`, `K(L,γ)`, `2L³/3` all live in the window half-width. **No statement below mixes
a `δ`-range with an `L`-range**, and no `δ`-coordinate claim is transported into an `L`-claim.

## 6. What this dictionary is, and is not

It **is**: an exact one-dimensional reparameterisation of the real part of a reflected zero pair.
The partition state space is one-dimensional (the manuscript's §2.1), so all five coordinates
carry **identical information**; they differ only in functional form.

It is **not**, and nothing below asserts: that RH is a binary-partition physical system; that
Fisher, Poincaré or Compound metrics govern zeta zeros; that hyperbolic geometry proves RH; that
the manuscript's equipartition-differential results transfer to zero distributions; or that
special relativity and RH share dynamics. The exact available statement is only that
`(½+δ, ½−δ)` is a conserved binary partition of the real parts of a reflected pair. Everything
beyond that is proved or refused in `GENERATOR_AUDIT.md`.

Per the manuscript's own epistemic tiers (§11), the content of this file is **Tier 1 (exact)**:
identities and a dictionary, nothing more.
