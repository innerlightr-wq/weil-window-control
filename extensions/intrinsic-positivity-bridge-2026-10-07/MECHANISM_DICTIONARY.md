# MECHANISM_DICTIONARY.md

> **NOTE 2026-10-07 — [`SEMILOCAL_LITERATURE_UPDATE.md`](SEMILOCAL_LITERATURE_UPDATE.md).**
> Semi-local Sonin spaces and the relevant structural identifications are established in later work
> (CCM arXiv:2310.18423v2, 2024). The present audit did not locate the specific positive-trace
> comparison with sufficient discrepancy control needed for our application. **This is a gap in the
> verified transfer, not an absence of semi-local geometry.** Row 4's phrase *"open semi-locally"*
> below **remains correct** and is now narrower in meaning: open = **the identity and the estimate**
> (CC2021 Thm 4.6 / Thm 6.11 analogues), **not** the spaces. Two entries are refined, inline:
> the citation of the Ramanujan deposit (see `SOURCES_AND_SCOPE.md`) and the Gram-form row of §2,
> which now records that a **non-isometric** map does not carry a Gram form to a Gram form.
> §3's dictionary (`L = ½log2`, `e^{2L} = 2`, prime activation) and every numerical entry are
> **unchanged** — adding a prime to the finite set `S` is a different operation from enlarging `L`.

## 1. The Ramanujan comparison, in one table

| step | Ramanujan (`thales-ramanujan-quarter-gamma`, Zenodo `10.5281/zenodo.20672985`) | Weil (this repo + Connes–Consani) | status |
|---|---|---|---|
| 1 | analytic factorization (Clausen: a `₂F₁` square equals a `₃F₂`) | arithmetic quadratic form `Q_ζ(g) = Σ_ρ \|ĝ(γ_ρ)\|²` via the explicit formula | both **established** |
| 2 | parameter response (the Euler derivative; quarter-integer lattice `{1/4,1/2,3/4}`) | independently constructed positive Hilbert-space expression `Tr(ϑ(g)Sϑ(g)*) = ‖ϑ(g)S‖²_HS ≥ 0` | both **established**, by different means |
| 3 | independently supplied modular/period identities | explicit discrepancy `E(f)` of Theorem 4.6, identified independently (not by subtraction) | both **established** |
| 4 | exact arithmetic evaluation of `1/π` | an inequality sufficient for the desired positivity statement | **established only at the archimedean place and only on `supp g ⊂ [2^{-1/2},2^{1/2}]`** (CC2021 Thm 1 / 6.11); **open semi-locally** — *clarified 2026-10-07: the semi-local **spaces** and their identification are established (CCM2024 Def. 4.5, Props. 4.6–4.7, Thm 4.6); what is open is the semi-local **identity with an independently defined `E_S`** and the **bound on `E_S`**, which CCM2024 §1 calls a "strategy" it "expect[s]" to work* |

**Does anything actually transfer? No.** The Ramanujan side's load-bearing input is a
hypergeometric/modular identity supplied from outside the factorization; the Weil side's is a
Hilbert-space Gram square. There is no common operator, no common estimate, and no map between the
quarter-integer lattice and anything in the Sonin picture. The repository's own README says the
shared lattice "**is not a mechanism**" and that its contribution is "organizational and
diagnostic". **The analogy has served its organizational purpose and is retired here**: the rest of
this audit is native Sonin/Weil.

Specifically **not** imported: Thales constants, universal thresholds, any "prime Laplacian", and
the Clausen/Legendre material (Euler reflection vs Legendre's relation remains a distinction
internal to that paper).

## 2. Objects kept distinct

| object | what it is | what it does **not** give |
|---|---|---|
| a holomorphic square `H(z)²` | a function that is a square | **does not force real zeros** |
| a norm square `‖Bh‖²` | nonnegative by construction | no information about any *other* form |
| a Gram form `A S A*` | PSD because `= (AS)(AS)*` | needs **no** invariance of `ran S` under `A` |
| the image of an orthogonal projection under a **bounded invertible, non-isometric** map `J` — *added 2026-10-07* | `JPJ^{-1}` is **idempotent**; `JPJ*` is **positive** | **neither is both**, so neither is a Sonin compression with a free sign: `JPJ^{-1}` is not self-adjoint (a trace through it has no sign) and `JPJ*` is not a projection. CCM2024's `θ_S` is exactly such a `J` (their footnote 2: "hilbertian" = topological vector space structure), with critical-line multiplier `∏_{p∈S\{∞}} L_p(½+is)^{-1}` |
| a determinant / Legendre-type identity | a nonvanishing relation | **nonzero determinant ≠ positivity** |
| dissipation under a time evolution | requires a specified dynamical system | **parameter differentiation is not time evolution** |

Also kept in view: **off-line zero quartets retain the reflection symmetry** `ρ ↦ 1−ρ`, `ρ ↦ ρ̄`, so
symmetry of a configuration is no evidence of being on the line.

## 3. The exact dictionary: multiplicative ↔ additive ↔ our form

**Coordinates.** Connes–Consani work on `R*_+` with `d*x = dx/x` and test functions `g ∈ C_c^∞(R*_+)`,
`g*(x) = \overline{g(x^{-1})}`. The repository works in the additive logarithmic coordinate
`u = log x`. The scaling representation `ϑ` becomes translation; Haar measure `d*x` becomes `du`.

**Support.** `g` supported in `[2^{-1/2}, 2^{1/2}]` ⟺ `u ∈ [−½log2, ½log2]`, i.e.

```
      L_CC  =  (1/2) log 2  =  0.34657359...
```

Their Lemma 6.10 is stated on exactly `H = L²(I)`, `I = [−½log2, ½log2]` — the same interval.
The autocorrelation `g*g*` then lives on `[−2L, 2L] = [−log2, log2]`. **This factor of two is
derived from the convolution, not imported**: the repository's own prime term is
`von_mangoldt_terms(2L)` in `src/zeta_window.py`, i.e. `n` with `log n ≤ 2L`.

**Which form is which.** From the explicit formula, with `W_∞ := −W_R`,

```
      Q_zeta(g)  =  Pole(g)  +  W_inf(g*g*)  -  sum_p W_p(g*g*) .
```

So **`W_∞` is *not* the repository's whole `Arch + Pole − log π − Prime`**. The correspondence is

| Connes–Consani | repository (`src/zeta_window.py`) |
|---|---|
| `W_∞(g*g*)` | `Arch − log π · I` |
| — | `Pole` (even sector: `2F(i/2)² ≥ 0`), separate |
| `Σ_p W_p(g*g*)` | `Prime` |

The constant `log 4π + γ` in `W_R` splits as `log 4 + log π + γ`; the repository carries `−log π`
as its own block, which is exactly the `§7` "recorder split" convention already documented there
(`A = Pole + Arch − log π·I`, `2M = Prime`).

**Prime activation.** `W_p` vanishes on test functions supported in the **open** interval
`(p^{-1}, p)`. Hence (verified exactly in `checks.py` §1):

| `L` | `2L` | `e^{2L}` | active `n` |
|---|---|---|---|
| `½log2 = 0.346574` | `log 2` | `2` (exactly) | **none** — `n = 2` sits *on* the boundary, and `g*g*` has support in the open `(1/2,2)` |
| `0.55` | `1.10` | `3.004` | `2, 3` |
| `0.80` | `1.60` | `4.953` | `2, 3, 4` |
| `1.00` | `2.00` | `7.389` | `2, 3, 4, 5, 7` |

> **The Connes–Consani window is archimedean-only.** The first prime activates strictly beyond
> `L = ½ log 2`. Zhu's certified window `L ≤ 0.8`, which the manuscript already uses, **does**
> contain active primes (`n = 2, 3, 4`).

**Test functions, kept apart.** `g` on `R*_+` (multiplicative, their primary object); `h = g*g*`
(its autocorrelation, where the Weil functional is evaluated); and the repository's `f`, real and
even on `[−L,L]` in the additive coordinate. These are three different objects and the support
conventions differ by the convolution factor of two.
