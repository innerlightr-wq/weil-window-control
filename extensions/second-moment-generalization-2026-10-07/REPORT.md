# Reciprocal spectral perturbations and the window second moment

Exploration, 2026-10-07. Untracked; nothing in the manuscript, `notes/`, `data/` or any other
tracked file was modified. Repository at HEAD `12fde24`, clean, local == committed.

**Verdict: INTEGRATE.** A correct projected-second-moment statement (Theorem A below) sharpens
Remark 10 of the note by saying exactly when a norm bound *is* the response, but its machinery is
classical and the rank-two identity underneath it is published. No nontrivial extension survives.

---

## 0. Baseline fixed before generalizing

| | |
|---|---|
| Repo | `12fde24`, clean, in sync with `origin/main`; local manuscript == committed (the Retraction-18 revision) |
| `~/Downloads` | holds **both** PDFs: `weil_window_control.pdf` (md5 `e8cbd789…` = Zenodo v1 `23195177`) and `weil_window_control_revised_2026-10-07.pdf` (md5 `a29f9670…` = Zenodo `23212998` = `paper/weil_window_control.pdf`) |
| Audit reused | `extensions/openai-quasi-rh-audit-2026-10-07/` — specifically its separation of a **relative** perturbation estimate from an **independent positivity input** |

Four corrections carried in, all confirmed by computation here:

1. **The coefficient is `/6`, not `/3`.** For a single reciprocal real pair, `a = log rho`,
   `lambda_min(T_R) = -binom(R+2,3) a^2 + O_R(a^4) = -R(R+1)(R+2)a^2/6 + O_R(a^4)`.
   Verified exactly for `R = 1,2,4,8,16,32` (`s2_baseline.py` §2.3): `2 V_mu = binom(R+2,3)` with
   `V_mu = R(R+1)(R+2)/12` as an exact rational.
2. **`K(L,gamma)` is a norm, not an eigenvalue.** Theorem 8 is a *conditional lower bound*; there
   is no universal equality `lambda_min = -4K delta^2`. Quantified in §5 below.
3. **The exponents 2 and 3 are in different variables** — displacement `a` and window size `L`.
   Nothing here is a quadratic equation becoming cubic.
4. **No positivity assumption is dropped when a window is extended.** Every statement below is a
   statement about `lambda_min` of an explicit finite matrix, never about the ζ form.

---

## 1. The rank-two baseline — correct, and classical

For a finite positive measure `mu` on `R` of bounded support, on `L^2(mu)`,
`(T_a f)(x) = int 2 cosh(a(x-y)) f(y) dmu(y)`. Since
`2cosh(a(x-y)) = e^{ax}e^{-ay} + e^{-ax}e^{ay}`, writing `u = e^{ax}`, `w = e^{-ax}` gives
**`T_a = u (x) w + w (x) u`**, of rank at most two, and the standard 2x2 reduction gives

> `lambda_pm(a) = <u,w> +/- sqrt(<u,u><w,w>) = M +/- sqrt( (int e^{2ax} dmu)(int e^{-2ax} dmu) )`.

**This is published.** Alfakih, *On yielding and jointly yielding entries of Euclidean distance
matrices*, [arXiv:1609.07055](https://arxiv.org/pdf/1609.07055), **§2.2, Proposition 2.2, p. 7**,
verbatim: *"Let `a` and `b` be two nonzero, nonparallel vectors in `R^r`, `r >= 2`, and let
`Psi = ab^T + ba^T`. Then `Psi` has exactly one positive eigenvalue `lambda_1` and one negative
eigenvalue `lambda_r`, where `lambda_1 = a^T b + ||a|| ||b||` and `lambda_r = a^T b - ||a|| ||b||`."*
Same 2x2 proof. **So §1 is calibration, not a theorem.**

**Expansion.** With `X ~ mu/M`, odd cumulants cancel between `phi(2a)` and `phi(-2a)`, giving

> `lambda_-(a) = -2 V_mu a^2 - (2/3) V_4 a^4 + O(a^6)`,  `V_mu = int (x-xbar)^2 dmu`, `V_4 = int (x-xbar)^4 dmu`.

Both terms are `<= 0`. Verified to 40+ digits on five measures (`s2_baseline.py` §2.2): the residual
`(lambda_- + 2V_mu a^2)/a^4` converges to `-(2/3)V_4` in every case.

**Why a *centred* moment appears even off-centre.** The kernel depends only on `x-y`, so translating
all nodes leaves every matrix entry — hence `lambda_-` — unchanged. The coefficient must therefore be
translation invariant, and the only second moment that is, is the centred one. Confirmed numerically
on `[10,12]` (`s2_baseline.py`, "translated interval"): coefficient `2 V_mu` with `V_mu` centred.

**Degenerate case, stated separately.** By Cauchy–Schwarz
`sqrt(<u,u><w,w>) >= |<u,w>| = M`, with equality iff `u || w` in `L^2(mu)`, i.e. iff `e^{2ax}` is
`mu`-a.e. constant. Hence **`lambda_- < 0` strictly iff `a != 0` and `supp mu` has at least two
points.** For a one-point measure `T_a` is rank *one* with spectrum `{2M}`: the "`-`" branch is not
an eigenvalue at all, and `V_mu = 0` consistently. (This is why the one-point row of §2.1 is a
*degeneracy*, not a mismatch.)

**Explicit remainder.** With `R_mu = max |x - xbar|` on `supp mu` (`<= diam`),
`C - 1 = E[cosh(2aX) - 1] <= 2a^2 m_2 cosh(2aR_mu)` and `sqrt(C^2-S^2) <= C` give

> **`lambda_-(a) >= -2 V_mu a^2 cosh(2 a R_mu)`.**

Checked at `a = 0.01, 0.1, 0.5` on every non-degenerate case. **The control parameter is `a R_mu`,
not `a`** — a fixed-window `O(a^4)` is *not* uniform in a growing window, exactly as cautioned.

**Calibrations.** (a) counting measure on `{0..R}` -> `binom(R+2,3)` exactly. (b) Lebesgue on
`[-L,L]` -> `lambda_-(a) = 2L - sinh(2aL)/a` exactly, with continuous extension `0` at `a=0`, and
series `-(4/3)L^3 a^2 - (4/15)L^5 a^4 - ...` matching `-2V_mu = -4L^3/3`. (c) nonuniform weights,
irregular nodes, translated interval, one point — all pass.

*Representation note:* we always work in the Euclidean orthonormal frame
`T_ij = sqrt(w_i w_j) 2cosh(a(x_i-x_j))`. This is **not** Toeplitz unless the nodes are equispaced
*and* the weights constant, and we never call it Toeplitz otherwise.

---

## 2. What is geometric and what is normalization

With mass `M ~ L^s` and *normalized* variance `V_mu/M ~ L^2`, the coefficient is
`2 V_mu = 2 M (V_mu/M) ~ L^{s+2}`. Measured log-log slopes (`s3_scaling.py` §3.2):

| construction | mass | predicted `s+2` | measured slope |
|---|---|---|---|
| Lebesgue on `[-L,L]` | `~L^1` | 3 | **3.00000** |
| probability-normalized on `[-L,L]` | `~L^0` | 2 | **2.00000** |
| 5 fixed nodes stretched by `L` | `~L^0` | 2 | **2.00000** |
| Lebesgue `x L^2` | `~L^3` | 5 | **5.00000** |

> **Cubic window growth is not a law; it is the case `s = 1`.** It requires the mass to grow
> linearly with the window. Probability normalization gives `L^2` instead.

Both of the note's settings have `s = 1`: counting measure on `{0..R}` has mass `R+1` (coefficient
`~R^3/6`, measured `coeff/R^3 -> 1/6`), and Lebesgue on `[-L,L]` has mass `2L` (coefficient
`4L^3/3`). That — not a universal principle — is what the dictionary's shared
`(window length)^3/6` records.

*Measure scaling vs inner-product scaling.* `mu -> c mu` scales the map **and** the inner product,
so `lambda -> c lambda` and `V_mu -> c V_mu`, consistently (checked: ratio exactly 7 for `c=7`).
Multiplying only the inner product by `c` with the linear map fixed instead gives `lambda -> lambda/c`.
Different operations; not interchangeable.

*`N` is not `L`.* At fixed `L`, refining the quadrature `N = 25..400` converges to
`2L - sinh(2aL)/a` and changes no exponent (`s3_scaling.py` §3.4). The `L`-power comes from the
measure.

---

## 3. THE ONE THEOREM: a projected second moment

> **Theorem A (projected second-moment response).** Let `mu` be a finite positive measure on `R`
> with finite support of `n` points and weights `w_i > 0`; work in `L^2(mu)` (dimension `n`). Let
> `1` and `x` denote the constant and coordinate functions. Let `G = G^* >= 0` be fixed and put
> `A_a = G + T_a`, so `A_0 = G + 2|1><1| >= 0`. Let `E = ker A_0` and suppose `E != {0}`. Then:
>
> 1. `E` is contained in the mean-zero subspace: `v in E => <1,v> = 0` and `G v = 0`.
> 2. `T_a = T_0 + a^2 Q + O(a^4)` with `Q = |x^2><1| + |1><x^2| - 2|x><x|`, and the compression of
>    `Q` to `E` is exactly
>    **`P_E Q P_E = -2 |P_E x><P_E x|`.**
> 3. Consequently, if `P_E x != 0`,
>    **`lambda_min(A_a) = -2 a^2 ||P_E x||^2 + O(a^4)`**,
>    the `O(a^4)` depending on the spectral gap of `A_0` on `E^perp` and on `a R_mu` as in §1.
> 4. If `P_E x = 0` the quadratic term vanishes identically.
> 5. If instead `A_0 > 0` with lowest eigenspace `E_0 ⊥ 1`, the same compression gives
>    `lambda_min(A_a) = lambda_0 - 2a^2||P_{E_0}x||^2 + O(a^4)`, which stays **positive** for
>    `a^2 < lambda_0 / (2||P_{E_0}x||^2)`.
>
> *Proof of (2):* for `v in E`, `<v,Qv> = <v,x^2><1,v> + <v,1><x^2,v> - 2|<x,v>|^2 = -2|<x,v>|^2`,
> the first two terms vanishing by (1). (3) is first-order degenerate perturbation theory in the
> parameter `a^2` (Rellich/Kato). ∎

**When does the coefficient equal `V_mu`?** Exactly when `P_E x = x - xbar 1`, i.e. when the
background imposes no constraint beyond mean-zero — in particular when `G = 0`.
**Otherwise the background strictly reduces it to a projected moment.**

Verified (`s4_background.py`, 4 nodes, exact rational projections, Decimal eigenvalues):

| case | `dim E` | `||P_E x||^2` | predicted coeff | measured `(λ-λ₀)/a²` at `a=1e-4` |
|---|---|---|---|---|
| (A) `G = 0` | 3 | `10 = V_mu` | `-20` | `-20.0000002` |
| (B) `G = 3|v><v|`, `v` generic | 2 | `49/5` | `-19.6` | `-19.6000002` |
| (C) `G = 3|x_c><x_c|` | 2 | `0` | `0` | `-8.0e-17` (response is `-0.8 a^6`) |
| (D) `G = 5 P_{1^perp}`, `A_0 > 0` | — | `10` | `-20` | `-20.0000002`, `λ_min = 4.99999980 > 0` |

Case (C) is the sharpest check: the quadratic *and* quartic terms vanish and the response is `a^6`.
The quartic cancellation is explicit — the first-order `a^4` term of `T_a` contributes
`+(1/2)|<x^2,v>|^2 a^4` and second-order perturbation theory in `a^2 Q` contributes
`-(1/2)|<x^2,v>|^2 a^4`.

**Status.** Classical machinery (Rellich/Kato degenerate perturbation theory) plus a two-line
compression. Correct and useful; not a discovery.

**Infinite dimensions: NOT attempted.** Theorem A is a finite-matrix statement. Extending it needs
form domains, closedness/lower bounds and a spectral-gap hypothesis for the *unbounded* Weil form;
a finite-matrix proof is not a theorem about it. No such extension is claimed.

---

## 4. Oscillatory extension, and what it says about `K`

For the real conjugation-closed quartet `4cos(gamma(x-y))cosh(a(x-y)) = 4 Re cosh((a+i gamma)(x-y))`:
`A_0 = 4(|c><c| + |s><s|) >= 0` with `c = cos(gamma x)`, `s = sin(gamma x)`, and the `a^2` term has
kernel `2(x-y)^2 cos(gamma(x-y))`. For `v ⊥ c,s` the bare-`c(y)`/`s(y)` terms die and

> **`P_E Q_gamma P_E = -4[ |P_E(xc)><P_E(xc)| + |P_E(xs)><P_E(xs)| ]`.**

Confirmed to 8 digits (`s5_oscillatory.py` §5.1): predicted `-1.26283476` vs measured `-1.26283474`
(`L=1, gamma=3`); predicted `-5.81875010` vs measured `-5.81874999` (`L=1.6, gamma=5`).

**A phase that is merely unitary conjugation is a control, not a result** — it moves `|F'|` around
without changing the spectrum. What matters here is that the *real* quartet is rank 4 and its
compression is a **sum of two projected moments**, one per parity.

**`K(L,gamma) = L^3 k(gamma L)`,** `k(z) = int_{-1}^{1} t^2 sin^2(zt) dt`. Measured (`§5.2`):

| `z = gamma L` | 0.01 | 0.1 | 1 | 10 | 200 |
|---|---|---|---|---|---|
| `k(z)` | 4.000e-5 | 3.990e-3 | 0.3141 | 0.2859 | 0.3355 |
| `k(z)/z^2` | 0.39999 | 0.39905 | — | — | — |

> **`k(z) -> 1/3` as `z -> infinity`**, so `K ~ L^3/3`: the cubic law is the **large-`gamma L`**
> regime. **`k(z) ~ (2/5) z^2` as `z -> 0`**, so small `gamma L` gives **`K ~ (2/5) gamma^2 L^5`**,
> not `L^3`. And `k(0) = 0` exactly: at `gamma = 0` the quartet degenerates to a double pair and the
> sensitivity vanishes. `k` is not monotone (`k(2) = 0.581 > 1/3`), matching the note's observation
> that the `sin^2` correction is non-monotone in `L`. **Replacing `sin^2` by `1/2` is legitimate
> only in the averaging regime `gamma L >> 1`.**

**The (a)/(b)/(c) distinction, made quantitative.** In the real-even sector `x cos(gamma x)` is odd
and `x sin(gamma x)` is even, so only the latter contributes, and its *unprojected* norm squared is
exactly `K(L,gamma)`. But `<x sin(gamma x), cos(gamma x)> != 0`, so the projection strictly reduces
it. **Theorem A therefore says precisely when `4K` is the response and when it is only a bound:
the response is `4||P_E(x sin(gamma x))||^2 <= 4K`, with equality iff `x sin(gamma x) ⊥ range(A_0)`.**

**Honest negative.** In this toy model the reduction is tiny — `||P_E(x sin)||^2 / K` is
`0.9966, 0.9963, 0.9993, 0.9997, 0.9991` at `L = 0.8, 1.0, 1.3, 1.6, 2.0`. The note's measured
`C/4K` is `0.419, 0.635, 0.802, 0.963, 0.947`. **So this model does NOT explain the note's gap.**
The gap must come from the real baseline — archimedean, pole and prime terms and all remaining
zeros — not from the on-circle quartet alone. The qualitative lesson (a norm identity for
`F'(gamma)` bounds, and does not compute, the ground-state response) stands; the quantitative
explanation does not. **Any bridge to the ζ theorem must still carry its remaining-zero hypothesis.**

---

## 5. Why "symmetric spectral problem" is not enough

| | |
|---|---|
| (a) | `A(a) = diag(a,-a)`: `spec A(a) = spec A(-a)` as a **set**, yet `A(a) != A(-a)`. The analytic branches `±a` are smooth; their **ordered minimum** `-|a|` is only Lipschitz — linear, not quadratic. Spectral symmetry ≠ even family; ordered min ≠ analytic branch. |
| (b) | `A(a) = diag(1+a^2, 2)`: even, analytic, `lambda_min` **increases** quadratically. |
| (c) | Vanishing quadratic term with higher-order response: §3 case (C), measured `-0.8 a^6`. Also trivially `diag(-a^4,1)`. |
| (d) | Strictly positive baseline: §3 case (D), `lambda_min(A_0) = 5`, so positivity survives all `|a| < 0.5`. Expanding about a positive baseline ≠ creating a negative eigenvalue. |
| (e) | Same support, different mass/variance: Lebesgue on `[-1,1]` gives coeff `1.333`; the same support with mass `x3` gives `4.000`; two atoms at `±1` of mass 2 give `4.000`, of mass 20 give `40.0`. Support length fixes neither the coefficient nor its `L`-power. |
| (f) | Reparametrisation: `a = c delta` multiplies the coefficient by `c^2` (the note's `a = delta log q`); `a = L delta` **adds 2** to the apparent window power (`L^3 -> L^5`). Fix physical parameters before comparing exponents. |

---

## 6. Prior art, with locations

| Item | Status |
|---|---|
| Rank-two eigenvalues `a^T b ± ||a|| ||b||` | **FOUND, classical.** Alfakih, arXiv:1609.07055, §2.2 **Proposition 2.2, p. 7**, quoted in §1 above; same 2x2 proof. Verified by reading the PDF text. |
| Rellich's theorem (analytic Hermitian families) | **CONFIRMED as the classical antecedent** via Barbarino–Noferini, arXiv:2211.15539, abstract: *"generalizes the celebrated theorem of Rellich for matrix-valued functions that are analytic and Hermitian on the real line."* The **introduction's exact statement NOT CHECKED** — the abstract page did not expose it. |
| Greenbaum–Li–Overton, arXiv:1903.00785 | **Read (abstract).** First-order only, general (non-Hermitian) matrices, simple eigenvalues; multiple eigenvalues only referenced. **Does not contain Theorem A.** |
| Kato, *Perturbation Theory for Linear Operators* (degenerate Rayleigh–Schrödinger) | **NOT CHECKED** — not available here. This is the standard reference for Theorem A(3) and I do not claim more than that. |
| Paley–Wiener derivative-evaluation norm | **NOT FOUND in sources checked.** Searches returned the reproducing kernel `sin(eta(u-v))/(pi(u-v))` but no explicit derivative-evaluation norm. `K(L,gamma) = ||u sin(gamma u)||^2` is the Riesz-representer computation and is classical in substance. **"Not found in sources checked" is not novelty.** |
| Moment matrices, reciprocal Toeplitz perturbations | **NOT CHECKED** systematically. |

Fisher-information/variance readings were **not** pursued: no exact statistical model was supplied,
and `V_mu` here is an unnormalized measure moment, not a variance of a probability law unless `M=1`.

The prime-comb and coupling-constant experiments were **not** rerun.

---

## 7. Answers

**1. What explains the quadratic dependence?** `T_a = u (x) w + w (x) u` with `u = e^{ax}`,
`w = e^{-ax}`, and `lambda_- = M - sqrt(<u,u><w,w>)`. At `a = 0`, `u = w = 1` and Cauchy–Schwarz is
an equality; `a` tilts `u` and `w` apart, and the Cauchy–Schwarz *defect* is second order in the
tilt. The coefficient is the centred second moment because the kernel depends only on `x-y`, so the
answer must be translation invariant.

**2. What explains cubic window scaling, and how does normalization change it?** It is
`coeff = 2V_mu ~ L^{s+2}` where `L^s` is the mass. Cubic is the single case `s = 1` — Lebesgue
measure on `[-L,L]` (mass `2L`) and counting measure on `{0..R}` (mass `R+1`). Probability
normalization gives `L^2`; a fixed node set stretched gives `L^2`; mass `~L^3` gives `L^5`.
**Normalization-dependent, not universal.**

**3. When does the coefficient measure the actual response rather than only a bound?** Theorem A:
the response is `2||P_E x||^2` where `E = ker A_0`. It equals the full moment `V_mu` iff the
background imposes nothing beyond mean-zero; it is strictly smaller when the background constrains
`x`; it vanishes when `P_E x = 0`; and it does not create a negative eigenvalue at all when
`A_0 > 0` and `a` is below the gap. In the oscillatory case the response is
`4||P_E(x sin(gamma x))||^2 <= 4K(L,gamma)`, with equality iff `x sin(gamma x) ⊥ range(A_0)`.

**4. Does anything useful and apparently new survive beyond the rank-two identity?** Useful: yes —
Theorem A and its oscillatory analogue give the exact criterion separating "norm bound" from
"response", which is what Remark 10 of the note asserts qualitatively. New: **no.** The rank-two
identity is published (located above), and Theorem A is degenerate perturbation theory plus a
two-line compression. The one quantitative hope — that the projection explains the note's
`C/4K` gap — was **tested and failed** (§4).

---

## 8. Verdict

**INTEGRATE.** The correct broader statement is the *projected* second moment, with cubic window
growth demoted to the normalization case `s = 1` and the `K ~ L^3` law demoted to the regime
`gamma L >> 1`. It strengthens the note's Remark 10 by making the (a)/(b)/(c) distinction exact and
checkable. It is largely classical, rests on a published rank-two identity, and does not explain the
measured saturation gap. There is no credible route from here to a nontrivial new theorem, so this
is not PURSUE; and because a usable sharpening of Remark 10 does survive, it is not CLOSE.

**Central question, answered.** There *is* a useful projected-second-moment theorem for reciprocal
spectral perturbations (Theorem A), and cubic window growth *is* a specified normalization-dependent
case (`mass ~ L^1`) rather than a universal `delta^2 L^3` principle.

## Reproduce

```
python3 s2_baseline.py       # rank-two identity, expansion, remainder, calibrations
python3 s3_scaling.py        # the L^{s+2} law, measure vs inner product, N-convergence
python3 s4_background.py     # Theorem A, cases (A)-(D)
python3 s5_oscillatory.py    # quartet compression; K = L^3 k(gamma L); projection test
python3 s6_counterexamples.py
```

`kernel2.py` carries two independent routes — the closed form and a general Jacobi eigensolver that
assumes nothing about rank — in 60-digit `Decimal` (no `mpmath` on this machine). **Stable digits
are not a rigorous enclosure**, and nothing here is offered as one.

**Two tooling bugs were found and fixed during this work, both recorded in the scripts:** the first
Jacobi implementation updated the `(p,q)` block in place from already-rotated columns, silently
breaking trace preservation (caught by a PSD control); and slicing `str(Decimal)` for display drops
the exponent, which made a `1e-25` eigenvalue read as order 1 and briefly produced a wrong reading
of case (C).
