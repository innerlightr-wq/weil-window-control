# CLAIMS.md — every claim, classified

**Classes.** `PROVED` = an analytic or exact-rational derivation given in `DERIVATION.md`
(sub-labelled `exact/Q` where it is checked by exact rational arithmetic, which is stronger
than a numerical check). `CERTIFIED` = an Arb ball enclosure with a reproducible
certificate (`src/certify.py`, 400 bits). `NUMERICAL` = high-precision computation with a
precision-doubling or scale-free residual check, **uncertified**. `CONJECTURAL` = stated,
not established. `REFUTED` = shown false by an explicit counterexample.

No row is `CERTIFIED` unless an interval enclosure exists for *that* statement. No row is
`PROVED` on the strength of a computation.

---

## Part 2 — the mean-gap diagnostic

| # | claim | class | where |
|---|---|---|---|
| M1 | `A == 1/2` identically; `G^2 == A M_H`; hence `M_H = 2G^2` and the mean triple is a function of the **single** scalar `D = A - M_H = 2 delta^2` | **PROVED** (exact/Q) | §1.1 |
| M2 | `A - M_H = D`; `A - G = 1/2 - sqrt(1/4 - D/2)`, strictly increasing and convex in `D` | **PROVED** (exact/Q) | §1.3 |
| M3 | `G - M_H` is strictly concave in `D`, vanishes at `D = 0` (equality) **and** at `D = 1/2` (complete dominance), with maximum `1/8` attained exactly at `D = 3/8`, where `G = 1/4`; hence **two-to-one** onto `(0, 1/8)` | **PROVED** (exact/Q) | §1.3 |
| M4 | along any evolution: `(A-M_H)' = D'`, `(A-G)' = D'/(4G)`, `(G-M_H)' = D'(1 - 1/(4G))`; the last **flips sign** at `D = 3/8` | **PROVED** | §1.4 |
| M5 | `M_H'' = 1 - 4 M_H` in the arc coordinate `delta = (1/2) sin theta`, and it is **not** dynamics: every `c_0 + c_1 cos 2theta` satisfies `f'' = -4(f - c_0)` | **PROVED** + negative control | §1.5 |
| M6 | for `\|delta\| > 1/2` the partition reading fails (`G` imaginary, `M_H < 0`); `sup D^zeta = 1/2`, `sup D^DH = 1.6025374835598227777878195 > 1/2` | **PROVED** (exact/Q) / **CERTIFIED** (DH value) | §1.6, certify C4 |
| M7 | `tau_c(F) = Delta^2/2 = sup_rho D(rho) = sup_rho (A - M_H)(rho)`: **de Bruijn's threshold is exactly the maximal arithmetic–harmonic mean gap**. For DH it is `1.602537483559822777787819`, reproducing the parent's committed value to all 22 digits | **PROVED** (dictionary identity) / **CERTIFIED** (DH value) | §1.7 |
| M8 | the means are **not new coordinates**: the committed partition audit's `h` *is* `G` and its `R` *is* `1/(2 M_H)`; its verdict MISMATCH and its "monotonicity is coordinate-free" finding apply verbatim to `A - G` and `A - M_H` | **PROVED** (exact/Q re-check) | §1.2 |

## Part 3 — the specified evolution

| # | claim | class | where |
|---|---|---|---|
| D1 | the dictionary: `H_0^zeta = xi(1/2+iz/2)/8`, `H_0^DH = Lambda_F(1/2+iz/2)/2`, `d_tau H = -d_z^2 H`, `delta = -Im z/2`, `D = (Im z)^2/2`, `kappa` from the Gauss sum | **NUMERICAL** (residuals `<= 1.3e-40` relative, dps 90) | §0.3 |
| F1 | `z'(tau) = H''(z)/H'(z)` and `D'(tau) = Im(z) * Im(H''/H')` on a simple branch | **PROVED** (exact chain rule) | §2.1 |
| F2 | `H_tau` even of order 1 ⟹ `h_tau(w) = H_tau(sqrt w)` has order `1/2` ⟹ **genus 0** ⟹ absolutely convergent canonical product with **no** exponential factor; `H'/H(z) = 2z sum_k 1/(z^2-a_k^2)` | **PROVED** | §3.1–3.2 |
| F3 | `H''/H'(a_j) = 2 sum_{k != j} 1/(a_j - z_k)` over the full symmetric zero list, in the **paired** order (the unpaired sum is not absolutely convergent) | **PROVED** | §3.3 |
| F4 | `D'(tau) = -2 y_j sum_{k != j} (y_j - y_k)/\|z_j - z_k\|^2`, **absolutely convergent** given a bounded zero strip | **PROVED** | §3.4 |
| F5 | own-orbit term `= -1 - y_j^2/\|z_j\|^2` (`Re z_j != 0`) or exactly `-1` (`Re z_j = 0`); the conjugate partner contributes exactly `-1` and the partner `-conj z_j` exactly `0` | **PROVED** | §3.5 |
| F6 | a zero contributes **positively** to `D'` only if it lies strictly further from the critical line **on the same side** | **PROVED** | §3.5 |
| F7 | if `\|Im z_j\|` is maximal over the zeros then `D'(tau) <= -1`; with a margin `eps` below the supremum, `D' <= -1 + 2\|y_j\| eps S_j` | **PROVED** | §3.5 |
| F8 | (F7) is the differential form of de Bruijn's bound and integrates to `Delta_tau^2 <= max(Delta_0^2 - 2 tau, 0)` | **PROVED as a differential statement; the envelope argument through collisions and the attainment of the supremum are NOT carried out** — no new proof of de Bruijn's theorem is claimed | §3.5 |
| F8b | **(S1d)** a single competing zero can lift `D_j'` above `-1` only if `R = \|Re z_j - Re z_k\| < Delta_tau/sqrt3`, **sharp**, attained exactly at `(y_j,y_k,h) = (2Delta/3, Delta, Delta/3)`; crude `N`-competitor form `R < Delta_tau sqrt(N/2)` | **PROVED** (exact/Q); **SUPERSEDES** the `Delta/sqrt2` and `Delta sqrt(N)/2` of the first pass — the latter was **wrong** | AUDIT §1.1 |
| F8c | **(S1d')** the sharper **per-zero** form: `R < y_j`, i.e. `\|gamma_j - gamma_k\| < \|delta_j\|`, necessary in **all** cases (AM–GM covers the strip-truncated branch). Necessary, **not** sufficient (7 of 19 grid points at `R = 0.9 y_j` admit no competitor) | **PROVED**; grid-verified `19x19x200`, zero violations | AUDIT §1.1 |
| F8e | **single vs several:** F8b and F8c bound **one** competitor's contribution. For `N` competitors the necessary condition relaxes to `R < min(N y_j, Delta_tau sqrt(N/2))` — linearly and like `sqrt N`. **No conclusion here rests on `N = 1`:** the toy counterexample has `N = 1`, and every DH census target is excluded even by the `N = 3` form | **PROVED** + **NUMERICAL** | AUDIT §1.1 |
| F8d | consequences: the minimal counterexample is confluent in the ordinate (`R = 0`); no DH census quartet can be one (competitors outside the per-zero cutoff by **189, 1086, 406**); for `zeta` the mean spacing crosses the `(S1d)` threshold at `gamma* = 2 pi exp(4 pi sqrt2) = 328290294.5`, Lehmer-type pairs being inside far earlier | **PROVED** + **NUMERICAL** | AUDIT §1.1 |
| F9 | `Escr(tau) := sum_j D_j = (1/2) sum_j (Im z_j)^2` satisfies `Escr' = -2 sum_{j<k} (y_j-y_k)^2/\|z_j-z_k\|^2 <= 0`, with equality iff all `Im z_j = 0` | **PROVED** exact/Q for polynomial flows; **PROVED conditional on `m < infinity`** for `H_tau`. **CORRECTED:** the first pass gave the hypothesis as `Escr < infinity`, which is not sufficient and with `m = infinity` is *incompatible* | §3.6, §4.1, AUDIT §1.3 |
| F10 | each conjugate pair contributes exactly `1` to the dissipation, so `Escr' <= -m(tau) <= -2` while `Escr > 0`, hence `Lambda <= Escr(0)/2` | **PROVED** conditional on `Escr(0) < infinity` **and** absolute continuity of `Escr` on `[0,Lambda]` (stated as a hypothesis, not derived — it is the step where the infinite case differs from the polynomial one) | §3.7 |
| F11 | `Escr(0)/2 >= Delta_0^2/2` always, with equality iff the non-real zeros are a single purely imaginary pair — so **de Bruijn's bound is always at least as strong as (F10)** | **PROVED** | §3.7 |
| F12 | collisions: a real simple branch stays real; a branch leaves the axis only through a collision; the symmetric functions of a colliding cluster are analytic through it, the individual `D_j` need not be | **PROVED** | §3.8 |

## Part 4 — the controls

| # | claim | class | where |
|---|---|---|---|
| C1 | **universal monotonicity of the individual `D` is false**: for `P_0 = (z^2+1)(z^2+4)`, `D'(0) = +1/3` at `z = i` while `D'(0) = -19/3` at the extremal `z = 2i` | **REFUTED** (exact/Q, three independent derivations) | §4.1 |
| C2 | exact threshold for `P_0 = (z^2+a^2)(z^2+b^2)`: inner-pair `D'(0) = 4b^2/(a^2-b^2) - 1 > 0` iff `a^2 < 5b^2`, and `= 0` exactly at `a^2 = 5b^2` | **PROVED** (exact/Q) | §4.1 |
| C3 | minimality: a single symmetry orbit always has `D' < 0`, so two orbits (degree 4) is the minimal counterexample | **PROVED** | §4.1 |
| C4 | toy collision structure: `tau_* = (10-sqrt46)/24`, `tau_** = (10+sqrt46)/24`; `Escr` continuous, `Escr'` **jumps** `-12 -> -8`; `Escr(tau_**) = 0` exactly | **PROVED** (exact/Q) | §4.1 |
| C5 | `sum_j Re(z_j^2)` increases at exactly `2n(n-1)` for **every** polynomial heat flow — monotone, universal, and **not** a displacement diagnostic (it grows even when all zeros are real) | **PROVED** (exact/Q) | §4.1 |
| H1 | DH: `D'(0) < 0` at all four census orbits — `-1.16044463430610822178`, `-1.04435944844102968031`, `-1.01207352431659836488`, `-1.12778399346959232977` | **CERTIFIED** (Arb 400 bits, radii `<= 4.5e-32`), **as a function of the supplied `z_1`** | §4.2 |
| H2 | DH: the cross-orbit remainder at `O_1` is `<= -0.1604316745137291191935` | **CERTIFIED** | §4.2 |
| H3 | DH: the partial sums of (F4) over an independently computed zero list converge to the exact cross-orbit term; remainder `0.14 %` of `\|D'\|` at cutoff `\|Im s\| < 260` | **NUMERICAL** (measured residual, explicitly **not** a bound) | §4.2 |
| H4 | DH: an independent scan finds **122** on-line ordinates below `T = 200`, and `122 + 2x4 = 130` — the parent extension's contour count and reconciliation, reproduced from separate code | **NUMERICAL** | §4.2 |
| H5 | DH: the off-line quartet at `gamma_1` sits in a zero gap of `3.05` mean spacings, so a mean-spacing density heuristic **overestimates** the drag by a factor `8.12` (`-1.3026848` predicted vs `-0.16043167` actual) and is not usable as a dissipation lower bound | **NUMERICAL** | §4.2 |
| H6 | global extremality of `O_1` is **not** established and no finite census can establish it; for `O_2, O_3, O_4` the extremality hypothesis is **false already inside the census**, so their negative `D'` is a numerical fact, not a forced one | **PROVED** (logical) + **NUMERICAL** | §4.2 |
| Z1 | `zeta`: a real branch stays real — `\|Im z\| ~ 1e-97` and `\|Im z'\| ~ 1e-98` against `\|Re z'\| ~ 0.58`, from an initial guess perturbed by `1e-3`; hence `D == 0` and `D' == 0` identically | **PROVED** + **NUMERICAL** | §4.3 |
| Z2 | `zeta`: the first backward collisions of three adjacent low pairs, located by solving `H_tau = H'_tau = 0` | **NUMERICAL** (scale-free residuals; see `results/part4_zeta.txt`) | §4.3 |
| Z3 | `zeta`: the emergent conjugate pair is continued below the collision; `Im z ~ sqrt(2(tau_c-tau))`, so `D ~ (tau_c - tau)` and `D' ~ -1` | **NUMERICAL** | §4.3 |
| Z4 | the (Z3) pair is **not** extremal at `tau ~ -1.28`, so (F7) does **not** apply to it: the exact two-zero model puts the Lehmer pair's collision scale at `-7.1e-4` against `-1.04` for the closest pair among the first twenty zeros — a factor `1459.3753` — so zeros leave the axis in order of increasing gap and the low-height pair is one of the last | **PROVED** (the two-zero model is exact) + **NUMERICAL** (the gaps) | `src/zeta_lehmer.py` |
| Z5 | **in neither control is the sign-forcing hypothesis verifiable**: DH's `O_1` is extremal only within a finite census and `O_2..O_4` provably are not; `zeta`'s tracked pair provably is not. Extremality quantifies over all zeros, so no finite computation can establish it | **PROVED** (logical) | REPORT §4.3 |

## Part 5 — the alternative evolution

| # | claim | class | where |
|---|---|---|---|
| A1 | the `kappa`-family is **rejected**: `Lambda_{F_kappa}(s) = Lambda_{F_kappa}(1-s)` holds **only** at `kappa = kappa_DH` (the root-number defect is `8.1e-51` there and `1.9e-3` at `kappa_DH + 10^-3`), `zeta` occurs at no `kappa`, and no inference returns to `zeta` | **PROVED** (the root-number identity) + **NUMERICAL** | `src/alt_evolution.py` |
| A2 | no justified alternative evolution is available in this repository; the branch is closed without a second candidate | **PROVED** by exhaustion of the stated preference + the parent audits' prior findings | REPORT §7 |

## Part 6 — noncircularity

| # | claim | class | where |
|---|---|---|---|
| N1 | **nothing here distinguishes `zeta` from Davenport–Heilbronn.** Every hypothesis used — evenness, reality, double-exponential decay, bounded strip, simple zero, canonical product — DH satisfies, and the DH control instantiates every conclusion, with all four certified quartets moving *toward* the line at `tau = 0` | **PROVED** (the control is the witness) | REPORT §8 |
| N2 | the one hypothesis in the inventory that DH **provably fails** is `m < infinity`, finitely many non-real zeros | **PROVED** given Davenport–Heilbronn 1936 (infinitude), a **bibliographic** statement not verified here. **Upgraded from CONJECTURAL** by AUDIT §1.3: the corrected hypothesis needs only the 1936 infinitude theorem, not Levinson's count and not a summability guess | AUDIT §1.3, §3 |
| N6 | **eleven of the fourteen constraints in the inventory are satisfied by DH**, which has certified off-line zeros and a certified negative Weil form; exactly **two** are primitive separators (`K12` Euler product, `K13` finiteness), one with its route closed by the parent Phase H and one attached to a conclusion proved weaker than de Bruijn's | **PROVED** (the control is the witness) | AUDIT §2, §3 |
| N7 | the flow block `{K3,K4,K5,K6}` is **two** independent facts, not four: `K4 ==> K3`, `K6` is weaker than `K3`, and `K5` is not implied by `K4` | **PROVED** | AUDIT §2 |
| N8 | the one window-block combination that works (`Q_rep,delta <= Q_add,delta`, so the certificate covers the count-preserving convention) is **PRE-EXISTING** — the parent's `cold_certificate.json` carries `Q_add_minus_Q_replacement = 3.410319615341804e-25` and `2 F(14)^2` reproduces it to 16 digits | **PROVED**, and **not new** | AUDIT §2 |
| **C1** | **CERTIFIED: `H_tau` has a nonreal zero at `tau = 4/25 = 0.16` exactly, hence `Lambda_DH > 0.16`.** Arb 700 bits, 512 arcs; disc centre `171.2502 - 0.1823i` (exact dyadic), radius `1/16`; `\|Im centre\| - r = 0.1198 > 0`; `min\|H_tau\| >= 5.756765507e-31`; winding number **exactly 1** by the half-plane method (net advance 4, 0 non-adjacent steps); negative control on a disc shifted `0.6` returns winding **0**. Rigorous tails `T_20 <= 9.954e-47`, `A_10 <= 6.296e-47` | **CERTIFIED** (T2-C). Improves the published certified `Lambda_DH >= 6.8268e-29` by `2.3e27`; bracket `[0.16, 1.6025…]`, a factor of `10.02` | AUDIT §4b |
| C1a | the representation `H_tau = (1/2) sum (tau/4)^m Lambda^{(2m)}(s)/m!` that C1 rests on, cross-checked at `tau != 0` against the defining integral by oscillation-resolving quadrature sharing no code: **relative** deviations `5.9e-31` to `9.6e-45` at five points (the first is an absolute `3.3e-64` on a value of size `5.6e-34`), quadrature precision doubling `1.6e-63` to `7.8e-92`, also relative; a wrong-constant control `(tau/2)^m` disagrees by `18026` | **NUMERICAL**, with a discriminating control | AUDIT §4b |
| C1b | the two inputs to C1's tail bounds are gated, not quoted: `S(1) = 1.3236800759485 <= 2` justifies `\|Phi_F\| <= 8e^{3u}e^{-pi e^{4u}/5}`, and `(log f)'(3) <= -409024.8675 <= -1` justifies `int_3^inf f <= f(3)` | **CERTIFIED** (Arb) | AUDIT §4b |
| C1c | `W = net/4` **exactly** (not up to an error): abutting balls, `\|e_i\| < pi/2` per arc, telescoping with `H(2pi) = H(0)` | **PROVED** | AUDIT §4b |
| C2 | the certificate needs **no** trajectory estimate, **no** isolation hypothesis and nothing uniform over an interval — one parameter, one contour, one integer | **PROVED** (by construction) | AUDIT §4b |
| L1 | **proposed lemma, not proved, and NO LONGER NEEDED (superseded numerically by C1):** sign identity (equality, so an *upper* bound on `\|D'\|`) + the unconditional shell count `=> Lambda >= y_1^3/(3(A y_1 + B))`. Measured on DH `O_1`: `Lambda_DH >~ 0.01366592` against the published certified `6.8268e-29`, a factor `2.0018e+26`, narrowing the bracket from 28 orders to `2.07`. Robust to inflating `C_0` by `10^6` | **CONJECTURAL** (plausibility estimate; the Jensen shell count is standard but not derived here, isolation checked only at `tau = 0`, branch simplicity for `tau > 0` unchecked, `C_0` measured) | AUDIT §4 |
| N3 | (F10)'s corollary is nevertheless **strictly weaker** than de Bruijn's classical bound (F11), so N2 buys nothing | **PROVED** | §3.7 |
| N4 | the whole route is **forward-only**: every statement propagates reality toward increasing `tau` and says nothing at `tau = 0`, which is where the arithmetic and RH live | **PROVED** (the direction of F7/F9) | REPORT §8 |
| N5 | `D` decreasing is **not** a positivity mechanism: at `tau = 0` the DH form is certifiably negative (`-0.679071346441480590986...`, parent extension) **while** all four of its certified off-line quartets have `D'(0) < 0` | **CERTIFIED** on both sides (the negative witness is the parent's, reused) | REPORT §8 |

---

## Corrections made by `AUDIT.md`

| what | class before | class after |
|---|---|---|
| proximity cutoff constant `Delta/sqrt2` | PROVED (valid) | **superseded** by the sharp `Delta/sqrt3` |
| `N`-competitor form `Delta sqrt(N)/2` | PROVED | **REFUTED** — contradicts its own `N = 1` case; correct form `Delta sqrt(N/2)` |
| (S2) hypothesis `Escr < infinity` | PROVED conditional | **REFUTED as stated** — with `m = infinity` the identity is incompatible; correct hypothesis `m < infinity` |
| the per-zero cutoff `R < y_j` | not stated | **PROVED** (available and missed) |
| the cross-check column header shortened to "series vs integral" | ambiguous | **labelled "relative"**; the script's own transcript always said "relative deviation", the ambiguity was introduced in transcription only |
| the `N = 1` cutoff quoted without the `N > 1` relaxation | stated unqualified | **qualified**: `R < min(N y_j, Delta sqrt(N/2))`; the `N = 1` form is not necessary for several competitors |
| the §4 trajectory lemma as the route to a better `Lambda_DH` | proposed | **superseded** by the fixed-parameter certificate C1, which is certified and 12x stronger |
| `Delta_DH/sqrt3` printed as `1.0336096` | — | **corrected to `1.0336142`**; caught by `tests/verify_numbers.py`, not by reading — the second arithmetic slip that verifier has found |
| `Q_rep <= Q_add` as a new combination | believed new | **pre-existing** in the parent extension |
| DH fails `sum_rho D < infinity` | CONJECTURAL | **upgraded to PROVED** given the 1936 theorem, for the corrected hypothesis |

## Withdrawn or not claimed

* **Not claimed:** a new proof of de Bruijn's theorem (F8 is the differential statement only).
* **Not claimed:** any improvement on the de Bruijn–Newman upper bound for `Lambda`, for
  `zeta` or for DH. (F10) is weaker than de Bruijn by (F11).
* **Not claimed:** any novelty for the zero-velocity equation `z' = H''/H'`, the pair
  dissipation of (F9), or the Hadamard/genus argument. These are textbook; the parent audit
  already records `z' = H''/H'` as such. What is assembled here is their application to the
  mean-gap diagnostic and the counterexample C1.
* **Not claimed:** that `Escr` is an "entropy". It is a sum of squared displacements with a
  proved monotonicity and an explicit dissipation; no thermodynamic reading is attached and
  none is needed.
* **Not claimed:** anything about where any zero of `zeta` is.
* **Explicitly reported as a failed first attempt:** the dictionary check `d_tau H = -d_z^2 H`
  at `z = 171.4` **failed** at dps 40 with a relative residual of `6.0e-10`. It was a
  quadrature-resolution artefact, not a failed identity; recorded in `DERIVATION.md` §0.3
  rather than silently fixed.
* **Caught by `tests/verify_numbers.py`, not by reading:** this file quoted `O_4`'s
  certified `D'(0)` as `-1.12778399346959232978` at 20 decimals. The Arb midpoint is
  `-1.1277839934695923297705678...`, whose 21st significant digit is `0`, so it rounds
  **down** to `...232977`. Corrected. The 30-digit value quoted in `REPORT.md` and
  `DERIVATION.md` was right; only this one truncation was wrong. That is exactly the class
  of error the parent extensions' note verifiers were built to catch, and it is why one
  exists here.
* **Explicitly reported as a wrong estimate of my own:** a mean-spacing heuristic predicted
  the DH cross-orbit term as `-1.3026848`; the true value is `-0.16043167`. A first, sloppier
  version of this very sentence quoted `-0.58`, which was an uncomputed intermediate; the
  heuristic is now evaluated in `src/dh_interaction.py` rather than estimated by hand. The heuristic was wrong
  because the quartet sits in a `3.05x` zero gap (H5). It is recorded because the same
  heuristic is what a dissipation lower bound would need.
