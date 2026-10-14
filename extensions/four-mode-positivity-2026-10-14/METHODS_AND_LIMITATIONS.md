# METHODS_AND_LIMITATIONS.md

Every statement in this extension, sorted by what kind of thing it is. Read §1
as the authoritative scope of the whole directory.

---

## 1. Evidence tiering

### (a) Standard facts, no novelty claimed

| | |
|---|---|
| Weil's explicit formula in `psi`-form; Weil positivity ⟺ RH | classical |
| Positivity ⟺ `Q_zeta(f) >= 0` for **all** admissible `f` | classical |
| Rank-one threshold: for `H > 0`, `H - 4bb^T` has a negative direction iff `4b^T H^{-1} b > 1` | classical — Schur complement / Sherman–Morrison / Birman–Schwinger |
| Sylvester's criterion; Sylvester's law of inertia; exact rational `LDL^T` | classical |
| Gershgorin / row-sum bound `\|x^T Delta x\| <= (max_j sum_l rad_{jl})\|\|x\|\|^2` | classical |
| Polarization of a quadratic form; Hermite–Hadamard; integration by parts; integral/midpoint comparison of monotone and convex series | classical |
| AM–GM (`2uv <= tau u^2 + v^2/tau`) | classical |
| Outward-rounded exact rational interval arithmetic | standard technique |
| `sinh(t)/t` increasing | elementary |

### (b) Previously published in this repository, imported unchanged

The entire certified path except the tail bound: the explicit formula
implementation, the interval class and every proved transcendental series, the
canonical autocorrelation reduction `C(u) = sum_k[A_k sin(w_k u) + B_k(U-u)cos(w_k u)]`,
the pole / archimedean / `log pi` / prime-power terms, the prime support cutoff
(`log 5 > 2L`, so exactly `n = 2,3,4`), the `gamma_E` elimination
`H_M - gamma_E = ln(M+1) - T`, and the exact quartet closed forms. Imported by
`scripts/weil_core.py` from
[`extensions/planted-witness-closure-2026-10-08/verify_certificate.py`](../planted-witness-closure-2026-10-08/verify_certificate.py),
sha256 `c49858a1…`, **checked at import time** — a mismatch raises rather than
silently running on something else. Also published: the certified interval
`[1/4, 49/100]` and the point bound at `delta = 2/5`, both of which this
extension *recomputes* rather than cites (assertion groups 3–4 of
`verify_new_certificate.py`).

### (c) New analytic derivations — proved here

| | where | status |
|---|---|---|
| Proposition T: `O(1/M^2)` archimedean tail, two integrations by parts | [`ANALYTIC_TAIL.md`](ANALYTIC_TAIL.md) §3 | proved, with rigorous remainder |
| Lemma E: endpoint identities `C(U)=0`, `C''(U)=0`, `C'(U)=-C'(0)`, `C'(0)=0`, `C''(0)=-\|\|f'\|\|^2=-sum_k w_k^2 x_k^2` | §2 | proved; each also evaluated numerically and charged, not substituted |
| Proposition T3: `O(1/M^3)` archimedean tail, three integrations by parts, collapsing to an exact diagonal quadratic correction | §4 | proved, with rigorous remainder |
| Lemma S1: `O(1/M^2)` midpoint enclosure of `S_1` by convexity | §4 | proved |
| Lemma T: `O(1/M^3)` enclosure of the `gamma_E` series `T` | §5 | proved |
| Cell loss bound via AM–GM, and the derivative bounds `Ahat, Bhat` | [`NON_DETECTION.md`](NON_DETECTION.md) §2 | proved, elementary |

All dependencies of these are present in this directory or in the published
files it imports; nothing is quoted from an unavailable source.

### (d) Newly certified computational theorems

| | statement | verifier |
|---|---|---|
| **Theorem A** | `Q_add,delta(f_x) < 0` for every `delta in [9/100, 1/4]`, one frozen rational witness, 69 closed cells; union with the inherited certificate gives `[9/100, 49/100]` | `verify_new_certificate.py` |
| **Theorem O1** | `G >= m0 I` on the four-mode subspace, `m0 >= 1.659599e-10 > 0`, **unconditional** | `verify_route_o.py` `O1`,`O2`,`O2b` |
| **Theorem O2** | `M(delta) >= 0` — no negative witness in this subspace — for every `delta in [0, 3831/64000]`, 6072 closed cells | `verify_route_o.py` `O3`,`O4` |

Computer-assisted with ordinary runtime dependencies (CPython, its
`int`/`Fraction` implementation, the OS). **Not** proof-assistant formal
verification; nothing is formalised in Lean or any other system; **not**
externally peer-reviewed.

### (e) Numerical diagnostics — not certified

* `Theta(delta)` and the inertia of `M(delta)` in
  [`NON_DETECTION.md`](NON_DETECTION.md) §3, computed at the **midpoint** of the
  interval matrix. The onset bracket `(27/400, 137/2000) = (0.0675, 0.0685)` is a
  crossing bracket of that midpoint computation, **not a certified sharp
  threshold**.
* The spectrum of `Gt_mid` (`1.746e-10, 2.711e-06, 3.846e-03, 1.528`) and the
  condition number `8.75e+09`: these *are* exact rational inertia brackets of
  `Gt_mid`, but `Gt_mid` is the midpoint of the enclosure of `G`, so only the
  smallest, via Theorem O1, carries a certified consequence for `G` itself.

### (f) Withdrawn

* The onset bracket `(17/200, 9/100) = (0.085, 0.09)` and the whole `Theta`
  table of an earlier round of this work: they were computed from `G_core`, the
  **tail-free** matrix, not from `G`. `lambda_min(G_core) = 5.722182e-06` versus
  `lambda_min(G) = 1.746061e-10`. Both the numbers and the bracket are
  withdrawn; the corrected diagnostic is in
  [`NON_DETECTION.md`](NON_DETECTION.md) §3 and the diagnosis is in
  [`FOUR_MODE_POSITIVITY.md`](FOUR_MODE_POSITIVITY.md) §4.1. **The detection
  certificate never depended on it** (it is a direct sign evaluation with the
  tail included, no matrix), and the error was in the safe direction.

### (g) Open, recorded as ongoing work rather than as results

* Closing the gap `(0.05986, 0.09)` between certified non-detection and
  certified detection. Theorem O2's run was budget-limited; the target was
  `0.065` and it reached `0.05986`. A matter of compute at the present cell
  widths, not of a missing idea — but it was not done, so nothing about it is
  claimed, and the numerical onset bracket sits **inside** that uncertified gap.
* A certified (not midpoint) statement about the onset itself, which would need
  monotonicity of `Theta` or a rigorous enclosure of the crossing.
* Porting the sharper tail into the independent implementation, so that the left
  endpoint `delta = 9/100` is corroborated by two implementations rather than
  one.
* Any statement uniform in `L`, in `gamma`, or in the dimension.

## 2. Duplication check against the author's existing work

Done **before** the mathematics, and it changed the deliverable. The repository
and the local checkout were searched for `inverse-Gram`, `rank-one`,
`Sherman–Morrison`, `Schur complement`, `detection threshold`,
`near-null overlap`, `stiffness`, `node constraint`, `baseline-weighted`,
`H^{-1}`.

**A closely related earlier attempt exists in the author's local notes** (an
unpublished working directory, not part of this repository and not reproduced
here). Its scope, recovered rather than assumed:

| | that attempt | this extension |
|---|---|---|
| dimension | `N = 12, 16, 20, 24` modes | `N = 4` modes |
| displacement regime | infinitesimal, `delta_* ~ 4.3e-7 … 5.9e-7` | finite, `delta ~ 0.06 … 0.49` |
| perturbation | **rank one**, `H_t = H_0 - 4t vv*`, `t = delta^2` (the Taylor/response model) | **exact finite quartet**, `M = G + 4a(delta)a(delta)^T - 4b(delta)b(delta)^T`, both vectors `delta`-dependent |
| baseline | `H_0 = Q_0 = M_N + 2\|c><c\|` (the **replacement** convention) | `H = G + 4aa^T` (the **additive** convention) |
| weighting tested | `S_v = v^* H_0^{-1} v` against the unweighted near-null overlap | `Theta = 4 b^T H^{-1} b` against a fixed vector's `4(b.x)^2/(x^T H x)` |
| positivity input | **external**: Zhu's `lambda_min(M_N) >= lambda*(0.8) >= 8.9e-18`, quoted with its normalisation and range | **none** — Theorem O1 is unconditional |
| verdict | **null/mismatch**: the two agree to `1.12%`, because `lambda_2/lambda_1 = 5.4e5` makes the near-null space effectively one-dimensional, so "the weighting has nothing to reweight" | — |

**Consequences, which are the reason this extension looks the way it does.**

1. **The baseline-weighted matrix criterion is not new.** It already existed in
   the author's own notes and is classical anyway (Proposition 1). It is
   **reused, not renamed**, and no novelty is claimed for it. The new content
   therefore had to be in the *certificates*, not in the criterion — which is
   what §1(c), §1(d) deliver.
2. The earlier null verdict was reached in the **near-singular,
   high-dimensional, infinitesimal-`delta`** regime for the **rank-one** model
   against the **replacement** baseline. Its own stated reason
   (`lambda_2/lambda_1 ≈ 5.4e5` at `N >= 12`) is a large-`N` spectral property.
   At `N = 4` the measured ratio is `2.711e-06/1.746e-10 ≈ 1.55e+04`, so the
   same explanation applies in weaker form — and indeed the three-way comparison
   at `N = 4` found no advantage for the weighting either. That is a **negative
   result, reported as such**, not a result of this extension.
3. The earlier notes also record, and this extension inherits: that the rank-one
   threshold is classical with no novelty claimed; that `arXiv:2010.15102` was
   **not** consulted and is not relied on; and the separation of the "near-null"
   and "signal-dominant reoptimised" regimes of
   [`notes/proofs.md`](../../notes/proofs.md) §B.5, with the finding that a
   `lambda = 0` threshold cannot reproduce the reoptimised `C(L)` (shortfall
   `2.51e3`).

## 3. Reproduction environment

Python 3.14.7, CPython standard library **only** — no `numpy`, no `mpmath`, no
`sympy`, no `scipy`, and none installed on the machine where this was produced.
Nothing in this directory imports a third party package. Measured on one core:

| | runtime | peak RSS |
|---|---|---|
| `verify_new_certificate.py` (69 cells) | ~17 min | 16 MiB |
| `verify_route_o.py` (`O1`–`O5`, 6072 cells) | ~45 min | ~20 MiB |
| `scripts/route_o1.py` (Theorem O1) | ~5 min | 17 MiB |
| `scripts/indep_crosscheck.py` | ~8 min | 16 MiB |
| `tests/*.py` | < 2 min each | < 20 MiB |

The first run of anything recomputes the `x`-independent archimedean sums
(~4 min) and caches them in `results/arch_sums_cache.json`; the cache key
carries `M`, the grid bit depth, `L`, `gamma`, the dimension, the basis string
and the **sha256 of the published verifier**, so it cannot go stale across a
change of input. Delete it for a cold run. See [`REPRODUCE.md`](REPRODUCE.md).

## 4. Limitations, collected

1. **One window, one height, one basis, one dimension.** `L = 4/5`,
   `gamma = 14`, four `EVEN_D` modes. Nothing is uniform in any of these.
2. **The quartet is hypothetical.** Four zeros at `1/2 ± delta ± 14i` are
   **added** to the true zero multiset and none is removed. `gamma = 14` is a
   prescribed rational height, **not an observed ordinate of `zeta`**, and no
   zero of `zeta` is located, constrained or verified anywhere here.
3. **Four-mode positivity is not Weil positivity.** Theorem O1 is a
   finite-dimensional, measure-zero slice of the condition equivalent to RH, and
   its sign is the one RH predicts. **Nothing here proves, weakens, or
   conditionally establishes the Riemann hypothesis.**
4. **Non-detection in a finite basis is not positivity.** Theorem O2 says no
   witness exists *in this four-dimensional space*. A larger basis could detect
   at smaller `delta`.
5. **The onset is not certified** and the gap `(0.05986, 0.09)` is open.
6. **The left endpoint `9/100` rests on one implementation.** The independent
   route's budget is ~`90x` too loose there.
7. **Computer-assisted, not formally verified, not peer-reviewed.** CPython's
   integer and `Fraction` arithmetic and the OS are trusted. No proof assistant
   was used.
8. **`Q_zeta > 0` on one vector is not positive definiteness** — that needed
   Theorem O1, and Theorem O1 is still only about this subspace.
