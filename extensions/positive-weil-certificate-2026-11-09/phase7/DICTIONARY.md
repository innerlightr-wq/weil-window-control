# Normalization dictionary and theorem-to-gap table

**Paper under audit.** Ronnie Andrews, Alain Connes, Antonio La Gatta,
*The uncertainty operator, the Weil ground state and the zeros of zeta: a
numerical study and conditional positivity of the Weil quadratic form*,
arXiv:2611.13129**v1**, 9 Nov 2026 (19:13:52 UTC), 57 pp.

* local pinned copy `notes/2611.13129v1.pdf`
* `sha256 = 6a48849465a22705540d9a5ad246156babcb73bdb8d0f0571fb07a68c943fe67`
* arXiv abstract page checked 2026-11-12: **v1 is the only version**. No later
  version to diff.

Repository side: `weil-window-control`, branch
`positive-weil-certificate-2026-11-09`, commit
`7515396c7f2aa17b53d565a4f11f8e937ee164dc`, subtree
`extensions/positive-weil-certificate-2026-11-09/`.

---

## 1. The dictionary — PROVED, then verified numerically

### 1.1 Variables, window, measure

| paper | repository | status |
|---|---|---|
| `u in R*_+`, `d*u = du/u` | `x in R`, `dx` | `x = t = log u` |
| window `supp g subset [lambda^-1, lambda]` | `supp f subset [-L, L]` | **`lambda = e^L`**, i.e. `L = log lambda` |
| `L^2(R*_+, d*u)` norm | `L^2(dx)` norm | identical: `d*u = dt` |
| `iota: u -> 1/u`, sectors `QW^ev / QW^odd` | parity in `x`, even/odd sectors | identical: `iota` is `t -> -t` |
| `F(x) = int g(u) h(ux) d*u` (eq. 6) | `C_{jk}(s) = int f_j(t) f_k(t-s) dt` | `F(e^s) = C(s)`, the additive autocorrelation |
| `lambda^2 = 5` (Table 1 first row) | `L = (log 5)/2 = 0.8047189562...` | — |
| `lambda^2 = e^{8/5} = 4.95303242439511480365` | **`L = 4/5`** (this repository's certified window) | — |

`lambda = exp(L)` is therefore **confirmed, not assumed**. Numerically
`L = 4/5` gives `lambda^2 = 4.953032424395114804`, which sits **0.94 % below**
the paper's smallest tabulated window `lambda^2 = 5`.

### 1.2 The form — term by term

Writing `R(t) = f(e^t) = f(t)` and `C(s) = int R(t)R(t+s) dt`
(so `C(0) = ||f||^2`, `C` even, `supp C subset [-2L, 2L]`):

| paper term | value in the additive variable | repository term |
|---|---|---|
| `W_{0,2}(F) = int F(x)(x^{1/2}+x^{-1/2}) d*x` | `2 (int R e^{u/2})(int R e^{-u/2})` | `Pole = 2 s F(i/2)^2`, `s = +1` even / `-1` odd |
| `W_p(F) = log p sum_m p^{-m/2}[F(p^m)+F(p^-m)]` | `sum_{n=p^m} (2 Lambda(n)/sqrt n) C(log n)` | `Prime_{jk}` — **identical, and finite/exact** (`log n < 2L`) |
| `-W_R(F)`, eq. (8) | Fourier multiplier `m(tau) - c_0` | `Arch - LogPi`, multiplier `Re psi(1/4+i tau/2) - log pi` |

**The archimedean identity, proved.** With
`m(tau) = int_0^inf (1 - cos(tau s)) e^{s/2}/sinh s ds` and
`c_0 = log 4pi + gamma_E + int_0^inf (e^{s/2}-1)/sinh s ds` (the paper's own
Lemma 8.13 decomposition),

```
    m(tau) - c_0 = int_0^inf [1 - e^{s/2} cos(tau s)]/sinh s ds - log 4pi - gamma_E .
```

The integral representation `psi(z) = -gamma + int_0^inf (e^{-u} - e^{-zu})/(1-e^{-u}) du`
at `z = 1/4 + i tau/2`, substituting `u = 2s` and using `1 - e^{-2s} = 2 e^{-s} sinh s`, gives

```
    Re psi(1/4 + i tau/2) = -gamma + int_0^inf [e^{-s} - e^{s/2} cos(tau s)]/sinh s ds .
```

Subtracting, the difference is `int_0^inf (1 - e^{-s})/sinh s ds - log 4`, and
`(1-e^{-s})/sinh s = 2 e^{-s}/(1+e^{-s})` integrates to `2 log 2 = log 4` exactly.
Hence

```
    m(tau) - c_0  =  Re psi(1/4 + i tau/2) - log pi      for every real tau.
```

**Conclusion (dictionary theorem).** For `f` supported in `[-L,L]` in the form
domain, with `g(u) = f(log u)` and `lambda = e^L`,

```
    QW_lambda(g, g)  =  Q(f)      and      ||g||^2_{L^2(d*u)} = ||f||^2_{L^2(dx)} ,
```

with no factor of 2, no sign convention left free, and the two pole signs
matching sector by sector (for odd `R`, `int R e^{u/2} = - int R e^{-u/2}`, so
`W_{0,2} = -2 (int R e^{-u/2})^2`, which is the repository's `s = -1`).
Consequently `epsilon_lambda(1) = lambda_Q^{even}(L)` and
`epsilon^{odd}_lambda(1) = lambda_Q^{odd}(L)` with `lambda = e^L`.

### 1.3 Numerical verification of the dictionary

`work/dictionary_check.py` implements eqs. (6), (7), (8) **directly from the
paper** — multiplicative variable, `log 4pi + gamma_E`, the `1/(x - x^{-1})`
kernel — with no code shared with the repository, and evaluates them on the
repository's own Legendre test functions at `L = 4/5`. Against
`phase5/certs/conventions.json`:

| test function | paper `Psi(F)` (independent) | repo `Q(f)` (certified ball) | rel. diff |
|---|---|---|---|
| `phi_0` | `0.080274054959269611213691659381317` | `[0.0802740549592696112136916593813 +/- 1.70e-32]` | `2.1e-31` |
| `phi_2` | `0.389336213205013655883159773665532` | `[0.389336213205013655883159773666 +/- 4.69e-31]` | `1.2e-30` |
| `phi_0 - phi_2` | `0.117009437226589934779775252494441` | `[0.117009437226589934779775252494 +/- 4.42e-31]` | `3.8e-30` |

Agreement to ~30 digits, i.e. to the width of the serialized repository balls.
The pole moments `int f e^{-u/2}` also match `conventions.json`'s
`pole_moment` field to 30 digits. **The dictionary is confirmed.**

### 1.4 What the dictionary does NOT identify

* The paper's **multiplicative window parameter** `lambda` is the window, not a
  cutoff. The repository's `T` (frequency cutoff, 200/400/800) and its
  `lambda_Q = min Rayleigh quotient` are different objects: `T` indexes the
  *comparison form* `R_T <= Q` inside a fixed window, and has no counterpart in
  the paper. The paper's truncation parameter `N` (number of log-Fourier modes)
  is the analogue of the repository's Legendre degree `N`, not of `T`.
* `gamma = 2 pi lambda^2` (prolate bandwidth) is unrelated to `gamma_E`.
* The paper's `epsilon_lambda(j)` are eigenvalues; the repository's `a`, `d`,
  `kappa`, `mu` are Schur-block constants. Only `mu` and `U` bracket an
  eigenvalue.
* The paper's `delta_n(lambda) = 1 - |chi_n(lambda)|` are **uncertainty-operator**
  rungs (Theorem 2.1), not Weil eigenvalues. They are what Hypothesis 3.4
  compares the Weil spectrum against.

### 1.5 One documentation defect found in the repository (not a result defect)

`DERIVATION.md` §1 prints the archimedean term as

```
    Arch_{jk} = -gamma C(0) + sum_{m>=0}[ C(0)/(m+1) - 2 int_0^{2L} C e^{-(m+1/2)u} du ]
```

As printed this series **diverges**: the `m`-th term behaves like
`C(0)/m - 2C(0)/m = -C(0)/m`. The correct grouping, which follows from eq. (8)
via `1/sinh s = 2 sum_{m>=0} e^{-(2m+1)s}`, is

```
    -W_R = -(log 4pi + gamma_E) C(0)
           + 2 sum_{m>=0}[ C(0)/(2m+1) - int_0^{2L} C(u) e^{-(2m+1/2)u} du ] ,
```

whose terms are `O(1/m^2)`. The constant differs from the printed one by exactly
`log 4` — the same `log 4` identified in §1.2 above. **Nothing computed depends
on the printed formula**: the certified pipeline uses the frequency form
`H(t) = Re psi(1/4+it/2) - log pi` (`phase2/src/baseline.py::H`,
`phase2/src/block.py`), and the closed-form time-domain evaluator uses
`zeta_window.arch_closed`. This is a transcription error in prose only, and it
should be corrected in `DERIVATION.md` before the note is deposited.

---

## 2. Theorem-to-gap table

Page/equation references are to v1. "Testable here" means: with the
repository's existing certified machinery, at a modest window.

| result | exact hypotheses | what it establishes | what it does NOT give | testable here |
|---|---|---|---|---|
| **Thm 2.1** (p. 9) `U_lambda = 2 - P_lambda - hat P_lambda` | none beyond Halmos' two-projection theorem | spectrum of `U_lambda` below 1 is exactly the simple `delta_n(lambda)`, `n` even; complement is the Sonin space with `U = 2` | nothing about `QW` | yes — and used: it is what makes `delta_n` *exact* rather than asymptotic |
| **eq. (17), §2.5** Fuchs ladder | fixed `n`, `gamma -> inf`, error `1+O(gamma^{-1})` | `delta_n ~ A_n gamma^{n+1/2} e^{-2gamma}` | **not uniform in `n`**; the paper itself warns the comparison "must use the exact `D_j`" (§4.2 item 1): at `lambda^2 = 5` Fuchs is off by `1.0` in `D_2` | yes — Appendix A reproduced, see §3 below |
| **Thm 3.1/3.2** (pp. 13–14) upper bound `epsilon_lambda <= C delta_4 (log lambda)^3` | Gevrey repair of [11,16]; `sigma_{4K} <= sigma_0 < 1` | an **upper** bound by exhibiting a test vector | no lower bound of any kind | not needed — the repository already has a *better*, certified upper bound at `L = 4/5` |
| **Hyp. 3.4**, eq. (32) (p. 14) | — | `c delta_{4j} <= epsilon_lambda(j) <= C delta_{4j} log lambda` | **the lower halves are the hypothesis.** The paper is explicit (Fact 4.1): because `QW^N_lambda` is a *compression*, min–max gives `epsilon^N_lambda(1) >= epsilon_lambda(1)`, so **every computed ratio is an upper bound** and no computation in the paper can test a lower half | **YES — and this is where this repository is the only available instrument.** See §3 |
| **§4.2–4.4, Tables 1–4** | `N`-truncation, Arb inverse iteration | `epsilon_lambda(j)/delta_{4j}` upper bounds, `1.40 ... 11.2` for `lambda^2 <= 50`; interlacing; the `~lambda^2` count | nothing certified from below; radii "meaningless" at `lambda^2 = 100` (their §4.1), midpoint iteration thereafter | comparison, yes |
| **Fact 5.3** (p. 23) prolate conditioning | numerical, `N`-truncated | after removing `lambda^2-1` (even) / `lambda^2-2` (odd) prolate directions and the pole form, `0.27 <= mu_lambda <= 0.78` for `lambda^2 <= 50`; `mu_lambda ~ 2/lambda` | not certified; not a lower bound on `QW_lambda` | in principle, but needs the `E`-map / Gevrey construction — **not** in this repository |
| **Lemma 8.13** (p. 47) | `f` real, in the form domain | `Lambda(f) <= QW_lambda(f,f) + C_0 lambda ||f||^2` | `C_0` not made explicit | its *proof* is what supplies the multiplier identity used in §1.2 |
| **Thm 8.14** (p. 47) | **(H2') at scale `sigma`**, `sigma > 1`, plus the choice of `n_0` | `QW_lambda >= (log sigma)/2 ||f||^2` on `W_lambda(sigma,n_0)`, of codimension `<= 2 + 2(n_0+1)|M|` | **conditional on (H2'), which is semilocal and open.** Remark 8.15(b): in the *archimedean* case all hypotheses are known — but that case is `S = {infinity}`, no primes, and settles nothing semilocal | no — (H2') is not a finite computation |
| **Lemma 8.16** (p. 48) Schur complement | `Q` bounded below; `dim V < inf`; `B*v` in the dual of the form domain of `D`; **`D >= d > 0` on all of `W`** | `epsilon_lambda >= min(0, lambda_min(S))`, `S = A - ||D^{-1/2}B*||^2`; and `Q >= 0 iff S >= 0` | nothing without `D >= d > 0`, which is exactly what Thm 8.14 (hence (H2')) supplies | **structurally identical to the repository's own chain** — see §4 |
| **Prop. 8.17** (p. 49) monotonicity | — | `epsilon_lambda` non-increasing in `lambda`; and `epsilon_lambda >= -eta(lambda)` **for all large `lambda`** with `eta -> 0` implies `QW_lambda >= 0` for every `lambda`, i.e. RH | the hypothesis is an *asymptotic family* statement. **No finite set of windows supplies it**, however many and however rigorous | the monotonicity half is usable immediately (see §4); the implication is not |
| **(R1)** (p. 49) | — | `Rad(QW)` is a convolution ideal — proved there in two lines | — | n/a |
| **(R2)** (p. 49) | uniform-in-index extension of [11] up to `tau`, plus a Gram lower bound beyond the Hermite regime | would give `alpha, beta = O(tau^{1/2} lambda^{C'})` in (73) | **"we leave this open"** — the paper's own words | no |
| **Prop. 8.18** (p. 50) | (H2'), (R1), (R2), `tau = lambda^{-A}`, **and** `lambda_min(S|_{V'}) >= -eta(lambda)`, `eta -> 0`, for every large `lambda` | RH | three open hypotheses plus an asymptotic family statement about the transition-and-pole block. Fuchs' fixed-`n` asymptotics is explicitly *not* used there, but uniformity at transition indices growing with `lambda` is still required | no |
| **§8.8 item 3** | — | asks for exactly "a numerical audit of the three bounds (73) at a modest cutoff, with the Gram conditioning, the pole block and the mixed terms separated and with rigorous error bounds" | — | **this is the paper's own request, and the closest match to this repository's machinery** |

### Classification of the paper's contribution to this project

| class | content |
|---|---|
| **usable theorem** | Thm 2.1 (exact rungs) and Appendix A (Bouwkamp): these give `delta_n(lambda)` exactly, independently of the Weil form. Prop. 8.17's **monotonicity half** (`epsilon_lambda` non-increasing) is unconditional and immediately usable. |
| **conditional reduction** | Thm 8.14 + Lemma 8.16 + Prop. 8.18: a reduction of RH to the sign of a Schur complement on an `O(A sigma^2 log lambda) + 2`-dimensional block — conditional on (H2'), (R1), (R2). Not usable as an input here; valuable as a *target description*. |
| **numerical inspiration** | Fact 5.3 (prolate conditioning leaves `mu_lambda ~ 2/lambda`), Tables 2/4 (the ratios), §4.4 (the counts). These told us where the smallness lives. |
| **no transfer** | §6, §7 (ground state, Mellin zeros, the three distances): a different target (zero locations), with no path to a window positivity constant. Thm 3.1/3.2's upper bound is weaker than what this repository already certifies at `L = 4/5`. |

**Critically: the paper provides no lower bound on `epsilon_lambda` at any
`lambda`, and by its own min–max remark its numerics cannot provide one.**
