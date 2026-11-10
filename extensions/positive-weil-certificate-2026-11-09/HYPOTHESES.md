# HYPOTHESES.md — input ledger and control table

Build: `extensions/positive-weil-certificate-2026-11-09`, dated 2026-11-09.
Repository HEAD at branch point: `ad78c59a6ca7d2a376c86220e583c5ae94e4962b`.

Status vocabulary, used strictly:

| status | meaning |
|---|---|
| **PROVED** | written argument here, with its assumptions named |
| **IMPORTED** | theorem taken from a named source, located by number |
| **CERTIFIED** | valid Arb-ball / exact-rational certificate *including* analytic errors |
| **NUMERICAL** | precision-stable but not certified |
| **UNRESOLVED** | no valid conclusion |

---

## 1. Positivity core (coefficient-sign condition only)

Nothing in this block uses continuation, the functional equation, gamma factors,
conductor, growth, or remainder estimates.

| # | hypothesis | source | applies to | used in | status |
|---|---|---|---|---|---|
| K1 | `log F(s) = sum_{n>=2} c_n n^{-s}` converges absolutely for `Re s > 1` | definition | each control object separately | the whole core | **PROVED** per object |
| K2 | `c_n >= 0` | per object | zeta; `zeta*L(chi)`; `zeta(2s-1)` | sign of `sum c_n n^{-sigma} P(t log n)` | **PROVED** (A3, C1, D1) |
| K3 | `P(theta) = 3 + 4cos theta + cos 2theta = 2(1+cos theta)^2 >= 0` | exact factorization in `Q[cos theta]` | all real `theta` | the 3-4-1 inequality | **PROVED** (A1, A2) |
| K4 | `sum_n c_n n^{-sigma} P(t log n) >= 0` for `sigma > 1` | K1+K2+K3 | same | `zeta(s)^3 |zeta(s+it)|^4 |zeta(s+2it)| >= 1` | **PROVED** |

## 2. Analytic wrapper (declared separately, as the brief requires)

| # | hypothesis | source | used in | status |
|---|---|---|---|---|
| W1 | meromorphic continuation of `zeta` past `Re s = 1` | classical | boundary theorem | **IMPORTED** |
| W2 | simple pole at `s = 1`, residue 1 | classical | the `(sigma-1)^{-3}` factor | **IMPORTED** |
| W3 | regularity of `zeta` at `1 + 2 i t_0` for `t_0 != 0` | W1 + W2 | the third factor stays bounded | **IMPORTED** |
| W4 | explicit-formula / Weil pairing in the window normalization of `src/zeta_window.py` | `notes/stage3_assembly.md` | every matrix entry | **IMPORTED** (repo) |
| W5 | `Re psi(1/4 + it/2)` series `-gamma + sum_n [1/(n+1) - (n+1/4)/((n+1/4)^2+t^2/4)]` | classical digamma series | gate monotonicity G1–G3 | **IMPORTED** |

## 3. The finite-to-infinite reduction

| # | hypothesis | source | used in | status |
|---|---|---|---|---|
| R1 | one-stroke reduction `Q(f) >= 2F(i/2)^2 + (1/pi)int_0^{T#}(Psi_L-beta)|F|^2 + beta||f||^2` whenever `Psi_L >= beta` on `[T#,inf)` | Zhu arXiv:2608.24827**v2**, Thm 1.1 | would convert a finite PSD check into whole-window positivity | **IMPORTED, NOT EXECUTED HERE** |
| R2 | `sup_t P_L(t) = A_L`, optimal by Weyl equidistribution | Zhu, Lemma 3.2 | the comb envelope | **IMPORTED** |
| R3 | `t -> Re psi(1/4+it/2)` strictly increasing on `t>0` | here, from W5 | sharpened gate `beta_sharp` | **PROVED** (G1–G3) |
| R4 | leading `N x N` Legendre block PSD + super-exponential tail `eps_D, eps_B` | Zhu, Thm 1.1 | the actual L=0.8 lower bound | **NOT REPRODUCED** — see SOURCE_AUDIT §3 |

**Consequence.** This build does **not** establish any whole-window lower bound.
`m_L` is **UNRESOLVED** at every `L` here. What is reproduced is the reduction's
*scalar gate* and its *cost*, not its matrix step.

## 4. Control table (Phase A)

| control | object | coefficient core | verdict | where |
|---|---|---|---|---|
| **A** | `zeta` | `c_{p^k} = 1/k`, else 0 | **PASSES**: `c_n >= 0`, `b_n = c_n log n = Lambda(n)` | A3, A4 |
| **B** | Davenport–Heilbronn, `a_n` periodic mod 5 `= (1,k,-k,-1,0)`, `k = [sqrt(10-2sqrt5)-2]/(sqrt5-1)` | `b_3 = -k log 3 < 0` | **FAILS**, as it must | A5–A9 |
| **C** | real nonprincipal `chi` | `chi(p)^k/k` can be `< 0`; `(1+chi(p)^k)/k >= 0` | unpaired **FAILS**, paired `zeta*L(chi)` **PASSES** | C1–C3 |
| **D** | `F_2(s) = zeta(2s-1)` | `c_{p^{2k}} = p^k/k >= 0` | **PASSES the core, FAILS the 1/2-line normalization** | D1–D4 |

Control **B** notes, as the brief requires: the coefficient hypothesis failing for DH
does **not** mean every positive-polynomial sum on DH must be negative. The regression
control for an actual negative Weil witness is the repository's certified
`Q_add,delta < 0` on `delta in [9/100, 49/100]`, re-run here (Phase B regression).

Control **D** is the decisive one for scope: it **passes** K1–K4 and still has its
completed symmetry centred at `3/4`. So the coefficient core is **not sufficient** for
a conclusion about the line `1/2`. It is *not* a counterexample to RH, and *not* a
counterexample to a correctly normalized L-function statement.

Control **C** note: the pairing passing K2 does **not** exclude an exceptional real
zero of `L(s,chi)`. Nothing here bears on that.

## 5. Statement F control

An ingredient claimed to separate `zeta` from DH must (i) hold for `zeta`, (ii) fail
for DH, and (iii) enter an estimate that controls zero geometry.

| ingredient | (i) zeta | (ii) fails for DH | (iii) controls zero geometry | passes F? |
|---|---|---|---|---|
| coefficient positivity `c_n >= 0` | yes (A3) | yes, `c_3 = -k < 0` (A8) | **no** — control **D** has it and sits at `3/4` | **NO** |
| comb mass `A_L` finiteness | yes | yes (DH has no Euler product, `A_L` is not even defined the same way) | only through `T#`, a *cost*, not a zero bound | **NO** |

Passing F would not be sufficient for RH or for universal Weil positivity in any case.

## 5a. Phase 2 ledger (baseline at L = 4/5)

| # | item | status |
|---|---|---|
| B1 | `Psi_L(t) >= beta := H(T)-A_L` for `|t| >= T` | **PROVED + CERTIFIED** |
| B2 | `Q(f) >= R_T(f)` (split, exterior discarded) | **PROVED** |
| B3 | `F_n(t) = i^n sqrt(2L(2n+1)) j_n(tL)`; `sum (2n+1) j_n^2 = 1` | **PROVED** |
| B4 | `g_N <= (2LT/pi) S_N(TL)` | **CERTIFIED** |
| B5 | Schur-type block bound | **PROVED** |
| B6 | `a = 1.03e-17` lower-bounds the least eigenvalue of the `N x N` Legendre compression of `R_T`, both sectors | **CERTIFIED** (Rump, `||I-G||` = `2.35e-25` / `4.59e-33`) |
| B7 | **whole-window `Q_zeta(f) >= 1e-17 ||f||^2` at `L = 4/5`** | **CERTIFIED**, margin `1.5554e-19` |
| B8 | parity decomposition covers arbitrary complex admissible `f` | **PROVED** (pole term `= 2C^2 - 2S^2` decouples; `Q(u+iv) = Q(u)+Q(v)`) |
| B21 | precision is OWNED: scope restored on return and exception, caches keyed by precision, wide/non-finite = UNRESOLVED | **PROVED + TESTED** (13 checks) |
| B22 | inputs built outside a precision scope carry the ambient radius and cannot be recovered | **DEMONSTRATED** (3 instances this run) |
| B23 | `mu_odd = 2.3e-17` certifies for the WHOLE odd sector | **CERTIFIED**, margin `1.2116e-18` |
| B24 | `inf_even <= U = 1.901382911145464406030000e-17 < 2.3e-17 <= inf_odd`; the even sector determines `lambda_Q` | **CERTIFIED** (separation `3.986170888545355939700000e-18`); endpoint corrected 2026-11-10, R13 |
| B25 | `R_200(f_400) = 1.24275215e-17 +/- 8.10e-26` by direct scalar quadrature | **CERTIFIED** |
| B26 | `beta_800 = 1.90476107093077689 > beta_400` | **CERTIFIED** |
| B27 | `e_800/e_400` = `0.498` (`f_400`), `0.48788` (`f_200`); but `e_400/e_200` = `0.711` on `f_400` | **CERTIFIED** — ratios are NOT constant, so no decay law is claimed |
| B28 | `R_800(f_400) = 1.6681e-17`, `R_800(f_200) = 3.01603e-17` | **CERTIFIED** |
| B29 | `lambda_R(800) <= 1.6681e-17`, so the whole `T=800` route is capped at `+18.3%` over `1.41e-17` | **CERTIFIED** |
| B15 | cutoff identity: `R_{T2,b2}-R_{T1,b1} >= 0` for `T2>=T1`, `b2>=b1` | **PROVED** |
| B16 | `beta_400 = 1.211613695058065 > beta_200` | **CERTIFIED** |
| B17 | `T=400`, `N=440` (degrees 0..439, 220/parity): `a=1.42e-17` even, `2.4e-17` odd | **CERTIFIED** |
| B18 | **whole-window `Q_zeta(f) >= 1.41e-17 ||f||^2`** at `L=4/5` | **CERTIFIED** |
| B19 | compact time-domain `Q`: `Arch(g)` + exact `c = -g'(0+) = f(L)^2`; inner GL exact | **PROVED + CERTIFIED**, validated on `phi_0` to 21 digits |
| B20 | `q(f_400)` enclosed, midpoint `1.90138291114546440e-17`, upper endpoint `1.901382911145464406030000e-17`; quote as `1.9014e-17` | **CERTIFIED** — best upper bound for `lambda_Q` here; see R13 |
| B10 | `R_{T,beta}` is monotone non-decreasing in `beta`; `A(b2)-A(b1) = (b2-b1)(I-K_T)`, `0 <= K_T <= I` | **PROVED + CERTIFIED** (residual `1.6e-353`) |
| B11 | `lambda_R <= lambda_Q`, `lambda_R <= lambda_A`; no relation between `lambda_A` and `lambda_Q` | **PROVED** |
| B12 | `r_f = 1.031022030245990e-17` on the frozen even-block minimizer | **CERTIFIED** (upper bound for `lambda_R`, `lambda_A`) |
| B13 | `e_f >= 1.530372869310224e-17` on that vector (band `[200,600]`) | **CERTIFIED** lower bound |
| B14 | `q_f` (variational upper bound for `lambda_Q` from that vector) | **NOT COMPUTED**; `q_f >= 2.561e-17` does not improve `2.112992612e-17` |
| B9 | quadrature remainder (Bernstein ellipse, `b_e = 0.30`, `rho = 5`, order 64) folded into the entries **once** | **CERTIFIED**, max entry radius `4.99e-60` |
| L-A | Lemma A, finite-measure energy control | **PROVED** (standard) |
| L-B | Lemma B, uniform local occupation | **PROVED** (standard) |
| P1 | an occupation bound consumed in a certificate | **NONE** — `E` empty at the gate; penalty unpayable at a smaller cutoff |

## 5b. Retracted

| # | statement | status |
|---|---|---|
| X1 | `meas{t in [T,2T] : P_L(t) > theta A_L} <= c T^{1-eta}` | **REFUTED** — positive density by continuous Kronecker–Weyl; see `RETRACTIONS.md` R1 |
| X3 | the printed midpoint `1.90138291114546440e-17` is an upper bound for `lambda_Q` | **FALSE, corrected** — it is the ball midpoint; `RETRACTIONS.md` R13 |
| X2 | "no bound on `meas(E)` can control `int_E |F|^2`" | **FALSE, withdrawn** — Lemmas A and B do exactly that; only the *global density* version is true. `RETRACTIONS.md` R6 |

## 6. Explicitly NOT assumed

- nonnegative log-coefficients imply RH — refuted by control D;
- positivity methods universally stop at `1 - c/log t`;
- Vinogradov–Korobov does not use positivity;
- absence of a witness proves positivity;
- a determinant/Schur reformulation is a new arithmetic estimate;
- an added quartet is an actual zeta zero configuration;
- the remaining-zero form is positive;
- the near-resonant set of the prime comb is sparse in the power-saving sense
  (**refuted**, X1);
- clearing the envelope gate implies the finite reduction is positive.
