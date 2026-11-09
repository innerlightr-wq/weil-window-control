# REPORT.md

Build dated 2026-11-09, branch `positive-weil-certificate-2026-11-09` off
`ad78c59`. All numbers below are reproduced by `verify.py`.

> **PHASE 2 (2026-11-09).** Baseline whole-window positivity at `L = 4/5`:
> **UNRESOLVED**, with the gap now isolated to a single quantity. The exterior
> gate, the function-space complement and the coupling are all **CERTIFIED** and
> are *not* the obstruction: at `N = 220` the chain closes as soon as
> `lam_A >= 1.0146e-28`, where `lam_A` is the least eigenvalue of the `N x N`
> Legendre compression. That matrix was **not assembled** here. See §7.
> A further claim is corrected: `RETRACTIONS.md` **R6**.
>
> **Revised 2026-11-09 after external review. See [`RETRACTIONS.md`](RETRACTIONS.md).**
> The "next step" of the first version — a power-saving sparsity bound on the
> near-resonant set — is **FALSE and fully retracted** (R1): by continuous
> Kronecker–Weyl that set has *positive density*, confirmed here to four
> significant figures. Three further statements are corrected: the exact
> multiplier applies to `T_1` only (R2), the 69 recomputed cells cover
> `[9/100, 1/4]` with `[1/4, 49/100]` imported (R3), and a `923 x 923` assembly
> is an *attempted* certificate, not a guaranteed one (R4).

---

## 1. The strongest conclusion, and its domain

> **Certified.** For `L = 81/100` and the even sector of `L^2(-L,L)`,
> ```
>     lambda_min(81/100)  <=  6.565025426e-18 ,
> ```
> with the enclosure `[6.56502542589e-18 +/- 4.92e-31]`, from a single frozen
> rational test vector (denominator `10^60`, 28 even cosine modes) evaluated in
> Arb ball arithmetic against the exact closed-form Weil kernel. At `L = 4/5`
> the same pipeline gives `lambda_min(4/5) <= 2.112992612e-17`.

**Domain.** This is an **upper** bound on the window infimum. It is
unconditional — one admissible test function suffices, no completeness, no
truncation of the function space, no appeal to zeros or to RH. It says **nothing
about positivity**.

> **Proved.** `t -> Re psi(1/4 + it/2)` is strictly increasing on `t > 0`
> (termwise from the digamma series), hence for every `T#`
> ```
>     Psi_L(t) >= beta_sharp(T#) := Re psi(1/4 + i T#/2) - log pi - A_L ,   t >= T# ,
> ```
> which dominates the imported `beta_paper(T#) = log(T#/2pi) - 1/T# - A_L`.

> **Certified, exact.** Crossing `L_5 = (log 5)/2 = 0.804718956…`, the comb mass
> jumps by **exactly** `2 log 5 / sqrt 5 = 1.43952503111…`, so the `T#` demanded
> by the one-stroke reduction multiplies by **exactly** `exp(2 log 5/sqrt 5) =
> 4.2186915938…`.

**No whole-window lower bound is established at any `L`.** `m_L` is
**UNRESOLVED**.

## 2. Controls

| control | verdict | why |
|---|---|---|
| **A** zeta | **passes** the coefficient core | `c_{p^k}=1/k >= 0`, `b_n = Lambda(n)` |
| **B** Davenport–Heilbronn | **fails**, as required | `b_3 = -kappa log 3 < 0` with `kappa = 0.284079…> 0`; re-derived here from the Dirichlet inverse, together with `b_6 = (1+kappa^2) log 6 = 1.936356…` |
| **C** real `chi` | unpaired **fails**, `zeta·L(chi)` **passes** | `chi(p)^k/k` can be `-1`; `(1+chi(p)^k)/k in {0,1/k,2/k} >= 0` |
| **D** `zeta(2s-1)` | **passes the core and still sits at 3/4** | `c_{p^{2k}} = p^k/k >= 0`; zeros map `rho -> (1+rho)/2` |
| **regression** repo's negative witness | **passes** | `verify_new_certificate.py`, exit 0. **Recomputed**: 69 cells covering `[9/100, 1/4]`, worst unit-`L2` margin `1.276748e-06`. **Imported**: the published certificate on `[1/4, 49/100]`, of which the join at `1/4` and one interior point `Q_add,2/5 <= -0.01451076575` were recomputed. The union is certified only by combining the two (R3). |

**Control D is the load-bearing one.** It satisfies every hypothesis of the
positivity core and its completed symmetry is centred at `3/4`. So
**coefficient positivity is not sufficient** for a conclusion about the line
`1/2`. It is not a counterexample to RH, nor to a correctly normalized
L-function statement. Control B failing does **not** make every positive-polynomial
sum on DH negative; that is why the certified negative witness is kept as a
separate regression.

**Statement F.** Neither candidate separator passes: coefficient positivity
holds for `zeta`, fails for DH — but control D has it too and sits at `3/4`, so
it does not enter any estimate controlling zero geometry. Comb-mass finiteness
controls a *cost*, not a zero location.

## 3. Does `L = 0.81` close? No. The deficit, measured.

It does not close, and the obstruction is **not** that some inequality was found
false. It is that the reduction's own gate moves discontinuously.

| quantity | `L = 4/5` | `L = 81/100` | factor |
|---|---|---|---|
| active prime powers (`log n < 2L`) | `{2,3,4}` | `{2,3,4,5}` | +1 |
| comb mass `A_L` | `2.94197352522` | `4.38149855633` | `+2log5/sqrt5` exactly |
| `T_1 = 2pi e^{A_L}` | `119.086556188` | `502.389453522` | `4.218691593787` — **exact**, and only here (R2) |
| exact gate `T_gate` | `119.086906076` | `502.389536459` | `4.218679895` (**not** exact, R2) |
| paper's root `T_paper` | `120.082403914` | `503.388460911` | `4.192025180` (**not** exact, R2) |
| `T#` at the paper's safety ratio `1.6655` | `200` | `~838` | `4.19` |
| modes `N ~ e L T#/2` | `~217` | `~923` | `4.24` |
| assembly cost `~ N^2 T#` | `1` | **`~75.5`** | — |
| best certified upper bound here | `2.1130e-17` | `6.5650e-18` | `0.311` |

So at `L = 81/100` the same route must do roughly **75x the assembly work** to
certify a margin that is **at most `6.565e-18`** — below the `8.9e-18` already
certified at `L = 0.8` by the external source. (The monotone decrease itself is
trivial: a larger window is a larger space, so `lambda_min` cannot increase. What
is new is the *certified value*: the target a `0.81` certificate must clear.)

**Which inequality fails.** None is violated. What fails is *feasibility*: the
gate needs `T# > T_gate(81/100) = 502.3895…`, and the Legendre/Bessel block that
Theorem 1.1 then requires has `N ~ 923` at `>= 50` digits. That assembly was not
performed here, and the external reproducibility materials could not be obtained
(`SOURCE_AUDIT.md` §3), so even the `L = 0.8` matrix step is **unverified** in
this build.

**What quantitative improvement would suffice.** Either
(a) **attempt** the `~923 x 923` assembly with certified tails `eps_D, eps_B` —
roughly 75x the published compute. This is an *attempt*, not a guarantee (R4):
clearing the envelope gate is necessary but does not establish that the finite
reduction stays positive once the tail errors and the complement bound are
included. Or (b) lower the gate — **closed to pointwise methods**, since
`sup_t P_L = A_L` is attained at `t=0` and approached at arbitrarily large `t`
by Weyl equidistribution, so no pointwise envelope beats `A_L`.

## 4. Was any new arithmetic estimate obtained?

**No.** The one thing proved here that is not in the imported sources —
`beta_sharp >= beta_paper` from digamma monotonicity — moves a **constant**:
`+0.0050` at `T# = 200` (`0.83%` off the admissible `T#` at `L = 4/5`, `0.20%`
at `L = 81/100`). It does not touch `A_L`, and therefore does not touch the
doubly exponential threshold `T_1 = 2pi e^{A_L}` with `A_L = (4+o(1))e^L`. It is
a sharpening of bookkeeping, not an arithmetic estimate, and it is reported as
such.

The Schur/block route of the brief's §5 was **not** pursued as new content: a
prior audit in this repository (`channel-compatibility-exploration-2026-10-24`)
already established that such a reformulation "only reorganizes the existing
positivity problem", and that any new estimate must bound `E^T S^{-1} E` from
**analytic** data. Nothing here does that, so nothing is claimed.

## 5. Limitations

1. **Whole-window status UNRESOLVED at every `L`.** The finite matrices give
   upper bounds only; the infinite-dimensional complement (defect (ii)) is
   untouched. Note that defect (i), analytic truncation per entry, is *absent*
   in this normalization — the prime sum is finite by support and the
   archimedean term is closed-form — so (i) was never the obstacle, and bounding
   it would not have discharged (ii).
2. **The external `L=0.8` lower bound is imported, not verified.** Its matrix
   step was not re-run; its reproducibility materials were not located.
3. **`beta_sharp` is proved but unused.** No certificate in this build consumes
   it, because the matrix step it would feed was not performed.
4. **The second implementation is not an independent derivation.** `certify.py`
   re-implements the *same* closed form as `src/zeta_window.py` in a different
   arithmetic backend. Agreement (V3a, max difference 0 at `N=6`) tests the
   arithmetic and the coding, not the derivation.
5. **`lambda_min` here is the even sector.** The odd sector is assembled and
   checked for finiteness (V4) but is not part of the certified bounds.
6. No quantitative zero-free region; only the boundary theorem.

## 6. Single best next step — REVISED

**Complete the `L = 0.8` lower-bound reproduction before scaling anything.**

The pipeline built here has never produced a whole-window lower bound. Until it
has produced one — by reproducing the external `L = 0.8` result end to end,
including the Legendre/Bessel assembly and the tail constants `eps_D, eps_B`
that `SOURCE_AUDIT.md` §3 records as unobtained — scaling it 75x to `L = 0.81`
would demonstrate nothing. That reproduction is also the only way to test
whether this pipeline *can* certify positivity at all.

**The previous "next step" is retracted.** It proposed a power-saving bound
`meas E_{L,theta}(T) <= c T^{1-eta}` on the near-resonant set. That is **false**
at every fixed `L`: the set has *positive density*, equal to the Haar measure of
an open neighbourhood of the origin in the prime-phase torus, by continuous
Kronecker–Weyl. Measured here, the density converges to that Haar measure to
four significant figures over four decades of `T` (`RETRACTIONS.md` R1).

**The research question that survives**, correctly posed: density alone cannot
control `int_E |F|^2`, because `g(u) = f(u) cos(omega u)` keeps the support and
the parity of `f` while translating its Fourier profile onto any chosen
resonance. So the question is not how small `E` is, but

> how much Fourier energy an admissible **window-supported** `f` can place on
> the unfavourable set `E`, and whether the favourable region absorbs it.

That needs (i) effective equidistribution with the correct **positive-density**
main term plus an error, (ii) a concentration bound for window-supported `f` on
`E`, and (iii) an explicit inequality combining them. Nothing here establishes
any of the three, and the failure of the retracted estimate does not cap other
extensions of the detector.


---

# 7. PHASE 2 — baseline at L = 4/5, and the local-occupation pilot

## 7.1 Status

> **Whole-window positivity at `L = 0.8`: UNRESOLVED**, because the `N x N`
> finite block was not assembled or certified. Everything else in the chain is
> certified, and the residual requirement is `lam_A >= 1.0146e-28` at `N = 220`.

## 7.2 What is now certified (the reduction, derived here)

Derived in the repository normalization, not imported (`phase2/src/baseline.py`,
`DERIVATION.md` §5). With `H(t) = Re psi(1/4+it/2) - log pi`, `beta := H(T) - A_L`:

| step | content | status |
|---|---|---|
| 1 gate | `H` strictly increasing (G1–G3) and `P_L <= A_L` give `Psi_L(t) >= beta` for all `|t| >= T`. At `T = 200`, `beta = 0.51846573324385`. **The exterior set is empty.** | **PROVED + CERTIFIED** |
| 2 split | `Q(f) >= R_T(f) = beta||f||^2 + 2 s F(i/2)^2 + (1/2pi)int_{-T}^{T}(Psi_L-beta)|F|^2` | **PROVED** |
| 3 basis | Legendre: `F_n(t) = i^n sqrt(2L(2n+1)) j_n(tL)`, with `sum_n (2n+1) j_n(x)^2 = 1` | **PROVED** |
| 4 complement | `g_N <= (2LT/pi) S_N(TL)`, `S_N(x) = sum_{n>=N}(2n+1)j_n(x)^2` | **CERTIFIED** |
| 5 block bound | `lam_min(R_T) >= (lam_A + lam_D - sqrt((lam_A-lam_D)^2 + 4 kappa^2))/2` | **PROVED** |

Certified constants at `L = 4/5`, `T = 200` (`certs/baseline_constants.json`),
`M = sup_{|t|<=T}|Psi_L - beta| <= 8.83262267769`, `|p| = sqrt(2 sinh L) = 1.33274602396`:

| `N` | `g_N` | `|p_perp|` | `lam_D` | `kappa` | **`lam_A` required** |
|---|---|---|---|---|---|
| 200 | `2.2755e-17` | `3.23e-515` | `0.51846573` | `4.2133e-8` | `3.4240e-15` |
| **220** | `6.7424e-31` | `1.17e-575` | `0.51846573` | `7.2527e-15` | **`1.0146e-28`** |
| 250 | `2.4380e-55` | `8.87e-668` | `0.51846573` | `4.3612e-27` | `3.6685e-53` |

**The complement and the coupling are not the obstruction.** They were the parts
most likely to fail, and they are certified to be `1e-28` and smaller, against a
window minimum believed to be `~9e-18` — ten orders of margin.

## 7.3 What is missing, exactly

Only `lam_A`: the least eigenvalue of the `N x N` compression of `R_T` to the
first `N` Legendre modes, with

```
   (R_T)_{nm} = beta delta_{nm} + 2 s p_n p_m
                + (2L/pi) sqrt((2n+1)(2m+1)) (-1)^{(n-m)/2}
                  int_0^{T} (Psi_L(t) - beta) j_n(tL) j_m(tL) dt .
```

**Blocker class: resource / implementation, not mathematics.** The integrand is
analytic in the strip `|Im t| < 1/2` (the nearest digamma pole is at `t = i/2`),
so a panelled Gauss–Legendre rule converges geometrically; with 800 panels of
width `1/4` and 32 nodes the Bernstein-ellipse error is below `1e-50`. What was
not built here is (i) certified Gauss–Legendre nodes and weights in Arb,
(ii) the rigorous panel error bound, and (iii) a certified `LDL^T` positivity
check of the resulting `110 x 110` even-sector matrix at the `1e-28` level. Each
is standard; together they are a substantial build.

Measured costs that informed this: a single entry via Arb's rigorous
`acb.integral` over `[0,200]` took **27–32 s**, so per-entry rigorous
integration is `~50 h` for one sector and was rejected in favour of a fixed rule
with an error bound, which was not completed.

**Not attempted, deliberately:** `L = 0.81`, per the accepted sequencing (R5).

## 7.4 Local-occupation pilot

**No occupation bound was consumed in any certificate, and none can be at this
margin.** Two findings, both decisive:

1. **At the gate cutoff the unfavourable set is empty.** `beta := H(T) - A_L`
   and `H` increasing give `Psi_L(t) >= beta` for every `|t| >= T`, so
   `E = {}` and `b q = 0`. There is nothing for an occupation bound to buy at
   the baseline.
2. **For a smaller cutoff the penalty cannot be paid.** `b <= M = 8.8326`, and
   the margin that must absorb `b q` is `lam_A ~ 1e-17`. So `q <= 1.13e-18` is
   required. Lemma A needs `meas(E) <= 3.93e-18`; Lemma B needs
   `eps <= 1e-18` in *every* cell. A Logvinenko–Sereda/Kovrijkine thick-set
   bound has the same shape and cannot reach `1e-18` on a set of positive
   density either.

The lemmas are correct and informative in general (`q_B = 0.09…0.26` against
measured ratios `0.05…0.10`); they are simply ~18 orders of magnitude from what
this particular comparison needs. That is a property of the margin, not of the
lemmas.

## 7.5 Phase 2 next step

**Build the certified quadrature and the `LDL^T` check.** Concretely: certified
Gauss–Legendre nodes/weights in Arb, the Bernstein-ellipse panel error bound for
an integrand analytic in `|Im t| < 1/2`, assembly of the `110 x 110` even-sector
matrix and the odd-sector counterpart, and an interval `LDL^T` factorization of
`R_T - m I`. The certified constants above then convert any `m >= 1.0146e-28`
into a whole-window lower bound at `L = 4/5`.

---

# 8. PHASE 3 — the finite blocks, assembled and certified

> ## **CERTIFIED: `Q_zeta(f) >= 1e-17 ||f||_2^2` for every admissible `f` supported in `[-4/5, 4/5]`, in both parity sectors, hence for arbitrary complex tests.**

Obtained independently. The external `8.9e-18` was never an input; it is a
comparison point only, and this build's own certified **upper** bound
`2.112992612e-17` brackets the result from the other side.

## 8.1 The certified chain

| quantity | value | what it bounds |
|---|---|---|
| `beta = H(T) - A_L` | `0.5184657332438467199873536` | the scalar symbol: `Psi_L(t) >= beta` for all `|t| >= T = 200`; the unfavourable exterior set is **empty** |
| `M` | `8.832622677693132757724221` | `sup_{|t|<=T}|Psi_L - beta|`, from monotonicity of `H`, no covering |
| `S_N(TL)` | `6.619351600478776e-33` | `sum_{n>=220}(2n+1) j_n(160)^2` |
| `g_N` | `6.742416174588454e-31` | Fourier energy fraction in `[-T,T]` of `V_N^perp` — a property of the **basis complement** |
| `\|p_perp\|` | `1.169181247302579e-575` | pole-vector tail |
| `kappa` | `7.252662296635726e-15` | the coupling `\|\|B\|\|` |
| `kappa^2` | `5.260111038904140e-29` | — |
| **`a`** | **`1.03e-17`** | least eigenvalue of the **exact compression of `R_T` to `V_N`** — not of `Q`, not of the infinite operator |
| `d` | `0.5184657332438467` (even); odd differs by `< 1.6e-396` | `R_T` on `V_N^perp` |
| `kappa^2/d` | `1.014553267772122e-28` | the Schur penalty (break-even point, **not** the bound) |
| **`mu`** | **`1e-17`** | the certified whole-window constant |

Final test, per sector, by three **strict interval signs** rather than the
cancellation-prone `2x2` formula: `a > mu`, `d > mu`,
`(a-mu)(d-mu) > kappa^2`, with margin `1.55539719921e-19`.

`N = 220` means Legendre degrees `0..219`: **110 modes per parity**. The odd
sector carries the pole term with `s = -1` (subtracted in `d`, not dropped).

**Parity and complex tests.** The form decouples by parity: `C_{nm} = 0` for
`n + m` odd because `v_n v_m` is then odd on the symmetric interval, and the
pole term is `2 F(i/2)F(-i/2) = 2C^2 - 2S^2` with `C = int f cosh(u/2)`,
`S = int f sinh(u/2)` — the even and odd parts separately. For complex
`f = u + iv` the Hermitian form gives `Q(f) = Q(u) + Q(v)`, and each real part
splits into its parities. So the two sector bounds cover the full class, and the
whole-window constant is the smaller, here equal in both.

## 8.2 Assembly

Panelled Gauss–Legendre: `800` panels of half-length `1/8` on `[0,200]`,
**order 64**, nodes and weights from FLINT's `arb.legendre_p_root(n,k,weight=True)`.
Quadrature remainder by the Bernstein-ellipse bound with imaginary semi-axis
`b_e = 0.30` (strictly inside the pole-free strip `|Im z| < 1/2`), `rho = 5.0`,
folded into the entries as a widening **once**. Max entry radius `4.99e-60`.

Holomorphic extensions, as required: `H_hol(z) = [psi(1/4+iz/2)+psi(1/4-iz/2)]/2
- log pi` (`Re psi(1/4+iz/2)` is **not** holomorphic), and the basis product
through the entire `v_n(z) = sqrt(2L(2n+1)) j_n(Lz)`.

Positive definiteness by **Rump's preconditioned criterion**
`||I - R^{-T}(A - aI)R^{-1}||_inf < 1`, giving `2.3470257e-25` (even) and
`4.5876981e-33` (odd).

Elapsed `1343 s`, peak RSS `205.5 MiB`.

## 8.3 Four defects found and fixed during this run

1. **Wrong pole phases.** The quadratic form is assembled in `y`-coordinates
   (`y_n = (-1)^{floor(n/2)} x_n`), so the pole vector needs the same phase.
   Omitting it is invisible at degree 0 (`sigma_0 = +1`) and wrong from degree 2
   on; it produced a genuinely **negative pivot** `-0.3949` and would have been
   reported as indefiniteness. Caught by comparing `|F(t)|^2` and `F(i/2)`
   against direct integration — now check X1.
2. **`acb.digamma` returns NaN** on balls whose real part has radius `>~ 0.1`,
   even strictly inside the right half-plane. Replaced by the explicit
   recurrence `psi(z) = psi(z+K) - sum_{k<K} 1/(z+k)` with adaptive `K`.
3. **`acb.bessel_j` is useless on wide balls** — width `~1e32` at `|w| ~ 80`,
   not fixable by precision (checked at 900 and 3000 bits). Replaced by the
   uniform bound `|v_n(z)| <= sqrt(2L) e^{L|Im z|}` from the Legendre integral
   representation, valid for all `n`.
4. **Naive interval `LDL^T` is unusable here**: entry radii `~1e-38` amplify by
   about `6.5` per pivot, reaching `O(1)` by step 47. Preconditioning replaces
   that by a factor of order `cond(R)`.

Also: FLINT's adaptive integrator **silently swallowed** an exception from the
`1/sqrt(w)` representation and returned a wrong value rather than failing. The
cross-check is therefore run only away from the origin, with the
non-analyticity protocol wired explicitly.

## 8.4 What this is and is not

* It **is** a whole-window lower bound at `L = 4/5`, consuming the full analytic
  reduction in the repository normalization (explicit formula, pole,
  archimedean, prime comb) — not merely nonnegative coefficients of `log F`.
* It is **not** RH, and it says nothing about any other window. There is no
  propagation theorem here, and Zhu's Theorem 1.4 argues the route is capped.
* The repository's explicit-formula normalization (`W4`) remains an **imported**
  hypothesis; this build verified it only through the independent cross-check of
  `Q[0,0]` against the repository's own closed form (agreement to 20 digits plus
  a consistent exterior remainder).
* No new arithmetic estimate was obtained.

---

# 9. PHASE 4 — audit of the certificate, and the measured reduction loss

**The lower certificate survives the audit.** Nothing in §8 is retracted. Two
*interpretive* claims made alongside it are corrected (`RETRACTIONS.md` R11, R12),
and one new diagnostic changes how the result should be read.

## 9.1 Headline table

| quantity | value | status |
|---|---|---|
| `lambda_Q` lower bound | **`1e-17`** | **CERTIFIED** (§8), re-audited here |
| `lambda_Q` upper bound | **`2.112992612e-17`** | **CERTIFIED** (phase 1 artifact), retained unchanged |
| two-sided enclosure | `1e-17 <= lambda_Q <= 2.112992612e-17` | **CERTIFIED**, width `1.1129926e-17` |
| `r_f` on the frozen even-block minimizer | `1.031022030245990e-17` | **CERTIFIED** — upper bound for `lambda_R`, and for `lambda_A` in this block |
| `e_f` on the same vector | **`>= 1.530372869310224e-17`** | **CERTIFIED** lower bound, band `[200,600]` only |
| `e_f / r_f` | **`>= 1.4843`** | **CERTIFIED** |
| `q_f = r_f + e_f` | `>= 2.561394900e-17` | derived; **does not improve** the upper bound, which is retained |
| exception-handling audit | characterisation **withdrawn**, R12 | resolved |

## 9.2 The three spectral quantities, kept apart

```
   lambda_Q = inf Q(f)/||f||^2        over the admissible whole window
   lambda_R = inf R_{T,beta}(f)/||f||^2   over the same window
   lambda_A = min eigenvalue of the finite compression of R in ONE sector
```

`Q >= R` pointwise and an infimum over a subspace is `>=` the infimum over the
whole space, so

```
   lambda_R <= lambda_Q        and        lambda_R <= lambda_A ,
```

with **no** inequality between `lambda_A` and `lambda_Q`. The numerical even-block
value `1.0310220e-17` estimates **`lambda_A`**; §8 implied it was comparable to
`lambda_Q`, which it is not. The certified `r_f = 1.031022030245990e-17` is the
same quantity done properly.

## 9.3 The measured loss, and what it means

On the minimizer of the even block — the vector most favourable to the reduction —

```
   r_f = 1.031022030245990e-17 ,      e_f >= 1.530372869310224e-17 ,
```

so **the reduction discards at least 1.48 times as much as it keeps**, and that
is from the band `[200,600]` alone; the true `e_f` is larger. Consequently
`q_f >= 2.561e-17`, which is **above** the certified upper bound
`2.112992612e-17`, so this vector is a poor test function for `Q` and the earlier
upper bound stands.

Read correctly: `e_f` is the loss **on this one vector**, not the gap between
`lambda_R` and `lambda_Q`. But it does show the frequency-envelope reduction is
lossy at exactly the place one would hope it were tight, which is the useful
diagnostic this run was for.

## 9.4 The preconditioned positive-definiteness proof, stated exactly

With `B = A - aI` Hermitian and `R` the (frozen, treated as exact input) upper
triangular Cholesky factor of the midpoint of `B`, set `C = R^{-1}` and

```
   G = C^* B C = R^{-T} (A - aI) R^{-1} ,
```

computed by two explicit forward substitutions in ball arithmetic. The verified
inequality is

```
   || I - G ||_infinity < 1 ,
```

giving `2.3470257e-25` (even) and `4.5876981e-33` (odd) at `a = 1.03e-17`.

*Why this proves `A - aI > 0`.* `||I-G||_inf < 1` puts every eigenvalue of the
Hermitian `G` in `(0,2)`, so `G > 0`. If `C` were singular then `G = C^*BC`
would be singular, contradicting `G > 0`; so `C` is nonsingular. Then for any
`x != 0`, writing `x = Cy` with `y = C^{-1}x != 0`, `x^*Bx = y^*Gy > 0`. Hence
`B > 0`, i.e. `lambda_min(A) > a`. This is a **congruence**, not a one-sided
inverse residual — the latter would establish nonsingularity only.

## 9.5 Exception handling — resolved, and a claim withdrawn

`acb.integral` is used in **exactly one** place, `rigorous_entry`, called only by
the verifier's spot-check; `block.assemble` never calls it (grep of the function
body, check A5d/A5e). The certificate therefore never depended on that path.

The minimal reproducer (`phase3/src/callback_repro.py`) shows that on this
version every tested path is safe: a raising callback propagates as `ValueError`,
an invalid operation returns `nan`, and the documented non-finite protocol works.
The real residual issue is that an escaping Python exception surfaces as a
`SystemError` only when control returns to the interpreter — a **binding-level**
matter in the callback adapter, since FLINT's C core has no knowledge of Python
exceptions. **My earlier claim that the integrator returned a silently wrong
value is withdrawn** (R12): the wrong numbers came from a run whose matrix side
still carried the pole-phase bug, at a different `T`, and I conflated interleaved
output from two runs.

The wrapper is hardened regardless: the callback records any failure in a flag,
and the integration is rejected on that flag rather than on the returned value.

## 9.6 Next experiment, from the measured loss

`e_f/r_f >= 1.48` says the exterior discard — not the finite block — is where the
reduction is losing. The cheap, bounded follow-up is to **re-run the existing
pipeline at a larger `T` at fixed `L = 4/5`**, which shrinks the discarded
exterior while leaving the gate valid (`beta = H(T) - A_L` only grows with `T`,
and by R11 that is favourable), at a cost of roughly `N ~ e L T / 2` more modes.
Measuring `r_f` and `e_f` at `T = 400` against the values above would say
directly whether the enclosure width `1.11e-17` is reducible by this route or
whether the loss is intrinsic.

Explicitly **not** recommended: lowering `beta` (refuted, R11), `L = 0.81`, or
any new occupation/out-of-band programme.

---

# 10. PHASE 5 — T = 400, and the original form in the time domain

> ## **IMPROVED BRACKET: `1.41e-17 <= lambda_Q <= 1.90138291114546440e-17`** at `L = 4/5`
>
> **Both** endpoints improved. Lower raised from `1e-17` by the `T = 400`
> certificate (+41%); upper lowered from `2.112992612e-17` by the Rayleigh
> quotient of `Q` at the new `T = 400` minimizer, evaluated through the compact
> time-domain route. Width `4.9138e-18`, down from `1.11299e-17` —
> **56% narrower**.

## 10.1 Why `T = 200` could not be pushed further

`lambda_R(200) <= r_200(f_200) = 1.0310220302459896493e-17` is a certified
*ceiling* on the unchanged comparison: a lower bound justified solely by
`Q >= R_200` can never exceed it. Against the achieved `1e-17` that left only
**3.1022%** of room. This is a limit on the comparison, **not** on `lambda_Q`.

## 10.2 The cutoff comparison (proved)

With `Psi_L >= beta_i` on `|t| > T_i`, `T_2 >= T_1`, `beta_2 >= beta_1`:

```
   R_{T2,beta2}(f) - R_{T1,beta1}(f)
     = (1/2pi) int_{T1<|t|<=T2} (Psi_L - beta_1)|F|^2
     + ((beta_2-beta_1)/2pi) int_{|t|>T2} |F|^2      >= 0 ,
```

both terms nonnegative. So `lambda_R(T_2) >= lambda_R(T_1)` **for the exact
forms**. It does not follow that a certified numerical bound improves, and it
does **not** order finite matrices of different sizes: changing `N` changes the
compression, so only matrices in the same subspace are comparable.

Certified: `beta_400 = H(400) - A_L = 1.211613695058065`, versus
`beta_200 = 0.518465733243847`; the increase is `0.693147961814` (`~ log 2`),
and `beta_400 > beta_200` strictly by the proved monotonicity of `H`.

## 10.3 The `T = 400` certificate

`N = 440` was chosen from the **actual** inequalities, not the heuristic
`e L T / 2 = 434.9`: at `N = 440`, `S_N(320) = 9.26e-63`, `g_N = 1.886e-60`,
`kappa = 1.308343771e-29`, so the chain closes for any `a > 1.41e-58`.
Degrees `0..439`, **220 modes per parity**. 800 panels of half-length `1/4`,
Gauss-Legendre order 64.

| sector | `a` (certified) | `d` | `kappa^2/d` | `mu` | margin |
|---|---|---|---|---|---|
| even | `1.42e-17` | `1.2116136950581` | `1.413e-58` | `1.41e-17` | `1.21161369506e-19` |
| odd | `2.4e-17` | `1.2116136950581` | `1.413e-58` | `1.47e-17` | `1.12680073640e-17` |

Rump residuals `||I-G||_inf`: `1.56889e-15` (even, at `a = 1.3e-17`),
`4.35257e-22` (odd, at `a = 2.4e-17`). **Whole-window `mu = min = 1.41e-17`**;
the even sector binds. Assembly `3532 s`, total `3581 s`, peak RSS `393.4 MiB`.

The even block's own minimizer gives `R_400(f_400) = 1.4330724047355850e-17`, so
`lambda_A(400, even)` is pinned in `[1.42e-17, 1.4331e-17]` and the `T = 400`
route has only about **1.6%** of room left above the certified `1.41e-17`.

## 10.4 The original form in the time domain

`Q` is now evaluated directly on the geometric side, with **no long oscillatory
frequency integral**, via the autocorrelation `g` and the compact archimedean
identity (Zhu v2 Lemma 2.5) in the repository normalization:

```
   Arch(g) = -[gamma_E + log pi + log(1-e^{-4L})] g(0)
             + int_0^{2L} 2[e^{-2x}g(0) - e^{-x/2}g(x)]/[1-e^{-2x}] dx ,
   Q(f)    = 2 s F(i/2)^2 + Arch(g) - sum_{log n < 2L} (2 Lambda(n)/sqrt n) g(log n) .
```

Two things make this cheap and rigorous. `g` is a polynomial of degree
`<= 2 deg f` on `[0,2L]`, so the inner integral is computed by Gauss-Legendre
with `deg f + 2` nodes and is **exact** — no inner quadrature error at all. And
the removable endpoint is handled by an exact identity rather than a limit
estimate: for **even** `f`,

```
   g'(x) = -f(L) f(L-x) + int_{-L}^{L-x} f'(s+x) f(s) ds ,
```

and at `x = 0` the integral is `(1/2)[f^2]_{-L}^{L} = 0`, so
**`c = -g'(0+) = f(L)^2`** exactly, and the integrand limit is
`c - (3/2) g(0)`. The outer integrand is analytic in `|Im x| < pi`; panels are
geometric near `0` so no Bernstein ellipse encloses the origin, and the
`[0, 1e-40]` head is bounded by an explicit constant.

**Validation.** On `f = phi_0` (the constant) the compact route gives
`0.080274054959269611214 +/- 3.09e-22`, matching the repository's independent
closed-form `Q[0,0] = 0.080274054959269611214` to all 21 digits.

**Disclosure.** This is a *second representation* of the same form, not an
independent implementation: it shares the Arb library, `block.Symbol` for the
comb data and `block.pole_vector`.

## 10.5 The 2x2 table (all entries CERTIFIED)

| | `R_200` | `R_400` | `Q` (compact route) |
|---|---|---|---|
| `f_200` | `1.0310220302459896493e-17` | `2.3520530467310079487e-17` | `3.64859175290407702e-17` |
| `f_400` | not evaluated (see below) | `1.4330724047355850300e-17` | **`1.90138291114546440e-17`** |

Same-vector discarded contributions, both certified:

```
   e_200(f_200) = Q - R_200 = 2.617569722658087374e-17
   e_400(f_200) = Q - R_400 = 1.29653870617306907e-17
```

The discard on the old vector **halves** when the cutoff doubles, exactly the
sign the cutoff identity predicts, and `R_400(f_200) > R_200(f_200)` as it must.
The earlier band-only bound `e_200 >= 1.5303728693e-17` is confirmed: it
captured `58.5%` of the true tail.

`R_200(f_400)` was **not** evaluated: `f_400` has degrees up to `438` and the
`T = 200` block only spans degrees `0..218`, so truncating it would change the
function. Evaluating `R_200` there would need a separate assembly on the larger
basis and was not worth the cost for a diagnostic.

Same-vector discards on the two minimizers:

```
   e_200(f_200) = 2.617569722658087374e-17        e/r = 2.539
   e_400(f_400) = 4.6831050640987937e-18          e/r = 0.327
```

`q(f_200) = 3.64859175290407702e-17` is **above** the old upper bound
`2.112992612e-17`, so that vector does not improve it — the expected outcome.
But `q(f_400) = 1.90138291114546440e-17` is **below** it, so the upper endpoint
moves. A trial function's `Q` value bounds `lambda_Q` above only when it is
smaller than the incumbent; here `f_400` is, and `f_200` is not.

## 10.6 Status, exactly

| | old | new | status |
|---|---|---|---|
| whole-window lower bound | `1e-17` | **`1.41e-17`** | CERTIFIED, `T=400`, both sectors |
| best original-`Q` upper bound | `2.112992612e-17` | **`1.90138291114546440e-17`** | CERTIFIED, Rayleigh at `f_400` |
| bracket width | `1.11299e-17` | **`4.9138e-18`** | 56% narrower |
| `e(f_200)` | `>= 1.5304e-17` (band) | `= 2.617569722658087374e-17` (exact) | CERTIFIED |
| `e_400(f_200)` | — | `1.29653870617306907e-17` | CERTIFIED, halved |
| `e_400(f_400)` | — | `4.6831050640987937e-18` | CERTIFIED |
| complement `g_N` / coupling `kappa` at `T=400` | — | `1.886e-60` / `1.308343771e-29` | CERTIFIED |
| even/odd margins at `mu` | — | `1.21161e-19` / `1.12680e-17` | CERTIFIED |
| time / memory | — | `3581 s`, `393.4 MiB` (assembly); `1543 s` (cross table) | measured |
| unresolved | — | `R_200(f_400)`; whether `a` can exceed `1.42e-17` | see below |

**Status: improved bracket.** Not an intrinsic-barrier claim anywhere, and
nothing here bears on RH, zero locations, larger support, or any arithmetic
estimate.

## 10.7 Next

The `T = 400` route is nearly exhausted: `lambda_A(400, even)` is pinned in
`[1.42e-17, 1.4331e-17]`, so at most about **1.6%** more is available from this
cutoff, and the even sector binds. The residual gap is now between
`lambda_R(400) <= 1.4331e-17` and `lambda_Q <= 1.90138e-17`, i.e. the
reduction's own loss, which `e_400(f_400) = 4.68e-18` measures directly.
`T = 800` would be the same experiment again at roughly 4x the assembly cost;
whether it is worth running should be decided from whether `e` keeps halving.
