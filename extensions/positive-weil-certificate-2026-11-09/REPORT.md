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

> **2026-11-10, FINAL FOR THIS BUILD (`RETRACTIONS.md` R14 + addendum).**
> Published bracket: **`1.41e-17 <= lambda_Q(4/5) <= 1.674e-17`**, width
> `2.64e-18`. Upper endpoint `1.67312553407068186807726837736e-17` from the
> CLOSED-FORM route (digamma / Lerch Phi / Hurwitz zeta, no quadrature) at
> **160 even cosine modes** (`phase5/src/upper_highN.py`). This BEATS the
> withdrawn time-domain endpoint `1.9014e-17` by 12.0%. Against Zhu's certified
> `[8.9e-18, 2.27e-17]` the width is 80.9% smaller (factor 5.23), and his
> numerical reference `1.656e-17` lies inside. Parity margin `6.2687e-18`.
>
> **ERRATUM, 2026-11-10 (`RETRACTIONS.md` R14) — THE UPPER BOUND CHANGED.**
> The compact time-domain evaluator omits the TRUNCATION error of its outer
> panelled Gauss-Legendre rule (it bounds only the `[0,eta]` head), so every
> value it produced is **withdrawn as a certificate**: `q(f_400)`, `q(f_200)`,
> the endpoint `1.9014e-17`, and the discards `e_T = q - r_T` with their ratios.
> They are retained as MEASUREMENTS. The published bracket is now
> **`1.41e-17 <= lambda_Q(4/5) <= 2.1130e-17`**, the upper endpoint coming from
> the phase-1 CLOSED-FORM route (`src/certify.py`: digamma / Lerch Phi /
> Hurwitz zeta, no quadrature anywhere), replayed under owned precision in
> `phase5/src/replay_upper.py`. The parity separation survives on that witness
> with margin `1.870073882785092e-18`. **Unaffected:** the lower bound, the
> odd-sector bound, every constant in the chain, the scalar `R_T` values (whose
> routine DOES add a certified remainder) and the ceilings built from them.
>
> **ERRATUM, 2026-11-10 (`RETRACTIONS.md` R13).** Earlier revisions of this file
> printed the bracket's upper endpoint as `1.90138291114546440e-17`. That is the
> **midpoint** of the `q(f_400)` enclosure, not an upper bound for it, and
> `phase5/certs/parity.json` separately recorded a looser endpoint
> `1.9013829111454646892e-17` produced by a 53-bit re-parse. The strict endpoint
> is `1.901382911145464406030000e-17`. Both errors were conservative and nothing
> was withdrawn *by R13*. **Note, however, that R14 above later withdrew that
> endpoint entirely**: the quantity it bounds is itself uncertified, so the
> instruction to "quote `1.9014e-17`" given here is superseded — the bracket to
> quote is `1.41e-17 <= lambda_Q(4/5) <= 1.674e-17`. A technical note
> consolidating this build is in [`note/`](note/).

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

> ## ~~**IMPROVED BRACKET: `1.41e-17 <= lambda_Q <= 1.9014e-17`**~~ — **SUPERSEDED**
>
> **This phase's upper endpoint is WITHDRAWN (`RETRACTIONS.md` R14): it came from
> the compact time-domain route, whose outer quadrature remainder is unbounded.
> The current bracket is `1.41e-17 <= lambda_Q(4/5) <= 1.674e-17`. The lower
> endpoint below stands. Section retained as the record of what this phase did.**
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
But `q(f_400)` has midpoint `1.90138291114546440e-17` and endpoint
`1.901382911145464406030000e-17`, both **below** it, so the upper endpoint
moves. A trial function's `Q` value bounds `lambda_Q` above only when it is
smaller than the incumbent; here `f_400` is, and `f_200` is not.

## 10.6 Status, exactly

| | old | new | status |
|---|---|---|---|
| whole-window lower bound | `1e-17` | **`1.41e-17`** | CERTIFIED, `T=400`, both sectors |
| best original-`Q` upper bound | `2.112992612e-17` | ~~`1.9014e-17`~~ | **WITHDRAWN, R14** — compact route, remainder unbounded. Current: `1.674e-17`, closed-form route |
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

---

# 11. PHASE 6 — precision ownership, parity separation, scalar pilots

> **SUPERSEDED (`RETRACTIONS.md` R14): the upper endpoint quoted here is withdrawn;
> the current bracket is `1.41e-17 <= lambda_Q(4/5) <= 1.674e-17`.** As written at
> the time: bracket unchanged and re-verified at `1.9014e-17`.
> New: a certified **parity separation** identifying the even sector as the one
> that determines `lambda_Q`, a repaired precision discipline, the missing
> `R_200(f_400)` entry, and a measured cost basis for the `T = 800` decision.

## 11.1 Precision ownership, repaired at the source

`flint.ctx.prec` is a process global. Routines set it on entry and never
restored it, so a callee silently changed its caller's context. **This defect
appeared four times in this project, three of them inside this session**, and
each time it produced a plausible-looking wrong constant rather than an error:

| where | symptom |
|---|---|
| phase 2 `M_sup` | wide-ball digamma returned NaN, the comparison skipped those cells, `M` came out `1831` instead of `8.83` — and was **used** |
| phase 4 `verify_phase3` | `validate_conventions` reset to 300 bits; the complement constants came back `[+/- 7.4e+132]` |
| phase 5 `run_block400` | `verify_pd_preconditioned` left 1400 bits; the EVEN sector's `sector_dk` ran at 900 and gave `d = [+/- 2.8e+183]` |
| phase 6, three times in this run's own code | `arb(4)/5` built at 53 bits → integral radius `1e105`; coefficients built at 53 bits → `1.8e-15`; the reference constant `arb("0.08027...")` built at 53 bits → a meaningless comparison |

The last three are the important lesson: scoping alone is not enough. **An input
constructed outside the scope carries the ambient precision's radius, and
entering a higher-precision scope afterwards does not recover it.** The rule
adopted is therefore *own your inputs as well as your scope*:

* `phase5/src/prec.py` — `owned(bits)` wraps `ctx.workprec` and **asserts** the
  caller's context was restored (on normal return and on exception), `owns()` is
  the decorator form and refuses an implicit precision, `PrecCache` keys every
  entry by precision so a 300-bit value can never be served to a 2000-bit
  caller, and `require_tight` turns a wide or non-finite enclosure into an
  explicit **UNRESOLVED**, never a midpoint read.
* `phase5/src/constants.py` — every precision-dependent certificate constant
  (`beta, M, S_N, g, kappa, d_even, d_odd`) is built inside `owned(prec)` and
  cached by `(L, T, N, tail_extra, prec)`.
* `phase5/src/scalar.py` — rebuilds `L`, the `Symbol` and the coefficient vector
  **inside** the scope from exact rational data, and checks `beta` for tightness.

`phase5/test_prec.py` (13 checks): nesting, restoration after a controlled
exception, rejection of an invalid precision, identical enclosures from ambient
53/300/1200/2000, order-independence, precision-keyed caching, fail-closed on
wide and non-finite values, and — the regression that matters —
`constants(..., prec=200)` is now **refused** rather than silently wrong.

**Why the cached matrices did not need reassembly.** `block.assemble` never
*reads* an ambient precision: it sets `prec_bessel` on entry and `prec_mat`
around each flush, so every entry was produced at a precision the assembly
chose. (It leaked that precision outward, which is what corrupted the constants
computed afterwards, but it never inherited one.) Entries are stored as Arb
balls with explicit radii and re-read by `arb(str)`, which is exact. Only the
dependent constants were affected, and those are recomputed here at an owned
precision. The `T = 400` certificate was replayed through the repaired path:
even `||I-G||_inf = 1.4973903e-14` at `a = 1.42e-17`, odd `4.3535838e-22` at
`a = 2.4e-17`, both PD, `mu_even = 1.41e-17` retained with margin
`1.21161369506e-19`.

## 11.2 Certified parity separation

Re-evaluated from the exact certified data, not the rounded report values:

```
   odd sector:   a = 2.4e-17 ,  d = 1.2116136950580645495 ,  kappa = 1.308343771e-29
   mu_odd = 2.3e-17 :  a > mu,  d > mu,  (a-mu)(d-mu) > kappa^2
                       margin = 1.21161369506e-18
```

so `mu_odd = 2.3e-17` certifies — the `1.47e-17` used in §10 was unnecessarily
conservative (an artefact of the ladder I searched, not of the data). With
`U` the exact upper endpoint of the `q(f_400)` enclosure,

```
   U = 1.901382911145464406030000e-17 ,
   mu_odd - U = 3.986170888545355939700000e-18 > 0 ,
```

and since `f_400` is **even**,

```
   inf_even Q/||f||^2  <=  U  <  mu_odd  <=  inf_odd Q/||f||^2 .
```

With the orthogonal parity decomposition this **identifies the even infimum as
the whole-window infimum**. It does *not* prove simplicity, uniqueness, or that
the infimum is attained.

**Corrected 2026-11-10, after obtaining Zhu's full text.** Saying only that
"parity separation is not new" was too weak. `arXiv:2608.24827v2`
**Theorem 6.2** certifies the separation at support `1.6`, i.e. at `L = 4/5`,
with `8.206e-15 <= lam_odd_1 <= 2.347e-14` — an odd-sector constant about
**357x stronger** than the `2.3e-17` obtained here — together with
`lam_even_2 >= 2.085e-12`, and concludes that the ground state is **simple and
even**. What is recorded above is only that this pipeline's own certified
numbers reproduce the ordering. The global lower bound
remains `1.41e-17`, governed by the even sector.

## 11.3 Original-form validation above degree 0

The `phi_0` check cannot see a degree-2 phase error (`sigma_0 = +1`), which is
exactly how the earlier pole-phase defect survived. Added (`phase5/src/conventions.py`,
13 checks, all passing), with basis coordinates stated and the pole moment computed
**directly from the physical-space function** rather than via `block.pole_vector`:

| case | pole moment (direct integral) | agrees with helper |
|---|---|---|
| `phi_0` | `1.298912903748447` | yes, to `1e-25` |
| `phi_2` | `0.03051622422538102` | yes |
| `phi_0 - phi_2` | `1.268396679523066` | yes |

The prime shifts are re-enumerated independently of `block.Symbol` and agree.
`Q(phi_0)` still matches the repository's exact `Q[0,0]` (difference
`-3.0834062e-22`).

**Exactness, stated precisely.** For `f` of degree `d` the integrand
`f(t) f(t-x)` has degree `2d` in the integration variable, so an `m`-node
Gauss-Legendre rule is exact when `2m-1 >= 2d`, i.e. `m >= d+1`;
`compact.autocorr` uses `m = d+2`. For `d = 438` that is `440` nodes. Exactness
means **no polynomial quadrature truncation error in exact arithmetic** — node,
weight, coefficient and rounding enclosures are still carried as ball radii, and
the outer archimedean integral keeps its own certified Bernstein remainder.

All decisive numbers are published as midpoint-radius enclosures in
`phase5/certs/`; printed decimals are summaries, not exact values.

## 11.4 The completed scalar table

`R_200(f_400)` was the missing entry. Evaluated as a **direct scalar integral**
on the full degree-438 function — no truncation, no enlarged matrix:

```
   R_200(f_400) = 1.24275215e-17  +/- 8.10e-26     [CERTIFIED]
                  (quadrature remainder <= 7.75e-26, 200 panels, order 64, 933 s)
```

| | `R_200` | `R_400` | `Q` |
|---|---|---|---|
| `f_200` | `1.0310220302459896493e-17` | `2.3520530467310079487e-17` | `3.64859175290407702e-17` |
| `f_400` | **`1.24275215e-17`** | `1.4330724047355850300e-17` | `1.90138291114546440e-17` |

Cutoff monotonicity now checked on **each fixed function separately**, as
required: `R_400(f_200) > R_200(f_200)` and `R_400(f_400) > R_200(f_400)`.
The scalar route was validated against the matrix on `f_200`:
`1.03102203e-17 +/- 1.95e-26`, agreeing to `1e-24`.

## 11.5 What `q_f - r_f` is, and is not

For `f_400`, `q_f - r_f = 4.6831050640987937e-18` is **the discarded-term
identity evaluated on that function**, within its certified enclosure. It is
**not** `lambda_Q - lambda_R400`. Both `lambda_R400 <= r_f` and
`lambda_Q <= q_f` are upper bounds, and two upper bounds do not make their
difference an equality for the difference of the infima.

The valid fixed-reduction ceiling to retain: `lambda_R400 <= r_f = 1.4330724e-17`,
which leaves about **1.636%** of headroom above the certified `1.41e-17` for any
certificate based solely on the unchanged `R_400`.

## 11.6 `T = 800`: scalar pilot, measured cost, and the decision

`beta_800 = H(800) - A_L = 1.90476107093077689 > beta_400 = 1.2116136950580645`,
the increase again `0.693147...` (`~ log 2`) — certified, and the ordering the
cutoff identity needs.

**Pilot status: COMPLETE** (12634 s total; it finished after §11 was first
written, and this section is the finished result).

| | `R_200` | `R_400` | `R_800` | `Q` |
|---|---|---|---|---|
| `f_200` | `1.0310220302459896e-17` | `2.3520530467310079e-17` | `3.01603e-17` | `3.6485917529040770e-17` |
| `f_400` | `1.24275215e-17` | `1.4330724047355850e-17` | `1.6681e-17` | `1.9013829111454644e-17` |

Strictly increasing in `T` along **each** fixed function, as the cutoff identity
requires. Discarded contributions and their ratios, by rigorous interval
division:

```
   f_400 :  e_200 = 6.586308e-18   e_400 = 4.6831050640988e-18   e_800 = 2.333e-18
            e_400/e_200 = 0.711          e_800/e_400 = 0.498 +/- 2.36e-4
   f_200 :  e_400 = 1.2965387061731e-17  e_800 = 6.3256e-18
            e_800/e_400 = 0.48788 +/- 9.07e-6
```

**There is no decay law here.** On `f_400` the successive ratios are `0.711`
then `0.498` — not a constant rate. Two functions happening to give `~0.49` over
the same doubling is two samples of one step, not a rate in `T`, and comparing
differently optimised vectors would establish even less. `R_800` on a fixed
vector is a variational **upper** bound for `lambda_R800`; it yields no
whole-window lower bound and does not predict what a reoptimised `f_800` does.

**The decisive number is the ceiling.** `f_400` is admissible, so

```
   lambda_R(800)  <=  R_800(f_400)  =  1.6681e-17 ,
```

and any certificate resting on `Q >= R_800` is bounded by that. Against the
current lower endpoint `1.41e-17` the entire `T = 800` route is therefore worth

```
   at most a factor 1.1830 , i.e. +18.3% ,
```

and the reoptimised `lambda_R(800)` can only be smaller than this ceiling, so the
realisable gain is less. This replaces the `~35%` figure estimated in §10.7 from
`lambda_Q` alone, which was an overestimate: it bounded the gap to `lambda_Q`,
not to what the `T = 800` *comparison* can reach.

That overrun is itself the useful measurement: **the per-Bessel cost grows with
the argument**, so node count alone understates the scaling. Measured basis:

| stage | measurement |
|---|---|
| assembly, `T=200`, `N=220`, `M = 800x64` | `1338 s` |
| assembly, `T=400`, `N=440`, same `M` | `3532 s` — `2.64x` for `2x N` at fixed `M` |
| factorization (Rump, both sectors + checkpoint) at `N=440` | under `50 s` |
| scalar leg, `T=200`, 12800 nodes x 440 Bessel | `933 s` |
| scalar leg, `T=800`, 51200 nodes x 440 Bessel | `> 8700 s` (incomplete) — `> 2.4x` the node-count scaling |

For a full `T = 800` matrix (`N ~ 880`, `M ~ 1600x64`): Bessel count rises `4x`
and per-Bessel cost roughly `2.5x` from the larger argument, so the assembly is
**about `10x` the `T = 400` run, i.e. of order 10 hours** — not `4x`.
Accumulation (`N^2 M`) and factorization (`N^3`) each scale by `8x`, but both
were small at `T = 400` (factorization under 50 s), so they do not dominate.
All of this is an estimate anchored on measurements, not a guarantee.

> **Recommendation: CONSOLIDATE at `T = 400`.** The pilot has now landed and it
> argues against scaling, more sharply than the `lambda_Q` estimate did. The
> `T = 400` route is `~1.6%` from exhausted, so further gain needs a larger
> cutoff — but the whole `T = 800` route is capped at `+18.3%`
> (`lambda_R(800) <= 1.6681e-17`), and the realisable gain is strictly less,
> against a measured assembly cost of order **10 hours**. The bracket would at
> best narrow from `[1.41, 1.9014]e-17` to about `[1.67, 1.9014]e-17`, i.e. from
> width `4.91e-18` to `2.33e-18`, and only if the reoptimised `f_800` reached the
> ceiling exactly.
>
> If a `T = 800` attempt is nonetheless wanted, it should be **separately
> scoped** with that ceiling stated in advance, so the result is judged against
> `+18.3%` rather than against an open-ended expectation.

Neither the incomplete pilot nor any failed certificate here is an intrinsic
barrier. Nothing in this run bears on RH, zero locations, or larger support.

---

# 12. PHASE 7 (2026-11-12) — constants sharpened from the frozen blocks, conditioned coercivity, and the Connes-paper dictionary

No matrix was reassembled. `phase4/certs/block400_{even,odd}.json` were consumed
as committed; only the Rump shift and a rank-`K` rational bump were applied.
Full record, scripts, logs and artifacts:
`~/scratch/rh-connes-continuation-2026-11-12/notes/CONTINUATION_REPORT_2026-11-12.md`,
whose §10 records the corrections made after same-day external review.

## 12.1 The sharpened unconditioned constants

| quantity | committed | **2026-11-12** | factor | Rump `||I-G||_inf` |
|---|---|---|---|---|
| even-sector `a` at `T=400` | `1.42e-17` (`mu = 1.41e-17`) | **`1.43300048828e-17`** (`mu = 1.43300048827856e-17`) | `1.0163` on `mu` | `1.12012e-12` |
| odd-sector `a` at `T=400` | `2.4e-17` (`mu = 1.47e-17`) | **`1.3857e-14`** (`mu = 1.38569999999861e-14`) | **`577`** | `5.65580e-18` |

`a` is the Rump-certified shift on the finite block; `mu` is the whole-window
constant that `certify_final` then certifies, and it is **strictly below** `a`
by construction — the Schur test requires `(a-mu)(d-mu) > kappa^2`. Every `mu`
below is published as a short exact decimal with its margin; see §12.4.

New whole-window bracket at `L = 4/5`, the **even** sector still binding:

```
    1.43300048827856e-17  <=  lambda_Q(4/5)  <=  1.67312553407068186807726837736e-17
```

short quotable form `1.4330e-17 <= lambda_Q(4/5) <= 1.674e-17` (the upper
endpoint rounds **up**: `1.6731e-17` would not be a bound);
width `2.4013e-18`, down from `2.64e-18`. Parity separation margin
`1.38403e-14` (was `6.2687e-18`), ratio `mu_odd/U = 828`. The odd constant now
**exceeds** Zhu's Thm 6.2 value `8.206e-15` by a factor `1.689`; Zhu's
simplicity and evenness results are untouched and remain his.

**The `T = 400` even route is exhausted.** The certified floor
`1.43300048828e-17` and the independently measured ceiling
`R_400(f_400) = 1.4330724047355850e-17` (`phase4/certs/cross.json`) pin
`lambda_min(A_even, T=400)` to a relative width of `5.0e-5` — not the `~1.6%`
of §10.3, which measured the gap from the *committed* `mu` rather than from
`lambda_min(A)`. Confirmed a third way: inverse iteration on the midpoint block
gives `lambda_1(A_even) = 1.43307240474e-17`, matching `R_400(f_400)` to twelve
digits, and `lambda_1(A_odd) = 1.38574029573e-14`. §11.6's `T = 800` headroom is
correspondingly `+16.4%`, not `+18.3%`; **CONSOLIDATE at `T = 400`** stands,
strengthened.

The tests are not vacuous: the same routine **fails** just above each certified
shift (`1.43308715820e-17` even, `2.4576e-14` odd).

## 12.2 Certified conditioned coercivity, **by parity sector**

For any explicit vectors `v_i` and any `c >= 0`,
`A + c sum_i v_i v_i^T >= a I` implies `<A f, f> >= a ||f||^2` for every `f`
with `<v_i, f> = 0`, since the rank-`K` bump vanishes there. One Rump test on a
rank-`K` **exact-rational** perturbation of a frozen sector block therefore
certifies conditioned coercivity; by Courant-Fischer it also certifies
`lambda_{K+1}(A) >= a`. The `v_i` were chosen by inverse iteration on the
midpoint matrix — a step needing no rigour, since only the interval test is a
certificate.

| sector | conditions | certified `a` | `||I-G||_inf` | sector `mu` | also |
|---|---|---|---|---|---|
| **even** | 4 | `0.61755371093750` | `2.6216530e-32` | **`0.617553710936882`** | `lambda_5(A_even) >= 0.61755371` |
| **even** | 5 | `0.85349273681641` | `2.1380206e-25` | **`0.853492736815556`** | `lambda_6(A_even) >= 0.85349274` |
| **odd** | 3 | `0.10680580139160` | `2.1597066e-37` | **`0.106805801391493`** | `lambda_4(A_odd) >= 0.10680580` |
| **odd** | 4 | `0.78413009643555` | `9.2955279e-29` | **`0.784130096434765`** | `lambda_5(A_odd) >= 0.78413010` |

The complement constant (`d = beta - M g_N` even, `beta - M g_N - 2|p_perp|^2`
odd, the pole term carried) and the coupling
`kappa = 1.3083437712853882590e-29` are **unchanged**; only `a` changed.

**Combined over both sectors** (the form decouples by parity, and complex `f`
splits as `Q(u) + Q(v)` with each part splitting again):

| conditions | certified lower bound on `Q_zeta`, **both sectors** |
|---|---|
| 4 even + 3 odd = 7 | `Q_zeta(f) >= 0.106805801391493 ||f||^2` |
| 4 even + 4 odd = 8 | `Q_zeta(f) >= 0.617553710936882 ||f||^2` |

So the smallness of `lambda_Q(4/5)` is carried by **four even and four odd
explicit directions out of 440**. Conditioning on them raises the certified
constant by sixteen orders of magnitude. **Four conditions on the combined space
cannot do this** — min-max forbids it, since the odd sector carries three
minuscule eigenvalues of its own; the per-sector counts matter and are stated
above.

**Scope.** This proves that *some* explicit subspace of the stated dimension
works. It does **not** identify that subspace with the `E`-images of the
zero-integral prolate combinations (arXiv:2611.13129v1 Fact 5.1), which is
untested here. Removing directions changes the minimum over what remains; it
does not raise `lambda_Q(4/5)`, which stays in the §12.1 bracket. And it does
not escape the gate: `d` still needs `T > T_gate(L) = 2 pi e^{A_L}(1+o(1))`
with `A_L ~ 4 e^L`, so the route remains doubly exponentially capped in `L`.

## 12.3 The dictionary to arXiv:2611.13129v1 (Andrews-Connes-La Gatta)

`lambda = e^L` and `QW_lambda(g,g) = Q(f)` for `g(u) = f(log u)` — **proved**,
the two archimedean normalizations differing by exactly
`int_0^inf (1-e^{-s})/sinh s ds = 2 log 2`, and verified to ~30 digits by an
independent implementation of their eqs. (6)(7)(8) against
`phase5/certs/conventions.json` (`phi_0`, `phi_2`, `phi_0 - phi_2`), pole
moments included. Hence `eps_lambda(1) = lambda_Q^{even}(L)` and
`eps^odd_lambda(1) = lambda_Q^{odd}(L)` with `lambda = e^L`; our window
`L = 4/5` is `lambda^2 = 4.95303242439511480365`, just below their smallest
tabulated `lambda^2 = 5`.

Their eigenvalues are eigenvalues of compressions, so by min-max they bound the
limiting eigenvalues only from above — the authors state this (their §4.1
"Convergence in `N`", Fact 4.1) and treat their values as numerical evidence
with that limitation. What their method cannot do is **certify** the
infinite-dimensional lower bound. Combining §12.1 with the exact Theorem 2.1
rungs at `lambda = e^{4/5}` (Bouwkamp's method, App. A, validated against their
printed 30-digit `chi_0`, `chi_4` at `lambda^2 = 5`, and reproducing their
`-log delta_4 = 41.185`, `-log delta_6 = 33.833`, `-log delta_10 = 21.261`):

```
    delta_4 = 2.2391127080494187067e-18     delta_6  = 3.4175764349345855425e-15
    delta_8 = 2.4465397668381e-12           delta_10 = 9.3660403962553e-10

    6.3998587 <= eps_lambda(1)/delta_4(lambda) <= 7.4722703        (even)
                 eps^odd_lambda(1)/delta_6(lambda) >= 4.0546277    (odd)
```

against their `7.06` and `4.84` at `lambda^2 = 5`, both upper bounds. This is an
explicit certified enclosure at a stated window with a reproducible
certificate — **no priority is claimed**, and `delta_4`, `delta_6` are
high-precision and validated but **not yet interval enclosures**, which is the
one uncertified input and the gating item for every ratio quoted here.

Independent corroborations that fell out of the eigenvector step:

* the counts below `10^{-2}`: **exactly four** even (both bounds certified, see
  §12.5) and **at most three** odd — the paper's Table 3 at `lambda^2 = 5`
  records `4` and `3`;
* the ratios `lambda_j(A)/delta_{4j} = 6.400, 2.954, 2.734, 1.511` (even,
  `j = 1..4`) and `lambda_j(A)/delta_{4j+2} = 4.055, 2.590, 2.634` (odd) all sit
  in or just below the paper's ranges. **These are NUMERICAL**: `lambda_j(A)` are
  eigenvalues of a compression of `R_T`, and since `Q >= R_T` they bound
  `eps_lambda(j)` from neither side for `j >= 2`;
* their sharp test eq. (33) gives
  `log(delta_8/delta_4) - log(lambda_2/lambda_1) = 0.773044` against their
  `0.771` at `lambda^2 = 5` and their reported range `0.53-0.83` over
  `5 <= lambda^2 <= 300`;
* Fact 5.3 item 3 predicts that the odd sector, after `lambda^2 - 2 = 3`
  conditions, retains "one transitional eigenvalue `~ 0.08-0.10`" which the
  pole condition removes, "multiplying the minimum by a factor between `2.1`
  and `9` (`lambda^2 = 5`)". We certify that transitional eigenvalue at
  `>= 0.106805801391493` and the factor at `7.342`.

The cross-sector test shows the one discrepancy and it is accounted for:
`log(lambda_1(A_odd)/lambda_1(A_even)) = 6.8742` against exact
`log(delta_6/delta_4) = 7.3306`, shortfall `0.4564` versus the paper's `0.378`;
the excess `0.079` is exactly `log` of the ratio of the two sectors' `R_T`
deficits (`0.9065` even, `0.8378` odd).

Classification of the paper's contribution: Theorem 2.1 and Appendix A are
**usable** (exact rungs); Proposition 8.17's monotonicity half is
**unconditional and usable**; Theorem 8.14 + Lemma 8.16 + Proposition 8.18 are a
**conditional reduction** resting on (H2') — a hypothesis about smoothed
semilocal pair vectors — plus (R1), (R2) (the latter explicitly left open by the
authors) and an asymptotic-family hypothesis that no finite set of windows can
supply; §6-§7 (ground state, Mellin zeros) give **no transfer**, and their
Theorem 3.1/3.2 upper bound is weaker than what is already certified here.
**Nothing in the paper lower-bounds `eps_lambda` at any `lambda`.**

Theorem 8.14 is nevertheless the paper's real contribution to this project,
because it delivers the complement constant `d = (log sigma)/2` with no `T` and
no `A_L` — i.e. it is a gate-free route to the one estimate that makes this
repository's chain doubly exponentially expensive in `L`. That, and not its
numerics, is what would scale.

## 12.4 The Schur margin, published explicitly

`certify_final` requires the three strict signs `a > mu`, `d > mu`,
`(a-mu)(d-mu) > kappa^2`, so `mu` must sit strictly below `a` by at least the
coupling penalty. Its supremum satisfies `(a-mu*)(d-mu*) = kappa^2`, i.e.
`a - mu* = kappa^2/d` to leading order:

```
    a - mu*  =  1.41280e-58   at L = 4/5        (kappa^2/d = 1.4127964e-58)
             =  4.96559e-57   at lambda_-       (kappa^2/d = 4.9655891e-57)
```

reproducing `kappa^2/d` to six digits in both cases — a check on the Schur
formula itself. A `mu` published at that distance from `a` cannot be serialized:
its margin is `1e-41` of `a` relatively, far below any reasonable print width,
so a round-tripped value would be a ball straddling `a` and would **fail**
replay. Each `mu` above is therefore

```
    mu  =  a * (1 - 10^-12),  truncated DOWN to 15 significant decimal digits,
```

giving up `10^-12` of `a` relatively — about `10^30` times the true margin — so
it certifies with large room while remaining an exact short decimal that
replays bit-for-bit. All eight certificates, their margins, and the fact that
`certify_final` **fails** at `mu = a` in every case, are recorded in
`~/scratch/rh-connes-continuation-2026-11-12/certs/schur_margins.json` and
audited by `work/margin_audit.py`, which exits non-zero unless both directions
hold.

## 12.5 The even-sector count, certified both ways

The conditioned certificates bound the counts from **above** only:
`eps_lambda(5) >= 0.6175537109` gives at most four even eigenvalues below
`10^{-2}`, and `eps^odd_lambda(4) >= 0.1068058014` at most three odd. For the
lower side, `src/certify.py::Win(4,5).matrix(N)` is the Weil form on the first
`N` even cosine modes in an **orthonormal** basis, in closed form through
digamma / Lerch `Phi` / Hurwitz `zeta` as Arb balls — a **compression of `Q`**,
so by min-max `lambda_j(M) >= eps_lambda(j)` and its eigenvalues are *upper*
bounds. For an explicit 4-dimensional `W` of frozen rationals,

```
    eps_lambda(4)  <=  lambda_4(M)  <=  lambda_max(M|_W)  <  c ,
```

certified by the `4x4` test `c G - B > 0`, `G = W^T W`, `B = W^T M W`. At
`N = 160`, `prec = 2920`:

```
    eps_lambda(4)  <=  0.00311839891820  <  10^-2        CERTIFIED
```

Hence **exactly four even eigenvalues of `Q` below `10^{-2}` at `L = 4/5`**, both
bounds certified. `Win` implements the even basis only, so in the odd sector only
**at most three** is certified; the observed `lambda_j(A_odd)` are numerical.

The same compression gives certified upper bounds `eps_lambda(1..5) <=
1.6731255e-17, 8.4227418e-12, 6.4517875e-7, 0.0031183989, 0.62402164`; running
the even conditioning at `K = 1, 2, 3` would supply matching lower bounds and
bracket the first five even eigenvalues individually. Left to a later version.

---

# 13. PHASE 7b (2026-11-12) — a second window, and what it tests

Same machinery at a second window, to make the comparison with
arXiv:2611.13129v1's `lambda^2 = 5` row as close as this code allows.

## 13.1 The window

`L~ = 80471895621705018/100000000000000000`, strictly below `L* = (log 5)/2` by
`7.3e-18`, so `lambda_-^2 = 4.999999999999999927`. The point of choosing it
below `L*` is that `log n < 2 L~` keeps `n = 5` out of the prime comb, so
`A_L~ = 2.941973525223620455504` is **identical** to `A_{4/5}`,
`T_gate = 119.086906076` is unchanged, and `T = 400` remains admissible — the
run costs one assembly (`3678 s`, peak RSS `393.3 MiB`) rather than the ten
hours an `A_L` jump would force. (`A_L` jumps by `2 log 5/sqrt 5` for
`L > L*`, not at `L*`: at `L*` itself the strict inequality still excludes
`n = 5`, and in the time domain that term carries `C(2L) = 0`. So `L*` is itself
a cheap window — see §13.4.)

## 13.2 The certificates

Constants through `phase5/src/constants.py` at `prec = 2000`, all passing
`require_tight`; blocks `220x220`, max radius `1.36e-48`;
`beta_400 = 1.211613695058064549472`, `kappa = 7.756529961646111179258e-29`,
`d_even = d_odd = 1.211613695058064549472`.

| sector | certified `a` = `mu` (lower endpoint) | Rump `||I-G||_inf` |
|---|---|---|
| even | `8.273438682556152343750000e-18` | `5.05672988450e-11` |
| odd | `7.885009765625000000000000e-15` | `4.01037156749e-16` |

## 13.3 What it tests

**Monotonicity, rigorously.** `eps_lambda` is non-increasing in `lambda`
(arXiv:2611.13129v1 Prop. 8.17, unconditional), so with
`4.95303 < lambda_-^2`:

```
    mu_even(lambda_-) = 8.2734e-18  <=  eps(lambda_-)  <=  eps(e^{4/5})  <=  1.674e-17
```

which holds. This is the only falsification test the pair of runs supports: a
lower bound at the *larger* window exceeding our own certified upper bound at
the *smaller* one would contradict Prop. 8.17 together with
`phase5/certs/upper_highN.json`. It does not. **A lower bound at `lambda_-`
exceeding an external upper bound at exactly `lambda^2 = 5` would not be
contradictory** — the smaller window legitimately has the larger minimum — and
an earlier draft of this section wrongly asserted otherwise.

**Two-window scaling — a DIAGNOSTIC, not a certified scaling result.** The two
certified lower bounds are **certificate floors, not eigenvalues** (each sits at
roughly `90.6 %` of its `eps` in the even sector, and the two happen to sit at
similar fractions, which is an observation), and the rungs below are
**high-precision values, not interval enclosures**. The agreement is evidence
that the two runs and the dictionary are mutually consistent; nothing is
inferred from it:

| | `mu(lambda^2=5.000)/mu(lambda^2=4.953)` | exact rung ratio | rel. diff |
|---|---|---|---|
| even | `0.577351` | `delta_4(5)/delta_4(4.953) = 0.579867` | **`0.43 %`** |
| odd | `0.569027` | `delta_6(5)/delta_6(4.953) = 0.592859` | `4.0 %` |

The even agreement to `0.43 %` is consistent with the exponent structure
`-log delta_4 = 4 pi lambda^2 - 9 log lambda - l_4`, and uses none of the
paper's eigenvalues — but it compares floors with uncertified rungs and is a
diagnostic only.

**Matched-window ratios**, against the paper's own printed
`delta_4(sqrt 5) = 1.2983880087183625604e-18`,
`delta_6(sqrt 5) = 2.02613988982e-15`:

```
    eps_lambda(1)/delta_4      ~>  6.3720849      (paper, N=960: 7.06)
    eps^odd_lambda(1)/delta_6  ~>  3.8916413      (paper, N=960: 4.84)
```

**These are NUMERICAL COMPARISONS, not established inequalities.** The
denominators are evaluated at `lambda^2 = 5` while the numerators are certified
at `lambda_-^2 = 5 - 7.3e-17`, and they are high-precision values rather than
outward enclosures. Consistency only. The paper's values are **rounded printed numbers**, not
outward enclosures, so they bound nothing; and the residual window mismatch
(`7.3e-18` in `L`, about `8.4e-16` relative in `eps` by the derivative
estimate) is an **uncertified** transport.

## 13.4 What is still needed for a genuine head-to-head

1. An `arb`-parameterized window. `block.Symbol`, `baseline.{comb_terms, A_L,
   M_sup}` and `block.assemble` all take a rational `L_num/L_den`, and
   `L* = (log 5)/2` is irrational. Re-parameterizing is a change to
   sign-critical code and must be validated by re-running the full 176-check
   sequence at `L = 4/5` through the new signature and requiring bit-identical
   certificates before use. Once done, the matched run costs one assembly.
2. An outward enclosure for the paper's side, rather than a printed decimal.
3. A quadrature-free **upper** bound at `lambda_-` (the
   `phase5/src/upper_highN.py` route needs no gate). By monotonicity it would
   bound `eps(lambda)` for every `lambda >= lambda_-`, including
   `lambda^2 = 5`, rigorously — the one piece that would let the comparison
   bite from this side.

Three defects in the new driver code for this phase, all precision ownership or
search logic, all failures to certify rather than false certifications, are
recorded as `RETRACTIONS.md` R15 with the rule they imply.

---

# 14. PHASE 7c (2026-11-12) — a two-sided bracket at `lambda_-`, and what transports

The quadrature-free upper route (`phase5/src/upper_highN.py`'s two-stage
structure, with `src/certify.py::Win` and `L` parameterized) run at `L~`.
`N = 160` even cosine modes; exploratory `mpmath` stage locates the vector and
certifies nothing; the frozen exact rationals are enclosed in Arb inside
`owned(prec = 2920)`. **No quadrature anywhere**, so there is no R14-type
remainder to bound.

```
    [CERTIFIED]  Q(f_160)/||f_160||^2   = [9.47476951083049032849709e-18 +/- 8.10e-43]
    [CERTIFIED]  outward upper endpoint =  9.47476951083049032849709e-18
    safe short decimal (rounds UP)      =  9.475e-18
```

## 14.1 The bracket

```
    8.27343868254787e-18  <=  eps(lambda_-)  <=  9.475e-18 ,
        lambda_-^2 = 4.999999999999999927 ,  even sector
```

width `1.2013e-18`, ratio `1.1452`. The second window this repository brackets
from both sides.

| window | certified lower | certified upper | width |
|---|---|---|---|
| `lambda^2 = 4.95303242439511480365` (`L = 4/5`) | `1.43300048827856e-17` | `1.674e-17` | `2.40e-18` |
| `lambda^2 = 4.999999999999999927` (`L~`) | `8.27343868254787e-18` | `9.475e-18` | `1.20e-18` |

## 14.2 The paper's value lies inside it

arXiv:2611.13129v1 reports `eps_lambda(1)` at `lambda^2 = 5`, `N = 960`, as
`-log eps = 39.231`, i.e. `9.1663e-18` (equivalently `7.06 x delta_4`):

```
    8.2734e-18   <   9.1663e-18   <   9.4748e-18          (74.3 % of the bracket)
```

A **reconstructed finite-compression estimate** — the lowest eigenvalue of their
`N = 960` log-Fourier compression, recovered from a printed `-log eps` — falls
inside a two-sided certificate at a window differing by `7.3e-18` in `L`. Their
basis, normalization and codebase are all different from ours, which is what
makes the agreement informative; but this is **numerical consistency only**, not
an inequality between certified objects, and their figure is **not** an
independently known value of the infinite-dimensional minimum.

## 14.3 What transports upward, and what does not

By Prop. 8.17 (`eps_lambda` non-increasing in `lambda`, unconditional), the
**upper** bound transports to every larger window:

```
    eps_lambda  <=  9.475e-18     for every lambda >= lambda_- ,
                                  INCLUDING lambda^2 = 5 exactly.
```

A rigorous upper bound at the paper's own window, with no parameter transfer
and no reliance on their numerics; `3.37 %` weaker than their reconstructed
`9.1663e-18`, which is the right side to be on.

**The lower bound does not transport upward.** `mu(lambda_-) = 8.2734e-18` says
nothing about `eps` at `lambda^2 = 5` or beyond — monotonicity runs the other
way. A certified *lower* bound at exactly `lambda^2 = 5` still needs the
`arb`-window reparameterization of §13.4 item 1.

## 14.4 Ratio comparison, two-sided — NUMERICAL

With `delta_4(sqrt 5) = 1.2983880087183625604e-18`:

```
    eps_lambda(1)/delta_4(lambda)  in  [ 6.3720849 , 7.2973329 ]        (paper: 7.06)
```

two-sided now, and containing `7.06`. **A numerical comparison, not an
established enclosure**: `delta_4` is high-precision rather than
interval-enclosed, and it is evaluated at `lambda^2 = 5` while the bracket is
certified at `lambda_-^2 = 5 - 7.3e-17`. The `L = 4/5` pair
`[6.3998587, 7.4722703]` has the matching window but the same missing
enclosure. Interval-enclosing `delta_{4j}` remains the gating item for every
ratio in §§12-14.

A last diagnostic: the two certified **upper** bounds scale as
`U(lambda_-)/U(e^{4/5}) = 0.56629` against the rung ratio
`delta_4(5)/delta_4(4.953) = 0.57987`, a `2.3 %` agreement — same status as
§13.3's `0.43 %`, i.e. a diagnostic.

---

# 15. PHASE 7d (2026-11-12) — a rigorous comparison of the two windows' minima

With `eps_0 = eps(e^{4/5})` and `eps_- = eps(lambda_-)`, pairing the lower
endpoint at each window with the upper endpoint at the other gives

```
    mu_- / U_0   <=   eps_- / eps_0   <=   U_- / mu_0 .
```

This uses **only certified endpoints** — no prolate denominator, no parameter
transfer, no monotonicity transport, nothing from arXiv:2611.13129v1. Taking all
four inputs as the exact decimals that are packaged (published `mu`, which round
down; published `safe_decimal` upper bounds, which round up):

```
    mu_-/U_0 = 8.27343868254787e-18 / 1.674e-17  = 0.494231701466...
    U_-/mu_0 = 9.475e-18 / 1.43300048827856e-17  = 0.661200053838...
```

hence, outward-rounded,

```
    0.494  <=  eps(lambda_-) / eps(e^{4/5})  <=  0.662 ,
```

i.e. **the true spectral minimum falls by between `33.8 %` and `50.6 %`**
between `lambda^2 = 4.95303242439511480365` and
`lambda^2 = 4.999999999999999927`. Consuming the full-precision endpoints
instead tightens it to `[0.494490, 0.661184]`, same 3-decimal outward form.

This is a statement about the **true minima** of two infinite-dimensional forms,
with every input certified and every division taken in the safe direction. It is
stronger than the floor-vs-rung diagnostic of §13.3, and it is the first
quantitative statement here relating the window Weil form at two windows.

**It does not establish any scaling law.** The rung ratio
`delta_4(5)/delta_4(4.953) = 0.579867` lies inside the bracket, at `51.2 %`
across. That is a **numerical rung-ratio reference**, not a prediction of
Hypothesis 3.4: the hypothesis compares `eps_lambda(j)` with `delta_{4j}` up to
constants `c <= C` and a factor `log lambda`, and says nothing about the ratio of
two Weil minima equalling the ratio of their rungs — that needs the extra
assumption that the proportionality factor is the same at both windows.
Narrowing the bracket below `0.058` (each window's own bracket ratio below about
`1.05`, against `1.1676` at `L = 4/5` and `1.1452` at `lambda_-`, i.e. roughly a
threefold narrowing) would test **that stronger approximation**, not the
hypothesis.

## 15.1 Two qualifications, retained

* **The comparison with `9.1663e-18` is numerical consistency only.** That
  figure is a **reconstructed finite-compression estimate** — the lowest
  eigenvalue of the paper's `N = 960` log-Fourier compression, recovered from a
  printed `-log eps = 39.231` — lying inside our bracket at a nearby window. It
  is **not** an independently known value of the infinite-dimensional minimum,
  and it is never used as a bound.
* **At exactly `lambda^2 = 5` there is only the transported upper bound.**
  `eps_lambda <= 9.475e-18` for every `lambda >= lambda_-` (Prop. 8.17). The
  lower bound runs the other way: `eps_lambda >= 8.27343868254787e-18` for every
  `lambda <= lambda_-`, i.e. on the whole window range
  `0 < L <= 0.80471895621705018` — an interval, not a point — but it does
  **not** extend past that endpoint, and in particular not to `lambda^2 = 5`.

## 15.2 Release acceptance

`~/scratch/rh-connes-continuation-2026-11-12/work/acceptance.py` applies this
patch to a worktree at `7515396` and re-verifies **every** release claim from
the packaged artifact files rather than from values held in memory: the 6 patched
files and all 176 repository checks on the patched tree; the 8 Schur margins
re-certified from their own stored `a`, `d`, `kappa`, `mu` (with `mu = a`
required to fail); the 4 conditioned certificates re-certified by rebuilding the
rank-`K` bump from the stored frozen rational vectors; the 2 `lambda_-` blocks
from their packaged checkpoints; the `lambda_-` upper bound from its stored
exact rationals; the two-window ratio from the four packaged decimals; and the
even-sector count certificate. Stage **B2** separately re-runs the two
**improved** unconditioned `L = 4/5` block bounds
(`A - 1.43300048828e-17 I > 0`, `A - 1.3857e-14 I > 0`) from the frozen blocks,
which stage A does **not** cover: `phase4/verify_t400.py` tests the committed
`1.42e-17` / `2.4e-17` shifts, and `certify_final` assumes `A >= a I` rather than
checking it. It reverts the patch, checks the worktree clean, and exits non-zero
on any failure, so it is usable as a release gate.
**Result: PASS, 42 checks, 0 failed.**
