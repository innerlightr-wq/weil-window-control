# REPORT.md

Build dated 2026-11-09, branch `positive-weil-certificate-2026-11-09` off
`ad78c59`. All numbers below are reproduced by `verify.py`.

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
