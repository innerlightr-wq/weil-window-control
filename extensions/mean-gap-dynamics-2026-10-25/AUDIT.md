# AUDIT.md — how the available constraints combine

A single bounded audit, run after [`REPORT.md`](REPORT.md) and under its stopping rule.
Four tasks: resolve the outstanding claims; map the constraints; apply the
Davenport–Heilbronn control; propose at most one intermediate lemma.

Scripts: `src/audit1_claims.py`, `src/audit2_constraints.py`, `src/audit4_lemma.py`.
Transcripts: `results/audit1.txt`, `results/audit2.txt`, `results/audit4.txt`.
Data: `data/audit1.json`, `data/audit2.json`, `data/audit4.json`.

---

## 0. Verdict, and the stopping rule applied

**The stopping rule fires for the RH branch, and does not fire for one bounded
Davenport–Heilbronn deliverable.** Both halves, stated plainly:

**Stop the RH branch.** The audit maps fourteen constraints. **Eleven of them are
satisfied by Davenport–Heilbronn**, which simultaneously has four reconciled off-line
quartets and a certified *negative* finite-window Weil form. So no combination drawn from
those eleven can exclude off-line zeros — not as a matter of difficulty, but because a
counterexample satisfying the whole combination is already certified inside this
repository. Exactly **two** constraints in the inventory are primitive separators: the
Euler product, whose one available route the parent Phase H closed, and finiteness of the
off-line zero count, which this audit newly identifies as the correct hypothesis for the
energy identity and which is attached to a conclusion *proved weaker* than de Bruijn's.
Everything else the audit produced in the RH direction is a familiar heat-flow bound or a
relabelling. That is the stopping rule's first clause, met.

**One thing survived, and it has now been certified — and it is not about RH.** The audit does identify a new quantitative
interaction between two independently established constraints — this extension's pointwise
sign identity (an *equality*, so it bounds `|D'|` from above as well as below) and the
unconditional shell-counting bound that the parent extensions already use for `N(T)`. It
yields a **lower** bound on `Lambda`, which is the wrong direction for RH and vacuous for
`zeta`. But applied to Davenport–Heilbronn it gives, as a plausibility estimate backed by
measured ingredients,

```
Lambda_DH  >~  0.0137        against the published certified   Lambda_DH >= 6.8268e-29 ,
```

a factor of `2.0 x 10^26`, narrowing the published bracket from **28 orders of magnitude to
about 2**.

**And then the cheaper route was taken instead, and it certifies more.** Rather than prove
that trajectory lemma, §4b builds a **single fixed-parameter certificate**: `H_tau` has a
nonreal zero at `tau = 4/25` exactly, by an Arb winding-number computation on one contour.
So
```
              Lambda_DH  >  0.16     -- CERTIFIED, not estimated,
```
a factor of `2.3 x 10^27` on the published lower bound, leaving the bracket
`[0.16, 1.6025…]` — **a factor of 10**. This needs no trajectory estimate, no isolation
hypothesis, and nothing uniform over an interval, so the two qualifications against §4 do
not arise. The §4 lemma is now **not needed** for the numerical improvement and remains a
separate theoretical question.

**Three of my own claims needed correction, one of them substantively.** §1.

---

## 1. Task 1 — resolving the outstanding claims

### 1.1 The proximity cutoff: valid, not sharp, and its `N`-competitor form was wrong

`src/audit1_claims.py` §1a, exact over `Q`.

**What was claimed** (`REPORT.md` §3.5b): a single competitor can lift `D_j'` above `-1`
only if `R = |Re z_j - Re z_k| < Delta_tau/sqrt2`; with `N` competitors,
`R < Delta_tau sqrt(N)/2`.

**Finding 1 — the `N`-form is wrong.** Each competitor contributes at most
`Delta^2/(2R^2)`, so `N` of them need `N Delta^2/(2R^2) > 1`, i.e.
`R < Delta sqrt(N/2)` — **not** `Delta sqrt(N)/2`. The stated form gives `Delta/2` at
`N = 1`, contradicting the `Delta/sqrt2` the same paragraph states. Corrected.

**Finding 2 — the single-competitor constant is not sharp.** Maximising
`2(1-u)u/(s^2+u^2)` (`s = R/Delta`, `u = h/Delta`) exactly: the stationary point gives
`s^2 = u^2/(1-2u)`, and substituting into `2u(1-u) = s^2+u^2` gives `3u^2 - 4u + 1 = 0`,
hence `u = 1/3` and `s^2 = 1/3`.

> **Sharp cutoff.** `R < Delta_tau/sqrt3 = 0.57735026919 Delta_tau`, attained exactly at
> `(y_j, y_k, h) = (2Delta/3, Delta, Delta/3)`, where the contribution is **exactly 1**.
> The previously stated `Delta/sqrt2` is valid but 22.5 % too loose
> (`phi(1/sqrt2) = 0.732050807569 < 1`).

**Finding 3 — a sharper per-zero form, not previously noticed.** Do *not* optimise over
`y_j`. For fixed `y_j` the unconstrained maximum of `2 y_j h/(R^2+h^2)` over `h > 0` is at
`h = R`, with value `y_j/R`. So

> **`R < y_j` is necessary**, i.e. in the Dirichlet variable
> ```
>          | gamma_j - gamma_k |  <  | delta_j |  .
> ```
> If the optimum `h = R` is cut off by the strip (`h = Delta - y_j < R`), then
> `2 y_j h > R^2 + h^2 >= 2Rh` by AM–GM, so `y_j > R` again — **so it is necessary in all
> cases**, with no hypothesis beyond the strip. Verified on a `19 x 19 x 200` grid: zero
> violations. It is necessary and **not** sufficient (7 of 19 grid points have
> `R = 0.9 y_j` and still admit no competitor).

For a **single** competitor the audited condition is therefore the conjunction
`R < min(y_j, Delta_tau/sqrt3)`, of which the first factor is new and is the binding one for
a shallow zero.

> **Single versus several interactions — the distinction, stated explicitly.** Every bound
> above concerns **one** competitor's contribution exceeding `1`. Several smaller
> contributions can sum past `1` while no single one does, so the single-competitor
> conditions are **not necessary for `N > 1`**. With `N` competitors all at ordinate
> distance `>= R`, each contributes at most `min(y_j/R, Delta^2/(2R^2))`, so
> ```
>        R  <  min( N y_j ,  Delta_tau sqrt(N/2) )        [necessary, N competitors]
> ```
> Both factors relax — linearly and like `sqrt N` — so a configuration with many shallow
> competitors is **not** excluded by the `N = 1` form. **No conclusion in this extension
> rests on the `N = 1` case:** the toy counterexample has `N = 1`, and every DH census
> target is excluded even by the `N = 3` form (nearest competitor `R = 20.4` or `56.9`
> against thresholds `0.446` to `1.851`). Recorded because the `N = 1` form is the
> quotable one and the easiest to over-apply.

Re-applied:

| | previous form | per-zero form |
|---|---|---|
| toy `(z^2+1)(z^2+4)`, `y_j = 1`, `R = 0` | satisfied | satisfied |
| DH `O_2` vs `O_1` | outside by `45x` | outside by **`189x`** |
| DH `O_3` vs `O_1` | outside by `128x` | outside by **`1086x`** |
| DH `O_4` vs `O_1` | outside by `144x` | outside by **`406x`** |

### 1.2 The strip bound is not the transition threshold

`src/audit1_claims.py` §1b. Three quantities were in play and are now kept apart:

| | meaning | status |
|---|---|---|
| `Lambda` | the de Bruijn–Newman constant: the **actual** transition | not computable from the zero set at `tau = 0` |
| `tau_c = Delta_0^2/2 = sup_rho D(rho)` | de Bruijn's **sufficient** threshold; the identity of `DERIVATION.md` §1.7 | an **upper** bound for `Lambda` |
| `Escr(0)/2` | the (S2) bound | also an upper bound, **weaker** |

Always `Lambda <= tau_c <= Escr(0)/2`, and **both can be strict.** Exactly, on the toy:

```
Lambda = (10 + sqrt46)/24 = 0.6992637492968862 ,   tau_c = 2 ,   Escr(0)/2 = 5/2
```

so `tau_c/Lambda = 2.8601511`. For Davenport–Heilbronn, `tau_c^DH = 1.602537483559822…`
**is the upper endpoint of the parent extension's own certified bracket**
`[6.8268e-29, 1.6025…]`; the identity names that endpoint and does not narrow the bracket
at all.

*Audit finding:* `REPORT.md` never asserts `tau_c = Lambda`, and calling `tau_c` "de
Bruijn's threshold" is correct — but the distinction was only implicit. It is now explicit,
with the factor-2.86 toy gap as a standing reminder.

### 1.3 The energy identity's hypothesis — a real error

`src/audit1_claims.py` §1c.

**What was claimed:** (S2) holds "whenever `Escr < infinity`". **That is wrong, and not
merely incomplete.** Every conjugate pair `{z_j, conj z_j}` with `y_j != 0` contributes
**exactly 1** to the dissipation, so with `m(tau)` the number of non-real zeros

```
Diss(tau) >= m(tau)/2 ,        hence     m = infinity  ==>  Diss = + infinity .
```

But `Escr = sum_j D_j` can be **finite** with `m = infinity` — take `D_j = 2^{-j}`,
exhibited in the script. So in that regime the identity is not unproved but
**incompatible**: a finite left side and an infinite right side.

> **Corrected hypothesis.** (S2) and its corollary require
> ```
>          m(tau) < infinity        (FINITELY MANY non-real zeros) ,
> ```
> not `Escr(tau) < infinity`. Under `m < infinity` everything is automatic: `Escr` is a
> finite sum, the double sum converges absolutely (pairs of real zeros contribute `0`;
> each pair with a non-real member is bounded by `2Delta/|z_j-z_k|^2` and
> `sum_k |z_j-z_k|^{-2} < infinity` by order 1), and the rearrangement is legitimate.

**The correction strengthens the filter rather than weakening it.** `REPORT.md` §6(d)
recorded `sum_rho D(rho) < infinity` as a hypothesis DH is "not known to satisfy", resting
on Levinson's `>> T` count *and* on a guess about summability, and explicitly declined to
claim the filter was passed. The corrected hypothesis `m < infinity` is failed by DH on the
strength of the **original 1936 theorem alone** — Davenport and Heilbronn proved `F` has
infinitely many zeros off the critical line. No summability guess and no second source is
needed. (Bibliographic, not verified here, exactly as the parent extension records it.)

So `K13` becomes a genuine separator — and §3 shows it is attached to a useless conclusion.

---

## 2. Task 2 — the constraint map

Full table with sources in `results/audit2.txt`; `data/audit2.json` carries it as data.
The three means count as **one** diagnostic, as the brief requires.

| id | constraint | restricts | freedom remaining | DH |
|---|---|---|---|---|
| `K0` | the mean-gap diagnostic (**one** scalar `D = 2delta^2`) | nothing — `A = 1/2` *is* the functional equation, `G^2 = A M_H` | all of `delta` | SAT |
| `K1` | functional equation / reflection | zeros in orbits | `delta` free in `(-1/2,1/2)` | SAT |
| `K2` | bounded zero strip | `sup D <= 1/2` (`zeta`), `1.6025…` (DH) | the whole interior | SAT |
| `K3` | de Bruijn forward bound `Lambda <= sup D(0)` | the parameter where reality is guaranteed | everything at `tau = 0` | SAT |
| `K4` | extremal sign rule `D' <= -1` at the deepest zero | the deepest zero's motion | every non-extremal zero (**refuted** exactly) | SAT |
| `K5` | proximity cutoff `\|gamma_j-gamma_k\| < \|delta_j\|`, `R < Delta/sqrt3` | *which* zeros can have growing `D` | configurations with a nearby deeper zero | SAT |
| `K6` | energy dissipation; `Lambda <= Escr(0)/2` | the aggregate | — (weaker than `K3`) | FAIL, via `K13` |
| `K7` | `delta^2`–`L^3` law (relative; Retraction 18) | how much an off-line zero can **cost** the form (a *lower* bound) | the sign of `lambda*(L)` for `L > 0.8` | SAT |
| `K8` | certified `lambda*(0.8) >= 8.9e-18` | the window form is positive there | nothing about zero locations | SAT-N |
| `K9` | certified four-mode detection on `[9/100,49/100]`; certified **non**-detection on `[0,0.05986]` | whether one 4-mode window sees a quartet | `delta in (0.05986,0.09)`; all other windows | SAT |
| `K10` | moderate-gap budget: depth only via `kappa`; trace `delta`-blind; off-line pair at ordinate distance `0`; log-additive | the ARTZ budget's `delta`-sensitivity | — | SAT |
| `K11` | confluent cluster bound; `NG-2`, `NG-2'`, `NG-3` | detectability under competing ordinates | — | SAT |
| `K12` | **Euler product / multiplicativity of `b_n`** | the arithmetic side's structure | — | **FAIL** |
| `K13` | **finitely many off-line zeros (`m < infinity`)** | enables `K6` | — | **FAIL** |

### What adds information beyond the others

**The flow block `{K3,K4,K5,K6}` is not four facts — it is two.**

* `K4 ==> K3`: integrating `D' <= -1` at the deepest zero *is* the strip bound.
* `K6 ==> Lambda <= Escr(0)/2`, which this extension **proves** is never better than
  `K3`'s (`Escr(0)/2 >= Delta_0^2/2` always), and whose hypothesis is `K13`.
* `K5` is **not** implied by `K4`. `K4` is an envelope statement about one zero; `K5` is a
  per-zero statement about every zero, and says something where `K4` is silent.

> **The flow block reduces to `K4` (= de Bruijn) and `K5` (new).**

**The window block `{K7,K8,K9}` has no usable new internal combination.**

* `K7 + K8` give `lambda_min(Q_delta) >= 8.9e-18 - c(L)delta^2` for `L <= 0.8`: a *lower*
  bound. Excluding an off-line zero needs an *upper* one. Wrong direction.
* `K8 + K9` **does** combine, and it is **already in the repository.** The added-quartet
  and count-preserving forms differ by exactly `2F(gamma)^2 >= 0`, so
  `Q_rep,delta <= Q_add,delta` and the certificate covers the count-preserving convention —
  which answers the parent's own "the quartet is artificial" caveat. I derived this in the
  audit and then found the parent had already computed it: `cold_certificate.json` carries
  the field `Q_add_minus_Q_replacement = 3.410319615341804e-25`, and
  `2 F(14)^2 = 3.410319615341804e-25` reproduces it to 16 digits (ratio `1.0`). It is
  `2.35e-23` of the certified margin — negligible. Stated in
  `planted-witness-closure-2026-10-08/CERTIFICATE_THEOREM.md` and
  `planted-witness-consolidation-2026-10-08/VERIFICATION_REPORT.md`. **Not new.**
  It still does not reach RH: concluding anything about `zeta` needs `Q_rest >= 0`, which
  is Theorem 8's hypothesis and is exactly what is unavailable.

**The detectability block `{K10,K11}` touches the flow block in exactly one place.**
`K11`'s `NG-2'` (margin exactly `0` at multiplicity `>= 2`) and `K5` describe the *same*
configuration from two sides: an off-line zero whose competitor shares its ordinate. It is
a coincidence of configuration, not a new estimate, and it is why the minimal
counterexample of `REPORT.md` §4.1 is confluent in the ordinate.

**`K0` adds nothing to anything.** Its only non-vacuous consequence is the dictionary
identity `tau_c = sup D`, which names `K3`'s threshold without narrowing it (§1.2).

---

## 3. Task 3 — the Davenport–Heilbronn filter

```
DH SATISFIES (11 of 14):  K0, K1, K2, K3, K4, K5, K7, K8, K9, K10, K11
DH FAILS      (3 of 14):  K6, K12, K13     -- and K6 fails only BECAUSE its
                                              hypothesis is K13, so the PRIMITIVE
                                              separators number TWO: K12 and K13.
```

`K8` is marked `SAT-N`, not `SAT`: DH's own window form is numerically positive at small
`L` (`lambda_min = 4.03e-5` at `L = 1.0`, parent level-2 scan) with no certificate, so
small-window positivity is satisfied by DH as far as anyone has measured, and is certainly
not a separator.

And DH simultaneously has: a certified **negative** finite-window Weil form
(`-0.679071346441480590986377586093…`, Arb radius `1.4e-116`); four reconciled off-line
quartets below `T = 200`, the deepest at `delta = 0.3085`; and (this extension, certified)
all four with `D'(0) < 0`.

> **Therefore no combination drawn from `{K0,K1,K2,K3,K4,K5,K7,K8,K9,K10,K11}` can, by
> itself, exclude off-line zeros.** DH satisfies every one and has off-line zeros. This is
> the parent source-flow audit's *filter*, applied to the full inventory rather than to one
> theorem's hypotheses.

Any future argument must use `K12` or `K13` **essentially**. Status of each:

* **`K12` (Euler product).** The parent Phase H showed that splitting the arithmetic side
  by the support of `b_n` does **not** separate "has off-line zeros" from "has no Euler
  product": `lambda_min` of the truncation is negative at 17 of 18 configurations,
  including windows where the true form is positive. *The route is closed; the property is
  not.*
* **`K13` (finitely many off-line zeros).** Newly identified here as the correct hypothesis
  for `K6` (§1.3), and DH fails it on the 1936 theorem alone. But `K6`'s conclusion is
  **proved weaker** than `K3`'s. **`K13` is a genuine separator attached to a useless
  conclusion** — which is, precisely, why the audit stops the RH branch rather than
  continuing it.

---

## 4. Task 4 — the one proposed lemma

Not an RH lemma. It is proposed because the stopping rule's second clause is met — a new
quantitative interaction between two independently established constraints — and because
its target is a **published number** that the parent extension itself calls structurally
weak.

### The interaction

* **(A)** This extension's pointwise sign identity (T1, proved):
  `D'(tau) = -2 y_j sum_{k != j} (y_j - y_k)/|z_j - z_k|^2`. It is an **equality**, so it
  bounds `|D'|` from **above** as well as below. Every previous use took only the lower
  side.
* **(B)** The unconditional shell count `#{k : |z_k - w| <= 1} <= C_0 log(2+|w|)`, the exact
  analogue for `H_tau` of the `N(T+1) - N(T) = O(log T)` estimate the parent extensions
  already use for `N(T)`. It follows from Jensen's inequality and the order-1 growth bound
  `|H_tau(z)| <= exp(C|z| log|z|)` already derived in `DERIVATION.md` §3.1.

An **upper** bound on `|D'|` is a **lower** bound on how long a branch takes to reach the
critical line, hence a lower bound on `Lambda`.

### Statement

> **Proposed Lemma.** Let `H_tau` be as in the dictionary of `DERIVATION.md` §0, with all
> zeros in `|Im z| <= Delta_0` and the shell count (B) holding with constant `C_0`,
> uniformly for `tau in [0, tau_0]`. Let `z_1(tau)` be the branch of a simple non-real zero
> with `z_1(0) = 2gamma_1 - 2i delta_1`, and suppose that on `[0, tau_0]`
>
> * the branch stays **simple**,
> * `|Re z_1(tau) - Re z_1(0)| <= 1`, and
> * **(isolation)** no zero outside the orbit `{+-z_1, +-conj z_1}` lies within
>   `2|Im z_1(tau)|` of `z_1(tau)`.
>
> Then, writing `y = Im z_1`,
> ```
>   |D_1'(tau)|  <=  A + B/|y| ,
>     A = 1 + y^2/|z_1|^2 + 2 C_0 Delta_0 log(2+|z_1|) ,     B = C_0 Delta_0 ,
> ```
> and consequently `D_1(tau) > 0` — the branch has **not** reached the critical line — for
> every
> ```
>   tau  <  y_1^3 / ( 3 ( A y_1 + B ) ) .
> ```
> Hence `Lambda >= min( tau_0 , y_1^3/(3(A y_1 + B)) )`.

The proof steps are the sign identity (proved), the Jensen shell count (standard, **not**
derived here), and one integration of `|y'| <= A/|y| + B/|y|^2`.

### Plausibility, every ingredient measured

`src/audit4_lemma.py`, on the Davenport–Heilbronn quartet `O_1` and the 354-point zero list
computed by `src/dh_interaction.py`.

| ingredient | measured | status |
|---|---|---|
| isolation at `tau = 0` | nearest non-orbit zero at `3.94481176273`, conjugate partner at `2\|y_1\| = 1.23406872983` | **holds, margin `3.19659x`** |
| shell constant `C_0` | `0.38792761` (max over four sample ordinates) | measured on a finite list, order of magnitude only |
| `S' = sum'_k \|z_1-z_k\|^{-2}` | measured `0.20851027` vs hypothesis bound `1.8753802` | slack `8.99419x` |
| `\|D'(0)\|` | certified `1.1604446343061082218` vs bound `9.286639392` | slack `8.00266x` |
| `A`, `B` | `8.161100685`, `0.6944960613` | from the hypothesis, not the measurement |
| **`tau*`** | **`>= 0.01366592`** | plausibility estimate |

```
published certified   Lambda_DH >= 6.8268151773112564447e-29
this estimate         Lambda_DH >~ 0.01366592              -- factor 2.0018e+26
published upper       Lambda_DH <= 1.602537483559822777787819
bracket               28 orders of magnitude  ->  about 2.07
```

Robust to the one measured constant: inflating `C_0` by `10^6` still leaves a factor
`2.2e+20` over the published bound (the estimate scales like `1/C_0`).

### What the lemma is not

* **Not about RH.** A *lower* bound on `Lambda` is evidence only in the direction of RH
  being false, and is vacuous for `zeta`, where no off-line zero is known. It cannot be run
  on `zeta` at all.
* **Not arithmetic.** It uses no Euler product and no prime coefficient, so by §3 it does
  not distinguish `zeta` from DH — and it does not need to, because its target is DH's own
  constant.
* **Not proved.** The Jensen shell count is standard but not derived here; the isolation
  hypothesis is checked only at `tau = 0`; and the branch staying simple with bounded real
  part is **not** checked for `tau > 0`. The constant `C_0` is measured on a finite list; a
  theorem needs the Jensen constant, which will be larger — which is what the robustness
  sweep is for.
* **Not a bracket improvement until written.** The estimate is `>~`, not `>=`.

### Noncircularity

Does it distinguish `zeta` from DH? **No**, and it does not claim to. Does a distinguishing
arithmetic assumption enter? **None.** Does it assume real zeros, Weil positivity, or any
RH-equivalent? **No** — it assumes the opposite, that an off-line zero exists, and is
vacuous otherwise. Does it only establish eventual approach? **It establishes the reverse:
a floor on how slowly the approach must happen.** What additional estimate would constrain
the original parameter value? For the DH bracket's *upper* end, nothing here helps; that is
`K2` and `sigma_0`.

---

## 4b. CERTIFIED: the fixed-parameter route, carried out

**Recommendation taken.** Rather than prove the trajectory lemma of §4, a **single
fixed-parameter certificate** was built — and it certifies far more than the lemma would
have estimated. `src/dh_fixed_tau.py` (navigation), `src/certify_fixed_tau.py`
(certificate), `results/certificate_fixed_tau.txt`, `data/certificate_fixed_tau.json`.

### What is certified

> **Theorem (computer-assisted, Arb at 700 bits).** The Davenport–Heilbronn heat family
> `H_tau` has a **nonreal** zero at `tau = 4/25 = 0.16` exactly. Hence
> ```
>                     Lambda_DH  >  0.16 .
> ```

| | |
|---|---|
| parameter | `tau = 4/25`, exact rational |
| disc | centre `171.2502 - 0.1823 i` (exact dyadic), radius `1/16` exact |
| `\|Im centre\| - r` | `0.11980000000000` — the closed disc is **disjoint from the real axis** |
| `min \|H_tau\|` on `dD` | `>= 5.756765507e-31` (plain decimal, rounded **down**) |
| worst relative enclosure radius | `0.0121658` over 512 arcs |
| winding number | **exactly 1** (net half-plane advance `4`, zero non-adjacent steps) |
| negative control | the same machinery on a disc shifted `0.6` in `Re z` returns winding **`0`** |

### Why this avoids every weakness of §4 and of the parent's own bound

* **No trajectory estimate, no isolation hypothesis, nothing uniform over an interval.**
  One parameter value, one contour, one integer. The two qualifications raised against §4 —
  single-versus-several interactions, and uniformity of the shell count and isolation over
  the deformation interval — **do not arise at all**.
* **No perturbation comparison with `H_0`.** The parent extension's `Lambda_DH >= 6.83e-29`
  came from `|H_tau - H_0| <= |tau| C` with `C` bounding `int x^2|Phi_F|e^{tau_0 x^2+Mx}dx`,
  which replaces `cos(zx)` by `e^{Mx}` and so **discards every bit of oscillatory
  cancellation** — the reason that bound is `10^{-29}` and not `10^{-1}`. Here `H_tau` is
  evaluated *directly*, from
  ```
    H_tau(z) = (1/2) sum_{m>=0} (tau/4)^m Lambda_F^{(2m)}(s)/m! ,     s = 1/2 + i z/2 ,
  ```
  (because `d_tau H = -d_z^2 H` and `d_z^2 = -(1/4) d_s^2`), so all the oscillation stays
  inside `Lambda_F`, which Arb encloses rigorously.
* **Two rigorous tails, both majorised by non-oscillatory positive integrals**: the
  `tau`-series truncation `T_M <= 9.9540082e-47` at `M = 20`, and the arc expansion
  `A_K <= 6.2959781e-47` at `K = 10`. Both are computed by rigorous Arb integration on
  `[0,3]` plus an explicit decay bound beyond.
* **The special functions are evaluated at exact POINTS, never on wide balls.** Feeding a
  ball of radius `3.1e-4` into the order-43 Hurwitz-zeta series inflates the enclosure of
  `H_tau` to `+/- 1.9e-28`, swamping its true size `3.8e-31`. Expanding about exact arc
  midpoints and bounding the arc tail separately does not.
* **The winding number is computed by half-planes, not by summing arguments.** Summing 512
  increments each enclosed to a relative `1.2e-2` accumulates a radius of order
  `2 sum(rel)`, which is *scale-invariant* in the arc count — the first attempt returned
  `+/- 3.69`, enclosing many integers. Placing each enclosure in one of the four open
  half-planes gives an **exact** integer with no accumulation.

### The representation is cross-checked independently at `tau != 0`

This is the one place where a mis-derived constant in (*) would not have been caught: the
internal `M`-consistency of the series holds for *any* constant, the stage-0 dictionary
audit checks `d_tau H = -d_z^2 H` through the **quadrature** route only and `H_0 = Lambda_F/2`
at `tau = 0` only, and the winding gates are self-consistent whatever `H` is.
`src/check_taylor_vs_integral.py` closes it against the defining integral
`int_0^inf e^{tau u^2} Phi_F(u) cos(zu) du`, by oscillation-resolving mpmath quadrature
that shares no code with the series and calls nothing beyond `exp` and `cos`:

| `z` | `tau` | `\|H\|` | series vs integral, **relative** | quadrature dps 90 → 140, **relative** |
|---|---|---|---|---|
| `171.2502 - 0.1823i` (the disc centre) | `4/25` | `5.6058091e-34` | `5.926e-31` | `1.593e-63` |
| `171.2502 - 0.1823i` | `1/10` | `3.7851744e-30` | `4.525e-39` | `2.081e-67` |
| `171.3100 - 0.2500i` | `4/25` | `1.1821328e-30` | `3.913e-34` | `7.985e-67` |
| `20.0 - 0.4000i` | `4/25` | `0.0021313049` | `1.109e-44` | `3.16e-92` |
| `5.0 + 0.3000i` | `1/5` | `0.37388813` | `9.632e-45` | `7.817e-92` |

Both deviation columns are **relative** (`|series - integral|/|integral|`, and likewise for
the dps 90 → 140 comparison), as `src/check_taylor_vs_integral.py:93-94` computes them. On
the headline row that is an **absolute** agreement of `3.3e-64` on a value of size
`5.6e-34`, i.e. the quadrature resolves the ~34 digits of cancellation and then another 30.

**And the check discriminates.** The same series with `(tau/2)^m` in place of `(tau/4)^m`
disagrees with the integral by a relative `18026` — so this is a real test of the constant,
not a tautology.

*How the gap was found:* a background job launched early in the session and killed before
it finished later reported completion, having reached only the `tau = 0` case. `tau = 0` is
the trivial one (the series collapses to its `m = 0` term, already audited at stage 0), and
noticing that the interesting case had never run is what prompted this check.

### The two inputs to the tail bounds are gated, not quoted

* `|Phi_F(u)| <= 8 e^{3u} e^{-pi e^{4u}/5}` on `u >= 0`: since
  `Phi_F = 4e^{3x} sum n a_n e^{-pi n^2 T/5}` with `T = e^{4x} >= 1` and `|a_n| <= 1`, one has
  `|Phi_F| <= 4 e^{3x} e^{-pi T/5} S(T)` with `S(T) = sum_{n>=1} n e^{-pi(n^2-1)T/5}`
  **decreasing** in `T`; `S(1) = 1.3236800759485`, so `4 S(1) = 5.294720304 <= 8`. The
  parent's constant is valid with `1.51x` of slack.
* `int_3^inf f <= f(3)` for every integrand used: `(log f)'(3) <= -409024.8675 <= -1` at the
  largest power `p = 52`, and it only decreases because `(4pi/5)e^{12} = 409047.4056`
  dominates.

### Why `net/4` is the winding number *exactly*

If arc `i`'s enclosure lies in the open half-plane of direction `alpha_i = h_i pi/2`, a
continuous branch has `arg H = alpha_i + e_i` with `|e_i| < pi/2` there. The `N` balls
**abut** (ball `i` is centred at `theta_i = 2 pi i/N` with half-width `Delta theta/2`, so they
tile the circle), and at a shared endpoint the same `H` value gives
`alpha_{i+1} - alpha_i = e_i - e_{i+1} in (-pi, pi)` — consistent with a step of `0` or `+-1`,
and the step is uniquely forced when the current index is not in the next arc's set.
Telescoping once round with `H(2pi) = H(0)` gives `net * pi/2 = 2 pi W`, so

```
                        W = net/4      exactly, not up to an error.
```

### Two of my own mistakes on the way, both logged in the script

1. **`K = 4` was not enough.** `A_4 ~ 2e-23` swamped `H_tau ~ 1e-30` by seven orders and
   every arc enclosure contained `0` — the same cancellation-loss disease as the parent's
   bound. `A_K ~ (rho u*)^{K+1}/(K+1)!` with `u* ~ 0.6`, so `K = 10` costs four more series
   coefficients and buys 29 orders. Raising the arc count instead would have needed ~61000
   arcs.
2. **Summing argument increments does not work**, for the scale-invariance reason above.

### The resulting bracket

```
previous :  6.8268151773112564447e-29  <=  Lambda_DH  <=  1.602537483559822777787819
                                                     span: 28 orders of magnitude
now      :                       0.16  <   Lambda_DH  <=  1.602537483559822777787819
                                                     span: a factor of 10.02
```

An improvement of `2.3 x 10^27` on the lower end, **certified, not estimated** — and it
overshoots the §4 lemma's own plausibility estimate (`>~ 0.0137`) by a factor of 12. The
navigation run shows why: the branch stays nonreal up to `tau ~ 0.161` and only then
reaches the axis, so the trajectory bound was conservative by an order of magnitude.

### Scope, stated precisely

* **Nothing here bears on RH.** For `zeta`, RH requires `Lambda_zeta <= 0` and
  Rodgers–Tao give `Lambda_zeta >= 0`; a positive lower bound derived from an **existing**
  off-line zero of a **different** function says nothing about `zeta`, where no off-line
  zero is known and the construction cannot be run at all.
* `0.16` is **not** `Lambda_DH`. It is a certified strict lower bound. The bracket is a
  factor of 10 wide and this does not locate the transition.
* The certificate is exact-rational/Arb ball arithmetic run by CPython and python-flint. It
  is **not** proof-assistant verification and has **not** been peer-reviewed.
* The disc centre came from navigation, but **nothing depends on the centre being near a
  zero**: the winding number is computed from `H_tau` on the contour alone, and the
  negative control confirms the machinery returns `0` where there is no zero.
* The §4 trajectory lemma is **not** proved and is now **not needed** for the numerical
  improvement. It remains a separate theoretical question, and §4's caveats stand.

## 5. Corrections made to the earlier documents

Append-only, following the parent extensions' convention. Each was found by this audit, not
by review.

1. **The `N`-competitor cutoff `R < Delta sqrt(N)/2` was WRONG** (`REPORT.md` §3.5b,
   `DERIVATION.md` §3.5b, `CLAIMS.md` F8b). Correct crude form `R < Delta sqrt(N/2)`;
   the stated version contradicted its own `N = 1` case. §1.1.
2. **The single-competitor constant `Delta/sqrt2` was not sharp.** Sharp value
   `Delta/sqrt3`, attained exactly at `(2Delta/3, Delta, Delta/3)`. The old constant is
   valid, so no conclusion drawn from it was wrong — only loose. §1.1.
3. **`Escr < infinity` was the wrong hypothesis for (S2)** (`DERIVATION.md` §3.6,
   `CLAIMS.md` F9/F10, `REPORT.md` §3.6/§6(d)). It is not merely insufficient: with
   `m = infinity` the identity is incompatible. Correct hypothesis `m < infinity`. §1.3.
4. **A new per-zero cutoff `|gamma_j-gamma_k| < |delta_j|` was available and missed.** Not
   an error, but the previous documents reported the weaker form only. §1.1.
5. **One combination I believed new is pre-existing** (`Q_rep <= Q_add`). Recorded in §2,
   with the parent's committed field reproduced to 16 digits.
6. **A comparison in `src/audit4_lemma.py` §3 failed on its first run** because it compared
   `S_1` (orbit terms included, `0.865`) against a bound for `S_1'` (orbit terms excluded,
   `0.532`). The bound was never wrong; the comparison was. Logged in the script itself.

---

## 6. What to do next

**Recommended: close the RH branch here.** §3 is the defensible stopping point the brief
asked for: an eleven-of-fourteen DH-satisfied inventory, with the two separators' status
each independently documented. Nothing further in the mean-gap or heat-flow direction is
warranted, and this audit does not propose any.

**Done, not optional:** §4b **is** the standalone sharpening, and it is certified rather
than proposed. What remains, if anything, is cosmetic: `tau` could be pushed from `0.16`
toward the observed branch limit `~0.161`, which changes nothing. The bracket is now a
factor of 10 wide and narrowing it further means locating `Lambda_DH`, a different and much
harder problem.

**The §4 trajectory lemma is no longer needed** and is not proposed. It would still be a
theoretical contribution — a general rate bound rather than a single parameter value — but
it is now strictly weaker than what §4b certifies on this object, and its two open
qualifications (several interactions; uniformity of the shell count and isolation over an
interval) are real work for no numerical gain.

**Label discipline:** §4b is a Davenport–Heilbronn result. It is not RH work, and the
parent extensions' own non-claims apply to it verbatim.
