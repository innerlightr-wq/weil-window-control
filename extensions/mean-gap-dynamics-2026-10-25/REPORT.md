# Mean-gap dynamics under the de Bruijn–Newman flow

A bounded investigation of whether the arithmetic/geometric/harmonic mean gaps of the
reflected-coordinate partition of a zero carry information about zero displacement beyond
`delta^2`, and whether any specified evolution forces the sign of their change.

**Extension of [`weil-window-control`](../../README.md). No DOI, not part of the paper, not
peer-reviewed.** Derivations and assumptions: [`DERIVATION.md`](DERIVATION.md). Every claim
with its class: [`CLAIMS.md`](CLAIMS.md). Reproduction: [`REPRODUCE.md`](REPRODUCE.md).

> **Read [`AUDIT.md`](AUDIT.md) after this file.** It also carries the one certified
> result this line of work produced — `Lambda_DH > 0.16`, a factor `2.3e27` on the
> published bound — which is a Davenport–Heilbronn fact and **not** RH progress.
> A later audit of how these constraints
> combine **corrected three claims below** — the proximity cutoff's sharpness and its
> `N`-competitor form, and the hypothesis of the energy identity — and applied an explicit
> stopping rule. Each correction is marked in place.

---

## 0. Verdict in plain language

**The mean gap is a re-expression of `2 delta^2` and nothing more.** The arithmetic mean is
identically `1/2` — that is the functional equation, not information — and the geometric and
harmonic means are algebraically the same quantity (`G^2 = A M_H`). So all three means, and
all three gaps between them, are functions of the one number `D = 2 delta^2`. One of the
three gaps, `G - M_H`, is *worse* than `D`: it is two-to-one, rising to a maximum of exactly
`1/8` at `D = 3/8` and falling back to zero at complete dominance, so it cannot tell a
shallow off-line zero from a deep one.

**There is exactly one place where the mean language is not a relabelling,** and it is an
identity, not a mechanism: in this repository's normalisation, **de Bruijn's threshold is
exactly the largest mean gap over the zero set**,

```
tau_c  =  Delta^2/2  =  sup_rho D(rho)  =  sup_rho ( A - M_H )(rho) .
```

For `zeta` that is `1/2`; for Davenport–Heilbronn it is `1.6025374835598227777878195`,
which reproduces the parent extension's `tau_DH` to all 22 committed digits.

**The specified evolution does force the sign — for one zero, and that zero is the wrong
one.** Along the heat flow, the exact chain rule gives `D'(tau) = Im(z) Im(H''/H')`, and the
canonical product turns this into a clean sum over the other zeros. Reading its signs:

* a zero's own symmetry orbit always pushes `D` **down**, by at least `1` (the conjugate
  partner contributes exactly `-1`);
* another zero can push `D` **up** only if it is strictly further off the critical line on
  the same side;
* so if a zero is at the *maximum* distance from the line, its `D` must decrease, at rate
  at least `1`.

That last statement integrates to exactly de Bruijn's strip bound and nothing better. And
for any zero that is **not** extremal the sign is not forced: it is **refuted**, by an exact
rational counterexample in a degree-4 polynomial heat flow, where the inner zero of
`(z^2+1)(z^2+4)` has `D'(0) = +1/3` — its mean gap *grows* — while the extremal zero has
`D'(0) = -19/3`.

Reading the same sum quantitatively gives the one new statement here: a competing zero can
lift `D_j'` above `-1` **only if it shares the ordinate to within the target's own
displacement**, `|gamma_j - gamma_k| < |delta_j|` — and, after optimising, only if
`|gamma_j - gamma_k| < delta_max/sqrt3`, which is sharp. (Both sharpened in
[`AUDIT.md`](AUDIT.md) §1.1; the first pass gave the looser `delta_max/sqrt2` and a wrong
`N`-competitor form.) Necessary, not sufficient — but cheap, and it both explains why no
Davenport–Heilbronn census quartet can be a counterexample (every competitor is outside by a
factor of **189 or more**) and shows that the minimal counterexample is *confluent in the
ordinate*, the same configuration the repository already identifies as the detector's
structural blind spot.

**Nothing here distinguishes `zeta` from Davenport–Heilbronn.** Every hypothesis used is one
DH satisfies, and the DH control instantiates every conclusion: all four of its certified
off-line quartets have `D'(0) < 0` at `tau = 0` — between `-1.012` and `-1.161`, interval-
certified — at the very window where the parent extension certifies its Weil form is
**negative** (`-0.679071346441480590986...`). So **`D` decreasing is not a positivity
mechanism**: the one object in this repository with certified off-line zeros has its mean
gaps shrinking and its Weil form certifiably indefinite at the same parameter.

**One genuinely monotone functional does exist, and it is weaker than what is already known.**
The *total* mean gap `Escr = sum_rho D(rho)` has an exact dissipation identity and is
non-increasing. It is a properly defined functional, its monotonicity is proved, and it is
rigorously the sum of the displacement diagnostic — the three things the brief demands
before such an object may be mentioned. Its unconditional corollary is an upper bound on
`Lambda` that is **provably never better** than de Bruijn's. So it is real and it is not
progress.

**And the whole route is forward-only.** Every statement here propagates reality toward
*increasing* `tau` and says nothing at `tau = 0`, which is where the arithmetic and RH live.
The exact budget identity `Escr(0) = 2 int_0^Lambda (dissipation) dtau` makes that explicit:
showing `Escr(0) = 0` — which *is* RH — needs either `Lambda <= 0`, which is RH, or vanishing
dissipation, which is the configuration to be excluded.

**Verdict: diagnostic re-expression — with three decisive obstructions, one exact dictionary
identity, and one new necessary condition.** The obstructions: (i) the mean triple collapses
to the single scalar `D`, and one of its gaps carries strictly less; (ii) individual
monotonicity of `D` is refuted exactly, and the one case where the sign *is* forced integrates
to de Bruijn's classical bound and no further; (iii) the Davenport–Heilbronn control satisfies
every hypothesis used and has its mean gaps shrinking while its Weil form is certifiably
negative, so `D` decreasing carries no positivity content. **No computational pattern and no
change of variables here is progress toward RH, and nothing here bears on RH at all.**

---

## 1. State of the repository, and the overlap this work has with it

### 1.1 A finding about the checkout, recorded first

The local working tree of `weil-window-control` was on branch **`projection-formula`** at
`16c9f0d`, which is **not an ancestor of `origin/main`** and is missing most of the work
this brief refers to. `origin/main` was at `dfdabab`. The authoritative state was read from
`origin/main` via `git archive` into a scratch directory; **no branch was switched, nothing
was stashed, reset, committed or pushed, and the two modified tracked files
(`README.md`, `extensions/resolution-inheritance/src/wform.py`) and the untracked
`extensions/density-stiffness-audit/` were left exactly as found.**

This directory is therefore **untracked** in a working tree that is behind the remote. Before
it is committed, the checkout should be brought to `origin/main` and this directory moved
across unchanged; the scripts have no dependency on any other path in the repository.

### 1.2 What was read, and treated as authoritative

| file (on `origin/main`) | what it fixes for this work |
|---|---|
| `README.md`, `REPORT.md` | tiering conventions, non-claims, the `delta^2`–`L^3` law's status |
| `notes/proofs.md` §§B.1–B.4 | the `delta^2` law, `c(L) = 8L^3(1/3+1/sqrt5)`, Theorem B and **Retraction 18** (Theorem 7 is a *relative* estimate) |
| `extensions/rh-partition-coordinate-audit-2026-10-07/` | **the principal overlap** — see §1.3 |
| `extensions/davenport-heilbronn-2026-10-23/` | the certified negative witness, the census, the four quartets, `sigma_0`, `kappa` |
| `extensions/source-flow-detector-audit-2026-10-24/` | **the heat-family normalisation**, `tau_DH`, the `Lambda_DH` bracket, and the obstruction instrument ("the filter") |
| `extensions/moderate-gap-spectral-bridge-2026-10-21/06_CONFLUENT_BLIND_SPOTS.md` | the confluent cluster bound, `NG-2`, `NG-2'`, `NG-3`; the off-line pair's ordinate gap is exactly `0` |
| `extensions/dirichlet-windows-2026-10-23/` | the real-pair parity split; the `gamma = 0` blindness is one sector |

The "projected-rank counterexample" of the brief is on the unmerged `projection-formula`
branch, `notes/projection_sensitivity.md` Remark D.3 and §D′ (the rank is `2`, not `1`, away
from `theta = 0`; `binom(R+2,3)` is the `theta = 0` member only). It is recorded as read; it
plays no role here, because this investigation touches no window form.

### 1.3 Overlap, stated before any new claim

`extensions/rh-partition-coordinate-audit-2026-10-07` already contains the partition
dictionary, with verdict **MISMATCH**, and **its coordinates already are these means**:

```
its altitude   h = sqrt(p q)   IS   the geometric mean G
its response   R = 1/(4 p q)   IS   1/(2 M_H)
its asymmetry  chi = 2|delta| ,  and  4 h^2 + chi^2 = 1   (re-checked here, exact over Q)
```

So **the arithmetic and harmonic means are not new coordinates**, and that audit's two
findings — *monotonicity is coordinate-free, so no monotonicity theorem can be obtained by
changing coordinates*, and *convexity is coordinate-dependent and therefore not a statement
about the Weil form* — apply verbatim to `A - M_H` and `A - G`.

Three things here are not in that audit:

1. the **non-monotone** gap `G - M_H` (its dictionary has no non-monotone member) and its
   exact two-to-one structure;
2. everything about the **evolution** — it audited coordinates on a fixed configuration and
   explicitly disclaimed dynamics;
3. the identity `tau_c = sup D`.

The `delta^2`–`L^3` law is also untouched: this investigation never evaluates a window form,
so neither `c(L)`, `K(L,gamma)`, Theorem 7's relative status, nor Retraction 18 is affected.

---

## 2. The diagnostic, fixed (Part 2)

`src/means.py`, exact over `Q`, **ALL PASS**. Full statements in
[`DERIVATION.md`](DERIVATION.md) §1. The content:

```
delta = beta - 1/2 ,  D = 2 delta^2 ;   p = beta ,  q = 1 - beta
A = 1/2 ,   G = sqrt(1/4 - delta^2) ,   M_H = 1/2 - 2 delta^2 ,   D = A - M_H
```

* **`A == 1/2` identically** — this is the functional equation, carrying no information.
* **`G^2 == A M_H`** (`AM . HM = GM^2`), so `M_H = 2 G^2`: `G` and `M_H` are the same
  quantity algebraically.
* Therefore the whole triple is a function of the **single scalar** `D`.
* `A - M_H = D`; `A - G = 1/2 - sqrt(1/4 - D/2)` is strictly increasing and convex in `D`,
  hence a monotone reparameterisation carrying exactly `D`'s information.
* `G - M_H = sqrt(1/4 - D/2) - 1/2 + D` is **not** monotone. It is strictly concave, vanishes
  at `D = 0` (**equality**, `p = q = 1/2`, the RH point) and again at `D = 1/2` (**complete
  dominance**, `beta in {0,1}`, `G = M_H = 0`), with maximum exactly `1/8` at exactly
  `D = 3/8` — where `G = 1/4` exactly, since `sqrt(1/4 - 3/16) = 1/4`. Both the critical
  point and the critical value are rational, so this is not a numerical optimum.
* Along **any** evolution `tau |-> D(tau)`:
  `(A-M_H)' = D'`, `(A-G)' = D'/(4G)`, and `(G-M_H)' = D'(1 - 1/(4G))`, the last **flipping
  sign at `D = 3/8`**.
* The arc-coordinate identity `M_H'' = 1 - 4 M_H` (with `delta = (1/2)sin theta`) is true and
  **empty**: `theta` is a chart on a fixed configuration, and *every* `c_0 + c_1 cos 2theta`
  satisfies `f'' = -4(f - c_0)`. Negative control run at three `theta` for three `(c_0,c_1)`.
  **No attractor, no restoring law, no dynamics.**
* Outside `0 <= beta <= 1` the partition reading fails (`pq < 0`, `G` imaginary, `M_H < 0`)
  while `D` is unaffected. The ranges matter: `sup D^zeta = 1/2` exactly, but
  `sup D^DH = 2(sigma_0 - 1/2)^2 = 1.6025374835598227777878195 > 1/2` (interval-certified),
  so **on the one object here with certified off-line zeros the partition interpretation is
  unavailable on part of the admissible range.**

### The one non-vacuous consequence

`tau_c = Delta^2/2` and `sup_rho D(rho) = Delta^2/2` are literally the same number, so

> **de Bruijn's threshold is exactly the supremum of the arithmetic–harmonic mean gap over
> the zero set.** "All zeros of `H_tau` are real for `tau >= Delta^2/2`" reads: *all zeros
> become real once `tau` exceeds the initial maximal mean gap.*

This is a dictionary identity. It is worth recording because it is the only place where the
mean language names something, and it is **not** a mechanism: the direction of the
implication is unchanged by renaming its threshold.

---

## 3. The evolution, derived independently (Part 3)

### 3.1 The dictionary, verified before calculating

`src/check_dictionary.py`, dps 90, **ALL PASS**: `H_0^zeta(z) = xi(1/2+iz/2)/8` to `9.5e-89`
relative; `d_tau H = -d_z^2 H` for **both** kernels to `1.3e-40` relative;
`H_0^DH(z) = Lambda_F(1/2+iz/2)/2` to `2.0e-90`; `Lambda_F(s) = Lambda_F(1-s)` to `7.1e-58`;
`kappa` from the **Gauss sum** against its closed form to `3.9e-62` and against the parent's
committed 40 digits to `1.2e-41`; and `z = 2 gamma - 2i delta`, so `delta = -Im z/2` and
`D = (Im z)^2/2`, exactly as the brief states.

> **Recorded because it would otherwise look like a result.** The PDE check at `z = 171.4`
> first **FAILED**, with a relative residual of `6.0e-10`. The identity is fine and the
> kernel is fine: at `|z| ~ 171` the integrand oscillates ~170 times per unit and `mp.quad`
> on a coarse node list returned ten digits, while `H_0` itself is of size `1e-30` because of
> the gamma factor. With two quadrature nodes per oscillation and a finite-difference step
> matched to the precision, the residual is `1.3e-40`. A check that fails for want of
> resolution must not be recorded as a failed identity.

### 3.2 The exact chain rule

Differentiating `H_tau(z(tau)) = 0` and using `d_tau H = -H''`:

```
z'(tau) = H''(z)/H'(z) ,          D'(tau) = Im(z) * Im( H''(z)/H'(z) ) .
```

Both reproduce the brief's formulas, signs and factors. They are **exact chain-rule
identities and nothing more**, valid on a simple branch and undefined where `H' = 0`.
Verified two independent ways at the DH orbit `O_1` — an analytic `s`-chain-rule route using
Hurwitz-zeta and digamma derivatives, and a quadrature route differentiating `cos(zu)` under
the integral sign, sharing no special-function code — agreeing to `1.07e-46`, with
precision-doubling agreement `2.8e-47`.

### 3.3 The canonical product, done in the variable where it is legitimate

`H_tau` is **even** (`Phi` real and even) and real on `R`. Its growth is
`exp(C|z| log|z|)` because `Phi` decays double-exponentially, so `H_tau` has order `1` — and
in `z` the genus is `1`, where Hadamard would leave an unknown factor `e^{cz}`. Passing to
`w = z^2`, `h_tau(w) := H_tau(sqrt w)` has order `1/2 < 1`, hence **genus 0**, hence

```
h_tau(w) = h_tau(0) prod_k (1 - w/w_k) ,     absolutely convergent, NO exponential factor,
H'/H (z) = 2z sum_k 1/(z^2 - a_k^2) ,        w_k = a_k^2 .
```

That is the normalisation the brief asks for, and it is why the squared variable is used.
At a simple zero this gives

```
H''/H' (a_j) = 2 sum_{k != j} 1/(a_j - z_k)
```

over the full symmetric zero list **in the paired order** — the unpaired sum is not
absolutely convergent and the pairing is not optional.

### 3.4 The sign identity, and what it does and does not control

> **Identity.** At a simple zero `z_j = x_j + i y_j` of `H_tau`,
> ```
>      D'(tau)  =  -2 y_j  sum_{k != j}  (y_j - y_k) / | z_j - z_k |^2 ,
> ```
> and if the zeros lie in a bounded strip — one of de Bruijn's own hypotheses, true for
> `zeta` and for DH — **this sum converges absolutely.**

Classifying it exactly as the brief asks:

| | |
|---|---|
| **an exact chain-rule identity** | `D' = Im(z) Im(H''/H')`, and the interaction form above |
| **a sufficient condition that merely restates `D' <= 0`** | "all terms non-positive", i.e. "no zero lies further off the line on the same side" — which for the extremal zero is automatic and otherwise is the thing in question |
| **a condition derivable independently of the arithmetic** | `\|y_j\|` maximal `=> D' <= -1`. It uses only evenness, decay and a bounded strip. **It is the differential form of de Bruijn's bound**, and it applies to Davenport–Heilbronn verbatim |
| **a numerical trend** | all four DH census quartets have `D'(0) < 0`, including the three for which the extremality hypothesis is **false** |

The signed structure, read off term by term:

* the conjugate partner `conj z_j` contributes **exactly `-1`**;
* the partner `-conj z_j` contributes **exactly `0`**;
* the mirror `-z_j` contributes `-y_j^2/|z_j|^2 <= 0`;
* every zero with `|y_k| <= |y_j|` on either side contributes `<= 0` — in particular every
  real zero;
* a zero contributes **positively only if it is strictly further off the line on the same
  side**.

Hence the own-orbit term is `<= -1` always, and **the extremal zero has `D' <= -1`.**
Integrating that gives `Delta_tau^2 <= max(Delta_0^2 - 2tau, 0)` — exactly the classical
threshold of §2. **No new proof of de Bruijn's theorem is claimed**: the envelope argument
through collisions and the attainment of the supremum are not carried out here.

### 3.5 Collisions

Three exact statements (`DERIVATION.md` §3.8): a real simple branch stays real, because
`H''/H'` is real there, so **a branch can leave the axis only through a collision**; the
collision locus is `H_tau = H'_tau = 0`, solved by Newton using `d_tau H = -H''` and
`d_tau H' = -H'''` so that `H'` is never divided by at the collision itself; and the
elementary symmetric functions of a colliding cluster extend analytically through it while
the individual `D_j` need not. The toy model makes the last point explicit and exhibits the
consequence: `Escr` is continuous through the collision but `Escr'` **jumps**, from `-12` to
`-8`, because two heights have merged and their pair term has fallen to zero.

### 3.5b A necessary ordinate-separation condition — the one new cheap statement

`src/separation_bound.py`, `DERIVATION.md` §3.5b. Bounding a single competitor's positive
contribution by `Delta_tau^2/(2R^2)` with `R = |Re z_j - Re z_k|` (using
`max_{0<a<b} 2a(b-a) = b^2/2`):

> **A single competing zero can lift `D_j'` above `-1` only if**
> ```
>    R = |Re z_j - Re z_k|  <  min( y_j , Delta_tau/sqrt3 ) ,
> ```
> **i.e. `|gamma_j - gamma_k| < min(|delta_j|, delta_max/sqrt3)`** — a competitor must share
> the ordinate to within the target's own displacement. With `N` competitors all at distance
> `>= R`, the crude form is `R < Delta_tau sqrt(N/2)`. Necessary, **not** sufficient.
>
> **Amended by [`AUDIT.md`](AUDIT.md) §1.1.** The first version gave `Delta_tau/sqrt2`
> (valid, 22.5 % loose; sharp is `1/sqrt3`, attained at `(2Delta/3, Delta, Delta/3)`), an
> `N`-form `Delta sqrt(N)/2` that was **wrong**, and missed the per-zero factor `R < y_j`
> entirely. No conclusion changes — the old constant was valid, so everything it excluded
> stays excluded — but the per-zero factor is much stronger and is now the binding one.

Three readings:

1. **The minimal counterexample is confluent in the ordinate.** In §4.1 the competitor sits
   at `R = 0` exactly. That is the configuration the repository's
   `06_CONFLUENT_BLIND_SPOTS.md` identifies as structural — an off-line zero's
   functional-equation partner is always at ordinate distance `0`. Here the partner supplies
   the `-1`, and a **second, deeper orbit at the same ordinate** overturns it.
2. **No DH census quartet can be a counterexample.** Under the per-zero form the deepest
   competitor of each of `O_2, O_3, O_4` is `O_1`, outside `|gamma_j-gamma_k| < |delta_j|`
   by factors **189, 1086, 406**. The four certified `D'(0) < 0` are now *explained*.
3. **For `zeta` it is not vacuous, and that is as far as it goes.** The threshold is
   `|gamma_j - gamma_k| < 1/(2 sqrt 2) = 0.3535533906`; the mean spacing reaches it at
   `gamma* = 2 pi exp(4 pi sqrt 2) = 328290294.5`, and Lehmer-type pairs are inside it far
   earlier (`6709/6710`: `|d gamma| = 0.0377`). So it **localises** where non-monotonicity
   could occur and rules nothing out. **No claim is made that it occurs.**

Like everything else here it uses no arithmetic, so it does not distinguish `zeta` from DH,
and it adds nothing at `tau = 0` beyond the strip width.

### 3.6 The one monotone functional, with its hypothesis and its price

Symmetrising the sign identity over pairs:

> **`Escr(tau) := sum_j D_j(tau) = (1/2) sum_j (Im z_j)^2`**, the total mean gap, satisfies
> ```
>     Escr'(tau) = -2 sum_{j<k} (y_j - y_k)^2 / |z_j - z_k|^2  <=  0 ,
> ```
> whenever the number `m(tau)` of non-real zeros is **finite**, with equality iff every
> `Im z_j = 0`.
>
> **Amended by [`AUDIT.md`](AUDIT.md) §1.3 — a real error.** The first version said
> "whenever `Escr < infinity`". That is not sufficient and in the infinite case is
> *incompatible*: each conjugate pair contributes exactly `1` to the dissipation, so
> `m = infinity` forces `Diss = +infinity` while `Escr` can be finite.

This is a functional that is explicitly defined, whose monotonicity is **proved**, and which
is rigorously the **sum of the zero-displacement diagnostic** — so it meets the brief's bar
for being mentioned at all. No thermodynamic reading is attached, and the word "entropy" is
not used for it.

Its price. Each conjugate pair contributes exactly `1` to the dissipation, so
`Escr' <= -m(tau) <= -2` while `Escr > 0`, giving `Lambda <= Escr(0)/2` — provided `Escr` is
absolutely continuous on `[0, Lambda]`, which is a stated **hypothesis**, not a consequence,
and is the one place the infinite case differs from the polynomial one. But every non-real
zero has a conjugate partner with the same `D`, so `Escr(0)/2 >= Delta_0^2/2`:

> **de Bruijn's classical bound is always at least as strong, and strictly stronger as soon
> as there are four or more non-real zeros.** On the toy: `Delta_0^2/2 = 2`,
> `Escr(0)/2 = 5/2`, truth `0.699`.

And the obstruction is now exactly locatable. The only **unconditional** lower bounds on the
dissipation are the two terms that symmetry supplies for free: the conjugate-partner term,
exactly `1` per non-real pair, and the mirror term `y_j^2/|z_j|^2`, which is `1.3e-5` at the
DH quartet and negligible at any large height. Anything more needs a lower bound on
`sum_{j<k}(y_j-y_k)^2/|z_j-z_k|^2` over the remaining pairs — i.e. a guarantee that off-line
zeros are not clustered at equal imaginary part near one another, which is a statement about
the configuration to be excluded. §4.2 measures how badly the obvious heuristic fails to
supply it.

---

## 4. Controls (Part 4)

### 4.1 Finite polynomial heat flow — the counterexample, exact over `Q`

`src/toy.py`, **ALL PASS**. `P_tau = exp(-tau d_z^2)P_0` is a finite sum, so
`d_tau P = -P''` holds identically in `Q[z,tau]` (exact residual `0`) and the whole of §3
holds with a finite product, no convergence question and no tail.

> **Counterexample.** `P_0(z) = (z^2+1)(z^2+4)`:
> ```
>   z = 2i  (extremal) :  D'(0) = -19/3  <= -1      -- as the sign lemma requires
>   z = i   (inner)    :  D'(0) = +1/3   >  0       -- the mean gap GROWS
> ```
> **Universal monotonicity of an individual `D` is refuted, not merely unproved.**

Three independent derivations agree exactly: the interaction identity over `Q(i)`; the closed
form `P_tau = z^4 + (5-12tau)z^2 + (4-10tau+12tau^2)` with
`D_in = (5 - 12tau - sqrt(96tau^2-80tau+9))/4`, whose derivative at `0` is rational because
the surd is `sqrt 9 = 3`; and root-finding at 14 nodes agreeing to `1e-30`. `D_in` increases
across the **whole** branch, `0.5 -> 0.827` — not a single-point artefact.

The threshold is exact: for `(z^2+a^2)(z^2+b^2)` with `0 < b < a`, the inner pair has
`D'(0) = 4b^2/(a^2-b^2) - 1`, so `D'(0) > 0` iff `a^2 < 5b^2`, and `= 0` exactly at
`a^2 = 5b^2`. Checked at five exact rational points including the boundary.

**Minimality:** a single symmetry orbit can never do this — for a lone quartet the own-orbit
term is the whole sum and `D' = -1 - y^2/|z|^2 < 0`. Two orbits, degree 4, is minimal.

**Negative control on the "entropy" idea.** The other obvious quadratic functional,
`sum_j Re(z_j^2) = sum_j (x_j^2 - y_j^2)`, increases at **exactly `2n(n-1)`** for every
polynomial heat flow of degree `n` (because `d_tau P` has no `z^{n-1}` term, so `e_1` is
conserved and `e_2' = -n(n-1)`). It is monotone, universally — and it is **not** a
displacement diagnostic: it grows even when every zero is real and `D == 0` throughout.
Monotone is cheap. *Monotone and rigorously tied to `D`* is what §3.6 supplies, and §3.6 is
what turns out to be weak.

### 4.2 Davenport–Heilbronn — certified off-line zeros at `tau = 0`

`src/dh_flow.py`, `src/dh_interaction.py`, `src/certify.py`. `D'(0)` comes from `H''/H'`, so
**the infinite zero sum is never truncated.**

| orbit | `delta` | `D` | certified `D'(0)` (Arb, 400 bits) | own-orbit term | cross-orbit rest |
|---|---|---|---|---|---|
| `O_1` | `0.3085171825` | `0.1903657037` | `-1.16044463430610822178313161420` | `-1.0000129598` | `-0.1604316745` |
| `O_2` | `0.1508300806` | `0.0454994264` | `-1.04435944844102968030569764300` | `-1.0000017455` | `-0.0443577029` |
| `O_3` | `0.0743560504` | `0.0110576445` | `-1.01207352431659836487600009985` | `-1.0000001995` | `-0.0120733248` |
| `O_4` | `0.2242576946` | `0.1005830272` | `-1.12778399346959232977056785846` | `-1.0000016107` | `-0.1277823828` |

All strictly negative, ball radii `<= 4.5e-32`, with one-sided upper bounds printed as plain
decimals rounded **up**. Each is certified *as an analytic function of the supplied `z_1`*
(`H'(z_1) != 0` is part of the certificate); that `z_1` encloses an actual zero rests on the
parent census, which is **T2**. That is the parent extension's own convention and it is not
upgraded here. **Respecting the distinction the brief insists on:** the fully certified
object in the parent work is the *negative Weil witness* `Q_F(g_*)/\|g_*\|^2`, which does not
depend on any zero location; the *quartet decomposition* and these `D'` values are certified
only as functions of numerical zero coordinates.

**The interaction identity checked with its remainder.** An independent scan of the critical
line (sign changes of the real `Lambda_F(1/2+it)`, bisected to `1e-35`) gives 169 on-line
ordinates below `T = 260`; summing the identity over them and the four quartets:

| cutoff `\|Im s\| <` | zeros used | partial cross-orbit sum | remainder vs exact | rel. to `\|D'\|` |
|---|---|---|---|---|
| 25 | 16 | `-0.00046489` | `-0.15997` | 13.8 % |
| 100 | 104 | `-0.14965664` | `-0.010775` | 0.93 % |
| 200 | 256 | `-0.15823179` | `-0.0021999` | 0.19 % |
| 260 | 350 | `-0.15877279` | `-0.0016589` | 0.14 % |

The remainder shrinks monotonically, as an absolutely convergent tail must. It is a
**measured residual, not a bound**, and nothing rests on it.

*An independent reproduction.* The scan finds **122** on-line ordinates below `T = 200`, and
`122 + 2x4 = 130` — the parent extension's contour count and its reconciliation, obtained
here from separate code.

**A wrong estimate of my own, and why it matters.** A mean-spacing density heuristic predicts
the cross-orbit contribution as `-1.3026848`; the truth is `-0.16043167` — an **overestimate by
a factor 8.12**, and the prediction is larger in magnitude than the exact conjugate-partner
term `-1`, so a valid version of it would roughly double the decay rate. The reason is
structural: the nearest on-line ordinates to `gamma_1 = 85.69934849` are `83.1088410847` and
`87.6474763326`, a gap of `4.5386` against a mean spacing of `1.4881` — **3.05 mean
spacings**. The off-line quartet sits in a large zero gap. So the obvious density heuristic
is **not** usable as the dissipation lower bound that §3.6 identifies as the missing
estimate. That is the single most useful negative fact this control produced.

**Extremality, and what the sign lemma actually covers here.** `O_1` has the largest
`|delta|` in the census, so *within the census* the lemma applies to it. But the census
covers `0 < Im s < 200` while the proved strip allows `|delta| <= sigma_0 - 1/2 = 0.8951`,
nearly three times `|delta_1|`. **Global extremality of `O_1` is not established and no
finite census can establish it.** For `O_2, O_3, O_4` the hypothesis is false already inside
the census, so their negative `D'` is a numerical fact, not a forced one.

### 4.3 `zeta` — no RH assumption, no planted family

`src/zeta_flow.py`; measured output in `results/part4_zeta.txt`.

* **(Z1)** From an initial guess perturbed off the axis by `1e-3`, Newton returns `|Im z|` to
  `~3e-57` and `|Im z'|` to `~5e-58`, against `|Re z'| = 0.5796` — so the reality is
  *measured*, not imposed by the arithmetic (the iterate is kept complex and re-perturbed at
  every step). A real branch is pinned to the axis: `D == 0` and `D' == 0` identically, and
  tracking a real `zeta` zero forward in `tau` can **never** produce an off-line zero.
  *"The zero stayed real" is evidence of nothing* — the brief's warning, instantiated.
* **(Z2)** A branch leaves the axis only through a collision; the first backward collisions
  of three adjacent low pairs are located by solving `H_tau = H'_tau = 0`, with **scale-free**
  residuals (`|H|/|H''|`, `|H'|/|H''|`) because the raw `|H|` is `~e^{-pi gamma/4}` and an
  absolute threshold on it is meaningless. *The first version of this gate used an absolute
  threshold and reported three genuine double zeros as UNRESOLVED.*
* **(Z3)** The emergent conjugate pair is continued below the collision at
  `tau_c = -1.2755439545316`, `z_c = 154.053285221206571`. The local quadratic model gives
  `Im z ~ sqrt(2(tau_c-tau))`, hence `D ~ (tau_c - tau)` and `D' -> -1`. Measured:

  | `tau_c - tau` | `Im z` | `D` | `D'` | `\|H\|/\|H''\|` |
  |---|---|---|---|---|
  | `0.0001` | `0.014142259864659` | `0.00010000175703977` | `-1.00003514107` | `1.9e-25` |
  | `0.001` | `0.044725288797559` | `0.0010001757290125` | `-1.00035148585` | `2.9e-25` |
  | `0.01` | `0.14154573803933` | `0.010017597978549` | `-1.00352238631` | `1.6e-25` |
  | `0.1` | `0.45118798544583` | `0.10178529910533` | `-1.03599366929` | `1.6e-25` |
  | `0.4` | `0.92738204282596` | `0.43001872667803` | `-1.15519264123` | `2.1e-25` |
  | `0.8` | `1.3629476299572` | `0.92881312100299` | `-1.34537668298` | `2.5e-25` |

  `D` tracks `tau_c - tau` to four digits at the smallest step, and `D' -> -1` as it must.
  Newton converges in 4–7 steps; the break tolerance is set **above** the measured
  quadrature noise floor (`~1e-24` in `|H|/|H''|`), because a tighter one is unreachable —
  the first version used `1e-32` and burned its whole iteration budget at the floor.
* **One further observation.** At `tau_c - tau = 0.8` the pair has `D = 0.9288 > 1/2`, i.e.
  `|delta| > 1/2`: the zero has left `0 <= beta <= 1` altogether, so the binary-partition
  reading of the means is **unavailable** there while `D` remains perfectly well defined.
  **The partition interpretation fails along the flow even for `zeta`** — §2's range remark,
  instantiated.

**A correction to the natural reading of (Z3)** (`src/zeta_lehmer.py`). It is tempting to say
the emergent pair is the only non-real pair, hence extremal, hence the sign lemma *forces*
`D' <= -1`. **That is wrong, and cheaply so.** Adjacent `zeta` gaps shrink like
`2 pi/log(gamma/2 pi)`, and the exact two-zero model puts a collision at `tau ~ -(gap in z)^2/8`,
so close pairs at large height leave the axis **far earlier**:

| pair | `gamma_n` | gap in `gamma` | gap in `z` | collision scale `tau_c ~` |
|---|---|---|---|---|
| 19/20 | `75.70469069908393` | `1.44014937` | `2.88029874` | `-1.0370151` |
| 13/14 | `59.34704400260235` | `1.484734522` | `2.969469044` | `-1.1022183` |
| 6708/6709 | `7004.043723499329` | `1.019142676` | `2.038285351` | `-0.5193259` |
| **6709/6710** (Lehmer) | `7005.062866174921` | `0.03769849773` | `0.07539699545` | **`-0.00071058837`** |

The Lehmer pair's scale is smaller by a factor of exactly **1459.3753**. So under the backward flow
zeros leave the axis in order of increasing gap, and the low-height pair of (Z3) is one of the
**last** to leave: at `tau ~ -1.28` it is certainly **not** extremal, and the sign lemma does
**not** apply to it. Any `D' <= -1` observed there is a numerical fact *consistent with* the
lemma, not an instance of it.

> **So the position is uniform across both controls: in neither is the one sign-forcing
> hypothesis verifiable.** On the DH side `O_1` is extremal only within a finite census and
> `O_2, O_3, O_4` are provably not extremal; on the `zeta` side the tracked pair provably is
> not. This is not a defect of the computation — extremality is a statement about *all* zeros,
> so no finite computation can establish it, and that is precisely why the forced-sign
> statement integrates only to the global strip bound. A finite scan establishes no global
> monotonicity, and none is claimed.

---

## 5. The one alternative evolution considered (Part 5)

The brief's stated preference is a family already in the repository whose dependence on prime
coefficients can be tracked. There is exactly one: the **`kappa`-family**

```
F_kappa(s) = [ (1 - i kappa) L(s,chi) + (1 + i kappa) L(s,chi-bar) ] / 2 ,
             coefficients  a = (1, kappa, -kappa, -1, 0)  periodic mod 5,
```

with `kappa` a genuine arithmetic parameter, `tan(arg eps_chi / 2)` from the Gauss sum.
`src/alt_evolution.py` tests the question that decides it.

| the brief's criterion | answer |
|---|---|
| what varies | `kappa`, Gauss-sum derived — a real arithmetic parameter |
| what is preserved | the Dirichlet series, its periodicity mod 5, its abscissa — **and nothing else** |
| does `zeta` occur at a specified parameter? | **NO.** `zeta` is in the family at no `kappa` (conductor 1 vs 5; `zeta` has a pole) |
| does it preserve useful arithmetic relations? | **NO.** `Lambda_{F_kappa}(s) = Lambda_{F_kappa}(1-s)` holds **only** at `kappa_DH`. The root-number defect `(1-i kappa)eps_chi - (1+i kappa)` is `8.1e-51` there and `1.9e-3` already at `kappa_DH + 10^-3`; at `kappa = 0.2` the functional equation fails by `31 %`. Off the DH point there is no completed function, no explicit formula, no spectral side, and no finite-window Weil form at all |
| zeros of a deformed function, an approximant, or a detector polynomial? | zeros of a genuine Dirichlet series — but one with no functional equation, so not the spectrum of any Weil form |
| what inference could return to `zeta`? | **NONE.** There is no parameter path from the family to `zeta`, so no limit argument exists |

**Three of five criteria fail outright. The branch is closed, and no second alternative is
pursued** — as the brief directs when no justified alternative is available. The repository's
own prior findings already cover the obvious remainder: the source-flow audit §7 records that
*no appropriate prime-side representation of the deformed family has been established* at
`tau != 0`, and the parent DH extension §15 (Phase H) records that splitting the arithmetic
side by the support of `b_n` does **not** separate "has off-line zeros" from "has no Euler
product" (`lambda_min` of the truncation is negative at 17 of 18 configurations, including
windows where the true form is positive). Neither is reopened.

---

## 6. Noncircularity (Part 6)

Applied to each claim that could look promising.

**(a) The mean-gap collapse (M1–M3).** Does it distinguish `zeta` from DH? **No** — it is
algebra about `beta`. Does a distinguishing arithmetic assumption enter? **None does.** Does
it assume an RH-equivalent? **No.** It is an obstruction, not a route, and it is reported as
one.

**(b) The chain rule and the sign identity (F1, F4).** No arithmetic enters. Every hypothesis
— evenness, reality, double-exponential decay, a bounded strip, a simple zero — **Davenport–
Heilbronn satisfies**, with a positive kernel into the bargain (parent audit §6b). This is
the parent audit's **filter**, applied to my own statements, and they do not pass it.

**(c) The forced sign at the extremal zero (F7).** Not RH-equivalent as a hypothesis —
"this zero is extremal" is weaker than RH. But it establishes **only eventual approach**:
integrated, it is de Bruijn's forward bound `Delta_tau^2 <= Delta_0^2 - 2tau`, which says
nothing at `tau = 0`. And `tau = 0` is where the arithmetic and RH live. *What additional
estimate would constrain the original parameter value?* A lower bound on the pair dissipation
not presupposing where the zeros are — and §4.2 shows the obvious density heuristic
overestimates it by `8.12x` in the one case that can be measured, precisely because off-line
zeros sit in zero gaps. **This is the same obstruction the parent source-flow audit §7
records, and this investigation does not move it.**

**(c') The separation bound (F8b), the one new statement.** Does it distinguish `zeta` from
DH? **No** — it is read off the same sign identity and uses only the strip half-width. Does a
distinguishing arithmetic assumption enter? **None.** Does it assume an RH-equivalent? **No**,
and notably it is *necessary*, not sufficient, so it cannot be mistaken for one. Does it only
establish eventual approach? It establishes nothing about approach at all; it localises where
a *failure* of approach could live. What would make it useful? A matching **sufficient**
condition, i.e. a lower bound on the positive contribution rather than an upper one — which
is the dissipation lower bound of (c) again, in another guise. So it survives its own
noncircularity test and is still only a filter.

**(d) The monotone functional (F9, F10).** An equivalent reformulation is not automatically
circular, so this deserves the careful version. `Escr` is properly defined, its monotonicity
is proved, and it is rigorously the sum of `D` — it is not a redescription of the target. Its
hypothesis `sum_rho D(rho) < infinity` is the **only** hypothesis found in this whole
investigation that DH is **not known** to satisfy: DH has infinitely many off-line zeros
(Davenport–Heilbronn 1936; `>> T` up to height `T`, Levinson 1972 — **bibliographic
statements, neither source read here**), and the four census displacements show no decay. So
it is a difference *in kind* from the hypotheses of (b), all of which DH provably satisfies.

**Amended by [`AUDIT.md`](AUDIT.md) §1.3, and the amendment strengthens this paragraph.**
The correct hypothesis is not `sum_rho D(rho) < infinity` but `m < infinity`, finitely many
non-real zeros. DH fails *that* on the strength of the **original 1936 infinitude theorem
alone** — no summability guess and no second source. So the filter **is** passed: this is
the one hypothesis in the whole inventory that provably excludes the control. And it still
does not matter, because **(F11) proves (F10)'s corollary is never better than de Bruijn's
bound.** A genuine separator attached to a useless conclusion is the sharpest way to state
why this branch stops; see `AUDIT.md` §3.

**(e) Is the route forward-only?** **Yes, provably.** Every statement here propagates reality
toward increasing `tau`. Given `Escr(0) < infinity` (so that `Escr` is defined) and
`Escr(Lambda) = 0` (forward preservation), the budget identity
`Escr(0) = 2 int_0^Lambda (dissipation) dtau` is exact, and it makes the circularity
explicit: showing `Escr(0) = 0` — which *is* RH — requires either
`Lambda <= 0`, which is RH, or vanishing dissipation, which is the configuration to be
excluded.

**(f) The decisive one — what the DH control actually shows.** At `tau = 0`, DH's window form
is certifiably **negative** (`-0.679071346441480590986377586093...`, parent extension, Arb
radius `1.4e-116`), **while** all four of its certified off-line quartets have `D'(0) < 0`.

> **`D` decreasing is not a positivity mechanism.** The one object in this repository with
> certified off-line zeros has its mean gaps shrinking and its Weil form certifiably
> indefinite at the same parameter. Any proposal that reads "the mean gap decreases" as
> evidence of positivity has a certified counterexample inside this repository.

---

## 7. The five questions, answered directly

**1. Does the mean gap reveal anything beyond squared zero displacement?**
**No.** `A == 1/2` is the functional equation and `G^2 = A M_H`, so the triple is a function
of the single scalar `D = 2 delta^2`. `A - M_H` *is* `D`; `A - G` is a strictly monotone
convex reparameterisation of it; and `G - M_H` reveals strictly **less**, being two-to-one
onto `(0, 1/8)` with its maximum exactly `1/8` at exactly `D = 3/8`. The one non-vacuous
consequence of the renaming is the identity `tau_c = sup_rho (A - M_H)(rho)`.

**2. Does the specified evolution force its sign of change in any nontrivial setting?**
**Yes, in exactly one, and it is the classical one.** At a zero of maximal distance from the
critical line, `D'(tau) <= -1`, forced, with the conjugate partner alone contributing exactly
`-1`. That integrates to de Bruijn's threshold `Delta^2/2` and to nothing better. For any
**non-extremal** zero the sign is **not** forced — it is refuted, by an exact rational
counterexample where an individual `D` grows across a whole branch. Individual monotonicity
must not be inferred from the global strip bound, and here it provably cannot be.

**3. Is any resulting constraint specific to zeta arithmetic?**
**No.** Every hypothesis used is one Davenport–Heilbronn satisfies; the DH control
instantiates every conclusion; and the derivations never touch an Euler product, a prime
coefficient, or any multiplicative relation. The one hypothesis DH is not *known* to satisfy
is `sum_rho D(rho) < infinity`, and its corollary is provably weaker than what is already
known, so even that buys nothing.

**4. What exactly fails in the Davenport–Heilbronn control?**
Nothing fails numerically — and that is the failure. All four certified quartets move
**toward** the line (`D'(0)` between `-1.012` and `-1.161`, interval-certified) at exactly
the parameter where DH's window form is certifiably negative. So the diagnostic behaves
*well* on the object where the conclusion one would want is *false*: `D` decreasing carries
no positivity content. Three further specifics: the extremality hypothesis — the only
hypothesis that forces the sign — holds for `O_1` only within a finite census and is **false
for the other three** (and on the `zeta` side it is provably false for the pair that can be
tracked, §4.3, so it is unverifiable in *both* controls); `sup D^DH = 1.6025 > 1/2`, so the
binary-partition interpretation is unavailable on part of the DH strip; and the obvious
density heuristic for the dissipation overestimates it by `8.12x` there, because the off-line
quartet sits in a `3.05x` zero gap. One thing the control does supply positively: the
separation bound of §3.5b *explains* why no census quartet can be a counterexample — every
competitor is outside the per-zero cutoff `|gamma_j - gamma_k| < |delta_j|` by a factor of
**189 to 1086** (`AUDIT.md` §1.1) — so the four negative `D'(0)` are no longer merely
observed.

**5. Is there one finite, worthwhile next lemma — or is this only a diagnostic
re-expression?**
**It is a diagnostic re-expression, with three obstructions, one identity and one necessary
condition worth keeping (§0).**
The candidate next lemma is already proved here: the dissipation identity
`Escr' = -2 sum_{j<k}(y_j-y_k)^2/|z_j-z_k|^2`, a genuine monotone functional rigorously tied
to `D`. I am **not** proposing it as a next step, because (F11) proves its only unconditional
corollary is weaker than de Bruijn's bound, and the estimate that would strengthen it — a
dissipation lower bound — is measurably unavailable from zero density (§4.2) and amounts to
assuming off-line zeros are not clustered, which is the configuration in question. What is
worth keeping is four lines: `tau_c = sup_rho (A - M_H)(rho)`; the sign rule *a zero's
displacement can grow only if another zero lies further off the line on the same side*; the
quantitative form of that rule — *only if it shares the ordinate to within `|delta_j|`, and
to within `delta_max/sqrt3`* (§3.5b, sharpened in [`AUDIT.md`](AUDIT.md) §1.1), which is the
one new statement here and is cheap to test against any configuration; and the refutation of
individual monotonicity. Each closes a door; none opens one.

---

## 8. Non-claims

1. **Nothing here proves or advances the Riemann hypothesis**, and nothing here is about
   where any zero of `zeta` is.
2. **No new proof of de Bruijn's theorem is claimed.** (F7) is the differential statement;
   the envelope argument through collisions is not carried out.
3. **No improvement on any bound for `Lambda`** — for `zeta` or for Davenport–Heilbronn. The
   bound (F10) is proved *weaker* than the classical one.
4. **No novelty is claimed** for `z' = H''/H'`, for the pair-dissipation identity, or for the
   Hadamard/genus argument; all are textbook, and the parent audit already records the first
   as such. What is assembled here is their application to the mean-gap diagnostic, the
   counterexample, and the `tau_c = sup D` identity.
5. **"Dynamics" means a specified parameter evolution.** No physical time and no mechanism
   generating `zeta` is claimed or implied. `tau` is not time.
6. **The word "entropy" is not used** for `Escr` or anything else, and no thermodynamic
   reading is attached.
7. **A monotone functional is not an RH mechanism.** `sum_j Re(z_j^2)` is universally
   monotone and completely uninformative; that is in the record as a negative control.
8. **A finite scan establishes no global statement**, and starting from a real zero that
   stays real is not evidence that every zero is real — (Z1) is recorded precisely to make
   that concrete.
9. **The certified `D'(0)` values are certified as functions of the supplied zero
   coordinates.** The parent's fully unconditional object is its negative Weil witness, which
   depends on no zero location; that distinction is preserved and not blurred.
10. **The repository is unchanged outside this directory.** No manuscript, note, data file,
    certificate or licence was touched; no branch was switched; nothing was committed or
    pushed.

## 9. Repository effect, and software

**Effect: none.** Every byte of this work is inside
`extensions/mean-gap-dynamics-2026-10-25/`. No tracked file anywhere in the repository was
modified — **including the root `README.md`**, which is left for the author to decide on,
since the working tree is behind `origin/main` (§1.1). The pre-existing modifications and
the untracked `extensions/density-stiffness-audit/` are untouched.

```
Python 3.12.2 (CPython)      mpmath 1.3.0       numpy 1.26.4 (not used on any reported path)
python-flint 0.9.0 (Arb)     -- installed into the ambient environment for the certification;
                                it is additive, matches the version the parent DH extension
                                declares, and changes no pinned behaviour.
```

No randomness is used anywhere. All substantive loops live in script files. Every reported
number is produced by a script in `src/` and recorded in `results/` and `data/`.
