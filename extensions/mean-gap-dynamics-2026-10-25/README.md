# Mean-gap dynamics under the de Bruijn–Newman flow

Extension of [`weil-window-control`](../../README.md). **Carries no DOI, is not part of the
paper, and has not been peer-reviewed.** Findings, non-claims and limitations:
[`REPORT.md`](REPORT.md). Every identity with its assumptions and normalisation:
[`DERIVATION.md`](DERIVATION.md). Every claim with its class:
[`CLAIMS.md`](CLAIMS.md). Reproduction: [`REPRODUCE.md`](REPRODUCE.md).

---

## The question

For a zero `rho = beta + i gamma` put `delta = beta - 1/2` and `D = 2 delta^2`. On
`0 <= beta <= 1` the reflected-coordinate partition `p = beta`, `q = 1 - beta` has arithmetic
mean `A = 1/2`, geometric mean `G = sqrt(1/4 - delta^2)` and harmonic mean
`M_H = 1/2 - 2 delta^2`, so `D = A - M_H`. Do the gaps between these means say anything about
zero displacement beyond `delta^2`, and does any specified evolution force the sign of their
change?

## The answers

**No, and almost never.**

* `A == 1/2` **is** the functional equation and `G^2 = A M_H`, so the whole mean triple is a
  function of the single scalar `D = 2 delta^2`. `A - M_H` *is* `D`; `A - G` is a monotone
  reparameterisation; and `G - M_H` carries strictly **less** — it is two-to-one, with
  maximum exactly `1/8` at exactly `D = 3/8`, vanishing both at equality (`D = 0`) and at
  complete dominance (`D = 1/2`).
* One identity is worth keeping: in this repository's normalisation **de Bruijn's threshold
  is exactly the largest mean gap**, `tau_c = Delta^2/2 = sup_rho (A - M_H)(rho)`.
* Along the heat flow `d_tau H = -d_z^2 H` the exact chain rule gives
  `D' = Im(z) Im(H''/H')`, and the canonical product (taken in `w = z^2`, where the genus is
  `0`) turns it into
  ```
      D'(tau)  =  -2 y_j  sum_{k != j}  (y_j - y_k)/|z_j - z_k|^2 ,
  ```
  absolutely convergent given a bounded zero strip. A zero's displacement can **grow** only
  if another zero lies strictly further off the critical line on the same side; its conjugate
  partner always contributes exactly `-1`. So the **extremal** zero obeys `D' <= -1`, which
  integrates to de Bruijn's bound and nothing better.
* For a **non-extremal** zero the sign is **refuted**, exactly over `Q`: for
  `P_0(z) = (z^2+1)(z^2+4)` under the polynomial heat flow, `D'(0) = +1/3` at `z = i` while
  `D'(0) = -19/3` at the extremal `z = 2i`.
* Quantitatively, a competing zero can lift `D_j'` above `-1` **only if it shares the
  ordinate to within the target's own displacement**, `|gamma_j - gamma_k| < |delta_j|`
  (and, after optimising, to within `delta_max/sqrt3`, which is sharp) — necessary, not
  sufficient, and the one new statement here that is cheap to apply. The minimal
  counterexample is confluent in the ordinate; no Davenport–Heilbronn census quartet comes
  within a factor of 189 of the threshold.
* **Nothing distinguishes `zeta` from Davenport–Heilbronn.** All four of DH's certified
  off-line quartets have `D'(0) < 0` — interval-certified, between `-1.012` and `-1.161` — at
  the very parameter where the parent extension certifies its Weil form is **negative**. So
  **`D` decreasing is not a positivity mechanism.**

**Nothing here proves or advances the Riemann hypothesis.** Full non-claims in
[`REPORT.md`](REPORT.md) §8.

## A certified Davenport–Heilbronn result

The follow-up audit produced one certified number, and it is **not** about RH:

> **`H_tau` has a nonreal zero at `tau = 4/25 = 0.16` exactly, so `Lambda_DH > 0.16`.**
> Arb at 700 bits, 512 arcs, winding number exactly 1 by the half-plane method, with a
> negative control. This improves the published certified `Lambda_DH >= 6.8268e-29` by a
> factor of `2.3 x 10^27`, leaving the bracket `[0.16, 1.602537…]` — **a factor of 10**
> instead of 28 orders of magnitude.

It needs no trajectory estimate, no isolation hypothesis and nothing uniform over a
parameter interval: one parameter value, one contour, one integer. See
[`AUDIT.md`](AUDIT.md) §4b. **For `zeta`, RH needs `Lambda <= 0` and Rodgers–Tao give
`Lambda >= 0`; a positive lower bound from an existing off-line zero of a different
function says nothing about `zeta`.**

## The follow-up audit

[`AUDIT.md`](AUDIT.md) audits how these constraints combine. It **corrected three claims
above** (the proximity cutoff's sharp constant, its `N`-competitor form — which was wrong —
and the hypothesis of the energy identity), mapped fourteen constraints and found **eleven of
them satisfied by Davenport–Heilbronn**, and applied an explicit stopping rule: **the RH
branch stops there.** One bounded Davenport–Heilbronn deliverable survives — a plausibility
estimate `Lambda_DH >~ 0.0137` against the published certified `6.83e-29`, a factor `2e26` —
and it has since been **certified** at the much stronger value `Lambda_DH > 0.16` by the
fixed-parameter route above. Explicitly not RH work.

## A note on the checkout

The local working tree was on branch `projection-formula` at `16c9f0d`, **not an ancestor of
`origin/main`** (`dfdabab`) and missing most of the work this builds on. The authoritative
state was read from `origin/main` via `git archive` into a scratch directory. **No branch was
switched, nothing was stashed, reset, committed or pushed**, and the pre-existing
modifications and untracked paths were left exactly as found. This directory is untracked and
depends on no other path in the repository. See [`REPORT.md`](REPORT.md) §1.1.

## Layout

```
REPORT.md       findings in plain language, then exact statements and limitations
AUDIT.md        the follow-up constraint audit: corrections, the map, the DH filter, one lemma
DERIVATION.md   every identity, with assumptions, tiers and the normalisation audit
CLAIMS.md       the classification table: proved / certified / numerical / conjectural / refuted
REPRODUCE.md    dependency order and runtimes
src/
  kernels.py          the zeta and DH kernels, H_tau and its z-derivatives, the dictionary
  check_dictionary.py Stage 0: verify the normalisation BEFORE any dynamics  (dps 90)
  means.py            Part 2: the mean-gap algebra, exact over Q
  toy.py              Control 1: the polynomial heat flow and the exact counterexample
  dh_flow.py          Control 2: D'(0) at the four DH census orbits, two routes
  dh_interaction.py   the interaction identity with its measured remainder
  zeta_flow.py        Control 3: real branches, collisions, the emergent pair
  zeta_lehmer.py      Control 3 addendum: why the tracked pair is NOT the extremal zero
  separation_bound.py the necessary ordinate-separation condition (S1d)
  audit1_claims.py    AUDIT Task 1: cutoff sharpness, Lambda vs tau_c, the (S2) hypothesis
  audit2_constraints.py AUDIT Tasks 2-3: the constraint map and the DH filter
  audit4_lemma.py     AUDIT Task 4: the one proposed lemma, measured
  heat_taylor.py      H_tau from the s-derivatives of Lambda_F, with its two rigorous tails
  dh_fixed_tau.py     navigation: the O_1 branch at fixed tau, and numerical winding numbers
  check_taylor_vs_integral.py  the certificate's representation, checked at tau != 0
                      against the defining integral, with a wrong-constant control
  certify_fixed_tau.py CERTIFICATE: a nonreal zero at tau = 4/25, so Lambda_DH > 0.16
  alt_evolution.py    Part 5: the kappa-family, evaluated and rejected
  certify.py          Arb (400-bit) certification of the representative statements
  gq.py               exact Gaussian rationals, for the toy model
tests/
  regression.py       expected values written into the test, plus a negative control
  verify_numbers.py   every headline number in the Markdown, checked against data/
results/              console transcripts, one per script
data/                 JSON results and the certificates
run_suite.sh          the whole thing in dependency order; non-zero exit on any failure
MANIFEST.sha256       sha256 of every file here; `shasum -a 256 -c` verifies clean
```

Every script asserts its own inequalities and residuals and exits non-zero if any gate
fails; none prints a stored label.

## Licence

As the parent repository: **MIT** for code under `src/` and `tests/`, **CC BY 4.0** for the
Markdown, `data/` and `results/`.
