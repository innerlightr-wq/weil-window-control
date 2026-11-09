# positive-weil-certificate-2026-11-09

First build of a positive finite-window Weil certificate pipeline.

**Scope in one line.** This build reproduces a reference detector and a control
table, re-runs the repository's certified negative witness as a regression,
reproduces the *scalar gate* of a published finite-to-infinite reduction and
sharpens its constant, and measures what the next window `L = 81/100` costs.
**It establishes no whole-window lower bound at any `L`.** `m_L` is
**UNRESOLVED** throughout. Nothing here is progress toward RH.

---

## Current result

| claim | status |
|---|---|
| `P(theta)=3+4cos+cos2 = 2(1+cos)^2`; `zeta(s)^3|zeta(s+it)|^4|zeta(s+2it)| >= 1` for `Re s>1`; no zeros on `Re s=1`, `t != 0` | **PROVED** (wrapper inputs W1–W3 named) |
| `t -> Re psi(1/4+it/2)` strictly increasing; hence `Psi_L(t) >= beta_sharp(T#)` for `t >= T#`, with `beta_sharp >= beta_paper` | **PROVED** (G1–G3) |
| exact gate `T_gate(4/5) = 119.086906076…`, `T_gate(81/100) = 502.389536459…` | **CERTIFIED** (Arb, certified bisection) |
| `A_L(4/5) = 2.94197352522…`, `beta*(T#=200) = 0.5134667749…` | **CERTIFIED**, reproduces the external 2.9420 / 0.5135 |
| crossing `L_5=(log 5)/2`, `A_L` jumps by **exactly** `2 log 5/sqrt 5`, so the required `T#` multiplies by **exactly** `exp(2 log 5/sqrt 5) = 4.2186915938…` | **CERTIFIED** |
| `lambda_min(4/5) <= 2.11e-17`, `lambda_min(81/100) <= 7.2e-18` | **CERTIFIED** upper bounds (frozen rational vectors); see `certs/` |
| whole-window lower bound `Q_zeta(f) >= m_L ||f||^2` at `L = 4/5` or `81/100` | **UNRESOLVED** |

An **upper** bound needs one test function and is unconditional. A **lower**
bound needs the whole window. The two are never merged in the reporting, and
`verify.py` check V6 fails if they are.

## Reproduction

```sh
cd extensions/positive-weil-certificate-2026-11-09
P="PYTHONPATH=src:../../src ../../.venv/bin/python"

eval $P src/phase_a.py        # reference detector + control table   (~2 s)
eval $P src/phase_b.py        # the reduction's scalar gate          (~10 s)
eval $P src/phase_c.py        # activation, cost jump, Gram bounds   (~1 min)
eval $P src/upper_bound.py    # frozen vectors -> certified bounds   (~12 min)
eval $P verify.py             # COLD VERIFIER, nonzero exit on failure

# regression control, the repository's existing certified negative witness:
( cd ../four-mode-positivity-2026-10-14 && ../../.venv/bin/python verify_new_certificate.py )
```

Environment: `.venv` Python 3.14.7, `mpmath 1.3.0`, `python-flint 0.9.0` (Arb).
Every sign-critical path is Arb ball arithmetic or exact rationals; mpmath is used
only to *locate* candidate vectors, which are then frozen and certified.

## Layout

```
HYPOTHESES.md    input ledger, control table, Statement-F filter
DERIVATION.md    normalization, reference inequality, the gate proved here
SOURCE_AUDIT.md  exact paper versions by md5/DOI, and what was unavailable
REPORT.md        result, limitations, the next missing estimate
src/wcore.py     exact Chebyshev + Arb comb mass / prime activation
src/phase_a.py   reference detector and controls A-D
src/phase_b.py   one-stroke gate: A_L, T_1, T_gate, beta_sharp vs beta_paper
src/phase_c.py   L=81/100: activation, cost factor, finite Gram upper bounds
src/certify.py   SECOND implementation of the Weil form, in Arb balls
src/upper_bound.py  frozen rational vectors -> certified upper bounds
certs/           frozen inputs and machine-readable results
data/            phase_c.json
verify.py        cold verifier
```

## What this is not

* not a whole-window positivity certificate;
* not a reproduction of the external `L = 0.8` lower bound — the matrix step of
  that reduction was **not** re-run (`SOURCE_AUDIT.md` §3);
* not a new arithmetic estimate: the one sharpening proved here (`beta_sharp`)
  moves a constant by under 1%, not the doubly exponential threshold;
* not a quantitative zero-free region;
* not a priority claim, and no manuscript is packaged.
