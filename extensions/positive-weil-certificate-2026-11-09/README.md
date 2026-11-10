# positive-weil-certificate-2026-11-09

First build of a positive finite-window Weil certificate pipeline.

> # **CERTIFIED: `1.41e-17 <= lambda_Q(4/5) <= 1.674e-17`**
>
> Width `2.64e-18`. **Lower** endpoint: a `T = 400` frequency-cutoff certificate,
> both parity sectors (`REPORT.md` §10). **Upper** endpoint
> `1.67312553407068186807726837736e-17`: the Rayleigh quotient of `Q` at a frozen
> rational witness over **160 even cosine modes**, evaluated in **closed form**
> through digamma, Lerch `Phi` and Hurwitz `zeta` — a route that performs **no
> quadrature at all** (`phase5/src/upper_highN.py`).
>
> Parity separation, on that witness: `1.674e-17 < 2.3e-17`, margin
> `6.2687e-18`, so the **EVEN** sector determines `lambda_Q`. The separation
> itself is not new — Zhu's Thm 6.2 has it at this window with an odd constant
> `357x` stronger.
>
> ### ⚠ The compact time-domain route is WITHDRAWN (`RETRACTIONS.md` R14)
>
> An earlier headline here read `1.41e-17 <= lambda_Q <= 1.9014e-17`, with the
> upper endpoint from the compact time-domain evaluator. **That evaluator omits
> the truncation error of its outer panelled Gauss–Legendre rule** — it bounds
> only the integrable head — so its published radii are incomplete and every
> value it produced is withdrawn **as a certificate** and retained as a
> **measurement**: `q(f_400) = 1.90138291114546440e-17`, `q(f_200)`, and the
> discards `e_T = q - r_T` with their ratios.
>
> **Unaffected:** the lower bound and the odd-sector bound, whose interior
> quadrature carries a certified Bernstein-ellipse remainder folded into the
> matrix entries exactly once; `beta`, `M`, `kappa`, `g_N`, `S_N`; and the scalar
> `R_T` values, whose routine *does* add a certified remainder and publishes it.
>
> The replacement endpoint above is **better** than the withdrawn one by `12.0%`,
> so the defect cost nothing in the bracket. A rigorous outer remainder is still
> **open** and is no longer on the critical path
> (`phase5/src/outer_convergence.py`).
>
> An earlier erratum (`R13`) corrected two mis-quoted endpoints: an enclosure
> **midpoint** quoted as a bound, and a `53`-bit re-parse. Both were conservative.
>
> **A technical note consolidating this build is in [`note/`](note/)**, with a
> claim-to-artifact index, a rounding-direction checker, a Zenodo package and a
> reproducibility archive. It carries **no DOI**: nothing has been deposited.
>
> Superseded headline (kept for provenance):
>
> Lower bound `Q_zeta(f) >= 1e-17 ||f||_2^2` for every admissible `f` supported in
> `[-4/5, 4/5]`, both parity sectors, hence arbitrary complex tests; upper bound from
> the phase-1 frozen-vector artifact. Obtained independently; the external `8.9e-18`
> was never an input. **Audited in phase 4 and it stands** (`REPORT.md` §9); two
> interpretive claims were corrected (`RETRACTIONS.md` R11, R12).
> **This is not RH** and says nothing about any other window.
>
> Measured there: on the even block's own minimizer the reduction discards
> `e_f >= 1.5304e-17` while keeping `r_f = 1.0310e-17` — a loss ratio `>= 1.48`.
>
> **PHASE 2 (2026-11-09).** Baseline whole-window positivity at `L = 4/5` is
> **UNRESOLVED**, but the gap is now a single number: the gate, complement and
> coupling are **CERTIFIED**, and at `N = 220` the chain closes as soon as the
> `N x N` Legendre block has `lam_A >= 1.0146e-28`. That block was not
> assembled. Phase 2 lives in [`phase2/`](phase2/); see `REPORT.md` §7.
>
> **Revised after external review, 2026-11-09.** One claim is **retracted**:
> the proposed sparsity ("power saving") target for the near-resonant set is
> **false** — that set has positive density. Three statements are corrected.
> See [`RETRACTIONS.md`](RETRACTIONS.md).

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
| crossing `L_5=(log 5)/2`, `A_L` jumps by **exactly** `2 log 5/sqrt 5`, so `T_1 = 2pi e^{A_L}` multiplies by **exactly** `exp(2 log 5/sqrt 5) = 4.218691593787…` (this is exact for `T_1` **only** — the digamma gate gives `4.218679895`, the `-1/T` gate `4.192025180`) | **CERTIFIED** |
| `lambda_min(4/5) <= 1.674e-17` (160 cosine modes), `lambda_min(81/100) <= 7.2e-18` | **CERTIFIED** upper bounds (frozen rational vectors), closed-form route, no quadrature; see `certs/` and `phase5/certs/upper_highN.json` |
| gate `beta = H(T)-A_L = 0.51846573324385` at `T=200`; exterior set **empty** | **CERTIFIED** |
| complement `g_220 = 6.7424e-31`, coupling `kappa = 7.2527e-15`, `lam_D = 0.51846573` | **CERTIFIED** |
| Lemmas A and B (finite-measure and uniform-local-occupation energy control) | **PROVED** (standard tools) |
| **whole-window lower bound `Q_zeta(f) >= 1e-17 ||f||^2` at `L = 4/5`** | **CERTIFIED** (phase 3): `a = 1.03e-17`, `d = 0.51846573`, `kappa = 7.2527e-15`, margin `1.5554e-19`, both sectors |
| whole-window lower bound at `L = 81/100` | **not attempted** (sequencing, R5) |

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
eval $P src/refutation.py     # refutation of the retracted sparsity target (~3 min)
eval $P verify.py             # phase 1 COLD VERIFIER, nonzero exit on failure

# phase 2:
PH="PYTHONPATH=phase2/src:src ../../.venv/bin/python"
eval $PH phase2/src/lemmas.py      # Lemmas A,B + informative tests     (~40 s)
eval $PH phase2/src/baseline.py    # certified gate/complement/coupling (~45 s)
eval $PH phase2/src/pilot.py       # local-occupation feasibility       (~5 s)
eval $PH phase2/verify_phase2.py   # phase 2 COLD VERIFIER              (~45 s)
eval $PH phase2/src/run_block.py   # ASSEMBLE + certify the blocks      (~22 min)
eval $PH phase2/verify_phase3.py   # phase 3 COLD VERIFIER (20 checks)  (~3 min)
eval $PH phase3/src/audit.py         # beta monotonicity + K_T identity (~2 min)
eval $PH phase3/src/callback_repro.py # minimal callback reproducer      (~1 min)
eval $PH phase3/src/minimizer.py     # frozen minimizer, r_f and e_f    (~7 min)
eval $PH phase3/verify_audit.py      # phase 4 AUDIT VERIFIER (16 checks)(~2 min)

# phase 5 (T = 400 and the time-domain original form):
PH4="PYTHONPATH=phase2/src:phase4/src ../../.venv/bin/python"
eval $PH4 phase4/src/cutoff.py       # cutoff identity, beta_400          (~1 min)
eval $PH4 phase4/src/eval_f200.py    # Q(f_200) by the compact route      (~7 min)
eval $PH4 phase4/src/run_block400.py # ASSEMBLE + certify at T=400        (~60 min)
eval $PH4 phase4/src/cross.py        # the 2x2 table, Q(f_400)            (~26 min)
eval $PH4 phase4/verify_t400.py      # phase 5 COLD VERIFIER (18 checks)  (~5 min)

# phase 6 (precision ownership, parity separation, scalar pilots):
PH5="PYTHONPATH=phase2/src:phase4/src:phase5/src ../../.venv/bin/python"
eval $PH5 phase5/test_prec.py        # precision-ownership tests (13)     (~1 min)
eval $PH5 phase5/src/parity.py       # T=400 replay + parity separation   (~2 min)
eval $PH5 phase5/src/conventions.py  # degree-2 convention checks (13)    (~1 min)
eval $PH5 phase5/src/pilot800.py     # R_200(f_400) + T=800 scalar pilot  (~3.5 h)

# the CURRENT upper bound (closed form, no quadrature) and its replays:
PH6="PYTHONPATH=src:phase2/src:phase5/src:../../src ../../.venv/bin/python"
eval $PH6 phase5/src/upper_highN.py 160       # locate + certify, 160 modes  (~11 min)
eval $PH6 phase5/src/upper_highN.py --verify  # re-certify from the rationals (~25 s)
eval $PH6 phase5/src/replay_upper.py          # 28-mode fallback replay       (~1 s)
eval $PH5 phase5/src/outer_convergence.py     # R14 diagnostic, NOT a bound   (~26 min)
eval $PH5 phase4/src/emit_certificate_T400.py  # rewrite certificate_T400   (~1 min)

# the technical note's own checker (rounding direction, scope, parity labels):
PYTHONPATH=phase2/src ../../.venv/bin/python note/consistency_check.py  # 47 checks

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
RETRACTIONS.md   corrections after review (R1 density target; R6 measure/modulation)
note/            technical note, claim index, Zenodo package, archive
phase5/          precision ownership, parity separation, scalar pilots
phase4/          T=400 certificate; compact time-domain evaluation of Q
phase3/          audit: beta monotonicity, callback reproducer, frozen minimizer
phase2/          baseline reduction, Lemmas A/B, occupation pilot, verifiers
  src/block.py   Legendre assembly, certified GL quadrature, Rump PD verification
  certs/         frozen matrix enclosures + whole_window_certificate.json
src/wcore.py     exact Chebyshev + Arb comb mass / prime activation
src/phase_a.py   reference detector and controls A-D
src/phase_b.py   one-stroke gate: A_L, T_1, T_gate, beta_sharp vs beta_paper
src/phase_c.py   L=81/100: activation, cost factor, finite Gram upper bounds
src/certify.py   SECOND implementation of the Weil form, in Arb balls
src/upper_bound.py  frozen rational vectors -> certified upper bounds
src/refutation.py   Kronecker-Weyl refutation of the retracted sparsity target
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
* not a priority claim. A manuscript **is** packaged in [`note/`](note/) as a
  technical note / preprint; it is not peer-reviewed, carries no DOI, and
  claims none of the above.
