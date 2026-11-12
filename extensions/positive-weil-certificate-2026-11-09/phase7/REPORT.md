# RH / Weil-positivity continuation — 12 November 2026

Investigation directory `~/scratch/rh-connes-continuation-2026-11-12/`.
Read-only detached worktree of the research branch in `wt/`. No repository file
was modified, no branch moved, nothing committed, pushed, deposited or sent.

---

## 0. Plain language

The repository holds a rigorous, machine-checked statement about one fixed
window: the Weil quadratic form, restricted to test functions supported in
`[-4/5, 4/5]`, is bounded below by a specific tiny positive number. The new
Connes-and-coauthors paper studies the *same* operator — written in a
multiplicative variable instead of an additive one — over a whole family of
windows, and conjectures that its smallest eigenvalue is comparable to a
quantity coming from a completely different operator (the "uncertainty
operator", whose eigenvalues are classical prolate-spheroidal quantities,
computable exactly).

Three things came out of today.

1. **The two normalizations are the same form.** Their window parameter
   `lambda` is exactly `e^L` where `L` is ours, and after an archimedean
   identity that works out to a clean `log 4`, their quadratic form equals ours
   on the nose — no stray factor of two, no sign left free. I proved this and
   then checked it numerically: an independent implementation of their formulas
   reproduces the repository's own certified values to about thirty digits.

2. **Our certificate turns out to be the one instrument that can test the
   unproven half of their central hypothesis.** Their numerics compute
   eigenvalues of truncated matrices, which — as they say themselves — can only
   ever give *upper* bounds on the true eigenvalue. The lower half of their
   Hypothesis 3.4 is therefore untested by any computation in the paper. Our
   pipeline produces exactly the missing side. Computing their prolate rung
   `delta_4` exactly at our window and dividing, the certified bracket becomes
   a certified bracket for *their* ratio: between 6.40 and 7.47, against the
   7.06 they report at the nearest window they tabulate. That is a rigorous
   instance of the hypothesis's lower half.

3. **Two certified numbers improved today for free.** The matrices were never
   reassembled; the frozen enclosures already on disk simply had more in them
   than was extracted. The even sector moved up 1.6 %, which turns out to
   exhaust the `T = 400` route completely. The odd sector moved up by a factor
   of **577**, overtaking the best previously published odd bound for this
   window by a factor 1.69.

What did **not** happen: nothing here proves the Riemann Hypothesis, says
anything about any larger window, or closes any of the paper's three open
hypotheses. The route the prompt proposed — reorganizing the finite block along
prolate directions — was tested and found **inapplicable**, for a measurable
reason given in §5.

---

## 1. Authoritative local state

### 1.1 Git facts (recorded, not inferred from timestamps)

| item | value |
|---|---|
| repository | `/home/elias/GitHub/weil-window-control` |
| checked-out branch | `dh-note-ziad-citation-2026-11-10` @ `f510c13d2bc2bc1e9e3256397afa81bdac078070` (**= `origin/`**) |
| research branch | `positive-weil-certificate-2026-11-09` @ `7515396c7f2aa17b53d565a4f11f8e937ee164dc` (**= `origin/`**, byte-identical) |
| `main` | `ad78c59a6ca7d2a376c86220e583c5ae94e4962b` (= `origin/main`) — matches the prompt; its root README is **not** the positivity status |
| worktrees before today | one (the main checkout) |
| worktree added today | `~/scratch/rh-connes-continuation-2026-11-12/wt` @ detached `7515396` |
| dirty tracked files | `README.md` only (+56 lines, the DH-note section on the dh-note branch) |

### 1.2 Is local ahead of GitHub?

**For the positivity work: no.** All 106 files of
`extensions/positive-weil-certificate-2026-11-09/` exist only inside git, at
`7515396`, identical to `origin`. The 19 files physically present in that
directory in the working tree are **build products only**: 4 LaTeX
`aux/log/out`, `weil-L45-certificate-1.2.0.tar.gz`, and 14 `__pycache__`
`.pyc` files. No uncommitted source, certificate or note revision.

**Elsewhere: yes, four directories are ahead** — untracked, present in **no**
branch, tag or remote:

| directory | date | subject | bears on today? |
|---|---|---|---|
| `extensions/odd-sector-residual-margin-2026-11-02` | 3 Nov | **artificial finite background** `Q_fin`, planted quartets, `kappa_odd`, synthetic `Gamma = {1,3,7/2}`; its own overview says "everything below is a finite-model result … yields no RH or GRH consequence" | **No.** This is the *artificial added-quartet* object. Its "odd sector" is the odd subspace of a synthetic model, **not** the odd parity sector of the zeta window form. Kept strictly distinct, as the brief requires. |
| `extensions/spectral-isometry-transfer-2026-11-02` | 2 Nov | audit of arXiv:2610.40470v1 — verdict **NO TRANSFER**, with a proved limitation | no |
| `extensions/suzuki-function-field-window-2026-11-02` | 2 Nov | audit of arXiv:2610.41115v1 — exact operator dictionary, **no arithmetic transfer** | no |
| `extensions/newtonian-weil-response-2026-10-08` | 8 Oct | earlier response computation | no |

No file dated 9–12 November exists outside the research branch except the
LaTeX build products above and the Connes PDF in `~/Downloads`. **No prior work
on arXiv:2611.13129 exists on this machine** (grep for `2611.13129`,
`La Gatta`, `uncertainty operator` over `~/Downloads ~/scratch ~/Documents`
returns only the prompt file itself).

### 1.3 Verification actually run

The documented **fast check** (`note/README.md`, "verify everything from the
frozen artifacts"), in the pinned venv (`python-flint 0.9.0`, `mpmath 1.3.0`,
Python 3.14.7). Log: `logs/fastverify.log`.

| verifier | exit | checks |
|---|---|---|
| `verify.py` | 0 | all pass |
| `phase2/verify_phase2.py` | 0 | all pass |
| `phase2/verify_phase3.py` | 0 | all pass |
| `phase3/verify_audit.py` | 0 | all pass |
| `phase4/verify_t400.py` | 0 | all pass |
| `phase5/test_prec.py` | 0 | all pass |
| `note/consistency_check.py` | 0 | 52 pass, 0 fail |
| **total** | **all 0** | **176 `[PASS]`, 0 `[FAIL]`** |

Two notes on what this is:

* It is a **frozen-artifact replay, not a cold rebuild.** No matrix is
  reassembled. The one piece of genuine recomputation is `verify_t400` check
  T4, which re-runs Rump's preconditioned positive-definiteness test on the
  frozen `220 x 220` enclosures and re-runs `certify_final`; `beta`, `M`,
  `S_N`, `g`, `kappa`, `d` are recomputed from `baseline` rather than read.
* `note/README.md` advertises **166** checks; the sequence actually emits
  **176**. The note's table is stale by ten checks (all passing). Minor
  documentation defect.

**Regression control** (`extensions/four-mode-positivity-2026-10-14/verify_new_certificate.py`):
exit 0, `NEW CERTIFICATE: ALL ASSERTIONS PASS` — 69 cells of certified
Davenport–Heilbronn **negativity** `Q_add,delta(f_x) < 0` on `delta in [9/100, 1/4]`,
*and* `Q_zeta(frozen) > 0` recomputed (`lo = 5.124239331e-03`). The DH negative
witnesses and actual-zeta positivity remain distinct objects, as required.

### 1.4 Evidence-backed status table

Verification levels: **M** manuscript assertion only · **D** saved certificate
data · **R** successful replay of a frozen artifact · **I** independently
regenerated here today · **C** conditional theorem · **N** numerical
observation.

| # | claim | precise scope | file / commit | verifier | level | unresolved dependency |
|---|---|---|---|---|---|---|
| 1 | `Q(f) >= 1.41e-17 ||f||^2` | every admissible `f` supported in `[-4/5,4/5]`, both parity sectors, hence arbitrary complex tests | `phase4/certs/certificate_T400.json:whole_window_mu` @ `7515396` | `phase4/verify_t400.py` T4; `consistency_check` T1–T8 | **R** | superseded today by #2 |
| 2 | **`Q(f) >= 1.4330004882800e-17 ||f||^2`** | same scope | frozen `phase4/certs/block400_even.json` @ `7515396`; new shift in `work/sharpened_T400.json` | `work/push_shift.py`, `work/finalize.py` (Rump PD, `||I-G|| = 1.12e-12 < 1`), then `certify_final` | **I** | none beyond the imported normalization `W4` |
| 3 | `Q(f) >= 1.3857e-14 ||f||^2` on the **odd** sector | odd parity only, `L = 4/5` | same frozen odd block | same | **I** | none; **577x** the committed `2.3e-17`, **1.69x** Zhu Thm 6.2's `8.206e-15` |
| 4 | `lambda_Q(4/5) <= 1.67312553407068186807726837736e-17` | one frozen rational witness, 160 even cosine modes, closed form, **no quadrature** | `phase5/certs/upper_highN.json:outward_upper_endpoint` | `phase5/src/upper_highN.py --verify`; `consistency_check` U3–U6d | **R** | none |
| 5 | `lambda_min(A_even, T=400)` is pinned to rel. width `5.0e-5` | the `440`-mode compression of `R_400` | floor #2 vs ceiling `R400_f400 = 1.4330724047355850e-17` in `phase4/certs/cross.json` | two independent routes (Rump PD; closed-form Rayleigh quotient at the frozen minimizer) | **I** | — |
| 6 | `lambda = e^L`, and `QW_lambda(g,g) = Q(f)` exactly | the paper's eqs (6)(7)(8) vs this repository's form | `notes/DICTIONARY.md` §1; `work/dictionary_check.py` | analytic proof + 30-digit agreement on `phi_0`, `phi_2`, `phi_0-phi_2` against `phase5/certs/conventions.json` | **I** | none |
| 7 | exact rungs `delta_4(e^{4/5}) = 2.2391127080494187067e-18`, `delta_6 = 3.4175764349345855425e-15` | Theorem 2.1 rungs of `U_lambda`, not Weil eigenvalues | `work/prolate_rungs.py` | matches the paper's printed 30-digit `chi_0`, `chi_4` at `lambda^2=5` to `6.4e-30` / `1.6e-30`; stable under `K: 140->280` and `dps: 80->110`; `gamma->0` gives `mu_n = n(n+1)` | **I (high-precision, not interval)** | an interval enclosure of `delta_4` is the one uncertified input in #8 |
| 8 | **`6.3998587 <= epsilon_lambda(1)/delta_4(lambda) <= 7.4722703`** at `lambda = e^{4/5}` | even sector; a certified instance of the **lower** half of Hypothesis 3.4 | #2, #4, #6, #7 | derived in `work/finalize.py` | **I**, conditional on #7's precision | #7 not yet an interval |
| 9 | `epsilon^{odd}_lambda(1)/delta_6(lambda) >= 4.0546277` | odd sector, same `lambda` | #3, #6, #7 | same | **I**, same caveat | same |
| 10 | gate: `T_gate(4/5) = 119.086906076…`; `Psi_L >= beta` on `|t| >= T`, exterior set empty | `L = 4/5`, `T = 200, 400, 800` | `phase4/certs/cutoff.json`, `phase5/certs/pilot800.json` | `verify_t400` T1; `consistency_check` G1–G3 | **R** | — |
| 11 | `q(f_400) = 1.90138291114546440e-17`, `q(f_200)`, `e_T`, discard ratios | **MEASUREMENTS, not certificates** | `phase4/certs/cross.json`, `phase5/certs/pilot800.json` | — | **N** | R14's outer Gauss–Legendre remainder, still open |
| 12 | `lambda_R(800) <= R_800(f_400) = 1.6681e-17`; the whole `T=800` route is capped at `+18.3%` over `1.41e-17` | a **ceiling** on the comparison, not an achieved bound | `phase5/certs/pilot800.json` | `consistency_check` H3–H4 | **D/N** | recompute the headroom against #2: see §6 |
| 13 | even-sector derived fields of `phase4/certs/finite_block_certificate_T400.json` are corrupt (`d = [+/- 2.75e+183]`, `kappa = nan`, `final = null`) | provenance record only | that file | disclosed in `CLAIM_INDEX.md`; no claim reads it | **D** | preserve unaltered |
| 14 | Thm 8.14 / Lemma 8.16 / Prop 8.18 of arXiv:2611.13129v1 | reduction of RH to a Schur block of dim `O(A sigma^2 log lambda)+2` | `notes/DICTIONARY.md` §2 | reading | **C** | (H2'), (R1), (R2), and an asymptotic-family hypothesis |

### 1.5 Contradictions resolved by inspecting producers, not prose

* **`CLAIM_INDEX.md` §5 vs §7 on the discard ratios.** §5 labels them
  "**MEASUREMENTS**"; §7 says they "are certified *enclosures*, but no decay law
  is inferred from them and they enter no bound." **§5 is right and §7 is
  wrong as written.** The ratios are quotients of `e_T = q - r_T`, and `q` came
  from `phase4/src/compact.py::arch_compact`, whose returned radius omits the
  truncation error of its outer panelled Gauss–Legendre rule (R14). An
  arithmetic-only radius on a quantity whose rule error is unbounded is not an
  enclosure. The ratios are measurements; §7's clause "are certified
  enclosures" should be struck. (`consistency_check.py` D1–D3 check the
  *arithmetic* of the division, which is correct — that is not the same as
  certifying the inputs.)
* **`README.md` `kappa = 7.2527e-15` vs `CLAIM_INDEX` `kappa = 1.3083e-29`.**
  Not a contradiction: the first is `T = 200`, the second `T = 400`.
* **`README.md`'s "~1.6 % of room left" at `T = 400`.** Superseded today:
  the real figure is `5.0e-5` (§4.2). The `1.6%` was the gap from the
  *committed* `1.41e-17` to the ceiling, not from `lambda_min(A)`.

---

## 2. The corrected baseline, preserved

Carried forward unchanged from the research branch, with the corrections
restated so they cannot be lost:

* **`R14` stands.** The time-domain compact evaluator
  (`phase4/src/compact.py::arch_compact`) returns `const + tot + head`, where
  `head` bounds only the `[0, eta]` piece; **no term bounds the truncation
  error of the panelled Gauss–Legendre rule**. Ball arithmetic cannot supply
  it: the gap between a Gauss sum and its integral is a property of the rule,
  not of the arithmetic. Withdrawn **as certificates**, retained as
  measurements: `q(f_400) = 1.90138291114546440e-17`, `q(f_200)`,
  `U = 1.901382911145464406030000e-17`, the bracket
  `1.41e-17 <= lambda_Q <= 1.9014e-17`, and every `e_T = q - r_T` with its
  ratios. The alternative upper route (`phase5/src/outer_convergence.py`)
  **does not repair it** — it records a refinement difference (`+/- 2.71e-35`),
  explicitly not a bound, and explains why the naive Bernstein route fails at
  `deg g = 877`.
* **`R13` stands.** A ball midpoint quoted as the bracket's upper endpoint, and
  a 53-bit re-parse of a serialized ball. Both conservative; both type errors
  that precision cannot detect. A lower bound is consumed through `.lower()`,
  an upper bound through `.upper()`, never through a midpoint or a printed
  decimal. **This rule was applied to every new number below** — see
  `work/finalize.py`, which re-certifies each sharpened shift and then prints
  the lower endpoints.
* **`phase4/certs/certificate_T400.json` is the authoritative whole-window
  record**; `phase4/certs/finite_block_certificate_T400.json`'s even-sector
  derived fields are corrupt and retained unaltered for provenance only.
* **`T = 800` has scalar pilot information, not a certificate.** Its `18.3 %`
  figure is a *ceiling* on what the comparison can reach, not an achieved
  improvement.
* **Credit.** The finite-to-infinite reduction is Zhu's (arXiv:2608.24827v2).
  Checked against `SOURCE_AUDIT.md` and `note/README.md` v1.0.2: Zhu's
  Theorem 6.2 certifies the parity **separation at this very window** with
  `lambda_odd_1 >= 8.206e-15`, and establishes that the ground state is simple
  and even. **That odd constant was 357x stronger than this repository's
  `2.3e-17`. After §4 it is 1.69x *weaker* than the `1.3857e-14` certified
  today** — the first time this repository's odd-sector constant has passed the
  published one. The simplicity and evenness results remain Zhu's alone and are
  not reproduced here. (Zhu's *matrix* step was never re-run here;
  `SOURCE_AUDIT.md` §3.)
* `A_L` optimality for any pointwise comb envelope (R2), and the refutation of
  the density/"power saving" target (R1, positive density by continuous
  Kronecker–Weyl) both stand and were not revisited.

---

## 3. The paper: what it is and what it is not

See `notes/DICTIONARY.md` for the full theorem-to-gap table with page and
equation references. The three distinctions the brief asks for, confirmed:

* **Theorem 8.14 (p. 47) is conditional on (H2') at scale `sigma`** — a
  semilocal localization hypothesis. Remark 8.15(b) says that in the
  *archimedean* case (`S = {infinity}`, `lambda^2 < 2.25`, no prime enters) the
  stronger form of (H2') is known and no absorption is needed. That settles
  nothing semilocal: the archimedean case is precisely the regime with no
  arithmetic in it.
* **Lemma 8.16 (p. 48) is a Schur reduction, available only once `D >= d > 0`
  on the *whole* complement `W`** and `B*v` lies in the dual of `D`'s form
  domain. Its coercivity input is Theorem 8.14, hence (H2').
* **Proposition 8.17 (p. 49)** has two halves. The monotonicity
  (`epsilon_lambda` non-increasing in `lambda`, because the form on a larger
  window restricts to the form on a smaller one) is **unconditional and
  usable**. The implication — `epsilon_lambda >= -eta(lambda)` for **all large
  `lambda`** with `eta -> 0` gives RH — is a statement about an asymptotic
  family. **No finite collection of window computations, however rigorous,
  supplies it.** Our result at one `lambda` is not a step along that hypothesis.
* **Proposition 8.18 (p. 50)** additionally needs (R1) (proved in two lines),
  **(R2) (the authors write "we leave this open")**, Gram/uniform estimates, and
  the sign of the Schur complement on the transition-and-pole block for every
  large `lambda`. Assuming positivity of the conditioned form to conclude
  positivity is explicitly *not* what they do — but `V'` is where all the
  arithmetic is compressed, and its sign is open.
* **Fuchs' asymptotics (eq. 17) are fixed-`n`.** The paper is careful about
  this and says the comparison "must use the exact `D_j`" (§4.2 item 1): at
  `lambda^2 = 5` Fuchs is off by `1.0` in `D_2` in absolute terms. §8.6 states
  explicitly that "Fuchs' fixed-`n` asymptotics plays no role here". Everything
  below therefore uses the **exact** rungs from Theorem 2.1 / Appendix A, never
  eq. (17).
* **Their recovery of many zeta zeros (§7) and the real zeros of approximating
  transforms do not establish the global convergence their argument needs**, and
  are not used here.

---

## 4. What was tested, and what came out

### 4.1 The dictionary — PASSED (§1.4 row 6, `notes/DICTIONARY.md` §1)

Derived before any numerical comparison, as the brief requires. Analytic proof
plus a fully independent reimplementation of the paper's eqs. (6)(7)(8)
agreeing with the repository's certified values to ~30 digits on three test
functions. The archimedean identity reduces to
`int_0^inf (1-e^{-s})/sinh s ds = 2 log 2`.

Byproduct: a **divergent series** as printed in `DERIVATION.md` §1 (prose only;
nothing computed uses it). Diagnosis and the correct grouping in
`notes/DICTIONARY.md` §1.5.

### 4.2 Sharpening the frozen `T = 400` blocks — PASSED, and it closes the route

The brief asked whether a reorganization of the finite block improves the
certificate or exposes its limiting estimate. The honest first question is
whether the *existing* block had been read out tightly, and it had not. Nothing
was reassembled: `phase4/certs/block400_{even,odd}.json` were consumed as
given, and the Rump shift was pushed by doubling and then 14 bisections
(`work/push_shift.py`, log `logs/push_shift.log`, 4.2 s per test).

| sector | committed `a` | certified today | factor | Rump `||I-G||` at the new shift | sector `mu` after the Schur chain |
|---|---|---|---|---|---|
| even | `1.42e-17` | **`1.4330004882800e-17`** | `1.0092` | `1.12012e-12` | `1.4330004882800e-17` |
| odd | `2.4e-17` | **`1.3857e-14`** | **`577.4`** | `5.65580e-18` | `1.3857e-14` |

Both re-certified a second time in `work/finalize.py`, and both quoted through
the **lower endpoint** of the resulting ball. `kappa^2/d = 1.4127964e-58`, so
the Schur penalty is `58` orders below the shift and `mu = a` to within
`1e-58`.

**Built-in negative control.** The test is not vacuously passing: the same
routine **fails** just above each certified shift — at `1.43308715820e-17`
(even) and `2.4576e-14` (odd), by midpoint Cholesky breakdown. The bisection is
tight to `~5e-5` relative in both sectors. This is the control appropriate to
*this* decomposition (an upper edge on the shift ladder); the separate planted
negative witness is §1.3's regression control.

**New whole-window result at `L = 4/5`** (even sector binds, as before):

```
    1.4330004882800e-17  <=  lambda_Q(4/5)  <=  1.67312553407068186807726837736e-17
```

width `2.4013e-18` (was `2.64e-18`, `-9.1 %`); lower endpoint up `1.63 %`.
Parity separation `mu_odd / U = 828`, margin `1.38403e-14` (was `6.27e-18`).

**Why this closes the `T = 400` route.** Two independent computations now meet:
the certified floor `1.43300048828e-17` (Rump PD on the assembled block) and
the independently measured ceiling `R_400(f_400) = 1.4330724047355850e-17`
(closed-form Rayleigh quotient at the frozen `T=400` minimizer,
`phase4/certs/cross.json`). They bracket `lambda_min(A_even, T=400)` to a
**relative width of `5.0e-5`**. The `T = 400` comparison is therefore
*exhausted*, not merely "about 1.6 % from exhausted" as `REPORT.md` §10.3 says.

**Consequence for the `T = 800` headroom.** `lambda_R(800) <= 1.6681e-17`
(`phase5/certs/pilot800.json`). Against the new lower endpoint the entire
`T = 800` route is worth at most a factor `1.6681/1.43300 = 1.164`, i.e.
**+16.4 %**, not `+18.3 %` — and strictly less in practice, against a measured
assembly cost of order ten hours. The `CONSOLIDATE at T = 400` recommendation
of `REPORT.md` §11.6 is strengthened, not weakened.

### 4.3 The exact prolate rungs — PASSED (§1.4 row 7)

Appendix A of the paper (Bouwkamp's method) reimplemented in `mpmath`:
symmetric tridiagonal eigenproblem (75) on even Legendre indices,
`chi_n = lambda sqrt 2 d_0 / phi_n(0)` from eq. (78). Validation ladder, all
passing:

* `gamma -> 0` gives `mu_n = 0, 6, 20, 42 = n(n+1)`;
* at `lambda^2 = 5`, `chi_0` and `chi_4` match the paper's own printed
  30-digit values to `6.4e-30` and `1.6e-30`;
* `-log delta_4 = 41.18540817` and `-log delta_6 = 33.83264394` at
  `lambda^2 = 5` reproduce their table's `D_1 = 41.185` and `33.833`;
* stable to all 20 printed digits under `K: 140 -> 280` and `dps: 80 -> 110`.

At our window `lambda = e^{4/5}` (`lambda^2 = 4.95303242439511480365`):

```
    delta_4 = 2.2391127080494187067e-18      (-log = 40.640452)
    delta_6 = 3.4175764349345855425e-15      (-log = 33.309845)
```

These are high-precision and validated, **not interval-certified** — the single
uncertified input to §4.4.

### 4.4 The result: a certified instance of Hypothesis 3.4's lower half

The paper's Fact 4.1 is explicit that because `QW^N_lambda` is a *compression*,
min–max gives `epsilon^N_lambda(1) >= epsilon_lambda(1)`, so **every ratio the
paper computes is an upper bound for the true one**. The lower halves of
Hypothesis 3.4 (eq. 32) are consequently untested by any computation in the
paper, and cannot be tested by that method. Combining §4.2 with §4.3 and the
dictionary:

| quantity at `lambda = e^{4/5}`, `lambda^2 = 4.953032` | certified here | paper, nearest tabulated window `lambda^2 = 5`, `N = 960` |
|---|---|---|
| `epsilon_lambda(1)` (even) | `[1.4330004883e-17, 1.6731255341e-17]` | `9.1663e-18` (upper bound) |
| **`epsilon_lambda(1)/delta_4(lambda)`** | **`[6.3998587, 7.4722703]`** | `7.06` (upper bound, still decreasing in `N`) |
| `epsilon^{odd}_lambda(1)` | `>= 1.3857e-14` | `9.8065e-15` (upper bound) |
| **`epsilon^{odd}_lambda(1)/delta_6(lambda)`** | **`>= 4.0546277`** | `4.84` (upper bound) |

Transporting the paper's ratio to our `lambda` at face value gives
`7.06 x delta_4 = 1.58081e-17`, which lies at **61.6 %** of our certified
bracket — inside it, comfortably. In the odd sector `4.84 x delta_6 =
1.65411e-14`, and our certified lower bound is `0.838` of that. **Nothing is
contradicted in either direction, and the two sides are mutually confirming:**
their number is the first independent check that this repository's certificate
is not an artifact of its own reduction, and our certificate is the first
rigorous check that their eigenvalue is not an artifact of truncation from
below.

**Stated as a result.** At `lambda = e^{4/5}` the lower half of Hypothesis 3.4
holds in the even sector with `c = 6.39`, and in the odd sector with
`c = 4.05` — certified, modulo an interval enclosure for `delta_4`, `delta_6`.
This is one point, not a family; it is **not** evidence for the hypothesis
asymptotically, and `Remark 3.5` of the paper makes clear the hypothesis is
strictly stronger than window positivity.

### 4.5 Error budget for the tested decomposition

The decomposition is `H = V_N (+) V_N^perp` with `V_N` the first `N = 440`
Legendre modes (`220` per parity), `A` the exact compression of `R_400`, `D`
its restriction to `V_N^perp`, `B` the coupling; `S = A - B D^{-1} B*`.
`D >= d I` is established on the **full** relevant complement before any
inverse is taken (`DERIVATION.md` §5 steps 4–5), and `B*v` pairs with `D`'s
form domain through the Cauchy–Schwarz surrogate `|F(t)|^2 <= 2L S_N(tL)`.
Positivity of `A` and `D` separately is **not** the argument: the three strict
interval signs `a > mu`, `d > mu`, `(a-mu)(d-mu) > kappa^2` are.

| budget line | value | where it enters | margin against the certified shift `a ~ 1.43e-17` |
|---|---|---|---|
| **Gram normalization** | **exactly zero** | `phi_n(u) = sqrt((2n+1)/(2L)) P_n(u/L)` is orthonormal by construction; the Gram matrix is `I`, not an approximation | — |
| finite-block entries | max radius `1.240e-48` (even), `1.340e-48` (odd), Bernstein-ellipse quadrature remainder (`b_e = 0.30`, `rho = 5.0`) folded in **once** | `A` | `31` orders |
| prime sum | **exactly zero truncation** — `C_{jk}` is supported in `|u| <= 2L`, so only `log n < 2L` contributes | `A` | — |
| archimedean sum | **exactly zero truncation** — closed form in digamma / Lerch `Phi` / Hurwitz `zeta` | `A` | — |
| prolate approximations | **none in the certified chain.** The prolate quantities `delta_4`, `delta_6` appear only in the *comparison* of §4.4, where the error is `mpmath` precision (validated: stable under `K: 140->280`, `dps: 80->110`; matches the paper's 30 digits) | §4.4 only | the one uncertified input |
| pole terms | `|p|^2 = 2 sinh L` exactly; tail `|p_perp| = 7.9934e-1282` | `kappa`, and `d` in the odd sector where `s = -1` is **carried, not dropped** | `1265` orders |
| mixed / coupling block | `kappa = 1.3083437712853882590e-29` | `kappa` | `12` orders |
| omitted modes (complement) | `S_N(320) = 9.2601e-63`, `g_N = 1.8864e-60` | `d = beta - M g_N` | `43` orders |
| **inverse / Schur penalty** | `kappa^2/d = 1.4127963647515877234e-58` | the `(a-mu)(d-mu) > kappa^2` test | **`41` orders** |
| preconditioned PD residual | `||I - R^{-T}(A-aI)R^{-1}||_inf = 1.12012e-12` (even), `5.65580e-18` (odd) | Rump's criterion, `< 1` required | `12` / `18` orders |
| **exterior discard** | `lambda_Q - lambda_R(400)`: `<= 2.4013e-18` certified; `~1.4774e-18` (`9.3 %`) against the paper's transported value | the only term of order the result itself | **this is the binding estimate** |

Every sign-critical line is Arb ball arithmetic or exact rationals. `mpmath`
locates candidate vectors, which are frozen as exact rationals and then
certified; it appears in no pass/fail decision. High precision alone is treated
as a diagnostic (`phase5/test_prec.py` enforces precision *ownership*: a call
without an explicit `prec=` is a `TypeError`, and `constants()` at 200 bits is
refused rather than silently wrong).

---

## 5. Why the proposed prolate-preconditioning route is inapplicable — measured, not asserted

The working candidate in the brief was to use prolate / arithmetic-map
directions to organize or precondition the finite block while keeping rigorous
complement and coupling control. It was tested and it cannot help, for three
reasons, all now quantitative:

1. **A basis change cannot move the number being certified.** `a` is the least
   eigenvalue of the *exact compression* of `R_T` to `V_N` — a basis-independent
   quantity. Reorganizing `V_N` along prolate directions changes the matrix, not
   its spectrum. (The brief says this; §4.2 confirms there is nothing left to
   gain anyway.)
2. **The finite block is no longer the loss.** After §4.2 the certified floor
   and the independently measured ceiling pin `lambda_min(A_even, T=400)` to a
   relative width of `5.0e-5`. There is no room in `A` to recover.
3. **Conditioning was never the binding constraint.** Rump's residual at the
   certified shift is `1.12e-12`, twelve orders below the `< 1` threshold, and
   the preconditioner already replaces the naive interval-`LDL^T` radius growth
   (a factor `~6.5` per pivot, `O(1)` by step 47) with a factor of order
   `cond(R) ~ 1e9`. Better-conditioned directions would improve a number that
   is not close to failing.

Where the loss actually sits: **the exterior discard**. The reduction throws
away `(1/2pi) int_{|t| >= T}(Psi_L - beta)|F|^2 >= 0` entirely. At `T = 400`
that costs about `9.3 %` of `lambda_Q` (`~1.48e-18`), against `5.0e-5` for the
block. Prolate conditioning addresses the block; the discard is untouched by
any choice of basis inside `V_N`.

There is one place where the prolate structure is genuinely the right tool, and
it is the paper's own §8.8 item 3: an audit of the three bounds (73) with the
Gram conditioning, the pole block and the mixed terms separated, with rigorous
error bounds. That is a **different** decomposition — `V_low (+) V_tr (+) C rho`
with `V_low` the deep-ladder conditioning vectors — and it needs the `E`-map
(Gevrey repair, alias sums, Poisson defect) which this repository does not
implement. Note also that it would require exactly the Gram lower bound the
authors leave open (§8.8 item 2), whereas the Legendre decomposition used here
has Gram `= I` **exactly**. Scoping that is a separate project, not a session.

---

## 6. What remains uncertified

| # | item | status |
|---|---|---|
| 1 | `delta_4(lambda)`, `delta_6(lambda)` | high-precision + validated against the paper's own 30 printed digits and stable under `K`/`dps` doubling, but **not an interval enclosure**. The single uncertified input to §4.4. |
| 2 | the outer Gauss–Legendre truncation remainder of `phase4/src/compact.py` (R14) | **open**. Off the critical path (the certified upper endpoint is quadrature-free), but it means `q(f_T)`, `e_T` and all discard ratios stay measurements. |
| 3 | `lambda_Q(4/5)` itself, between `1.4330004883e-17` and `1.6731255341e-17` | a `16.8 %` window. The paper's transported value `1.58081e-17` sits inside it but is not a bound. |
| 4 | the `T = 800` whole-window certificate | not attempted; capped at `+16.4 %` over today's lower endpoint, at order ten hours' assembly. |
| 5 | `L = 81/100` and every other window | not attempted. `A_L` jumps by exactly `2 log 5 / sqrt 5` above `L_5 = (log 5)/2`, taking `T_gate` from `119.09` to `502.39`. |
| 6 | the paper's (H2'), (R1-as-used), (R2) | open in the paper; (R2) explicitly left open by the authors. |
| 7 | Fact 5.3 (prolate conditioning leaves `mu_lambda ~ 2/lambda`) | numerical in the paper, `N`-truncated, **not certified anywhere**. This is §7's target. |
| 8 | `W4`, the repository's explicit-formula normalization | still an **imported** hypothesis, cross-checked only through `Q[0,0]` — but now also, independently, through the 30-digit agreement of §4.1 with a reimplementation of a *different* published normalization. That is a genuine strengthening of `W4`'s evidence, not a proof of it. |

**In flight at the time of writing** (`logs/run_block_lam5.log`,
`certs/certificate_T400_lam5.json`): a certified whole-window bound at
`L~ = 80471895621705018/100000000000000000`, which is strictly below
`L* = (log 5)/2` and within `1.3e-17` of it, giving
`lambda^2 = 4.999999999999999927`. The pre-checks passed: the prime comb is
`{2,3,4}` with `n = 5` excluded (so `A_L~ = A_{4/5} = 2.941973525223620455503909`
exactly, and `T = 400 > T_gate = 119.09` remains admissible, which is why this
costs one assembly rather than ten hours), and `Psi_L(400) >= beta_400`. This
is the **direct head-to-head** with the paper's `lambda^2 = 5` row
(`epsilon_lambda(1) = 9.1663e-18`, `epsilon^odd = 9.8065e-15`, both `N = 960`
upper bounds), using their own printed `delta_4 = 1.2983880087183625604e-18`
and `delta_6 = 2.02613988982e-15`. Expected `~8.3e-18` (even) and `~8.2e-15`
(odd) if today's ratios carry over. **Falsification: if the certified lower
bound at `L~` exceeds `9.1663e-18` (even) or `9.8065e-15` (odd), then one of
the two computations is wrong**, since by min–max their truncated eigenvalue is
an upper bound for the same quantity and the transport over `1.3e-17` in `L` is
nil. That is a real contradiction test, and the first one available between
these two programmes.

---

## 7. The next target as it stood mid-session

> **Superseded during the session.** The coercivity half of Lemma F below was
> found to be executable with no prolate construction at all, and it is now
> **certified** — see §8. What remains of Lemma F is the *identification* of the
> conditioning space with the prolate span; the revised target is §9. This
> section is kept because its dependency analysis and stop criteria are what
> made §8 possible, and because its falsification criteria still apply to §9.


**Lemma F (certified prolate conditioning at one window).** Let `L = 4/5`,
`lambda = e^L`, `T = 400`, `N = 440`, and let `A_even` be the frozen `220 x 220`
even-sector enclosure of the compression of `R_T` already committed at `7515396`
(`phase4/certs/block400_even.json`). Let `P_4 subset V_N^{even}` be the span of

* the projections onto `V_N^{even}` of the `E`-images of the four zero-integral
  prolate combinations `f_k^{ev} = I_0 h_{4k,lambda} - I_{4k} h_{0,lambda}`,
  `k = 1,...,4`, where `I_n = int h_{n,lambda}` and
  `E(f)(u) = u^{1/2} sum_{1 <= n < lambda^2} f(nu)` (a sum of at most four terms
  at this `lambda`, with the integral split at `u = lambda/n`), together with
* the pole direction `rho^{ev} = (u^{1/2} + u^{-1/2})|_window`, equivalently the
  linear functional `l(f) = M f(1/2)`.

Then
```
    mu_lambda := inf { Q(f)/||f||^2 : f admissible, f perp P_4 }   >=   1/2 .
```

**Why this is the target.** It is the one statement for which (i) the paper
supplies the prediction (Fact 5.3: `0.27 <= mu_lambda <= 0.78` for
`lambda^2 <= 50`, `lambda mu_lambda in [1.6, 2.5]`, so `~0.78` at
`lambda^2 ~ 5`), (ii) the paper explicitly asks for a certified version (§8.8
item 4, "sharpening the conditioning of Theorem 8.14 to the `lambda^2-1`,
`lambda^2-2` vectors of Fact 5.3, and understanding the law `mu_lambda ~ 2/lambda`";
and item 3, "a numerical audit of the three bounds (73) at a modest cutoff …
with rigorous error bounds"), and (iii) **this repository can do it with no new
assembly.** Removing five directions out of 220 should lift the certified
constant from `1.43e-17` to order `1/2` — sixteen orders of magnitude — which
is the only mechanism in sight that is not doubly-exponentially capped.

**Why it is not just another decimal.** Theorem 8.14 obtains a coercivity
`d = (log sigma)/2` on a conditioned space only under the open semilocal
hypothesis (H2'). Lemma F would be an **unconditional** instance of the same
input at one window, certified. With it, Lemma 8.16 applies with
`V = P_4 (+) C rho`, `d = 1/2`, and whole-window positivity at `L = 4/5`
collapses to the sign of a `5 x 5` Schur complement — Proposition 8.18's exact
shape, with none of its hypotheses, at one `lambda`.

**Dependencies, each already discharged or in hand.**

1. The prolate functions `h_{n,lambda}` and the exact rungs: **done today**,
   validated against the paper's printed values (`work/prolate_rungs.py`). The
   Legendre coefficients `d^{(n)}` are what the recurrence needs; no library
   routine, no `pro_ang1`.
2. The `E`-map on the window: finite alias sum, at most four terms, no Gevrey
   repair needed — the repair in [11,16] exists to put `f` in `S_0^{ev}`, which
   matters for `K_lambda in Rad(QW)` exactly; here the conditioning vectors need
   only be *explicit*, not radical, because `mu_lambda` is a Rayleigh quotient
   over their orthogonal complement, not a radical statement.
3. **No Gram lower bound is required** — the Legendre basis is orthonormal
   exactly, so the five conditions are imposed by exact Gram–Schmidt inside
   `V_N^{even}` with Gram `= I`. This sidesteps the paper's own open item (§8.8
   item 2), which is a Gram bound beyond the Hermite regime.
4. Complement and coupling: **unchanged**. `d = beta - M g_N = 1.2116...`,
   `kappa = 1.3083e-29`, `g_N = 1.8864e-60` are properties of `V_N^perp`, not of
   the conditioning inside `V_N`, so the certified chain of `DERIVATION.md` §5
   applies verbatim — and `d = 1.21 > 1/2` already.
5. Certification: `block.verify_pd_preconditioned` on the `215 x 215`
   compression, at a shift of order `1/2` instead of `1e-17`. The radii
   (`<= 1.34e-48`) are 48 orders below that shift, so the ill-conditioning that
   forced Rump's criterion at `1e-17` is **absent** at this scale; even a naive
   interval Cholesky should succeed.

**Falsification / stop criteria, declared in advance.**

* **Falsified** if the certified `mu_lambda` is below `10^{-2}` — i.e. a
  minuscule eigenvalue survives the five conditions. That would contradict
  Fact 5.3 item 2 ("in the even sector the `lambda^2-1` prolate vectors alone
  already give the first non-minuscule eigenvalue, to 4 digits") at this
  `lambda`, and would mean the identification of the minuscule eigenspace with
  the prolate span fails where the dictionary is otherwise exact. Report as a
  negative result against the paper's Fact 5.3, not as a failure of the code.
* **Stop** if the four `E`-image vectors, after projection onto
  `V_N^{even}`, have a wrong-parity component above `10^{-6}` relative: the
  paper measures `10^{-14}` for small `k`, so anything larger means the alias
  sum or the split at `u = lambda/n` is mis-implemented, and the vectors must be
  validated against their Table 5 (`sin angle(p_1, eta_lambda(1)) = 3.4e-4` at
  `lambda^2 = 5`) before proceeding.
* **Stop** if `mu_lambda` lands in `[10^{-2}, 10^{-1}]` rather than near
  `0.78`: that is the odd sector's "transitional eigenvalue" behaviour
  (Fact 5.3 item 3) appearing in the even sector, and means the pole condition
  is doing more work than the prolate ones. Diagnose before claiming.
* **Not a target:** narrowing `[1.4330e-17, 1.674e-17]` further, or any
  `T >= 800` assembly. Both were measured today and both are capped
  (`5.0e-5` and `+16.4 %` respectively).

### Second-order items, in priority order

1. **Interval-enclose `delta_4`, `delta_6`** (Bouwkamp eigenproblem in Arb with
   a certified Legendre-coefficient tail, plus an enclosure of
   `phi_n(0) = sum d_k bar P_k(0)`). This is the only uncertified input to
   today's headline and it is a contained piece of work.
2. Finish and read out the `lambda^2 = 5` head-to-head now running (§6).
3. Correct `DERIVATION.md` §1's divergent archimedean series, strike
   `CLAIM_INDEX.md` §7's "are certified enclosures" clause for the discard
   ratios, update `note/README.md`'s check count `166 -> 176`, and record the
   two sharpened constants. **Proposed as a patch only; nothing was edited.**

---

## 8. ADDENDUM — §7's target was executable today, and it is now certified

While the `lambda^2 = 5` assembly ran, it became clear that Lemma F's
**coercivity** half needs no prolate construction and no `E`-map at all. The
mechanism: for **any** vectors `v_1,...,v_K` and any `c >= 0`,

```
        A + c * sum_i v_i v_i^T  >=  a I
   ==>  <A f, f> >= a ||f||^2    for every f with <v_i, f> = 0 for all i,
```

because the rank-`K` bump contributes exactly `0` on that subspace. So a single
Rump test on a rank-`K` perturbation of the **already-frozen** matrix certifies
conditioned coercivity, for whatever subspace the `v_i` happen to span. Their
provenance needs no rigour; only the interval test does. By Courant–Fischer the
same test certifies `lambda_{K+1}(A) >= a`.

### 8.1 Choosing the subspace (NUMERICAL, no rigour claimed)

`work/lowvecs.py`: block inverse iteration with a Rayleigh–Ritz step on the
**midpoint** of `phase4/certs/block400_even.json` (`mpmath`, `dps 50`,
Cholesky once, 14 iterations). Log `logs/lowvecs.log`. The five lowest
eigenvalues of the `220 x 220` even block at `L = 4/5`, `T = 400`:

| `j` | `lambda_j(A_even)` | residual `||A w - lam w||` |
|---|---|---|
| 1 | `1.43307240474e-17` | `2.8e-51` |
| 2 | `7.22796052796e-12` | `2.7e-51` |
| 3 | `5.59834111275e-7` | `2.4e-51` |
| 4 | `0.00286337178089` | `3.1e-36` |
| 5 | `0.617583280625` | `2.6e-3` (least converged) |

`lambda_1 = 1.43307240474e-17` reproduces `R_400(f_400) =
1.4330724047355850e-17` from `phase4/certs/cross.json` to twelve digits — a
**third** independent route to the same number (Rump floor, closed-form
Rayleigh quotient at the frozen minimizer, and now inverse iteration), and the
certified floor of §4.2 sits `5.0e-5` below it.

**Exactly four even eigenvalues lie below `10^{-2}`.** The paper's Table 3 at
`lambda^2 = 5` records exactly `4` for `QW_lambda^{ev}`, and Fact 4.2 item 1
predicts `lambda^2` of them per class. Independent agreement.

### 8.2 The eigenvalue ladder against the exact rungs — an unplanned second test

`work/rungs_ladder.py` computes the exact rungs at `lambda = e^{4/5}` from
Theorem 2.1 / Appendix A (not Fuchs):

```
    delta_4  = 2.2391127080494e-18     delta_8  = 2.4465397668381e-12
    delta_12 = 2.0478074953012e-7      delta_16 = 0.0018956321010246
```

| `j` | `lambda_j(A_even)/delta_{4j}(lambda)` | paper, `lambda^2 = 5` |
|---|---|---|
| 1 | `6.400` | `7.06` (Table 2) |
| 2 | `2.954` | "between `3.3` (`lambda^2=5`) and `6.2` (`lambda^2=50`)" (§4.2 item 2) |
| 3 | `2.734` | — |
| 4 | `1.511` | Table 4 range at `lambda^2=5`: `1.77 ... 7.06` |

All four sit in or just below the paper's range, consistently — which is
exactly right, since `lambda_j(A)` are eigenvalues of a compression of
`R_T <= Q` (running about `9 %` low, as `lambda_1` does), while the paper's
values are `N`-truncated upper bounds still decreasing in `N`.

And the paper's own sharp test, eq. (33):

```
    log(lambda_2/lambda_1)        = 13.131062
    log(delta_8/delta_4)  (exact) = 13.904106
    shortfall                     =  0.773044
```

against the paper's `0.771` at `lambda^2 = 5` and its reported range
`0.53 ... 0.83` over `5 <= lambda^2 <= 300` (§4.2 item 3, "with no trend in
`lambda`"). **Agreement to `0.002` in a quantity of size `13.9`, from a
different basis, a different normalization and a different codebase.** This was
not planned; it fell out of the eigenvector step.

### 8.3 The certificate

`work/condition.py`, log `logs/condition.log`, artifact
`certs/conditioned_coercivity_L45.json`. The five Ritz vectors are frozen as
**exact rationals with denominator `10^34`**; the bump constant is `c = 4`; the
shift is bisected 20 times in `[0, 2]`; the surviving shift is re-certified; the
Schur chain then runs with the **unchanged** `d = beta - M g_N = 1.2116136950581`
and `kappa = 1.3083437712853882590e-29`.

| conditions | certified `a` | Rump `||I-G||_inf` | whole-window conditioned `mu` | also certifies |
|---|---|---|---|---|
| `K = 4` | `0.61755371093750` | `2.6216530e-32` | **`0.6175537109375`** | `lambda_5(A_even) >= 0.6175537109` |
| `K = 5` | `0.85349273681641` | `2.1380206e-25` | **`0.8534927368164`** | `lambda_6(A_even) >= 0.8534927368` |

**The result, stated exactly.** There are four explicit frozen rational vectors
`v_1,...,v_4` in the `220`-mode even Legendre coordinate basis of
`extensions/positive-weil-certificate-2026-11-09` (recorded in
`certs/conditioned_coercivity_L45.json`) such that every admissible **even**
`f` supported in `[-4/5, 4/5]`, in the form domain, with
`<v_i, f> = 0` for `i = 1,...,4`, satisfies

```
    Q_zeta(f)  >=  0.6175537109375 * ||f||_2^2 ,
```

and with a fifth direction the constant rises to `0.8534927368164`. This is
**unconditional** and certified in Arb ball arithmetic throughout, with the
rank-`K` bump built from exact rationals.

Against `1.4330004882800e-17` for the unconditioned window: removing four
directions out of `220` lifts the certified constant by **sixteen orders of
magnitude**. In the paper's scale-free form, `lambda mu_lambda = 2.2255409285 x
0.8534927368 = 1.8995`, against Fact 5.3's "`1.6 <= lambda mu_lambda <= 2.5`
throughout" and the law `mu_lambda ~ 2/lambda` (`2/lambda = 0.8987` here,
against our certified `0.8535` — within `5 %`).

### 8.4 Exactly what this does and does not establish

* It **is** an unconditional, certified instance of the *input* that
  Theorem 8.14 obtains only under the open semilocal hypothesis (H2'): a
  coercivity constant of order `1` on an explicit finite-codimension subspace.
  With it, Lemma 8.16 applies at this window with `d = 0.6175`, and
  whole-window positivity at `L = 4/5` reduces to the sign of a `4 x 4` Schur
  complement — Proposition 8.18's exact shape, with none of its hypotheses, at
  one `lambda`.
* It **is** a certified counterpart to Fact 5.3 items 1–2, whose values in the
  paper are numerical, `N`-truncated and uncertified. Ours is a lower bound,
  theirs an upper bound; the two are complementary, as in §4.4.
* It is **not** an identification of the conditioning space with the prolate
  span. The `v_i` are approximate eigenvectors chosen numerically; they are
  *not* the `E`-images of the zero-integral prolate combinations
  `f_k^{ev} = I_0 h_{4k,lambda} - I_{4k} h_{0,lambda}`. Fact 5.1's claim
  (`sin angle(p_k, eta_lambda(k)) = 3.4e-4` at `lambda^2 = 5`) is **untested
  here** and is now the remaining half of §7's target.
* It says nothing about any other window, and nothing about RH. The reduction
  to a `4 x 4` block is a reduction at one `lambda`; Proposition 8.17's
  implication needs the family.
* The constant `0.6176` is not sharp: `K = 4` gives `lambda_5(A_even)`, and
  `A_even` is a compression of `R_400 <= Q`, so the true conditioned minimum of
  `Q` is somewhat larger.

---

## 9. The single best next mathematical target — revised

### 9.1 The strategic finding this session produced

The certified chain has two halves, and they are now sharply unequal.

* **The block half is cheap and unconditional.** §8 lifts the certified
  constant on the finite block from `1.43e-17` to `0.6176` by removing four
  explicit directions, at a cost of one Rump test on an already-frozen matrix.
  No gate, no assembly, no hypothesis.
* **The complement half still forces the gate.** `d = beta - M g_N` with
  `beta = H(T) - A_L` requires `T > T_gate(L)`, and
  `T_gate(L) = 2 pi e^{A_L}(1 + o(1))` with `A_L = sum_{log n < 2L} 2 Lambda(n)/sqrt n ~ 4 e^L`
  by Chebyshev. So the whole apparatus — **including §8's conditioned
  result** — is **doubly exponentially capped in `L`**, and `R2` already says
  `A_L` is optimal for any *pointwise* comb envelope while `R1` closed the
  density escape (positive density, continuous Kronecker–Weyl).

So the binding obstruction is **not** the finite block, **not** its
conditioning, and **not** the exterior discard (which §4.5 shows costs `9.3 %`,
a bounded factor). It is that `d > 0` is obtained from a pointwise gate.
**Theorem 8.14 is exactly a gate-free route to `d`** — `d = (log sigma)/2`, with
no `T` and no `A_L` — and that is the paper's real contribution to this
project. It is conditional on (H2').

This reframes what the repository should want from the paper: not its
numerics, not its upper bounds, but **(H2')**.

### 9.2 Target — Lemma P (prolate identification), the remaining half of §7

**Lemma P.** At `lambda = e^{4/5}`, let `h_{n,lambda}` be the prolate functions
of eq. (12) with coefficients from the Bouwkamp eigenproblem (75) — already
computed and validated in `work/prolate_rungs.py` — and let
`f_k^{ev} = I_0 h_{4k,lambda} - I_{4k} h_{0,lambda}`, `I_n = int h_{n,lambda}`,
be the zero-integral combinations of eq. (38). Let
`p_k = Proj_{V_N^{even}} E(f_k^{ev})` with
`E(f)(u) = u^{1/2} sum_{1 <= n < lambda^2} f(nu)` (at most **four** terms at this
`lambda`, split at `u = lambda/n`), Gram–Schmidt in order `k = 1,...,4`. Then

```
    sin angle( span{p_1,...,p_4},  span{v_1,...,v_4} )  <=  10^{-2} ,
```

where `v_1,...,v_4` are the frozen rational vectors of
`certs/conditioned_coercivity_L45.json`.

**Why this, and not something else.** §8.3 certifies the *constant* but proves
only that *some* explicit 4-dimensional subspace works. What makes
Theorem 8.14 worth having, by the authors' own Remark 8.15(c), is that it
describes the conditioning vectors **a priori**. Lemma P supplies that
description at this window. With §8.3 and Lemma P together, the statement
becomes: *at `lambda = e^{4/5}`, `Q_zeta >= 0.6176 ||f||^2` on the orthogonal
complement of the `E`-images of the four zero-integral prolate combinations and
the pole direction* — Proposition 8.18's exact shape, instantiated,
unconditionally, with the block explicitly identified. That is the strongest
statement this machinery can reach at a single window, and it is reachable.

**Dependencies, and why each is already discharged.**

1. Prolate coefficients `d^{(n)}` and the exact rungs — **done and validated**
   today against the paper's own printed 30-digit `chi_0`, `chi_4`
   (`work/prolate_rungs.py`).
2. `I_n = int_R h_{n,lambda}`: exact, no quadrature, from
   `int_{-1}^{1} bar P_k = sqrt 2 delta_{k0}` (the paper's App. A item 2) —
   so `I_n` is a single coefficient times a constant.
3. `P_k(z)` at the arguments `nu/lambda` the alias sum needs: three-term
   recurrence, stable on `[-1,1]`, evaluated directly (App. A item 3). No
   interpolation, no `pro_ang1`.
4. The projection onto `V_N^{even}`: inner products
   `int_{lambda^{-1}}^{lambda} E(f)(u) phi_m(log u) d*u` by certified
   Gauss–Legendre on each piece between consecutive breakpoints `lambda/n`,
   reusing `block.py`'s panelled rule with its Bernstein-ellipse remainder.
   The integrand is a finite sum of Legendre polynomials in `nu/lambda` times
   `u^{1/2}` — analytic on each piece, so the existing remainder applies.
5. **No Gevrey repair is needed.** The repair of [11,16] exists to place `f` in
   `S_0^{ev}` so that `K_lambda in Rad(QW)` *exactly*. Lemma P is a statement
   about angles between subspaces, and §8.3's certificate holds for whatever
   subspace the vectors span — so the conditioning vectors need only be
   explicit, not radical. **No Gram lower bound is needed either**: the Legendre
   basis is orthonormal exactly, which sidesteps the paper's open §8.8 item 2.
6. Certification: the subspace angle is a `4 x 4` SVD of `V^T P` with both
   bases interval-orthonormalised — small, and at a scale (`1e-2`) forty
   orders above the matrix radii.

**Falsification / stop criteria, declared in advance.**

* **Falsified** if `sin angle > 10^{-1}`. Then Fact 5.1 fails at this `lambda`
  in this normalization, and since §4.4, §8.1 and §8.2 have already confirmed
  the *eigenvalue* half of the dictionary to three independent tests, the
  failure would be located squarely in the **eigenvector** half — a genuine
  negative result about the paper's Fact 5.1, to be reported as such and not as
  a code failure.
* **Stop before certifying** if the wrong-parity component of any
  `E(f_k^{ev})` exceeds `10^{-6}` relative. The paper measures `10^{-14}` for
  small `k`, so anything larger means the alias sum or the split at
  `u = lambda/n` is mis-implemented. Validate against their Table 5
  (`sin angle(p_1, eta_lambda(1)) = 3.4e-4` at `lambda^2 = 5`) first.
* **Stop** if the certified angle lands in `[10^{-2}, 10^{-1}]`: that is the
  plunge-region behaviour the paper reports only near `k ~ lambda^2`, not at
  `k = 1`, and would need diagnosis before any claim.
* **Not targets.** Narrowing `[1.4330e-17, 1.674e-17]`; any `T >= 800`
  assembly; any window with `L > (log 5)/2` (where `A_L` jumps by exactly
  `2 log 5/sqrt 5` and `T_gate` goes from `119.09` to `502.39`). All three were
  measured and all three are capped.

### 9.3 Second-order items, in priority order

1. **Interval-enclose `delta_{4j}`** (`j = 1,...,4`) — Bouwkamp eigenproblem in
   Arb with a certified Legendre-coefficient tail and an enclosure of
   `phi_n(0) = sum_k d_k bar P_k(0)`. This is the only uncertified input to
   §4.4, §8.2 and §8.3's comparisons.
2. Read out the `lambda^2 = 5` head-to-head (§6) when the assembly lands.
3. Propose, as a patch and nothing more: `DERIVATION.md` §1's divergent
   archimedean series; `CLAIM_INDEX.md` §7's "are certified enclosures" clause
   on the discard ratios; `note/README.md`'s check count `166 -> 176`; and the
   three new constants (`1.4330004882800e-17`, `1.3857e-14`, `0.6175537109375`).

---

## 10. CORRECTIONS AND EXTENSION after external review (same day)

Four points were raised. **Three are accepted and one is a mis-framing of the
mechanism.** One of the accepted points (parity scope) led to a strictly
stronger result, below.

### 10.1 Parity scope — ACCEPTED; the claim was even-sector, and is now both

**The reviewer is right that the scope had to be stated, and right about the
min–max consequence.** The computation was always even-sector only
(`work/condition.py` reads `phase4/certs/block400_even.json` and uses
`d_even = beta - M g_N`), and §8.3 of this report does say "every admissible
**even** `f`". But **the proposed patch text (`REPORT.md` §12.2) omitted the
parity label**, and so did the summary table. That is a genuine defect and is
corrected in the regenerated patch.

The reviewer's min–max argument is correct and sharp: the odd sector also
carries minuscule eigenvalues, so four conditions on the **combined** space
cannot leave `0.61755`. Quantified today — the odd block was frozen too, so
this cost one run:

**Odd sector, frozen `phase4/certs/block400_odd.json`, inverse iteration
(`work/lowvecs_odd.py`, `logs/lowvecs_odd.log`):**

| `j` | `lambda_j(A_odd)` | residual |
|---|---|---|
| 1 | `1.38574029573e-14` | `5.2e-51` |
| 2 | `2.42537508217e-9` | `2.4e-51` |
| 3 | `6.88937205273e-5` | `5.0e-51` |
| 4 | `0.106807075206` | `1.1e-13` |
| 5 | `0.785079621892` | `1.4e-2` |

`lambda_1 = 1.38574029573e-14` confirms §4.2's certified `1.3857e-14` is tight.
**Exactly three** odd eigenvalues lie below `10^{-2}` — the paper's Table 3 at
`lambda^2 = 5` records exactly `3` for `QW_lambda^{odd}`, against `4` for the
even sector. Both counts now match independently.

**Odd certificates** (`work/condition_odd.py`, `logs/condition_odd.log`,
`certs/conditioned_coercivity_L45_odd.json`), same bump mechanism, with the
odd-sector `d = beta - M g_N - 2|p_perp|^2` (the pole term **carried**, not
dropped):

| conditions | certified `a` | Rump `||I-G||_inf` | whole-window odd `mu` |
|---|---|---|---|
| `K = 3` | `0.10680580139160` | `2.1597066e-37` | **`0.1068058013916016`** |
| `K = 4` | `0.78413009643555` | `9.2955279e-29` | **`0.7841300964355469`** |

**Corrected statement, both sectors.** With the two sets of frozen rational
directions together:

| conditions | certified lower bound on `Q_zeta` over **both** sectors |
|---|---|
| `4` even `+ 3` odd `= 7` | `Q_zeta(f) >= 0.1068058013916016 ||f||^2` |
| `4` even `+ 4` odd `= 8` | `Q_zeta(f) >= 0.6175537109375 ||f||^2` |

(the minimum of the two sector constants, since the form decouples by parity
and arbitrary complex `f` splits as `Q(u) + Q(v)` with each part splitting
again by parity).

**This reproduces the paper's Fact 5.3 item 3, which is a subtle prediction.**
At `lambda^2 = 5` the counts are `lambda^2 - 1 = 4` even and `lambda^2 - 2 = 3`
odd prolate conditions; Fact 5.3 item 3 then says the odd sector is left with
"one transitional eigenvalue `eps^odd_lambda(lambda^2-1) ~ 0.08-0.10`, the
partner of the plunge rung", which "the pole condition removes, multiplying the
minimum by a factor between `2.1` (`lambda^2 = 100, 150`) and `9`
(`lambda^2 = 5`)". We find that transitional eigenvalue at
**`0.1068058013916016`**, certified, and a fourth odd condition multiplies the
minimum by **`7.342`** — inside `[2.1, 9]` and near the `lambda^2 = 5` end, as
it should be at `lambda^2 = 4.953`. In the scale-free form,
`lambda mu_lambda = 1.374` with 8 conditions (Fact 5.3: `1.6 <= lambda mu_lambda <= 2.5`;
ours is below because `A` is a compression of `R_T <= Q`, see §10.5).

### 10.2 "The lower half is untested" — ACCEPTED, rephrased

The reviewer is right that this was too strong as rhetoric. The paper does
discuss convergence in `N` explicitly (its §4.1 "Convergence in `N`", its
Fact 4.1), and states that its values overestimate the limiting eigenvalues.
The precise statement, which is what §4.4 and §8 actually rest on:

> The paper's computations provide **numerical evidence** for both halves of
> Hypothesis 3.4, with the limitation — which the authors state — that an
> eigenvalue of a compression bounds the limiting eigenvalue only from above.
> They therefore **do not certify** the infinite-dimensional lower bound, and no
> computation of that kind can.

Every occurrence of "untested" is replaced by that formulation in the
regenerated patch and above. The substance of §4.4 and §8 is unchanged: our
certificates supply the side their method cannot.

### 10.3 "First rigorous two-sided value" — ACCEPTED, withdrawn

Dropped. I did not check prior work, and the reviewer is right that an existing
certified Weil bound combined with a certified prolate rung already yields a
ratio enclosure in principle. The defensible claim, and the only one now made:

> **an explicit certified enclosure** `6.3998587 <= eps_lambda(1)/delta_4(lambda) <= 7.4722703`
> at `lambda = e^{4/5}`, with a reproducible certificate, subject to `delta_4`
> being enclosed rather than high-precision.

No priority is claimed, consistently with the extension's own README.

### 10.4 The `lambda^2 = 5` comparison — ACCEPTED; my falsification criterion was WRONG

The reviewer is right and this is the one place where I stated something
invalid. I wrote (§6) that a certified lower bound at
`lambda_-^2 = 4.999999999999999927` exceeding their upper bound at
`lambda^2 = 5` would mean one of the two computations is wrong. **That does not
follow.** By Proposition 8.17, `eps_lambda` is non-increasing in `lambda`, so
`eps(lambda_-) >= eps(sqrt 5)`: the smaller window legitimately has the larger
minimum. My hedge ("the transport over `1.3e-17` in `L` is nil") rests on a
derivative estimate `d(-log eps)/dL ~ 8 pi e^{2L} - 9 ~ 115` read off the Fuchs
exponent — a heuristic, not a certified parameter-transfer bound.

Two further points of the reviewer's, both correct: matched parameters are what
the comparison needs, and their `9.1663e-18` is a **rounded printed number**
(reconstructed from `-log eps = 39.231` and from `7.06 x delta_4`), not an
outward enclosure, so it cannot be used as a bound at all.

**Matched parameters are not reachable with this code.** `L* = (log 5)/2` is
irrational, and `block.Symbol`, `baseline.A_L`, `baseline.comb_terms`,
`baseline.M_sup` and `block.assemble` are all parameterized by a *rational*
`L_num/L_den`. Running at exactly `L*` requires re-parameterizing the certified
assembly path to take an `arb`, which is a change to sign-critical code and is
not something to do quickly. (Worth noting for later: `L*` itself *is* an
admissible cheap window, because `log n < 2L*` excludes `n = 5` strictly, so
`A_{L*} = 2.941973525...` and `T_gate = 119.09` — the `A_L` jump happens for
`L > L*`, not at it. So the matched run would still cost one assembly, not ten
hours, once the code takes an irrational `L`.)

**Corrected falsification criterion for the run now in flight.** What is
rigorous is monotonicity plus our own certified upper bound at the smaller
window: since `4.95303 < lambda_-^2`,

```
    mu(lambda_-)  <=  eps(lambda_-)  <=  eps(e^{4/5})  <=  1.67312553407068186807726837736e-17 ,
```

so **`mu(lambda_-) > 1.6731255341e-17` would be a genuine contradiction** —
against Proposition 8.17 together with `phase5/certs/upper_highN.json`. That is
weaker than what I claimed, and it is what the run can actually test.

Everything else the run gives is **consistency**, not falsification: the
comparison with `~9.17e-18` is informative only modulo the uncertified
transport above. And one thing it gives rigorously in the useful direction:
by monotonicity a certified **upper** bound at `lambda_-` would upper-bound
`eps(lambda)` for every `lambda >= lambda_-`, including `lambda^2 = 5` — so the
cheap quadrature-free upper route is the piece worth adding there.

### 10.5 One point I do not accept: the angle is not needed at all

The reviewer writes that "a small `L^2` angle alone does not automatically
transfer the coercivity bound: the form errors and mixed terms also need
control", and that Lemma P's endpoint should be coercivity on the complement of
the explicit prolate space rather than an angle.

**The second half is right and is now the target (§11). The first half
describes a route I am not taking.** The bump mechanism of §8 needs no angle
and no perturbation theory:

```
      A + c * sum_i p_i p_i^T  >=  a I       (one Rump test)
  ==> <A f, f> >= a ||f||^2  for every f with <p_i, f> = 0,
```

for **any** explicit vectors `p_i` whatsoever — in particular for the prolate
`E`-images themselves. So the coercivity statement on the complement of
`V_prol` is obtained by running exactly the certified test of §8.3 with
`v_i := p_i`, with no transfer step, no form-error budget for the substitution
and no mixed terms beyond the `kappa` already in the chain. The subspace angle
against the numerically-chosen `v_i`, and against the paper's Table 5
(`sin angle(p_1, eta_lambda(1)) = 3.4e-4` at `lambda^2 = 5`), becomes a
**diagnostic** that the `E`-map was implemented correctly — not a step in the
proof. That is a strictly cleaner route than the one the reviewer's caution
addresses, and it is why §11 states the endpoint as the coercivity inequality
directly.

The reviewer's framing point stands, though, and I have adopted it: the
deliverable is `f perp V_prol => Q(f) >= c_prol ||f||^2` with `c_prol > 0`
explicit, not an angle.

### 10.6 Accepted without reservation

* "Removing those directions changes the minimum over the remaining space; it
  does not raise the minimum of the original form." Correct, and §8.4 says so.
  `lambda_Q(4/5)` remains in `[1.4330004882800e-17, 1.6731255341e-17]`;
  nothing in §8 or §10 moves it.
* Theorem 8.14 uses **smoothed semilocal pair vectors** and assumes (H2').
  Identifying our subspace with ordinary prolate `E`-images gives a concrete
  local connection and leaves (H2'), (R1), (R2) untouched. Stated in §8.4 and
  kept.
* The prolate denominator must be enclosed, not merely high-precision. It is
  §9.3 item 1 and remains the gating item for every ratio claim.

### 10.7 A further corroboration, with its discrepancy explained

The odd-sector ladder against the exact odd rungs at `lambda = e^{4/5}`
(`delta_6 = 3.4175764349346e-15`, `delta_10 = 9.3660403962553e-10`,
`delta_14 = 2.6157962275826e-5`; the `lambda^2 = 5` values
`-log delta_6 = 33.832643944` and `-log delta_10 = 21.261448843` reproduce the
paper's printed `33.833` and `21.261`):

| `j` | `lambda_j(A_odd)/delta_{4j+2}` | paper, `lambda^2 = 5` |
|---|---|---|
| 1 | `4.055` | `4.84` (Table 4, range `2.92-4.84`) |
| 2 | `2.590` | — |
| 3 | `2.634` | — |

The cross-sector test is where a discrepancy shows and it is fully accounted
for:

```
    log(lambda_1(A_odd)/lambda_1(A_even)) = 6.8742
    exact log(delta_6/delta_4)            = 7.3306
    shortfall                             = 0.4564
```

against the paper's `0.378` at `lambda^2 = 5` and its reported range
`0.29-0.39`. The excess is **exactly** the sector-dependent `R_T` deficit:
`lambda_1(A_even)` is `0.9065` of the paper's transported `eps_lambda(1)` and
`lambda_1(A_odd)` is `0.8378` of `eps^odd_lambda(1)`, and
`log(0.9065/0.8378) = 0.0789`, giving `0.378 + 0.079 = 0.457` against the
observed `0.4564`. So our shortfall is the paper's shortfall plus a known,
quantified artifact of comparing compressions of `R_T` rather than eigenvalues
of `Q` — not a disagreement.

---

## 11. Target, restated as the reviewer asks — and in the reviewer's order

### 11.1 Lemma P' (coercivity off the explicit prolate space)

Supersedes §9.2. Let `lambda = e^{4/5}`. In **each parity sector separately**,
let `V_prol` be defined explicitly as follows: with `h_{n,lambda}` the prolate
functions of eq. (12) (Bouwkamp coefficients from (75), validated in
`work/prolate_rungs.py`), `I_n = int_R h_{n,lambda}` exact from
`int_{-1}^1 bar P_k = sqrt 2 delta_{k0}`, and

```
    f_k^{ev}  = I_0 h_{4k,lambda}   - I_{4k}   h_{0,lambda} ,   k = 1,...,4
    f_k^{odd} = I_2 h_{4k+2,lambda} - I_{4k+2} h_{2,lambda} ,   k = 1,2,3
    E(f)(u)   = u^{1/2} sum_{1 <= n < lambda^2} f(nu)      (four terms here)
```

put `p_k = Proj_{V_N^{sector}} E(f_k^{sector})`, Gram-Schmidt in order
`k = 1, 2, ...`, together with the pole direction
`rho^{ev/odd} = (u^{1/2} +/- u^{-1/2})|_window`. **Then**

```
    f perp V_prol   ==>   Q_zeta(f) >= c_prol ||f||^2 ,   c_prol > 0 explicit,
```

certified, in each sector, with `c_prol` expected near `0.62` (even, 4
conditions) and near `0.78` (odd, 4 conditions including the pole direction).

**How it is proved — no angle, no transfer.** The `p_k` are frozen as exact
rationals and fed straight into the certified test of §8.3: one Rump test on
`A + c sum_k p_k p_k^T - a I` in the frozen sector block, then the unchanged
Schur chain (`d`, `kappa` untouched). The implication is exact because the
rank-`K` bump vanishes on `{f : <p_k, f> = 0}`. The subspace angle against the
`v_i` of §8.3 and against the paper's Table 5 is computed only as a
**diagnostic** that the `E`-map is right.

**Why it is worth doing even though §8/§10.1 already give the constants.** What
Theorem 8.14 contributes, by the authors' own Remark 8.15(c), is the *a priori
description* of the conditioning vectors. §8.3 and §10.1 certify constants on
subspaces found by inverse iteration — existence, not description. Lemma P'
replaces "some explicit 4-dimensional subspace" by "the `E`-images of the
zero-integral prolate combinations and the pole direction", which is a
statement that can be *written down at any `lambda`* and is therefore the only
part of this machinery with a chance of being uniform.

**Falsification / stop criteria.**

* **Falsified** if no `a >= 10^{-2}` certifies with `v_i := p_k`. Since §4.4,
  §8.1, §8.2, §10.1 and §10.7 have now confirmed the *eigenvalue* half of the
  dictionary by five independent tests, such a failure would be located in the
  **eigenvector** half — a negative result about Fact 5.1 at this `lambda`, to
  be reported as such.
* **Stop before certifying** if the wrong-parity component of any
  `E(f_k^{sector})` exceeds `10^{-6}` relative (the paper measures `10^{-14}`
  for small `k`), or if the diagnostic angle to `eta_lambda(1)` disagrees with
  their Table 5 value `3.4e-4` at `lambda^2 = 5` by more than an order of
  magnitude. Either means the alias sum or the split at `u = lambda/n` is
  mis-implemented.
* **Not a target.** Narrowing `[1.4330e-17, 1.674e-17]`; any `T >= 800`
  assembly; any `L > (log 5)/2`.

### 11.2 The order of work, adopted from the review

1. **Finish the matched-parameter comparison.** Concretely: re-parameterize
   `block.Symbol`, `baseline.{comb_terms, A_L, M_sup}` and `block.assemble` to
   accept an `arb` window (or an exact algebraic/symbolic `L`) so that
   `L* = (log 5)/2` can be used exactly. This is a change to sign-critical
   code and needs its own validation pass: re-run the full 176-check sequence
   at `L = 4/5` through the new signature and require bit-identical
   certificates before using it at `L*`. Note `A_{L*} = 2.9419735252...` and
   `T_gate(L*) = 119.09`, so the matched run costs one assembly.
   Until then, the `lambda_-` run stands as a bracket at `lambda_-`, with the
   corrected criterion of §10.4 — and a quadrature-free **upper** bound at
   `lambda_-` should be added, since by monotonicity it bounds `eps(lambda)`
   for every `lambda >= lambda_-`, including `lambda^2 = 5`, rigorously.
2. **Verify the sector and denominator certificates.** Sector: done in §10.1
   (both sectors, with the counts `4` and `3` matching the paper's Table 3, and
   the combined-space statement at `7` and `8` conditions). Denominator:
   interval-enclose `delta_{4j}` and `delta_{4j+2}` for `j <= 4` — Bouwkamp
   eigenproblem in Arb with a certified Legendre-coefficient tail and an
   enclosure of `phi_n(0) = sum_k d_k bar P_k(0)`. **This is now the gating
   item**: every ratio claim in §4.4, §8.2 and §10.7 is conditional on it.
3. **Prove positivity on the complement of the explicit prolate space** —
   Lemma P', §11.1.

### 11.3 What this tells us about the larger-window question

The reviewer's closing question is whether the structural description can
replace the expensive complement estimate at larger windows. The honest answer
from §9.1 is: **not by itself.** The complement constant `d = beta - M g_N`
comes from the gate `beta = H(T) - A_L`, so it needs `T > T_gate(L) = 2 pi e^{A_L}(1+o(1))`
with `A_L ~ 4 e^L` — doubly exponential in `L`, and `R2` says `A_L` is optimal
for any pointwise comb envelope while `R1` closed the density escape. Lemma P'
describes the *block*; it does not touch `d`.

What would replace the expensive complement estimate is exactly
**Theorem 8.14**, which delivers `d = (log sigma)/2` with no `T` and no `A_L` —
and is conditional on (H2'). So the value of Lemma P' is not that it scales:
it is that it identifies, at one window and unconditionally, the finite block
that Theorem 8.14 describes, so that the remaining gap between this repository
and the paper's programme is isolated to a single named hypothesis about the
*complement*, rather than spread across the block, the conditioning and the
complement together. That isolation is the result worth having.

---

## 12. The `lambda^2 = 5` run: THREE defects in my own scripts, then the result

The run completed (assembly `3678 s`, peak RSS `393.3 MiB`). Its first pass
produced almost nothing usable, for three reasons, **all mine**, and two of
them instances of the precision-ownership class this repository already
documented in `REPORT.md` §11.1 and `RETRACTIONS.md` R13(b). They are recorded
here in full because the repository ships a repaired path specifically to
prevent them and I bypassed it.

### 12.1 Defect A — `d` and `kappa` computed outside an owned scope

`work/run_block_lam5.py::sector_dk` recomputes `S_tail`, `g`, `pole_vec_tail`,
`kappa` and `d` inline. It is called for the **even** sector before any
`verify_pd_preconditioned` call has set `ctx.prec`, so it ran at whatever
precision `block.assemble` left behind. Result, from `logs/run_block_lam5.log`:

```
    -- even: d=[+/- 3.09e+230] kappa=nan
```

`d` a ball of radius `3e+230` centred on nothing and `kappa = nan`, so
`certify_final` could not close and the even sector reported `mu >= 0` —
**vacuous**. The odd sector, called *after* the even sector's Rump tests had set
`ctx.prec = 1400`, came out clean (`d = 1.21161369506`,
`kappa = 7.756529962e-29`). A precision bug that depends on call order, which is
exactly the signature §11.1 describes.

**This is the same corruption pattern as the disclosed artifact defect in
`phase4/certs/finite_block_certificate_T400.json`** (even-sector
`d = [+/- 2.75e+183]`, `kappa = nan`, `final = null`). The repository's response
to that was `phase5/src/constants.py`: `constants(L_num, L_den, T, N, prec)`
computes every one of these inside `owned(prec)` and runs `require_tight` on
`beta`, `M`, `kappa`, `d_even`, `d_odd`. **Had I called it, `require_tight`
would have raised rather than writing a vacuous certificate.** I wrote a fresh
`sector_dk` instead, and reproduced the defect the repaired path exists to stop.

### 12.2 Defect B — a one-sided shift search

`push()` doubles upward from a start guess and breaks on the first failure, with
no downward search. The odd start guess was `8e-15`; the true
`lambda_min(A_odd)` is just below it, at `~7.9e-15`. So the odd sector certified
**nothing at all** — not because anything was wrong with the matrix, but because
the search could only go the wrong way. A bad guess should cost a few seconds,
not the whole sector.

### 12.3 Defect C — R13(b), verbatim, in my repair script

The repair (`work/finalize_lam5.py`) loads the checkpointed blocks and re-runs
the certification. Its first attempt failed to certify **any** shift, down to
`7.5e-30`. Cause: `constants()` exits its `owned(prec)` scope and restores
`ctx.prec` to the ambient `53` bits, and I then parsed the serialized matrix
entries with `arb(v)` at that precision. Measured:

```
    entry "[0.0795639938261656777682... +/- 5.12e-49]"
      parsed at   53 bits -> radius 7.18e-19
      parsed at 1400 bits -> radius 5.12e-49
    inflation factor 1.4e30
```

Entry radii of `7e-19` against a `lambda_min` of `8.3e-18` leave nothing to
certify. R13(b) states this in one sentence — *"an input constructed outside the
scope carries the ambient precision's radius; re-parsing a serialized ball is
exactly such an input"* — and I walked into it anyway. Fixed by setting
`ctx.prec` explicitly before the parse.

**Why this matters beyond my scripts.** Three for three, the sign-critical
failures today were precision ownership, not mathematics. The repository's
existing countermeasures (`owned`, `require_tight`, `phase5/test_prec.py`'s
"calling without `prec=` is a `TypeError`") are the right ones, and the lesson
is narrow and actionable: **any new driver must go through
`phase5/src/constants.py`, and any re-parse of a serialized ball must name its
precision.** That belongs in the extension's documentation; it is in the
regenerated patch.

### 12.4 The result, after repair

`work/finalize_lam5.py`, log `logs/finalize_lam5.log`, artifact
`certs/certificate_T400_lam5_repaired.json`. No reassembly: the blocks
checkpointed by the first run (`certs/block400_lam5_{even,odd}.json`, `220x220`,
max radius `1.36e-48`, no `nan`/`inf`) were reused. Constants from
`constants(prec=2000)`, all passing `require_tight`.

Window: `L~ = 80471895621705018/100000000000000000`, strictly below
`L* = (log 5)/2` by `7.3e-18`, with the prime comb `{2,3,4}` (`n = 5` excluded),
`A_L~ = 2.941973525223620455504` identical to `A_{4/5}`, `beta_400 =
1.211613695058064549472`, `T = 400 > T_gate = 119.09`. So
`lambda_-^2 = 4.999999999999999927`.

| sector | certified `a` | Rump `||I-G||_inf` | `mu` (lower endpoint) |
|---|---|---|---|
| even | `8.273438682556152343750000e-18` | `5.05672988450e-11` | **`8.273438682556152343750000e-18`** |
| odd | `7.885009765625000000000000e-15` | `4.01037156749e-16` | **`7.885009765625000000000000e-15`** |

The first run's even `a = 8.2734375e-18` was itself a **valid** pass (
`verify_pd_preconditioned` sets its own precision); the repaired bisection is
larger by a factor `1.0000001`. Only the chain to the whole window was lost.

### 12.5 What the run does and does not test

**The rigorous criterion of §10.4 — passed, no contradiction.** Monotonicity
(Prop. 8.17) plus our own certified upper bound at the smaller window:

```
    mu_even(lambda_-) = 8.2734e-18  <=  eps(lambda_-)  <=  eps(e^{4/5})  <=  1.674e-17   OK
    mu_odd(lambda_-)  = 7.8850e-15  <=  eps^odd(lambda_-)  <=  eps^odd(e^{4/5})           OK
```

both hold with room. Nothing is falsified.

**An unplanned and sharper consistency test.** The two certified lower bounds at
two windows should scale like the exact rungs if the dictionary's `lambda`
dependence is right:

| | `mu(lambda^2=5.000)/mu(lambda^2=4.953)` | exact rung ratio | rel. diff |
|---|---|---|---|
| even | `0.577351` | `delta_4(5)/delta_4(4.953) = 0.579867` | **`0.43 %`** |
| odd | `0.569027` | `delta_6(5)/delta_6(4.953) = 0.592859` | `4.0 %` |

The even agreement to `0.43 %` is a two-window test of the Fuchs exponent
structure `-log delta_4 = 4 pi lambda^2 - 9 log lambda - l_4` through two
*certified* bounds, with no appeal to the paper's eigenvalues at all. This is
the first test today that uses two windows rather than one.

**The matched-window ratios**, now at `lambda^2 = 5 - 7.3e-17`, against the
paper's own printed `delta_4(sqrt 5) = 1.2983880087183625604e-18` and
`delta_6(sqrt 5) = 2.02613988982e-15`:

```
    eps_lambda(1)/delta_4(lambda)      >=  6.3720849      (paper, N=960: 7.06)
    eps^odd_lambda(1)/delta_6(lambda)  >=  3.8916413      (paper, N=960: 4.84)
```

The window mismatch is now `7.3e-18` in `L`, i.e. a relative effect on `eps` of
about `8.4e-16` by the derivative estimate — still an **uncertified** transport,
but no longer a 1 % effect. The paper's `7.06` and `4.84` remain **rounded
printed numbers**, not outward enclosures, so they bound nothing; the comparison
is consistency only, in the direction §10.4 allows.

**What is still missing for a genuine head-to-head**, unchanged from §11.2
item 1: an `arb`-parameterized window so that `L*` can be used exactly, and an
outward enclosure of the paper's side. Plus the quadrature-free **upper** bound
at `lambda_-`, which by monotonicity would bound `eps(lambda)` for every
`lambda >= lambda_-` — including `lambda^2 = 5` — rigorously, and is the one
piece that would make the comparison bite from our side.

---

## 13. CERTIFICATE REPAIR after second review — the `a = mu` defect

Four points were raised. **All four are accepted.** The first is a genuine
artifact defect and is the reason this section exists.

### 13.1 The defect: every published `mu` was unreplayable

The reviewer read `phase2/src/block.py::certify_final` and observed that with
`a = mu` the third test `(a-mu)(d-mu) > kappa^2` must fail against a nonzero
coupling, so a Rump pass on the finite block alone does not complete the step.
Checked, and the diagnosis is right in a way that is worse than cosmetic.

What the drivers did: bisect `mu` upward over `[0, a]` for 200-400 steps,
accepting only values where `certify_final(...)["ok"]` is true. That converges
to the supremum `mu*`, which satisfies `(a-mu*)(d-mu*) = kappa^2`, i.e.

```
    a - mu*  =  kappa^2/d  (to leading order)  =  1.41e-58  at L = 4/5
                                               =  4.97e-57  at lambda_-
```

**The Schur step was therefore genuinely performed, at full precision.** The
failure was in writing it down: `mu` was serialized with `.str(25)`, whose
truncation radius is `~3e-46` — about `2.1e12` times **larger** than the margin
it has to express. The published `mu` is consequently a ball straddling `a`, and
on replay `bool(a > mu)` is **False**:

```
    a  = [1.43300048828000000000000000000e-17 +/- 3e-51]
    mu = [1.4330004882800000000000000000e-17  +/- 3.01e-46]      <- as published
    a > mu  False      (a-mu)(d-mu) > kappa^2  False      ok = False
```

So the artifacts as written did not support the claim. Measured, confirmed, and
repaired — not argued away.

### 13.2 The repair: `mu` as a short exact decimal, with the margin shown

`work/margin_audit.py` locates the true supremum by bisection, then publishes

```
    mu  =  a * (1 - 10^-12),  truncated DOWN to 15 significant decimal digits,
```

a uniform rule. It gives up `10^-12` of `a` relatively — about `10^30` times the
true margin `kappa^2/d ~ 10^-41` of `a` — so it certifies with enormous room
while remaining an exact short decimal that a reader can quote and a replay can
parse exactly. Artifact: `certs/schur_margins.json`.

| certificate | `a` | **published `mu`** | `a - mu` | `certify_final` margin |
|---|---|---|---|---|
| `L=4/5` even, unconditioned | `1.43300048828e-17` | `1.43300048827856e-17` | `1.44e-29` | `1.74472372088e-29` |
| `L=4/5` odd, unconditioned | `1.3857e-14` | `1.38569999999861e-14` | `1.39e-26` | `1.68414303613e-26` |
| `L=4/5` even, 4 conditions | `0.61755371093750` | `0.617553710936882` | `6.18e-13` | `3.67129070187e-13` |
| `L=4/5` even, 5 conditions | `0.85349273681641` | `0.853492736815556` | `8.54e-13` | `3.05835298339e-13` |
| `L=4/5` odd, 3 conditions | `0.10680580139160` | `0.106805801391493` | `1.07e-13` | `1.18214444622e-13` |
| `L=4/5` odd, 4 conditions | `0.78413009643555` | `0.784130096434765` | `7.85e-13` | `3.35574624919e-13` |
| `lambda_-` even | `8.27343868255615234375e-18` | `8.27343868254787e-18` | `8.28e-30` | `1.00350011147e-29` |
| `lambda_-` odd | `7.885009765625e-15` | `7.88500976561711e-15` | `7.89e-27` | `9.55963205401e-27` |

All eight pass all three strict signs, and in all eight **`certify_final` FAILS
at `mu = a`** (margin `-1.71176342386e-58` at `L = 4/5`,
`-6.01637570459e-57` at `lambda_-`) — so the margin is doing real work, exactly
as the reviewer said it must. The audit asserts both directions and exits
non-zero otherwise.

**An internal confirmation fell out.** The bisected `a - mu*` came to
`1.41280e-58` at `L = 4/5` and `4.96559e-57` at `lambda_-`, against
`kappa^2/d = 1.4127964e-58` and `4.9655891e-57`. The break-even of the Schur
test reproduces `kappa^2/d` to six digits in both cases — a check on the
formula, not just on the arithmetic.

**Every constant in this report is unchanged** except for the `10^-12` relative
haircut on the published `mu`. The bracket at `L = 4/5` reads
`1.43300048827856e-17 <= lambda_Q(4/5) <= 1.674e-17`, and the short quotable
form is `1.4330e-17 <= lambda_Q(4/5) <= 1.674e-17`.

### 13.3 `1.6731e-17` was not a valid upper bound — R13(a), again, and mine

Also accepted, and this one is embarrassing in the same specific way the
repository already retracted once. I wrote the short upper endpoint as
`1.6731e-17` in four places. The certified outward endpoint is
`1.67312553407068186807726837736e-17`, and

```
    1.6731e-17  >=  1.67312553407068...e-17    is FALSE.
```

A decimal rounded **down** cannot bound a quantity from above. This is R13(a)
verbatim. The repository's own artifact already carries the correct short form
in `phase5/certs/upper_highN.json:safe_decimal = "1.674e-17"`, enforced by
`note/consistency_check.py` U4-U5, and I failed to use it. All four occurrences
are corrected to `1.674e-17`, and every short bound now in this report has been
re-checked for direction in ball arithmetic:

| short form | role | dominates / is dominated |
|---|---|---|
| `1.674e-17` | upper | `>= 1.673125534...e-17` TRUE |
| `1.6731255341e-17` | upper, full precision | TRUE |
| `1.4330e-17` | lower | `<= 1.43300048827856e-17` TRUE |

### 13.4 Relabelled: the two-window scaling is a diagnostic

Accepted. §12.5 and the patch's §13.3 compared `mu(lambda^2=5.000)/mu(lambda^2=4.953)`
with `delta_4(5)/delta_4(4.953)` and called it "a two-window test of the
exponent structure through two certified bounds". That overstates it twice
over: the numerators are **certificate floors, not eigenvalues** (each sits at
about `90.6 %` of its `eps` in the even sector, and the two happen to sit at
similar fractions — which is an observation, not a theorem), and the
denominators are **high-precision rung values, not interval enclosures**. The
agreement (`0.43 %` even, `4.0 %` odd) is a **diagnostic** that the two runs
and the dictionary are mutually consistent. It is not a certified
eigenvalue-scaling result, and nothing is inferred from it.

### 13.5 Relabelled: the `lambda_-` ratio inequalities are numerical comparisons

Accepted. The displayed

```
    eps_lambda(1)/delta_4  >=  6.3720849 ,      eps^odd_lambda(1)/delta_6  >=  3.8916413
```

are **not** established as inequalities. Two reasons, both the reviewer's: the
denominators `delta_4(sqrt 5)`, `delta_6(sqrt 5)` are evaluated at `lambda^2 = 5`
while the numerators are certified at `lambda_-^2 = 5 - 7.3e-17`, and they are
high-precision values rather than outward enclosures. They are **numerical
comparisons** until the enclosures exist and the windows match. The same caveat
applies to the `L = 4/5` pair `[6.3998587, 7.4722703]` on the enclosure count
only — there the window does match exactly, so the single missing ingredient is
the interval enclosure of `delta_4`.

### 13.6 Order of work, as instructed

1. repair or document the `a - mu` margin — **done**, §13.1-13.2,
   `certs/schur_margins.json`;
2. run the existing whole-window verifier — §13.7;
3. the inexpensive quadrature-free upper-bound check — §14;
4. consolidate the manuscript and reproducibility package — the regenerated
   `notes/proposed_2026-11-12.patch` and `notes/REPRODUCE.md`.

The exact-window reparameterization and the prolate-denominator enclosure stay
**outside version 1**, as instructed; they are recorded as the follow-on items
in §11.2.

### 13.7 The whole-window verifier, re-run

`logs/fastverify2.log`. All seven verifiers plus the upper-bound
re-certification, in the pinned venv, after every worktree revert:

| verifier | exit |
|---|---|
| `verify.py`, `verify_phase2.py`, `verify_phase3.py`, `verify_audit.py` | 0 |
| `phase4/verify_t400.py` | 0 |
| `phase5/test_prec.py`, `note/consistency_check.py` | 0 |
| `phase5/src/upper_highN.py --verify` | 0 — re-certified from exact rationals: `[1.673125534070681868077e-17 +/- 2.69e-39]` |

**176 `[PASS]`, 0 `[FAIL]`, every exit code 0.** The worktree is byte-identical
to `7515396`; nothing I did to it leaked into the verification.

---

## 14. A genuine two-sided bracket at `lambda_-`

The inexpensive quadrature-free upper bound, run at `L~` with the same two-stage
structure as `phase5/src/upper_highN.py` and the same closed-form evaluator
`src/certify.py::Win`, with `L` parameterized (`work/upper_lam5.py`, log
`logs/upper_lam5.log`, artifact `certs/upper_lam5.json`). `N = 160` even cosine
modes, exploratory `mpmath` stage to locate the vector (`647 s`, certifies
nothing), then the frozen exact rationals enclosed in Arb inside
`owned(prec = 2920)`. **No quadrature anywhere.**

```
    [CERTIFIED]  Q(f_160)/||f_160||^2   = [9.47476951083049032849709e-18 +/- 8.10e-43]
    [CERTIFIED]  outward upper endpoint =  9.47476951083049032849709e-18
    safe short decimal (rounds UP)      =  9.475e-18         dominates: TRUE
```

### 14.1 At `lambda_-` itself — a two-sided bracket

```
    8.27343868254787e-18  <=  eps(lambda_-)  <=  9.475e-18 ,
        lambda_-^2 = 4.999999999999999927 ,  even sector
```

width `1.2013e-18`, ratio `1.1452`. This is the second window at which this
repository brackets the window Weil form from both sides.

### 14.2 The paper's value lies inside it

The paper reports `eps_lambda(1)` at `lambda^2 = 5`, `N = 960`, as
`-log eps = 39.231`, i.e. `9.1663e-18` (equivalently `7.06 x delta_4`).

```
    8.2734e-18   <   9.1663e-18   <   9.4748e-18
```

It sits at **`74.3 %`** of our certified bracket. This is the cleanest statement
of the day: an independently computed eigenvalue, from a different
normalization, a different basis (log-Fourier vs Legendre/cosine) and a
different codebase, falls inside a two-sided certificate at a window differing
by `7.3e-18` in `L`. It is a consistency statement, not an inequality between
certified objects — their figure is a reconstructed printed decimal, not an
enclosure, and the windows are not identical.

### 14.3 What transports upward, and what does not

By Prop. 8.17 (`eps_lambda` non-increasing in `lambda`, unconditional), the
upper bound **does** transport to every larger window:

```
    eps_lambda  <=  9.475e-18     for every lambda >= lambda_- ,
                                  INCLUDING lambda^2 = 5 exactly.
```

That is a rigorous upper bound at the paper's own window, with no parameter
transfer and no reliance on their numerics. Against their reconstructed
`9.1663e-18` it is `3.37 %` weaker — the right side to be on, since theirs is
an `N`-truncated value and ours is a certified enclosure at a slightly smaller
window.

**The lower bound does not transport upward.** `mu(lambda_-) = 8.2734e-18` says
nothing about `eps` at `lambda^2 = 5` or beyond; monotonicity runs the other
way. Getting a certified *lower* bound at exactly `lambda^2 = 5` still requires
the `arb`-window reparameterization of §11.2 item 1.

Summary of the two windows:

| window | certified lower | certified upper | width |
|---|---|---|---|
| `lambda^2 = 4.95303242439511480365` (`L = 4/5`) | `1.43300048827856e-17` | `1.674e-17` | `2.40e-18` |
| `lambda^2 = 4.999999999999999927` (`L~`) | `8.27343868254787e-18` | `9.475e-18` | `1.20e-18` |

### 14.4 Ratio comparisons, two-sided — still NUMERICAL

With `delta_4(sqrt 5) = 1.2983880087183625604e-18`:

```
    eps_lambda(1)/delta_4(lambda)  in  [ 6.3720849 , 7.2973329 ]        (paper: 7.06)
```

two-sided now, and containing `7.06`. **Labelled a numerical comparison**, per
§13.5: `delta_4` is high-precision rather than interval-enclosed, and it is
evaluated at `lambda^2 = 5` while the bracket is certified at
`lambda_-^2 = 5 - 7.3e-17`. The `L = 4/5` pair `[6.3998587, 7.4722703]` has the
matching window but the same missing enclosure.

A last diagnostic: the two certified **upper** bounds scale as
`U(lambda_-)/U(e^{4/5}) = 0.56629` against the rung ratio
`delta_4(5)/delta_4(4.953) = 0.57987` — `2.3 %`, the same order as the
floors' `0.43 %` of §13.4 and, like it, a **diagnostic only**.

---

## 15. The two-window ratio — a rigorous comparison of the true minima

Pointed out in review, verified here. With `eps_0 = eps(e^{4/5})` and
`eps_- = eps(lambda_-)`, pairing the lower endpoint at each window with the
upper endpoint at the other gives

```
    mu_- / U_0   <=   eps_- / eps_0   <=   U_- / mu_0 .
```

This uses **only certified endpoints**: no prolate denominator, no parameter
transfer, no monotonicity transport, nothing from the paper. `work/window_ratio.py`
computes it in exact rational arithmetic and cross-checks in Arb balls.

Taking all four inputs as the **exact decimals that are actually packaged** —
the published `mu` (which round down) and the published `safe_decimal` upper
bounds `1.674e-17` and `9.475e-18` (which round up):

```
    mu_-/U_0 = 8.27343868254787e-18 / 1.674e-17   = 0.494231701466...
    U_-/mu_0 = 9.475e-18 / 1.43300048827856e-17   = 0.661200053838...
```

so, outward-rounded,

```
    0.494  <=  eps(lambda_-) / eps(e^{4/5})  <=  0.662 ,
```

i.e. **the true spectral minimum falls by between `33.8 %` and `50.6 %`
between `lambda^2 = 4.95303242439511480365` and
`lambda^2 = 4.999999999999999927`.**

Consuming the full-precision endpoints instead of the safe decimals tightens
this slightly, to `[0.494490, 0.661184]`, with the same outward 3-decimal form.

### 15.1 Why this is stronger than §13.4

§13.4's `0.43 %` compared two **certificate floors** against two **uncertified
rung values**, and is a diagnostic. This is a statement about the **true
minima** of the two infinite-dimensional forms, with every input certified and
every division taken in the safe direction. It is the first quantitative
statement in this project relating the window Weil form at two different
windows.

### 15.2 What it does not establish

It does **not** establish the prolate scaling law. The law predicts
`delta_4(5)/delta_4(4.953) = 0.579867`, which lies inside the bracket — at
`51.2 %` of the way across it — so the two are **consistent**, and that is all.
Testing the law to `+/- 5 %` would need the ratio bracket narrower than
`0.058`, i.e. each window's own bracket ratio below about `1.05`; they are
currently `1.1676` at `L = 4/5` and `1.1452` at `lambda_-`. So the tightening
that would make this a test of the law is roughly a threefold narrowing at both
windows — a concrete target, and a measurement of how far away it is.

### 15.3 Two qualifications retained in the release

Both adopted verbatim from the review.

* **The comparison with `9.1663e-18` is numerical consistency only.** That
  figure is a **reconstructed finite-compression estimate** — the lowest
  eigenvalue of the paper's `N = 960` log-Fourier compression, recovered from a
  printed `-log eps = 39.231` — lying inside our bracket at a nearby window. It
  is **not** an independently known value of the infinite-dimensional minimum,
  and §14.2's phrase "an independently computed eigenvalue" is corrected to say
  so. Nothing in this report treats it as a bound.
* **At exactly `lambda^2 = 5` there is only the transported upper bound.**
  `eps_lambda <= 9.475e-18` holds for every `lambda >= lambda_-` by Prop. 8.17.
  The positive lower bound runs the other way: `eps_lambda >= 8.27343868254787e-18`
  for every `lambda <= lambda_-`, i.e. on the whole range of windows
  `0 < L <= 0.80471895621705018` — an interval, not a point, which is worth
  stating — but it **does not extend past that endpoint**, and in particular not
  to `lambda^2 = 5`.

---

## 16. RELEASE ACCEPTANCE

`work/acceptance.py`, log `logs/acceptance.log`. The patch is applied to the
worktree and **every release claim is re-verified from the packaged artifact
files**, never from values held in memory by the run that produced them. The
patch is reverted at the end and the worktree checked clean.

| stage | what is re-verified | result |
|---|---|---|
| **A** | the patch applies to `7515396` and touches exactly 6 files; all seven repository verifiers plus `upper_highN --verify` run on the **patched** tree | 11 checks pass, **176 `[PASS]`** |
| **B** | all 8 Schur margins re-certified from `certs/schur_margins.json`'s own `a`, `d`, `kappa`, `mu` strings — and `mu = a` must still **fail** | 8 pass |
| **C** | the 4 conditioned certificates re-certified by rebuilding the rank-`K` bump from the frozen rational vectors **stored in the artifacts** (denominator `10^34`, `c = 4`) and re-running Rump | 4 pass |
| **D** | the 2 `lambda_-` block certificates re-certified from the packaged `220x220` checkpoints at the stored `a` | 2 pass |
| **E** | the `lambda_-` upper bound re-certified from its stored exact rationals; stored endpoint and `safe_decimal` both dominate | 3 pass |
| **F** | the two-window ratio recomputed from the four packaged exact decimals; the prolate law shown inside the bracket | 6 pass |
| — | worktree restored to `7515396` | pass |

**ACCEPTANCE: PASS — 35 checks, 0 failed.** The acceptance script exits non-zero on
any failure, so it is usable as a release gate.

Note on `D`: the re-certified Rump residuals (`1.73e-8` even, `9.12e-11` odd)
differ from the original run's (`5.06e-11`, `4.01e-16`) because the
preconditioner is recomputed from the re-parsed blocks. Both are far below the
`< 1` threshold, so both certify; the residual is a property of the
preconditioner, not of the claim.

### 16.1 Release scope, stated once

Two rigorous window brackets:

```
    1.43300048827856e-17  <=  lambda_Q(4/5)      <=  1.674e-17
    8.27343868254787e-18  <=  lambda_Q(L~)       <=  9.475e-18
```

plus, on `0 < L <= 0.80471895621705018`, `lambda_Q(L) >= 8.27343868254787e-18`
by monotonicity; the parity results (`4` even and `3` odd eigenvalues below
`10^{-2}`, matching the paper's Table 3 at `lambda^2 = 5`); the conditioned
coercivity results per sector (`4` even directions give `>= 0.617553710936882`,
`4` even `+ 4` odd give the same over both sectors); the proved normalization
dictionary `lambda = e^L`; the rigorous two-window ratio
`0.494 <= eps_-/eps_0 <= 0.662`; explicit reproduction instructions; and a
documented account of the failure modes (`RETRACTIONS.md` R15, four defects,
all failures to certify rather than false certifications).

**Out of version 1**, as instructed: the irrational-window reparameterization
and the interval enclosure of `delta_{4j}`. Every ratio against `delta_{4j}`
is labelled a numerical comparison until the second of those lands.

---

## 17. THREE CORRECTIONS before freezing version 1

All three accepted. Point 3 was a genuine gap in the release gate; point 1 was an
overstatement that turned out to be **repairable into the stronger claim** in the
even sector; point 2 is a mislabelling.

### 17.1 Point 3 — the gate did not cover the two improved `L = 4/5` shifts. FIXED

The reviewer is right, and this was the serious one. Acceptance stages C and D
re-ran `verify_pd_preconditioned` for the four conditioned certificates and the
two `lambda_-` blocks, and stage B re-ran the **scalar** `certify_final` for all
eight. But `certify_final` *assumes* `A >= a I`; it does not check it. And stage
A's `phase4/verify_t400.py` tests the shifts in `phase4/certs/certificate_T400.json`,
which are the **committed** `1.42e-17` and `2.4e-17` — not the improved
`1.43300048828e-17` and `1.3857e-14` that the headline bracket rests on. So the
two matrix bounds underpinning the release's lower endpoint were verified once,
earlier in the session, and never by the gate.

**Added as acceptance stage B2**, from the frozen blocks:

```
    [PASS] L=4/5 even: A - 1.43300048828e-17 I is PD   ||I-G|| = 1.1201203e-12
    [PASS] L=4/5 odd : A - 1.3857e-14 I is PD          ||I-G|| = 5.6557966e-18
    [PASS] each improved shift exceeds the committed one
```

Acceptance is now **42 checks, 0 failed** (`logs/acceptance.log`). The stage
header in the script states explicitly that stage A does not cover these two.

### 17.2 Point 1 — "exactly" vs "at most". Corrected, and the even sector upgraded

The reviewer's reasoning is exactly right. What the conditioned certificates give
is an **upper** bound on the count: `eps_lambda(5) >= 0.6175537109` certifies at
most four even eigenvalues below `10^{-2}`, and `eps^odd_lambda(4) >= 0.1068058014`
at most three odd. The counts I had quoted as "exactly four and exactly three"
came from `lambda_j(A)`, the eigenvalues of a **compression of `R_T`**, and since
`Q >= R_T` those bound `eps_lambda(j)` from **neither** side for `j >= 2`. That
claim was unsupported.

**For the even sector this is now certified both ways** (`work/counts_even.py`,
`certs/counts_even_L45.json`). `Win(4,5).matrix(N)` is the Weil form on the first
`N` even cosine modes in an **orthonormal** basis, in closed form
(digamma / Lerch `Phi` / Hurwitz `zeta`) as Arb balls — a **compression of `Q`**,
so by min–max `lambda_j(M) >= eps_lambda(j)`: its eigenvalues are *upper* bounds.
Then for an explicit 4-dimensional `W` of frozen rationals,

```
    eps_lambda(4)  <=  lambda_4(M)  <=  lambda_max(M|_W) ,
```

and `lambda_max(M|_W) < c` is certified by the `4x4` test `c G - B > 0` with
`G = W^T W`, `B = W^T M W`. At `N = 160`, `prec = 2920` (entry radii `~1.9e-868`;
they underflow `float`, which is why an earlier log line printed `0.0e0`):

```
    eps_lambda(4)  <=  0.00311839891820  <  10^-2        CERTIFIED
```

so **at least** four. With at most four from the conditioned certificate:

> **EXACTLY four even eigenvalues of `Q` below `10^{-2}` at `L = 4/5`** —
> both bounds certified.

**For the odd sector only "at most three" is certified.** `Win` implements the
even cosine basis only (`src/certify.py` docstring: "Basis (EVEN sector)"), so
there is no odd compression of `Q` available and no "at least three". The release
says **at most three**, consistent with the paper's Table 3 entry of `3`, and the
observed `lambda_j(A_odd)` values are labelled separately as numerical.

Corrected statements, replacing every earlier "exactly four / exactly three":

| claim | status |
|---|---|
| at most 4 even eigenvalues of `Q` below `10^{-2}` | **CERTIFIED** (`eps_lambda(5) >= 0.6175537109`) |
| at least 4 even eigenvalues of `Q` below `10^{-2}` | **CERTIFIED** (`eps_lambda(4) <= 0.0031183989`) |
| hence **exactly 4** even | **CERTIFIED** |
| at most 3 odd eigenvalues of `Q` below `10^{-2}` | **CERTIFIED** (`eps^odd_lambda(4) >= 0.1068058014`) |
| at least 3 odd | **not certified** — no odd compression of `Q` here |
| `lambda_j(A)` for `j >= 2` as eigenvalue estimates | **NUMERICAL**, observed in a compression of `R_T`; bounds `eps_lambda(j)` from neither side |

A by-product worth recording, not pursued: the same `Win` compression gives
certified upper bounds `eps_lambda(1..5) <= 1.6731255e-17, 8.4227418e-12,
6.4517875e-7, 0.0031183989, 0.62402164`, and running the even conditioning at
`K = 1, 2, 3` would give matching certified lower bounds — bracketing the first
five even eigenvalues individually. That is a one-command follow-on and is left
to v1.1, per the instruction to freeze scope.

### 17.3 Point 2 — `0.579867` is a numerical rung-ratio reference

Accepted. Hypothesis 3.4 compares `eps_lambda(j)` with `delta_{4j}(lambda)` up to
constants `c <= C` and a factor `log lambda`. It does **not** assert that the
ratio of two Weil minima equals the ratio of their rungs; that would need the
proportionality factor to be the same at both windows, which is an extra
assumption the hypothesis does not make. So `delta_4(5)/delta_4(4.953) = 0.579867`
is a **numerical rung-ratio reference**, not a prediction of the paper, and §15.2's
phrasing ("the law predicts") is wrong.

What stands unchanged: the certified interval

```
    0.494  <=  eps(lambda_-)/eps(e^{4/5})  <=  0.662
```

is independent of any of this. The rung ratio happens to lie inside it. A future
`+/- 5 %` comparison would be testing the **stronger approximation** that the
proportionality factor is window-independent — not Hypothesis 3.4 itself — and
should be described that way.
