
---

## R15 — precision ownership: a rule for NEW drivers, after three repeats on 2026-11-12

Not a retraction: nothing published here is withdrawn. This records a rule,
because the defect of `REPORT.md` §11.1 and of `R13(a)/(b)` recurred **four
times in one session** in new driver code written outside the repaired path, and each
time it produced a silently vacuous or silently unprovable certificate rather
than an error.

| # | what happened | consequence |
|---|---|---|
| A | a fresh `sector_dk` recomputed `S_tail`, `g`, `pole_vec_tail`, `kappa`, `d` inline, called before any `verify_pd_preconditioned` had set `ctx.prec` | `d = [+/- 3.09e+230]`, `kappa = nan` in the **even** sector only (the odd sector, called after a Rump test, was clean) — so `certify_final` could not close and the sector reported a **vacuous** `mu >= 0`. Exactly the corruption pattern of the disclosed defect in `phase4/certs/finite_block_certificate_T400.json`. |
| B | a shift search that doubles upward from a start guess and breaks on first failure | a start guess `1.4%` above the true `lambda_min` certified **nothing**, with no error |
| C | serialized matrix entries re-parsed with `arb(v)` after `constants()` had exited its `owned(prec)` scope and restored `ctx.prec = 53` | every entry radius inflated from `5.12e-49` to `7.18e-19`, a factor `1.4e30`; against a `lambda_min` of `8.3e-18` **no shift certified at all**, down to `7.5e-30` |
| D | `mu` bisected to the supremum (`a - mu* = kappa^2/d ~ 1e-58`) and then serialized with `.str(25)`, whose truncation radius is `~3e-46` | the published `mu` is a ball **straddling** `a`, so on replay `bool(a > mu)` is False and `certify_final` refuses. The Schur step had been performed correctly at full precision; the **artifact** did not support the claim |

### The rule

1. **Every new driver computes its certificate constants through
   `phase5/src/constants.py::constants(L_num, L_den, T, N, prec=...)`**, never
   inline. That function runs inside `owned(prec)` and applies `require_tight`
   to `beta`, `M`, `kappa`, `d_even`, `d_odd`. In case A it would have **raised**
   instead of writing a vacuous record.
2. **Any re-parse of a serialized ball names its precision**, by setting
   `ctx.prec` immediately before the parse or by parsing inside `owned(prec)`.
   `constants()` restores the ambient precision on exit, so code after it is at
   `53` bits unless it says otherwise. This is `R13(b)` applied to matrix
   deserialization rather than to a scalar endpoint.
3. **A shift search is two-sided.** A one-sided ladder cannot distinguish "the
   matrix is not positive definite at this shift" from "my guess was too high",
   and the two have opposite meanings.
4. **A published constant must be serializable with room to spare.** A `mu`
   bisected to the Schur supremum has margin `kappa^2/d`, which is `~1e-41` of
   `a` relatively and cannot survive any print width. Publish
   `mu = a*(1 - 10^-12)` truncated DOWN to 15 significant digits: it keeps
   `10^30` times the needed margin, loses `10^-12` of the constant, and replays
   exactly. A certificate that cannot be re-verified from its own artifact is
   not a certificate.

### What is unaffected

Nothing certified in phases 1-6 changes. Defects A-C occurred only in new
investigation drivers (`~/scratch/rh-connes-continuation-2026-11-12/work/`), and
all three were failures to *certify*, never false certifications: in each case
`verify_pd_preconditioned` sets its own precision internally, so the shifts it
reported as PASS were valid. Defect A's `a` was a valid pass whose chain was
lost; B and C certified nothing; D's underlying certification was sound and only
its serialization was unusable, now repaired in `certs/schur_margins.json`. The repaired results are in
`certs/certificate_T400_lam5_repaired.json` of that directory.
