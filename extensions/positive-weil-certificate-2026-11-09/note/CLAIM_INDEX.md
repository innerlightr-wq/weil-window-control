# CLAIM_INDEX.md — every claim in the note, mapped to its artifact

Each row gives: the claim as it appears in
`weil_two_window_positivity.tex`, the artifact file that holds it, the field
inside that file, the working precision at which it was produced, the published
enclosure, and the verifier that re-derives it.

All paths are relative to `extensions/positive-weil-certificate-2026-11-09/`.
Printed decimals in the note are **summaries**; the enclosures below are
authoritative. A serialized enclosure is a ball that *contains* the quantity, so
it must be consumed through the correct endpoint: lower endpoint for a lower
bound, upper endpoint for an upper bound. (`note/consistency_check.py` check
`U1` enforces exactly this for the one place it matters.)

---

## 0. Version 2.0.0 additions (phase 7)

Every row below is replayed by `phase7/acceptance.py`. `a` is the Rump-certified
shift on the finite block; `mu` is the whole-window constant that
`certify_final` then certifies, and is **strictly below** `a` because the Schur
test requires `(a-mu)(d-mu) > kappa^2`. Published `mu` values are short exact
decimals at `a(1 - 10^-12)`, which keeps about `10^30` times the true margin
`kappa^2/d` and survives serialization; a `mu` at the supremum cannot be
serialized (see `RETRACTIONS.md` R15(d)).

| claim | file | field | prec | value | verifier |
|---|---|---|---|---|---|
| `Q(f) >= 1.43300048827856e-17 ‖f‖²`, `L=4/5` both sectors | `phase7/certs/schur_margins.json` | `L45_even_uncond.mu_quotable` | 2000 | `a = 1.43300048828e-17`, margin `1.74472372088e-29` | `acceptance.py` B2, B |
| odd sector, `L=4/5` | same | `L45_odd_uncond.mu_quotable` | 2000 | `1.38569999999861e-14`, margin `1.68414303613e-26` | `acceptance.py` B2, B |
| `Q(f) >= 8.27343868254787e-18 ‖f‖²` at `L~` | `phase7/certs/certificate_T400_lam5_repaired.json` + `schur_margins.json` | `sectors.even.a`, `lam5_even.mu_quotable` | 2000 | `a = 8.27343868255615234375e-18`, `‖I-G‖ = 5.05672988450e-11` | `acceptance.py` D, B |
| odd sector at `L~` | same | `lam5_odd.mu_quotable` | 2000 | `7.88500976561711e-15`, `‖I-G‖ = 4.01037156749e-16` | `acceptance.py` D, B |
| `lambda_Q(L~) <= 9.47476951083049032849709e-18` | `phase7/certs/upper_lam5.json` | `outward_upper_endpoint` | 2920 | safe decimal `9.475e-18`; closed form, **no quadrature** | `acceptance.py` E |
| conditioned, even, 4 directions | `phase7/certs/conditioned_coercivity_L45.json` | `results.4` | 1400 | `a = 0.61755371093750`, `mu = 0.617553710936882`, `‖I-G‖ = 2.6216530e-32` | `acceptance.py` C, B |
| conditioned, even, 5 directions | same | `results.5` | 1400 | `mu = 0.853492736815556`, `‖I-G‖ = 2.1380206e-25` | `acceptance.py` C, B |
| conditioned, odd, 3 directions | `phase7/certs/conditioned_coercivity_L45_odd.json` | `results.3` | 1400 | `mu = 0.106805801391493`, `‖I-G‖ = 2.1597066e-37` | `acceptance.py` C, B |
| conditioned, odd, 4 directions | same | `results.4` | 1400 | `mu = 0.784130096434765`, `‖I-G‖ = 9.2955279e-29` | `acceptance.py` C, B |
| conditioning vectors (explicit rationals, denominator `10^34`, bump `c = 4`) | both files | `vectors` | exact | `5` even, `6` odd vectors of length `220` | `acceptance.py` C |
| `eps_even(4) <= 0.00311839891820 < 0.003119 < 1e-2` | `phase7/certs/counts_even_L45.json` | `certified_lambda4_upper` | 2920 | `4x4` test `cG - B > 0`, `G = WᵀW` verified PD | `acceptance.py` G |
| hence **exactly four** even eigenvalues below `1e-2` at `L = 4/5` | `counts_even_L45.json` + `conditioned_coercivity_L45.json` | — | — | at most 4 from `mu_5 >= 0.6175`, at least 4 from `c < 1e-2` | `acceptance.py` G |
| odd count | — | — | — | **at most three** only; no odd compression of `Q` is available here | — |
| `0.494 <= lambda_Q(L~)/lambda_Q(4/5) <= 0.662` | the four packaged endpoints | — | exact | `mu_-/U_0 = 0.494231701466`, `U_-/mu_0 = 0.661200053838` | `acceptance.py` F |
| the exact prolate rungs used for comparison | `phase7/src/prolate_rungs.py` | — | 80–110 dps | **high-precision, NOT interval enclosures** | — |

> **Directed decimals.** The gate checks the direction of every short bound and
> explicitly **rejects** `1.6731e-17` as an upper form for `U_0` and `0.6176` as
> a lower form for `0.617553710936882`. Accepted short forms: `1.674e-17`,
> `9.475e-18`, `0.003119` (upper); `1.4330e-17`, `8.2734e-18`, `0.6175`
> (lower).

> **Numerical, not certified.** Every ratio against a prolate rung
> `delta_{4j}(lambda)`; the counts read off a compression of `R_T` (since
> `Q >= R_T`, those bound the eigenvalues of `Q` from neither side for `j >= 2`);
> and the placement of the `9.1663e-18` value reconstructed from
> arXiv:2611.13129v1's printed `-log eps` inside the `L~` bracket.

## 1. Theorem (i) — the whole-window lower bound

| claim | file | field | prec | enclosure | verifier |
|---|---|---|---|---|---|
| `Q(f) >= 1.41e-17 ‖f‖²`, both sectors — **superseded by §0** | `phase4/certs/certificate_T400.json` | `whole_window_mu` | 2000 | `[1.410000000000000000000000e-17 +/- 1e-46]` | `phase4/verify_t400.py`, `note/consistency_check.py` T1–T8 |
| even sector `a` | same | `sectors.even.a` | 2000 (emitted); PD check at 1400 | `[1.420000000000000000000000e-17 +/- 1e-46]` | `phase5/src/parity.py` replay |
| even sector margin | same | `sectors.even.margin` | 2000 | `[1.21161369506e-19 +/- 1.94e-31]` | `note/consistency_check.py` T6-even |
| odd sector `a` | same | `sectors.odd.a` | 2000 (emitted); PD check at 1400 | `[2.400000000000000000000000e-17 +/- 1e-46]` | `phase5/src/parity.py` replay |
| `beta_400 = H(400) − A_L` | same | `beta` | 2000 | `[1.211613695058064549471502 +/- 1.89e-26]` | `phase4/verify_t400.py` |
| `M = sup_{\|t\|<=T}\|Psi_L − beta\|` | same | `M` | 2000 | `[9.525770639507350587208369 +/- 9.37e-26]` | `phase4/verify_t400.py` |
| `kappa` (coupling) | same | `kappa` | 2000 | `[1.308343771285388259029018e-29 +/- 1.31e-54]` | `phase4/verify_t400.py` |
| `kappa²/d` (Schur penalty) | same | `sectors.*.kappa2_over_d` | 2000 | `[1.412796364751587723429843e-58 +/- 3.25e-83]` | `phase4/verify_t400.py` |
| `g_N` (complement) | `phase4/certs/certificate_T400.json` | `g` | 2000 | `[1.886442193224755508462446e-60 +/- 2.48e-85]` | `phase4/verify_t400.py` |
| `S_N(320)` | `phase4/certs/finite_block_certificate_T400.json` | `sectors.odd.S_tail` | 2000 | `[9.260051461963608020971092e-63 +/- 3.67e-88]` | `phase4/verify_t400.py` |
| `\|p_perp\|` (pole tail) | same | `sectors.*.p_perp` | 2000 | `[7.993372152149564201076486e-1282 +/- 2.23e-1307]` | `phase4/verify_t400.py` |
| Rump residual, even | `phase5/logs/parity.log` | — | 1400 | `[1.4973903e-14 +/- 3.61e-22]` | `phase5/src/parity.py` |
| Rump residual, odd | same | — | 1400 | `[4.3535838e-22 +/- 2.65e-31]` | `phase5/src/parity.py` |
| assembly parameters (`N`, panels, order, `b_e`, `rho`, remainder) | `phase4/certs/finite_block_certificate_T400.json` | `meta.*` | — | see file | `phase4/verify_t400.py` |
| max entry radius after folding the remainder | `phase2/certs/block_even.json` | all entries | 384 | `< 4.99e-60` at `T=200` | `phase2/verify_phase3.py` X2b |

> **Provenance note on `certificate_T400.json`.** This file was originally
> written ad hoc and had no committed producer. It now has one:
> `phase4/src/emit_certificate_T400.py` rebuilds every derived field through the
> owned-precision path of `phase5/src/constants.py`, re-runs the three strict
> Schur signs for both sectors, and refuses to write if either fails. The
> regenerated file matches the original byte for byte except for one serialised
> radius (`3e-46` instead of `1e-46` on the odd `a`), both valid enclosures.
> Nothing decisive ever rested on the file: `phase4/verify_t400.py` recomputes
> `beta`, `M`, `S_N`, `g`, `kappa` and `d` from `baseline` and reads only the
> proposed shift `a` and the proposed constant `mu` from it, then re-runs the
> preconditioned PD check and `certify_final` itself.

> **Known artifact defect, disclosed.** In
> `phase4/certs/finite_block_certificate_T400.json` the **even** sector's
> `d`, `kappa`, `g`, `S_tail` and `final` fields are corrupt
> (`d = [+/- 2.75e+183]`, `kappa = nan`, `final = null`). That record was
> written by the run that carried the precision-ownership leak described in
> `REPORT.md` §11.1. It is **retained unaltered for provenance** and is not
> used by any claim. The authoritative even-sector record is
> `phase4/certs/certificate_T400.json`, recomputed through the repaired path and
> replayed in `phase5/certs/parity.json`. `phase4/verify_t400.py` reads the
> authoritative file.

## 2. Theorem (ii) — the certified upper bound

> **STATUS CHANGE, 2026-11-10 (`RETRACTIONS.md` R14).** The time-domain route of
> `phase4/src/compact.py` omits the truncation error of its outer panelled
> Gauss–Legendre rule. Every value it produced is **withdrawn as a certificate**
> and retained as a **measurement**. The certified endpoint now comes from the
> phase-1 closed-form route, which performs no quadrature.

| claim | file | field | prec | enclosure | verifier |
|---|---|---|---|---|---|
| **`Q(f_160)/‖f_160‖²`, certified** | `phase5/certs/upper_highN.json` | `enclosure` | 2920 | `[1.67312553407068186807726837736e-17 +/- 4.41e-47]` | `phase5/src/upper_highN.py` (re-certifies from the stored rationals); `note/consistency_check.py` U3–U6d |
| the 28-mode replay, retained as an independent check of the same route | `phase5/certs/upper_fallback.json` | `enclosure` | 1600 (also at 1200) | `[2.112992611721490779453e-17 +/- 1.9e-354]` | `phase5/src/replay_upper.py` (11 checks); `note/consistency_check.py` U1 |
| outward upper endpoint `U` | `phase5/certs/upper_highN.json` | `outward_upper_endpoint` | 2920 | `1.67312553407068186807726837736e-17` | same |
| safe quoted decimal | same | `safe_decimal` | — | `1.674e-17` | `note/consistency_check.py` U4–U5 |
| frozen witness `f_160` (160 even cosine modes) | `phase5/certs/upper_highN.json` | `x` | exact | rationals, denominator `10^180` | `phase5/src/upper_highN.py` |
| route is quadrature-free and independent of `compact.py` | `phase5/certs/upper_fallback.json` | `route`, `independent_of` | — | checked by source inspection | `phase5/src/replay_upper.py` U1–U2 |
| `q(f_400)` — **MEASUREMENT, not a certificate** | `phase4/certs/cross.json` | `Q_f400` | 1100 | `[1.90138291114546440e-17 +/- 6.03e-35]`, radius **incomplete** | — |
| `q(f_200)` — **MEASUREMENT, not a certificate** | `phase4/certs/q_f200.json` | `q` | 1100 | `[3.64859175290407702e-17 +/- 4.14e-35]`, radius **incomplete** | — |
| `‖f_200‖²` | `phase4/certs/q_f200.json` | `norm_sq` | 1100 | `[2.275244443650965360970550 +/- 8.18e-26]` | `phase4/verify_t400.py` |
| `c = −g'(0+) = f(L)²` | same | `c`, `f_at_L` | 1100 | `c = [4.538954705978261920675496e-15 +/- 3.26e-40]` | `phase4/verify_t400.py` |
| `Q(phi_0)` vs the repository closed form | `phase5/certs/conventions.json` | `phi_0.Q` | 900 | `[0.0802740549592696112136916593813 +/- 1.70e-32]` | `phase5/src/conventions.py` |
| `Q(phi_2)`, `Q(phi_0 − phi_2)` and direct pole moments | same | `phi_2.*`, `phi_0 - phi_2.*` | 900 | see file | `phase5/src/conventions.py` |

## 3. Theorem (iii) — the parity separation

| claim | file | field | prec | enclosure | verifier |
|---|---|---|---|---|---|
| `mu_odd = 2.3e-17` certifies | `phase5/certs/parity.json` | `odd.mu_odd`, `odd.ok` | 2000 | margin `[1.21161369506e-18 +/- 1.94e-30]` | `note/consistency_check.py` P1–P2, P7 |
| separation `mu_odd − U` | derived in the checker from the certified endpoint | — | 600 | `6.268744659293181e-18` — on the 160-mode witness | `note/consistency_check.py` P3a–P4 |
| the witness is even | `phase5/certs/upper_highN.json` | `N` = 160 even cosine modes | exact | — | `note/consistency_check.py` P6 |
| even sector binds | `phase5/certs/parity.json` | `whole_window_lower` | — | `"1.41e-17 (even sector binds)"` | `note/consistency_check.py` T8, P5 |

## 4. The gate, and monotonicity in `T`

| claim | file | field | prec | enclosure | verifier |
|---|---|---|---|---|---|
| `A_L(4/5)` | `phase4/certs/cutoff.json` | `A_L` | 1200 | `[2.941973525223620455503909 +/- 3.91e-25]` | `phase4/src/cutoff.py` |
| `beta_200` | same | `beta_200` | 1200 | `[0.5184657332438467199873536 +/- 2.62e-26]` | `phase2/verify_phase2.py` |
| `beta_400` | same | `beta_400` | 1200 | `[1.211613695058064549471502 +/- 1.89e-26]` | `phase4/verify_t400.py` |
| `beta_800` | `phase5/certs/pilot800.json` | `beta_800` | 1600 | `[1.904761070930776888621206 +/- 1.68e-25]` | `note/consistency_check.py` G1–G3 |
| increment `beta_400 − beta_200` | `phase4/certs/cutoff.json` | `beta_increase` | 1200 | `[0.6931479618142178294841484 +/- 7.24e-27]` | `note/consistency_check.py` G3 |
| beta-monotonicity identity residual | `phase3/certs/beta_monotonicity.json` | `identity_max_residual` | 700 | `[5.2292989671653341567e-143 +/- 4.54e-163]` | `phase3/verify_audit.py` |
| `K_T` diagonal range | same | `K_diag_min`, `K_diag_max` | 700 | `[0.94208886376952044561 …]`, `[0.99801326000383564212 …]` | `phase3/verify_audit.py` |

## 5. The scalar table and the discards

| claim | file | field | prec | enclosure | verifier |
|---|---|---|---|---|---|
| `R_200(f_200)` | `phase4/certs/cross.json` | `R200_f200` | 1100 | `[1.031022030245989649334614e-17 +/- 4.88e-42]` | `phase4/verify_t400.py` |
| `R_400(f_200)` | same | `R400_f200` | 1100 | `[2.352053046731007948698226e-17 +/- 3.28e-42]` | `phase4/verify_t400.py` |
| `R_800(f_200)` | `phase5/certs/pilot800.json` | `R_800(f_200)` | 1600 | `[3.01603e-17 +/- 5.30e-23]` | `note/consistency_check.py` C1 |
| `R_200(f_400)` | same | `R_200(f_400)` | 1600 | `[1.24275215e-17 +/- 8.10e-26]` | `note/consistency_check.py` C1 |
| `R_400(f_400)` | `phase4/certs/cross.json` | `R400_f400` | 1100 | `[1.433072404735585029967874e-17 +/- 4.94e-42]` | `phase4/verify_t400.py` |
| `R_800(f_400)` | `phase5/certs/pilot800.json` | `R_800(f_400)` | 1600 | `[1.6681e-17 +/- 4.60e-22]` | `note/consistency_check.py` C1 |
| quadrature remainders for the scalar legs | same | `*_quad_err` | 1600 | see file | — |
| `e_400(f_400)` — **MEASUREMENT** (depends on the withdrawn `q`) | `phase4/certs/cross.json` | `e400_f400` | 1100 | `[4.6831050640987937e-18 +/- 6.00e-35]` | `note/consistency_check.py` S3–S4 |
| discard ratios — **MEASUREMENTS** | `phase5/certs/pilot800.json` | `ratio_f_400`, `ratio_f_200` | 1600 | `[0.498 +/- 2.36e-4]`, `[0.48788 +/- 9.07e-6]` | `note/consistency_check.py` D1–D3 |
| ratios at the `200 -> 400` step | recomputed | — | 600 | `[0.7110365 +/- 1.61e-8]` (`f_400`), `[0.49532156 +/- 4.69e-9]` (`f_200`) | `note/consistency_check.py` D3 |

## 6. Ceilings, costs, and the superseded baseline

| claim | file | field | prec | value | verifier |
|---|---|---|---|---|---|
| ceiling `≤ 1.0310221e-17`, headroom `3.1023%` at `T = 200` | outward endpoint of `R200_f200` | — | 600 | endpoint `1.0310220302e-17` | `note/consistency_check.py` H1 |
| ceiling `≤ 1.4330725e-17`, headroom `1.6364%` at `T = 400` | outward endpoint of `R400_f400` | — | 600 | endpoint `1.4330724047e-17` | `note/consistency_check.py` H2 |
| ceiling `≤ 1.6681461e-17`, headroom `18.309%` at `T = 800` | outward endpoint of `R_800(f_400)` | — | 600 | endpoint `1.6681460e-17`, factor `1.1831` | `note/consistency_check.py` H3–H4 |
| assembly times and peak RSS | `phase2/certs/finite_block_certificate.json`, `phase4/certs/finite_block_certificate_T400.json` | `elapsed_s`, `peak_rss_mib` | — | `1343.1 s / 205.5 MiB`, `3580.8 s / 393.4 MiB` | measured, not certified |
| scalar-leg times | `phase5/logs/pilot800.log` | — | — | `933 s`, `8049 s`, `3651 s` (total `12634 s`) | measured, not certified |
| superseded `T = 200` certificate | `phase2/certs/whole_window_certificate.json` | `mu`, `status` | 1400 | `1e-17`, `CERTIFIED` | `phase2/verify_phase3.py`; `note/consistency_check.py` B1–B4 |

## 7. Claims that are measurements, not certificates

These appear in the note and are deliberately **not** certified claims. They are
reported as measured facts about this machine and this implementation, and no
certificate depends on them:

* every runtime and peak-memory figure;
* the extrapolated "order of ten hours" for a full `T = 800` assembly;
* the `2.64x` assembly growth and the per-Bessel cost observation;
* the discard ratios (see Section 5) — **measurements, not enclosures.** They
  are quotients of `e_T = q - r_T`, and `q` comes from the time-domain
  evaluator whose outer quadrature remainder is unbounded (`RETRACTIONS.md`
  R14); an arithmetic-only radius on such a quantity is not an enclosure. The
  checker's `D1`-`D3` verify the *arithmetic of the division*, which is correct;
  that is not a certification of the inputs. No decay law is inferred from them
  and they enter no bound.
