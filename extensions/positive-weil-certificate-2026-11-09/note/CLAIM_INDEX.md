# CLAIM_INDEX.md — every claim in the note, mapped to its artifact

Each row gives: the claim as it appears in
`weil_L45_certified_positivity.tex`, the artifact file that holds it, the field
inside that file, the working precision at which it was produced, the published
enclosure, and the verifier that re-derives it.

All paths are relative to `extensions/positive-weil-certificate-2026-11-09/`.
Printed decimals in the note are **summaries**; the enclosures below are
authoritative. A serialized enclosure is a ball that *contains* the quantity, so
it must be consumed through the correct endpoint: lower endpoint for a lower
bound, upper endpoint for an upper bound. (`note/consistency_check.py` check
`U1` enforces exactly this for the one place it matters.)

---

## 1. Theorem (i) — the whole-window lower bound

| claim | file | field | prec | enclosure | verifier |
|---|---|---|---|---|---|
| `Q(f) >= 1.41e-17 ‖f‖²`, both sectors | `phase4/certs/certificate_T400.json` | `whole_window_mu` | 2000 | `[1.410000000000000000000000e-17 +/- 1e-46]` | `phase4/verify_t400.py`, `note/consistency_check.py` T1–T8 |
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

## 2. Theorem (ii) — the original-form upper bound

| claim | file | field | prec | enclosure | verifier |
|---|---|---|---|---|---|
| `q(f_400) = Q(f_400)/‖f_400‖²` | `phase4/certs/cross.json` | `Q_f400` | 1100 | `[1.90138291114546440e-17 +/- 6.03e-35]` | `phase4/verify_t400.py`; `note/consistency_check.py` U1–U6 |
| strict upper endpoint `U` | `phase5/certs/parity.json` | `U_upper_endpoint` | 2000 | `[1.901382911145464406030000e-17 +/- 2.04e-43]` | `note/consistency_check.py` U1 |
| safe quoted decimal | same | `U_safe_decimal` | — | `1.9014e-17` | `note/consistency_check.py` U3–U5 |
| frozen vector `f_400` (220 even modes, degrees 0..438) | `phase4/certs/cross.json` | `f400_y` | exact | rationals, denominator `10^40` | `phase4/verify_t400.py` |
| frozen vector `f_200` (110 even modes, degrees 0..218) | `phase3/certs/frozen_minimizer.json` | `y` | exact | rationals, denominator `10^40` | `phase3/verify_audit.py` |
| `‖f_200‖²` | `phase4/certs/q_f200.json` | `norm_sq` | 1100 | `[2.275244443650965360970550 +/- 8.18e-26]` | `phase4/verify_t400.py` |
| `q(f_200)` | `phase4/certs/q_f200.json` | `q` | 1100 | `[3.64859175290407702e-17 +/- 4.14e-35]` | `phase4/verify_t400.py` |
| pole / Arch / comb split for `f_200` | same | `pole`, `arch`, `prime` | 1100 | see file | `phase4/verify_t400.py` |
| `c = −g'(0+) = f(L)²` | same | `c`, `f_at_L` | 1100 | `c = [4.538954705978261920675496e-15 +/- 3.26e-40]` | `phase4/verify_t400.py` |
| `Q(phi_0)` vs the repository closed form | `phase5/certs/conventions.json` | `phi_0.Q` | 900 | `[0.0802740549592696112136916593813 +/- 1.70e-32]` | `phase5/src/conventions.py` |
| `Q(phi_2)`, `Q(phi_0 − phi_2)` and direct pole moments | same | `phi_2.*`, `phi_0 - phi_2.*` | 900 | see file | `phase5/src/conventions.py` |

## 3. Theorem (iii) — the parity separation

| claim | file | field | prec | enclosure | verifier |
|---|---|---|---|---|---|
| `mu_odd = 2.3e-17` certifies | `phase5/certs/parity.json` | `odd.mu_odd`, `odd.ok` | 2000 | margin `[1.21161369506e-18 +/- 1.94e-30]` | `note/consistency_check.py` P1–P2, P7 |
| separation `mu_odd − U` | same | `separation` | 2000 | `[3.986170888545355939700000e-18 +/- 2.04e-43]` | `note/consistency_check.py` P3a–P4 |
| `f_400` is even | `phase4/certs/cross.json` | `f400_y` (220 even modes) | exact | — | `note/consistency_check.py` P6 |
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
| `e_400(f_400)` | `phase4/certs/cross.json` | `e400_f400` | 1100 | `[4.6831050640987937e-18 +/- 6.00e-35]` | `note/consistency_check.py` S3–S4 |
| discard ratios | `phase5/certs/pilot800.json` | `ratio_f_400`, `ratio_f_200` | 1600 | `[0.498 +/- 2.36e-4]`, `[0.48788 +/- 9.07e-6]` | `note/consistency_check.py` D1–D3 |
| ratios at the `200 -> 400` step | recomputed | — | 600 | `[0.7110365 +/- 1.61e-8]` (`f_400`), `[0.49532156 +/- 4.69e-9]` (`f_200`) | `note/consistency_check.py` D3 |

## 6. Ceilings, costs, and the superseded baseline

| claim | file | field | prec | value | verifier |
|---|---|---|---|---|---|
| `3.102%` headroom at `T = 200` | derived from `R200_f200` | — | 600 | `[3.1022030 +/- 2.46e-8]%` | `note/consistency_check.py` H1 |
| `1.636%` headroom at `T = 400` | derived from `R400_f400` | — | 600 | `[1.6363408 +/- 3.87e-8]%` | `note/consistency_check.py` H2 |
| `18.31%` ceiling at `T = 800` | derived from `R_800(f_400)` | — | 600 | `[18.30 +/- 8.23e-3]%`, factor `1.1830` | `note/consistency_check.py` H3–H4 |
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
* the discard ratios (see Section 5) — they are certified *enclosures*, but no
  decay law is inferred from them and they enter no bound.
