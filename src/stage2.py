#!/usr/bin/env python3
"""Stage 2: teeth. Planted off-circle Frobenius configurations.

Deliverables per config:
  - first R where lam_min(T_R) < 0 (Weil positivity breaks),
  - whether the ground-state zeros still lie on |z|=1 (H1),
  - what those zeros converge to as R grows,
  - whether simplicity or evenness of the ground state fails anywhere.
"""
import sys, os, csv, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from planted import t_from_blocks, betas, predicted_signature
from exact_inertia import rational_t, exact_inertia
from analysis import analyse_window
from weilform import toeplitz, inertia

RHO = '1.3'
# rational surrogates: cos(0.7)=0.764842 ~ 19/25 (phi=0.701369);
#                      cos(2.0)=-0.416147 ~ -21/50 (phi=2.004589).
# These make every t(n) rational, so the inertia is certified exactly over Q.
COS07, COS20 = '19/25', '-21/50'

CONFIGS = [
    # (tag, description, blocks)
    ('C0-online', 'CONTROL: on-line, angles +-0.7, +-2.0 (g_eff=2, "RH true")',
     [('online', 0.7), ('online', 2.0)]),
    ('C1-g1real', 'g=1 tooth: single real FE pair rho=1.3  (|a| = 2.069 sqrt q > 2 sqrt q)',
     [('realpair', RHO)]),
    ('C2-quartet', 'g=2 tooth: one off-line quartet rho=1.3, phi=0.7',
     [('quartet', RHO, 0.7)]),
    ('C3-twoquartets', 'literal note config: quartets rho=1.3 at phi=0.7 and phi=2.0 (g_eff=4)',
     [('quartet', RHO, 0.7), ('quartet', RHO, 2.0)]),
    ('C4-mixed', 'mixed: off-line quartet (1.3, 0.7) + on-line pair at 2.0 (g_eff=3)',
     [('quartet', RHO, 0.7), ('online', 2.0)]),
]

# exactly-rational twins (cos phi rational) for T1 certification
RAT = {
    'C0-online':      [('online', COS07), ('online', COS20)],
    'C1-g1real':      [('realpair', '13/10')],
    'C2-quartet':     [('quartet', '13/10', COS07)],
    'C3-twoquartets': [('quartet', '13/10', COS07), ('quartet', '13/10', COS20)],
    'C4-mixed':       [('quartet', '13/10', COS07), ('online', COS20)],
}

# how far off the circle -> how soon positivity breaks
SWEEP = ['3/2', '13/10', '11/10', '21/20', '101/100', '1003/1000', '10001/10000']


def run_config(tag, desc, blocks, Rmax, dps, rows, zrows, rat_blocks=None):
    t, geff = t_from_blocks(blocks, Rmax)
    tq = None
    if rat_blocks is not None:
        tq, _ = rational_t(rat_blocks, Rmax)
    bs = betas(blocks)
    print(f"\n=== {tag} ===\n  {desc}")
    print(f"  g_eff={geff}   |beta_j| = {[mp.nstr(abs(b), 6) for b in bs]}")
    print(f"  t(0..{min(Rmax,8)}) = {[mp.nstr(x, 8) for x in t[:9]]}")
    first_neg = None
    prev_zang = None
    for R in range(0, Rmax + 1):
        row, v, zs = analyse_window(t, R, [mp.arg(b) for b in bs], dps)
        row.update(config=tag, g_eff=geff, rho=RHO)
        ps = predicted_signature(blocks, R)
        row['sig_pred'] = (f"({ps[0]},{ps[1]},{ps[2]})" if ps else '')
        row['sig_matches_pred'] = (row['inertia'] == row['sig_pred']) if ps else None
        if first_neg is None and row['lam_min_sign'] == 'neg' and row['prec_stable']:
            first_neg = R
        row['first_neg_so_far'] = first_neg
        if tq is not None:
            ex = exact_inertia(tq, R)
            row['exact_inertia_Q'] = (f"({ex[0]},{ex[1]})" if ex else 'singular')
            row['exact_lam_min_neg'] = ('' if ex is None else str(ex[0] >= 1))
            row['exact_agrees'] = ('' if ex is None else
                                   str((ex[0] >= 1) == (row['lam_min_sign'] == 'neg')))
        else:
            row['exact_inertia_Q'] = row['exact_lam_min_neg'] = row['exact_agrees'] = ''
        zang = sorted(float(mp.arg(z)) for z in zs)
        row['zero_angles'] = ';'.join('%.10f' % a for a in zang)
        row['zero_angle_drift'] = ''
        if prev_zang is not None and len(zang) == len(prev_zang) + 1:
            pass
        prev_zang = zang
        rows.append(row)
        for z in sorted(zs, key=lambda z: mp.arg(z)):
            zrows.append(dict(config=tag, R=R, zero_angle=mp.nstr(mp.arg(z), 14),
                              abs_zero=mp.nstr(abs(z), 14)))
        print(f"  R={R:2d} n={R+1:2d} lam_min={row['lam_min']:>16s}[{row['lam_min_sign']:>4s}]"
              f" err={row['lam_min_err_est']:>9s} certain={str(row['sign_certain']):5s}"
              f" f64={row['lam_min_float64']:>11s} simple={str(row['lam_min_simple']):5s}"
              f" par={row['parity']:6s} H1:max||z|-1|={row['H1_max_abs_zero_minus_1']:>10s}"
              f" inert={row['inertia']:>10s} pred={row['sig_pred']:>10s}"
              f" exactQ={row['exact_inertia_Q']:>10s}")
    print(f"  --> FIRST R with lam_min < 0 (precision-stable): {first_neg}")
    bad = [r['R'] for r in rows if r['config'] == tag and r['H1_applicable']
           and not r['H1_holds']]
    print(f"  --> H1 (simple lam_min => all zeros on |z|=1): "
          f"{'HOLDS at every applicable R' if not bad else 'FAILS at R=' + str(bad)}")
    nonsimple = [r['R'] for r in rows if r['config'] == tag and not r['lam_min_simple']]
    nonparity = [r['R'] for r in rows if r['config'] == tag and r['parity'] == 'mixed']
    print(f"  --> lam_min NOT simple at R = {nonsimple or 'nowhere'}")
    print(f"  --> ground state NOT of definite parity at R = {nonparity or 'nowhere'}")
    return first_neg


def run_sweep(Rmax, dps, srows):
    print("\n\n######## rho-sweep: how far off-circle vs. when positivity breaks ########")
    for rho in SWEEP:
        for kind, blocks in (('realpair', [('realpair', rho)]),
                             ('quartet0.7', [('quartet', rho, 0.7)])):
 
            t, geff = t_from_blocks(blocks, Rmax)
            rb = [('realpair', rho)] if kind == 'realpair' else [('quartet', rho, COS07)]
            tq, _ = rational_t(rb, Rmax)
            exact_first = None
            for R in range(0, Rmax + 1):
                ex = exact_inertia(tq, R)
                if ex is not None and ex[0] >= 1:
                    exact_first = R
                    break
            first_neg, lam_at_break = None, None
            lams = []
            for R in range(0, Rmax + 1):
                row, _, _ = analyse_window(t, R, [mp.mpf(0)], dps)
                lams.append(row['lam_min'])
                if first_neg is None and row['lam_min_sign'] == 'neg' and row['prec_stable']:
                    first_neg, lam_at_break = R, row['lam_min']
            srows.append(dict(kind=kind, rho=rho, first_neg_R=first_neg,
                              first_neg_R_EXACT_over_Q=exact_first,
                              lam_min_at_break=lam_at_break,
                              lam_min_at_Rmax=lams[-1], Rmax=Rmax, dps=dps))
            print(f"  {kind:12s} rho={rho:>7s}  first R with lam_min<0: "
                  f"mpmath={str(first_neg):>5s} exact/Q={str(exact_first):>5s}"
                  f"   lam_min there = {str(lam_at_break):>20s}"
                  f"   lam_min(R={Rmax}) = {lams[-1]}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dps', type=int, default=50)
    ap.add_argument('--Rmax', type=int, default=20)
    ap.add_argument('--out', default=os.path.join(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))), 'results'))
    a = ap.parse_args()
    mp.mp.dps = a.dps
    print(f"mpmath working precision: {mp.mp.dps} decimal digits; Rmax={a.Rmax}")
    print("NOTE: planted t(n) depends only on the normalised radii rho = |alpha|/sqrt(q),")
    print("      so q drops out of Stage 2 entirely (see src/planted.py).")
    rows, zrows, srows = [], [], []
    summary = {}
    for tag, desc, blocks in CONFIGS:
        summary[tag] = run_config(tag, desc, blocks, a.Rmax, a.dps, rows, zrows,
                                  RAT.get(tag))
    run_sweep(a.Rmax, a.dps, srows)
    os.makedirs(a.out, exist_ok=True)
    for name, data in (('stage2_spectra.csv', rows), ('stage2_zeros.csv', zrows),
                       ('stage2_rho_sweep.csv', srows)):
        with open(os.path.join(a.out, name), 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(data[0].keys()))
            w.writeheader(); w.writerows(data)
        print(f"wrote {os.path.join(a.out, name)}  ({len(data)} rows)")
    print("\n######## first-negative-R summary ########")
    for k, v in summary.items():
        print(f"  {k:16s} -> {v}")


if __name__ == '__main__':
    main()
