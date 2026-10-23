#!/usr/bin/env python3
r"""Phase 4d: the crossover criterion, assembled from the committed CSVs only.

Phase 4c measures |F'(gamma_1)|^2/K at the PERTURBED minimiser and at the UNPERTURBED
ground state.  Where the two agree the minimiser did not move and C is a property of the
ground state; where they differ it re-optimised.  The proposed criterion is that the
minimiser re-optimises exactly when the available gain exceeds the cost of rotating out of
the ground state:

    R(L, delta) = 4 delta^2 K(L, gamma_1)  /  lambda_2(Q)        re-optimises iff R > 1,

the numerator being the most the delta^2 term can pay (Theorem 8's coefficient at the
extremal) and the denominator the price of the first available direction.

No new computation: lambda_2 comes from data/phase3_gateB.csv, K and the two saturations
from data/phase4c_decomposition.csv.  Pure arithmetic, so it is cheap and auditable.
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp                                                      # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    mp.mp.dps = 30
    gb = list(csv.DictReader(open(os.path.join(HERE, 'data', 'phase3_gateB.csv'))))
    dc = list(csv.DictReader(open(os.path.join(HERE, 'data',
                                               'phase4c_decomposition.csv'))))
    rows, bad = [], 0
    print(f"  {'L':>5s} {'function':9s} {'K':>10s} {'4d^2 K':>11s} {'lambda_2':>13s} "
          f"{'R':>12s} {'|Fp|^2/K':>11s} {'gs value':>12s} {'ratio':>10s} "
          f"{'moved?':>7s} {'R>1?':>6s} {'agree':>6s}")
    for r in dc:
        L, lab, de = r['L'], r['label'], mp.mpf(r['delta'])
        m = [g for g in gb if g['label'] == lab and g['L'] == L
             and g['sector'] == 'even_d' and g['N'] == '16']
        if not m:
            continue
        lam2 = mp.mpf(m[0]['lam_2'])
        K = mp.mpf(r['K'])
        gain = 4 * de ** 2 * K
        R = gain / lam2
        a, b = mp.mpf(r['Fp2_over_K']), mp.mpf(r['Fp_gs2_over_K'])
        ratio = a / b
        moved = r['minimiser_moved'] == 'yes'
        pred = bool(R > 1)
        ok = (moved == pred)
        if not ok:
            bad += 1
        print(f"  {L:>5s} {lab:9s} {mp.nstr(K, 6):>10s} {mp.nstr(gain, 5):>11s} "
              f"{mp.nstr(lam2, 6):>13s} {mp.nstr(R, 5):>12s} {mp.nstr(a, 6):>11s} "
              f"{mp.nstr(b, 6):>12s} {mp.nstr(ratio, 5):>10s} "
              f"{('YES' if moved else 'no'):>7s} {str(pred):>6s} "
              f"{('ok' if ok else 'MISS'):>6s}")
        rows.append(dict(L=L, label=lab, delta=r['delta'], K=mp.nstr(K, 10),
                         gain_4d2K=mp.nstr(gain, 10), lam_2=mp.nstr(lam2, 10),
                         R=mp.nstr(R, 10), Fp2_over_K=mp.nstr(a, 10),
                         Fp_gs2_over_K=mp.nstr(b, 10), ratio=mp.nstr(ratio, 8),
                         minimiser_moved=r['minimiser_moved'],
                         R_gt_1=str(pred), agree=('ok' if ok else 'MISS')))
    mv = [mp.mpf(x['R']) for x in rows if x['minimiser_moved'] == 'yes']
    fz = [mp.mpf(x['R']) for x in rows if x['minimiser_moved'] != 'yes']
    print(f"\n  cells: {len(rows)}   disagreements with R > 1: {bad}")
    if mv and fz:
        print(f"  smallest R among MOVED  : {mp.nstr(min(mv), 6)}")
        print(f"  largest  R among frozen : {mp.nstr(max(fz), 6)}")
        print(f"  the two classes are separated by a factor "
              f"{mp.nstr(min(mv) / max(fz), 6)} with the threshold R = 1 inside the gap: "
              f"{bool(max(fz) < 1 < min(mv))}")
    with open(os.path.join(HERE, 'data', 'phase4d_crossover.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print('  wrote data/phase4d_crossover.csv')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
