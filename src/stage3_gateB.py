#!/usr/bin/env python3
"""Gate B: lambda_min of the L=0.8 window form vs Zhu's certified [8.9e-18, 2.27e-17]."""
import sys, os, csv, argparse, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from zeta_window import build_matrix
from weilform import eigsym, parity_of

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ZHU_LO, ZHU_HI = mp.mpf('8.9e-18'), mp.mpf('2.27e-17')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', default='0.8')
    ap.add_argument('--N', default='2,4,6,8,10,12,14,16,18,20')
    ap.add_argument('--dps', default='40,70')
    ap.add_argument('--out', default='stage3_gateB_lambda_min.csv')
    a = ap.parse_args()
    Ns = [int(x) for x in a.N.split(',')]
    rows = []
    print(f"Gate B: L={a.L}; Zhu certified lambda* in [{ZHU_LO}, {ZHU_HI}]")
    for dps in [int(x) for x in a.dps.split(',')]:
        print(f"\n-- working precision {dps} decimal digits --")
        prev = None
        for N in Ns:
            mp.mp.dps = dps
            t0 = time.time()
            M = build_matrix(mp.mpf(a.L), N, dps)
            vals, vecs = eigsym(M)
            lam, lam2 = vals[0], (vals[1] if N > 1 else mp.inf)
            par, pres = parity_of(vecs[0])
            # the basis is all-even by construction; 'parity' here is the COEFFICIENT
            # flip symmetry, which is not a symmetry of this problem -- recorded only
            # to show the ground state is not a basis artefact.
            mono = '-' if prev is None else ('OK' if lam <= prev * (1 + mp.mpf(10) ** -20)
                                             + mp.mpf(10) ** (-(dps - 8)) else 'VIOLATED')
            inz = bool(ZHU_LO <= lam <= ZHU_HI)
            print(f"   N={N:3d} lam_min={mp.nstr(lam,10):>16s} sign={'+' if lam>0 else '-'}"
                  f" lam_2={mp.nstr(lam2,8):>12s} gap/lam_min={mp.nstr(lam2/lam,6):>10s}"
                  f" inZhu={str(inz):5s} monotone={mono:8s} ({time.time()-t0:.0f}s)")
            rows.append(dict(L=a.L, N=N, dps=dps, lam_min=mp.nstr(lam, 12),
                             lam_min_sign=('pos' if lam > 0 else 'neg'),
                             lam_2=mp.nstr(lam2, 10),
                             ratio_lam2_over_lam1=mp.nstr(lam2 / lam, 8),
                             in_zhu_interval=inz, monotone=mono,
                             coeff_flip_parity=par))
            prev = lam
    with open(os.path.join(ROOT, 'results', a.out), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote results/{a.out}")


if __name__ == '__main__':
    main()
