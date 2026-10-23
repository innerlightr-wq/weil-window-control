#!/usr/bin/env python3
r"""Phase 3 Gate B: lambda_min(N) for L(s,chi) in two independent complete even bases.

  EVEN_D : cos((2k+1) pi x / 2L)     vanishes at +-L
  EVEN   : cos(k pi x / L)           contains the constant, jumps at +-L
Both are complete in L^2_even[-L, L], so lambda_min(N) is non-increasing in N and the two
must converge to the same lambda*_chi(L).  Agreement between them is the cross-validation.

A zeta control is recomputed at every L in the same run.  At L = 0.8 it must land in Zhu's
certified interval [8.9e-18, 2.27e-17]; that is the external anchor, and it is the test
that a sector or sign bug in the chi code has not also broken the shared code path.

EVERY chi-side lambda_min is a VARIATIONAL UPPER BOUND for lambda*_chi(L) and is never used
as a lower bound (REPORT.md non-claim 8).  No certified positivity input exists on the chi
side.
"""
import argparse
import csv
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp                                                      # noqa: E402
import dirichlet_window as dw                                            # noqa: E402
from zeta_window import build_matrix, EVEN, EVEN_D                       # noqa: E402
from weilform import eigsym                                              # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ZHU_LO, ZHU_HI = mp.mpf('8.9e-18'), mp.mpf('2.27e-17')


def sweep(build, Ns, dps):
    prev, out = None, []
    for N in Ns:
        t0 = time.time()
        M = build(N)
        vals, _ = eigsym(M)
        lam, lam2 = vals[0], vals[1]
        mono = '-' if prev is None else ('OK' if lam <= prev + mp.mpf(10) ** (-(dps - 8))
                                         else 'VIOLATED')
        out.append((N, lam, lam2, mono, time.time() - t0))
        prev = lam
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--Ls', default='0.5,0.8,1.0,1.3')
    ap.add_argument('--Ns', default='4,6,8,10,12,14,16')
    ap.add_argument('--dps', type=int, default=50)
    ap.add_argument('--out', default='phase3_gateB')
    a = ap.parse_args()
    mp.mp.dps = a.dps
    Ns = [int(x) for x in a.Ns.split(',')]
    rows, summary = [], {}

    for Lx in a.Ls.split(','):
        L = mp.mpf(Lx)
        print("=" * 100)
        print(f"L = {Lx}   (dps = {a.dps})")
        print("=" * 100)
        for sector in (EVEN_D, EVEN):
            # ---- zeta control, same run, same code path
            res = sweep(lambda N: build_matrix(L, N, a.dps, sector), Ns, a.dps)
            N, lam, lam2, mono, dt = res[-1]
            inz = bool(ZHU_LO <= lam <= ZHU_HI)
            print(f"  zeta    {sector:7s} N={N:3d} lam_min={mp.nstr(lam, 10):>16s} "
                  f"lam_2={mp.nstr(lam2, 8):>12s} mono={mono:9s} "
                  f"inZhu={str(inz) if abs(float(L) - 0.8) < 1e-12 else '-':5s}")
            for N_, lam_, lam2_, mono_, dt_ in res:
                rows.append(dict(label='zeta', d='', q=1, a_par=0, L=Lx, sector=sector,
                                 N=N_, dps=a.dps, lam_min=mp.nstr(lam_, 14),
                                 sign=('pos' if lam_ > 0 else 'neg'),
                                 lam_2=mp.nstr(lam2_, 10), monotone=mono_,
                                 seconds=round(dt_, 1)))
            summary.setdefault('zeta', {})[f'{Lx}|{sector}'] = dict(
                lam_min=mp.nstr(lam, 14), N=N, monotone_all=all(
                    r[3] in ('OK', '-') for r in res),
                in_zhu=(inz if abs(float(L) - 0.8) < 1e-12 else None))
            lam_zeta = lam
            # ---- the five characters
            for d in dw.DISCRIMINANTS:
                lab = dw.CHARS[d]['label']
                res = sweep(lambda N, d=d: dw.build_matrix_chi(L, range(N), a.dps, d,
                                                               sector), Ns, a.dps)
                N, lam, lam2, mono, dt = res[-1]
                ratio = lam / lam_zeta if lam_zeta != 0 else mp.nan
                print(f"  {lab:7s} {sector:7s} N={N:3d} lam_min={mp.nstr(lam, 10):>16s} "
                      f"lam_2={mp.nstr(lam2, 8):>12s} mono={mono:9s} "
                      f"lam/lam_zeta={mp.nstr(ratio, 6)}")
                for N_, lam_, lam2_, mono_, dt_ in res:
                    rows.append(dict(label=lab, d=d, q=dw.CHARS[d]['q'],
                                     a_par=dw.CHARS[d]['a'], L=Lx, sector=sector, N=N_,
                                     dps=a.dps, lam_min=mp.nstr(lam_, 14),
                                     sign=('pos' if lam_ > 0 else 'neg'),
                                     lam_2=mp.nstr(lam2_, 10), monotone=mono_,
                                     seconds=round(dt_, 1)))
                summary.setdefault(lab, {})[f'{Lx}|{sector}'] = dict(
                    lam_min=mp.nstr(lam, 14), N=N,
                    monotone_all=all(r[3] in ('OK', '-') for r in res),
                    ratio_to_zeta=mp.nstr(ratio, 8))
        print()
    with open(os.path.join(HERE, 'data', f'{a.out}.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(HERE, 'results', f'{a.out}.json'), 'w') as fh:
        json.dump(summary, fh, indent=1)
    print(f'wrote data/{a.out}.csv and results/{a.out}.json')


if __name__ == '__main__':
    main()
