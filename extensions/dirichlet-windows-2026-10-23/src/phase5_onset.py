#!/usr/bin/env python3
r"""Phase 5b: at what L does a planted REAL zero at 1/2 + delta first produce a negative
Rayleigh quotient?

Phase 5a establishes that the response is rank one, PSD on the even sector and NSD on the
odd sector, so the question only has content in the ODD sector, where

    Q^odd - 2 a a^T  is indefinite   <=>   Theta(L, delta) := 2 a^T (Q^odd)^{-1} a > 1

(Schur complement, where Q^odd is positive definite on the subspace).  This script scans L
at fixed basis dimension and records Theta and lambda_min, with an explicit LDL witness
wherever lambda_min < 0, so each detection is a one-sided certificate: a single trial
function with Q < 0.

Theta scales like delta^2 at fixed L (checked), so the table is reported for each delta
separately rather than rescaled.

This is a statement about the DETECTOR (REPORT.md non-claim 2).
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
import planting as pl                                                    # noqa: E402
from zeta_window import build_matrix, ODD, EVEN_D                        # noqa: E402
from weilform import eigsym                                             # noqa: E402
from ldl import ldl, first_negative_pivot, witness, rayleigh             # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--Ls', default='0.5,0.8,1.1,1.4,1.7,2.0,2.3,2.6,3.0')
    ap.add_argument('--deltas', default='0.1,0.25,0.375,0.45,0.49')
    ap.add_argument('--N', type=int, default=12)
    ap.add_argument('--sector', default='odd')
    ap.add_argument('--out', default='phase5_onset')
    a = ap.parse_args()
    deltas = [mp.mpf(x) for x in a.deltas.split(',')]
    sector = a.sector
    rows = []
    labels = ['zeta'] + [dw.CHARS[d]['label'] for d in dw.DISCRIMINANTS]
    print(f"sector = {sector}, N = {a.N}; Theta > 1 <=> indefinite\n")
    hdr = f"  {'L':>5s} {'dps':>4s} {'function':9s} {'lam_min(Q)':>15s} " + \
          "".join(f"{'Th(' + mp.nstr(d, 4) + ')':>13s}" for d in deltas) + \
          f" {'first det':>10s}"
    print(hdr)
    for Lx in a.Ls.split(','):
        L = mp.mpf(Lx)
        dps = max(60, int(30 + 40 * float(Lx)))
        mp.mp.dps = dps
        for lab in labels:
            t0 = time.time()
            if lab == 'zeta':
                M = build_matrix(L, a.N, dps, sector)
            else:
                d = [k for k in dw.DISCRIMINANTS if dw.CHARS[k]['label'] == lab][0]
                M = dw.build_matrix_chi(L, range(a.N), dps, d, sector)
            lam0 = eigsym(M)[0][0]
            ths, dets = [], []
            for de in deltas:
                B, Ap = pl.real_block(L, range(a.N), de, sector)
                try:
                    x = mp.lu_solve(M, mp.matrix(Ap))
                    th = 2 * sum(Ap[i] * x[i] for i in range(a.N))
                except Exception:                                       # noqa: BLE001
                    th = mp.nan
                Z = M + B
                lam = eigsym(Z)[0][0]
                wit, ray = '', ''
                try:
                    Lm, D = ldl(Z)
                    j = first_negative_pivot(D)
                    if j is not None:
                        v = witness(Lm, j)
                        r, _, _ = rayleigh(Z, v)
                        wit, ray = str(j + 1), mp.nstr(r, 10)
                except ZeroDivisionError:
                    pass
                ths.append(th)
                dets.append(lam < 0)
                rows.append(dict(label=lab, sector=sector, L=Lx, N=a.N, dps=dps,
                                 delta=mp.nstr(de, 8),
                                 lam_unperturbed=mp.nstr(lam0, 12),
                                 lam_planted=mp.nstr(lam, 12),
                                 theta=mp.nstr(th, 10),
                                 theta_gt_1=str(bool(th > 1)),
                                 detected=('yes' if lam < 0 else 'no'),
                                 witness_support=wit, min_rayleigh=ray,
                                 rank1_eig=mp.nstr(-2 * sum(y ** 2 for y in Ap), 10)
                                 if sector == ODD else
                                 mp.nstr(2 * sum(y ** 2 for y in Ap), 10)))
            first = next((mp.nstr(d, 4) for d, ok in zip(deltas, dets) if ok), '-')
            print(f"  {Lx:>5s} {dps:>4d} {lab:9s} {mp.nstr(lam0, 12):>15s} "
                  + "".join(f"{mp.nstr(t, 7):>13s}" for t in ths)
                  + f" {first:>10s}  ({time.time() - t0:.0f}s)")
        print()
    with open(os.path.join(HERE, 'data', f'{a.out}_{sector}.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(HERE, 'results', f'{a.out}_{sector}.json'), 'w') as fh:
        json.dump(rows, fh, indent=1)
    # onset summary: smallest scanned L with a detection, per (function, delta)
    print("ONSET SUMMARY: smallest scanned L with lambda_min < 0")
    print(f"  {'function':9s} " + "".join(f"{'d=' + mp.nstr(d, 4):>12s}"
                                          for d in deltas))
    summ = {}
    for lab in labels:
        cells = []
        for de in deltas:
            ok = [r for r in rows if r['label'] == lab
                  and r['delta'] == mp.nstr(de, 8) and r['detected'] == 'yes']
            cells.append(min((float(r['L']) for r in ok), default=None))
        summ[lab] = {mp.nstr(d, 4): c for d, c in zip(deltas, cells)}
        print(f"  {lab:9s} " + "".join(f"{(str(c) if c is not None else 'none'):>12s}"
                                       for c in cells))
    with open(os.path.join(HERE, 'results', f'{a.out}_{sector}_summary.json'), 'w') as fh:
        json.dump(summ, fh, indent=1)
    print(f'wrote data/{a.out}_{sector}.csv')


if __name__ == '__main__':
    main()
