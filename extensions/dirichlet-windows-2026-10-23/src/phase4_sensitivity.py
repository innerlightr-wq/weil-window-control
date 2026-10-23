#!/usr/bin/env python3
r"""Phase 4: does the sensitivity law transfer to L(s, chi) with only the baseline changed?

Theorem 8's delta^2 coefficient is 4 K(L, gamma), K(L,gamma) = int u^2 sin^2(gamma u) du.
K is the squared norm of the derivative-evaluation functional f -> F'(gamma) on the
normalised admissible sector, and that sector does not know which L-function it is being
tested against -- so 4K should transfer UNCHANGED, and only lambda*_chi(L) should move.

This script measures, at the first ordinate gamma_1(chi) of each character:

   C(L) = [ lambda_min(Q_0) - lambda_min(Q_delta) ] / delta^2

for BOTH counterfactuals of notes/proofs.md B.1, kept distinct (paper S10):

   added quartet : Q_delta = Q + quartet(gamma_1, delta),   Q_0 = Q + 2 pair(gamma_1)
   moved pair    : Q_delta = Q - pair + quartet,            Q_0 = Q + pair(gamma_1)

and reports C/(4K(L, gamma_1)) next to the zeta ratio computed in the same run.

The parent repository's C_measured = -lambda_min(Q_delta)/delta^2 omits the baseline
because for zeta lambda_min(Q_0) ~ 1e-17.  That is NOT true for chi (there is no pole term
to cancel against), so the baseline-subtracted C is the quantity compared, and the raw
parent-convention value is reported alongside so the two can be told apart.
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
from zeta_window import build_matrix, EVEN_D                             # noqa: E402
from weilform import eigsym                                              # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NTAB = {'0.8': 16, '1.0': 16, '1.3': 20, '1.6': 22}


def lmin(M):
    return eigsym(M)[0][0]


def run_one(label, M, L, g1, deltas, N, sector=EVEN_D):
    idx = range(N)
    P = pl.pair_block(L, idx, g1, sector)
    base = {'added': M + 2 * P, 'moved': M + P}
    lam0 = {k: lmin(v) for k, v in base.items()}
    K = pl.K_closed(L, g1)
    out = []
    for de in deltas:
        Qd = pl.quartet_block(L, idx, g1, de, sector)
        for conv, Z in (('added', M + Qd), ('moved', M - P + Qd)):
            lam = lmin(Z)
            C_sub = (lam0[conv] - lam) / de ** 2
            C_raw = -lam / de ** 2
            out.append(dict(label=label, L=mp.nstr(L, 6), N=N, gamma_1=mp.nstr(g1, 14),
                            delta=mp.nstr(de, 8), convention=conv,
                            lam_unperturbed=mp.nstr(lmin(M), 12),
                            lam_base0=mp.nstr(lam0[conv], 12),
                            lam_delta=mp.nstr(lam, 12),
                            C_sub=mp.nstr(C_sub, 10), C_raw=mp.nstr(C_raw, 10),
                            K=mp.nstr(K, 12), fourK=mp.nstr(4 * K, 12),
                            ratio_C_over_4K=mp.nstr(C_sub / (4 * K), 8),
                            detected=('yes' if lam < 0 else 'no')))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--Ls', default='0.8,1.0,1.3')
    ap.add_argument('--deltas', default='0.25,0.375,0.02')
    ap.add_argument('--out', default='phase4_sensitivity')
    a = ap.parse_args()
    deltas = [mp.mpf(x) for x in a.deltas.split(',')]
    rows = []
    for Lx in a.Ls.split(','):
        L = mp.mpf(Lx)
        N = NTAB[Lx]
        dps = max(60, int(40 + 50 * float(Lx)))
        mp.mp.dps = dps
        print("=" * 110)
        print(f"L = {Lx}, N = {N}, dps = {dps}, sector = even_d")
        print("=" * 110)
        print(f"  {'function':8s} {'gamma_1':>12s} {'4K(L,g1)':>12s} {'delta':>8s} "
              f"{'conv':>6s} {'lam_0':>14s} {'lam_delta':>14s} {'C (sub)':>13s} "
              f"{'C/4K':>9s} {'det':>4s}")
        t0 = time.time()
        g1z = mp.im(mp.zetazero(1))
        Mz = build_matrix(L, N, dps, EVEN_D)
        rz = run_one('zeta', Mz, L, g1z, deltas, N)
        rows += rz
        for r in rz:
            print(f"  {'zeta':8s} {r['gamma_1']:>12s} {r['fourK']:>12s} {r['delta']:>8s} "
                  f"{r['convention']:>6s} {r['lam_base0']:>14s} {r['lam_delta']:>14s} "
                  f"{r['C_sub']:>13s} {r['ratio_C_over_4K']:>9s} {r['detected']:>4s}")
        for d in dw.DISCRIMINANTS:
            lab = dw.CHARS[d]['label']
            gs = [mp.mpf(x) for x in
                  open(os.path.join(HERE, 'data', f'zeros_{lab}.txt'))
                  if not x.startswith('#')]
            M = dw.build_matrix_chi(L, range(N), dps, d, EVEN_D)
            rr = run_one(lab, M, L, gs[0], deltas, N)
            rows += rr
            for r in rr:
                print(f"  {lab:8s} {r['gamma_1']:>12s} {r['fourK']:>12s} "
                      f"{r['delta']:>8s} {r['convention']:>6s} {r['lam_base0']:>14s} "
                      f"{r['lam_delta']:>14s} {r['C_sub']:>13s} "
                      f"{r['ratio_C_over_4K']:>9s} {r['detected']:>4s}")
        print(f"  ({time.time() - t0:.0f}s)\n")
    with open(os.path.join(HERE, 'data', f'{a.out}.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(HERE, 'results', f'{a.out}.json'), 'w') as fh:
        json.dump(rows, fh, indent=1)
    print(f'wrote data/{a.out}.csv')


if __name__ == '__main__':
    main()
