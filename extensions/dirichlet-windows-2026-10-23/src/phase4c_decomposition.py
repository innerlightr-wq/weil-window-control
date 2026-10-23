#!/usr/bin/env python3
r"""Phase 4c: WHY is C_chi/(4K) chi-dependent when 4K is not?

Theorem 8's coefficient is 4K(L,gamma) = 4 sup |F'(gamma)|^2, an extremal over a function
class that does not mention the L-function -- so it transfers unchanged.  Phase 4 finds that
the MEASURED coefficient C_chi = [lam_min(Q_0) - lam_min(Q_delta)]/delta^2 saturates it very
differently from character to character.  The parent repository's note (B.5) identifies the
mechanism for zeta: C is a property of the RE-OPTIMISED minimiser of Q_delta, and

    C = 4 ( F'(gamma)^2 + F(gamma) F''(gamma) )   at that minimiser,

with F(gamma) ~ 0 for zeta because the window minimiser nearly annihilates the leading
zeros.  This script measures the three pieces separately for each character:

    |F'|^2 / K            saturation of the DERIVATIVE extremal alone  (<= 1 by Cauchy-Schwarz)
    4 F F'' / C           how much of C is the cross term, i.e. the failure to annihilate
    F(gamma_1)^2          the annihilation residual itself

and, decisively, the SAME saturation at the UNPERTURBED ground state of Q_chi:

    |F'_gs|^2 / K         what the coefficient would be with no re-optimisation at all

If the two agree, the measured C is a property of the ground state and gamma_1(chi) and the
minimiser did not move; if they differ, the re-optimisation is doing the work.  That is the
one number that separates "the baseline changed" from "the detector behaves differently".
"""
import argparse
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp                                                      # noqa: E402
import dirichlet_window as dw                                            # noqa: E402
import planting as pl                                                    # noqa: E402
from zeta_window import build_matrix, F, norm2, EVEN_D                   # noqa: E402
from weilform import eigsym                                             # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NTAB = {'0.8': 16, '1.0': 16, '1.3': 20}


def Fderivs(v, L, g, N, sector=EVEN_D):
    L = mp.mpf(L)
    nrm = [mp.sqrt(norm2(k, L, sector)) for k in range(N)]

    def fun(t):
        return sum(v[k] * F(k, L, t, sector) / nrm[k] for k in range(N))
    return fun(g), mp.diff(fun, g), mp.diff(fun, g, 2)


def unit(v):
    n = mp.sqrt(sum(x ** 2 for x in v))
    return [x / n for x in v]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--Ls', default='0.8,1.0,1.3')
    ap.add_argument('--delta', default='0.02')
    a = ap.parse_args()
    de = mp.mpf(a.delta)
    rows = []
    for Lx in a.Ls.split(','):
        L = mp.mpf(Lx)
        N = NTAB[Lx]
        dps = max(60, int(40 + 50 * float(Lx)))
        mp.mp.dps = dps
        print("=" * 120)
        print(f"L = {Lx}, N = {N}, dps = {dps}, delta = {a.delta}, moved-pair at gamma_1")
        print("=" * 120)
        print(f"  {'function':9s} {'gamma_1':>9s} {'lam_0':>13s} {'C_sub':>12s} "
              f"{'4(Fp^2+F Fpp)':>14s} {'rel':>8s} {'K':>10s} {'|Fp|^2/K':>10s} "
              f"{'4F Fpp/C':>10s} {'F(g1)^2':>12s} {'C/4K':>9s} "
              f"{'|Fp_gs|^2/K':>12s} {'moved?':>8s}")
        jobs = [('zeta', mp.im(mp.zetazero(1)),
                 build_matrix(L, N, dps, EVEN_D))]
        for d in dw.DISCRIMINANTS:
            lab = dw.CHARS[d]['label']
            gs = [mp.mpf(x) for x in
                  open(os.path.join(HERE, 'data', f'zeros_{lab}.txt'))
                  if not x.startswith('#')]
            jobs.append((lab, gs[0],
                         dw.build_matrix_chi(L, range(N), dps, d, EVEN_D)))
        for lab, g1, M in jobs:
            P = pl.pair_block(L, range(N), g1, EVEN_D)
            Q0 = M + P                                  # moved-pair baseline (delta = 0)
            Qd = M - P + pl.quartet_block(L, range(N), g1, de, EVEN_D)
            v0 = eigsym(Q0)
            lam0 = v0[0][0]
            gs0 = unit(v0[1][0])
            G0, G1, G2 = Fderivs(gs0, L, g1, N)
            vals, vecs = eigsym(Qd)
            lam = vals[0]
            v = unit(vecs[0])
            C = (lam0 - lam) / de ** 2
            F0, F1, F2 = Fderivs(v, L, g1, N)
            pred = 4 * (F1 ** 2 + F0 * F2)
            K = pl.K_closed(L, g1)
            rel = abs(C - pred) / abs(C)
            print(f"  {lab:9s} {mp.nstr(g1, 7):>9s} {mp.nstr(lam0, 8):>13s} "
                  f"{mp.nstr(C, 8):>12s} {mp.nstr(pred, 8):>14s} {mp.nstr(rel, 3):>8s} "
                  f"{mp.nstr(K, 7):>10s} {mp.nstr(F1 ** 2 / K, 6):>10s} "
                  f"{mp.nstr(4 * F0 * F2 / C, 6):>10s} {mp.nstr(F0 ** 2, 6):>12s} "
                  f"{mp.nstr(C / (4 * K), 6):>9s} {mp.nstr(G1 ** 2 / K, 6):>12s} "
                  f"{('YES' if abs(G1 ** 2 - F1 ** 2) > abs(F1 ** 2) / 10 else 'no'):>8s}")
            rows.append(dict(label=lab, L=Lx, N=N, dps=dps, delta=a.delta,
                             gamma_1=mp.nstr(g1, 14), lam_0=mp.nstr(lam0, 12),
                             C_sub=mp.nstr(C, 12), pred_4Fp2_plus_4FFpp=mp.nstr(pred, 12),
                             rel_err=mp.nstr(rel, 6), K=mp.nstr(K, 12),
                             Fp2_over_K=mp.nstr(F1 ** 2 / K, 10),
                             cross_frac=mp.nstr(4 * F0 * F2 / C, 10),
                             F_at_g1_sq=mp.nstr(F0 ** 2, 10),
                             C_over_4K=mp.nstr(C / (4 * K), 10),
                             Fp_gs2_over_K=mp.nstr(G1 ** 2 / K, 10),
                             F_gs_at_g1_sq=mp.nstr(G0 ** 2, 10),
                             minimiser_moved=('yes' if abs(G1 ** 2 - F1 ** 2)
                                              > abs(F1 ** 2) / 10 else 'no')))
        print()
    with open(os.path.join(HERE, 'data', 'phase4c_decomposition.csv'), 'w',
              newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(HERE, 'results', 'phase4c_decomposition.json'), 'w') as fh:
        json.dump(rows, fh, indent=1)
    print('wrote data/phase4c_decomposition.csv')


if __name__ == '__main__':
    main()
