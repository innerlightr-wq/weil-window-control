#!/usr/bin/env python3
r"""Phase 5: the Siegel corner -- a planted REAL zero at beta = 1/2 + delta, gamma = 0.

A real zero beta = 1/2 + delta forces the partner 1 - beta = 1/2 - delta, so in the
z-variable (rho = 1/2 + i z) the configuration is the pair z = -+ i delta: purely imaginary
ordinates, i.e. exactly gamma = 0.  Its contribution to the window form is

    2 h(i delta) = 2 F(i delta) F(-i delta),

and F(i delta) = int f(x) e^{-delta x} dx, so

    even f :  F(-i delta) = +F(i delta)  ->  + 2 A_delta(f)^2,  A = int f cosh(delta u) du
    odd  f :  F(-i delta) = -F(i delta)  ->  - 2 S_delta(f)^2,  S = int f sinh(delta u) du

The perturbation is therefore RANK ONE and of a DEFINITE SIGN that depends on the parity of
the test function -- positive semidefinite on the even sector, negative semidefinite on the
odd sector.  This is the same mechanism as the pole term of zeta (z = -+ i/2), and
real_block at delta = 1/2 reproduces the parent repository's pole-term matrix exactly,
sector sign included (src/selftest_planting.py, check C).

Consequences measured here:

  EVEN sector.  lambda_min(Q + 2 a a^T) >= lambda_min(Q) for every L, every N, every delta.
  The detector cannot register a real zero at all.  This is stronger than "the delta^2
  coefficient 4K(L,0) vanishes": the response is positive to ALL orders in delta.

  ODD sector.  Q - 2 a a^T is indefinite iff  Theta := 2 a^T Q^{-1} a > 1  (Schur
  complement, valid where Q is positive definite on the subspace).  Theta is the reported
  diagnostic; it is non-decreasing in N on nested bases, so the first N and L at which it
  crosses 1 is the detection onset.  lambda_min and an LDL witness freeze it.

This is a statement about the DETECTOR, not about L(s, chi) (REPORT.md non-claim 2).
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
from zeta_window import build_matrix, EVEN_D, ODD                        # noqa: E402
from weilform import eigsym                                              # noqa: E402
from ldl import ldl, first_negative_pivot, witness, rayleigh             # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def theta(M, a):
    """2 a^T M^{-1} a, or None if M is singular / not usable."""
    try:
        x = mp.lu_solve(M, mp.matrix(a))
    except Exception:                                                   # noqa: BLE001
        return None
    return 2 * sum(a[i] * x[i] for i in range(len(a)))


def one(label, M, L, delta, N, sector):
    idx = range(N)
    B, Ap = pl.real_block(L, idx, delta, sector)
    lam0 = eigsym(M)[0][0]
    Z = M + B
    lam1 = eigsym(Z)[0][0]
    th = theta(M, Ap)
    # one-sided certificate: an explicit trial vector with negative Rayleigh quotient
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
    return dict(label=label, sector=sector, L=mp.nstr(L, 6), N=N,
                delta=mp.nstr(delta, 8),
                lam_unperturbed=mp.nstr(lam0, 12), lam_planted=mp.nstr(lam1, 12),
                delta_lam=mp.nstr(lam1 - lam0, 12),
                response_sign=('>=0' if lam1 >= lam0 else '<0'),
                rank1_eig=mp.nstr((-2 if sector == ODD else 2)
                                  * sum(x ** 2 for x in Ap), 10),
                theta=('' if th is None else mp.nstr(th, 10)),
                theta_gt_1=('' if th is None else str(bool(th > 1))),
                detected=('yes' if lam1 < 0 else 'no'),
                witness_support=wit, min_rayleigh=ray)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--Ls', default='0.5,0.8,1.0,1.3,1.6,2.0')
    ap.add_argument('--deltas', default='0.1,0.25,0.375,0.45,0.49')
    ap.add_argument('--Ns', default='8,12,16,20')
    ap.add_argument('--sectors', default='even_d,odd')
    ap.add_argument('--out', default='phase5_siegel')
    a = ap.parse_args()
    deltas = [mp.mpf(x) for x in a.deltas.split(',')]
    Ns = [int(x) for x in a.Ns.split(',')]
    rows = []
    for Lx in a.Ls.split(','):
        L = mp.mpf(Lx)
        dps = max(60, int(40 + 50 * float(Lx)))
        mp.mp.dps = dps
        print("=" * 118)
        print(f"L = {Lx}  dps = {dps}")
        print("=" * 118)
        for sector in a.sectors.split(','):
            for N in Ns:
                t0 = time.time()
                Mz = build_matrix(L, N, dps, sector)
                mats = [('zeta', Mz)]
                for d in dw.DISCRIMINANTS:
                    mats.append((dw.CHARS[d]['label'],
                                 dw.build_matrix_chi(L, range(N), dps, d, sector)))
                for label, M in mats:
                    for de in deltas:
                        r = one(label, M, L, de, N, sector)
                        rows.append(r)
                        print(f"  {sector:7s} N={N:3d} {label:8s} d={r['delta']:>7s} "
                              f"lam_0={r['lam_unperturbed']:>14s} "
                              f"lam_planted={r['lam_planted']:>14s} "
                              f"resp={r['response_sign']:>4s} "
                              f"Theta={r['theta']:>13s} Th>1={r['theta_gt_1']:>5s} "
                              f"det={r['detected']}")
                print(f"  ({sector}, N={N}: {time.time() - t0:.0f}s)")
        print()
    with open(os.path.join(HERE, 'data', f'{a.out}.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(HERE, 'results', f'{a.out}.json'), 'w') as fh:
        json.dump(rows, fh, indent=1)
    # structural assertion, over the whole table
    bad = [r for r in rows if r['sector'] != ODD and r['response_sign'] != '>=0']
    print(f"EVEN-sector monotonicity violations (must be 0): {len(bad)}")
    if bad:
        print(bad[:5])
    nd = [r for r in rows if r['detected'] == 'yes']
    print(f"rows with a detection (lambda_min < 0): {len(nd)}")
    for r in nd[:40]:
        print('   ', r['label'], r['sector'], 'L=' + r['L'], 'N=' + str(r['N']),
              'delta=' + r['delta'], 'lam=' + r['lam_planted'],
              'Theta=' + r['theta'], 'witness_dim=' + r['witness_support'])
    print(f'wrote data/{a.out}.csv')


if __name__ == '__main__':
    main()
