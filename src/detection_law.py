#!/usr/bin/env python3
"""Closed form for lam_min of the real off-line pair, and the detection threshold.

T1 derivation.  For the real functional-equation pair {rho, 1/rho} the window data is
t(n) = rho^n + rho^-n, so with u_i = rho^i, w_i = rho^-i  (i = 0..R),

    T_R = u w^T + w u^T .

On span{u, w} this is the 2x2 pencil with Gram G = [[|u|^2, u.w],[u.w, |w|^2]];
the two nonzero eigenvalues are  lam_+- = (u.w) +- |u| |w| .  Since u.w = R+1,

    (EXACT)   lam_min(T_R) = (R+1) - |u| |w| ,
              |u|^2 = (rho^{2R+2} - 1)/(rho^2 - 1),  |w|^2 = (rho^{-2R-2} - 1)/(rho^-2 - 1).

By Cauchy-Schwarz |u||w| >= u.w = R+1 with equality iff u || w iff rho = 1.  So
lam_min < 0 for EVERY R >= 1 as soon as rho != 1: the first window that can see
anything already sees the violation.  Asymptotically |u||w| ~ rho^{R+2}/(rho^2-1), so
lam_min ~ -rho^{R+2}/(rho^2-1): exponential divergence at rate rho.
"""
import sys, os, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from weilform import toeplitz, eigsym

mp.mp.dps = 80


def lam_min_closed(rho, R):
    rho = mp.mpf(rho)
    if rho == 1:
        return mp.mpf(0)
    u2 = (rho ** (2 * R + 2) - 1) / (rho ** 2 - 1)
    w2 = (rho ** (-2 * R - 2) - 1) / (rho ** -2 - 1)
    return (R + 1) - mp.sqrt(u2 * w2)


def verify():
    print("verify closed form against direct eigensolve (dps=80):")
    worst = mp.mpf(0)
    for rho in ['3/2', '13/10', '11/10', '101/100']:
        for R in range(1, 15):
            t = [mp.mpf(rho) ** n + mp.mpf(rho) ** (-n) for n in range(R + 1)]
            direct = eigsym(toeplitz(t, R))[0][0]
            closed = lam_min_closed(rho, R)
            worst = max(worst, abs(direct - closed) / max(abs(direct), mp.mpf(1)))
    print(f"  max relative discrepancy over rho in (1.5,1.3,1.1,1.01), R=1..14: "
          f"{mp.nstr(worst, 4)}")
    return worst


def threshold_table(path):
    """For each off-line distance eps (rho = 1+eps), the smallest R with
    |lam_min| above a given noise floor."""
    floors = [mp.mpf(10) ** -16, mp.mpf(10) ** -17, mp.mpf(10) ** -30, mp.mpf(10) ** -50]
    rows = []
    print("\nR needed for |lam_min| to clear a noise floor, rho = 1 + eps")
    print(f"{'eps':>10s} {'|lam|(R=1)':>14s} {'|lam|(R=10)':>14s} {'|lam|(R=100)':>14s}"
          + ''.join(f"{'R>'+mp.nstr(f,2):>14s}" for f in floors))
    for k in range(1, 13):
        eps = mp.mpf(10) ** (-k)
        rho = 1 + eps
        need = []
        for f in floors:
            # |lam_min| is strictly increasing in R, so binary search is exact
            lo, hi = 1, 1
            while hi <= 2 ** 22 and abs(lam_min_closed(rho, hi)) <= f:
                lo, hi = hi, hi * 2
            found = None
            if abs(lam_min_closed(rho, hi)) > f:
                while lo < hi:
                    mid = (lo + hi) // 2
                    if abs(lam_min_closed(rho, mid)) > f:
                        hi = mid
                    else:
                        lo = mid + 1
                found = lo
            need.append(found)
        rows.append(dict(eps=mp.nstr(eps, 3), rho=mp.nstr(rho, 14),
                         lam_R1=mp.nstr(lam_min_closed(rho, 1), 6),
                         lam_R10=mp.nstr(lam_min_closed(rho, 10), 6),
                         lam_R100=mp.nstr(lam_min_closed(rho, 100), 6),
                         **{f'R_for_floor_1e-{int(-mp.log10(f))}': need[i]
                            for i, f in enumerate(floors)}))
        print(f"{mp.nstr(eps,3):>10s} {mp.nstr(abs(lam_min_closed(rho,1)),6):>14s}"
              f" {mp.nstr(abs(lam_min_closed(rho,10)),6):>14s}"
              f" {mp.nstr(abs(lam_min_closed(rho,100)),6):>14s}"
              + ''.join(f"{str(n):>14s}" for n in need))
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"wrote {path}")
    # small-eps law at fixed R
    print("\nsmall-eps scaling at fixed R:  lam_min ~ -C_R eps^2 ?")
    for R in (1, 2, 3, 4, 5, 8, 16, 32):
        e1, e2 = mp.mpf(10) ** -6, mp.mpf(10) ** -7
        l1, l2 = abs(lam_min_closed(1 + e1, R)), abs(lam_min_closed(1 + e2, R))
        p = mp.log(l1 / l2) / mp.log(e1 / e2)
        binom = mp.mpf(R * (R + 1) * (R + 2)) / 6
        print(f"  R={R:3d}  C_R = {mp.nstr(l1/e1**2, 10):>16s}   binom(R+2,3) = {binom}"
              f"   ratio = {mp.nstr((l1/e1**2)/binom, 10)}   exponent = {mp.nstr(p, 8)}")


if __name__ == '__main__':
    verify()
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       'results', 'stage2_detection_law.csv')
    threshold_table(out)
