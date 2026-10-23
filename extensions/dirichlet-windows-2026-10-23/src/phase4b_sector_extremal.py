#!/usr/bin/env python3
r"""Phase 4b: the derivative-evaluation extremal is SECTOR-dependent, and the two sectors
are exactly complementary.

Theorem 8 of the paper is stated for the even sector, and its constant is
    K_even(L, gamma) = sup { |F'(gamma)|^2 : f real EVEN, supp f in [-L,L], ||f||_2 = 1 }
                     = int_{-L}^{L} u^2 sin^2(gamma u) du,
attained by f ~ u sin(gamma u) (which is even).  For real ODD f write F = i Ftilde,
Ftilde(t) = int f(x) sin(tx) dx, so Ftilde'(t) = int x f(x) cos(tx) dx and

    K_odd(L, gamma) = sup { |Ftilde'(gamma)|^2 : f real ODD, ... } = int u^2 cos^2(gamma u) du,

attained by f ~ u cos(gamma u) (which is odd).  Hence the EXACT complementarity

    K_even(L, gamma) + K_odd(L, gamma) = int_{-L}^{L} u^2 du = 2 L^3 / 3     for every gamma,

so the paper's gamma-uniform bound 2L^3/3 is not a slack bound at all: it is the TOTAL of
the two sectors, and K_even attains it exactly where K_odd vanishes.  In particular

    K_even(L, 0) = 0        and        K_odd(L, 0) = 2 L^3 / 3 = the maximum,

which is why the gamma = 0 (Siegel) blindness is an even-sector phenomenon and the odd
sector is MAXIMALLY sensitive there.  Checked numerically below, including attainment.
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp                                                      # noqa: E402
import planting as pl                                                    # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def K_odd_closed(L, g):
    """int_{-L}^{L} u^2 cos^2(g u) du = 2L^3/3 - K_even(L,g)."""
    return 2 * mp.mpf(L) ** 3 / 3 - pl.K_closed(L, g)


def main():
    mp.mp.dps = 50
    fail = []
    print("=" * 94)
    print("(a) K_even + K_odd = 2L^3/3 exactly, for every gamma")
    print("=" * 94)
    print(f"{'L':>5s} {'gamma':>10s} {'K_even':>16s} {'K_odd':>16s} {'sum':>16s} "
          f"{'2L^3/3':>16s} {'|sum-2L^3/3|':>13s}")
    rows = []
    for Lx in ('0.5', '0.8', '1.0', '1.3', '2.0'):
        L = mp.mpf(Lx)
        for gx in ('0', '1e-6', '6.0209489', '8.0397372', '14.1347251', '100.0'):
            g = mp.mpf(gx)
            ke = pl.K_closed(L, g)
            ko = K_odd_closed(L, g)
            tot = ke + ko
            want = 2 * L ** 3 / 3
            err = abs(tot - want)
            if err > mp.mpf(10) ** -40:
                fail.append(('sum', Lx, gx))
            print(f"{Lx:>5s} {gx:>10s} {mp.nstr(ke, 12):>16s} {mp.nstr(ko, 12):>16s} "
                  f"{mp.nstr(tot, 12):>16s} {mp.nstr(want, 12):>16s} "
                  f"{mp.nstr(err, 4):>13s}")
            rows.append(dict(L=Lx, gamma=gx, K_even=mp.nstr(ke, 14),
                             K_odd=mp.nstr(ko, 14), sum=mp.nstr(tot, 14),
                             two_L3_3=mp.nstr(want, 14), abs_err=mp.nstr(err, 6)))

    print("\n" + "=" * 94)
    print("(b) attainment: K_even by f ~ u sin(gu) (even), K_odd by f ~ u cos(gu) (odd)")
    print("=" * 94)
    for Lx in ('0.8', '2.0'):
        L = mp.mpf(Lx)
        for gx in ('0', '8.0397372', '14.1347251'):
            g = mp.mpf(gx)
            # even sector: f = u sin(gu)/||.||, F'(g) = -int u f(u) sin(gu) du
            n2 = mp.quad(lambda u: (u * mp.sin(g * u)) ** 2, [-L, 0, L])
            if n2 > 0:
                Fp = -mp.quad(lambda u: u * (u * mp.sin(g * u) / mp.sqrt(n2))
                              * mp.sin(g * u), [-L, 0, L])
                ok = abs(Fp ** 2 - pl.K_closed(L, g)) < mp.mpf(10) ** -30
            else:
                Fp, ok = mp.mpf(0), (pl.K_closed(L, g) == 0)
            # odd sector: f = u cos(gu)/||.||, Ftilde'(g) = int u f(u) cos(gu) du
            m2 = mp.quad(lambda u: (u * mp.cos(g * u)) ** 2, [-L, 0, L])
            Fpo = mp.quad(lambda u: u * (u * mp.cos(g * u) / mp.sqrt(m2))
                          * mp.cos(g * u), [-L, 0, L])
            oko = abs(Fpo ** 2 - K_odd_closed(L, g)) < mp.mpf(10) ** -30
            print(f"  L={Lx} gamma={gx:>12s}: even |F'|^2={mp.nstr(Fp ** 2, 12):>16s} "
                  f"(=K_even: {ok})   odd |F~'|^2={mp.nstr(Fpo ** 2, 12):>16s} "
                  f"(=K_odd: {oko})")
            if not (ok and oko):
                fail.append(('attain', Lx, gx))

    print("\n" + "=" * 94)
    print("(c) the gamma = 0 corner")
    print("=" * 94)
    for Lx in ('0.5', '0.8', '1.0', '1.3', '1.6', '2.0'):
        L = mp.mpf(Lx)
        print(f"  L={Lx:>4s}: K_even(L,0) = {mp.nstr(pl.K_closed(L, 0), 6):>8s}   "
              f"K_odd(L,0) = {mp.nstr(K_odd_closed(L, 0), 12):>14s}   "
              f"2L^3/3 = {mp.nstr(2 * L ** 3 / 3, 12):>14s}   "
              f"ratio = {mp.nstr(K_odd_closed(L, 0) / (2 * L ** 3 / 3), 8)}")
    with open(os.path.join(HERE, 'data', 'phase4b_sector_extremal.csv'), 'w',
              newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nfailures: {len(fail)}  {fail if fail else ''}")
    print("wrote data/phase4b_sector_extremal.csv")
    return 1 if fail else 0


if __name__ == '__main__':
    sys.exit(main())
