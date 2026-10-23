#!/usr/bin/env python3
r"""Regression tests for the planting blocks (Phases 4 and 5).

  (A) F_at_imag(k, L, 1/2, sector) == the parent repository's F_at_half, every sector;
  (B) F_at_imag against mp.quad of the defining integral, generic delta;
  (C) real_block(L, ., 1/2, sector) == the parent's POLE-TERM matrix (2 * pole_sign * P_j
      P_k / ||.||), every sector.  The parent's pole term is validated end-to-end by its
      zeta Gate A, so this identity anchors the Phase-5 block -- including its SECTOR SIGN
      -- to already-checked code;
  (D) real_block is rank one, and PSD on the even sectors / NSD on the odd sector, exactly;
  (E) quartet_block at delta = 0 equals 2 * pair_block (the quartet degenerates to a double
      on-line pair);
  (F) real_block equals the gamma -> 0 limit of quartet_block / 2 ... stated precisely:
      quartet_block(gamma=0, delta) == 2 * real_block(delta).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp                                                      # noqa: E402
import planting as pl                                                    # noqa: E402
from zeta_window import EVEN, EVEN_D, ODD, F_at_half, norm2, pole_sign   # noqa: E402

FAIL = []


def check(name, ok, detail=''):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ''))
    if not ok:
        FAIL.append(name)


def main():
    mp.mp.dps = 50
    TOL = mp.mpf(10) ** -40
    Ls = ('0.5', '0.8', '1.3')
    sectors = (EVEN_D, EVEN, ODD)

    print("(A) F_at_imag(., ., 1/2, .) == parent F_at_half")
    worst = mp.mpf(0)
    for sector in sectors:
        for Lx in Ls:
            L = mp.mpf(Lx)
            for k in range(6):
                worst = max(worst, abs(pl.F_at_imag(k, L, mp.mpf(1) / 2, sector)
                                       - F_at_half(k, L, sector)))
    check("max abs diff over 3 sectors x 3 L x 6 modes", worst < TOL,
          f"{mp.nstr(worst, 6)}")

    print("(B) F_at_imag vs quadrature at generic delta")
    worst = mp.mpf(0)
    for sector in sectors:
        L = mp.mpf('0.8')
        for k in range(4):
            w = pl.freq(k, L, sector) if hasattr(pl, 'freq') else None
            for dl in ('0.1', '0.37', '0.49', '-0.37'):
                dd = mp.mpf(dl)
                from zeta_window import freq
                ww = freq(k, L, sector)
                fun = ((lambda x: mp.sin(ww * x) * mp.e ** (-dd * x)) if sector == ODD
                       else (lambda x: mp.cos(ww * x) * mp.e ** (-dd * x)))
                num = mp.quad(fun, [-L, 0, L])
                worst = max(worst, abs(num - pl.F_at_imag(k, L, dd, sector)))
    check("max abs diff vs mp.quad (4 deltas incl. negative)", worst < mp.mpf(10) ** -35,
          f"{mp.nstr(worst, 6)}")

    print("(C) real_block(., ., 1/2, .) == parent pole-term matrix  [anchors the sector sign]")
    worst = mp.mpf(0)
    detail = []
    for sector in sectors:
        for Lx in Ls:
            L = mp.mpf(Lx)
            idx = range(5)
            B, _ = pl.real_block(L, idx, mp.mpf(1) / 2, sector)
            nrm = [mp.sqrt(norm2(k, L, sector)) for k in idx]
            P = [F_at_half(k, L, sector) for k in idx]
            sg = pole_sign(sector)
            for i in range(5):
                for j in range(5):
                    want = 2 * sg * P[i] * P[j] / (nrm[i] * nrm[j])
                    worst = max(worst, abs(B[i, j] - want))
            detail.append(f"{sector}:sign={sg}")
    check("max abs diff vs 2*pole_sign*P_j P_k/||.||", worst < TOL,
          f"{mp.nstr(worst, 6)}  ({', '.join(detail[:3])})")

    print("(D) real_block is rank one, PSD on even sectors and NSD on the odd sector")
    for sector in sectors:
        L = mp.mpf('0.8')
        for dl in ('0.1', '0.25', '0.49'):
            B, Ap = pl.real_block(L, range(6), mp.mpf(dl), sector)
            ev = sorted(mp.eigsy(B)[0], key=lambda z: z)
            nz = [v for v in ev if abs(v) > mp.mpf(10) ** -30]
            nrm2 = sum(x ** 2 for x in Ap)
            want = (-2 if sector == ODD else 2) * nrm2
            ok = (len(nz) == 1 and abs(nz[0] - want) < mp.mpf(10) ** -30
                  and (nz[0] < 0 if sector == ODD else nz[0] > 0))
            check(f"{sector:7s} delta={dl}: one nonzero eigenvalue = "
                  f"{'-' if sector == ODD else '+'}2||A||^2", ok,
                  f"nonzero eigs {[mp.nstr(v, 8) for v in nz]}, 2||A||^2 = "
                  f"{mp.nstr(want, 8)}")

    print("(E) quartet_block(gamma, 0) == 2 * pair_block(gamma)")
    worst = mp.mpf(0)
    for sector in sectors:
        L = mp.mpf('0.8')
        for g in ('6.0', '14.13'):
            Q = pl.quartet_block(L, range(5), mp.mpf(g), mp.mpf(0), sector)
            P2 = pl.pair_block(L, range(5), mp.mpf(g), sector)
            for i in range(5):
                for j in range(5):
                    worst = max(worst, abs(Q[i, j] - 2 * P2[i, j]))
    check("max abs diff", worst < TOL, f"{mp.nstr(worst, 6)}")

    print("(F) quartet_block(0, delta) == 2 * real_block(delta)")
    worst = mp.mpf(0)
    for sector in sectors:
        L = mp.mpf('0.8')
        for dl in ('0.1', '0.37'):
            Q = pl.quartet_block(L, range(5), mp.mpf(0), mp.mpf(dl), sector)
            B, _ = pl.real_block(L, range(5), mp.mpf(dl), sector)
            for i in range(5):
                for j in range(5):
                    worst = max(worst, abs(Q[i, j] - 2 * B[i, j]))
    check("max abs diff", worst < TOL, f"{mp.nstr(worst, 6)}")

    print()
    if FAIL:
        print(f"FAILED: {len(FAIL)}: {FAIL}")
        return 1
    print("ALL PLANTING REGRESSIONS PASSED")
    return 0


if __name__ == '__main__':
    sys.exit(main())
