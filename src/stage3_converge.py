#!/usr/bin/env python3
"""Gate B, done properly: lambda*(L) from TWO independent complete even bases.

  EVEN   : cos(k pi x/L)        -- contains the constant, jumps at +-L
  EVEN_D : cos((2k+1)pi x/2L)   -- vanishes at +-L
Both are complete in L^2_even[-L,L], so lambda_min(N) decreases to the same lambda*(L).
Agreement between them is the cross-validation; Zhu's interval is the external check.
"""
import sys, os, csv, argparse, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from zeta_window import build_matrix, F, norm2, EVEN, EVEN_D, ODD
from weilform import eigsym

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ZHU_LO, ZHU_HI = mp.mpf('8.9e-18'), mp.mpf('2.27e-17')


def gateA(L, sector, N, Ks, allg, dps):
    M = build_matrix(mp.mpf(L), N, dps, sector)
    nrm = [mp.sqrt(norm2(k, mp.mpf(L), sector)) for k in range(N)]
    out = []
    Z = mp.zeros(N, N)
    done = 0
    for K in Ks:
        for g in allg[done:K]:
            Fv = [F(k, mp.mpf(L), g, sector) / nrm[k] for k in range(N)]
            for j in range(N):
                for k in range(N):
                    Z[j, k] += 2 * Fv[j] * Fv[k]
        done = K
        d = max(abs(M[j, k] - Z[j, k]) for j in range(N) for k in range(N))
        out.append((K, allg[K - 1], d))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', default='0.8')
    ap.add_argument('--N', default='4,6,8,10,12,14,16,18,20')
    ap.add_argument('--dps', type=int, default=50)
    ap.add_argument('--gateA-K', default='50,200,800')
    ap.add_argument('--zeros', type=int, default=800)
    a = ap.parse_args()
    L = mp.mpf(a.L)
    Ns = [int(x) for x in a.N.split(',')]
    Ks = [int(x) for x in a.gateA_K.split(',')]
    mp.mp.dps = a.dps

    print(f"L = {a.L}, dps = {a.dps};  Zhu certified lambda*(0.8) in [{ZHU_LO}, {ZHU_HI}]")
    print("\n--- GATE A: geometric side vs zeros side (K zeros) ---")
    t0 = time.time()
    allg = [mp.im(mp.zetazero(n)) for n in range(1, a.zeros + 1)]
    print(f"    ({a.zeros} zeros in {time.time()-t0:.0f}s)")
    arows = []
    for sector in (EVEN, EVEN_D, ODD):
        for K, g, d in gateA(a.L, sector, 4, Ks, allg, 30):
            print(f"    sector={sector:7s} K={K:4d} gamma_K={mp.nstr(g,8):>10s} "
                  f"max|geom-zeros| = {mp.nstr(d,6)}")
            arows.append(dict(L=a.L, sector=sector, N=4, K=K, gamma_K=mp.nstr(g, 10),
                              max_abs_diff=mp.nstr(d, 8)))

    print("\n--- GATE B: lambda_min(N) from two independent even bases ---")
    brows = []
    res = {}
    for sector in (EVEN_D, EVEN):
        prev = None
        for N in Ns:
            t0 = time.time()
            M = build_matrix(L, N, a.dps, sector)
            vals, vecs = eigsym(M)
            lam = vals[0]
            mono = '-' if prev is None else ('OK' if lam <= prev + mp.mpf(10) ** (-(a.dps - 8))
                                             else 'VIOLATED')
            inz = bool(ZHU_LO <= lam <= ZHU_HI)
            print(f"    {sector:7s} N={N:3d} lam_min={mp.nstr(lam,10):>16s} "
                  f"lam_2={mp.nstr(vals[1],8):>12s} inZhu={str(inz):5s} mono={mono:8s} "
                  f"({time.time()-t0:.0f}s)")
            brows.append(dict(L=a.L, sector=sector, N=N, dps=a.dps,
                              lam_min=mp.nstr(lam, 12),
                              lam_min_sign=('pos' if lam > 0 else 'neg'),
                              lam_2=mp.nstr(vals[1], 10), in_zhu_interval=inz,
                              monotone=mono))
            res.setdefault(sector, []).append((N, lam))
            prev = lam
    for name, rows_ in (('stage3_gateA_explicit_formula.csv', arows),
                        ('stage3_gateB_lambda_min.csv', brows)):
        with open(os.path.join(ROOT, 'data', name), 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows_[0].keys()))
            w.writeheader(); w.writerows(rows_)
        print(f"wrote data/{name}")


if __name__ == '__main__':
    main()
