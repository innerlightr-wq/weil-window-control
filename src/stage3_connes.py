#!/usr/bin/env python3
"""Stage 3b/3c/3d.

3b  Connes's experiment: support [1,13] <-> 2L = log 13 (so exactly the prime powers
    up to 13 enter the geometric side).  Ground state, zeros of its Fourier transform,
    compared with the first 50 zeta zeros.
3c  Ground-state simplicity and the EVEN/ODD gap for L in {0.5, 0.8, 1.0, log(13)/2}.
3d  Recorder split  A = Pole + Arch - LogPi,  2M = Prime,  S = A - 2M; inertia of each.
"""
import sys, os, csv, argparse, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from zeta_window import (build_matrix, build_parts, F, norm2, freq,
                         EVEN, EVEN_D, ODD)
from weilform import eigsym, inertia

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def ground_state(L, N, dps, sector=EVEN_D):
    M = build_matrix(mp.mpf(L), N, dps, sector)
    vals, vecs = eigsym(M)
    return vals, vecs[0], M


def Ffun(c, L, sector):
    nrm = [mp.sqrt(norm2(k, L, sector)) for k in range(len(c))]
    def f(t):
        return sum(c[k] * F(k, L, t, sector) / nrm[k] for k in range(len(c)))
    return f


def real_zeros(fun, tmax, nsamp=4000):
    """Sign-change bracketing + bisection/secant refinement on (0, tmax)."""
    zs = []
    prev_t = mp.mpf(0)
    prev_v = fun(prev_t)
    for i in range(1, nsamp + 1):
        t = mp.mpf(tmax) * i / nsamp
        v = fun(t)
        if prev_v == 0:
            zs.append(prev_t)
        elif v * prev_v < 0:
            try:
                r = mp.findroot(fun, (prev_t, t), solver='anderson',
                                tol=mp.mpf(10) ** (-(mp.mp.dps - 5)))
                zs.append(r)
            except Exception:
                pass
        prev_t, prev_v = t, v
    return zs


def part3b(dps, N, nzeros, rows):
    L = mp.log(13) / 2
    print(f"\n=== 3b: Connes's window, support [1,13]  ->  L = log(13)/2 = {mp.nstr(L,12)} ===")
    print(f"    (2L = log 13, so the geometric side carries exactly the prime powers <= 13)")
    M = build_matrix(L, N, dps, EVEN_D)
    vals, vecs = eigsym(M)
    lam, c = vals[0], vecs[0]
    print(f"    N={N}, dps={dps}: lam_min = {mp.nstr(lam,12)}  lam_2 = {mp.nstr(vals[1],8)}"
          f"  simple = {vals[1] > vals[0] * 1000 or vals[1] - vals[0] > mp.mpf(10)**(-(dps//2))}")
    g = [mp.im(mp.zetazero(n)) for n in range(1, nzeros + 1)]
    fun = Ffun(c, L, EVEN_D)
    zs = real_zeros(fun, g[-1] * mp.mpf('1.02'), nsamp=6000)
    print(f"    ground-state F has {len(zs)} real zeros in (0, {mp.nstr(g[-1]*mp.mpf('1.02'),6)})")
    print(f"    {'n':>3s} {'gamma_n (zeta)':>22s} {'nearest F-zero':>22s} {'abs err':>12s} {'rel err':>11s}")
    for n in range(1, nzeros + 1):
        gn = g[n - 1]
        if not zs:
            break
        z = min(zs, key=lambda z: abs(z - gn))
        err = abs(z - gn)
        if n <= 12 or n % 5 == 0 or n == nzeros:
            print(f"    {n:3d} {mp.nstr(gn,16):>22s} {mp.nstr(z,16):>22s} "
                  f"{mp.nstr(err,6):>12s} {mp.nstr(err/gn,5):>11s}")
        rows.append(dict(part='3b', L=mp.nstr(L, 12), N=N, dps=dps, n=n,
                         gamma_n=mp.nstr(gn, 20), F_zero=mp.nstr(z, 20),
                         abs_err=mp.nstr(err, 8), rel_err=mp.nstr(err / gn, 8)))
    return L, lam


def part3c(Ls, N, dps, rows):
    print(f"\n=== 3c: ground-state simplicity and the EVEN/ODD gap (N={N}, dps={dps}) ===")
    print(f"    {'L':>12s} {'lam_min even':>18s} {'lam_2 even':>14s} {'lam_min odd':>18s}"
          f" {'odd - even':>16s} {'global GS':>10s}")
    for Lx in Ls:
        L = mp.mpf(Lx) if not isinstance(Lx, str) else mp.mpf(Lx)
        ve, _ = eigsym(build_matrix(L, N, dps, EVEN_D))
        vo, _ = eigsym(build_matrix(L, N, dps, ODD))
        gap = vo[0] - ve[0]
        which = 'even' if ve[0] <= vo[0] else 'ODD'
        print(f"    {mp.nstr(L,8):>12s} {mp.nstr(ve[0],10):>18s} {mp.nstr(ve[1],8):>14s}"
              f" {mp.nstr(vo[0],10):>18s} {mp.nstr(gap,8):>16s} {which:>10s}")
        rows.append(dict(part='3c', L=mp.nstr(L, 10), N=N, dps=dps,
                         lam_min_even=mp.nstr(ve[0], 12), lam_2_even=mp.nstr(ve[1], 10),
                         lam_min_odd=mp.nstr(vo[0], 12), lam_2_odd=mp.nstr(vo[1], 10),
                         even_odd_gap=mp.nstr(gap, 10), global_ground_state=which,
                         even_simple=bool(ve[1] - ve[0] > abs(ve[0]) * 1000)))


def part3d(Ls, N, dps, rows):
    print(f"\n=== 3d: recorder split  A = Pole + Arch - LogPi,  2M = Prime,  S = A - 2M ===")
    print(f"    {'L':>10s} {'sector':>8s} {'inertia A':>14s} {'inertia M':>14s}"
          f" {'inertia S':>14s} {'A indef':>8s} {'M indef':>8s} {'S PSD':>7s}")
    for Lx in Ls:
        L = mp.mpf(Lx)
        for sector in (EVEN_D, ODD):
            Pole, Arch, LogPi, Prime, _ = build_parts(L, N, dps, sector)
            A = Pole + Arch - LogPi
            Mr = Prime / 2
            S = A - 2 * Mr
            iA, _, _ = inertia(A); iM, _, _ = inertia(Mr); iS, _, _ = inertia(S)
            resid = max(abs((A - 2 * Mr - S)[i, j]) for i in range(N) for j in range(N))
            print(f"    {mp.nstr(L,6):>10s} {sector:>8s} {str(iA):>14s} {str(iM):>14s}"
                  f" {str(iS):>14s} {str(iA[0]>0 and iA[2]>0):>8s}"
                  f" {str(iM[0]>0 and iM[2]>0):>8s} {str(iS[0]==0):>7s}")
            rows.append(dict(part='3d', L=mp.nstr(L, 10), sector=sector, N=N, dps=dps,
                             A_inertia=str(iA), M_inertia=str(iM), S_inertia=str(iS),
                             A_indefinite=(iA[0] > 0 and iA[2] > 0),
                             M_indefinite=(iM[0] > 0 and iM[2] > 0),
                             S_psd=(iS[0] == 0), split_residual=mp.nstr(resid, 4)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--N', type=int, default=16)
    ap.add_argument('--dps', type=int, default=50)
    ap.add_argument('--nzeros', type=int, default=50)
    ap.add_argument('--parts', default='b,c,d')
    a = ap.parse_args()
    mp.mp.dps = a.dps
    rows3b, rows3c, rows3d = [], [], []
    Ls = ['0.5', '0.8', '1.0', mp.nstr(mp.log(13) / 2, 20)]
    if 'b' in a.parts:
        part3b(a.dps, a.N, a.nzeros, rows3b)
    if 'c' in a.parts:
        part3c(Ls, a.N, a.dps, rows3c)
    if 'd' in a.parts:
        part3d(Ls, a.N, a.dps, rows3d)
    for name, data in (('stage3b_connes_zeros.csv', rows3b),
                       ('stage3c_even_odd.csv', rows3c),
                       ('stage3d_recorder_split.csv', rows3d)):
        if data:
            with open(os.path.join(ROOT, 'data', name), 'w', newline='') as fh:
                w = csv.DictWriter(fh, fieldnames=list(data[0].keys()))
                w.writeheader(); w.writerows(data)
            print(f"wrote data/{name}")


if __name__ == '__main__':
    main()
