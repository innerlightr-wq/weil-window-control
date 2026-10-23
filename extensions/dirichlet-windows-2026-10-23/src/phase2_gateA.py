#!/usr/bin/env python3
r"""Phase 2 Gate A: assembled geometric side vs sum over zeros, for every character.

For each chi and each L the assembled matrix M_chi (D1) is compared, entry by entry,
against

    Z_K[j,k] = sum_{n<=K} 2 F_j(gamma_n) F_k(gamma_n) / (||f_j|| ||f_k||),

the factor 2 accounting for the ordinate pair +-gamma_n (the zeros of L(s,chi) for real chi
are symmetric about the real axis).  The discrepancy must DECAY at the basis-predicted rate
gamma_K^{-3} (times a log); the script fits the exponent.

Three controls, all reported whether or not they pass:

  * zeta control: the parent repository's own Gate A, recomputed in the same run with the
    same code path, so a sector or sign bug in the chi code shows up as a chi/zeta
    disagreement rather than as a plausible-looking number;
  * sign-flip on the CONDUCTOR term alone (LogQ -> -LogQ);
  * sign-flip on the ARCHIMEDEAN term alone (Arch -> -Arch).
If a flip does not blow the discrepancy up, Gate A is reported as UNINFORMATIVE for that
term, not as passed.
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
from zeta_window import F, build_matrix, norm2, EVEN, EVEN_D, ODD        # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def zeros_side(L, sector, N, gammas, Ks):
    """max-norm of (geometric - zeros) at each K, cumulative."""
    L = mp.mpf(L)
    nrm = [mp.sqrt(norm2(k, L, sector)) for k in range(N)]
    Z = mp.zeros(N, N)
    out = {}
    done = 0
    for K in Ks:
        for g in gammas[done:K]:
            Fv = [F(k, L, g, sector) / nrm[k] for k in range(N)]
            for j in range(N):
                for k in range(N):
                    Z[j, k] += 2 * Fv[j] * Fv[k]
        done = K
        out[K] = mp.matrix(Z)
    return out


def maxdiff(M, Z, N):
    return max(abs(M[j, k] - Z[j, k]) for j in range(N) for k in range(N))


def fit_exponent(xs, ys):
    """least-squares slope of log y vs log x."""
    n = len(xs)
    lx = [mp.log(x) for x in xs]
    ly = [mp.log(y) for y in ys]
    mx, my = sum(lx) / n, sum(ly) / n
    num = sum((a - mx) * (b - my) for a, b in zip(lx, ly))
    den = sum((a - mx) ** 2 for a in lx)
    return num / den


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--Ls', default='0.5,0.8,1.0')
    ap.add_argument('--N', type=int, default=4)
    ap.add_argument('--Ks', default='20,40,80,160')
    ap.add_argument('--dps', type=int, default=40)
    ap.add_argument('--sectors', default='even_d,even,odd')
    ap.add_argument('--zeta-zeros', type=int, default=160)
    a = ap.parse_args()
    mp.mp.dps = a.dps
    Ks = [int(x) for x in a.Ks.split(',')]
    Ls = [mp.mpf(x) for x in a.Ls.split(',')]
    sectors = a.sectors.split(',')
    rows, summary = [], {}

    # ---------------------------------------------------------------- zeta control
    print("=" * 96)
    print("ZETA CONTROL (parent repository Gate A, same run, same code path)")
    print("=" * 96)
    t0 = time.time()
    zg = [mp.im(mp.zetazero(n)) for n in range(1, a.zeta_zeros + 1)]
    print(f"  {len(zg)} zeta zeros, gamma_last = {mp.nstr(zg[-1], 10)} "
          f"({time.time() - t0:.0f}s)")
    for sector in sectors:
        for L in Ls:
            M = build_matrix(L, a.N, a.dps, sector)
            Zs = zeros_side(L, sector, a.N, zg, Ks)
            ds = [maxdiff(M, Zs[K], a.N) for K in Ks]
            ex = fit_exponent([zg[K - 1] for K in Ks], ds)
            print(f"  zeta   sector={sector:7s} L={mp.nstr(L, 4):>5s} "
                  + "  ".join(f"K={K}:{mp.nstr(d, 4)}" for K, d in zip(Ks, ds))
                  + f"   fitted exponent {mp.nstr(ex, 5)}")
            rows.append(dict(which='zeta', d='', q=1, a_par=0, label='zeta',
                             sector=sector, L=mp.nstr(L, 6), N=a.N, flip='none',
                             **{f'disc_K{K}': mp.nstr(d, 8) for K, d in zip(Ks, ds)},
                             gamma_Klast=mp.nstr(zg[Ks[-1] - 1], 10),
                             fit_exponent=mp.nstr(ex, 8)))
            summary.setdefault('zeta', {})[f'{sector}|{mp.nstr(L, 4)}'] = dict(
                disc=[mp.nstr(x, 8) for x in ds], exponent=mp.nstr(ex, 6))

    # ---------------------------------------------------------------- characters
    for d in dw.DISCRIMINANTS:
        lab, q, apar = dw.CHARS[d]['label'], dw.CHARS[d]['q'], dw.CHARS[d]['a']
        path = os.path.join(HERE, 'data', f'zeros_{lab}.txt')
        gs = [mp.mpf(x) for x in open(path) if not x.startswith('#')]
        print("=" * 96)
        print(f"{lab}  (d={d}, q={q}, a={apar});  {len(gs)} ordinates, "
              f"gamma_1={mp.nstr(gs[0], 12)}, gamma_last={mp.nstr(gs[-1], 10)}")
        print("=" * 96)
        myKs = [K for K in Ks if K <= len(gs)]
        for sector in sectors:
            for L in Ls:
                Zs = zeros_side(L, sector, a.N, gs, myKs)
                line = {}
                for flip in (None, 'logq', 'arch'):
                    M = dw.build_matrix_chi(L, range(a.N), a.dps, d, sector, flip=flip)
                    ds = [maxdiff(M, Zs[K], a.N) for K in myKs]
                    ex = fit_exponent([gs[K - 1] for K in myKs], ds)
                    line[flip or 'none'] = (ds, ex)
                    rows.append(dict(which='chi', d=d, q=q, a_par=apar, label=lab,
                                     sector=sector, L=mp.nstr(L, 6), N=a.N,
                                     flip=(flip or 'none'),
                                     **{f'disc_K{K}': mp.nstr(dd, 8)
                                        for K, dd in zip(myKs, ds)},
                                     gamma_Klast=mp.nstr(gs[myKs[-1] - 1], 10),
                                     fit_exponent=mp.nstr(ex, 8)))
                ds0, ex0 = line['none']
                dq, _ = line['logq']
                da, _ = line['arch']
                ratio_q = dq[-1] / ds0[-1]
                ratio_a = da[-1] / ds0[-1]
                print(f"  {sector:7s} L={mp.nstr(L, 4):>5s}  "
                      + "  ".join(f"K={K}:{mp.nstr(x, 4)}" for K, x in zip(myKs, ds0))
                      + f"  | exp {mp.nstr(ex0, 5):>8s}"
                      + f"  | flip logq x{mp.nstr(ratio_q, 4):>9s}"
                      + f"  flip arch x{mp.nstr(ratio_a, 4):>9s}")
                summary.setdefault(lab, {})[f'{sector}|{mp.nstr(L, 4)}'] = dict(
                    disc=[mp.nstr(x, 8) for x in ds0], exponent=mp.nstr(ex0, 6),
                    flip_logq_ratio=mp.nstr(ratio_q, 6),
                    flip_arch_ratio=mp.nstr(ratio_a, 6))

    with open(os.path.join(HERE, 'data', 'phase2_gateA.csv'), 'w', newline='') as fh:
        keys = sorted({k for r in rows for k in r})
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(HERE, 'results', 'phase2_gateA.json'), 'w') as fh:
        json.dump(summary, fh, indent=1)
    print("\nwrote data/phase2_gateA.csv and results/phase2_gateA.json")


if __name__ == '__main__':
    main()
