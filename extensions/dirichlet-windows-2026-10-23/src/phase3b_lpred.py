#!/usr/bin/env python3
r"""Phase 3b: does the CONDUCTOR enter the resolution scale like HEIGHT?

The parent repository's detection onset tracks L_pred = (1/2) log(gamma_*/2pi), read off
"T*(L) = 2 pi e^{2L} first reaches gamma_*".  For L(s, chi) the zero-counting density is
(1/2pi) log(q gamma / 2pi), i.e. the ANALYTIC CONDUCTOR is q gamma, so the candidate is

    L_pred^(q) = (1/2) log(q gamma_1(chi) / 2 pi)   vs   L_pred^(1) = (1/2) log(gamma_1/2pi).

Operationalisation (the parent's, src/stage3e_certify.py): the window "resolves" gamma_1
when the MOVED-PAIR counterfactual at gamma_1 first makes the form indefinite.  We bisect
in L for the onset at fixed delta and fixed nested basis, then compare the two predictors.

The two predictors differ by (1/2) log q, which over q = 3..8 spans 0.549..1.040 -- a range
that happens to DOMINATE the spread of (1/2) log(gamma_1/2pi) across these five characters
(gamma_1 in roughly 6..9 gives -0.02..0.18).  So the test is well posed even on five
conductors, although the FIT is not: five points and two one-parameter predictors is T3 for
any coefficient, T2 only for the residuals reported.
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
from weilform import eigsym                                             # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def dps_for(L):
    return max(60, int(40 + 50 * float(L)))


def indefinite(builder, L, g1, delta, N):
    dps = dps_for(L)
    mp.mp.dps = dps
    M = builder(L, N, dps)
    lam0 = eigsym(M)[0][0]
    P = pl.pair_block(L, range(N), g1, EVEN_D)
    Q = pl.quartet_block(L, range(N), g1, delta, EVEN_D)
    lam = eigsym(M - P + Q)[0][0]
    floor = abs(lam0) if lam0 < 0 else mp.mpf(0)
    trust = (lam < 0) and (abs(lam) > 100 * floor if floor > 0 else True)
    return trust, lam, lam0


def onset(builder, g1, delta, N, lo=mp.mpf('0.3'), hi=mp.mpf('1.2'),
          tol=mp.mpf('0.005')):
    # The bracket is capped at 1.2 deliberately.  Phase 4 shows that at delta = 0.375
    # every one of the six functions is already detected at L = 1.0 and that chi_-7 and
    # chi_8 are not at L = 0.8, so every onset lies inside (0.3, 1.0].  Probing up to
    # L = 3.2 costs ~300 s per evaluation (dps 200, ~140 prime powers) and buys nothing;
    # a `None` return would mean the cap was wrong and is reported as such.
    ok_hi, lam_hi, _ = indefinite(builder, hi, g1, delta, N)
    if not ok_hi:
        return None, lam_hi
    ok_lo, _, _ = indefinite(builder, lo, g1, delta, N)
    if ok_lo:
        return lo, None                      # already detected at the bottom of the range
    while hi - lo > tol:
        mid = (lo + hi) / 2
        ok, _, _ = indefinite(builder, mid, g1, delta, N)
        if ok:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2, None


def fit(xs, ys):
    """least-squares y = alpha + beta x; returns (alpha, beta, max |resid|, rms)."""
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    beta = sxy / sxx
    alpha = my - beta * mx
    res = [y - (alpha + beta * x) for x, y in zip(xs, ys)]
    return alpha, beta, max(abs(r) for r in res), mp.sqrt(sum(r ** 2 for r in res) / n), res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--delta', default='0.375')
    ap.add_argument('--N', type=int, default=16)
    ap.add_argument('--out', default='phase3b_lpred')
    a = ap.parse_args()
    delta = mp.mpf(a.delta)
    rows = []
    mp.mp.dps = 40
    jobs = [('zeta', 1, 0, mp.im(mp.zetazero(1)),
             lambda L, N, dps: build_matrix(L, N, dps, EVEN_D))]
    for d in dw.DISCRIMINANTS:
        lab = dw.CHARS[d]['label']
        gs = [mp.mpf(x) for x in open(os.path.join(HERE, 'data', f'zeros_{lab}.txt'))
              if not x.startswith('#')]
        jobs.append((lab, dw.CHARS[d]['q'], dw.CHARS[d]['a'], gs[0],
                     (lambda dd: (lambda L, N, dps:
                                  dw.build_matrix_chi(L, range(N), dps, dd, EVEN_D)))(d)))
    print(f"delta = {a.delta}, N = {a.N}, sector = even_d, moved-pair counterfactual\n")
    print(f"  {'function':8s} {'q':>2s} {'a':>2s} {'gamma_1':>14s} {'L_onset':>9s} "
          f"{'Lpred(1)':>9s} {'Lpred(q)':>9s} {'onset-Lp1':>10s} {'onset-Lpq':>10s}")
    for lab, q, apar, g1, builder in jobs:
        t0 = time.time()
        Lon, why = onset(builder, g1, delta, a.N)
        lp1 = mp.log(g1 / (2 * mp.pi)) / 2
        lpq = mp.log(q * g1 / (2 * mp.pi)) / 2
        if Lon is None:
            print(f"  {lab:8s} {q:>2d} {apar:>2d} {mp.nstr(g1, 12):>14s} "
                  f"{'none<=3.2':>9s}")
            rows.append(dict(label=lab, q=q, a_par=apar, gamma_1=mp.nstr(g1, 14),
                             delta=a.delta, N=a.N, L_onset='',
                             L_pred_1=mp.nstr(lp1, 8), L_pred_q=mp.nstr(lpq, 8),
                             seconds=round(time.time() - t0, 1)))
            continue
        print(f"  {lab:8s} {q:>2d} {apar:>2d} {mp.nstr(g1, 12):>14s} "
              f"{mp.nstr(Lon, 6):>9s} {mp.nstr(lp1, 6):>9s} {mp.nstr(lpq, 6):>9s} "
              f"{mp.nstr(Lon - lp1, 5):>10s} {mp.nstr(Lon - lpq, 5):>10s} "
              f"({time.time() - t0:.0f}s)")
        rows.append(dict(label=lab, q=q, a_par=apar, gamma_1=mp.nstr(g1, 14),
                         delta=a.delta, N=a.N, L_onset=mp.nstr(Lon, 8),
                         L_pred_1=mp.nstr(lp1, 8), L_pred_q=mp.nstr(lpq, 8),
                         resid_1=mp.nstr(Lon - lp1, 8), resid_q=mp.nstr(Lon - lpq, 8),
                         seconds=round(time.time() - t0, 1)))
    chi_rows = [r for r in rows if r['label'] != 'zeta' and r['L_onset']]
    summary = {}
    if len(chi_rows) >= 3:
        ys = [mp.mpf(r['L_onset']) for r in chi_rows]
        for key, name in (('L_pred_1', 'Lpred(1) = 1/2 log(g1/2pi)'),
                          ('L_pred_q', 'Lpred(q) = 1/2 log(q g1/2pi)')):
            xs = [mp.mpf(r[key]) for r in chi_rows]
            al, be, mx, rms, res = fit(xs, ys)
            print(f"\n  fit L_onset = alpha + beta * {name}")
            print(f"    alpha = {mp.nstr(al, 6)}   beta = {mp.nstr(be, 6)}   "
                  f"max|resid| = {mp.nstr(mx, 4)}   rms = {mp.nstr(rms, 4)}")
            print("    residuals: " + ", ".join(
                f"{r['label']}:{mp.nstr(v, 4)}" for r, v in zip(chi_rows, res)))
            summary[key] = dict(alpha=mp.nstr(al, 8), beta=mp.nstr(be, 8),
                                max_resid=mp.nstr(mx, 6), rms=mp.nstr(rms, 6),
                                residuals={r['label']: mp.nstr(v, 6)
                                           for r, v in zip(chi_rows, res)})
        # and the no-free-parameter versions: is onset - Lpred constant?
        for key in ('resid_1', 'resid_q'):
            vs = [mp.mpf(r[key]) for r in chi_rows]
            m = sum(vs) / len(vs)
            sd = mp.sqrt(sum((v - m) ** 2 for v in vs) / len(vs))
            print(f"  offset-only, {key}: mean {mp.nstr(m, 6)}  sd {mp.nstr(sd, 4)}  "
                  f"spread {mp.nstr(max(vs) - min(vs), 4)}")
            summary[key] = dict(mean=mp.nstr(m, 8), sd=mp.nstr(sd, 6),
                                spread=mp.nstr(max(vs) - min(vs), 6))
    with open(os.path.join(HERE, 'data', f'{a.out}.csv'), 'w', newline='') as fh:
        keys = sorted({k for r in rows for k in r})
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(HERE, 'results', f'{a.out}.json'), 'w') as fh:
        json.dump(dict(rows=rows, fits=summary), fh, indent=1)
    print(f'\nwrote data/{a.out}.csv')


if __name__ == '__main__':
    main()
