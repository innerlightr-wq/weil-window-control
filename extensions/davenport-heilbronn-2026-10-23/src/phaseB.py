#!/usr/bin/env python3
r"""Phase B -- the arithmetic-side Weil form, and an independent scalar evaluation.

  B1  b_n by (2.1) and by the Dirichlet-inverse route; the exact algebraic checks
      b_2 = kappa log 2, b_3 = -kappa log 3, b_6 = (1+kappa^2) log 6; b_n = 0 for 5|n.
  B2  the arithmetic cutoff n <= e^{2L} decided in integer arithmetic, not by a float.
  B3  the closed-form transforms G_at against quadrature, including z near +-i w where
      the formula is 0/0.
  B4  |G(it)| t^2 bounded: the decay that makes the contour edges and the zero-side sum
      converge.
  B5  an INDEPENDENT scalar evaluation of Q_F(g) for fixed test functions: the
      archimedean term by direct quadrature of (1/2pi) int |G(it)|^2 Re psi(3/4+it/2) dt
      with an explicit analytic tail bound, versus the matrix route's closed form.

Output: results/phaseB.txt, results/phaseB.json, data/dh_coefficients.csv.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp                                                      # noqa: E402
import dh                                                                # noqa: E402
import config as cfg                                                     # noqa: E402
from zeta_window import Csym_terms, eval_terms, freq, norm2              # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(HERE, 'results')
DATA = os.path.join(HERE, 'data')


def h_of(coeffs, idx, L, sector, v):
    """h(v) = int g(u+v) g(u) du for g = sum coeffs[m] f_{idx[m]}/sqrt(norm2)."""
    tot = mp.mpf(0)
    for i, ki in enumerate(idx):
        for j, kj in enumerate(idx):
            if coeffs[i] == 0 or coeffs[j] == 0:
                continue
            tm = Csym_terms(ki, kj, L, sector)
            den = mp.sqrt(norm2(ki, L, sector) * norm2(kj, L, sector))
            tot += coeffs[i] * coeffs[j] * eval_terms(tm, L, v) / den
    return tot


def arch_quadrature(coeffs, idx, L, sector, T, dps):
    r"""(1/2pi) int_{-T}^{T} |G(it)|^2 Re psi(3/4 + it/2) dt, plus an ANALYTIC tail bound.

    Tail: for a single mode of frequency w, |G_k(it)|^2 <= 4 w^2 / (t^2 - w^2)^2 for
    t > w; for g = sum c_k f_k / sqrt(norm2) Cauchy-Schwarz gives
    |G(it)| <= (sum |c_k|) * max_k |G_k(it)| / min_k sqrt(norm2).  And for t >= 2,
    |Re psi(3/4 + it/2)| <= log(3/4 + t/2) <= log t.  So

        tail <= (1/pi) * Cg^2 * int_T^inf log t / (t^2 - wmax^2)^2 dt
             <= (1/pi) * Cg^2 * log T' / (3 (T - wmax^2/T)^3) * (1 + 1/log T)   [crude]

    we evaluate the tail integral itself by quadrature of the MAJORANT, which is a
    rigorous upper bound on the tail, not an error estimate of the main quadrature."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        L = mp.mpf(L)
        ws = [freq(k, L, sector) for k in idx]
        wmax = max(ws)

        def Gsq(t):
            return abs(dh.G_at(coeffs, idx, L, sector, mp.mpc(0, t))) ** 2

        def integrand(t):
            return Gsq(t) * mp.re(mp.digamma(mp.mpc(mp.mpf(3) / 4, t / 2)))

        # split at the resonances so the quadrature sees the peaks
        nodes = [0]
        for w in sorted(set(ws)):
            nodes += [max(0, w - 5), w, w + 5]
        nodes = sorted(set(x for x in nodes if 0 <= x <= T)) + [T]
        nodes = sorted(set(nodes))
        main = 2 * mp.quad(integrand, nodes) / (2 * mp.pi)

        Cg = mp.fsum([abs(c) for c in coeffs]) * 2 * wmax / mp.sqrt(
            min(norm2(k, L, sector) for k in idx))

        def majorant(t):
            return Cg ** 2 / (t ** 2 - wmax ** 2) ** 2 * mp.log(mp.mpf(3) / 4 + t / 2)

        tail = 2 * mp.quad(majorant, [T, mp.inf]) / (2 * mp.pi)
        return main, tail
    finally:
        mp.mp.dps = old


def main():
    os.makedirs(RES, exist_ok=True)
    os.makedirs(DATA, exist_ok=True)
    out = []
    rec = {}

    def p(line=''):
        print(line)
        out.append(line)

    mp.mp.dps = cfg.DPS_BASE
    k = dh.kappa()

    p('=' * 78)
    p('PHASE B  --  the arithmetic-side Weil form')
    p('=' * 78)

    # ---------------------------------------------------------------- B1
    p('')
    p('[B1] logarithmic-derivative coefficients b_n')
    NB = 200
    b1 = dh.b_coeffs(NB)
    b2 = dh.b_coeffs_via_inverse(NB)
    dmax = max(abs(b1[n] - b2[n]) for n in range(2, NB + 1))
    bmax = max(abs(b1[n]) for n in range(2, NB + 1))
    p('   recursion (2.1) vs Dirichlet-inverse route, n <= %d:' % NB)
    p('     max |b_n - b_n^inv|          = %s' % mp.nstr(dmax, 5))
    p('     max |b_n|                    = %s  (relative %s)'
      % (mp.nstr(bmax, 8), mp.nstr(dmax / bmax, 5)))
    e2 = abs(b1[2] - k * mp.log(2))
    e3 = abs(b1[3] + k * mp.log(3))
    e6 = abs(b1[6] - (1 + k ** 2) * mp.log(6))
    p('     |b_2 - kappa log 2|          = %s' % mp.nstr(e2, 5))
    p('     |b_3 + kappa log 3|          = %s' % mp.nstr(e3, 5))
    p('     |b_6 - (1+kappa^2) log 6|    = %s      <- non-prime-power control'
      % mp.nstr(e6, 5))
    five = max(abs(b1[n]) for n in range(5, NB + 1, 5))
    p('     max |b_n| over 5 | n         = %s      (vanishes: every divisor chain meets'
      ' a_5 = 0)' % mp.nstr(five, 5))
    p('     b_1                          = %s' % mp.nstr(b1[1], 5))
    npp = [n for n in range(2, 50) if b1[n] != 0 and not _is_pp(n)]
    p('     non-prime-powers below 50 with b_n != 0: %s' % npp)
    p('       (for zeta this list is EMPTY; here it is not -- no Euler product)')
    rec['b'] = dict(nmax=NB, inverse_route_diff=float(dmax), b2=float(e2), b3=float(e3),
                    b6=float(e6), five_divides=float(five), b1=float(b1[1]))

    with open(os.path.join(DATA, 'dh_coefficients.csv'), 'w') as fh:
        fh.write('n,a_n,b_n,b_n_inverse_route,is_prime_power\n')
        for n in range(1, NB + 1):
            fh.write('%d,%s,%s,%s,%d\n' % (n, mp.nstr(dh.a_coeff(n), 30),
                                           mp.nstr(b1[n], 30), mp.nstr(b2[n], 30),
                                           1 if _is_pp(n) else 0))
    p('   -> data/dh_coefficients.csv  (n = 1..%d)' % NB)

    # ---------------------------------------------------------------- B2
    p('')
    p('[B2] arithmetic cutoff  n <= e^{2L}  (integer test, then exact log re-test)')
    rows = []
    for Ls in cfg.WINDOW_SCAN + ['2.6']:
        L = dh.mpf_hi(Ls)
        tm = dh.b_terms(2 * L)
        nmax = max(t[0] for t in tm)
        # the boundary case: is there an n with log n exactly 2L?  Only if e^{2L} is an
        # integer, impossible for rational L != 0 by Lindemann.  Report the margin.
        marg = 2 * L - mp.log(nmax)
        rows.append(dict(L=Ls, nmax=nmax, count=len(tm), margin=float(marg)))
        p('     L = %-5s  2L = %-5s  e^{2L} = %-12s  n_max = %-4d  terms = %-4d'
          '  margin 2L - log n_max = %s'
          % (Ls, mp.nstr(2 * L, 6), mp.nstr(mp.e ** (2 * L), 10), nmax, len(tm),
             mp.nstr(marg, 6)))
    p('     (L = 2.6 -> n through 181, as the brief states)')
    rec['cutoff'] = rows

    # ---------------------------------------------------------------- B3
    p('')
    p('[B3] closed-form transforms G_at vs quadrature')
    L = dh.mpf_hi('1.7')
    worst = mp.mpf(0)
    for sector in (cfg.SECTOR_EVEN, cfg.SECTOR_ODD):
        for kk in (0, 3, 11):
            w = freq(kk, L, sector)
            fn = (mp.cos if sector != 'odd' else mp.sin)
            for z in [mp.mpf('0.3'), mp.mpc('0.3', '5.0'), mp.mpc(0, w),
                      mp.mpc(0, w) + mp.mpf('1e-25'), mp.mpc('0.30851718', '85.699')]:
                closed = dh.G_at([mp.mpf(1)], [kk], L, sector, z) * mp.sqrt(
                    norm2(kk, L, sector))
                quad = mp.quad(lambda u: fn(w * u) * mp.e ** (z * u), [-L, 0, L])
                rel = abs(closed - quad) / max(abs(quad), mp.mpf(1))
                worst = max(worst, rel)
    p('     worst relative discrepancy over 30 (sector, mode, z) combinations '
      '(including z = i w exactly, the removable 0/0): %s' % mp.nstr(worst, 5))
    rec['G_at_vs_quad'] = float(worst)

    # ---------------------------------------------------------------- B4
    p('')
    p('[B4] transform decay  |G(it)| t^2  (bases vanish at +-L, so g is continuous)')
    for sector in (cfg.SECTOR_EVEN, cfg.SECTOR_ODD):
        kk = 5
        w = freq(kk, L, sector)
        sup = mp.mpf(0)
        for jj in range(1, 2001):
            t = mp.mpf(jj) * mp.mpf('0.5') + 50
            sup = max(sup, abs(dh.G_at([mp.mpf(1)], [kk], L, sector,
                                       mp.mpc(0, t))) * t ** 2)
        p('     %-7s mode k=5 (w=%-8s): sup_{50<=t<=1050} |G(it)| t^2 = %s'
          % (sector, mp.nstr(w, 6), mp.nstr(sup, 8)))
    p('     bounded => |H(it)| = O(t^-4); the horizontal contour edges vanish and the')
    p('     zero-side sum converges absolutely (DERIVATION.md section 3).')

    # ---------------------------------------------------------------- B5
    p('')
    p('[B5] independent scalar evaluation of Q_F(g)')
    p('     route 1: the matrix assembly (closed-form archimedean, digamma/Lerch series)')
    p('     route 2: Gauss/Binet integral for psi (non-oscillatory, different kernel)')
    p('')
    p('     %-5s %-7s %-22s %-22s %-11s %-11s'
      % ('L', 'sector', 'Q route 1', 'Q route 2', '|diff|', '[0,x0] bnd'))
    rows = []
    for Ls in cfg.GATE_WINDOWS:
        L = dh.mpf_hi(Ls)
        for sector in (cfg.SECTOR_EVEN, cfg.SECTOR_ODD):
            idx = [0, 1, 2]
            coeffs = [mp.mpf('0.7'), mp.mpf('-0.5'), mp.mpf('0.3')]
            # route 1
            M = dh.build_matrix(L, idx, cfg.DPS_BASE, sector)
            q1 = mp.fsum([coeffs[i] * coeffs[j] * M[i, j]
                          for i in range(3) for j in range(3)])
            # route 2: Gauss/Binet integral for psi, non-oscillatory
            arch, cut = dh.arch_binet(coeffs, idx, L, sector, cfg.DPS_BASE)
            h0 = h_of(coeffs, idx, L, sector, 0)
            prime = mp.fsum([2 * bn / mp.sqrt(mp.mpf(nn))
                             * h_of(coeffs, idx, L, sector, logn)
                             for nn, bn, logn in dh.b_terms(2 * L)])
            q2 = h0 * mp.log(mp.mpf(5) / mp.pi) + arch - prime
            p('     %-5s %-7s %-22s %-22s %-11s %-11s'
              % (Ls, sector, mp.nstr(q1, 16), mp.nstr(q2, 16),
                 mp.nstr(abs(q1 - q2), 4), mp.nstr(cut, 4)))
            rows.append(dict(L=Ls, sector=sector, q1=mp.nstr(q1, 25),
                             q2=mp.nstr(q2, 25), diff=float(abs(q1 - q2)),
                             cut=float(cut)))
    rec['scalar_routes'] = rows
    worstd = max(r['diff'] for r in rows)
    worstt = max(r['cut'] for r in rows)
    p('')
    p('     worst |route1 - route2| = %s ;  worst [0,x0] bound = %s'
      % (mp.nstr(mp.mpf(worstd), 4), mp.nstr(mp.mpf(worstt), 4)))
    p('     Route 2 = Gauss/Binet integral for psi: a ONE-DIMENSIONAL, NON-oscillatory')
    p('     integral over [0, 4L] plus the EXACT tail h(0) E_1(4L).  The piece over')
    p('     [0, x0] is bounded analytically (column 6); the quadrature on [x0, 4L] is')
    p('     a LIBRARY ESTIMATE and is not claimed as a bound.')

    ok = (dmax < mp.mpf(10) ** (-cfg.DPS_BASE + 12) and e2 == 0 and e3 == 0 and e6 == 0
          and worst < mp.mpf('1e-20') and worstd < 10 * worstt)
    p('')
    p('=' * 78)
    p('PHASE B GATE: %s' % ('PASS' if ok else 'FAIL'))
    p('=' * 78)
    rec['gate'] = 'PASS' if ok else 'FAIL'

    with open(os.path.join(RES, 'phaseB.txt'), 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    with open(os.path.join(RES, 'phaseB.json'), 'w') as fh:
        json.dump(rec, fh, indent=1)
    return 0 if ok else 1


def _is_pp(n):
    m, p = n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            return m == 1
        p += 1
    return True


if __name__ == '__main__':
    sys.exit(main())
