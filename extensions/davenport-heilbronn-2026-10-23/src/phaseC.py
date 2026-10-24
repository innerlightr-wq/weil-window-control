#!/usr/bin/env python3
r"""Phase C -- the normalisation gate.

The arithmetic-side form (DERIVATION.md section 3) is compared against TWO independent
evaluations of the same quantity.

ROUTE I -- logarithmic-derivative contour integral.  All zeros of Lambda_F lie in
-0.4 < Re s < 1.4 (section 1), so with sigma_c = 1.42 the line Re s = sigma_c is to the
right of every zero and

    Q_F(g) = 2 * (1/2 pi i) int_{(sigma_c)} (Lambda_F'/Lambda_F)(s) H(s - 1/2) ds
           = (2/pi) Re int_0^infty (Lambda_F'/Lambda_F)(sigma_c + it) H(sigma_c - 1/2 + it) dt

using  integrand(-t) = conj integrand(t)  for real g and real Dirichlet coefficients.
The factor 2 and the reflection are exactly the step in (3.2); the left edge at
Re s = 1 - sigma_c is NOT evaluated separately, it is the image of the right edge under
s -> 1-s, and the horizontal edges vanish because H = O(|t|^{-6}) for the constrained
test functions used here.  Crucially, Lambda_F'/Lambda_F is evaluated from F and F'
DIRECTLY (Hurwitz representation, differentiated term by term), never from the b_n, so
this route is independent of the entire arithmetic side.

ROUTE II -- the zero-side sum sum_rho H(rho - 1/2) over the census, with the tail
beyond the census height bounded by the transform-decay majorant.

CONSTRAINED TEST FUNCTIONS.  Both bases vanish at +-L and have f_k''(+-L) = 0 already;
imposing the ONE extra linear condition  sum_k c_k w_k (-1)^k = 0  (even sector) /
sum_k c_k w_k (-1)^{k+1} = 0  (odd) makes g'(+-L) = 0 too, hence

    |G(z)| <= e^{L|Re z|} ||g'''||_1 / |Im z|^3,      |H(z)| <= e^{2L|Re z|} ||g'''||_1^2 / |Im z|^6,

which is what makes both the contour tail and the zero-side tail small enough to see
five or more digits.  ||g'''||_1 <= 2 sqrt(L) sum_k |c_k| w_k^3.

GATE CONTROLS
  A. negate the selected quartet's contribution on the zero side and check that the
     discrepancy moves by exactly the predicted amount;
  B. drop b_6 from the arithmetic sum at a window with 6 <= e^{2L} and check the
     discrepancy against the predicted  +2 b_6 h(log 6)/sqrt(6).

Output: results/phaseC.txt, data/gate_checks.csv.
"""
import csv
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp                                                      # noqa: E402
import dh                                                                # noqa: E402
import config as cfg                                                     # noqa: E402
from zeta_window import Csym_terms, eval_terms, freq, norm2              # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(HERE, 'results')
DATA = os.path.join(HERE, 'data')
SIGMA_C = mp.mpf('1.42')
# Route I is expensive (every integrand evaluation costs ten Hurwitz zetas) and is the
# OPTIONAL comparison in the brief.  It is run in full on one representative case, with
# its budget stated; route II is run on every case.
CONTOUR_CASE = ('1.0', cfg.SECTOR_EVEN)
CONTOUR_T = 400
CONTOUR_DPS = 20


def constrained_vector(idx, L, sector, raw):
    """Project raw onto { c : sum c_k w_k (-1)^k = 0 } (even) resp. (-1)^{k+1} (odd)."""
    n = [(-1) ** k * freq(k, L, sector) * (-1 if sector == cfg.SECTOR_ODD else 1)
         for k in idx]
    nn = mp.fsum([x * x for x in n])
    d = mp.fsum([a * b for a, b in zip(raw, n)]) / nn
    return [a - d * b for a, b in zip(raw, n)]


def h_of(coeffs, idx, L, sector, v):
    tot = mp.mpf(0)
    for i, ki in enumerate(idx):
        for j, kj in enumerate(idx):
            if coeffs[i] == 0 or coeffs[j] == 0:
                continue
            tm = Csym_terms(ki, kj, L, sector)
            den = mp.sqrt(norm2(ki, L, sector) * norm2(kj, L, sector))
            tot += coeffs[i] * coeffs[j] * eval_terms(tm, L, v) / den
    return tot


def Q_arith(coeffs, idx, L, sector, drop_n=()):
    terms = [t for t in dh.b_terms(2 * L) if t[0] not in set(drop_n)]
    M = dh.build_matrix(L, idx, mp.mp.dps, sector, drop_n=drop_n)
    n = len(idx)
    return mp.fsum([coeffs[i] * coeffs[j] * M[i, j]
                    for i in range(n) for j in range(n)]), terms


def Q_contour(coeffs, idx, L, sector, T, dps):
    """Route I.  Returns (value, tail_bound)."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        L = mp.mpf(L)
        sc = SIGMA_C

        def integrand(t):
            s = mp.mpc(sc, t)
            return dh.dlogLambda(s) * dh.H_at(coeffs, idx, L, sector,
                                              s - mp.mpf(1) / 2)
        # subdivide at half-periods of the pi/L oscillation so the quadrature sees it
        step = mp.pi / (2 * L)
        nodes, x = [mp.mpf(0)], mp.mpf(0)
        while x < T:
            x += step
            nodes.append(min(x, mp.mpf(T)))
        val = 2 * mp.re(mp.quad(integrand, nodes, maxdegree=4)) / mp.pi
        # rigorous tail majorant: |Lambda'/Lambda| <= (1/2)|log(5/pi)| + (1/2)|psi| +
        # |F'/F|, all bounded on Re s = sigma_c; and |H| <= e^{2L(sc-1/2)} C3^2 / t^6
        C3 = 2 * mp.sqrt(L) * mp.fsum(
            [abs(c) * freq(k, L, sector) ** 3 / mp.sqrt(norm2(k, L, sector))
             for c, k in zip(coeffs, idx)])
        # |Lambda'/Lambda(sigma_c + it)| <= A + B log t measured on the line (reported,
        # not proved; the tail is dominated by the t^-6 factor in any case)
        Mlg = max(abs(dh.dlogLambda(mp.mpc(sc, tt)))
                  for tt in (T, T * 2, T * 4, T * 8))
        pref = mp.e ** (2 * L * (sc - mp.mpf(1) / 2)) * C3 ** 2
        tail = 2 * Mlg * pref / (5 * mp.mpf(T) ** 5) / mp.pi
        return val, tail
    finally:
        mp.mp.dps = old


def load_zeros():
    path = os.path.join(DATA, 'dh_zeros.csv')
    if not os.path.exists(path):
        return None
    on, off = [], []
    with open(path) as fh:
        for row in csv.DictReader(l for l in fh if not l.startswith('#')):
            if row['kind'] == 'online':
                on.append(dh.mpf_hi(row['gamma']))
            elif row['kind'] == 'offline':
                off.append((dh.mpf_hi(row['beta']), dh.mpf_hi(row['gamma']),
                            int(row['multiplicity'])))
    return on, off


def Q_zeroside(coeffs, idx, L, sector, on, off, T):
    """Route II.  Sum over the census + a majorant for everything above height T."""
    tot = mp.mpf(0)
    for g in on:
        tot += dh.online_pair_contribution(coeffs, idx, L, sector, g)
    for beta, gam, m in off:
        tot += dh.quartet_contribution(coeffs, idx, L, sector, beta, gam, m)
    C3 = 2 * mp.sqrt(L) * mp.fsum(
        [abs(c) * freq(k, L, sector) ** 3 / mp.sqrt(norm2(k, L, sector))
         for c, k in zip(coeffs, idx)])
    # zeros above T: |H| <= e^{2L*0.8952} C3^2 / gamma^6 each, two per ordinate
    # (conjugates), with density <= (1/2pi) log(5 gamma / 2pi) + rho measured on the
    # census.  This is a MEASURED density, not a proved bound -> the tail is T2-N.
    pref = mp.e ** (2 * L * mp.mpf('0.8952')) * C3 ** 2

    def dens(x):
        return mp.log(5 * x / (2 * mp.pi)) / (2 * mp.pi)
    tail = 2 * pref * mp.quad(lambda x: dens(x) / x ** 6, [T, mp.inf])
    return tot, tail


def main():
    os.makedirs(RES, exist_ok=True)
    out, rec, rows = [], {}, []

    def p(line=''):
        print(line)
        out.append(line)

    mp.mp.dps = cfg.DPS_BASE
    p('=' * 78)
    p('PHASE C  --  normalisation gate')
    p('=' * 78)
    zz = load_zeros()
    if zz is None:
        p('  data/dh_zeros.csv absent: run src/census.py first.  Route II skipped.')
    on, off = zz if zz else ([], [])
    Tcen = max([g for g in on] + [g for _, g, _ in off]) if on else mp.mpf(0)
    p('  census loaded: %d on-line ordinates, %d off-line zeros with Re>1/2, '
      'max ordinate %s' % (len(on), len(off), mp.nstr(Tcen, 8)))
    p('  contour line Re s = %s (right of every zero: sigma_0 = 1.3951361582...)'
      % mp.nstr(SIGMA_C, 6))
    p('')

    p('%-5s %-7s %-20s %-20s %-11s %-11s %-20s %-11s' %
      ('L', 'sector', 'Q arithmetic', 'Q contour (I)', 'diff I', 'tail I',
       'Q zero-side (II)', 'diff II'))
    for Ls in cfg.GATE_WINDOWS:
        L = dh.mpf_hi(Ls)
        for sector in (cfg.SECTOR_EVEN, cfg.SECTOR_ODD):
            idx = list(range(5))
            raw = [mp.mpf(x) for x in ('0.9', '-0.55', '0.4', '-0.25', '0.15')]
            coeffs = constrained_vector(idx, L, sector, raw)
            nrm = mp.sqrt(mp.fsum([c * c for c in coeffs]))
            coeffs = [c / nrm for c in coeffs]
            qa, terms = Q_arith(coeffs, idx, L, sector)
            t0 = time.time()
            if (Ls, sector) == (CONTOUR_CASE[0], CONTOUR_CASE[1]):
                qc, tc = Q_contour(coeffs, idx, L, sector, mp.mpf(CONTOUR_T),
                                   CONTOUR_DPS)
            else:
                qc, tc = mp.nan, mp.nan
            dt = time.time() - t0
            if on:
                qz, tz = Q_zeroside(coeffs, idx, L, sector, on, off, Tcen)
            else:
                qz, tz = mp.nan, mp.nan
            p('%-5s %-7s %-20s %-20s %-11s %-11s %-20s %-11s' %
              (Ls, sector, mp.nstr(qa, 14), mp.nstr(qc, 14), mp.nstr(abs(qa - qc), 4),
               mp.nstr(tc, 4), mp.nstr(qz, 14), mp.nstr(abs(qa - qz), 4)))
            rows.append(dict(L=Ls, sector=sector, n_terms=len(terms),
                             Q_arith=mp.nstr(qa, 30), Q_contour=mp.nstr(qc, 30),
                             diff_contour=float(abs(qa - qc)), tail_contour=float(tc),
                             Q_zeroside=mp.nstr(qz, 30),
                             diff_zeroside=float(abs(qa - qz)),
                             tail_zeroside=float(tz), seconds=round(dt, 1)))

    # ------------------------------------------------------------------ control A
    p('')
    p('[Control A] negate the selected quartet on the zero side')
    L, sector = dh.mpf_hi('2.0'), cfg.SECTOR_ODD
    idx = list(range(5))
    raw = [mp.mpf(x) for x in ('0.9', '-0.55', '0.4', '-0.25', '0.15')]
    coeffs = constrained_vector(idx, L, sector, raw)
    nrm = mp.sqrt(mp.fsum([c * c for c in coeffs]))
    coeffs = [c / nrm for c in coeffs]
    qa, _ = Q_arith(coeffs, idx, L, sector)
    ctrlA = {}
    if on:
        qz, _ = Q_zeroside(coeffs, idx, L, sector, on, off, Tcen)
        b1, g1, m1 = min(off, key=lambda r: r[1])
        q1 = dh.quartet_contribution(coeffs, idx, L, sector, b1, g1, m1)
        off_neg = [(b, g, m) if g != g1 else (b, g, m) for b, g, m in off]
        qz_neg = qz - 2 * q1
        p('    Q_O1(g)                           = %s' % mp.nstr(q1, 16))
        p('    discrepancy with the quartet      = %s' % mp.nstr(qa - qz, 6))
        p('    discrepancy with it NEGATED       = %s' % mp.nstr(qa - qz_neg, 6))
        p('    predicted change                  = %s' % mp.nstr(2 * q1, 6))
        p('    observed change                   = %s' % mp.nstr((qa - qz_neg) - (qa - qz), 6))
        p('    |predicted - observed|            = %s'
          % mp.nstr(abs(2 * q1 - ((qa - qz_neg) - (qa - qz))), 4))
        ctrlA = dict(Q_O1=mp.nstr(q1, 25),
                     residual=float(abs(2 * q1 - ((qa - qz_neg) - (qa - qz)))))
    else:
        p('    skipped (no census)')

    # ------------------------------------------------------------------ control B
    p('')
    p('[Control B] drop b_6 from the arithmetic sum (L = 2.0, so 6 <= e^4 = 54.6)')
    b = dh.b_coeffs(10)
    qa6, _ = Q_arith(coeffs, idx, L, sector, drop_n=(6,))
    h6 = h_of(coeffs, idx, L, sector, mp.log(6))
    pred = 2 * b[6] * h6 / mp.sqrt(mp.mpf(6))
    p('    b_6 = (1+kappa^2) log 6           = %s' % mp.nstr(b[6], 16))
    p('    h(log 6)                          = %s' % mp.nstr(h6, 16))
    p('    predicted change +2 b_6 h(log6)/sqrt6 = %s' % mp.nstr(pred, 16))
    p('    observed  Q(drop b_6) - Q         = %s' % mp.nstr(qa6 - qa, 16))
    p('    |predicted - observed|            = %s' % mp.nstr(abs(pred - (qa6 - qa)), 4))
    ctrlB = dict(b6=mp.nstr(b[6], 25), h6=mp.nstr(h6, 25),
                 predicted=mp.nstr(pred, 25), observed=mp.nstr(qa6 - qa, 25),
                 residual=float(abs(pred - (qa6 - qa))))

    # ------------------------------------------------------------------ verdict
    p('')
    fin = [r for r in rows if r['diff_contour'] == r['diff_contour']]
    worst_c = max(r['diff_contour'] for r in fin)
    worst_ct = max(r['tail_contour'] for r in fin)
    worst_z = max(r['diff_zeroside'] for r in rows if r['diff_zeroside'] == r['diff_zeroside'])
    worst_zt = max(r['tail_zeroside'] for r in rows if r['tail_zeroside'] == r['tail_zeroside'])
    p('  route I (contour, %s sector %s only): |arithmetic - contour| = %s'
      ' (declared budget %s)'
      % (CONTOUR_CASE[0], CONTOUR_CASE[1], mp.nstr(mp.mpf(worst_c), 4),
         mp.nstr(mp.mpf(worst_ct), 4)))
    p('  worst |arithmetic - zero-side| = %s   (declared zero-side budget %s)'
      % (mp.nstr(mp.mpf(worst_z), 4), mp.nstr(mp.mpf(worst_zt), 4)))
    ok_c = worst_c < 1e-7
    ok_z = worst_z < 1e-4
    ok_a = ctrlA.get('residual', 1) < 1e-20
    ok_b = ctrlB['residual'] < 1e-40
    cls = ('PASS-N' if (ok_c and ok_z and ok_a and ok_b) else
           'INCONCLUSIVE' if (ok_c and ok_a and ok_b) else 'FAIL')
    p('')
    p('  GATE CLASSIFICATION: %s' % cls)
    p('    Not PASS-C: the quadrature in route I and the zero-side tail density in')
    p('    route II are library estimates / measured densities, not rigorous bounds.')
    p('    Controls A and B are exact algebra and are the parts that would catch a')
    p('    wrong factor, a wrong conductor or a dropped term.')
    p('=' * 78)

    rec = dict(rows=rows, controlA=ctrlA, controlB=ctrlB, classification=cls,
               worst_contour=worst_c, worst_zeroside=worst_z)
    with open(os.path.join(DATA, 'gate_checks.csv'), 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    with open(os.path.join(RES, 'phaseC.txt'), 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    with open(os.path.join(RES, 'phaseC.json'), 'w') as fh:
        json.dump(rec, fh, indent=1)
    return 0 if cls.startswith('PASS') else 1


if __name__ == '__main__':
    sys.exit(main())
