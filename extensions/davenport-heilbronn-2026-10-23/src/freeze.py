#!/usr/bin/env python3
r"""Choose and freeze the witness, from the scan in data/window_scan.csv.

TWO witnesses are frozen, because they answer two different questions and one of them
cannot carry a certificate.

  certificates/witness.json
      The CERTIFIABLE witness: the smallest scanned window at which the detected
      negativity is MACROSCOPIC, |lambda_min| >= config.MACRO_THRESHOLD, taking the
      deepest such (sector, level) at that window subject to a basis-size cap so that
      the Arb re-computation stays affordable.  This is the one src/certify.py encloses
      and the one every later phase uses.

  certificates/witness_smallest_L.json
      The witness at the SMALLEST scanned window with any detected negative value,
      whatever its size.  At that window the value is many orders of magnitude below
      the scale of the form, so it is reported with a precision ladder and is NOT
      certified.  It is kept because it, not the certifiable one, is the answer to
      "smallest tested window with a detected negative witness".

Neither is "the true first detection threshold" (non-claim 5), and no claim is made
about which sector detects first globally (non-claim 9).

Output: the two JSON files, results/freeze.txt.
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp                                                      # noqa: E402
import dh                                                                # noqa: E402
import config as cfg                                                     # noqa: E402
from phaseD import index_set, rayleigh, gram, GAMMA1, BETA1              # noqa: E402
from zeta_window import freq                                             # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES, DATA, CERT = (os.path.join(HERE, x) for x in ('results', 'data', 'certificates'))
DIG = 30


def build(row, out, p):
    Ls, sector, level = row['L'], row['sector'], int(row['level'])
    L = dh.mpf_hi(Ls)
    idx = index_set(L, sector, level - 1)
    assert len(idx) == int(row['N']), (len(idx), row['N'])
    mp.mp.dps = cfg.DPS_BASE
    M = dh.build_matrix(L, idx, cfg.DPS_BASE, sector)
    E, V = mp.eigsy(M)
    n = len(idx)
    lo = min(range(n), key=lambda i: E[i])
    v = [V[i, lo] for i in range(n)]
    nrm = mp.sqrt(mp.fsum([x * x for x in v]))
    v = [x / nrm for x in v]
    if v[max(range(n), key=lambda i: abs(v[i]))] < 0:
        v = [-x for x in v]
    vr = [mp.mpf(mp.nstr(x, DIG, strip_zeros=False)) for x in v]

    gm = gram(idx, L, sector)
    dev = max(abs(gm[i, j] - (1 if i == j else 0)) for i in range(n) for j in range(n))
    ge = mp.eigsy(gm, eigvals_only=True)
    p('    Gram: max|G - I| = %s, eigenvalues in [%s, %s], condition %s'
      % (mp.nstr(dev, 4), mp.nstr(min(ge), 14), mp.nstr(max(ge), 14),
         mp.nstr(max(ge) / min(ge), 8)))
    p('          (the basis is orthogonal with known norms, so G = I exactly; nothing')
    p('           was regularised and no direction was discarded)')

    vals = {}
    for dps in (cfg.DPS_BASE, cfg.DPS_HIGH, cfg.DPS_XHIGH):
        mp.mp.dps = dps
        Mx = dh.build_matrix(L, idx, dps, sector)
        q = rayleigh(Mx, vr)
        vals[dps] = q
        p('    dps = %-4d  Q_F(g)/||g||^2 = %s' % (dps, mp.nstr(q, 45)))
    mp.mp.dps = cfg.DPS_XHIGH
    p('    |val(%d) - val(%d)| = %-12s   |val(%d) - val(%d)| = %s'
      % (cfg.DPS_BASE, cfg.DPS_HIGH,
         mp.nstr(abs(vals[cfg.DPS_BASE] - vals[cfg.DPS_HIGH]), 4),
         cfg.DPS_HIGH, cfg.DPS_XHIGH,
         mp.nstr(abs(vals[cfg.DPS_HIGH] - vals[cfg.DPS_XHIGH]), 4)))
    p('    eigenvalue reported by the scan                 = %s' % row['lambda_min'])
    p('    |Rayleigh(rationalised v) - scan eigenvalue|    = %s'
      % mp.nstr(abs(vals[cfg.DPS_BASE] - mp.mpf(row['lambda_min'])), 4))
    mp.mp.dps = cfg.DPS_BASE
    return dict(
        L=Ls, sector=sector, level=level, idx=[int(k) for k in idx],
        frequencies=[mp.nstr(freq(k, L, sector), 40) for k in idx],
        coefficients=[mp.nstr(x, DIG) for x in vr],
        basis=('cos((2k+1) pi u / 2L)' if sector != cfg.SECTOR_ODD
               else 'sin((k+1) pi u / L)'),
        normalisation='g = sum_m c_m f_{idx[m]}(u)/sqrt(norm2_m); the reported scalar is'
                      ' Q_F(g)/||g||_2^2 with ||g||_2^2 = sum_m c_m^2',
        value_dps={str(k): mp.nstr(x, 50) for k, x in vals.items()},
        eigenvalue_scan=row['lambda_min'],
        gram_max_dev_from_I=float(dev),
        gram_condition=float(max(ge) / min(ge)),
        gamma1=GAMMA1, beta1=BETA1)


def main():
    mp.mp.dps = cfg.DPS_BASE
    out = []

    def p(line=''):
        print(line)
        out.append(line)

    with open(os.path.join(DATA, 'window_scan.csv')) as fh:
        rows = list(csv.DictReader(fh))
    neg = [r for r in rows if int(r['negative'])]
    p('=' * 92)
    p('FREEZE  --  selecting the witnesses from data/window_scan.csv')
    p('=' * 92)
    p('  %d of %d scanned (window, sector, level) triples have lambda_min < 0'
      % (len(neg), len(rows)))
    if not neg:
        p('  NO negative witness anywhere in the scan.')
        with open(os.path.join(RES, 'freeze.txt'), 'w') as fh:
            fh.write('\n'.join(out) + '\n')
        return 1

    Lmin = min(dh.mpf_hi(r['L']) for r in neg)
    small = [r for r in neg if dh.mpf_hi(r['L']) == Lmin]
    rs = min(small, key=lambda r: mp.mpf(r['lambda_min']))
    p('')
    p('  smallest scanned window with ANY detected negative value: L = %s' % rs['L'])
    for r in small:
        p('      %s level %d (N=%d): lambda_min = %s'
          % (r['sector'], int(r['level']), int(r['N']), r['lambda_min']))
    p('    -> freezing %s level %d' % (rs['sector'], int(rs['level'])))
    p('    CAUTION: |lambda_min| = %s here, against a form whose largest eigenvalue is'
      % mp.nstr(abs(mp.mpf(rs['lambda_min'])), 4))
    p('    of order 10.  The sign is reproducible (see the precision ladder below) but')
    p('    the value is far too small to build a ball certificate around at an')
    p('    affordable precision.  It is reported, not certified.')
    p('    This value is reported, not certified.')
    cs = build(rs, out, p)
    with open(os.path.join(CERT, 'witness_smallest_L.json'), 'w') as fh:
        json.dump(cs, fh, indent=1)
    p('    -> certificates/witness_smallest_L.json')

    macro = [r for r in neg
             if abs(mp.mpf(r['lambda_min'])) >= mp.mpf(cfg.MACRO_THRESHOLD)
             and int(r['N']) <= cfg.CERT_MAX_N]
    p('')
    if not macro:
        p('  no macroscopic witness within the basis-size cap N <= %d' % cfg.CERT_MAX_N)
        with open(os.path.join(RES, 'freeze.txt'), 'w') as fh:
            fh.write('\n'.join(out) + '\n')
        return 1
    Lm = min(dh.mpf_hi(r['L']) for r in macro)
    cand = [r for r in macro if dh.mpf_hi(r['L']) == Lm]
    rc = min(cand, key=lambda r: mp.mpf(r['lambda_min']))
    p('  smallest scanned window with a MACROSCOPIC negative value'
      ' (|lambda_min| >= %s, N <= %d): L = %s' % (cfg.MACRO_THRESHOLD, cfg.CERT_MAX_N, rc['L']))
    for r in cand:
        p('      %s level %d (N=%d): lambda_min = %s'
          % (r['sector'], int(r['level']), int(r['N']), r['lambda_min']))
    p('    -> freezing %s level %d as the CERTIFIABLE witness'
      % (rc['sector'], int(rc['level'])))
    cc = build(rc, out, p)
    with open(os.path.join(CERT, 'witness.json'), 'w') as fh:
        json.dump(cc, fh, indent=1)
    p('    -> certificates/witness.json')
    p('=' * 92)
    with open(os.path.join(RES, 'freeze.txt'), 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
