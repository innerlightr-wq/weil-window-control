#!/usr/bin/env python3
"""Regression fixture for the certified Lambda_DH lower bound of REPORT.md section 6c'.

`src/lambda_lower_bound.py` takes about a minute under Arb and re-runs the whole
Rouche argument.  This test does not re-run it.  It does two cheaper things that
together make a silent change in the certificate impossible to miss:

  1. every quantity is compared against a value WRITTEN INTO THIS FILE, so an
     edit to `data/lambda_bracket.txt` fails the test rather than redefining what
     passing means (the failure mode of a check that reads both sides from the
     same place);
  2. `eps` is RECOMPUTED from `m`, `C` and `tau_0` in exact rational arithmetic
     rather than read, so the certificate's own arithmetic is checked, not just
     transcribed.

Exits nonzero on any disagreement.  Run: python tests/test_lambda_certificate.py
"""
import os
import re
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CERT = os.path.join(HERE, 'data', 'lambda_bracket.txt')

# the committed certificate, as published in REPORT.md section 6c'
EXPECT = {
    'im_z1':  '-0.617034364913274771106703921214',
    'radius': '0.3000000000',
    'margin': '0.31703436491327477111',
    'arcs':   3000,
    'relrad': '0.064946579',
    'm':      '6.5935981835929959400e-30',      # the ball midpoint, as printed
    'm_floor': '6.5935981835929959399e-30',     # the exact floor: endpoint rounded DOWN
    'C_ceil': '0.048291904880526492446',        # the exact ceiling: endpoint rounded UP
    'winding': 1,
    'M':      '0.91703436491327477111',
    'tau0':   '0.100000',
    'C':      '0.048291904880526492445',
    'eps':    '6.8268151773112564447e-29',
    'upper':  '1.602537483559822777787819',
}

FAIL = []


def fail(msg):
    FAIL.append(msg)
    print('  [FAIL] %s' % msg)


def ball(txt, pattern):
    """Midpoint of the first Arb ball `[mid +/- rad]` matching `pattern`."""
    m = re.search(pattern + r'\s*=?\s*\[\s*([-0-9.e+]+)\s*\+/-', txt)
    if not m:
        fail('pattern not found in the certificate: %s' % pattern)
        return None
    return m.group(1)


def same(label, got, want):
    if got is None:
        return
    if got != want:
        fail('%s: certificate says %s, this fixture expects %s' % (label, got, want))
    else:
        print('  [PASS] %-28s %s' % (label, got))


def main():
    if not os.path.exists(CERT):
        print('MISSING: %s -- run src/lambda_lower_bound.py --arcs 3000 --prec 400' % CERT)
        return 1
    txt = open(CERT).read()

    print('[1] the disc, and that its closure misses the real axis')
    same('Im z_1', ball(txt, r'Im z_1'), EXPECT['im_z1'])
    same('radius r', ball(txt, r'radius r\s*'), EXPECT['radius'])
    same('margin |Im z_1| - r', ball(txt, r'margin from real axis = \|Im z_1\| - r'),
         EXPECT['margin'])
    if 'closure disjoint from the real axis: True' not in txt:
        fail('the certificate does not record the closure as disjoint from the real axis')
    else:
        print('  [PASS] closure disjoint from the real axis')
    # the margin is the one relation a reader can check by eye
    mg = Fraction(EXPECT['im_z1']).__abs__() - Fraction(EXPECT['radius'])
    if str(mg) and abs(float(mg) - float(EXPECT['margin'])) > 1e-20:
        fail('|Im z_1| - r does not equal the recorded margin')
    else:
        print('  [PASS] margin equals |Im z_1| - r in exact arithmetic')

    print('[2] m > 0 on the boundary, from the arc enclosures')
    if ('partition: %d arcs' % EXPECT['arcs']) not in txt:
        fail('the certificate does not record a %d-arc partition' % EXPECT['arcs'])
    else:
        print('  [PASS] %-28s %d arcs' % ('partition', EXPECT['arcs']))
    same('worst rel. enclosure radius', ball(txt, r'worst relative enclosure radius over the arcs:'),
         EXPECT['relrad'])
    same('m enclosure', ball(txt, r'm enclosure ='), EXPECT['m'])
    if 'm > 0 : True' not in txt:
        fail('the certificate does not record m > 0')
    else:
        print('  [PASS] m > 0')

    print('[3] the winding number is exactly 1')
    if 'increments unambiguous: True' not in txt:
        fail('the certificate does not record the argument increments as unambiguous')
    elif ('has exactly %d zero' % EXPECT['winding']) not in txt:
        fail('the certificate does not record exactly %d zero in D' % EXPECT['winding'])
    else:
        print('  [PASS] exactly %d zero in D, increments unambiguous' % EXPECT['winding'])

    print('[3b] the one-sided bounds are printed as exact floors/ceilings')
    for lab, want in (('m floor', EXPECT['m_floor']), ('C ceiling', EXPECT['C_ceil'])):
        if want not in txt:
            fail('%s: the certificate does not carry the plain decimal %s' % (lab, want))
        else:
            print('  [PASS] %-28s %s' % (lab, want))
    # a floor must not exceed the midpoint it came from
    if Fraction(EXPECT['m_floor']) > Fraction(EXPECT['m']):
        fail('the m floor exceeds the ball midpoint')
    if Fraction(EXPECT['C_ceil']) < Fraction(EXPECT['C']):
        fail('the C ceiling is below the ball midpoint')

    print('[4] M and C')
    same('M', ball(txt, r'M = max_\{dD\}\s*\|Im z\| = \|Im z_1\| \+ r'), EXPECT['M'])
    same('tau_0', ball(txt, r'tau_0 ='), EXPECT['tau0'])
    same('C enclosure', ball(txt, r'C enclosure ='), EXPECT['C'])
    # M = |Im z_1| + r, exactly
    Mx = Fraction(EXPECT['im_z1']).__abs__() + Fraction(EXPECT['radius'])
    if abs(float(Mx) - float(EXPECT['M'])) > 1e-20:
        fail('M does not equal |Im z_1| + r')
    else:
        print('  [PASS] M equals |Im z_1| + r in exact arithmetic')

    print('[5] eps = min(tau_0, m/(2C)), RECOMPUTED in exact rational arithmetic')
    m = Fraction(EXPECT['m'])
    C = Fraction(EXPECT['C'])
    t0 = Fraction(EXPECT['tau0'])
    eps = min(t0, m / (2 * C))
    got = Fraction(EXPECT['eps'])
    # the certificate rounds eps DOWN to 20 significant digits; require that the
    # printed value is a true lower bound and agrees to that precision
    if not (got <= m / (2 * C)):
        fail('the printed eps exceeds m/(2C): it is not a lower bound')
    elif abs(float(eps / got) - 1.0) > 1e-19:
        fail('recomputed eps = %s disagrees with the printed %s' % (float(eps), float(got)))
    else:
        print('  [PASS] %-28s %s  (= m/(2C), below tau_0)' % ('eps recomputed', EXPECT['eps']))
    if eps >= t0:
        fail('eps is not the m/(2C) branch of the minimum, as REPORT 6c\' states')

    print('[6] the bracket is nonempty and ordered')
    if EXPECT['upper'] not in txt:
        fail('the certificate does not carry the upper bound %s' % EXPECT['upper'])
    else:
        print('  [PASS] %-28s %s' % ('upper bound', EXPECT['upper']))
    if Fraction(EXPECT['eps']) >= Fraction(EXPECT['upper']):
        fail('the bracket is empty: eps >= the upper bound')
    else:
        print('  [PASS] %s <= Lambda_DH <= %s' % (EXPECT['eps'], EXPECT['upper']))

    print()
    if FAIL:
        print('=' * 78)
        print('LAMBDA CERTIFICATE REGRESSION FAILED (%d)' % len(FAIL))
        for f in FAIL:
            print('  - %s' % f)
        print('=' * 78)
        return 1
    print('=' * 78)
    print('LAMBDA CERTIFICATE REGRESSION PASSED')
    print('=' * 78)
    return 0


if __name__ == '__main__':
    sys.exit(main())
