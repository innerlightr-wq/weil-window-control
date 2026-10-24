#!/usr/bin/env python3
r"""Guard against RETRACTIONS.md R6: every decimal CONSTANT must be parsed at a fixed
high precision, never at the ambient mp.mp.dps.

Two checks:

  1. a behavioural one -- mp.mpf('1.8') built at dps 15 and at dps 160 differ in the
     17th digit, while dh.mpf_hi('1.8') does not depend on the ambient precision;
  2. a source one -- no file in src/ constructs a window length, a located zero or a
     perturbation parameter with a bare mp.mpf(<decimal string or CSV field>).
"""
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), 'src'))
import mpmath as mp                                                      # noqa: E402
import dh                                                                # noqa: E402

SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src')

# names that hold a decimal constant parsed from a string or a CSV field
# Only the constants that FEED THE FORM are precision-critical: a window length, the
# located zero, a perturbation parameter.  Contour coordinates, search step sizes and
# comparison tolerances are not -- a 1e-17 change in a tolerance or in the abscissa of a
# provably zero-free contour changes no reported value -- and are deliberately NOT
# flagged.
SUSPECT = re.compile(r"mp\.mpf\(\s*(Ls|c\['L'\]|row\['beta'\]"
                     r"|row\['gamma'\]|ts|BETA1|GAMMA1)\s*\)")
ALLOW = {'dh.py'}          # dh.py defines mpf_hi and may use mp.mpf internally

FAIL = []


def main():
    print('precision-constant guard (RETRACTIONS.md R6)')

    mp.mp.dps = 15
    lo = mp.mpf('1.8')
    hi_lo = dh.mpf_hi('1.8')
    mp.mp.dps = 200
    hi = mp.mpf('1.8')
    hi_hi = dh.mpf_hi('1.8')
    d_bare = abs(lo - hi)
    d_fix = abs(hi_lo - hi_hi)
    print('  |mp.mpf("1.8") at dps 15  -  at dps 200| = %s   (this is the bug)'
          % mp.nstr(d_bare, 6))
    print('  |dh.mpf_hi("1.8") at dps 15 - at dps 200| = %s   (ambient-independent)'
          % mp.nstr(d_fix, 6))
    if not (d_bare > mp.mpf('1e-20')):
        FAIL.append('the control did not reproduce the bare-mpf discrepancy')
    if d_fix != 0:
        FAIL.append('dh.mpf_hi is not ambient-independent')
    # and binary-exact decimals are unaffected either way
    mp.mp.dps = 15
    a = mp.mpf('2.0')
    mp.mp.dps = 200
    b = mp.mpf('2.0')
    print('  |mp.mpf("2.0") at dps 15  -  at dps 200| = %s   (binary-exact: unaffected)'
          % mp.nstr(abs(a - b), 6))
    if a != b:
        FAIL.append('2.0 is not binary-exact after all')

    print('  source scan:')
    for fn in sorted(os.listdir(SRC)):
        if not fn.endswith('.py') or fn in ALLOW:
            continue
        with open(os.path.join(SRC, fn)) as fh:
            for i, line in enumerate(fh, 1):
                if line.lstrip().startswith('#'):
                    continue
                m = SUSPECT.search(line)
                if m:
                    print('    %s:%d  %s' % (fn, i, line.strip()[:84]))
                    FAIL.append('%s:%d' % (fn, i))
    if not any(f for f in FAIL if ':' in f):
        print('    clean')

    print('')
    if FAIL:
        print('FAIL: %s' % FAIL)
        return 1
    print('PASS')
    return 0


if __name__ == '__main__':
    sys.exit(main())
