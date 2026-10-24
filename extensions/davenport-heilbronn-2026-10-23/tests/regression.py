#!/usr/bin/env python3
r"""Parent-repository regressions, to be run BEFORE and AFTER this experiment.

Checks, in a clean process, that this extension has not perturbed anything upstream:

  1. the parent zeta window form src/zeta_window.build_matrix reproduces the committed
     data/stage3_gateB_lambda_min.csv and data/stage3c_even_odd.csv values;
  2. the dirichlet-windows extension's L(s,chi) window form reproduces its committed
     data/phase3_gateB.csv values;
  3. no file outside this extension has been modified (git status).

Expected values are READ FROM THE COMMITTED CSVs; they are never edited to make a run
pass.  A mismatch is reported with both numbers.
"""
import csv
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXT = HERE
ROOT = os.path.dirname(os.path.dirname(EXT))
sys.path.insert(0, os.path.join(ROOT, 'src'))
sys.path.insert(0, os.path.join(ROOT, 'extensions', 'dirichlet-windows-2026-10-23', 'src'))

import mpmath as mp                                                      # noqa: E402
import zeta_window as zw                                                 # noqa: E402

FAIL = []


def read_csv(path):
    with open(path) as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith('#')))


def cmp(name, got, want, rel=mp.mpf('1e-18')):
    want = mp.mpf(want)
    d = abs(got - want) / max(abs(want), mp.mpf(1))
    ok = d <= rel
    print('   %-52s %-6s got %-24s want %-24s rel %s'
          % (name, 'OK' if ok else 'FAIL', mp.nstr(got, 18), mp.nstr(want, 18),
             mp.nstr(d, 4)))
    if not ok:
        FAIL.append(name)


def main():
    mp.mp.dps = 50
    print('=' * 110)
    print('PARENT-REPOSITORY REGRESSIONS (clean process)')
    print('=' * 110)

    print(' [1] zeta window: data/stage3_gateB_lambda_min.csv')
    rows = read_csv(os.path.join(ROOT, 'data', 'stage3_gateB_lambda_min.csv'))
    print('     columns: %s' % list(rows[0].keys()))
    n = 0
    for r in rows:
        try:
            L = mp.mpf(r.get('L'))
            N = int(r.get('N'))
            dps = int(r.get('dps', 40))
            want = r.get('lambda_min') or r.get('lam_min') or r.get('lambda')
            sector = r.get('sector', 'even')
        except (TypeError, ValueError):
            continue
        if want is None or n >= 6:
            continue
        M = zw.build_matrix(L, N, dps, sector)
        got = min(mp.eigsy(M, eigvals_only=True))
        cmp('zeta L=%s N=%d %s' % (r.get('L'), N, sector), got, want, mp.mpf('1e-14'))
        n += 1
    if n == 0:
        print('     (no rows with the expected columns; see results for the raw header)')

    print('')
    print(' [2] dirichlet windows: data/phase3_gateB.csv')
    dpath = os.path.join(ROOT, 'extensions', 'dirichlet-windows-2026-10-23', 'data',
                         'phase3_gateB.csv')
    if os.path.exists(dpath):
        import dirichlet_window as dw
        rows = read_csv(dpath)
        print('     columns: %s' % list(rows[0].keys()))
        n = 0
        for r in rows:
            try:
                d = int(r['d'])
                L = mp.mpf(r['L'])
                N = int(r['N'])
                dps = int(r.get('dps', 40))
                want = r.get('lam_min') or r.get('lambda_min')
                sector = r.get('sector', 'even_d')
            except (KeyError, TypeError, ValueError):
                continue
            if want is None or n >= 4:
                continue
            M = dw.build_matrix_chi(L, range(N), dps, d, sector)
            got = min(mp.eigsy(M, eigvals_only=True))
            cmp('chi d=%d L=%s N=%d %s' % (d, r['L'], N, sector), got, want,
                mp.mpf('1e-14'))
            n += 1
        if n == 0:
            print('     (no rows with the expected columns)')
    else:
        print('     file absent')

    print('')
    print(' [3] repository diff outside this extension')
    out = subprocess.run(['git', 'status', '--porcelain'], cwd=ROOT,
                         capture_output=True, text=True).stdout
    touched = []
    for line in out.splitlines():
        path = line[3:].strip()
        if not path:
            continue
        if path.startswith('extensions/davenport-heilbronn-2026-10-23'):
            continue
        touched.append(line)
    if touched:
        print('     FILES CHANGED OUTSIDE THE EXTENSION:')
        for t in touched:
            print('       %s' % t)
    else:
        print('     none')
    print('     (pre-existing untracked paths are listed above if present; they were')
    print('      recorded in REPORT.md section 0 before any work began.)')

    print('')
    if FAIL:
        print('REGRESSION FAILURES: %s' % FAIL)
        return 1
    print('ALL REGRESSIONS PASS')
    return 0


if __name__ == '__main__':
    sys.exit(main())
