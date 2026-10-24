#!/usr/bin/env python3
r"""Gate A, term-level sign-flip control (exact algebra).

The window matrix is  M_DH = Arch - LogPi + LogQ - Prime  (DERIVATION.md section 3).
Negating exactly one term must move Q(g) by exactly twice that term's contribution, with
the sign fixed by how the term enters.  Writing T(g) = c^T T c for each part,

    flip Arch  :  Q_flip - Q = -2 Arch(g)
    flip LogQ  :  Q_flip - Q = -2 LogQ(g)      [ LogQ(g) = log(5) h(0) ]
    flip Prime :  Q_flip - Q = +2 Prime(g)

These are identities of the assembly, so the residuals are at the precision floor; what
they test is that each term is present, enters with the stated sign, and carries the
stated conductor.  The LogQ row in particular is the control that would catch a wrong
conductor: it must equal -2 log(5) h(0) and nothing else.

This complements the two controls already in src/phaseC.py: control A negates the
selected quartet on the ZERO side (a sign flip of a spectral contribution), control B
omits b_6 from the arithmetic sum.  The present script flips the ARITHMETIC-side terms.

Output: data/gate_signflip.csv, results/gate_signflip.txt.
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp                                                      # noqa: E402
import dh                                                                # noqa: E402
import config as cfg                                                     # noqa: E402
from zeta_window import Csym_terms, eval_terms, norm2                    # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES, DATA = os.path.join(HERE, 'results'), os.path.join(HERE, 'data')


def quad(M, v):
    n = len(v)
    return mp.fsum([v[i] * v[j] * M[i, j] for i in range(n) for j in range(n)])


def h_at(coeffs, idx, L, sector, u):
    tot = mp.mpf(0)
    for i, ki in enumerate(idx):
        for j, kj in enumerate(idx):
            tm = Csym_terms(ki, kj, L, sector)
            den = mp.sqrt(norm2(ki, L, sector) * norm2(kj, L, sector))
            tot += coeffs[i] * coeffs[j] * eval_terms(tm, L, u) / den
    return tot


def main():
    mp.mp.dps = cfg.DPS_BASE
    out, rows = [], []

    def p(line=''):
        print(line)
        out.append(line)

    p('=' * 86)
    p('GATE A -- term-level sign-flip control (exact algebra)')
    p('=' * 86)
    p('  M_DH = Arch - LogPi + LogQ - Prime ;  flipping one term must move Q by exactly')
    p('  -2 Arch(g), -2 LogQ(g), or +2 Prime(g) respectively.')
    p('')
    p('  %-5s %-7s %-7s %-24s %-24s %-11s'
      % ('L', 'sector', 'flip', 'predicted change', 'observed change', '|difference|'))
    worst = mp.mpf(0)
    for Ls in cfg.GATE_WINDOWS:
        L = dh.mpf_hi(Ls)
        for sector in (cfg.SECTOR_EVEN, cfg.SECTOR_ODD):
            idx = list(range(5))
            v = [mp.mpf(x) for x in ('0.9', '-0.55', '0.4', '-0.25', '0.15')]
            nrm = mp.sqrt(mp.fsum([c * c for c in v]))
            v = [c / nrm for c in v]
            Arch, LogPi, LogQ, Prime, _ = dh.build_parts(L, idx, cfg.DPS_BASE, sector)
            n = len(idx)

            def asm(sA, sQ, sP):
                M = mp.zeros(n, n)
                for i in range(n):
                    for j in range(n):
                        M[i, j] = (sA * Arch[i, j] - LogPi[i, j]
                                   + sQ * LogQ[i, j] - sP * Prime[i, j])
                return M

            Q0 = quad(asm(1, 1, 1), v)
            for name, (sA, sQ, sP), pred in (
                    ('arch', (-1, 1, 1), -2 * quad(Arch, v)),
                    ('logq', (1, -1, 1), -2 * quad(LogQ, v)),
                    ('prime', (1, 1, -1), +2 * quad(Prime, v))):
                obs = quad(asm(sA, sQ, sP), v) - Q0
                d = abs(pred - obs)
                worst = max(worst, d)
                p('  %-5s %-7s %-7s %-24s %-24s %-11s'
                  % (Ls, sector, name, mp.nstr(pred, 18), mp.nstr(obs, 18),
                     mp.nstr(d, 4)))
                rows.append(dict(L=Ls, sector=sector, flip=name,
                                 predicted=mp.nstr(pred, 30),
                                 observed=mp.nstr(obs, 30), residual=float(d)))
            # the conductor identity, separately: LogQ(g) must be log(5) h(0)
            h0 = h_at(v, idx, L, sector, mp.mpf(0))
            lq = quad(LogQ, v)
            d = abs(lq - mp.log(mp.mpf(5)) * h0)
            worst = max(worst, d)
            p('  %-5s %-7s %-7s LogQ(g) - log(5) h(0) = %-32s  h(0) = %s'
              % (Ls, sector, 'cond', mp.nstr(d, 4), mp.nstr(h0, 16)))
            rows.append(dict(L=Ls, sector=sector, flip='conductor_identity',
                             predicted=mp.nstr(mp.log(mp.mpf(5)) * h0, 30),
                             observed=mp.nstr(lq, 30), residual=float(d)))
    p('')
    p('  worst residual over all rows: %s' % mp.nstr(worst, 4))
    ok = worst < mp.mpf('1e-50')
    p('  SIGN-FLIP CONTROL: %s' % ('PASS' if ok else 'FAIL'))
    p('=' * 86)
    with open(os.path.join(DATA, 'gate_signflip.csv'), 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    with open(os.path.join(RES, 'gate_signflip.txt'), 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
