#!/usr/bin/env python3
r"""Phase F -- masking, removal, and the finite-spectrum synthetic hybrid.

Two outputs are kept strictly apart:

  (A) FIXED-WITNESS contributions.  These are additive: removing an orbit changes
      Q_F(g_*) by exactly -Q_O(g_*), and the changes for different orbits add.
  (B) REOPTIMISED eigenvalues.  lambda_min of (M - M_O) is NOT lambda_min(M) - Q_O
      at anything, and differences of reoptimised eigenvalues are NOT an attribution.

SYNTHETIC HYBRID.  Keeping O_1 fixed, every OTHER enumerated off-line orbit is moved to
Re s = 1/2 at its original ordinate with multiplicity preserved: the quartet
{beta +- i gamma, 1-beta +- i gamma} becomes a pair of DOUBLE on-line zeros at
1/2 +- i gamma.  This is called a finite-spectrum hybrid and nothing else: the
unenumerated off-line zeros above the census height are untouched and unaccounted for,
so this is NOT "all other zeros moved to the line", and no analytic function is claimed
to have the modified spectrum.

Output: results/phaseF.txt, data/masking.csv.
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp                                                      # noqa: E402
import dh                                                                # noqa: E402
import orbitform                                                         # noqa: E402
import config as cfg                                                     # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES, DATA, CERT = (os.path.join(HERE, x) for x in ('results', 'data', 'certificates'))


def load_offline():
    path = os.path.join(DATA, 'dh_zeros.csv')
    off = []
    if not os.path.exists(path):
        return off
    with open(path) as fh:
        for row in csv.DictReader(l for l in fh if not l.startswith('#')):
            if row['kind'] == 'offline':
                off.append((dh.mpf_hi(row['beta']), dh.mpf_hi(row['gamma']),
                            int(row['multiplicity'])))
    off.sort(key=lambda r: r[1])
    return off


def main():
    out, rows = [], []

    def p(line=''):
        print(line)
        out.append(line)

    mp.mp.dps = cfg.DPS_BASE
    off = load_offline()
    p('=' * 96)
    p('PHASE F  --  masking and synthetic controls')
    p('=' * 96)
    if not off:
        p('  no census: run src/census.py first.')
        with open(os.path.join(RES, 'phaseF.txt'), 'w') as fh:
            fh.write('\n'.join(out) + '\n')
        return 1
    p('  enumerated off-line orbits (Re > 1/2) from the census:')
    for i, (b, g, m) in enumerate(off):
        p('    O_%d : beta = %-22s gamma = %-22s delta = %-20s mult %d'
          % (i + 1, mp.nstr(b, 18), mp.nstr(g, 18),
             mp.nstr(b - mp.mpf(1) / 2, 16), m))

    with open(os.path.join(CERT, 'witness.json')) as fh:
        c = json.load(fh)
    L, sector, idx = dh.mpf_hi(c['L']), c['sector'], c['idx']
    v = [mp.mpf(s) for s in c['coefficients']]
    n2 = mp.fsum([x * x for x in v])
    n = len(idx)
    M = dh.build_matrix(L, idx, cfg.DPS_BASE, sector)
    QF = mp.fsum([v[i] * v[j] * M[i, j] for i in range(n) for j in range(n)])
    EF = mp.eigsy(M, eigvals_only=True)
    p('')
    p('  frozen witness: L = %s, sector %s, N = %d' % (c['L'], sector, n))
    p('  Q_F(g_*)/||g||^2 = %s      lambda_min(Q_F) = %s'
      % (mp.nstr(QF / n2, 20), mp.nstr(min(EF), 20)))

    p('')
    p('(A) FIXED-WITNESS removal -- additive')
    p('    %-6s %-24s %-24s %-24s' %
      ('orbit', 'Q_O(g_*)/||g||^2', 'Q_F - Q_O, normalised', 'still negative?'))
    Os = []
    for i, (b, g, m) in enumerate(off):
        qo = dh.quartet_contribution(v, idx, L, sector, b, g, m)
        Os.append(qo)
        p('    O_%-4d %-24s %-24s %-24s'
          % (i + 1, mp.nstr(qo / n2, 18), mp.nstr((QF - qo) / n2, 18),
             'yes' if (QF - qo) < 0 else 'NO'))
        rows.append(dict(kind='fixed_witness_removal', orbit='O_%d' % (i + 1),
                         beta=mp.nstr(b, 30), gamma=mp.nstr(g, 30),
                         value=mp.nstr(qo / n2, 30),
                         after=mp.nstr((QF - qo) / n2, 30)))
    tot = mp.fsum(Os)
    p('    sum of all four                     = %s' % mp.nstr(tot / n2, 18))
    p('    Q_F - (all four), normalised        = %s' % mp.nstr((QF - tot) / n2, 18))
    p('    additivity check |sum - (Q_F - (Q_F - sum))| = %s'
      % mp.nstr(abs(tot - (QF - (QF - tot))), 4))

    p('')
    p('(B) REOPTIMISED -- NOT an additive attribution')
    p('    %-6s %-24s %-24s' % ('orbit', 'lambda_min(M - M_O)', 'change vs lambda_min(M)'))
    for i, (b, g, m) in enumerate(off):
        MO = orbitform.orbit_matrix(idx, L, sector, b, g, m)
        MR = mp.matrix(n, n)
        for a in range(n):
            for bb in range(n):
                MR[a, bb] = M[a, bb] - MO[a, bb]
        lam = min(mp.eigsy(MR, eigvals_only=True))
        p('    O_%-4d %-24s %-24s' % (i + 1, mp.nstr(lam, 18),
                                      mp.nstr(lam - min(EF), 18)))
        rows.append(dict(kind='reoptimised_removal', orbit='O_%d' % (i + 1),
                         beta=mp.nstr(b, 30), gamma=mp.nstr(g, 30),
                         value=mp.nstr(lam, 30),
                         after=mp.nstr(lam - min(EF), 30)))
    # remove O_1 only, then ask whether reoptimisation finds another witness
    MO1 = orbitform.orbit_matrix(idx, L, sector, off[0][0], off[0][1], off[0][2])
    MR1 = mp.matrix(n, n)
    for a in range(n):
        for bb in range(n):
            MR1[a, bb] = M[a, bb] - MO1[a, bb]
    lam1 = min(mp.eigsy(MR1, eigvals_only=True))
    p('')
    p('    removing O_1 alone: lambda_min = %s  -> reoptimisation %s another negative'
      ' direction' % (mp.nstr(lam1, 18), 'STILL FINDS' if lam1 < 0 else 'finds NO'))
    p('    (failure to find another witness would NOT restore whole-window positivity.)')

    # ------------------------------------------------------------- synthetic hybrid
    p('')
    p('(C) FINITE-SPECTRUM HYBRID: O_1 kept, O_2..O_%d moved to Re s = 1/2 at their'
      ' ordinates, multiplicity preserved (quartet -> two DOUBLE on-line zeros)' % len(off))
    MH = mp.matrix(n, n)
    for a in range(n):
        for bb in range(n):
            MH[a, bb] = M[a, bb]
    for b, g, m in off[1:]:
        MO = orbitform.orbit_matrix(idx, L, sector, b, g, m)
        MC = orbitform.orbit_matrix(idx, L, sector, mp.mpf(1) / 2, g, m)   # collapsed
        for a in range(n):
            for bb in range(n):
                MH[a, bb] += MC[a, bb] - MO[a, bb]
    QH = mp.fsum([v[i] * v[j] * MH[i, j] for i in range(n) for j in range(n)])
    lamH = min(mp.eigsy(MH, eigvals_only=True))
    p('    fixed witness : Q_hybrid(g_*)/||g||^2 = %-22s  (was %s)'
      % (mp.nstr(QH / n2, 18), mp.nstr(QF / n2, 18)))
    p('    change        : %-22s   sign: %s'
      % (mp.nstr((QH - QF) / n2, 18),
         'the other orbits were MASKING (removing them deepens the negative)'
         if QH < QF else
         'the other orbits were REINFORCING (removing them weakens the negative)'))
    p('    reoptimised   : lambda_min(hybrid) = %-22s  (was %s)'
      % (mp.nstr(lamH, 18), mp.nstr(min(EF), 18)))
    rows.append(dict(kind='synthetic_hybrid', orbit='O_2..O_%d -> line' % len(off),
                     beta='0.5', gamma='various',
                     value=mp.nstr(QH / n2, 30), after=mp.nstr(lamH, 30)))
    p('')
    p('    This is a FINITE-SPECTRUM HYBRID.  The unenumerated off-line zeros above the')
    p('    census height are untouched; no claim is made that "all other zeros were')
    p('    moved", and no Euler product or analytic function is attached to it.')

    with open(os.path.join(DATA, 'masking.csv'), 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    p('')
    p('  -> data/masking.csv')
    p('=' * 96)
    with open(os.path.join(RES, 'phaseF.txt'), 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
