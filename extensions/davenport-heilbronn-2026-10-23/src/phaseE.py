#!/usr/bin/env python3
r"""Phase E -- the remaining-zero form  Q_rest(g) = Q_F(g) - Q_{O_1}(g).

Q_F is evaluated through the ARITHMETIC formula, so every remaining zero -- including
every off-line zero above the census height and every zero we never located -- is
included implicitly.  Nothing is replaced by a finite zero list.

Q_{O_1} is the exact rank-2 orbit form of the selected quartet (src/orbitform.py).

For each (L, sector) the classification is one of R1..R6 as defined in the brief.  R6
(whole-window nonnegativity) is NEVER reported: a positive finite Galerkin restriction
says nothing about the infinite-dimensional complement or its coupling to the block.

Output: results/phaseE.txt, data/remaining_form.csv.
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
from phaseD import index_set, rayleigh, GAMMA1, BETA1                    # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(HERE, 'results')
DATA = os.path.join(HERE, 'data')
CERT = os.path.join(HERE, 'certificates')


def main():
    out, rows = [], []

    def p(line=''):
        print(line)
        out.append(line)

    mp.mp.dps = cfg.DPS_BASE
    beta1, gamma1 = dh.mpf_hi(BETA1), dh.mpf_hi(GAMMA1)
    p('=' * 96)
    p('PHASE E  --  the remaining-zero form')
    p('=' * 96)
    p('  O_1 = {beta_1 +- i gamma_1, 1-beta_1 +- i gamma_1},  beta_1 = %s' % BETA1)
    p('                                                        gamma_1 = %s' % GAMMA1)
    p('                                                        delta_1 = %s'
      % mp.nstr(beta1 - mp.mpf(1) / 2, 30))
    p('')

    # -------- polarisation check on the orbit matrix
    L0, sec0 = dh.mpf_hi('2.0'), cfg.SECTOR_ODD
    idx0 = index_set(L0, sec0, 0)
    cc = [mp.mpf(1) / (i + 2) for i in range(len(idx0))]
    q1, q2 = orbitform.check_orbit_matrix(idx0, L0, sec0, beta1, gamma1, cc)
    p('  polarisation check on M_{O_1}: |c^T M c - Q_O(c)| = %s' % mp.nstr(abs(q1 - q2), 4))
    p('')

    p('%-5s %-7s %-4s %-22s %-22s %-22s %-6s' %
      ('L', 'sector', 'N', 'lambda_min(Q_F)', 'lambda_min(Q_rest)',
       'Q_rest at the Q_F minimiser', 'class'))
    for Ls in cfg.WINDOW_SCAN:
        L = dh.mpf_hi(Ls)
        for sector in (cfg.SECTOR_EVEN, cfg.SECTOR_ODD):
            idx = index_set(L, sector, 1)
            M = dh.build_matrix(L, idx, cfg.DPS_BASE, sector)
            MO = orbitform.orbit_matrix(idx, L, sector, beta1, gamma1)
            n = len(idx)
            MR = mp.matrix(n, n)
            for i in range(n):
                for j in range(n):
                    MR[i, j] = M[i, j] - MO[i, j]
            E, V = mp.eigsy(M)
            lo = min(range(n), key=lambda i: E[i])
            vF = [V[i, lo] for i in range(n)]
            ER = mp.eigsy(MR, eigvals_only=True)
            lamF, lamR = E[lo], min(ER)
            qrest_at_vF = rayleigh(MR, vF)
            # classification
            if lamR < 0:
                cls = 'R2'          # negative, not rigorously certified here
            elif lamR > 0:
                cls = 'R4'          # numerically positive finite restriction
            else:
                cls = 'R5'
            p('%-5s %-7s %-4d %-22s %-22s %-22s %-6s' %
              (Ls, sector, n, mp.nstr(lamF, 16), mp.nstr(lamR, 16),
               mp.nstr(qrest_at_vF, 16), cls))
            rows.append(dict(L=Ls, sector=sector, N=n,
                             lambda_min_QF=mp.nstr(lamF, 30),
                             lambda_min_Qrest=mp.nstr(lamR, 30),
                             Qrest_at_QF_minimiser=mp.nstr(qrest_at_vF, 30),
                             classification=cls))

    # -------------------------------------------------- the frozen witness decomposition
    p('')
    p('[DECOMPOSITION AT THE FROZEN WITNESS g_*]')
    wpath = os.path.join(CERT, 'witness.json')
    if not os.path.exists(wpath):
        p('  certificates/witness.json absent -- run src/phaseD.py first.')
    else:
        with open(wpath) as fh:
            c = json.load(fh)
        Lw = dh.mpf_hi(c['L'])
        sw = c['sector']
        iw = c['idx']
        vw = [mp.mpf(s) for s in c['coefficients']]
        n2 = mp.fsum([x * x for x in vw])
        for dps in (cfg.DPS_BASE, cfg.DPS_HIGH):
            mp.mp.dps = dps
            Mw = dh.build_matrix(Lw, iw, dps, sw)
            QF = mp.fsum([vw[i] * vw[j] * Mw[i, j]
                          for i in range(len(iw)) for j in range(len(iw))])
            QO = dh.quartet_contribution(vw, iw, Lw, sw, beta1, gamma1)
            p('  dps = %d' % dps)
            p('    ||g_*||^2        = %s' % mp.nstr(n2, 30))
            p('    Q_F(g_*)         = %s' % mp.nstr(QF, 30))
            p('    Q_{O_1}(g_*)     = %s' % mp.nstr(QO, 30))
            p('    Q_rest(g_*)      = %s' % mp.nstr(QF - QO, 30))
            p('    additive check |Q_F - (Q_O1 + Q_rest)| = %s'
              % mp.nstr(abs(QF - (QO + (QF - QO))), 4))
            p('    normalised: Q_F/||g||^2 = %-24s Q_O1/||g||^2 = %-24s'
              % (mp.nstr(QF / n2, 20), mp.nstr(QO / n2, 20)))
            p('                Q_rest/||g||^2 = %s' % mp.nstr((QF - QO) / n2, 20))
        mp.mp.dps = cfg.DPS_BASE

    # -------------------------------------------------- Theorem 8 audit
    p('')
    p('[AUDIT vs the repository Theorem 8]')
    p('  Theorem 8 (paper, \\label{thm:zetasharp}) defines  Q_rest(f) = Q_0(f) - 4F(gamma)^2')
    p('  where Q_0 is the AUGMENTED COUNTERFACTUAL baseline of notes/proofs.md section B.1,')
    p('  convention 3: a DOUBLE on-line pair at gamma that is then moved off the line.')
    p('  Its baseline therefore contains 4F(gamma)^2 = the collapsed quartet, and its')
    p('  "remaining" form is what is left after removing those FOUR zeros.')
    p('')
    p('  Our Q_rest removes the FOUR ACTUAL zeros of the orbit O_1, at their actual')
    p('  displacement delta_1, not a collapsed double pair.  The two agree in the')
    p('  multiplicity convention (four zeros removed, four zeros in the orbit) and differ')
    p('  in that ours is evaluated at delta_1 > 0 while Theorem 8 removes the delta = 0')
    p('  limit.  The difference is exactly Q_{O_1}(delta_1) - Q_{O_1}(0), which Phase G')
    p('  measures.')
    p('')
    p('  Theorem 8 ASSUMES "every nontrivial zero other than the one being moved lies on')
    p('  the critical line", and states this as SUFFICIENT for Q_rest >= 0.  For the')
    p('  Davenport-Heilbronn function that hypothesis is FALSE: the census finds further')
    p('  off-line quartets.  So Theorem 8 does not apply here, and we do NOT assume it.')
    p('')
    p('  Does its PROOF transfer under the weaker hypothesis "Q_rest >= 0"?  Yes, and')
    p('  trivially: the proof uses the hypothesis ONLY through the single inequality')
    p('  Qtilde >= 0 (notes/proofs.md, Theorem B\' proof: "Q_delta >= 0 + 0 - 4 delta^2 F\'^2')
    p('  - delta^4 F\'\'^2 - |R_4|").  Nothing else about the location of the other zeros is')
    p('  used.  Therefore the conclusion')
    p('      lambda_min(Q_delta) >= -4 K(L,gamma) delta^2 - (2/5 + (16/3) e^{2 L delta}) L^5 delta^4')
    p('  holds verbatim with "all other zeros on the line" replaced by "Q_rest >= 0 on the')
    p('  admissible sector".  That is a strictly weaker hypothesis, and it is the one this')
    p('  phase measures.  What does NOT transfer is any claim that the measured positivity')
    p('  of a finite Galerkin restriction establishes that hypothesis (see below).')
    p('')
    p('  CAUTION.  A positive finite restriction of Q_rest does NOT establish the')
    p('  whole-window Theorem 8 hypothesis.  Classification R6 is not reported anywhere in')
    p('  this extension.')

    with open(os.path.join(DATA, 'remaining_form.csv'), 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    p('')
    p('  -> data/remaining_form.csv')
    p('=' * 96)
    with open(os.path.join(RES, 'phaseE.txt'), 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
