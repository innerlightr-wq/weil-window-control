#!/usr/bin/env python3
r"""Phase G -- sensitivity: the actual function versus a controlled deformation.

These are two different questions and are reported separately.

A.  ACTUAL-FUNCTION BEHAVIOUR.  The detected onset is compared with the heuristic
    reference scale  L_ref = (1/2) log(5 gamma_1 / 2 pi).  This is a REFERENCE SCALE,
    not a proved onset prediction.  Its provenance is a bandwidth heuristic whose
    source this extension did not verify (SOURCE_AUDIT.md E3); it is quoted because
    the brief asked for the comparison, and the comparison is reported with that caveat.

B.  CONTROLLED PERTURBATION.  Q_rest is held FIXED and the selected quartet's
    displacement is replaced by t*delta_1, multiplicities preserved:

        Q_t = Q_rest + Q_{O_1}(t delta_1, gamma_1),       t in [0, 1].

    At t = 0 the quartet has COLLAPSED to two double on-line zeros at 1/2 +- i gamma_1,
    not to a simple pair.  This is a SYNTHETIC deformation; the Davenport-Heilbronn
    function is not deformed, and nothing at t < 1 is an L-function.

    The COMPLETE quadratic coefficient, from DERIVATION.md (4.2) expanded in delta:
        P(delta) = m_0 + delta^2 m_2/2 + O(delta^4),  S(delta) = delta m_1 + O(delta^3),
        m_j = int g u^j e^{i gamma u} du,
        Q_O(delta) = 4|m_0|^2 + 4 delta^2 [ Re(m_0 conj m_2) - |m_1|^2 ] + O(delta^4).
    With F(t) = int g e^{itu} du one has m_0 = F(gamma), m_1 = -i F'(gamma),
    m_2 = -F''(gamma), so the coefficient is

        dQ/d(delta^2) = -4 [ |F'(gamma)|^2 + Re( F(gamma) conj F''(gamma) ) ],

    which for real even g is exactly the -4(F'^2 + F F'') of notes/proofs.md section B.1.
    The ISOLATED term -4|F'|^2 and the uniform bound -4K(L,gamma) are DIFFERENT
    quantities and all three are reported.

Output: results/phaseG.txt, data/sensitivity.csv.
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
from phaseD import index_set                                             # noqa: E402
from zeta_window import freq, norm2                                      # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES, DATA, CERT = (os.path.join(HERE, x) for x in ('results', 'data', 'certificates'))
GAMMA1 = '85.69934848537759217192926770894172903799'
BETA1 = '0.8085171824566373855533519606068441278507'


def K_even(L, g):
    return mp.quad(lambda u: u ** 2 * mp.sin(g * u) ** 2, [-L, 0, L])


def K_odd(L, g):
    return mp.quad(lambda u: u ** 2 * mp.cos(g * u) ** 2, [-L, 0, L])


def K_closed(L, g):
    """K_even(L,gamma) in closed form (paper eq. (K))."""
    return (L ** 3 / 3 - (L ** 2 * mp.sin(2 * g * L) / (2 * g)
                          + L * mp.cos(2 * g * L) / (2 * g ** 2)
                          - mp.sin(2 * g * L) / (4 * g ** 3)))


def moments(coeffs, idx, L, sector, gam):
    """m_j = int g u^j e^{i gamma u} du for j = 0,1,2, by quadrature on the basis."""
    def gfun(u):
        return mp.fsum([c * ((mp.sin(freq(k, L, sector) * u) if sector == cfg.SECTOR_ODD
                              else mp.cos(freq(k, L, sector) * u))
                             / mp.sqrt(norm2(k, L, sector)))
                        for c, k in zip(coeffs, idx)])
    return [mp.quad(lambda u: gfun(u) * u ** j * mp.e ** (1j * gam * u), [-L, 0, L])
            for j in (0, 1, 2)]


def main():
    out, rows = [], []

    def p(line=''):
        print(line)
        out.append(line)

    mp.mp.dps = cfg.DPS_BASE
    beta1, gam = dh.mpf_hi(BETA1), dh.mpf_hi(GAMMA1)
    d1 = beta1 - mp.mpf(1) / 2
    p('=' * 96)
    p('PHASE G  --  sensitivity')
    p('=' * 96)

    # ---------------------------------------------------------------- A
    Lref = mp.log(5 * gam / (2 * mp.pi)) / 2
    p('')
    p('(A) ACTUAL-FUNCTION BEHAVIOUR')
    p('    reference scale  L_ref = (1/2) log(5 gamma_1 / 2 pi) = %s' % mp.nstr(Lref, 16))
    p('    delta_1 = %s,  delta_1 * L at L = 2.0 : %s'
      % (mp.nstr(d1, 16), mp.nstr(d1 * 2, 10)))
    p('    This is a REFERENCE SCALE, not a predicted onset.  Its provenance (a')
    p('    bandwidth heuristic) was NOT verified; see SOURCE_AUDIT.md E3.  The measured')
    p('    onset is reported in Phase D and is compared here without any claim that the')
    p('    formula predicts it.')
    p('    Note delta_1 * L = %s at L = 2: the small-displacement approximation is NOT'
      % mp.nstr(d1 * 2, 8))
    p('    valid at the actual displacement, and no small-delta law is used in Phase D.')

    # ---------------------------------------------------------------- K identities
    p('')
    p('    K identities at gamma = gamma_1:')
    p('      %-6s %-22s %-22s %-22s %-14s' %
      ('L', 'K_even (quad)', 'K_even (closed)', 'K_odd (quad)', 'sum - 2L^3/3'))
    for Ls in ('1.0', '1.4', '1.8', '2.0', '2.6'):
        L = dh.mpf_hi(Ls)
        ke, kc, ko = K_even(L, gam), K_closed(L, gam), K_odd(L, gam)
        p('      %-6s %-22s %-22s %-22s %-14s'
          % (Ls, mp.nstr(ke, 16), mp.nstr(kc, 16), mp.nstr(ko, 16),
             mp.nstr(ke + ko - 2 * L ** 3 / 3, 4)))
        rows.append(dict(kind='K', L=Ls, sector='-', t='-',
                         a=mp.nstr(ke, 25), b=mp.nstr(ko, 25),
                         c=mp.nstr(ke + ko - 2 * L ** 3 / 3, 25), d=''))

    # ---------------------------------------------------------------- B
    p('')
    p('(B) CONTROLLED PERTURBATION  Q_t = Q_rest + Q_{O_1}(t delta_1, gamma_1)')
    with open(os.path.join(CERT, 'witness.json')) as fh:
        c = json.load(fh)
    L, sector, idx = dh.mpf_hi(c['L']), c['sector'], c['idx']
    v = [mp.mpf(s) for s in c['coefficients']]
    n2 = mp.fsum([x * x for x in v])
    n = len(idx)
    M = dh.build_matrix(L, idx, cfg.DPS_BASE, sector)
    MO1 = orbitform.orbit_matrix(idx, L, sector, beta1, gam)
    MR = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            MR[i, j] = M[i, j] - MO1[i, j]
    QR = mp.fsum([v[i] * v[j] * MR[i, j] for i in range(n) for j in range(n)]) / n2
    lamR = min(mp.eigsy(MR, eigvals_only=True))

    m0, m1, m2 = moments(v, idx, L, sector, gam)
    nv = mp.sqrt(n2)
    m0, m1, m2 = m0 / nv, m1 / nv, m2 / nv
    coef_full = 4 * (mp.re(m0 * mp.conj(m2)) - abs(m1) ** 2)
    coef_iso = -4 * abs(m1) ** 2
    Kv = K_even(L, gam) if sector != cfg.SECTOR_ODD else K_odd(L, gam)
    p('    frozen witness, L = %s, sector %s' % (c['L'], sector))
    p('    Q_rest(g_*)/||g||^2 = %-22s   lambda_min(Q_rest) = %s'
      % (mp.nstr(QR, 18), mp.nstr(lamR, 18)))
    p('    |m_0| = |F(gamma_1)|   = %s' % mp.nstr(abs(m0), 12))
    p("    |m_1| = |F'(gamma_1)|  = %s" % mp.nstr(abs(m1), 12))
    p("    |m_2| = |F''(gamma_1)| = %s" % mp.nstr(abs(m2), 12))
    p('    COMPLETE quadratic coefficient  4[Re(m0 conj m2) - |m1|^2] = %s'
      % mp.nstr(coef_full, 16))
    p("    ISOLATED derivative term        -4|F'(gamma_1)|^2          = %s"
      % mp.nstr(coef_iso, 16))
    p('    UNIFORM bound                   -4 K(L, gamma_1)           = %s'
      % mp.nstr(-4 * Kv, 16))
    p('    (these are three different numbers; the first is the actual expansion')
    p('     coefficient, the second drops the F F\'\' term, the third is a sup over the')
    p('     whole admissible sector and bounds only the second.)')
    p('')
    p('    %-6s %-24s %-24s %-24s %-14s' %
      ('t', 'Q_t(g_*)/||g||^2', 'quadratic model', 'residual', 'lambda_min(Q_t)'))
    q0 = None
    for ts in ('1.0', '0.8', '0.6', '0.4', '0.2', '0.1', '0.05', '0.0'):
        t = dh.mpf_hi(ts)
        MOt = orbitform.orbit_matrix(idx, L, sector,
                                     mp.mpf(1) / 2 + t * d1, gam)
        Mt = mp.matrix(n, n)
        for i in range(n):
            for j in range(n):
                Mt[i, j] = MR[i, j] + MOt[i, j]
        qt = mp.fsum([v[i] * v[j] * Mt[i, j]
                      for i in range(n) for j in range(n)]) / n2
        lamt = min(mp.eigsy(Mt, eigvals_only=True))
        if ts == '0.0':
            q0 = qt
        model = (q0 + coef_full * (t * d1) ** 2) if q0 is not None else mp.nan
        rows.append(dict(kind='perturbation', L=c['L'], sector=sector, t=ts,
                         a=mp.nstr(qt, 25), b=mp.nstr(lamt, 25), c='', d=''))
        p('    %-6s %-24s %-24s %-24s %-14s'
          % (ts, mp.nstr(qt, 18), '' if q0 is None else mp.nstr(model, 18),
             '' if q0 is None else mp.nstr(qt - model, 8), mp.nstr(lamt, 12)))
    # redo with q0 known, to print the model column properly
    p('')
    p('    (the t = 0 row is the collapsed configuration and fixes the model constant;')
    p('     rows are printed in decreasing t, so the model column is blank until t = 0.')
    p('     The reconciled table with the model is written to data/sensitivity.csv.)')
    p('')
    p('    t * delta_1 * L at t = 1 : %s  -- far outside any small-displacement regime,'
      % mp.nstr(d1 * L, 8))
    p('    so the quadratic model is reported WITH its residual and is not used as an')
    p('    approximation anywhere else in this extension.')

    # recompute the model residuals now that q0 is known, for the CSV
    for r in rows:
        if r['kind'] == 'perturbation':
            t = mp.mpf(r['t'])
            model = q0 + coef_full * (t * d1) ** 2
            r['c'] = mp.nstr(model, 25)
            r['d'] = mp.nstr(mp.mpf(r['a']) - model, 25)

    # ---------------------------------------------------------------- overlaps
    p('')
    p('(C) NORMALISED OVERLAPS with the isolated sensitivity optimiser')
    tgt = (lambda u: u * mp.sin(gam * u)) if sector != cfg.SECTOR_ODD \
        else (lambda u: u * mp.cos(gam * u))
    tname = 'u sin(gamma u)' if sector != cfg.SECTOR_ODD else 'u cos(gamma u)'
    tn = mp.sqrt(mp.quad(lambda u: tgt(u) ** 2, [-L, 0, L]))

    def ov(cf):
        nn = mp.sqrt(mp.fsum([x * x for x in cf]))
        s = mp.fsum([cf[m] / nn * mp.quad(
            lambda u: ((mp.sin(freq(k, L, sector) * u) if sector == cfg.SECTOR_ODD
                        else mp.cos(freq(k, L, sector) * u))
                       / mp.sqrt(norm2(k, L, sector))) * tgt(u), [-L, 0, L])
            for m, k in enumerate(idx)])
        return s / tn
    EF, VF = mp.eigsy(M)
    lo = min(range(n), key=lambda i: EF[i])
    vF = [VF[i, lo] for i in range(n)]
    p('    target (unit-normalised) : %s' % tname)
    p('    overlap with the frozen witness g_*            : %s' % mp.nstr(abs(ov(v)), 10))
    p('    overlap with the Q_F ground state              : %s' % mp.nstr(abs(ov(vF)), 10))
    ER, VR = mp.eigsy(MR)
    loR = min(range(n), key=lambda i: ER[i])
    p('    overlap with the Q_rest ground state           : %s'
      % mp.nstr(abs(ov([VR[i, loR] for i in range(n)])), 10))
    p('    A small overlap does not refute the sensitivity identity: u sin/cos(gamma u)')
    p('    optimises the ISOLATED functional f -> F\'(gamma), not the full form.')
    rows.append(dict(kind='overlap', L=c['L'], sector=sector, t=tname,
                     a=mp.nstr(abs(ov(v)), 25), b=mp.nstr(abs(ov(vF)), 25),
                     c=mp.nstr(abs(ov([VR[i, loR] for i in range(n)])), 25), d=''))

    with open(os.path.join(DATA, 'sensitivity.csv'), 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=['kind', 'L', 'sector', 't', 'a', 'b', 'c', 'd'])
        w.writeheader()
        for r in rows:
            w.writerow(r)
    p('')
    p('  -> data/sensitivity.csv   (columns a,b,c,d are kind-dependent; see header rows)')
    p('=' * 96)
    with open(os.path.join(RES, 'phaseG.txt'), 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
