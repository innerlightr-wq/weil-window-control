#!/usr/bin/env python3
r"""Phase A -- define and validate the actual Davenport--Heilbronn function.

Checks, in order:

  A1  kappa derived from the Gauss sum equals the closed form, and the root-number
      identity (1 - i kappa) eps_chi = 1 + i kappa holds exactly.
  A2  the three evaluation routes agree (A and B share mp.zeta; C does not; D is a
      Dirichlet-series sanity check at Re s = 4).
  A3  the completed functional equation Lambda_F(s) = Lambda_F(1-s) at 20 reproducibly
      chosen points, with a SCALE-AWARE residual, at two precisions.
  A4  entireness of the completion: F vanishes at the gamma poles s = -1, -3, -5, -7.
  A5  the zero-containing strip: |F(s) - 1| <= zeta(2) - 1 < 1 on Re s = 2, sampled, with
      the proof recorded in DERIVATION.md; reflected boundary Re s = -1.
  A6  refinement of the literature-reported off-line zero near 0.8085171824566374 +
      85.6993484853776i and its full symmetry orbit.

Output: results/phaseA.txt, results/phaseA.json, data/dh_offline_orbit.csv.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp                                                      # noqa: E402
import dh                                                                # noqa: E402
import config as cfg                                                     # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(HERE, 'results')
DATA = os.path.join(HERE, 'data')

SEED_BETA = '0.8085171824566374'
SEED_GAMMA = '85.6993484853776'

# 20 reproducible FE test points: fixed decimal strings, chosen to avoid the critical
# line (where the FE is trivial), the real axis, and the gamma poles.
FE_POINTS = [
    ('2.0', '1.0'), ('2.0', '7.5'), ('1.7', '0.3'), ('1.3', '12.25'),
    ('1.1', '30.5'), ('0.9', '3.75'), ('0.75', '21.0'), ('0.6', '5.5'),
    ('0.25', '9.125'), ('0.1', '2.5'), ('-0.3', '4.0'), ('-0.75', '17.5'),
    ('3.5', '0.5'), ('2.6', '45.25'), ('1.45', '85.699'), ('0.81', '85.7'),
    ('-0.5', '0.875'), ('0.35', '60.0'), ('1.9', '100.5'), ('0.45', '150.25'),
]


def fe_residual(sigma, t, dps):
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        s = mp.mpc(mp.mpf(sigma), mp.mpf(t))
        a = dh.Lambda_F(s)
        b = dh.Lambda_F(1 - s)
        scale = max(abs(a), abs(b))
        return abs(a - b), (abs(a - b) / scale if scale > 0 else mp.inf), scale
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

    p('=' * 78)
    p('PHASE A  --  the Davenport-Heilbronn function')
    p('=' * 78)

    # ---------------------------------------------------------------- A1
    p('')
    p('[A1] kappa from the Gauss sum')
    tau = dh.gauss_sum()
    eps = dh.root_number()
    k = dh.kappa()
    kc = dh.kappa_closed_form()
    p('  tau(chi)            = %s' % mp.nstr(tau, 30))
    p('  |tau| - sqrt(5)     = %s' % mp.nstr(abs(abs(tau) - mp.sqrt(5)), 5))
    p('  eps_chi             = %s' % mp.nstr(eps, 30))
    p('  |eps_chi| - 1       = %s' % mp.nstr(abs(abs(eps) - 1), 5))
    p('  kappa  (derived)    = %s' % mp.nstr(k, 40))
    p('  kappa  (closed form)= %s' % mp.nstr(kc, 40))
    p('  |difference|        = %s' % mp.nstr(abs(k - kc), 5))
    rni = abs((1 - 1j * k) * eps - (1 + 1j * k))
    p('  |(1-i k) eps - (1+i k)| = %s' % mp.nstr(rni, 5))
    rec['kappa'] = mp.nstr(k, 40)
    rec['kappa_closed_form_diff'] = float(abs(k - kc))
    rec['root_number_identity_residual'] = float(rni)
    rec['a_vector'] = [mp.nstr(x, 30) for x in dh.a_vector()]
    p('  (a_1..a_5)          = (%s)' % ', '.join(mp.nstr(x, 18) for x in dh.a_vector()))

    # ---------------------------------------------------------------- A2
    p('')
    p('[A2] evaluation routes   (A,B share mp.zeta; C is Euler-Maclaurin, independent)')
    pts = [('2.0', '0.0'), ('2.0', '1.0'), ('1.7', '5.3'), ('0.5', '14.1'),
           ('-0.3', '3.0'), ('0.8085171824566374', '85.6993484853776')]
    rows = []
    for sg, tt in pts:
        s = mp.mpc(mp.mpf(sg), mp.mpf(tt))
        A, B, C = dh.F_route_A(s), dh.F_route_B(s), dh.F_route_C(s)
        sc = max(abs(A), mp.mpf(10) ** (-cfg.DPS_BASE))
        p('  s = %-24s |A-B|/|A| = %-10s  |A-C|/|A| = %-10s' %
          ('%s%+sj' % (sg, tt), mp.nstr(abs(A - B) / sc, 4), mp.nstr(abs(A - C) / sc, 4)))
        rows.append(dict(sigma=sg, t=tt, rel_AB=float(abs(A - B) / sc),
                         rel_AC=float(abs(A - C) / sc)))
    sD = mp.mpf(4)
    D, B4 = dh.F_route_D(sD), dh.F_route_B(sD)
    p('  route D (direct series, Re s = 4, n <= 20000): |D-B| = %s   (tail bound %s)' %
      (mp.nstr(abs(D - B4), 4), mp.nstr(mp.mpf(20000) ** (-3) / 3, 4)))
    rec['routes'] = rows
    rec['route_D_diff'] = float(abs(D - B4))

    # ---------------------------------------------------------------- A3
    p('')
    p('[A3] completed functional equation  Lambda_F(s) = Lambda_F(1-s)')
    p('     %-20s %-13s %-13s %-13s' % ('s', '|Lam|', 'abs resid', 'rel resid'))
    fe = []
    worst_rel = mp.mpf(0)
    for sg, tt in FE_POINTS:
        ab, rel, sc = fe_residual(sg, tt, cfg.DPS_BASE)
        worst_rel = max(worst_rel, rel)
        p('     %-20s %-13s %-13s %-13s' %
          ('%s%+sj' % (sg, tt), mp.nstr(sc, 6), mp.nstr(ab, 4), mp.nstr(rel, 4)))
        fe.append(dict(sigma=sg, t=tt, scale=float(sc), abs=float(ab), rel=float(rel)))
    p('     worst relative residual at %d dps : %s' % (cfg.DPS_BASE, mp.nstr(worst_rel, 4)))
    worst_hi = mp.mpf(0)
    for sg, tt in FE_POINTS:
        _, rel, _ = fe_residual(sg, tt, cfg.DPS_HIGH)
        worst_hi = max(worst_hi, rel)
    p('     worst relative residual at %d dps : %s   (scales with precision)' %
      (cfg.DPS_HIGH, mp.nstr(worst_hi, 4)))
    rec['fe'] = fe
    rec['fe_worst_rel_base'] = float(worst_rel)
    rec['fe_worst_rel_high'] = float(worst_hi)

    # ---------------------------------------------------------------- A4
    p('')
    p('[A4] entireness: Gamma((s+1)/2) has poles at s = -1,-3,-5,... ; F must vanish')
    mp.mp.dps = cfg.DPS_BASE
    triv = []
    for m in (1, 3, 5, 7, 9):
        v = dh.F(mp.mpf(-m))
        p('     F(%-3d) = %-14s' % (-m, mp.nstr(v, 6)))
        triv.append(dict(s=-m, F=float(mp.re(v))))
    p('     (these are the trivial zeros; Lambda_F is entire and has none of them)')
    for m in (1, 3, 5):
        lv = dh.Lambda_F(mp.mpf(-m) + mp.mpf(10) ** (-20))
        p('     Lambda_F(%d + 1e-20) = %s   (finite: no pole)' % (-m, mp.nstr(lv, 8)))
    rec['trivial_zeros'] = triv

    # ---------------------------------------------------------------- A5
    p('')
    p('[A5] zero-containing strip')
    z2m1 = mp.zeta(2) - 1
    p('     bound |F(s)-1| <= sum_{n>=2} |a_n| n^-2 <= zeta(2)-1 = %s < 1  for Re s >= 2'
      % mp.nstr(z2m1, 12))
    mx = mp.mpf(0)
    for j in range(61):
        t = mp.mpf(j) * 5
        mx = max(mx, abs(dh.F(mp.mpc(2, t)) - 1))
    p('     sampled max |F(2+it)-1| over t = 0,5,...,300 : %s' % mp.nstr(mx, 8))
    p('     => no zeros with Re s >= 2; by Lambda_F(s)=Lambda_F(1-s), none with Re s <= -1.')
    p('     All zeros of Lambda_F lie in the open strip  -1 < Re s < 2.')
    rec['strip'] = dict(bound=float(z2m1), sampled_max=float(mx))

    # ---------------------------------------------------------------- A6
    p('')
    p('[A6] the selected off-line zero and its orbit')
    mp.mp.dps = cfg.DPS_HIGH
    seed = mp.mpc(mp.mpf(SEED_BETA), mp.mpf(SEED_GAMMA))
    p('     literature SEARCH SEED (not a certified root): %s' % mp.nstr(seed, 20))
    rho = mp.findroot(dh.F, seed, tol=mp.mpf(10) ** (-2 * cfg.DPS_HIGH + 20))
    beta, gamma = mp.re(rho), mp.im(rho)
    Fv = dh.F(rho)
    dF = mp.diff(dh.F, rho)
    p('     refined rho   = %s' % mp.nstr(rho, 40))
    p('     beta          = %s' % mp.nstr(beta, 40))
    p('     gamma         = %s' % mp.nstr(gamma, 40))
    p('     delta = b-1/2 = %s' % mp.nstr(beta - mp.mpf(1) / 2, 40))
    p('     |F(rho)|      = %s' % mp.nstr(abs(Fv), 6))
    p("     |F(rho)/F'(rho)| (Newton step) = %s" % mp.nstr(abs(Fv / dF), 6))
    p('     |seed - rho|  = %s' % mp.nstr(abs(seed - rho), 6))
    orbit = [rho, mp.conj(rho), 1 - rho, 1 - mp.conj(rho)]
    labels = ['beta+i*gamma', 'beta-i*gamma', '1-beta-i*gamma', '1-beta+i*gamma']
    rows = []
    for lab, z in zip(labels, orbit):
        fz = dh.F(z)
        lz = dh.Lambda_F(z)
        p('       %-16s %-46s |F| = %-10s |Lambda_F| = %s' %
          (lab, mp.nstr(z, 28), mp.nstr(abs(fz), 4), mp.nstr(abs(lz), 4)))
        rows.append(dict(label=lab, re=mp.nstr(mp.re(z), 40), im=mp.nstr(mp.im(z), 40),
                         absF=float(abs(fz))))
    # stability under precision
    mp.mp.dps = cfg.DPS_HIGH + 40
    rho2 = mp.findroot(dh.F, rho, tol=mp.mpf(10) ** (-2 * (cfg.DPS_HIGH + 20)))
    p('     reruns at %d dps: |rho(%d) - rho(%d)| = %s' %
      (cfg.DPS_HIGH + 40, cfg.DPS_HIGH, cfg.DPS_HIGH + 40, mp.nstr(abs(rho - rho2), 6)))
    rec['offline'] = dict(beta=mp.nstr(mp.re(rho2), 50), gamma=mp.nstr(mp.im(rho2), 50),
                          delta=mp.nstr(mp.re(rho2) - mp.mpf(1) / 2, 50),
                          absF=float(abs(dh.F(rho2))),
                          newton_step=float(abs(dh.F(rho2) / mp.diff(dh.F, rho2))),
                          seed_distance=float(abs(seed - rho2)),
                          precision_stability=float(abs(rho - rho2)),
                          orbit=rows)
    mp.mp.dps = cfg.DPS_HIGH
    with open(os.path.join(DATA, 'dh_offline_orbit.csv'), 'w') as fh:
        fh.write('label,re,im,absF\n')
        for r in rows:
            fh.write('%s,%s,%s,%.3e\n' % (r['label'], r['re'], r['im'], r['absF']))

    p('')
    p('=' * 78)
    ok = (abs(k - kc) < mp.mpf(10) ** (-cfg.DPS_BASE + 10) and rni < mp.mpf(10) ** (-40)
          and worst_rel < mp.mpf(10) ** (-cfg.DPS_BASE + 15))
    p('PHASE A GATE: %s' % ('PASS' if ok else 'FAIL'))
    rec['gate'] = 'PASS' if ok else 'FAIL'
    p('=' * 78)

    with open(os.path.join(RES, 'phaseA.txt'), 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    with open(os.path.join(RES, 'phaseA.json'), 'w') as fh:
        json.dump(rec, fh, indent=1)
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
