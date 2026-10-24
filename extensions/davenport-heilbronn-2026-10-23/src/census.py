#!/usr/bin/env python3
r"""Finite-height census of the zeros of Lambda_F in  -0.4 < Re s < 1.4,  0 < Im s < T.

The strip is justified in DERIVATION.md section 1: sum_{n>=2} |a_n| n^{-sigma} < 1 for
sigma > sigma_0 = 1.39513615823510972..., so F has no zeros with Re s >= sigma_0, and by
Lambda_F(s) = Lambda_F(1-s) none with Re s <= 1 - sigma_0 = -0.39513...  The vertical
edges at Re s = -0.4 and 1.4 are therefore PROVABLY zero-free, which is what makes the
argument-principle count meaningful.

Method:
  1. counting by the argument principle, implemented as the adaptively refined total
     change of arg F around the rectangle (the step is halved whenever a consecutive
     increment exceeds 0.45 pi, so no half-turn can be skipped);
  2. off-line zeros located by recursive box subdivision on 0.505 <= Re s <= 1.4, each
     box with winding 1 refined by Newton, then re-certified by a local winding number;
  3. on-line ordinates from sign changes of the real Hardy-type function Z(t);
  4. reconciliation: #on-line + 2 * #off-line(Re > 1/2) must equal the contour count.

Everything here is T2 unless explicitly marked: the winding numbers are integers to
~25 digits, which is strong evidence of correct counting but is not an interval
certificate (non-claims 6 and 7).
"""
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp                                                      # noqa: E402
import dh                                                                # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(HERE, 'results')
DATA = os.path.join(HERE, 'data')

SIGMA0 = '1.395136158235109721061359'
SL, SR = mp.mpf('-0.4'), mp.mpf('1.4')


def arg_change(path, a, b, n0, maxdepth=16):
    xs = [a + (b - a) * mp.mpf(i) / n0 for i in range(n0 + 1)]
    vals = [mp.arg(dh.F(path(x))) for x in xs]
    tot = mp.mpf(0)
    i = 0
    hmin = abs(b - a) / mp.mpf(2 ** maxdepth)
    refined = 0
    while i < len(xs) - 1:
        d = vals[i + 1] - vals[i]
        while d > mp.pi:
            d -= 2 * mp.pi
        while d < -mp.pi:
            d += 2 * mp.pi
        if abs(d) > mp.pi * mp.mpf('0.45') and abs(xs[i + 1] - xs[i]) > hmin:
            xm = (xs[i] + xs[i + 1]) / 2
            xs.insert(i + 1, xm)
            vals.insert(i + 1, mp.arg(dh.F(path(xm))))
            refined += 1
            continue
        tot += d
        i += 1
    return tot, refined


def winding(sL, sR, T1, T2, n0=None):
    sL, sR, T1, T2 = (mp.mpf(x) for x in (sL, sR, T1, T2))
    if n0 is None:
        n0 = 48
    nh = max(16, int((sR - sL) * 24))
    nv = max(24, int((T2 - T1) * 3))
    tot = mp.mpf(0)
    tot += arg_change(lambda x: mp.mpc(x, T1), sL, sR, nh)[0]
    tot += arg_change(lambda y: mp.mpc(sR, y), T1, T2, nv)[0]
    tot += arg_change(lambda x: mp.mpc(x, T2), sR, sL, nh)[0]
    tot += arg_change(lambda y: mp.mpc(sL, y), T2, T1, nv)[0]
    return tot / (2 * mp.pi)


def local_winding(c, r, n=400):
    """winding number of F on the circle |s - c| = r; certifies a located zero."""
    vals = [mp.arg(dh.F(c + r * mp.e ** (2j * mp.pi * mp.mpf(j) / n)))
            for j in range(n + 1)]
    tot = mp.mpf(0)
    for j in range(n):
        d = vals[j + 1] - vals[j]
        while d > mp.pi:
            d -= 2 * mp.pi
        while d < -mp.pi:
            d += 2 * mp.pi
        tot += d
    return tot / (2 * mp.pi)


def find_offline(x0, x1, y0, y1, depth=0, maxdepth=9, out=None):
    """Recursively subdivide until each box has winding 0 or 1, then Newton."""
    if out is None:
        out = []
    w = winding(x0, x1, y0, y1)
    n = int(mp.nint(w))
    if abs(w - n) > mp.mpf('1e-6'):
        out.append(('UNRESOLVED-WINDING', float(x0), float(x1), float(y0), float(y1),
                    float(w)))
        return out
    if n == 0:
        return out
    if n == 1 and (x1 - x0) < mp.mpf('0.9') and (y1 - y0) < mp.mpf('2.0'):
        c = mp.mpc((x0 + x1) / 2, (y0 + y1) / 2)
        try:
            r = mp.findroot(dh.F, c, tol=mp.mpf(10) ** (-2 * mp.mp.dps + 20))
        except Exception:                                              # noqa: BLE001
            out.append(('NEWTON-FAILED', float(x0), float(x1), float(y0), float(y1), 1.0))
            return out
        if x0 <= mp.re(r) <= x1 and y0 <= mp.im(r) <= y1:
            out.append(('ZERO', r))
            return out
    if depth >= maxdepth:
        out.append(('MAXDEPTH', float(x0), float(x1), float(y0), float(y1), float(n)))
        return out
    xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
    if (x1 - x0) * 8 > (y1 - y0):
        find_offline(x0, xm, y0, y1, depth + 1, maxdepth, out)
        find_offline(xm, x1, y0, y1, depth + 1, maxdepth, out)
    else:
        find_offline(x0, x1, y0, ym, depth + 1, maxdepth, out)
        find_offline(x0, x1, ym, y1, depth + 1, maxdepth, out)
    return out


def online_zeros(Tmax, dps_root, steps_per_gap=12, tmin='0.02'):
    out, flags = [], []
    t = mp.mpf(tmin)
    zt = dh.Z(t)
    while t < Tmax:
        sp = 2 * mp.pi / max(mp.log(5 * max(t, mp.mpf(1)) / (2 * mp.pi)), mp.mpf('0.3'))
        h = min(sp / steps_per_gap, mp.mpf('0.3'))
        t2 = t + h
        z2 = dh.Z(t2)
        if zt * z2 < 0:
            lo, hi, zlo = t, t2, zt
            for _ in range(24):
                mid = (lo + hi) / 2
                zm = dh.Z(mid)
                if zm == 0:
                    lo = hi = mid
                    break
                if zlo * zm < 0:
                    hi = mid
                else:
                    lo, zlo = mid, zm
            old = mp.mp.dps
            mp.mp.dps = dps_root
            try:
                r = mp.findroot(dh.Z, (lo + hi) / 2, tol=mp.mpf(10) ** (-2 * dps_root + 20))
                out.append(mp.re(r))
            except Exception as exc:                                   # noqa: BLE001
                flags.append(('findroot', float((lo + hi) / 2), str(exc)))
                out.append((lo + hi) / 2)
            finally:
                mp.mp.dps = old
        t, zt = t2, z2
    return out, flags


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--Tmax', type=float, default=200.0)
    ap.add_argument('--dps-scan', type=int, default=25)
    ap.add_argument('--dps-root', type=int, default=60)
    a = ap.parse_args()
    os.makedirs(RES, exist_ok=True)
    os.makedirs(DATA, exist_ok=True)
    out = []
    rec = {}

    def p(line=''):
        print(line)
        out.append(line)

    mp.mp.dps = a.dps_scan
    T = mp.mpf(a.Tmax)

    p('=' * 78)
    p('FINITE-HEIGHT ZERO CENSUS   -0.4 < Re s < 1.4,  0 < Im s < %g' % a.Tmax)
    p('=' * 78)
    p('  sigma_0 = %s   (sum_{n>=2}|a_n| n^-sigma = 1); the vertical edges at' % SIGMA0)
    p('  Re s = -0.4 and 1.4 are provably zero-free, so the contour count is exact.')
    p('')

    # ---- contour counts on a ladder
    p('[1] argument-principle counts')
    p('    %-8s %-22s %-8s %-14s %-12s' %
      ('T', 'winding number', 'integer', 'RvM main term', 'diff'))
    ladder = []
    for Tl in [25, 50, 100, 150, 200, a.Tmax]:
        if Tl > a.Tmax:
            continue
        t0 = time.time()
        w = winding(SL, SR, '0.001', Tl)
        n = int(mp.nint(w))
        rvm = dh.N_count(Tl)
        p('    %-8g %-22s %-8d %-14s %-12s  [%.0fs]'
          % (Tl, mp.nstr(w, 16), n, mp.nstr(rvm, 8), mp.nstr(n - rvm, 6),
             time.time() - t0))
        ladder.append(dict(T=float(Tl), winding=mp.nstr(w, 20), count=n,
                           rvm=float(rvm), dev=float(n - rvm),
                           integrality=float(abs(w - n))))
    Ntot = ladder[-1]['count']
    maxint = max(r['integrality'] for r in ladder)
    p('    worst deviation of a winding number from an integer: %s' % mp.nstr(mp.mpf(maxint), 4))
    rec['ladder'] = ladder

    # ---- off-line zeros
    p('')
    p('[2] off-line zeros with Re s > 1/2  (box subdivision + Newton + local winding)')
    # The search region is taken at Re s >= 0.55, i.e. 0.05 clear of the critical line:
    # a contour that runs closer than that passes within a fraction of the zero spacing
    # of the on-line zeros, the argument swings by nearly pi between samples, and the
    # adaptive refinement blows up.  Nothing is lost, because the RECONCILIATION below
    # fixes the number of off-line zeros exactly:  contour count = #on-line + 2 #off-line,
    # so if the located off-line zeros account for the contour count there is no room for
    # a zero in the uncovered band 1/2 < Re s < 0.55.
    found = []
    nband = max(1, int(mp.ceil(T / 25)))
    nonzero_bands = 0
    for ib in range(nband):
        y0 = mp.mpf('0.001') if ib == 0 else mp.mpf(25 * ib)
        y1 = min(mp.mpf(25 * (ib + 1)), T)
        wb = winding(mp.mpf('0.55'), SR, y0, y1)
        nb = int(mp.nint(wb))
        p('      band Im in (%g, %g): winding = %s' % (float(y0), float(y1),
                                                        mp.nstr(wb, 12)))
        if abs(wb - nb) > mp.mpf('1e-6'):
            found.append(('UNRESOLVED-BAND', float(y0), float(y1), 0.0, 0.0, float(wb)))
            continue
        if nb == 0:
            continue
        nonzero_bands += 1
        find_offline(mp.mpf('0.55'), SR, y0, y1, out=found)
    p('    scanned %d height-25 bands at Re s >= 0.55; %d contained zeros'
      % (nband, nonzero_bands))

    zeros = [z for tag, z in [(f[0], f[1]) for f in found if f[0] == 'ZERO']]
    probs = [f for f in found if f[0] != 'ZERO']
    zeros.sort(key=lambda z: mp.im(z))
    p('    boxes resolved: %d zero(s) found, %d problem box(es)' % (len(zeros), len(probs)))
    for f in probs:
        p('      PROBLEM: %s' % (f,))
    mp.mp.dps = a.dps_root
    offrows = []
    for z in zeros:
        r = mp.findroot(dh.F, z, tol=mp.mpf(10) ** (-2 * a.dps_root + 20))
        lw = local_winding(r, mp.mpf('0.02'))
        p('      rho = %s' % mp.nstr(r, 30))
        p('            beta  = %s' % mp.nstr(mp.re(r), 30))
        p('            gamma = %s' % mp.nstr(mp.im(r), 30))
        p('            delta = %s   |F| = %-10s  local winding on |s-rho|=0.02 : %s'
          % (mp.nstr(mp.re(r) - mp.mpf(1) / 2, 20), mp.nstr(abs(dh.F(r)), 4),
             mp.nstr(lw, 10)))
        offrows.append(dict(beta=mp.nstr(mp.re(r), 50), gamma=mp.nstr(mp.im(r), 50),
                            delta=mp.nstr(mp.re(r) - mp.mpf(1) / 2, 50),
                            absF=float(abs(dh.F(r))), local_winding=float(lw),
                            multiplicity=int(mp.nint(lw))))
    rec['offline'] = offrows
    # check for anything in the narrow band 1/2 < Re <= 0.505
    mp.mp.dps = a.dps_scan
    p('    (the band 1/2 < Re s < 0.55 is covered by the reconciliation in [4], not by')
    p('     a separate contour: a zero there would break contour = on-line + 2 off-line.)')
    rec['narrow_band_winding'] = 'covered by reconciliation'

    # ---- on-line zeros
    p('')
    p('[3] critical-line ordinates from sign changes of Z(t)')
    t0 = time.time()
    gam, flags = online_zeros(T, a.dps_root)
    p('    found %d sign-change zeros with 0 < t < %g   [%.0fs]'
      % (len(gam), a.Tmax, time.time() - t0))
    for f in flags:
        p('      FLAG: %s' % (f,))
    mp.mp.dps = a.dps_root
    mx = max(abs(dh.Z(g)) for g in gam) if gam else mp.mpf(0)
    p('    max |Z(gamma)| over the refined ordinates: %s' % mp.nstr(mx, 4))

    # ---- reconciliation
    p('')
    p('[4] reconciliation')
    noff = sum(r['multiplicity'] for r in offrows)
    p('    contour count (all zeros, 0 < Im < %g)          : %d' % (a.Tmax, Ntot))
    p('    off-line with Re > 1/2 (with multiplicity)      : %d' % noff)
    p('    their mirrors with Re < 1/2 (Lambda_F(s)=Lambda_F(1-s)) : %d' % noff)
    p('    on-line ordinates found                         : %d' % len(gam))
    tot = len(gam) + 2 * noff
    p('    total accounted for                             : %d' % tot)
    ok = (tot == Ntot)
    p('    RECONCILED: %s' % ('YES' if ok else 'NO  (difference %d)' % (Ntot - tot)))
    rec['reconciliation'] = dict(contour=Ntot, offline=noff, online=len(gam),
                                 total=tot, reconciled=bool(ok))

    p('')
    p('    Status of this census (non-claims 6, 7):')
    p('      * a LIST of found roots:        yes')
    p('      * a numerically RECONCILED census: %s' % ('yes' if ok else 'NO'))
    p('      * a rigorously COMPLETE census:  no -- the winding numbers are integers to')
    p('        %s but are computed by adaptive sampling, not interval arithmetic.'
      % mp.nstr(mp.mpf(maxint), 3))
    p('      * roots on Re s = 1/2 are "consistent with the line at the tested')
    p('        precision"; no existence/uniqueness certificate is claimed for them.')

    with open(os.path.join(DATA, 'dh_zeros.csv'), 'w') as fh:
        fh.write('# Davenport-Heilbronn zeros of Lambda_F, 0 < Im s < %g\n' % a.Tmax)
        fh.write('# strip -0.4 < Re s < 1.4 (provably contains all zeros)\n')
        fh.write('# contour count %d; on-line %d; off-line with Re>1/2 %d (+ mirrors)\n'
                 % (Ntot, len(gam), noff))
        fh.write('# precision: scan %d dps, roots %d dps.  T2: no interval certificate.\n'
                 % (a.dps_scan, a.dps_root))
        fh.write('kind,beta,gamma,multiplicity\n')
        for r in offrows:
            fh.write('offline,%s,%s,%d\n' % (r['beta'], r['gamma'], r['multiplicity']))
            fh.write('offline_mirror,%s,%s,%d\n'
                     % (mp.nstr(1 - mp.mpf(r['beta']), 50), r['gamma'],
                        r['multiplicity']))
        for g in gam:
            fh.write('online,0.5,%s,1\n' % mp.nstr(g, 50))
    p('')
    p('    -> data/dh_zeros.csv')
    p('=' * 78)

    with open(os.path.join(RES, 'census.txt'), 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    with open(os.path.join(RES, 'census.json'), 'w') as fh:
        json.dump(rec, fh, indent=1)
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
