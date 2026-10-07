#!/usr/bin/env python3
"""Boundary re-verification of t0 at high precision + a float64 safety-margin test.

python-flint/Arb is not installed in this environment, so interval arithmetic is replaced
by mpmath at 50 digits (an independent implementation path from the numba kernel) plus an
explicit float64 error budget. This is flagged in REPORT.md: the scan is T2 with a
controlled error model, NOT a formal interval certificate.
"""
import sys, os, csv, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import mpmath as mp
from symbol import Symbol
from scan import scan

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
mp.mp.dps = 50


def psi_mp(s, t):
    H = mp.re(mp.digamma(mp.mpf(1) / 4 + mp.mpc(0, 1) * mp.mpf(t) / 2)) - mp.log(mp.pi)
    return H - sum(mp.mpf(w) * mp.cos(mp.mpf(t) * mp.mpf(f))
                   for w, f in zip(s.weights, s.freqs))


if __name__ == '__main__':
    rows = list(csv.DictReader(open(os.path.join(ROOT, 'data', 'stage2_t0.csv'))))
    out = []
    print("%-6s %6s %14s %16s %16s %10s" %
          ("L", "beta", "t0_upper", "Psi(last dip mid)", "min Psi beyond t0", "verdict"))
    for r in rows:
        if float(r['beta']) != 0.0:
            continue
        L = float(r['L']); beta = 0.0
        s = Symbol(L)
        dips = np.load(os.path.join(ROOT, 'data', 'dips_L%s_b%s.npy' % (r['L'], r['beta'])))
        a, b = dips[-1]
        mid = 0.5 * (a + b)
        psi_mid = float(psi_mp(s, mid))          # must be < beta: a dip really is there
        # beyond t0_upper: sample densely in mpmath and confirm Psi > beta
        hi = float(r['scan_hi']); t0u = float(r['t0_upper'])
        ts = np.linspace(t0u + 1e-6, hi, 400)
        beyond = min(float(psi_mp(s, t)) for t in ts)
        ok = (psi_mid < beta) and (beyond > beta)
        print("%-6s %6.2f %14.6g %16.6f %16.6f %10s" %
              (r['L'], beta, t0u, psi_mid, beyond, "OK" if ok else "CHECK"))
        out.append(dict(L=L, beta=beta, t0_upper=t0u, psi_last_dip_mid=psi_mid,
                        min_psi_beyond_sampled=beyond, n_samples_beyond=len(ts), ok=bool(ok)))

    print("\nfloat64 safety-margin test: rerun with the certification threshold tightened")
    print("by 1e-12 (>> the verified 4e-15 evaluation error). t0 must not move.")
    for L in [0.9, 1.19, 1.4]:
        s = Symbol(L)
        r1 = scan(s, 0.0)
        r2 = scan(s, 1e-12)        # demand Psi >= 1e-12 instead of >= 0 to discard
        same = abs(r1['t0_upper'] - r2['t0_upper']) < 1e-6
        print("  L=%-5s t0(beta=0) = %.9g   t0(beta=1e-12) = %.9g   unchanged: %s"
              % (L, r1['t0_upper'], r2['t0_upper'], same))
        out.append(dict(L=L, margin_test=True, t0_beta0=r1['t0_upper'],
                        t0_beta_eps=r2['t0_upper'], unchanged=bool(same)))
    json.dump(out, open(os.path.join(ROOT, 'data', 'stage2_verification.json'), 'w'),
              indent=1, default=float)
    print("\nwrote data/stage2_verification.json")
