#!/usr/bin/env python3
r"""Stage 5: calibration on the exactly-solvable two-prime case.

For 0.549 < L < 0.693 the active set is {2,3} and

    Psi_L(t) = H(t) - a2 cos(t log 2) - a3 cos(t log 3),
    a2 = 2 log 2 / sqrt 2,  a3 = 2 log 3 / sqrt 3.

Here the beat picture is exact: the envelope of the two-frequency sum has period
2 pi / log(3/2) = 15.4959..., and the near-recurrences (times where BOTH cosines are
near 1, so P is near its maximum a2 + a3) are governed by the continued-fraction
convergents of log 3 / log 2.

We then switch primes on one at a time and measure how much dip structure the {2,3}
model still explains.
"""
import sys, os, json, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from fractions import Fraction
from symbol import Symbol, reduced_arg
from scan import scan

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
L2, L3 = np.log(2), np.log(3)


def convergents(x, n=10):
    out, a, p0, q0, p1, q1 = [], x, 1, 0, 0, 1
    for _ in range(n):
        ai = int(np.floor(a))
        p0, p1 = p1, ai * p1 + p0
        q0, q1 = q1, ai * q1 + q0
        out.append(Fraction(p1, q1))
        frac = a - ai
        if frac < 1e-15:
            break
        a = 1.0 / frac
    return out


def main():
    print("=" * 76)
    print("Stage 5(a): the {2,3} case -- beat structure, exactly")
    print("=" * 76)
    s = Symbol(0.6)
    a2, a3 = s.weights[0], s.weights[1]
    assert s.n.tolist() == [2, 3]
    env_T = 2 * np.pi / (L3 - L2)
    print("  a2 = %.9f   a3 = %.9f   A = a2+a3 = %.9f (= A_L %.9f)"
          % (a2, a3, a2 + a3, s.A))
    print("  envelope period 2 pi / log(3/2) = %.6f   (brief: approx 15.50)" % env_T)
    cf = convergents(L3 / L2, 8)
    print("  log3/log2 = %.9f ; convergents: %s"
          % (L3 / L2, ", ".join(str(c) for c in cf)))

    # exact alignment times: t log2 = 2 pi m  AND  t log3 = 2 pi n
    print("\n  best alignments t = 2 pi m / log 2 within the scan range, and P there:")
    r = scan(s, 0.0)
    hi = r['hi']
    print("    %4s %12s %12s %12s %10s" % ("m", "t", "P(t)", "A - P", "n/m vs log3/log2"))
    rows = []
    for m in range(1, int(hi * L2 / (2 * np.pi)) + 2):
        t = 2 * np.pi * m / L2
        if t > hi:
            break
        P = float(s.P(np.array([t]))[0])
        nn = t * L3 / (2 * np.pi)
        print("    %4d %12.6f %12.6f %12.6f %10.5f" % (m, t, P, s.A - P, nn / m))
        rows.append(dict(m=m, t=t, P=P, deficit=s.A - P, ratio=nn / m))

    # dip catalogue vs the beat prediction
    dips = r['dips']
    print("\n  %d dips in [2pi, %.4f]; checking each against the beat picture:" %
          (len(dips), hi))
    print("    %12s %12s %12s %12s" % ("dip start", "dip end", "min Psi", "nearest 2pi m/log2"))
    for a, b in dips:
        tm = np.linspace(a, b, 2001)
        psi = s.Psi(tm)
        k = int(np.argmin(psi))
        mnear = round(tm[k] * L2 / (2 * np.pi))
        print("    %12.6f %12.6f %12.6f %12.6f" %
              (a, b, psi[k], 2 * np.pi * mnear / L2))

    print("\n" + "=" * 76)
    print("Stage 5(b): how much of the dip structure the {2,3} beat still explains")
    print("=" * 76)
    print("""  The {2,3} DIP set is useless as a model past L=0.7: it is capped at
  T1(0.55) = %.2f, while the dips of L=1.6 run out to 4.7e6, so the covered
  fraction would be zero for trivial reasons.  The right comparison is with the
  {2,3} ALIGNMENT set, which is defined for all t:

      S_23(c) = { t : P_23(t) >= c (a2+a3) },   P_23 = a2 cos(t log2) + a3 cos(t log3).

  We report, for c = 1/2, the fraction of the real dip measure lying in S_23(c),
  against the density of S_23(c) itself (chance level, from the two-torus).  Their
  ratio is the enrichment: 1 = the {2,3} beat explains nothing beyond its own size.
""" % (2 * np.pi * np.exp(Symbol(0.55).A)))

    s23 = Symbol(0.6)
    a2, a3 = s23.weights[0], s23.weights[1]
    C = 0.5

    def P23(t):
        return a2 * np.cos(reduced_arg(t, L2)) + a3 * np.cos(reduced_arg(t, L3))

    # chance level: density of S_23 on the 2-torus (frequencies incommensurable)
    rng = np.random.default_rng(7)
    th = rng.uniform(0, 2 * np.pi, (2, 4_000_000))
    dens = float((a2 * np.cos(th[0]) + a3 * np.cos(th[1]) >= C * (a2 + a3)).mean())
    print("  chance density of S_23(1/2) = %.5f\n" % dens)

    print("  %-6s %-30s %8s %11s %9s %10s" %
          ("L", "active n", "n dips", "measure", "in S_23", "enrichment"))
    out = []
    for L in [0.6, 0.7, 0.9, 1.0, 1.19, 1.2, 1.4, 1.6]:
        sL = Symbol(L)
        rr = scan(sL, 0.0)
        tot = cov = 0.0
        for a, b in rr['dips']:
            n = max(64, int((b - a) / 0.01))
            tm = np.linspace(a, b, n)
            tot += (b - a)
            cov += (b - a) * float((P23(tm) >= C * (a2 + a3)).mean())
        frac = cov / tot if tot else float('nan')
        print("  %-6s %-30s %8d %11.4f %9.4f %10.3f" %
              (L, str(sL.n.tolist())[:30], rr['n_dips'], rr['dip_measure'],
               frac, frac / dens))
        out.append(dict(L=L, active=str(sL.n.tolist()), n_dips=rr['n_dips'],
                        measure=rr['dip_measure'], frac_in_S23=frac,
                        chance=dens, enrichment=frac / dens))

    with open(os.path.join(ROOT, 'data', 'stage5_twoprime.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
    json.dump(rows, open(os.path.join(ROOT, 'data', 'stage5_alignments.json'), 'w'),
              indent=1, default=float)
    print("\n  wrote data/stage5_twoprime.csv, data/stage5_alignments.json")


if __name__ == '__main__':
    main()
