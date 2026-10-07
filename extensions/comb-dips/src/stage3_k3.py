#!/usr/bin/env python3
r"""K3 protocol for the dip-measure excess.

The preregistered K3 requires, before calling anything structure:
  (1) three fresh surrogate seeds,
  (2) double precision on the affected L,
  (3) the alternative surrogate (log p +- 0.1 instead of +- 0.3).

Added on my own initiative, because it is the most likely mundane explanation:
  (4) a REGION SPLIT. Most dip measure sits at small t, where H(t) = Re psi - log pi is
      small and almost everything is a dip. That is a transient, not the equidistributed
      regime the null models describe, and the surrogate changes log 2 = 0.693 by up to
      +-0.3 -- a 43% change in the lowest frequency, which is exactly what controls the
      small-t structure. If the excess lives entirely at small t, it is an artefact of
      including a region where the null is not applicable.
"""
import sys, os, csv, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from symbol import Symbol
from scan import scan
from stage3_nulls import surrogate

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def measure_split(dips, tsplit):
    lo = sum(min(b, tsplit) - a for a, b in dips if a < tsplit)
    hi = sum(b - max(a, tsplit) for a, b in dips if b > tsplit)
    return lo, hi


def run(L, beta, n_sur=200, spread=0.3, seed=0):
    s = Symbol(L)
    r = scan(s, beta)
    tsplit = 2 * np.pi * np.exp(s.A / 2)          # geometric midpoint on the H scale
    rl, rh = measure_split(r['dips'], tsplit)
    rng = np.random.default_rng(seed)
    tot, los, his = [], [], []
    for _ in range(n_sur):
        sg = surrogate(s, rng, spread)
        q = scan(sg, beta)
        a, b = measure_split(q['dips'], tsplit)
        tot.append(q['dip_measure']); los.append(a); his.append(b)
    f = lambda v, x: float((np.asarray(v) >= x).mean())
    pct = lambda v: np.percentile(v, [2.5, 97.5])
    return dict(L=L, beta=beta, n_sur=n_sur, spread=spread, seed=seed, tsplit=tsplit,
                real_total=r['dip_measure'], real_lo=rl, real_hi=rh,
                band_total=pct(tot).tolist(), band_lo=pct(los).tolist(),
                band_hi=pct(his).tolist(),
                p_total=f(tot, r['dip_measure']), p_lo=f(los, rl), p_hi=f(his, rh))


if __name__ == '__main__':
    cells = [(1.0, 0.25), (1.6, 0.0), (1.6, 0.25), (1.2, 0.0), (1.4, 0.0)]
    out = []
    print("(1) three fresh seeds, and (3) the +-0.1 surrogate, on the affected cells")
    print("%5s %5s %8s %6s %11s %22s %8s" %
          ("L", "beta", "spread", "seed", "real meas", "surrogate 95%", "p(>=real)"))
    for L, beta in cells:
        for spread in (0.3, 0.1):
            for seed in (11, 22, 33):
                d = run(L, beta, 200, spread, seed)
                print("%5s %5.2f %8.1f %6d %11.3f [%9.3f,%9.3f] %8.3f" %
                      (L, beta, spread, seed, d['real_total'],
                       d['band_total'][0], d['band_total'][1], d['p_total']))
                out.append(d)
    print("\n(4) REGION SPLIT at t = 2 pi e^{A/2}: is the excess only at small t?")
    print("%5s %5s %10s %11s %22s %8s %11s %22s %8s" %
          ("L", "beta", "t split", "real low-t", "surrogate 95% low-t", "p", "real high-t",
           "surrogate 95% high-t", "p"))
    for L, beta in cells:
        d = [x for x in out if x['L'] == L and x['beta'] == beta
             and x['spread'] == 0.3 and x['seed'] == 11][0]
        print("%5s %5.2f %10.4g %11.3f [%9.3f,%9.3f] %8.3f %11.3f [%9.3f,%9.3f] %8.3f" %
              (L, beta, d['tsplit'], d['real_lo'], d['band_lo'][0], d['band_lo'][1],
               d['p_lo'], d['real_hi'], d['band_hi'][0], d['band_hi'][1], d['p_hi']))
    json.dump(out, open(os.path.join(ROOT, 'data', 'stage3_k3.json'), 'w'),
              indent=1, default=float)
    print("\nwrote data/stage3_k3.json")
