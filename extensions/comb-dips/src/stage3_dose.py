#!/usr/bin/env python3
r"""Dose-response: does the dip-measure excess fade as the surrogate spread drops below
the comma detunings?

Verdict rule (set before running):
  excess fades below the comma scales  -> EXPLAINED by comma-type near-commensurabilities
  excess persists well below all of them -> CANDIDATE ANOMALY (T2)

CONFOUND, controlled for. As spread -> 0 the surrogate tends to the real symbol, so the
excess must vanish for a trivial reason too. To separate that from the comma mechanism we
also measure OVERLAP: the fraction of real dip measure covered by surrogate dips. If the
dip patterns are still decorrelated at the spread where the excess vanishes, the vanishing
is informative; if the patterns have converged, it is trivial. Random-overlap baseline for
reference is (surrogate measure)/(range).
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from symbol import Symbol
from scan import scan
from stage3_nulls import surrogate
from stage3_k3 import measure_split

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SPREADS = [0.3, 0.1, 0.03, 0.01, 0.003]
CELLS = [(1.0, 0.25), (1.6, 0.0), (1.6, 0.25)]


def overlap_fraction(real_dips, sur_dips):
    """fraction of real dip measure lying inside surrogate dips"""
    tot = sum(b - a for a, b in real_dips)
    if tot == 0 or not sur_dips:
        return 0.0
    S = np.array(sur_dips)
    cov = 0.0
    for a, b in real_dips:
        j0 = np.searchsorted(S[:, 1], a)
        for c, d in S[j0:]:
            if c >= b:
                break
            lo, hi = max(a, c), min(b, d)
            if hi > lo:
                cov += hi - lo
    return cov / tot


def main(n_sur=100):
    out = []
    print("dose-response: real-minus-surrogate dip-measure excess vs surrogate spread")
    print("%5s %5s %8s %11s %11s %9s %9s %10s %10s" %
          ("L", "beta", "spread", "real meas", "sur median", "excess", "excess%",
           "p(>=real)", "overlap"))
    for L, beta in CELLS:
        s = Symbol(L)
        r = scan(s, beta)
        real = r['dip_measure']
        rng_len = r['hi'] - r['lo']
        tsplit = 2 * np.pi * np.exp(s.A / 2)
        real_lo_t, _ = measure_split(r['dips'], tsplit)
        for sp in SPREADS:
            rng = np.random.default_rng(777 + int(1000 * L) + int(100 * beta))
            ms, ov, mlo = [], [], []
            for _ in range(n_sur):
                sg = surrogate(s, rng, sp)
                q = scan(sg, beta)
                ms.append(q['dip_measure'])
                a_, _b = measure_split(q['dips'], tsplit)
                mlo.append(a_)
                ov.append(overlap_fraction(r['dips'], q['dips']))
            ms = np.array(ms); ov = np.array(ov); mlo = np.array(mlo)
            med = float(np.median(ms))
            exc = real - med
            p = float((ms >= real).mean())
            base = med / rng_len          # random-overlap baseline
            print("%5s %5.2f %8.3f %11.3f %11.3f %9.3f %9.2f %10.3f %10.4f" %
                  (L, beta, sp, real, med, exc, 100 * exc / med, p, float(ov.mean())))
            out.append(dict(L=L, beta=beta, spread=sp, n_sur=n_sur, real=real,
                            sur_median=med, sur_lo=float(np.percentile(ms, 2.5)),
                            sur_hi=float(np.percentile(ms, 97.5)),
                            excess=exc, excess_pct=100 * exc / med, p_ge_real=p,
                            overlap_mean=float(ov.mean()),
                            overlap_random_baseline=float(base),
                            real_low_t=real_lo_t, sur_low_t_median=float(np.median(mlo)),
                            excess_low_t=real_lo_t - float(np.median(mlo)),
                            p_low_t=float((mlo >= real_lo_t).mean())))
    json.dump(out, open(os.path.join(ROOT, 'data', 'stage3_dose.json'), 'w'),
              indent=1, default=float)
    print("\nwrote data/stage3_dose.json")


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 100)
