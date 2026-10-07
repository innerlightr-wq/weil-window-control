#!/usr/bin/env python3
r"""Stage 5(c) -- P4 against the right baseline.

Stage 5(b) measured the fraction of real dip measure lying in the {2,3} alignment
set S_23(1/2) and compared it with the unconditional density of S_23.  That ratio
(~3.8) is not the test: conditioning on being in a dip at all already forces the
heavy blocks to align.  The N1 model predicts that conditional probability exactly:

    Pr[ P_23 >= c (a2+a3)  |  P >= H(t) ]

on the torus, evaluated at the actual dip locations and weighted by dip length.
P4 holds iff the real fraction matches this, and fails iff the {2,3} block carries
measurably more (or less) than the null allows.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from symbol import Symbol, reduced_arg
from scan import scan
from stage3_nulls import surrogate

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
C = 0.5
NMC = 2_000_000
N_SUR = 20


def torus(sym, n=NMC, seed=0):
    """(P, P_23) sampled from the N1 torus law of this symbol."""
    rng = np.random.default_rng(seed)
    bl = {}
    for p, f, w in zip(sym.prime, sym.freqs, sym.weights):
        bl.setdefault(int(p), []).append((float(f), float(w)))
    P = np.zeros(n); P23 = np.zeros(n)
    A23 = 0.0
    for p, terms in bl.items():
        th = rng.uniform(0, 2 * np.pi, n)
        k0 = min(f for f, _ in terms)
        v = np.zeros(n)
        for f, w in terms:
            v += w * np.cos(round(f / k0) * th)
        P += v
        if p in (2, 3):
            P23 += v; A23 += sum(w for _, w in terms)
    return P, P23, A23


def predict(P, P23, A23, hs, lens):
    """N1 Pr[P23 >= C*A23 | P >= h], dip-length weighted over the given h's."""
    o = np.argsort(P); Ps = P[o]; hit = (P23[o] >= C * A23).astype(np.int64)
    suf_hit = np.concatenate([np.cumsum(hit[::-1])[::-1], [0]])
    suf_n = np.arange(Ps.size, -1, -1)
    i = np.searchsorted(Ps, hs, side='left')
    with np.errstate(invalid='ignore', divide='ignore'):
        q = np.where(suf_n[i] > 0, suf_hit[i] / np.maximum(suf_n[i], 1), np.nan)
    ok = np.isfinite(q)
    return float(np.sum(q[ok] * lens[ok]) / np.sum(lens[ok])), float(np.mean(suf_n[i] > 50))


def measure_in_S23(sym, dips, a2, a3, l2, l3):
    tot = cov = 0.0
    hs, lens = [], []
    for a, b in dips:
        n = max(64, int((b - a) / 0.01))
        tm = np.linspace(a, b, n)
        p23 = a2 * np.cos(reduced_arg(tm, l2)) + a3 * np.cos(reduced_arg(tm, l3))
        tot += (b - a); cov += (b - a) * float((p23 >= C * (a2 + a3)).mean())
        hs.append(float(sym.H(np.array([0.5 * (a + b)]))[0])); lens.append(b - a)
    return cov / tot, np.array(hs), np.array(lens)


def main():
    L2, L3 = np.log(2), np.log(3)
    rows = []
    print("  %-6s %9s %9s %9s %-18s %s" %
          ("L", "real", "N1 pred", "ratio", "N2 95% band", "P4"))
    for L in [0.6, 0.7, 0.9, 1.0, 1.19, 1.2, 1.4, 1.6]:
        sym = Symbol(L)
        a2, a3 = sym.weights[0], sym.weights[1]
        dips = scan(sym, 0.0)['dips']
        real, hs, lens = measure_in_S23(sym, dips, a2, a3, L2, L3)
        P, P23, A23 = torus(sym, seed=int(100 * L))
        pred, cov = predict(P, P23, A23, hs, lens)
        # N2: the same statistic on surrogates, using each surrogate's own 2- and
        # 3-blocks (the two heaviest), so the comparison is like for like
        rng = np.random.default_rng(99 + int(100 * L))
        vals = []
        for _ in range(N_SUR):
            s2 = surrogate(sym, rng, spread=0.3)
            s2.prime = sym.prime.copy(); s2.expo = sym.expo.copy()
            f2, f3 = s2.freqs[0], s2.freqs[1]
            d2 = scan(s2, 0.0)['dips']
            if d2:
                vals.append(measure_in_S23(s2, d2, a2, a3, f2, f3)[0])
        band = [float(np.quantile(vals, q)) for q in (0.025, 0.5, 0.975)]
        ok = band[0] <= real <= band[2]
        print("  %-6s %9.4f %9.4f %9.3f [%.4f, %.4f]   %s" %
              (L, real, pred, real / pred, band[0], band[2], "in band" if ok else "OUT"))
        rows.append(dict(L=L, real=real, n1_pred=pred, ratio=real / pred,
                         n1_support=cov, n2_band=band, in_band=bool(ok),
                         n_dips=len(dips)))
    json.dump(rows, open(os.path.join(ROOT, 'data', 'stage5_p4.json'), 'w'), indent=1)
    print("\nwrote data/stage5_p4.json")


if __name__ == '__main__':
    main()
