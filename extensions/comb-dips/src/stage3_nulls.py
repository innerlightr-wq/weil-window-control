#!/usr/bin/env python3
"""Stage 3: N1 (torus measure) predictions and N2 (random-frequency) surrogates."""
import sys, os, csv, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from symbol import Symbol, STAGE2_L
from scan import scan
import null_n1 as N1

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BETAS = [0.0, 0.25, 0.5]


def n1_table():
    rows = list(csv.DictReader(open(os.path.join(ROOT, 'data', 'stage2_t0.csv'))))
    out = []
    print("N1 -- torus-measure predictions vs the certified scan")
    print("%-6s %5s %12s %12s %8s %12s %12s %8s" %
          ("L", "beta", "t0 real", "t0 N1", "ratio", "meas real", "meas N1", "ratio"))
    for r in rows:
        L, beta = float(r['L']), float(r['beta'])
        s = Symbol(L)
        x, pdf = N1.tail_pdf(s)
        pi = N1.tail_prob(x, pdf)
        dips = np.load(os.path.join(ROOT, 'data', 'dips_L%s_b%s.npy' % (r['L'], r['beta'])))
        lens = dips[:, 1] - dips[:, 0]
        med = float(np.median(lens)) if lens.size else 0.0
        hi = float(r['scan_hi'])
        t0_pred = N1.predicted_t0(s, beta, pi, hi, med)
        t, dens, tail = N1.predict(s, beta, pi, hi)
        meas_pred = float(tail[0])
        t0_real = float(r['t0_upper']); meas_real = float(r['dip_measure'])
        print("%-6s %5.2f %12.5g %12.5g %8.3f %12.4f %12.4f %8.3f" %
              (L, beta, t0_real, t0_pred, t0_pred / t0_real, meas_real, meas_pred,
               meas_pred / meas_real))
        out.append(dict(L=L, beta=beta, t0_real=t0_real, t0_N1=t0_pred,
                        t0_ratio=t0_pred / t0_real, measure_real=meas_real,
                        measure_N1=meas_pred, measure_ratio=meas_pred / meas_real,
                        median_dip_len=med))
    with open(os.path.join(ROOT, 'data', 'stage3_n1.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
    return out


def surrogate(sym, rng, spread=0.3):
    """Replace each log p by log p + U(-spread, spread); prime powers keep k*(new log p)."""
    newlogp = {}
    for p in np.unique(sym.prime):
        newlogp[int(p)] = float(np.log(p) + rng.uniform(-spread, spread))
    freqs = np.array([sym.expo[i] * newlogp[int(sym.prime[i])] for i in range(len(sym.n))])
    return Symbol(sym.L, freqs=freqs, weights=sym.weights.copy(),
                  label=f"surrogate L={sym.L}")


def n2_table(n_sur=200, betas=(0.0,), spread=0.3, seed=0, Ls=None):
    Ls = Ls or STAGE2_L
    rows = list(csv.DictReader(open(os.path.join(ROOT, 'data', 'stage2_t0.csv'))))
    real = {(float(r['L']), float(r['beta'])): r for r in rows}
    out = []
    print("\nN2 -- %d random-frequency surrogates per L (log p +- %.1f)" % (n_sur, spread))
    print("%-6s %5s %12s %22s %8s %12s %22s %8s" %
          ("L", "beta", "t0 real", "t0 surrogate 95%", "in band", "meas real",
           "meas surrogate 95%", "in band"))
    for L in Ls:
        s = Symbol(L)
        for beta in betas:
            rng = np.random.default_rng(seed + int(1000 * L))
            t0s, ms, ns = [], [], []
            t_start = time.time()
            for _ in range(n_sur):
                sg = surrogate(s, rng, spread)
                r = scan(sg, beta)
                t0s.append(r['t0_upper']); ms.append(r['dip_measure']); ns.append(r['n_dips'])
            t0s = np.array(t0s); ms = np.array(ms); ns = np.array(ns)
            rr = real[(L, beta)]
            t0r = float(rr['t0_upper']); mr = float(rr['dip_measure']); nr = int(rr['n_dips'])
            lo_t, hi_t = np.percentile(t0s, [2.5, 97.5])
            lo_m, hi_m = np.percentile(ms, [2.5, 97.5])
            lo_n, hi_n = np.percentile(ns, [2.5, 97.5])
            in_t = lo_t <= t0r <= hi_t; in_m = lo_m <= mr <= hi_m; in_n = lo_n <= nr <= hi_n
            print("%-6s %5.2f %12.5g [%9.5g,%9.5g] %8s %12.3f [%9.3f,%9.3f] %8s" %
                  (L, beta, t0r, lo_t, hi_t, in_t, mr, lo_m, hi_m, in_m))
            out.append(dict(L=L, beta=beta, n_surrogates=n_sur, spread=spread, seed=seed,
                            t0_real=t0r, t0_lo=lo_t, t0_hi=hi_t, t0_in_band=bool(in_t),
                            t0_median=float(np.median(t0s)),
                            measure_real=mr, measure_lo=lo_m, measure_hi=hi_m,
                            measure_in_band=bool(in_m),
                            ndip_real=nr, ndip_lo=lo_n, ndip_hi=hi_n,
                            ndip_in_band=bool(in_n), runtime_s=time.time() - t_start))
    with open(os.path.join(ROOT, 'data', 'stage3_n2.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
    return out


if __name__ == '__main__':
    n1_table()
    n2_table(n_sur=int(sys.argv[1]) if len(sys.argv) > 1 else 200,
             betas=(0.0, 0.25, 0.5), Ls=[L for L in STAGE2_L if L not in (0.6, 0.8)])
