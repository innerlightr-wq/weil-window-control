#!/usr/bin/env python3
"""Stage 2: certified scan for t0(L,beta) over all L and beta."""
import sys, os, json, time, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from symbol import Symbol, STAGE2_L
from scan import scan

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BETAS = [0.0, 0.25, 0.5]

if __name__ == '__main__':
    rows = []
    print("%-6s %-34s %8s %6s %11s %7s %9s %12s %10s %9s %8s" %
          ("L", "active n", "A_L", "beta", "t0", "t0/T1", "t0/T*", "dip measure",
           "n dips", "bracket", "sec"))
    for L in STAGE2_L:
        s = Symbol(L)
        for beta in BETAS:
            t0 = time.time()
            r = scan(s, beta)
            el = time.time() - t0
            br = r['t0_upper'] - r['t0_lower']
            print("%-6s %-34s %8.4f %6.2f %11.5g %7.4f %9.4f %12.5f %10d %9.1e %8.1f" %
                  (L, str(s.n.tolist())[:34], s.A, beta, r['t0_upper'],
                   r['t0_upper'] / s.T1, r['t0_upper'] / s.Tstar,
                   r['dip_measure'], r['n_dips'], br, el))
            rows.append(dict(L=L, active_n=str(s.n.tolist()), n_active=len(s.n),
                             A=s.A, T1=s.T1, Tstar=s.Tstar, beta=beta,
                             scan_hi=r['hi'],
                             t0_lower=r['t0_lower'], t0_upper=r['t0_upper'],
                             t0_bracket=br,
                             t0_over_T1=r['t0_upper'] / s.T1,
                             t0_over_Tstar=r['t0_upper'] / s.Tstar,
                             dip_measure=r['dip_measure'],
                             measure_lower=r['measure_lower'],
                             measure_upper=r['measure_upper'],
                             n_dips=r['n_dips'], runtime_s=el,
                             h_boundary=r['h_boundary']))
            np.save(os.path.join(ROOT, 'data',
                                 'dips_L%s_b%s.npy' % (L, beta)), np.array(r['dips']))
    with open(os.path.join(ROOT, 'data', 'stage2_t0.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print("\nwrote data/stage2_t0.csv and per-(L,beta) dip lists")
