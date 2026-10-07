#!/usr/bin/env python3
"""comb-dips quick check: the Stage 1 gate, then recompute Stage 2 rows and
compare them against the committed data/stage2_t0.csv.  Exits non-zero on failure."""
import sys, os, csv, time, subprocess
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.dirname(HERE)
import numpy as np
from symbol import Symbol
from scan import scan

# rows to recompute: (L, beta).  L=1.19 is the Zhu-anchored case, L=1.6 the headline.
ROWS = [(1.19, 0.0), (1.6, 0.0)]
TOL_T0 = 1e-6          # absolute agreement of the certified bracket, in t
TOL_MEAS = 1e-6

def main():
    print("=" * 74)
    print("comb-dips quick check")
    print("=" * 74)

    print("\n[1/2] Stage 1 gate (A_L table vs Zhu; Zhu section 7 at L=1.19)\n")
    r = subprocess.run([sys.executable, os.path.join(HERE, 'stage1_validate.py')],
                       capture_output=True, text=True)
    tail = [l for l in r.stdout.splitlines() if 'GATE' in l or 'precision :' in l]
    print("\n".join("  " + l for l in tail))
    gate_ok = r.returncode == 0 and all('PASS' in l for l in tail) and len(tail) == 3
    if not gate_ok:
        print(r.stdout[-2000:]); print(r.stderr[-2000:])

    print("\n[2/2] Stage 2: recompute %d rows and compare with data/stage2_t0.csv\n"
          % len(ROWS))
    ref = {(float(x['L']), float(x['beta'])): x
           for x in csv.DictReader(open(os.path.join(ROOT, 'data', 'stage2_t0.csv')))}
    print("  %-5s %-6s %14s %14s %12s %10s %s"
          % ("L", "beta", "t0 (recomputed)", "t0 (committed)", "|diff|", "t0/T1", ""))
    scan_ok = True
    for L, beta in ROWS:
        sym = Symbol(L)
        t = time.time()
        s = scan(sym, beta)
        dt = time.time() - t
        want = ref[(L, beta)]
        d_t0 = abs(s['t0_lower'] - float(want['t0_lower']))
        d_m = abs(s['dip_measure'] - float(want['dip_measure']))
        ok = d_t0 < TOL_T0 and d_m < TOL_MEAS and s['n_dips'] == int(want['n_dips'])
        scan_ok &= ok
        print("  %-5s %-6s %14.4f %14.4f %12.2e %10.4f %s  (%.1fs)"
              % (L, beta, s['t0_lower'], float(want['t0_lower']), d_t0,
                 s['t0_lower'] / sym.T1, "OK" if ok else "MISMATCH", dt))
        if not ok:
            print("      measure %.9f vs %.9f (|diff| %.2e); n_dips %d vs %s"
                  % (s['dip_measure'], float(want['dip_measure']), d_m,
                     s['n_dips'], want['n_dips']))

    print("\n" + "=" * 74)
    print("  Stage 1 gate     : %s" % ("PASS" if gate_ok else "FAIL"))
    print("  Stage 2 agreement: %s" % ("PASS" if scan_ok else "FAIL"))
    print("=" * 74)
    return 0 if (gate_ok and scan_ok) else 1


if __name__ == '__main__':
    sys.exit(main())
