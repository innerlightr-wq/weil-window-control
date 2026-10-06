#!/usr/bin/env python3
"""Stage 3e, adversarial extension: does the DETECTION WINDOW depend on the HEIGHT
of the planted off-line zero, not just its depth?  (It does -- strongly.)"""
import sys, os, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from stage3e_teeth import zeros_matrix
from weilform import eigsym

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
mp.mp.dps = 70
K = 300
g = [mp.im(mp.zetazero(n)) for n in range(1, K + 1)]
Ls = ['0.4', '0.5', '0.6', '0.8', '1.0', '1.3', '1.6', '2.0', '2.5', '3.0']
heights = [('3.0', 'below gamma_1'), ('14.134725141734693', '= gamma_1'),
           ('18.0', 'between g1 and g2'), ('50.0', 'near gamma_10'),
           ('100.0', 'near gamma_29'), ('200.0', 'near gamma_79')]
delta = mp.mpf('0.1')
rows = []
print(f"planted depth delta = {delta} (Re rho = 0.6), K={K} true zeros, dps={mp.mp.dps}")
print(f"{'L':>6s} {'control':>14s}" + ''.join(f"{h[:9]:>13s}" for h, _ in heights))
first = {h: None for h, _ in heights}
for Lx in Ls:
    L = mp.mpf(Lx)
    N = max(12, int(mp.ceil(10 * L)))
    c = eigsym(zeros_matrix(L, N, g, []))[0][0]
    line = f"{Lx:>6s} {mp.nstr(c,5):>14s}"
    rec = dict(L=Lx, N=N, K=K, dps=mp.mp.dps, delta=str(delta),
               lam_min_control=mp.nstr(c, 10))
    for h, lab in heights:
        p = eigsym(zeros_matrix(L, N, g, [(mp.mpf(h), delta)]))[0][0]
        neg = p < 0 and abs(p) > 100 * (abs(c) if c < 0 else mp.mpf(10) ** (-mp.mp.dps + 10))
        if neg and first[h] is None:
            first[h] = Lx
        line += f"{('NEG ' if neg else '  + ') + mp.nstr(abs(p),3):>13s}"
        rec[f'lam_min_h{h}'] = mp.nstr(p, 10)
        rec[f'sign_h{h}'] = 'neg' if neg else 'pos'
    print(line)
    rows.append(rec)
print("\nfirst L at which lambda_min < 0, by planted height:")
for h, lab in heights:
    print(f"   gamma_* = {h:>20s} ({lab:<20s}) -> L = {first[h]}")
with open(os.path.join(ROOT, 'data', 'stage3e_height_sweep.csv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print("wrote data/stage3e_height_sweep.csv")
