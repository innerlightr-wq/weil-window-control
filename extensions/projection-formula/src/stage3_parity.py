#!/usr/bin/env python3
"""Stage 3 (H3): is the odd-sector floor high *because* odd f has F(0)=0?

In the EVEN basis f_k = cos(k pi x / L) one has int f_k = 0 for every k >= 1 and 2L for
k = 0, so the constraint F(0) = int f = 0 is exactly "delete the constant mode" (T1).
We compare the constrained-even floor with the odd floor.
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'src'))
import mpmath as mp
from zeta_window import build_matrix, F, norm2, EVEN, ODD
from weilform import eigsym

CELLS = [('0.5', 14, 50), ('0.6', 14, 60), ('0.8', 16, 70), ('1.0', 18, 80)]

rows = []
print('  %-5s %-4s %-16s %-16s %-16s %-11s %-11s' %
      ('L', 'N', 'even floor', 'constrained even', 'odd floor', 'constr/odd', 'even/odd'))
for Lx, N, dps in CELLS:
    mp.mp.dps = dps
    L = mp.mpf(Lx)
    Me = build_matrix(L, N, dps, EVEN)
    Mo = build_matrix(L, N, dps, ODD)
    # verify the claim int f_k = 0 for k>=1 (T1 check, numerical confirmation)
    nrm = [mp.sqrt(norm2(k, L, EVEN)) for k in range(N)]
    u = [F(k, L, mp.mpf(0), EVEN)/nrm[k] for k in range(N)]
    off = max(abs(u[k]) for k in range(1, N))
    # constrained even = delete the constant mode
    Mc = mp.matrix(N-1, N-1)
    for i in range(1, N):
        for j in range(1, N):
            Mc[i-1, j-1] = Me[i, j]
    le = eigsym(Me)[0][0]; lc = eigsym(Mc)[0][0]; lo = eigsym(Mo)[0][0]
    rows.append(dict(L=Lx, N=N, dps=dps, even=mp.nstr(le, 10), constrained=mp.nstr(lc, 10),
                     odd=mp.nstr(lo, 10), u_offdiag_max=mp.nstr(off, 4),
                     ratio_constr_odd=float(lc/lo), ratio_even_odd=float(le/lo)))
    print('  %-5s %-4d %-16s %-16s %-16s %-11.4g %-11.4g' %
          (Lx, N, mp.nstr(le, 8), mp.nstr(lc, 8), mp.nstr(lo, 8), float(lc/lo), float(le/lo)))
json.dump(rows, open('data/stage3_parity.json', 'w'), indent=1)
print()
print('  max |F_k(0)| over k>=1 (should be 0; confirms the constraint = drop mode 0):',
      ', '.join(r['u_offdiag_max'] for r in rows))
ok = all(1/3 <= r['ratio_constr_odd'] <= 3 for r in rows)
print('  K4: constrained-even within a factor 3 of odd at every L ?  %s' % ('YES' if ok else 'NO'))
for r in rows:
    print('     L=%-4s constrained/odd = %.4g %s' % (r['L'], r['ratio_constr_odd'],
          '' if 1/3 <= r['ratio_constr_odd'] <= 3 else '  <-- outside factor 3'))
