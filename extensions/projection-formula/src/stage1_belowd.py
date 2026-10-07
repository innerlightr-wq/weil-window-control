#!/usr/bin/env python3
"""Stage 1, R < d = 2g: no kernel. Does a threshold near-null space work?

For R < 2g the honest T_R is positive definite, so the perturbation is NON-degenerate and
first-order PT gives  lam_min(T^(a)) = lam_min(T) + a^2 <v_min, Q_2 v_min> + O(a^4),
with v_min the ground state -- a rank-one projection, not a kernel projection.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'src'))
import mpmath as mp
from weilform import toeplitz, eigsym
from stage1_ff import CURVES, curve_data, thetas_of, tvals, dt_exact

mp.mp.dps = 60
print('  %-6s %-3s %-4s %-7s %-16s %-16s %-11s %-14s'%('curve','g','R','a','lam_min(T^a)','prediction','rel.err','lam_min(T)'))
rows=[]
for tag, f, q, g in CURVES:
    d = 2*g
    A, p_all = curve_data(f, q, g, d+14)
    ths, _ = thetas_of(A, q, g); th = ths[0]
    t0 = tvals(p_all, q, g, d+14)
    for R in range(1, d):
        T = toeplitz(t0, R)
        vals, vecs = eigsym(T)
        lam0 = vals[0]; v = mp.matrix(vecs[0])
        # Q_2(c) = 2 Re(S_2 conj(S_0)) - 2|S_1|^2   (the a^2 form, full, no kernel assumption)
        S = lambda k: sum((mp.mpf(i)**k)*v[i]*mp.e**(1j*i*th) for i in range(R+1))
        S0,S1,S2 = S(0),S(1),S(2)
        Q2 = 2*mp.re(S2*mp.conj(S0)) - 2*abs(S1)**2
        for a in ('1e-2','1e-3','1e-4'):
            aa=mp.mpf(a)
            ta=[t0[n]+dt_exact(n,th,aa) for n in range(R+1)]
            lam=eigsym(toeplitz(ta,R))[0][0]
            pred=lam0+aa**2*Q2
            rel=abs(lam-pred)/abs(pred)
            rows.append(dict(tag=tag,g=g,R=R,a=a,rel=float(rel),lam0=mp.nstr(lam0,10),
                             S0=float(abs(S0)),Q2=mp.nstr(Q2,10)))
            print('  %-6s %-3d %-4d %-7s %-16s %-16s %-11.3e %-14s'%(
                tag,g,R,a,mp.nstr(lam,12),mp.nstr(pred,12),float(rel),mp.nstr(lam0,8)))
json.dump(rows,open('data/stage1_belowd.json','w'),indent=1)
print()
print('  |S_0(v_min)| (0 would mean the ground state already lies in the vanishing locus):')
for r in rows:
    if r['a']=='1e-4': print('    %-6s R=%-3d |S_0|=%.4f  lam_min(T)=%s'%(r['tag'],r['R'],r['S0'],r['lam0']))
