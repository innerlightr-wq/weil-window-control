#!/usr/bin/env python3
"""Stage 3 REVALIDATION for the revision.

(1) Convergence of the three floors in N (the original run used a single N per L).
(2) Comparison against Zhu arXiv:2608.24827 Table 1 (converged exact-assembly
    REFERENCE values, author's own words: "reference values, not certificates").
(3) The point-vs-interval test.  Zhu's stated mechanism (sec. 6, verbatim) is:
    "the asymmetry must instead reflect the constraint F_o(0) = 0, which denies odd
     test functions the zero-free frequency interval [0, gamma_1) in which even
     minimizers park most of their mass."
    Our control imposes the single linear condition F(0) = 0 on the even sector.
    That is strictly weaker than oddness: an even F with F(0)=0 has one zero at the
    origin but may still concentrate on (0, gamma_1), whereas an odd F vanishes at 0
    AND is antisymmetric, suppressing |F|^2 throughout a neighbourhood of 0.
    We therefore measure, for each minimiser, the mass fraction
        m = int_0^{gamma_1} |F|^2 dt  /  int_0^infty |F|^2 dt,  denominator = pi
    (Plancherel: F_k(t) = int f_k e^{itx} dx, ||f||_2 = 1, |F|^2 even => int_0^inf = pi).
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'src'))
import mpmath as mp
from zeta_window import build_matrix, F, norm2, EVEN, ODD
from weilform import eigsym

# Zhu arXiv:2608.24827, Table 1 (reference values, not certificates)
ZHU = {'0.5': ('9.34e-7', '1.94e-4', 208), '0.6': ('1.61e-9', '5.97e-7', 371),
       '0.7': ('4.18e-13', '2.56e-10', 612), '0.8': ('1.65e-17', '1.57e-14', 949),
       '0.9': ('4.14e-23', '7.22e-20', 1746), '1.0': ('5.88e-30', '1.49e-26', 2535)}
CELLS = [('0.5', 50, (14, 20, 26)), ('0.6', 60, (14, 20, 26)),
         ('0.8', 90, (16, 22, 28)), ('1.0', 120, (18, 24, 30))]

def floors(L, N, dps):
    Me = build_matrix(L, N, dps, EVEN); Mo = build_matrix(L, N, dps, ODD)
    Mc = mp.matrix(N-1, N-1)
    for i in range(1, N):
        for j in range(1, N): Mc[i-1, j-1] = Me[i, j]
    (ve, Ve), (vc, Vc), (vo, Vo) = eigsym(Me), eigsym(Mc), eigsym(Mo)
    return (ve[0], mp.matrix(Ve[0])), (vc[0], mp.matrix(Vc[0])), (vo[0], mp.matrix(Vo[0]))

print('(1)+(2) convergence in N, and comparison with Zhu Table 1\n')
print('  %-5s %-4s %-14s %-14s %-14s %-9s %-11s %-11s' %
      ('L','N','even floor','constr. even','odd floor','constr/odd','even/Zhu','odd/Zhu'))
conv = []
for Lx, dps, Ns in CELLS:
    mp.mp.dps = dps; L = mp.mpf(Lx)
    for N in Ns:
        (le,_), (lc,_), (lo,_) = floors(L, N, dps)
        ze, zo, _ = ZHU[Lx]
        print('  %-5s %-4d %-14s %-14s %-14s %-9.4g %-11.4g %-11.4g' %
              (Lx, N, mp.nstr(le,7), mp.nstr(lc,7), mp.nstr(lo,7),
               float(lc/lo), float(le/mp.mpf(ze)), float(lo/mp.mpf(zo))))
        conv.append(dict(L=Lx, N=N, dps=dps, even=mp.nstr(le,10), constr=mp.nstr(lc,10),
                         odd=mp.nstr(lo,10), ratio_constr_odd=float(lc/lo),
                         even_over_zhu=float(le/mp.mpf(ze)), odd_over_zhu=float(lo/mp.mpf(zo)),
                         ratio_even_odd=float(le/lo), zhu_ratio=ZHU[Lx][2]))
    print()

print('(3) point-vs-interval: mass fraction of |F|^2 on [0, gamma_1)\n')
print('  %-5s %-4s %-9s %-22s %-22s %-22s' %
      ('L','N','g1','even minimiser','constr-even minimiser','odd minimiser'))
mass = []
for Lx, dps, Ns in CELLS:
    mp.mp.dps = dps; L = mp.mpf(Lx); N = Ns[-1]
    g1 = mp.im(mp.zetazero(1))
    (le,ve), (lc,vc), (lo,vo) = floors(L, N, dps)
    nrm_e = [mp.sqrt(norm2(k, L, EVEN)) for k in range(N)]
    nrm_o = [mp.sqrt(norm2(k, L, ODD)) for k in range(N)]
    def Ffun(c, sector, nrm, drop0=False):
        off = 1 if drop0 else 0
        return lambda t: sum(c[k-off]*F(k, L, t, sector)/nrm[k] for k in range(off, N))
    out = []
    for tag, c, sector, nrm, drop0 in (('even', ve, EVEN, nrm_e, False),
                                       ('constr', vc, EVEN, nrm_e, True),
                                       ('odd', vo, ODD, nrm_o, False)):
        G = Ffun(c, sector, nrm, drop0)
        num = mp.quad(lambda t: G(t)**2, [0, g1])
        frac = num/mp.pi
        out.append(frac)
        mass.append(dict(L=Lx, N=N, which=tag, mass_0_g1=float(frac)))
    print('  %-5s %-4d %-9s %-22.10f %-22.10f %-22.10f' %
          (Lx, N, mp.nstr(g1,6), float(out[0]), float(out[1]), float(out[2])))
json.dump(dict(convergence=conv, mass=mass), open('data/stage3_revalidate.json','w'), indent=1)
