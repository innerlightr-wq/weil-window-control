import sys, time
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)                                            # this note's wform.py
sys.path.insert(0, os.path.join(_HERE, '..', '..', '..', 'src'))     # the project's modules
_DATA = os.path.join(_HERE, '..', 'data')
import mpmath as mp
import zeta_window as zw
import wform as W

print("Validation: generic evaluator vs the project's build_matrix (EVEN sector)\n")
print("  %-5s %-4s %-5s %14s %14s %12s" % ("L","N","dps","||Q_proj||_F","max|diff|","rel"))
for (L, N, dps) in [(0.5,4,30),(0.8,5,30),(0.8,5,40),(1.0,4,40)]:
    mp.mp.dps = dps + 10
    Lm = mp.mpf(L)
    prim = zw.von_mangoldt_terms(2*Lm)
    b = W.inner_basis(Lm, N, 'even')
    nrm = [mp.sqrt(W.l2(f,f)) for f in b]
    t = time.time()
    M = W.scale_matrix(W.Q_matrix(b, 2*Lm, prim), nrm)
    mp.mp.dps = dps
    R = zw.build_matrix(L, N, dps, zw.EVEN)
    d = max(abs(M[i,j]-R[i,j]) for i in range(N) for j in range(N))
    nf = W.fro(R)
    print("  %-5s %-4d %-5d %14s %14s %12s   (%.0fs)" %
          (L, N, dps, mp.nstr(nf,6), mp.nstr(d,4), mp.nstr(d/nf,4), time.time()-t))
