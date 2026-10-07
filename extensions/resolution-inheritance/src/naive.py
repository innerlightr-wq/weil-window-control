"""Control: what you get if you compare same-INDEX blocks of the project's own
non-nested bases.  cos(k pi x/L) and cos(k pi x/L') span different subspaces, so any
agreement or disagreement here is a basis artifact, not a statement about the form."""
import sys
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)                                            # this note's wform.py
sys.path.insert(0, os.path.join(_HERE, '..', '..', '..', 'src'))     # the project's modules
_DATA = os.path.join(_HERE, '..', 'data')
import mpmath as mp, zeta_window as zw, wform as W
print("  %-5s %-5s %-4s %-5s %12s" % ("L","L'","N","dps","E_old(naive)"))
for (L,Lp,N,dps) in [(0.5,0.8,4,30),(0.5,0.8,6,30),(0.5,0.8,8,30),(0.8,1.0,6,30),(0.8,1.6,6,30),(1.0,1.2,6,30)]:
    A = zw.build_matrix(L, N, dps, zw.EVEN)
    B = zw.build_matrix(Lp, N, dps, zw.EVEN)
    d = W.fro(mp.matrix([[B[i,j]-A[i,j] for j in range(N)] for i in range(N)]))/W.fro(A)
    print("  %-5s %-5s %-4d %-5d %12s" % (L,Lp,N,dps,mp.nstr(d,4)))
