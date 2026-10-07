"""Convergence of the coupling in basis size and precision, with a basis-free metric.

Frobenius ratio dilutes when you add weakly-coupled high modes, so also report
  sigma_cross = sup{ Q(f,g) : f in H_L, g in H_new, ||f||=||g||=1 }  (top singular value
of the cross block -- invariant under orthonormal changes of basis inside each block),
normalised by ||Q||_op.  If H_{L'} = H_L (+) N were a reducing decomposition this is 0."""
import sys, time, json
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)                                            # this note's wform.py
sys.path.insert(0, os.path.join(_HERE, '..', '..', '..', 'src'))     # the project's modules
_DATA = os.path.join(_HERE, '..', 'data')
import mpmath as mp, zeta_window as zw, wform as W, weilform as wf

def svmax(M, rows, cols):
    """top singular value of the sub-block, via sqrt(lambda_max(B^T B))."""
    B = mp.matrix([[M[i,j] for j in cols] for i in rows])
    G = B.T*B
    v,_ = wf.eigsym(G)
    return mp.sqrt(max(v[-1], mp.mpf(0)))

def run(L,Lp,Nin,Nnew,dps,kind):
    mp.mp.dps = dps+12
    L,Lp = mp.mpf(L), mp.mpf(Lp)
    if kind=='even': bi,bn = W.inner_basis(L,Nin,'even'), W.annulus_basis(L,Lp,Nnew,'even')
    else:            bi,bn = W.inner_basis(L,Nin,'even_d'), W.sine_annulus_basis(L,Lp,Nnew)
    b=bi+bn; n=len(b); nrm=[mp.sqrt(W.l2(f,f)) for f in b]
    M = W.scale_matrix(W.Q_matrix(b,2*Lp,zw.von_mangoldt_terms(2*Lp)), nrm)
    mp.mp.dps = dps
    old,new = range(Nin), range(Nin,n)
    ev,_ = wf.eigsym(M)
    opn = max(abs(ev[0]),abs(ev[-1]))
    fro_c = W.fro(M,old,new)/W.fro(M)
    sig = svmax(M,old,new)/opn
    mx  = max(abs(M[i,j]) for i in old for j in new)/opn
    return fro_c, sig, mx

print("  %-5s %-5s %-5s %-5s %-7s %12s %14s %14s %8s" %
      ("L","L'","Nin","Nnew","dps","E_cross(F)","sigma_x/||Q||","max|cross|","t"))
rows=[]
CASES=[(0.5,0.8,n,n,30,'even') for n in (4,6,8,12,16)] + \
      [(0.5,0.8,8,8,d,'even') for d in (20,30,50)] + \
      [(0.5,0.8,8,8,30,'even_d'),(0.8,1.6,8,8,30,'even'),(1.0,1.2,8,8,30,'even')]
for c in CASES:
    t=time.time(); f_,s_,m_ = run(*c)
    print("  %-5s %-5s %-5d %-5d %-7s %12s %14s %14s %7.0fs" %
          (c[0],c[1],c[2],c[3],"%d/%s"%(c[4],c[5]), mp.nstr(f_,4), mp.nstr(s_,4),
           mp.nstr(m_,4), time.time()-t))
    rows.append(dict(L=c[0],Lp=c[1],N=c[2],dps=c[4],kind=c[5],
                     fro=mp.nstr(f_,6),sig=mp.nstr(s_,6),mx=mp.nstr(m_,6)))
json.dump(rows,open(os.path.join(_DATA,'conv.json'),'w'),indent=1)
