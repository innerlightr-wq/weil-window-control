"""(a) Which PART of the Weil form carries the old/new coupling?
   (b) Does the delta^2 coefficient C(L) = 4 F'(gamma)^2 decompose additively?"""
import sys, time
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)                                            # this note's wform.py
sys.path.insert(0, os.path.join(_HERE, '..', '..', '..', 'src'))     # the project's modules
_DATA = os.path.join(_HERE, '..', 'data')
import mpmath as mp, zeta_window as zw, wform as W, weilform as wf

GAMMA = mp.mpf('14.134725141734693790457251983562470270784257115699')

def parts_entry(fp, gp, U, prim):
    C0 = W.csym(fp, gp, mp.mpf(0))
    pole = 2*W.p_half(fp)*W.p_half(gp)
    archv = W.arch(fp, gp, U)
    logpi = mp.log(mp.pi)*C0
    pr = sum(2*lam/mp.sqrt(mp.mpf(n))*W.csym(fp,gp,logn) for n,lam,logn in prim)
    return pole, archv, logpi, pr

def Fprime(fp):
    """F'(gamma) = -int x f(x) sin(gamma x) dx."""
    tot = mp.mpf(0)
    for (x0,a,p,q) in fp:
        tot += mp.quad(lambda x: x*mp.cos(a*(x-x0))*mp.sin(GAMMA*x), [p,q])
    return -tot

def run(L, Lp, Nin, Nnew, dps, kind='even'):
    mp.mp.dps = dps+12
    L, Lp = mp.mpf(L), mp.mpf(Lp)
    if kind=='even': bi,bn = W.inner_basis(L,Nin,'even'), W.annulus_basis(L,Lp,Nnew,'even')
    else:            bi,bn = W.inner_basis(L,Nin,'even_d'), W.sine_annulus_basis(L,Lp,Nnew)
    b = bi+bn; n = len(b)
    nrm = [mp.sqrt(W.l2(f,f)) for f in b]
    prim = zw.von_mangoldt_terms(2*Lp)
    P = [mp.zeros(n,n) for _ in range(4)]
    for i in range(n):
        for j in range(i,n):
            vals = parts_entry(b[i],b[j],2*Lp,prim)
            for k,v in enumerate(vals):
                P[k][i,j]=P[k][j,i]=v/(nrm[i]*nrm[j])
    M = P[0]+P[1]-P[2]-P[3]
    old, new = range(Nin), range(Nin,n)
    names = ['Pole','Arch','LogPi','Prime']
    print("\n  L=%s -> L'=%s   Nin=%d Nnew=%d dps=%d basis=%s" % (L,Lp,Nin,Nnew,dps,kind))
    print("    %-8s %16s %16s %14s" % ("part","||part||_F","||cross||_F","cross/part"))
    for k,nm in enumerate(names):
        a_, c_ = W.fro(P[k]), W.fro(P[k],old,new)
        print("    %-8s %16s %16s %14s" % (nm, mp.nstr(a_,6), mp.nstr(c_,6),
              mp.nstr(c_/a_,4) if a_>0 else "-"))
    print("    %-8s %16s %16s %14s" % ("TOTAL", mp.nstr(W.fro(M),6), mp.nstr(W.fro(M,old,new),6),
          mp.nstr(W.fro(M,old,new)/W.fro(M),4)))

    # ---- delta^2 additivity ----
    mp.mp.dps = dps
    g = [Fprime(f)/nr for f,nr in zip(b,nrm)]
    vals, vecs = wf.eigsym(M)
    v = vecs[0]
    SI = sum(v[i]*g[i] for i in old); SN = sum(v[j]*g[j] for j in new)
    wI = mp.sqrt(sum(v[i]**2 for i in old)); wN = mp.sqrt(sum(v[j]**2 for j in new))
    # ground state of Q_L alone, in the same inner basis
    mp.mp.dps = dps+12
    QL = W.scale_matrix(W.Q_matrix(bi,2*L,zw.von_mangoldt_terms(2*L)), nrm[:Nin])
    mp.mp.dps = dps
    vL_, vecL = wf.eigsym(QL); u = vecL[0]
    SL = sum(u[i]*g[i] for i in old)
    # angle between the inner part of the L' ground state and the L ground state
    ip = abs(sum(v[i]*u[i] for i in old))/wI if wI>0 else mp.mpf(0)
    C_Lp, C_L = 4*(SI+SN)**2, 4*SL**2
    print("    ground state: ||v_inner||=%s  ||v_new||=%s   |<v_in/|v_in|, gs(Q_L)>|=%s"
          % (mp.nstr(wI,6), mp.nstr(wN,6), mp.nstr(ip,6)))
    print("    C(L')=%s   C(L)=%s   4S_I^2=%s  4S_N^2=%s  cross 8 S_I S_N=%s"
          % (mp.nstr(C_Lp,6), mp.nstr(C_L,6), mp.nstr(4*SI**2,6), mp.nstr(4*SN**2,6),
             mp.nstr(8*SI*SN,6)))
    print("    additive?  C(L)+4S_N^2 = %s   vs C(L') = %s   rel err %s"
          % (mp.nstr(C_L+4*SN**2,6), mp.nstr(C_Lp,6),
             mp.nstr(abs(C_L+4*SN**2-C_Lp)/abs(C_Lp),4)))
    print("    CS bounds: 4L^3/3 = %s -> 4L'^3/3 = %s   (Taylor incr 4L^2 dL = %s)"
          % (mp.nstr(4*L**3/3,6), mp.nstr(4*Lp**3/3,6), mp.nstr(4*L**2*(Lp-L),6)))

for c in [(0.5,0.8,6,6,30,'even'), (0.8,1.0,6,6,30,'even'), (0.5,0.8,6,6,30,'even_d')]:
    t=time.time(); run(*c); print("    (%.0fs)" % (time.time()-t))
