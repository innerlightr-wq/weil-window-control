#!/usr/bin/env python3
"""Is the discarded Lemma-5 slack ADDITIVE over off-line pairs?  And why not?

Admissible ARTZ window (tau=1/200).  g pairs of equal depth at normalised ordinate
separation sep (sep >= 1 means NO local excess over the mean zero density, so these
configurations violate no counting constraint)."""
import math, sys
sys.path.insert(0, __file__.rsplit('/',1)[0])
import linalg as la
exec(open(__file__.rsplit('/',1)[0] + '/cluster_smooth.py').read().split('for y in (0.5')[0])

def pos_spec(offsets, y):
    g = len(offsets)
    B = [('g',x,y) for x in offsets] + [('h',x,y) for x in offsets]
    C = {}
    def RK(D,e):
        k=(round(D,9),round(e,9))
        if k not in C: C[k]=ReK(D,e)
        return C[k]
    def IK(D,e):
        k=('i',round(D,9),round(e,9))
        if k not in C: C[k]=ImK(D,e)
        return C[k]
    def ge(a,b):
        D=a[1]-b[1]; ya,yb=a[2],b[2]
        if a[0]=='g' and b[0]=='g': return 0.5*(RK(D,ya+yb)+RK(D,ya-yb))
        if a[0]=='h' and b[0]=='h': return 0.5*(RK(D,ya+yb)-RK(D,ya-yb))
        if a[0]=='g' and b[0]=='h': return -0.5*(IK(D,ya+yb)-IK(D,ya-yb))
        return -0.5*(IK(-D,ya+yb)-IK(-D,yb-ya))
    n=2*g
    M=[[ge(B[i],B[j]) for j in range(n)] for i in range(n)]
    q=la.spec_of_DM([2.0]*g+[-2.0]*g, M)
    qp=sorted([v for v in q if v>1e-9], reverse=True)
    return qp, M

print("=== slack growth in the number g of off-line pairs (admissible window) ===")
for y, sep in ((1.0, 2.0), (1.0, 3.0), (0.5, 2.0)):
    kap = ReK(0.0, 2*y); ideal = (kap-1)**2
    a = 2*math.pi*y
    print("\n  a = delta logT = %.4f, kappa = %.5f, (kappa-1)^2 = %.6g, sep = %g"
          % (a, kap, ideal, sep))
    print("  %-5s %-12s %-12s %-12s %-10s" % ("g","total/ideal","per pair","q_max/ideal^.5","frac in q_max"))
    data=[]
    for g in (1,2,4,8,16,32):
        qp, M = pos_spec([i*sep for i in range(g)], y)
        tot = sum((v-2.0)**2 for v in qp)
        top = (qp[0]-2.0)**2
        data.append((g, tot/ideal))
        print("  %-5d %-12.4f %-12.5f %-12.4f %-10.4f"
              % (g, tot/ideal, tot/ideal/g, qp[0]/math.sqrt(ideal), top/tot))
    # fit total/ideal = A + B ln g
    import statistics
    xs=[math.log(g) for g,_ in data[1:]]; ys=[v for _,v in data[1:]]
    n=len(xs); sx=sum(xs); sy=sum(ys); sxx=sum(x*x for x in xs); sxy=sum(x*v for x,v in zip(xs,ys))
    B=(n*sxy-sx*sy)/(n*sxx-sx*sx); A=(sy-B*sx)/n
    pred=[A+B*x for x in xs]
    res=max(abs(p-v) for p,v in zip(pred,ys))
    print("  log fit: total/ideal = %.4f + %.4f ln g   (max residual %.4f over g=2..32)"
          % (A, B, res))
    # linear fit for comparison
    xs2=[g for g,_ in data]; ys2=[v for _,v in data]
    n=len(xs2); sx=sum(xs2); sy=sum(ys2); sxx=sum(x*x for x in xs2); sxy=sum(x*v for x,v in zip(xs2,ys2))
    B2=(n*sxy-sx*sy)/(n*sxx-sx*sx); A2=(sy-B2*sx)/n
    res2=max(abs(A2+B2*x-v) for x,v in zip(xs2,ys2))
    print("  lin fit: total/ideal = %.4f + %.4f g        (max residual %.4f)" % (A2,B2,res2))

print("\n=== mechanism: the g-directions ALIGN as the depth grows ===")
print("  corr(g_z, g_z') = [ReK(D+2iy)+ReK(D)]/(1+kappa)  at separation D")
print("  %-8s %-10s %-10s %-10s %-10s" % ("a=2pi y","kappa","D=1","D=2","D=5"))
for y in (0.1, 0.25, 0.5, 1.0, 1.5):
    kap = ReK(0.0, 2*y)
    row = []
    for D in (1.0, 2.0, 5.0):
        row.append((ReK(D, 2*y) + ReK(D, 0.0))/(1+kap))
    print("  %-8.4f %-10.4f %-10.5f %-10.5f %-10.5f" % (2*math.pi*y, kap, *row))
print("\n  as a grows, eta^2(u) cosh(2 pi y u) concentrates at |u| -> 1/2, so all g_z")
print("  (which differ only by the phase e^{-2pi i u xi} sampled there) become collinear.")
