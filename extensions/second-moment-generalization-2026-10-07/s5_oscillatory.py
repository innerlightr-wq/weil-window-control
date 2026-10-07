"""Section 5: the real conjugation-closed quartet, and what it says about K(L,gamma)."""
import sys, math; sys.path.insert(0,'.')
from kernel2 import *
getcontext().prec = 50

def quartet_matrix(xs, ws, a, g):
    """4 cos(g(x-y)) cosh(a(x-y)) = 4 Re cosh((a+ig)(x-y)); rank <= 4."""
    n=len(xs)
    return [[(ws[i]*ws[j]).sqrt()*4*dcosh(a*(xs[i]-xs[j]))*D(str(math.cos(g*float(xs[i]-xs[j]))))
             for j in range(n)] for i in range(n)]

def midpoint(L,N):
    h=2*L/N; return [-L+h*(D(i)+D('0.5')) for i in range(N)], [h]*N

def ip(u,v,ws): return sum(a*b*w for a,b,w in zip(u,v,ws))
def proj_perp(u, basis, ws):
    w=list(u)
    B=[]
    for b in basis:
        c=list(b)
        for e in B:
            t=ip(c,e,ws)/ip(e,e,ws); c=[ci-t*ei for ci,ei in zip(c,e)]
        if ip(c,c,ws)>0: B.append(c)
    for e in B:
        t=ip(w,e,ws)/ip(e,e,ws); w=[wi-t*ei for wi,ei in zip(w,e)]
    return w

print("="*86)
print("5.1  quartet kernel: predicted compression  P_E Q P_E = -4[ |P_E(x c)><.| + |P_E(x s)><.| ]")
print("="*86)
print("  4cos(g d)cosh(a d) = 4cos(g d) + 2a^2 d^2 cos(g d) + O(a^4),  d = x-y.")
print("  A_0 = 4cos(g(x-y)) = 4(|c><c| + |s><s|) >= 0, rank 2, c=cos(gx), s=sin(gx).")
print("  For v perp c,s the terms with a bare c(y) or s(y) die, leaving -4(|<xc,v>|^2+|<xs,v>|^2).")
print()
for L, g, N in [(D('1.0'), 3.0, 120), (D('1.6'), 5.0, 160)]:
    xs, ws = midpoint(L, N)
    c=[D(str(math.cos(g*float(x)))) for x in xs]
    s=[D(str(math.sin(g*float(x)))) for x in xs]
    xc_=[x*ci for x,ci in zip(xs,c)]; xs_=[x*si for x,si in zip(xs,s)]
    Pxc=proj_perp(xc_,[c,s],ws); Pxs=proj_perp(xs_,[c,s],ws)
    # 2x2 Gram of the two projected vectors -> largest eigenvalue
    g11,g12,g22 = ip(Pxc,Pxc,ws), ip(Pxc,Pxs,ws), ip(Pxs,Pxs,ws)
    tr, det = g11+g22, g11*g22-g12*g12
    lam_big = (tr + (tr*tr-4*det).sqrt())/2
    pred = -4*lam_big
    lam0 = jacobi_eigenvalues(quartet_matrix(xs,ws,ZERO,g))[0]
    print(f"  L={L}, gamma={g}, N={N}:  lambda_min(A_0) = {float(lam0):.3e}")
    print(f"     ||P_E(x cos)||^2={float(g11):.8f}  ||P_E(x sin)||^2={float(g22):.8f}  "
          f"cross={float(g12):.2e}  -> predicted coeff = {float(pred):.8f}")
    for aa in ('1e-3','1e-4'):
        av=D(aa); lam=jacobi_eigenvalues(quartet_matrix(xs,ws,av,g))[0]
        print(f"     a={aa}:  (lam-lam0)/a^2 = {float((lam-lam0)/av**2):.8f}")

print()
print("="*86); print("5.2  K(L,gamma) = L^3 k(gamma L):  the cubic is the LARGE-gamma-L regime"); print("="*86)
def K(L,g):   # exact closed form, bracketed as in the TeX
    return L**3/3 - ( L**2*math.sin(2*g*L)/(2*g) + L*math.cos(2*g*L)/(2*g**2) - math.sin(2*g*L)/(4*g**3) )
def k(z):     # K = L^3 k(z), z = gamma L ; k(z) = int_{-1}^1 t^2 sin^2(z t) dt
    if z == 0: return 0.0
    return 1/3 - ( math.sin(2*z)/(2*z) + math.cos(2*z)/(2*z**2) - math.sin(2*z)/(4*z**3) )
print(f"  {'z=gamma*L':>12} {'k(z)':>14} {'k(z)/z^2':>14} {'(2/5) z^2':>14}   regime")
for z in (0.01,0.1,0.5,1,2,5,10,50,200):
    kk=k(z)
    print(f"  {z:>12} {kk:>14.8f} {kk/z**2:>14.8f} {0.4*z*z:>14.8f}   "
          f"{'small: k ~ (2/5)z^2 -> K ~ (2/5) gamma^2 L^5' if z<=0.1 else ('large: k -> 1/3 -> K ~ L^3/3' if z>=10 else 'intermediate')}")
print(f"  k(0) = 0 exactly: at gamma=0 the quartet degenerates to a double pair and K vanishes.")
print(f"  consistency: K(L,g) vs L^3 k(gL) at L=0.8,g=14.1347 -> {K(0.8,14.1347):.10f} vs {0.8**3*k(14.1347*0.8):.10f}")

print()
print("="*86); print("5.3  does the projection reduce the coefficient below the unprojected norm?"); print("="*86)
print("  In the real-EVEN sector, x cos(gx) is odd and x sin(gx) is even, so only x sin(gx)")
print("  survives; its unprojected norm^2 IS K(L,gamma). But <x sin(gx), cos(gx)> != 0, so")
print("  the projection strictly reduces it:")
print(f"  {'L':>6} {'K(L,g)':>14} {'||P_E(x sin)||^2':>18} {'ratio':>10}")
for Lf in (0.8,1.0,1.3,1.6,2.0):
    L=D(str(Lf)); g=14.1347; N=400
    xs,ws=midpoint(L,N)
    cv=[D(str(math.cos(g*float(x)))) for x in xs]
    xsv=[x*D(str(math.sin(g*float(x)))) for x in xs]
    P=proj_perp(xsv,[cv],ws)                      # even sector: baseline range is span{cos}
    num=float(ip(P,P,ws)); den=K(Lf,g)
    print(f"  {Lf:>6} {den:>14.8f} {num:>18.8f} {num/den:>10.5f}")
print("  => an UNPROJECTED norm identity for F'(gamma) is an upper bound on the response,")
print("     not the response. This is exactly the (a)/(b)/(c) distinction of Remark 10.")
