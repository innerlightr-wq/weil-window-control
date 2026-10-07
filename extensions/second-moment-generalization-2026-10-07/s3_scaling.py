"""Section 3: what is geometric and what is normalization."""
import sys; sys.path.insert(0,'.')
from kernel2 import *
getcontext().prec = 60
def dd(x): return D(str(x))

def coeff(xs, ws):            # |quadratic coefficient| = 2 V_mu
    return 2*centred_moment(xs, ws)
def midpoint(L, N):           # midpoint rule on [-L,L]: w = 2L/N
    h = 2*L/N
    return [-L + h*(D(i)+D('0.5')) for i in range(N)], [h]*N

print("="*84); print("3.1  five constructions, all with the SAME kernel 2cosh(a(x-y))"); print("="*84)
print(f"  {'construction':<42} {'M':>12} {'V/M':>12} {'V_mu':>14} {'coeff=2V':>14}")
rows=[]
for L in (1,2,4,8):
    Lv=D(L)
    # (i) Lebesgue on [-L,L]  (mass grows like L^1)
    xs,ws = midpoint(Lv, 400); rows.append((f"(i) Lebesgue [-L,L], L={L}", xs, ws))
    # (ii) probability-normalised on [-L,L]  (mass 1)
    xs2 = xs; ws2 = [w/(2*Lv) for w in ws]; rows.append((f"(ii) probability on [-L,L], L={L}", xs2, ws2))
    # (iii) fixed 5 nodes stretched by L
    xs3 = [Lv*D(t)/4 for t in (-4,-1,0,2,4)]; ws3=[ONE]*5; rows.append((f"(iii) 5 fixed nodes x L, L={L}", xs3, ws3))
    # (v) fixed profile (triangle) rescaled to support [-L,L], mass fixed = 1
    xs5,ws5 = midpoint(Lv,400)
    prof=[(ONE-abs(x)/Lv) for x in xs5]; s=sum(p*w for p,w in zip(prof,ws5))
    ws5=[p*w/s for p,w in zip(prof,ws5)]; rows.append((f"(v) triangle profile, mass 1, L={L}", xs5, ws5))
for name,xs,ws in rows:
    M=mass(xs,ws); V=centred_moment(xs,ws)
    print(f"  {name:<42} {str(M)[:12]:>12} {str(V/M)[:12]:>12} {str(V)[:14]:>14} {str(coeff(xs,ws))[:14]:>14}")

print()
print("  (iv) increasing nodes at FIXED spacing 1 (counting measure on {0..R}):")
for R in (4,8,16,32,64):
    xs=[D(i) for i in range(R+1)]; ws=[ONE]*(R+1)
    print(f"      R={R:<4} M={R+1:<5} V_mu={str(centred_moment(xs,ws))[:12]:<12} coeff={str(coeff(xs,ws))[:12]:<12} coeff/R^3={float(coeff(xs,ws)/D(R)**3):.6f}")

print()
print("="*84); print("3.2  the scaling law  coeff ~ L^(s+2)  when mass ~ L^s and V/M ~ L^2"); print("="*84)
import math
def slope(pairs):   # log-log slope of the last two points
    (x1,y1),(x2,y2)=pairs[-2],pairs[-1]
    return math.log(float(y2)/float(y1))/math.log(float(x2)/float(x1))
for label, build, s_expected in [
    ("Lebesgue  (mass ~ L^1)",      lambda L: midpoint(L,400),                                   1),
    ("probability (mass ~ L^0)",    lambda L: (midpoint(L,400)[0],[w/(2*L) for w in midpoint(L,400)[1]]), 0),
    ("5 fixed nodes (mass ~ L^0)",  lambda L: ([L*D(t)/4 for t in (-4,-1,0,2,4)],[ONE]*5),        0),
    ("Lebesgue x L^2 (mass ~ L^3)", lambda L: (midpoint(L,400)[0],[w*L*L for w in midpoint(L,400)[1]]), 3),
]:
    pts=[]
    for L in (2,4,8,16):
        xs,ws=build(D(L)); pts.append((D(L), coeff(xs,ws)))
    print(f"  {label:<30} measured log-log slope = {slope(pts):8.5f}   predicted s+2 = {s_expected+2}")

print()
print("="*84); print("3.3  measure scaling vs inner-product scaling (NOT the same thing)"); print("="*84)
xs,ws = midpoint(D(1),200); c=D(7)
l1 = lambda_pm_closed(xs,ws,D('0.01'))[1]
l2 = lambda_pm_closed(xs,[c*w for w in ws],D('0.01'))[1]
print(f"  mu -> 7mu : lambda_- {str(l1)[:20]} -> {str(l2)[:20]};  ratio = {float(l2/l1):.10f}  (expect 7)")
print("  Scaling the MEASURE scales both the map and the inner product: lambda -> c lambda,")
print("  and V_mu -> c V_mu, so coeff = 2V_mu is consistent. Multiplying only the inner")
print("  product by c, with the linear map fixed, instead gives lambda -> lambda/c. Different.")

print()
print("="*84); print("3.4  N-convergence at FIXED L does not create L-scaling"); print("="*84)
exact = lambda L,a: 2*L - dsinh(2*a*L)/a
for L in ('1','4'):
    Lv=D(L); a=D('0.001')
    print(f"  L={L}:  exact 2L-sinh(2aL)/a = {str(exact(Lv,a))[:22]}")
    for N in (25,50,100,200,400):
        xs,ws=midpoint(Lv,N)
        print(f"     N={N:<4} quadrature lambda_- = {str(lambda_pm_closed(xs,ws,a)[1])[:22]}")
print("  => N refines the quadrature at fixed L; the L-power comes from the measure, not N.")
