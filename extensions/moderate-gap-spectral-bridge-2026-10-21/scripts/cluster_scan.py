#!/usr/bin/env python3
"""How much budget does a CLUSTER of g off-line pairs cost?

sum_j (q_j-2)^2  depends only on the Gram of {g_l, h_l} for the off-line pairs
(P and the on-line zeros do not enter Q), so this is a self-contained computation.
Kernel: eta^2 = 1 on (-1/2,1/2)  (K = sinc; structural experiment, not an
admissible ARTZ window -- see 05_GAP_CONDITIONING.md)."""
import math, sys
sys.path.insert(0, __file__.rsplit('/',1)[0])
import linalg as la

def ReK(xi, et):
    r = xi*xi + et*et
    if r < 1e-24: return 1.0
    return (xi*math.sin(math.pi*xi)*math.cosh(math.pi*et)
            + et*math.cos(math.pi*xi)*math.sinh(math.pi*et))/(math.pi*r)
def ImK(xi, et):
    r = xi*xi + et*et
    if r < 1e-24: return 0.0
    return -(et*math.sin(math.pi*xi)*math.cosh(math.pi*et)
             - xi*math.cos(math.pi*xi)*math.sinh(math.pi*et))/(math.pi*r)
def kap(y): return 1.0 if abs(y)<1e-18 else math.sinh(2*math.pi*y)/(2*math.pi*y)

def slack(offsets, y):
    g = len(offsets)
    B = [('g',x,y) for x in offsets] + [('h',x,y) for x in offsets]
    def ge(a,b):
        D = a[1]-b[1]; ya, yb = a[2], b[2]
        if a[0]=='g' and b[0]=='g': return 0.5*(ReK(D,ya+yb)+ReK(D,ya-yb))
        if a[0]=='h' and b[0]=='h': return 0.5*(ReK(D,ya+yb)-ReK(D,ya-yb))
        if a[0]=='g' and b[0]=='h': return -0.5*(ImK(D,ya+yb)-ImK(D,ya-yb))
        return -0.5*(ImK(-D,ya+yb)-ImK(-D,yb-ya))
    n = 2*g
    M = [[ge(B[i],B[j]) for j in range(n)] for i in range(n)]
    D = [2.0]*g + [-2.0]*g
    q = la.spec_of_DM(D, M)
    qp = [v for v in q if v > 1e-9]
    return sum((v-2.0)**2 for v in qp), len(qp), sum(v*v for v in qp)

y = 1.0; ideal = (kap(y)-1)**2
print("eta^2=1 window, y=%.1f, kappa=%.4f, one-pair slack (kappa-1)^2 = %.6g\n" % (y, kap(y), ideal))
print("%-6s | %s" % ("sep", "  ".join("g=%-2d total/ideal" % g for g in (1,2,4,8,16,32))))
best = (1e9, None)
for sep in (0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.5,2.0,5.0,20.0):
    row = []
    for g in (1,2,4,8,16,32):
        s,_,_ = slack([i*sep for i in range(g)], y)
        row.append(s/ideal)
        if g >= 8 and s/ideal/g < best[0]: best = (s/ideal/g, (sep,g))
    print("%-6g |  " % sep + "  ".join("%-14.4f" % v for v in row))
print("\nworst (smallest) slack-per-pair / ideal : %.5f  at sep=%g, g=%d" % (best[0], *best[1]))
# does the cluster TOTAL grow with g at the worst separation?
sep = best[1][0]
print("\ncluster total at sep=%g as g grows (in units of (kappa-1)^2):" % sep)
for g in (1,2,4,8,16,32,64):
    s,npos,_ = slack([i*sep for i in range(g)], y)
    print("   g=%-3d  total=%-12.5f  per pair=%-12.6f  n_+=%d" % (g, s/ideal, s/ideal/g, npos))
