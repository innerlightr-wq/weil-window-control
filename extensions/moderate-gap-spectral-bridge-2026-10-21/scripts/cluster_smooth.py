#!/usr/bin/env python3
"""Same cluster experiment with an ADMISSIBLE ARTZ window.

eta_tau^2 = psi_tau^2 * A / A_tau,  A(x)=cos(sqrt2 x)/(sqrt2 sin(1/sqrt2)) on (-1/2,1/2),
psi_tau a C^infinity cutoff equal to 1 on |x| <= 1/2 - tau and 0 outside (-1/2,1/2).
This is exactly the family of Lamzouri Lemma 3.2 / ARTZ (9)-(11).

Question: for g off-line pairs of equal depth placed at ordinate separation s,
how does the discarded Lemma-5 slack  sum_j (q_j - 2)^2  scale with g?
"""
import math, sys
sys.path.insert(0, __file__.rsplit('/',1)[0])
import linalg as la

TAU = 1.0/200.0
def chi(t):
    if t <= 0: return 0.0
    if t >= 1: return 1.0
    a = math.exp(-1.0/t); b = math.exp(-1.0/(1.0-t))
    return a/(a+b)
def psi(x):
    ax = abs(x)
    if ax <= 0.5-TAU: return 1.0
    if ax >= 0.5: return 0.0
    return chi((0.5-ax)/TAU)
def A(x):
    return math.cos(math.sqrt(2)*x)/(math.sqrt(2)*math.sin(1/math.sqrt(2))) if abs(x) < 0.5 else 0.0
NQ = 24000
_h = 1.0/NQ
_u = [-0.5 + i*_h for i in range(NQ+1)]
_w = [(1.0 if i in (0,NQ) else (4.0 if i % 2 else 2.0))*_h/3.0 for i in range(NQ+1)]
_e2 = [psi(u)**2 * A(u) for u in _u]
_Atau = sum(w*e for w, e in zip(_w, _e2))
_e2 = [e/_Atau for e in _e2]          # now int eta^2 = 1

def ReK(xi, et):
    s = 0.0
    for u, w, e in zip(_u, _w, _e2):
        if e == 0.0: continue
        s += w*e*math.cos(2*math.pi*xi*u)*math.cosh(2*math.pi*et*u)
    return s
def ImK(xi, et):
    s = 0.0
    for u, w, e in zip(_u, _w, _e2):
        if e == 0.0: continue
        s += w*e*math.sin(2*math.pi*xi*u)*math.sinh(2*math.pi*et*u)
    return -s

def slack(offsets, y, cache={}):
    g = len(offsets)
    B = [('g',x,y) for x in offsets] + [('h',x,y) for x in offsets]
    def RK(D, e):
        k = (round(D,9), round(e,9))
        if k not in cache: cache[k] = ReK(D, e)
        return cache[k]
    def IK(D, e):
        k = ('i', round(D,9), round(e,9))
        if k not in cache: cache[k] = ImK(D, e)
        return cache[k]
    def ge(a,b):
        D = a[1]-b[1]; ya, yb = a[2], b[2]
        if a[0]=='g' and b[0]=='g': return 0.5*(RK(D,ya+yb)+RK(D,ya-yb))
        if a[0]=='h' and b[0]=='h': return 0.5*(RK(D,ya+yb)-RK(D,ya-yb))
        if a[0]=='g' and b[0]=='h': return -0.5*(IK(D,ya+yb)-IK(D,ya-yb))
        return -0.5*(IK(-D,ya+yb)-IK(-D,yb-ya))
    n = 2*g
    M = [[ge(B[i],B[j]) for j in range(n)] for i in range(n)]
    q = la.spec_of_DM([2.0]*g + [-2.0]*g, M)
    return sum((v-2.0)**2 for v in q if v > 1e-9)

for y in (0.5, 1.0):
    kap = ReK(0.0, 2*y); ideal = (kap-1)**2
    a_norm = 2*math.pi*y     # = delta log T
    print("\nADMISSIBLE window tau=1/200 : y=%.2f  (a = delta logT = %.3f)  kappa=%.5f  "
          "(kappa-1)^2=%.6g" % (y, a_norm, kap, ideal))
    print("%-6s | %s" % ("sep", "  ".join("g=%-2d" % g for g in (1,2,3,4,6,8,12,16))))
    worst = (1e9, None)
    for sep in (0.15,0.3,0.5,0.7,1.0,1.4,2.0,3.0,6.0):
        row = []
        for g in (1,2,3,4,6,8,12,16):
            s = slack([i*sep for i in range(g)], y)/ideal
            row.append(s)
            if g >= 4 and s/g < worst[0]: worst = (s/g, (sep,g,s))
        print("%-6g |  " % sep + "  ".join("%-7.3f" % v for v in row))
    print("   worst slack-per-pair / (kappa-1)^2 = %.5f  at sep=%g, g=%d (cluster total %.3f)"
          % (worst[0], worst[1][0], worst[1][1], worst[1][2]))
