"""Block inverse iteration on the MIDPOINT of the frozen T=400 even block, to
   choose a 5-dimensional conditioning subspace W.

   This step is NUMERICAL and needs no rigour: it only selects W.  The
   certificate is the interval test in work/condition.py, which holds for
   whatever W comes out of here (A + c*sum v_i v_i^T >= a I implies
   A(f) >= a||f||^2 for every f with <v_i,f> = 0).
"""
import json, time
from mpmath import mp, mpf, matrix, cholesky, sqrt, fabs

mp.dps = 50
raw = json.load(open("wt/extensions/positive-weil-certificate-2026-11-09/"
                     "phase4/certs/block400_odd.json"))
n = len(raw)
import re as _re
_NUM = _re.compile(r"^\[?\s*([+-]?\d+(?:\.\d*)?(?:[eE][+-]?\d+)?)")
def mid(s):            # "[x +/- r]" -> x ;  "[+/-r]" or "[+/- r]" -> 0
    h = s.split("+/-")[0].lstrip("[").strip()
    if not h: return mpf(0)
    m = _NUM.match(h)
    if not m: raise ValueError("unparsable entry %r" % s)
    return mpf(m.group(1))
A = matrix(n, n)
for i in range(n):
    for j in range(n):
        A[i, j] = mid(raw[i][j])
print("loaded %dx%d, diag[0]=%s diag[-1]=%s" % (n, n, mp.nstr(A[0,0],10), mp.nstr(A[n-1,n-1],10)))

t0 = time.time(); Lc = cholesky(A); print("cholesky %.0fs" % (time.time()-t0))
def solve(b):
    y = [mpf(0)]*n
    for i in range(n):
        y[i] = (b[i] - sum(Lc[i,k]*y[k] for k in range(i)))/Lc[i,i]
    x = [mpf(0)]*n
    for i in range(n-1, -1, -1):
        x[i] = (y[i] - sum(Lc[k,i]*x[k] for k in range(i+1, n)))/Lc[i,i]
    return x

K = 6
import random; random.seed(7)
V = [[mpf(random.random()-0.5) for _ in range(n)] for _ in range(K)]
def orth(V):
    out = []
    for v in V:
        for u in out:
            c = sum(a*b for a, b in zip(v, u))
            v = [a - c*b for a, b in zip(v, u)]
        nr = sqrt(sum(a*a for a in v))
        out.append([a/nr for a in v])
    return out
V = orth(V)
for it in range(14):
    V = orth([solve(v) for v in V])
    rq = [sum(v[i]*sum(A[i,j]*v[j] for j in range(n)) for i in range(n)) for v in V]
    rq = sorted(rq)
    print("  iter %2d  Rayleigh quotients: %s" % (it+1, ", ".join(mp.nstr(r,7) for r in rq)), flush=True)
# final: Rayleigh-Ritz in the span
G = matrix(K, K)
AV = [[sum(A[i,j]*v[j] for j in range(n)) for i in range(n)] for v in V]
for p in range(K):
    for q in range(K):
        G[p,q] = sum(V[p][i]*AV[q][i] for i in range(n))
from mpmath import eigsy
ev, Z = eigsy(G)
print("\nRayleigh-Ritz eigenvalues of the %d-dim span:" % K)
for k in range(K): print("   lam_%d ~ %s" % (k+1, mp.nstr(ev[k], 12)))
W = []
for k in range(K):
    w = [sum(Z[p,k]*V[p][i] for p in range(K)) for i in range(n)]
    nr = sqrt(sum(a*a for a in w)); W.append([a/nr for a in w])
# residual check (numerical only)
for k in range(K):
    Aw = [sum(A[i,j]*W[k][j] for j in range(n)) for i in range(n)]
    r = sqrt(sum((Aw[i]-ev[k]*W[k][i])**2 for i in range(n)))
    print("   ||A w_%d - lam_%d w_%d|| = %s" % (k+1,k+1,k+1, mp.nstr(r,6)))
json.dump({"ritz": [mp.nstr(e, 20) for e in ev],
           "vectors": [[mp.nstr(x, 34) for x in w] for w in W]},
          open("src/lowvecs_odd.json", "w"))
print("\nwrote work/lowvecs_odd.json")
