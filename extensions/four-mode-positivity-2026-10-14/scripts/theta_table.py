"""NUMERICAL DIAGNOSTIC (not a certificate): Theta(delta) = 4 b^T H^{-1} b with
H = G + 4 a a^T, evaluated at the midpoint of the certified interval matrix for
G (Gt3_mid, from the O(1/M^3) tail) and at the midpoints of a(delta), b(delta).
Proposition 1: M(delta) = H - 4 b b^T has a negative direction iff Theta > 1
(valid when H > 0, which exact LDL confirms at each delta below).

This supersedes an earlier table computed with the TAIL-FREE matrix G_core.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Q
import weil_core as C, ratlin as R
NDIM = C.NDIM
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
s1 = json.load(open(os.path.join(HERE, 'results', 'route_o_step1.json')))
G = [[Q(v) for v in row] for row in s1['Gt3_mid']]
Gc = C.G_core()
Gcm = [[(Gc[j][l].lo + Gc[j][l].hi)/2 for l in range(NDIM)] for j in range(NDIM)]
mid = lambda iv: (iv.lo + iv.hi)/2
print(f"lambda_min(G)      bracket {[float(x) for x in R.min_eig_bracket(G)]}")
print(f"lambda_min(G_core) bracket {[float(x) for x in R.min_eig_bracket(Gcm)]}")
print(f"the omitted term is -(3/4)||x||^2 S1 : S1 enclosure in step1 = "
      f"{[float(Q(t)) for t in s1['S1']]}")
print()
hdr = ['13/200','133/2000','27/400','137/2000','7/100','17/200','9/100','1/10','1/5','1/4']
print(f"{'delta':>12} {'Theta(G)':>14} {'inertia(M)':>12} {'Theta(G_core)':>15} {'H>0':>5}")
rows = []
for h in hdr:
    d = Q(h)
    a, b = C.ab_vectors(d)
    a1 = [mid(t) for t in a]; b1 = [mid(t) for t in b]
    for lab, Gm in (('G', G), ('Gcore', Gcm)):
        H = [[Gm[j][l] + 4*a1[j]*a1[l] for l in range(NDIM)] for j in range(NDIM)]
        piv = R.cholesky_pd(H)
        y = R.solve(H, b1); th = 4*R.dot(b1, y)
        if lab == 'G':
            thG = th; hpos = piv is not None and all(p > 0 for p in piv)
            M = [[H[j][l] - 4*b1[j]*b1[l] for l in range(NDIM)] for j in range(NDIM)]
            inert = R.inertia(M)
        else:
            thC = th
    print(f"{h:>12} {float(thG):>14.6f} {str(inert):>12} {float(thC):>15.6f} {str(hpos):>5}")
    rows.append(dict(delta=h, theta_G=str(thG), theta_Gcore=str(thC),
                     inertia=list(inert), H_pd=hpos))
json.dump(rows, open(os.path.join(HERE, 'results', 'theta_table.json'), 'w'), indent=1)
