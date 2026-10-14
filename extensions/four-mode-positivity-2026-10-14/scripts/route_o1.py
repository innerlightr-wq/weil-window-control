"""ROUTE O, step 1 (O(1/M^3) version): certify  G >= m0 I,  m0 > 0,
on the four-mode subspace, from the explicit formula alone.

Using tail3:   Q_zeta(x) = Qsharp3(x) +- R'(x),
   Qsharp3(x) = Qcore3(x) - (3/4)||x||^2 S1 + 2 S3 sum_k w_k^2 x_k^2
is EXACTLY a quadratic form with interval coefficients (C'(0) = 0, and
C''(0) = -sum w_k^2 x_k^2 exactly), and
   R'(x) <= 2 sup|C'''|(x) * S4bound  <=  rho3 ||x||^2 ,
with rho3 from the row-sum bound on the bilinear coefficient matrices.
Then  G >= Gt3_mid - (sigma3 + rho3) I  and lambda_min(Gt3_mid) is certified by
exact rational LDL.  Since 4 a a^T >= 0, the same m0 bounds H(delta) for EVERY delta.
"""
import sys, os, json, resource
sys.path.insert(0, os.path.dirname(__file__))
resource.setrlimit(resource.RLIMIT_AS, (4 << 30, 4 << 30))
from fractions import Fraction as Q
import weil_core as C
import tail3 as T3
import ratlin as R
V = C.V; Iv = C.Iv; NDIM = C.NDIM; U, M = C.U, C.M
e = lambda j: [Q(1) if i == j else Q(0) for i in range(NDIM)]

# ---- the sharp-core matrix: polarize Qsharp3 (quadratic part only) ----------
S1, S3 = T3.S1_sharp(), T3.S3_enclosure()

def Qsharp3(x):
    d = T3.Qcore3(x)
    C0 = d['C0']
    w2x2 = Iv(0)
    for k in range(NDIM):
        w = C.wk(k)
        w2x2 = w2x2 + w.sq() * Iv(Q(x[k]) * Q(x[k]))
    return d['core'] - Iv(Q(3, 4)) * Iv(C0) * S1 + Iv(2) * S3 * w2x2

diag = [Qsharp3(e(j)) for j in range(NDIM)]
Gt = [[None]*NDIM for _ in range(NDIM)]
for j in range(NDIM):
    Gt[j][j] = diag[j]
for j in range(NDIM):
    for l in range(j+1, NDIM):
        v = [Q(1) if i in (j, l) else Q(0) for i in range(NDIM)]
        Gt[j][l] = Gt[l][j] = (Qsharp3(v) - diag[j] - diag[l]) * Iv(Q(1, 2))
Gt_mid = [[(Gt[j][l].lo + Gt[j][l].hi)/2 for l in range(NDIM)] for j in range(NDIM)]
rad = [[(Gt[j][l].hi - Gt[j][l].lo)/2 for l in range(NDIM)] for j in range(NDIM)]
sigma = max(sum(rad[j]) for j in range(NDIM))
print('Gt3_mid  (matrix of Qsharp3):')
for j in range(NDIM):
    print('   [' + '  '.join(f'{float(Gt_mid[j][l]):+.14f}' for l in range(NDIM)) + ' ]')
print(f'  sigma3 = max row sum of radii = {float(sigma):.4e}')

# ---- rho3 : uniform bound on R'(x)/||x||^2 ---------------------------------
def ABcoef():
    def AB(x):
        C._set_X(x); return V.C_canonical()
    dA = [AB(e(j)) for j in range(NDIM)]
    al = [[[Iv(0)]*NDIM for _ in range(NDIM)] for _ in range(NDIM)]
    be = [[[Iv(0)]*NDIM for _ in range(NDIM)] for _ in range(NDIM)]
    for k in range(NDIM):
        for j in range(NDIM):
            al[k][j][j] = dA[j][0][k]; be[k][j][j] = dA[j][1][k]
    for j in range(NDIM):
        for l in range(j+1, NDIM):
            Aj, Bj = AB([Q(1) if i in (j, l) else Q(0) for i in range(NDIM)])
            for k in range(NDIM):
                al[k][j][l] = al[k][l][j] = (Aj[k]-dA[j][0][k]-dA[l][0][k])*Iv(Q(1,2))
                be[k][j][l] = be[k][l][j] = (Bj[k]-dA[j][1][k]-dA[l][1][k])*Iv(Q(1,2))
    return al, be

al, be = ABcoef()
absmax = lambda iv: max(abs(iv.lo), abs(iv.hi))
K3 = Iv(0)
for k in range(NDIM):
    w = C.wk(k)
    ra = max(sum(absmax(al[k][j][l]) for l in range(NDIM)) for j in range(NDIM))
    rb = max(sum(absmax(be[k][j][l]) for l in range(NDIM)) for j in range(NDIM))
    K3 = K3 + Iv(ra) * w * w.sq() + Iv(rb) * (Iv(3) * w.sq() + Iv(U) * w * w.sq())
rho3 = Iv(2) * K3 * T3.S4_bound()
print(f'  Kappa3 (uniform bound on sup|C\'\'\'|(x)/||x||^2) <= {float(K3.hi):.4f}')
print(f'  rho3 <= {float(rho3.hi):.4e}')
eps = Iv(sigma) + Iv(0, rho3.hi)
print(f'  sigma3 + rho3 <= {float(eps.hi):.4e}')
for lab, x in (('frozen', C.XFROZEN), ('e_0', e(0)), ('e_3', e(3))):
    d = T3.Q_zeta3(x); nsq = sum(Q(t)*Q(t) for t in x)
    pv = float(d['tail3_info']['Rmag'].hi / nsq)
    print(f'  check {lab:>7}: per-vector |R\'|/||x||^2 = {pv:.4e} <= rho3  '
          f'{pv <= float(rho3.hi)}')

# ---- certify lambda_min(Gt3_mid) ------------------------------------------
lo, hi = Q(0), Q(1)
for _ in range(160):
    t = (lo + hi)/2
    B = [[Gt_mid[j][l] - (t if j == l else 0) for l in range(NDIM)] for j in range(NDIM)]
    d = R.cholesky_pd(B)
    if d is not None and all(p > 0 for p in d): lo = t
    else: hi = t
print(f'\n  lambda_min(Gt3_mid) in [{float(lo):.6e}, {float(hi):.6e}]')
br = R.eig_brackets(Gt_mid)
print('  full spectrum (each eigenvalue bracketed by exact inertia bisection,')
print('  #{eig < t} = n_neg(Gt3_mid - t I); NO deflation):')
spec = [(a + b)/2 for a, b in br]
for a, b in br:
    print(f'   [{float(a):.9e}, {float(b):.9e}]  width {float(b-a):.2e}')
assert br[0][0] <= lo and hi <= br[0][1] or (br[0][0] <= hi and lo <= br[0][1]), \
    'the smallest inertia bracket must agree with the LDL bisection for lambda_min'
print(f'   condition number <= {float(br[-1][1]/br[0][0]):.3e}')

tcert = Q(int(float(lo)*10**14), 10**14)
Bc = [[Gt_mid[j][l] - (tcert if j == l else 0) for l in range(NDIM)] for j in range(NDIM)]
piv = R.cholesky_pd(Bc)
ok = piv is not None and all(p > 0 for p in piv)
m0 = Iv(tcert) - eps
print(f'\n  t_cert = {tcert} = {float(tcert):.6e}')
print(f'  LDL pivots of Gt3_mid - t_cert I: ' + ', '.join(f'{float(p):.4e}' for p in piv))
print(f'  all pivots > 0 : {ok}')
print(f'  => G >= m0 I with m0 >= {float(m0.lo):.6e}    POSITIVE: {m0.lo > 0}')
print(f'  => H(delta) >= m0 I for EVERY delta, since 4 a a^T >= 0')
json.dump(dict(Gt3_mid=[[str(v) for v in row] for row in Gt_mid],
               rad=[[str(v) for v in row] for row in rad], sigma3=str(sigma),
               Kappa3=[str(K3.lo), str(K3.hi)], rho3=[str(rho3.lo), str(rho3.hi)],
               eps=[str(eps.lo), str(eps.hi)], lam_min=[str(lo), str(hi)],
               spectrum=[str(s) for s in spec], t_cert=str(tcert),
               ldl_pivots=[str(p) for p in piv], m0_lo=str(m0.lo),
               S1=[str(S1.lo), str(S1.hi)], S3=[str(S3.lo), str(S3.hi)], M=M),
          open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'results', 'route_o_step1.json'), 'w'), indent=1)
print(f'\npeak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024:.1f} MiB')
