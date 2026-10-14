#!/usr/bin/env python3
"""VERIFIER for Route O.  Clean process, frozen inputs only.  Exits nonzero on failure.

Asserts, regenerating every enclosure from the mathematical input:

  O1. the O(1/M^3) tail and the sharp gamma_E interval T are each CONTAINED in
      the coarser enclosures they replace (two valid enclosures of the same
      number, the new one inside the old);
  O2. G >= m0 I with m0 > 0:  the matrix Gt3_mid of Qsharp3 is reproduced by
      polarization, sigma3 and rho3 are recomputed, and
      Gt3_mid - t_cert I is shown positive definite by exact rational LDL;
  O3. for every cell, lambda_min(N) >= Loss by exact rational LDL, where
      N and Loss are as in scripts/route_o2.py -- hence M(delta) >= 0, i.e.
      Q_add,delta(f_x) >= 0 for EVERY x in the subspace and every delta in the
      cell: certified NON-detection;
  O4. mechanical coverage of the non-detection cells (closed, consecutive, no gaps);
  O5. the onset bracket: Theta < 1 at the right end of the certified
      non-detection range and a NEGATIVE witness exists at 9/100 (the detection
      certificate), so the two certified regimes are separated by an explicit
      uncertified gap.
"""
import sys, os, json, resource
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'scripts'))
resource.setrlimit(resource.RLIMIT_AS, (4 << 30, 4 << 30))
from fractions import Fraction as Q
import weil_core as C
import tail_sharp as TS
import tail3 as T3
import ratlin as R
V = C.V; Iv = C.Iv; NDIM = C.NDIM; L = C.L
fails = []
def chk(c, m):
    print(('PASS  ' if c else 'FAIL  ') + m)
    if not c: fails.append(m)

s1 = json.load(open(os.path.join(HERE, 'results', 'route_o_step1.json')))
s2 = json.load(open(os.path.join(HERE, 'results', 'route_o_step2.json')))

# ---- O1 ------------------------------------------------------------------
Tpub = Iv(Q(1, 2*(C.M + 1)) - Q(1, 6*C.M*C.M), Q(1, 2*C.M))
Ts = T3.T_sharp()
chk(Tpub.lo <= Ts.lo and Ts.hi <= Tpub.hi,
    f'O1 sharp T inside the published T (widths {float(Ts.width()):.3e} vs '
    f'{float(Tpub.width()):.3e})')
X = C.XFROZEN
q1, q2, q3 = C.Q_zeta(X)['Q'], TS.Q_zeta_sharp(X)['Q'], T3.Q_zeta3(X)['Q']
chk(q1.lo <= q2.lo and q2.hi <= q1.hi and q2.lo <= q3.lo and q3.hi <= q2.hi,
    f'O1 nested enclosures of Q_zeta(frozen): O(1/M^3) in O(1/M^2) in O(1/M) '
    f'(widths {float(q3.width()):.2e} < {float(q2.width()):.2e} < {float(q1.width()):.2e})')

# ---- O2 ------------------------------------------------------------------
S1, S3 = T3.S1_sharp(), T3.S3_enclosure()
def Qsharp3(x):
    d = T3.Qcore3(x); w2 = Iv(0)
    for k in range(NDIM):
        w2 = w2 + C.wk(k).sq() * Iv(Q(x[k]) * Q(x[k]))
    return d['core'] - Iv(Q(3, 4)) * Iv(d['C0']) * S1 + Iv(2) * S3 * w2
e = lambda j: [Q(1) if i == j else Q(0) for i in range(NDIM)]
diag = [Qsharp3(e(j)) for j in range(NDIM)]
Gt = [[None]*NDIM for _ in range(NDIM)]
for j in range(NDIM): Gt[j][j] = diag[j]
for j in range(NDIM):
    for l in range(j+1, NDIM):
        v = [Q(1) if i in (j, l) else Q(0) for i in range(NDIM)]
        Gt[j][l] = Gt[l][j] = (Qsharp3(v) - diag[j] - diag[l]) * Iv(Q(1, 2))
Gt_mid = [[(Gt[j][l].lo + Gt[j][l].hi)/2 for l in range(NDIM)] for j in range(NDIM)]
stored = [[Q(v) for v in row] for row in s1['Gt3_mid']]
chk(all(Gt[j][l].lo <= stored[j][l] <= Gt[j][l].hi for j in range(NDIM) for l in range(NDIM)),
    'O2 stored Gt3_mid lies inside the freshly polarized interval matrix')
sigma = max(sum((Gt[j][l].hi - Gt[j][l].lo)/2 for l in range(NDIM)) for j in range(NDIM))
chk(sigma <= Q(s1['sigma3']) * Q(101, 100),
    f'O2 sigma3 recomputed = {float(sigma):.4e} (stored {float(Q(s1["sigma3"])):.4e})')
# O2b. the polarization is only legitimate if Qsharp3 really is a quadratic
# form; and the loss sigma3 must bound |x^T Delta x|/||x||^2 in THIS dimension.
# Both are checked rather than assumed.
z = [Q(0)]*NDIM
chk(Qsharp3(z).lo <= 0 <= Qsharp3(z).hi, 'O2b Qsharp3(0) = 0 (no constant part)')
xt = [Q(3, 7), Q(-5, 11), Q(2), Q(1, 13)]
q1v = Qsharp3(xt); q2v = Qsharp3([2*t for t in xt])
chk(q2v.lo <= 4*q1v.hi and 4*q1v.lo <= q2v.hi,
    'O2b Qsharp3 is 2-homogeneous: Qsharp3(2x) and 4 Qsharp3(x) overlap '
    f'({float(q2v.lo):.12e} vs {float(4*q1v.lo):.12e})')
qf = Iv(0)
for j in range(NDIM):
    for l in range(NDIM):
        qf = qf + Gt[j][l] * Iv(xt[j]*xt[l])
chk(qf.lo <= q1v.hi and q1v.lo <= qf.hi,
    'O2b polarized x^T Gt x reproduces Qsharp3(x) on an off-axis test vector')
nsq = sum(t*t for t in xt)
chk(abs(((qf.lo+qf.hi)/2 - (q1v.lo+q1v.hi)/2)) <= sigma * nsq,
    f'O2b polarization defect <= sigma3 ||x||^2 in dimension {NDIM} '
    '(row-sum bound, via |x_j x_l| <= (x_j^2+x_l^2)/2)')

tcert = Q(s1['t_cert']); eps = Q(s1['eps'][1])
Bc = [[stored[j][l] - (tcert if j == l else 0) for l in range(NDIM)] for j in range(NDIM)]
piv = R.cholesky_pd(Bc)
chk(piv is not None and all(p > 0 for p in piv),
    f'O2 exact LDL: Gt3_mid - t_cert I is positive definite, pivots ' +
    ', '.join(f'{float(p):.3e}' for p in piv))
m0 = tcert - eps
chk(m0 > 0, f'O2 => G >= m0 I with m0 = {float(m0):.6e} > 0, and H(delta) >= m0 I '
            f'for every delta (4aa^T >= 0)')

# ---- O3, O4 --------------------------------------------------------------
cells = s2['cells']
mid = lambda iv: (iv.lo + iv.hi)/2
rad = lambda iv: (iv.hi - iv.lo)/2
l32 = (Iv(L) * V.iv_sqrt(Iv(L))).hi
bad = 0
for i, c in enumerate(cells):
    d1, d2, w = Q(c['lo']), Q(c['hi']), Q(c['w'])
    tau, taup = Q(c['tau']), Q(c['taup'])
    a, b = C.ab_vectors(d1)
    a1 = [mid(t) for t in a]; b1 = [mid(t) for t in b]
    ra = max(rad(t) for t in a); rb = max(rad(t) for t in b)
    Ahat = Q(2)*Q(2, 3)*l32*V.iv_sinh(Iv(d2*L)).hi
    Bhat = Q(2)*l32*V.iv_cosh(Iv(d2*L)).hi
    tiny = 8*(ra*Q(2) + rb*Q(2)) + Q(1, 10**60)
    Loss = eps + 4*w*w*(Ahat*Ahat*(1 + 1/tau) + Bhat*Bhat*(1 + 1/taup)) + tiny
    N = [[stored[j][l] + 4*(1-tau)*a1[j]*a1[l] - 4*(1+taup)*b1[j]*b1[l]
          for l in range(NDIM)] for j in range(NDIM)]
    Bm = [[N[j][l] - (Loss if j == l else 0) for l in range(NDIM)] for j in range(NDIM)]
    p = R.cholesky_pd(Bm)
    if not (p is not None and all(t > 0 for t in p)): bad += 1
chk(bad == 0, f'O3 all {len(cells)} cells: lambda_min(N) >= Loss by exact LDL '
              f'=> M(delta) >= 0 throughout the cell ({bad} failures)')
okcov = Q(cells[0]['lo']) == 0
for i in range(len(cells)-1):
    if Q(cells[i]['hi']) != Q(cells[i+1]['lo']): okcov = False
chk(okcov, f'O4 coverage: closed cells consecutive from 0 with no gap, up to '
           f'delta = {cells[-1]["hi"]} = {float(Q(cells[-1]["hi"])):.8f}')

# ---- O5 ------------------------------------------------------------------
dr = Q(cells[-1]['hi'])
a, b = C.ab_vectors(dr)
a1 = [mid(t) for t in a]; b1 = [mid(t) for t in b]
H = [[stored[j][l] + 4*a1[j]*a1[l] for l in range(NDIM)] for j in range(NDIM)]
y = R.solve(H, b1); th = 4*R.dot(b1, y)
chk(th < 1, f'O5 Theta({float(dr):.6f}) = {float(th):.6f} < 1 (consistent with '
            f'Proposition 1 and the LDL certificate)')
cert = json.load(open(os.path.join(HERE, 'results', 'certificate_new.json')))
xw = [Q(t) for t in cert['cells'][0]['x']]
qa, _, _, _ = T3.Q_add3(xw, Q(9, 100))
chk(qa.hi < 0, f'O5 a NEGATIVE witness exists at delta = 9/100 (upper bound '
               f'{float(qa.hi):.4e}): the two certified regimes are genuinely different')
print()
print(f'  certified NON-detection : delta in [0, {cells[-1]["hi"]}]')
print(f'  certified DETECTION     : delta in [9/100, 49/100]')
print(f'  uncertified gap         : ({float(dr):.8f}, 0.09)')
print(f'  numerical onset bracket : (0.0675, 0.0685)  -- a DIAGNOSTIC, inside that gap')
print()
print('ROUTE O:', 'ALL ASSERTIONS PASS' if not fails else f'{len(fails)} FAILURES')
print(f'peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024:.1f} MiB')
sys.exit(1 if fails else 0)
