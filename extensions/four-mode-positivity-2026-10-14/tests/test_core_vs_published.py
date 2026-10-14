"""Validate weil_core against the PUBLISHED certificate numbers at the frozen vector.
Exits nonzero on any mismatch."""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))
from fractions import Fraction as Q
import weil_core as C
fails = []
def chk(c, m):
    print(('PASS  ' if c else 'FAIL  ') + m)
    if not c: fails.append(m)

t0 = time.time()
C.arch_sums(progress=True)
print(f'  arch_sums built in {time.time()-t0:.1f}s (M = {C.M})', file=sys.stderr)

X = C.XFROZEN
d = C.Q_zeta(X)
# published figures, CERTIFICATE_THEOREM.md section 5.5
pub = {'pole': (0.45715833231743186, 4.8e-77),
       'logpi': (1.1447298858493788, 5.9e-78),
       'prime': (0.51733913650311314, 1.7e-76)}
for k, (v, _) in pub.items():
    iv = d[k]
    chk(float(iv.lo) <= v <= float(iv.hi), f'{k}: published {v} inside recomputed {iv}')
chk(float(d['arch'].lo if 'arch' in d else d['arch_core'].lo + d['tail'].lo) <= 1.2100349293664769
    <= float(d['arch_core'].hi + d['tail'].hi),
    'arch: published [1.2100349293664769, 1.2101336300807928] inside recomputed')
QZ = d['Q']
chk(float(QZ.lo) <= 0.0051242393314170015 and 0.0052229400457328018 <= float(QZ.hi),
    f'Q_zeta: published [0.00512423933, 0.00522294005] inside recomputed '
    f'[{float(QZ.lo):.11f}, {float(QZ.hi):.11f}]')
chk(QZ.lo > 0, f'Q_zeta(frozen) > 0 : lo = {float(QZ.lo):.6e}')
chk(float(d['tail'].hi) <= 4.9350337755161665e-05*1.0000001,
    f'archimedean tail allowance = +-{float(d["tail"].hi):.9e} (published 4.935033776e-05)')

Qa, Ad, Bd, _ = C.Q_add(X, Q(2, 5))
chk(Qa.hi <= Q(-145107657498227, 10**16),
    f'point: Q_add,2/5(frozen) <= -0.0145107657498227 ; recomputed hi = {float(Qa.hi):.16f}')
chk(float(Qa.lo) <= -0.0146094664641386 + 1e-15,
    f'point: recomputed lo = {float(Qa.lo):.16f} (published -0.0146094664641386)')
Fg = C.F_at_gamma(X)
chk(abs(float(Fg.mid())) < 1e-12, f'node residual F(14) = {float(Fg.mid()):.6e} (published 4.129e-13)')

# G_core consistency: x^T G_core x must equal Qcore(x) for the frozen vector
G = C.G_core()
q1 = C.quad(G, X); q2 = C.Qcore(X)['core']
ov_lo, ov_hi = max(q1.lo, q2.lo), min(q1.hi, q2.hi)
chk(ov_lo <= ov_hi, f'polarization: x^T G_core x = {q1} overlaps Qcore(x) = {q2}')
# The enclosure width is dominated by the gamma_E-elimination interval
# T in [1/(2(M+1)) - 1/(6M^2), 1/(2M)], of width (2/3)/M^2 ~ 3.9e-11 per unit of
# C0 -- NOT by the 2^-260 grid.  It is 6 orders below the archimedean tail
# allowance 4.94e-5, so it never limits anything here.
Twidth = float(Q(1, 2*C.M) - Q(1, 2*(C.M+1)) + Q(1, 6*C.M*C.M))
chk(q2.width() <= Q(3, 2) * (Q(1, 2*C.M) - Q(1, 2*(C.M+1)) + Q(1, 6*C.M*C.M)),
    f'Qcore width {float(q2.width()):.3e} is the gamma_E term T (width {Twidth:.3e}) '
    f'times C0, not grid rounding')
chk(q1.width() < Q(1, 10**8),
    f'polarized x^T G_core x width {float(q1.width()):.3e} < 1e-8, i.e. 6 orders '
    f'below the archimedean tail allowance 4.94e-5')
# a,b vectors reproduce the scalar quartet
a, b = C.ab_vectors(Q(2, 5))
da, db = C.dot(a, X), C.dot(b, X)
chk(da.lo <= Ad.mid() <= da.hi and db.lo <= Bd.mid() <= db.hi,
    'a.x and b.x reproduce the scalar A_delta, B_delta at delta = 2/5')
# delta = 0 multiplicity limit: b = 0 and Q_add,0 = Q_zeta + 4 F(gamma)^2
a0, b0 = C.ab_vectors(Q(0))
chk(all(bb.lo == 0 == bb.hi for bb in b0), 'delta = 0: b = 0 exactly')
Q0, A0, B0, _ = C.Q_add(X, Q(0))
chk(abs(float((Q0 - (QZ + C.Iv(4) * Fg.sq())).mid())) < 1e-60,
    'delta = 0: Q_add,0 = Q_zeta + 4 F(gamma)^2 (multiplicity limit)')

print()
print('CORE VALIDATION:', 'ALL PASS' if not fails else f'{len(fails)} FAILURES')
import resource
print(f'peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024:.1f} MiB')
sys.exit(1 if fails else 0)
