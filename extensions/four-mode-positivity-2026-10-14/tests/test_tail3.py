"""Validate the O(1/M^3) tail and the sharper T against the two coarser routes."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))
from fractions import Fraction as Q
import weil_core as C
import tail_sharp as TS
import tail3 as T3
V = C.V
fails = []
def chk(c, m):
    print(('PASS  ' if c else 'FAIL  ') + m)
    if not c: fails.append(m)

Tpub = C.Iv(Q(1, 2*(C.M + 1)) - Q(1, 6*C.M*C.M), Q(1, 2*C.M))
Ts = T3.T_sharp()
chk(Tpub.lo <= Ts.lo and Ts.hi <= Tpub.hi,
    f'T: sharp [{float(Ts.lo):.18f}, {float(Ts.hi):.18f}] inside published '
    f'[{float(Tpub.lo):.18f}, {float(Tpub.hi):.18f}]')
chk(Ts.width() < Q(1, 10**14),
    f'T: width {float(Ts.width()):.3e} (published {float(Tpub.width()):.3e}), '
    f'factor {float(Tpub.width()/Ts.width()):.3g} narrower')

X = C.XFROZEN
for lab, x in (('frozen', X), ('e_0', [Q(1),Q(0),Q(0),Q(0)]), ('e_3', [Q(0),Q(0),Q(0),Q(1)])):
    d2 = TS.Q_zeta_sharp(x); d3 = T3.Q_zeta3(x); dp = C.Q_zeta(x)
    q2, q3, qp = d2['Q'], d3['Q'], dp['Q']
    chk(q2.lo <= q3.lo and q3.hi <= q2.hi,
        f'{lab}: O(1/M^3) enclosure inside the O(1/M^2) one '
        f'(widths {float(q3.width()):.3e} vs {float(q2.width()):.3e})')
    chk(qp.lo <= q3.lo and q3.hi <= qp.hi, f'{lab}: and inside the published O(1/M) one')
    hw = float(d3['tail3_info']['Rmag'].hi)
    chk(hw < 1e-11, f'{lab}: O(1/M^3) remainder |R\'| <= {hw:.4e}')

# the two endpoint terms that were shown to vanish identically
A, B = T3.Qcore3(X)['A'], T3.Qcore3(X)['B']
chk(abs(float(T3.dC_at_U(A, B).mid())) < 1e-60,
    f"C'(U) = {float(T3.dC_at_U(A,B).mid()):.2e} vanishes (cos(w_k U) = -1 and C'(0)=0)")
chk(abs(float(T3.ddC_at(A, B, C.Iv(Q(8,5))).mid())) < 1e-60,
    f"C''(U) = {float(T3.ddC_at(A,B,C.Iv(Q(8,5))).mid()):.2e} vanishes (sin(w_k U) = 0)")
# C''(0) = -sum w_k^2 x_k^2, checked against the canonical-form value
cc = T3.ddC_at(A, B, C.Iv(0))
pred = -sum(float(C.wk(k).mid())**2 * float(X[k])**2 for k in range(4))
chk(abs(float(cc.mid()) - pred) < 1e-9 * abs(pred),
    f"C''(0) = {float(cc.mid()):.12f} = -sum w_k^2 x_k^2 = {pred:.12f}  (= -||f'||^2)")
print()
print('TAIL3 VALIDATION:', 'ALL PASS' if not fails else f'{len(fails)} FAILURES')
sys.exit(1 if fails else 0)
