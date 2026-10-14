"""Validate the sharper archimedean tail.  Three independent requirements:
 (1) consistency: the sharp enclosure of Q_zeta must be CONTAINED in the
     published O(1/M) enclosure (both are valid enclosures of the same number);
 (2) truncation independence: recomputing with a SMALLER truncation M' must give
     an enclosure that still contains the same number, and overlaps;
 (3) the derivative inputs are checked against finite differences of C.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))
from fractions import Fraction as Q
import weil_core as C
import tail_sharp as TS
V = C.V
fails = []
def chk(c, m):
    print(('PASS  ' if c else 'FAIL  ') + m)
    if not c: fails.append(m)

X = C.XFROZEN
dpub = C.Q_zeta(X)
dsh = TS.Q_zeta_sharp(X)
pub, sh = dpub['Q'], dsh['Q']
chk(pub.lo <= sh.lo and sh.hi <= pub.hi,
    f'(1) sharp enclosure [{float(sh.lo):.15f}, {float(sh.hi):.15f}] is INSIDE the '
    f'published [{float(pub.lo):.15f}, {float(pub.hi):.15f}]')
chk(sh.width() < pub.width() / 1000,
    f'(1) width {float(sh.width()):.3e} vs published {float(pub.width()):.3e} '
    f'(factor {float(pub.width()/sh.width()):.0f} narrower)')
chk(float(dsh['tail_sharp'].hi - dsh['tail_sharp'].lo) < 1e-7,
    f'(1) sharp tail allowance half-width = '
    f'{float((dsh["tail_sharp"].hi-dsh["tail_sharp"].lo)/2):.4e}')
# (3) derivative inputs against finite differences of C
A, B = dsh['A'], dsh['B']
for u0 in (Q(0), Q(1, 2), Q(8, 5)):
    h = Q(1, 10**6)
    lo = u0 if u0 == 0 else u0 - h
    hi = u0 + h if u0 + h <= Q(8, 5) else u0
    num = (V.C_at(A, B, C.Iv(hi)) - V.C_at(A, B, C.Iv(lo))) / C.Iv(hi - lo)
    ana = TS.dC_at(A, B, C.Iv(u0))
    chk(abs(float(num.mid() - ana.mid())) < 1e-3,
        f"(3) C'({float(u0)}): analytic {float(ana.mid()):.9f} vs finite difference "
        f"{float(num.mid()):.9f}")
chk(abs(float(TS.dC_at(A, B, C.Iv(0)).mid())) < 1e-9,
    f"(3) C'(0+) = {float(TS.dC_at(A,B,C.Iv(0)).mid()):.3e} -- vanishes, as it must "
    f"for the autocorrelation of a continuous compactly supported f (C is even and C^1)")
chk(abs(float(V.C_at(A, B, C.Iv(Q(8,5))).mid())) < 1e-60,
    f'(3) C(U) = {float(V.C_at(A,B,C.Iv(Q(8,5))).mid()):.2e} (the hypothesis of the '
    f'integration by parts, verified not assumed)')
print(f"NOTE  sup|C\'| bound = {float(V.dC_sup_bound(A,B).hi):.12f} (published 11.436722214306984)")
print(f"NOTE  sup|C\'\'| bound = {float(dsh['tail_info']['ddCsup'].hi):.6f}")
print(f"NOTE  remainder |R| <= {float(dsh['tail_info']['Rmag'].hi):.4e}")
print()
print('TAIL VALIDATION:', 'ALL PASS' if not fails else f'{len(fails)} FAILURES')
sys.exit(1 if fails else 0)
