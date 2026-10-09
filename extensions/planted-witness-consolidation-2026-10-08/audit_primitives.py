#!/usr/bin/env python3
r"""Adversarial audit of the rational interval primitives of verify_certificate_consolidated.py.

Every test is a PROPERTY that must hold for a correct outward-rounded interval type, checked
against exact rational arithmetic or against a proved inequality -- never against libm.
Any failure prints FAIL and the script exits 1.
"""
import importlib.util, random, sys
from fractions import Fraction as Q

spec = importlib.util.spec_from_file_location('vc', 'verify_certificate_consolidated.py')
vc = importlib.util.module_from_spec(spec)
sys.modules['vc'] = vc
# the module executes main() under __main__ only
spec.loader.exec_module(vc)

Iv, CF = vc.Iv, vc.CertificationFailure
fails = []


def chk(name, cond, extra=''):
    if cond:
        print(f'  ok    {name}')
    else:
        print(f'  FAIL  {name}  {extra}')
        fails.append(name)


def encloses(iv, exact):
    return iv.lo <= exact <= iv.hi


random.seed(20261008)
RS = [Q(random.randint(-10**7, 10**7), random.randint(1, 10**7)) for _ in range(400)]


print('A. outward rounding and containment (the defining property)')
chk('lo <= hi for every literal', all(Iv(r).lo <= Iv(r).hi for r in RS))
chk('a rational literal is enclosed', all(encloses(Iv(r), r) for r in RS))
chk('negatives enclosed (rounding direction not flipped)',
    all(encloses(Iv(r), r) for r in RS if r < 0))


try:
    Iv(Q(16), Q(3))
    chk('empty interval rejected', False, 'no exception')
except CF:
    chk('empty interval rejected', True)

print('B. arithmetic is a sound enclosure')
bad = []
for i in range(0, len(RS) - 1, 2):
    p, q = RS[i], RS[i+1]
    P, Qi = Iv(p), Iv(q)
    for nm, iv, ex in (('add', P + Qi, p + q), ('sub', P - Qi, p - q),
                       ('mul', P * Qi, p * q), ('neg', -P, -p), ('sq', P.sq(), p * p)):
        if not encloses(iv, ex):
            bad.append((nm, p, q))
    if q != 0 and not encloses(P / Qi, p / q):
        bad.append(('div', p, q))
chk('add/sub/mul/neg/sq/div all enclose the exact value', not bad, bad[:3])

print('C. sign handling in mul and sq (all four quadrants)')
quad = [(Q(2), Q(3)), (Q(-2), Q(3)), (Q(2), Q(-3)), (Q(-2), Q(-3))]
chk('mul correct in all sign combinations',
    all(encloses(Iv(a) * Iv(b), a * b) for a, b in quad))
_m = Iv(Q(-1), Q(2)) * Iv(Q(-3), Q(5))
_prods = [Q(-1)*Q(-3), Q(-1)*Q(5), Q(2)*Q(-3), Q(2)*Q(5)]
chk('mul of straddling intervals encloses every endpoint product',
    _m.lo <= min(_prods) and _m.hi >= max(_prods), f'{float(_m.lo)},{float(_m.hi)}')
s = Iv(Q(-1), Q(2)).sq()
chk('sq of a straddling interval has lo <= 0 and hi >= 4', s.lo <= 0 and s.hi >= 4,
    f'{float(s.lo)},{float(s.hi)}')
chk('sq is never negative', Iv(Q(-7), Q(-3)).sq().lo >= 9 - Q(1, 10**9))

print('D. division by an interval containing zero must be refused')
for z in (Iv(Q(-1), Q(2)), Iv(Q(0)), Iv(Q(0), Q(5)), Iv(Q(-5), Q(0))):
    try:
        Iv(Q(1)) / z
        chk(f'1/{z!r} refused', False, 'no exception')
    except CF:
        chk(f'1/[{float(z.lo):.3g},{float(z.hi):.3g}] refused', True)

print('E. sqrt: brackets the true root, verified by squaring (no libm)')
bad = []
for r in RS:
    if r <= 0:
        continue
    s = vc.iv_sqrt(Iv(r))
    if not (s.lo >= 0 and s.lo * s.lo <= r <= s.hi * s.hi):
        bad.append(r)
chk('lo^2 <= x <= hi^2 for 400 random positives', not bad, bad[:3])
chk('sqrt(0) = 0', vc.iv_sqrt(Iv(Q(0))).lo == 0)
try:
    vc.iv_sqrt(Iv(Q(-1)))
    chk('sqrt of a negative refused', False, 'no exception')
except CF:
    chk('sqrt of a negative refused', True)

print('F. exp: proved functional equations, checked as interval statements')
e1 = vc.iv_exp(Iv(Q(1)))
chk('exp(1) in (2.718281828, 2.718281829)',
    e1.lo > Q(2718281828, 10**9) and e1.hi < Q(2718281829, 10**9))
chk('exp(0) = 1 exactly enclosed', encloses(vc.iv_exp(Iv(Q(0))), Q(1)))
bad = []
for t in [Q(1, 3), Q(-7, 5), Q(13, 2), Q(-64, 5), Q(-372, 5)]:
    a, b = vc.iv_exp(Iv(t)), vc.iv_exp(Iv(-t))
    pr = a * b
    if not (pr.lo <= 1 <= pr.hi):
        bad.append(('exp(t)exp(-t)=1', t))
    s2 = vc.iv_exp(Iv(t / 2)).sq()
    if not (s2.lo <= a.hi and a.lo <= s2.hi):
        bad.append(('exp(t/2)^2 = exp(t)', t))
chk('exp(t)exp(-t) encloses 1, and exp(t/2)^2 agrees with exp(t)', not bad, bad[:3])
chk('exp is increasing on the tested grid',
    all(vc.iv_exp(Iv(Q(i, 4))).lo <= vc.iv_exp(Iv(Q(i+1, 4))).hi for i in range(-40, 40)))

print('G. sin/cos: Pythagoras, parity, and range reduction at large argument')
bad = []
for t in [Q(0), Q(1, 7), Q(-3, 2), Q(56, 5), Q(224), Q(-224), Q(31415, 10**4)]:
    s, c = vc.iv_sin(Iv(t)), vc.iv_cos(Iv(t))
    p = s.sq() + c.sq()
    if not (p.lo <= 1 <= p.hi):
        bad.append(('sin^2+cos^2=1', t, float(p.lo), float(p.hi)))
    if not (-1 - Q(1, 10**60) <= s.lo and s.hi <= 1 + Q(1, 10**60)):
        bad.append(('|sin|<=1', t))
    sm, cm = vc.iv_sin(Iv(-t)), vc.iv_cos(Iv(-t))
    if not (sm.lo <= -s.lo and -s.hi <= sm.hi) or not (cm.lo <= c.hi and c.lo <= cm.hi):
        bad.append(('parity', t))
chk('sin^2+cos^2 encloses 1, |sin|,|cos| <= 1, parity holds (incl. t=224=gamma_max*L)',
    not bad, bad[:3])
# the primitive has NO range reduction; with nmax=1200 its proved truncation reaches
# |t| ~ 371.  Beyond that it must REFUSE, never return an unsound enclosure.
try:
    vc._sincos_rat(Q(1120))
    chk('out-of-range argument refused (no silent unsound result)', False, 'no exception')
except CF:
    chk('out-of-range argument t=1120 refused (no silent unsound result)', True)
chk('the largest argument the certificate uses, gamma*L = 224, is inside the range',
    vc._sincos_rat(Q(224)) is not None)
# double angle: a genuinely independent constraint
bad = []
for t in [Q(1, 3), Q(7, 2), Q(56, 5), Q(112)]:
    s, c = vc.iv_sin(Iv(t)), vc.iv_cos(Iv(t))
    s2, c2 = vc.iv_sin(Iv(2*t)), vc.iv_cos(Iv(2*t))
    lhs, rhs = Iv(Q(2)) * s * c, s2
    if not (lhs.lo <= rhs.hi and rhs.lo <= lhs.hi):
        bad.append(('sin2t', t))
    lhs, rhs = c.sq() - s.sq(), c2
    if not (lhs.lo <= rhs.hi and rhs.lo <= lhs.hi):
        bad.append(('cos2t', t))
chk('double-angle identities consistent', not bad, bad[:3])

print('H. the exact phase facts the reduction of C(u) depends on')
L = vc.L
for k in range(4):
    w = Q(2*k+1) * Q(5, 8)          # w_k / pi
    chk(f'w_{k} L / pi = {w*L} is a half-integer', (w * L * 2).denominator == 1)
for j in range(4):
    for k in range(4):
        p1 = (Q(2*j+1) - Q(2*k+1)) * Q(5, 8) * L
        p2 = (Q(2*j+1) + Q(2*k+1)) * Q(5, 8) * L
        if not (p1.denominator == 1 and p2.denominator == 1):
            fails.append(f'phase {j},{k}')
chk('every phase (w_j +- w_k)L is an INTEGER multiple of pi', 'phase 0,0' not in fails)
# and therefore sin of it is enclosed tightly around 0 by the implementation

print('I. logs and pi')
pi = vc.PI
chk('pi in (3.14159265358979323, 3.14159265358979324)',
    pi.lo > Q(314159265358979323, 10**17) and pi.hi < Q(314159265358979324, 10**17))
l2, l3 = vc.LOG2, vc.LOG3
chk('log2 in (0.69314718055994530, 0.69314718055994532)',
    l2.lo > Q(69314718055994530, 10**17) and l2.hi < Q(69314718055994532, 10**17))
e = vc.iv_exp(l2)
chk('exp(log 2) encloses 2', e.lo <= 2 <= e.hi)
e = vc.iv_exp(l3)
chk('exp(log 3) encloses 3', e.lo <= 3 <= e.hi)
s = l2 + l3
e = vc.iv_exp(s)
chk('exp(log2+log3) encloses 6', e.lo <= 6 <= e.hi)
lp = vc.LOGPI
e = vc.iv_exp(lp)
chk('exp(log pi) encloses pi', e.lo <= pi.hi and pi.lo <= e.hi)
chk('e^{8/5} < 5 and log 4 < 8/5  (prime-sum cutoff)',
    vc.iv_exp(Iv(Q(8, 5))).hi < 5 and (Iv(Q(2)) * l2).hi < Q(8, 5))

print('J. width monotonicity: a wider input never gives a narrower output')
w1 = (Iv(Q(1, 3)) * Iv(Q(1, 3))).hi - (Iv(Q(1, 3)) * Iv(Q(1, 3))).lo
t = Iv(Q(1, 3) - Q(1, 10**6), Q(1, 3) + Q(1, 10**6))
w2 = (t * t).hi - (t * t).lo
chk('widening the input widens the product', w2 > w1)

print()
if fails:
    print(f'AUDIT FAILED: {len(fails)} property/properties -> {fails}')
    sys.exit(1)
print('AUDIT PASSED: every interval-primitive property held.')
