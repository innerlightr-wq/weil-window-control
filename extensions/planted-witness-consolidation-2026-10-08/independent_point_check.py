#!/usr/bin/env python3
r"""INDEPENDENT point check of Q_add,2/5(f) < 0.   Standard library only.

SCOPE OF INDEPENDENCE (read this before calling it independent):

  * It does NOT read certificate.json, cold_certificate.json, any qzeta cache, or any
    numerical output of verify_certificate*.py.  It does not import that module.
  * It rebuilds the full geometric-side Q_zeta(f) and the exact quartet term from the
    frozen mathematical input (x, L, gamma) only.
  * The ARCHIMEDEAN TERM USES A DIFFERENT CLOSED FORM.  The first program sums
    131071 Laplace integrals and bounds the tail by [3/4 C(0) + 1/2 sup|C'|]/(M-3/4).
    This program evaluates the same series in closed form through digamma and trigamma
    at four complex points:

        Arch(f) = sum_k [ U B_k Re psi(z_k) + A_k Im psi(z_k) + (B_k/2) Re psi'(z_k) ]
                  - 2 sum_{m>=0} E_m I_m^(1) ,      z_k = 1/4 - i w_k/2 ,

    derived in ANALYTIC_AUDIT.md section 5.  The Euler constant cancels identically here,
    so no harmonic-number bound is used at all.  The error budget is completely different:
    a g-series tail for psi and a midpoint-rule bound for psi', not a 1/M tail.
  * SHARED, and documented as such: the frozen input; the classical explicit formula; the
    canonical (A_k, B_k) reduction of the autocorrelation (recoded here from the source
    C_terms, and cross-checked against C(0) = sum x_k^2 and C(2L) = 0); the idea of exact
    rational interval arithmetic (re-implemented below, not imported).
  * This is not external peer review.

TRUNCATIONS, chosen from the error budget BEFORE looking at the result:
    JREC = 50      digamma/trigamma recurrence depth   (Re z -> 50.25)
    JG   = 20000   terms of the g-series for psi        tail <= 0.5/(50+JG) = 2.49e-5
    psi' midpoint bound                                 <= 1/(12 (R-3/2)^3) = 7.2e-7
    MEXP = 24      terms of the E-series                tail <= 1e-28
  Target: total budget well below 1e-3, against an expected margin near 1.45e-2.
"""
from fractions import Fraction as Q
import json
import sys

GB = 200
GRID = 1 << GB


def fl(x):
    return Q(x.numerator * GRID // x.denominator, GRID)


def ce(x):
    return Q(-((-x.numerator) * GRID // x.denominator), GRID)


class R:
    """Real interval, exact rational endpoints, outward rounding. Re-implemented here."""
    __slots__ = ('a', 'b')

    def __init__(self, a, b=None):
        if b is None:
            b = a
        a, b = Q(a), Q(b)
        assert a <= b, f'empty interval [{a},{b}]'
        self.a, self.b = fl(a), ce(b)

    def __repr__(self):
        return f'[{float(self.a):.17g},{float(self.b):.17g}]'

    def __add__(s, o):
        o = R.mk(o)
        return R(s.a + o.a, s.b + o.b)
    __radd__ = __add__

    def __neg__(s):
        return R(-s.b, -s.a)

    def __sub__(s, o):
        return s + (-R.mk(o))

    def __rsub__(s, o):
        return R.mk(o) + (-s)

    def __mul__(s, o):
        o = R.mk(o)
        p = (s.a*o.a, s.a*o.b, s.b*o.a, s.b*o.b)
        return R(min(p), max(p))
    __rmul__ = __mul__

    def inv(s):
        assert not (s.a <= 0 <= s.b), f'inv of {s}'
        return R(Q(1)/s.b, Q(1)/s.a)

    def __truediv__(s, o):
        return s * R.mk(o).inv()

    def __rtruediv__(s, o):
        return R.mk(o) * s.inv()

    def sq(s):
        if s.a >= 0:
            return R(s.a*s.a, s.b*s.b)
        if s.b <= 0:
            return R(s.b*s.b, s.a*s.a)
        return R(0, max(s.a*s.a, s.b*s.b))

    def mag(s):
        return max(abs(s.a), abs(s.b))

    @staticmethod
    def mk(x):
        return x if isinstance(x, R) else R(x)


class C:
    """Complex interval as a pair of real intervals."""
    __slots__ = ('re', 'im')

    def __init__(self, re, im=0):
        self.re, self.im = R.mk(re), R.mk(im)

    def __repr__(self):
        return f'({self.re!r} + i{self.im!r})'

    def __add__(s, o):
        o = C.mk(o)
        return C(s.re + o.re, s.im + o.im)
    __radd__ = __add__

    def __neg__(s):
        return C(-s.re, -s.im)

    def __sub__(s, o):
        return s + (-C.mk(o))

    def __rsub__(s, o):
        return C.mk(o) + (-s)

    def __mul__(s, o):
        o = C.mk(o)
        return C(s.re*o.re - s.im*o.im, s.re*o.im + s.im*o.re)
    __rmul__ = __mul__

    def inv(s):
        d = s.re.sq() + s.im.sq()
        return C(s.re/d, -s.im/d)

    def __truediv__(s, o):
        return s * C.mk(o).inv()

    def __rtruediv__(s, o):
        return C.mk(o) * s.inv()

    def abs2(s):
        return s.re.sq() + s.im.sq()

    def widen(s, eps):
        """widen both parts by the rational radius eps (for enclosing a truncation error)"""
        e = R(-abs(eps), abs(eps))
        return C(s.re + e, s.im + e)

    @staticmethod
    def mk(x):
        return x if isinstance(x, C) else C(R.mk(x))


# ---------------------------------------------------------------- elementary functions
def r_exp(x, n=200):
    """exp of a real interval.  Range reduction t = 2^j u with |u| <= 1, then j squarings;
    negative arguments go through the reciprocal.  So the Taylor sum is ALWAYS evaluated at
    |u| <= 1 with positive terms -- no alternating cancellation and no huge rationals.
    Proved remainder |R_k| <= |u|^k/k! / (1 - |u|/(k+1))."""
    def small(u):
        a = abs(u)
        assert a <= 1
        k, f = 1, Q(1)
        while True:
            f *= k
            if k + 1 > a and a**k / f < Q(1, 1 << 240):
                break
            k += 1
            assert k <= n, 'r_exp truncation'
        s, term = Q(0), Q(1)
        for i in range(k):
            s += term
            term = term * u / (i + 1)
        rem = a**k / f / (1 - a/Q(k + 1))
        return R(s - rem, s + rem)

    def one(t):
        t = Q(t)
        neg = t < 0
        a = -t if neg else t
        j = 0
        while a > 1:
            a = a / 2
            j += 1
        v = small(a)
        for _ in range(j):
            v = v.sq()
        return v.inv() if neg else v

    lo, hi = one(x.a), one(x.b)
    return R(min(lo.a, hi.a), max(lo.b, hi.b))


def r_sincos(t, n=1200):
    """sin, cos of a RATIONAL t; |R_n| <= |t|^n/n! (Lagrange)."""
    t = Q(t)
    a = abs(t)
    k = 2
    f = Q(2)
    while a**k / f > Q(1, 1 << 240):
        k += 2
        f = f * (k - 1) * k
        assert k <= n, 'r_sincos truncation'
    rem = a**k / f
    S = Cc = Q(0)
    term = Q(1)
    for i in range(k):
        if i % 4 == 0:
            Cc += term
        elif i % 4 == 1:
            S += term
        elif i % 4 == 2:
            Cc -= term
        else:
            S -= term
        term = term * t / (i + 1)
    return R(S - rem, S + rem), R(Cc - rem, Cc + rem)


def r_sqrt(x):
    import math as _m
    assert x.a >= 0
    lo = Q(_m.isqrt(x.a.numerator * GRID * GRID // x.a.denominator), GRID)
    nn = -((-x.b.numerator) * GRID * GRID // x.b.denominator)
    r = _m.isqrt(nn)
    if r*r < nn:
        r += 1
    hi = Q(r, GRID)
    assert lo*lo <= x.a and hi*hi >= x.b, 'r_sqrt bracket'
    return R(lo, hi)


def _atanh_inv(n, K=140):
    s = Q(0)
    for k in range(K):
        s += Q(1, (2*k + 1) * n**(2*k + 1))
    rem = Q(1, (2*K + 1) * n**(2*K + 1)) / (1 - Q(1, n*n))
    return R(s, s + rem)


def _atan_inv(n, K=70):
    s = Q(0)
    for k in range(K):
        s += Q((-1)**k, (2*k + 1) * n**(2*k + 1))
    return R(s - Q(1, (2*K+1)*n**(2*K+1)), s + Q(1, (2*K+1)*n**(2*K+1)))


PI = R(16)*_atan_inv(5) - R(4)*_atan_inv(239)
LOG2 = R(2)*_atanh_inv(3)
LOG3 = LOG2 + R(2)*_atanh_inv(5)


def r_log(x, K=260):
    """log of a positive real interval, via 2^e * m and 2 atanh((m-1)/(m+1))."""
    assert x.a > 0
    def one(t):
        e = 0
        while t >= 2:
            t /= 2
            e += 1
        while t < 1:
            t *= 2
            e -= 1
        z = (t - 1)/(t + 1)
        s = Q(0)
        for k in range(K):
            s += z**(2*k + 1)/(2*k + 1)
        az = abs(z)
        rem = az**(2*K + 1)/(2*K + 1)/(1 - az*az)
        return R(e)*LOG2 + R(2*(s - rem), 2*(s + rem))
    a, b = one(x.a), one(x.b)
    return R(min(a.a, b.a), max(a.b, b.b))


LOGPI = r_log(PI)


def r_atan_small(t, K=160):
    """arctan of a real interval with |t| <= 1/2; alternating series, proved remainder."""
    assert t.mag() <= Q(1, 2)
    def one(v):
        s = Q(0)
        for k in range(K):
            s += Q((-1)**k) * v**(2*k + 1) / (2*k + 1)
        rem = abs(v)**(2*K + 1)/(2*K + 1)
        return R(s - rem, s + rem)
    a, b = one(t.a), one(t.b)
    return R(min(a.a, b.a), max(a.b, b.b))


def c_log(z):
    """principal log of a complex interval with Re > 0: (1/2)log|z|^2 + i atan(Im/Re)."""
    assert z.re.a > 0
    return C(R(Q(1, 2))*r_log(z.abs2()), r_atan_small(z.im/z.re))


# ------------------------------------------------------ digamma and trigamma, rigorously
JREC = 50
JG = 20000     # tail <= 0.5/(50+JG) = 2.50e-5, far below the expected margin
NG = 8              # terms of g(x) = sum_{n>=2} (-1)^n/(n x^n)


_GREM = {}


def _g(x, mlow):
    """g(x) = 1/x - log(1+1/x) = sum_{n>=2} (-1)^n/(n x^n) for a complex interval x.

    mlow is a POSITIVE INTEGER with mlow <= |x| (here |x| >= Re x, and Re x is known
    exactly), so the remainder bound costs only integer arithmetic:
        |sum_{n>NG}| <= sum_{n>NG} 1/(n |x|^n) <= 1/((NG+1) mlow^{NG+1}) * 1/(1 - 1/mlow).
    """
    inv = x.inv()
    p = inv*inv
    s = C(0)
    for n in range(2, NG + 1):
        s = s + C(R(Q((-1)**n, n)))*p
        p = p*inv
    rem = _GREM.get(mlow)
    if rem is None:
        rem = Q(1, (NG + 1) * mlow**(NG + 1)) * Q(mlow, mlow - 1)
        _GREM[mlow] = rem
    return s.widen(rem)


def digamma(z):
    r"""psi(z) = psi(z+JREC) - sum_{j<JREC} 1/(z+j),
        psi(w)  = log w - sum_{j>=0} g(w+j),
    the second identity being  sum_{j>=0}[1/(w+j) - log(1+1/(w+j))] = log w - psi(w),
    proved by telescoping (ANALYTIC_AUDIT.md section 5.2).  Tail after JG terms:
        |sum_{j>=JG} g(w+j)| <= (1/2)(1-1/Rw)^{-1} / (Rw + JG - 1),  Rw = Re w."""
    acc = C(0)
    for j in range(JREC):
        acc = acc + (z + C(R(j))).inv()
    w = z + C(R(JREC))
    tot = C(0)
    Rwi = JREC                      # integer lower bound on Re w (= JREC + 1/4)
    for j in range(JG):
        tot = tot + _g(w + C(R(j)), Rwi + j)
    Rw = w.re.a
    tail = Q(1, 2)/(1 - Q(1)/Q(Rw))/(Q(Rw) + JG - 1)
    psi_w = (c_log(w) - tot).widen(tail)
    return psi_w - acc


def trigamma(z):
    r"""psi'(z) = sum_{j<JREC} 1/(z+j)^2 + psi'(z+JREC),
        psi'(w) = 1/(w-1/2) + E,   |E| <= 1/(12 (Re w - 3/2)^3)
    by the midpoint rule on unit cells (ANALYTIC_AUDIT.md section 5.3)."""
    acc = C(0)
    for j in range(JREC):
        t = (z + C(R(j))).inv()
        acc = acc + t*t
    w = z + C(R(JREC))
    Rw = w.re.a
    E = Q(1, 12)/(Q(Rw) - Q(3, 2))**3
    return acc + (w - C(R(Q(1, 2)))).inv().widen(E)


# ------------------------------------------------------------------- the frozen problem
L = Q(4, 5)
U = 2*L
GAMMA = Q(14)
XNUM = [99815221908, -354082045264, 928280739889, 54384691304]
XDEN = 10**12
X = [Q(n, XDEN) for n in XNUM]
W = R(5)*PI/8                       # w_k = (2k+1) W


def wk(k):
    return R(2*k + 1)*W


def canonical():
    r"""Recoded from src/zeta_window.py C_terms / Csym_terms for sector even_d, divided by
    nrm_j nrm_k = L, and reduced using sin(beta) = 0, cos(beta) = +-1 (every beta is an
    integer multiple of pi).  Result: C(u) = sum_k [A_k sin(w_k u) + B_k (U-u) cos(w_k u)]."""
    A = [R(0)]*4
    B = [R(0)]*4
    A = [R(0) for _ in range(4)]
    B = [R(0) for _ in range(4)]
    for j in range(4):
        for k in range(4):
            xx = R(X[j]*X[k]/L)
            if j == k:
                B[k] = B[k] + xx*R(Q(1, 2))
                p2 = R(4*k + 2)*W
                A[k] = A[k] + xx/(R(2)*p2) + xx/(R(2)*p2)
            else:
                p1 = R(2*(j - k))*W
                p2 = R(2*j + 2*k + 2)*W
                cb1 = Q((-1)**((j - k) % 2))
                cb2 = Q((-1)**((j + k + 1) % 2))
                c1 = R(1)/(R(2)*p1)
                c2 = R(1)/(R(2)*p2)
                A[k] = A[k] + xx*(c1*R(cb1))
                A[j] = A[j] - xx*(c1*R(cb1))
                A[k] = A[k] - xx*(c2*R(cb2))
                A[j] = A[j] - xx*(c2*R(cb2))
    return A, B


SQL = r_sqrt(R(L))
SQ2 = r_sqrt(R(2))
SQ3 = r_sqrt(R(3))
SING, COSG = r_sincos(GAMMA*L)          # sin(56/5), cos(56/5)
MEXP = 24


def arch_via_digamma(A, B):
    r"""Arch(f) = sum_k [U B_k Re psi(z_k) + A_k Im psi(z_k) + (B_k/2) Re psi'(z_k)]
                  - 2 sum_{m>=0} E_m I_m^(1).
    z_k = 1/4 - i w_k/2.  The Euler constant cancels identically; no harmonic number and
    no 1/M tail appear."""
    tot = R(0)
    for k in range(4):
        w = wk(k)
        z = C(R(Q(1, 4)), -w/R(2))
        ps = digamma(z)
        pp = trigamma(z)
        tot = tot + R(U)*B[k]*ps.re + A[k]*ps.im + (B[k]/R(2))*pp.re
    # the E-series
    e_tot = R(0)
    for m in range(MEXP):
        s = Q(4*m + 1, 2)
        E = r_exp(R(-s*U))
        I1 = R(0)
        for k in range(4):
            w = wk(k)
            d = R(s*s) + w.sq()
            I1 = I1 + A[k]*w/d - B[k]*(R(s*s) - w.sq())/(d*d)
        e_tot = e_tot + E*I1
    # geometric tail of the E-series
    sM = Q(4*MEXP + 1, 2)
    bnd = R(0)
    for k in range(4):
        w = wk(k)
        bnd = bnd + R(A[k].mag())*w/R(sM*sM) + R(B[k].mag())/R(sM*sM)
    Etail = r_exp(R(-sM*U)).b * bnd.b / (1 - float(0))  # e^{-sU} decays by e^{-16U/5} each step
    geo = Q(1)/(1 - Q(1, 2))  # crude: ratio < 1/2 for sure since e^{-16*1.6/5} = e^{-5.12}
    tailmag = Q(r_exp(R(-sM*U)).b) * Q(bnd.b) * geo
    return tot - R(2)*e_tot + R(-tailmag, tailmag)


def pole_prime_logpi(A, B):
    """Independently recoded pole, prime and log-pi terms."""
    CH = (r_exp(R(L/2)) + r_exp(R(-L/2)))/R(2)
    t = R(0)
    for k in range(4):
        w = wk(k)
        p = R((-1)**k)*R(2)*w*CH/(w.sq() + R(Q(1, 4)))/SQL
        t = t + R(X[k])*p
    pole = R(2)*t.sq()

    def C_at(tt):
        s = R(0)
        for k in range(4):
            arg = wk(k)*tt
            sn, cs = r_sincos(arg.a) if arg.a == arg.b else (None, None)
            if sn is None:                      # interval argument: midpoint + radius
                m = (arg.a + arg.b)/2
                sn0, cs0 = r_sincos(m)
                rad = max(arg.b - m, m - arg.a)
                sn = R(max(sn0.a - rad, Q(-1)), min(sn0.b + rad, Q(1)))
                cs = R(max(cs0.a - rad, Q(-1)), min(cs0.b + rad, Q(1)))
            s = s + A[k]*sn + B[k]*(R(U) - tt)*cs
        return s
    prime = R(2)*(LOG2*(R(1)/SQ2)*C_at(LOG2) + LOG3*(R(1)/SQ3)*C_at(LOG3)
                  + LOG2*R(Q(1, 2))*C_at(R(2)*LOG2))
    C0 = sum((Q(x)*Q(x) for x in X), Q(0))
    return pole, prime, LOGPI*R(C0), C0


def quartet_complex_sinc(delta):
    r"""A_delta - i B_delta = F(gamma + i delta) = sum_k x_k [ sin((w_k-zeta)L)/(w_k-zeta)
        + sin((w_k+zeta)L)/(w_k+zeta) ] / sqrt(L),  zeta = gamma + i delta.

    A DIFFERENT algebraic route from the cosh/sinh split used by the first program:
    the complex-sinc form.  Endpoint values use
      sin((2k+1)pi/2 +- 56/5 +- i dL) = sin(X) cosh(Y) + i cos(X) sinh(Y),
      X = (2k+1)pi/2 +- 56/5,  Y = +- dL,
    and sin X = (-1)^k cos(56/5), cos X = -(-1)^k sin(+-56/5)."""
    d = Q(delta)
    dL = d*L
    chY = (r_exp(R(dL)) + r_exp(R(-dL)))/R(2)
    shY = (r_exp(R(dL)) - r_exp(R(-dL)))/R(2)
    tot = C(0)
    for k in range(4):
        w = wk(k)
        sg = R((-1)**k)
        # a = w_k - zeta : aL = (2k+1)pi/2 - 56/5 - i dL
        #   sin X = (-1)^k cos(56/5),  cos X = +(-1)^k sin(56/5)   (X = ..-56/5)
        sinaL = C(sg*COSG*chY, sg*SING*(-shY))
        den = C(w - R(GAMMA), R(-d))
        tot = tot + C(R(X[k]))*(sinaL/den)
        # a = w_k + zeta : aL = (2k+1)pi/2 + 56/5 + i dL
        #   sin X = (-1)^k cos(56/5),  cos X = -(-1)^k sin(56/5)
        sinaL2 = C(sg*COSG*chY, -sg*SING*shY)
        den2 = C(w + R(GAMMA), R(d))
        tot = tot + C(R(X[k]))*(sinaL2/den2)
    Fz = tot*C(R(1)/SQL)
    return Fz.re, -Fz.im          # A_delta, B_delta


def main():
    print('INDEPENDENT point check of Q_add,2/5(f) < 0')
    print('  truncations fixed in advance: JREC=%d  JG=%d  NG=%d  MEXP=%d  grid=2^-%d'
          % (JREC, JG, NG, MEXP, GB))
    A, B = canonical()
    C0 = sum((Q(x)*Q(x) for x in X), Q(0))
    # structural self-checks of the recoded reduction
    c0 = R(0)
    for k in range(4):
        c0 = c0 + B[k]*R(U)
    assert c0.a <= C0 <= c0.b, 'C(0) != sum x_k^2'
    print(f'  C(0) from A,B  = {c0!r}   vs sum x_k^2 = {float(C0):.18f}   OK')
    cU = R(0)
    for k in range(4):
        sn, cs = r_sincos(Q(0))
        arg = wk(k)*R(U)
        m = (arg.a + arg.b)/2
        sn0, _ = r_sincos(m)
        rad = max(arg.b - m, m - arg.a)
        cU = cU + A[k]*R(sn0.a - rad, sn0.b + rad)
    print(f'  C(2L) trig part = {cU!r}   ((U-u) factor vanishes at u=2L identically)')

    pole, prime, logpi, _ = pole_prime_logpi(A, B)
    ar = arch_via_digamma(A, B)
    Qz = pole + ar - logpi - prime
    print(f'  pole   = {pole!r}')
    print(f'  arch   = {ar!r}        <- digamma/trigamma closed form')
    print(f'  logpi  = {logpi!r}')
    print(f'  prime  = {prime!r}')
    print(f'  Q_zeta = {Qz!r}   width {float(Qz.b - Qz.a):.3e}')

    Ad, Bd = quartet_complex_sinc(Q(2, 5))
    Qa = Qz + R(4)*Ad.sq() - R(4)*Bd.sq()
    print(f'  A_2/5  = {Ad!r}   4A^2 upper = {float((R(4)*Ad.sq()).b):.9e}')
    print(f'  B_2/5  = {Bd!r}   4B^2 lower = {float((R(4)*Bd.sq()).a):.9e}')
    print(f'  Q_add,2/5 = {Qa!r}')
    neg = Qa.b < 0
    print(f'\n  END-TO-END NEGATIVE UPPER BOUND: {neg}   '
          f'upper = {float(Qa.b):.12f}   eta = {float(-Qa.b):.12f}')
    out = dict(scope_of_independence='see module docstring',
               truncations=dict(JREC=JREC, JG=JG, NG=NG, MEXP=MEXP, gridbits=GB),
               Q_zeta=dict(lo=str(Qz.a), hi=str(Qz.b)),
               pole=dict(lo=str(pole.a), hi=str(pole.b)),
               arch=dict(lo=str(ar.a), hi=str(ar.b)),
               logpi=dict(lo=str(logpi.a), hi=str(logpi.b)),
               prime=dict(lo=str(prime.a), hi=str(prime.b)),
               A_delta=dict(lo=str(Ad.a), hi=str(Ad.b)),
               B_delta=dict(lo=str(Bd.a), hi=str(Bd.b)),
               Q_add_2_5=dict(lo=str(Qa.a), hi=str(Qa.b)),
               negative_upper_bound=bool(neg), eta=str(-Qa.b), eta_f=float(-Qa.b))
    json.dump(out, open('independent_point_check.json', 'w'), indent=1)
    if not neg:
        raise SystemExit('INDEPENDENT CHECK DID NOT ESTABLISH A NEGATIVE UPPER BOUND')
    print('  wrote independent_point_check.json')


if __name__ == '__main__':
    main()
