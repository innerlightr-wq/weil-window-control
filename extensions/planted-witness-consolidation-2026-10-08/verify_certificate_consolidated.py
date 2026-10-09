#!/usr/bin/env python3
r"""Rigorous verifier for the added-quartet witness certificate.  Standard library only.

ARITHMETIC MODEL.  Every quantity is an interval with exact Fraction endpoints, rounded
OUTWARD onto a dyadic grid (denominator 2^GRIDBITS) after each operation.  There are no
float operations and no libm calls anywhere in the certified path: pi, log 2, log 3, exp,
sin, cos, sqrt are all produced by truncated series whose remainders are bounded by proved
inequalities stated at each function.  Floats appear only in the optional human-readable
printout and in the regression comparisons, which are clearly labelled and carry no weight.

The script FAILS (raises) rather than relaxing any bound.

Object certified (see CERTIFICATE_THEOREM.md):
    Q_add,delta(f) = Q_zeta(f) + 4 [ A_delta(f)^2 - B_delta(f)^2 ]
    A_delta(f) = int_{-L}^{L} f(u) cos(14 u) cosh(delta u) du
    B_delta(f) = int_{-L}^{L} f(u) sin(14 u) sinh(delta u) du
with L = 4/5, gamma = 14 exactly, f an explicit element of the 4-mode EVEN_D span.
"""
from fractions import Fraction as Q
import json
import sys

GRIDBITS = 260
GRID = 1 << GRIDBITS


class CertificationFailure(Exception):
    pass


def _fl(x):
    """round Fraction x DOWN onto the grid"""
    return Q(x.numerator * GRID // x.denominator, GRID)


def _ce(x):
    """round Fraction x UP onto the grid"""
    return Q(-((-x.numerator) * GRID // x.denominator), GRID)


class Iv:
    """A closed interval [lo, hi] of rationals; all operations round outward."""
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        if hi is None:
            hi = lo
        lo, hi = Q(lo), Q(hi)
        if lo > hi:
            raise CertificationFailure(f'empty interval [{lo}, {hi}]')
        self.lo, self.hi = _fl(lo), _ce(hi)

    # ---- basics
    def __repr__(self):
        return f'[{float(self.lo):.17g}, {float(self.hi):.17g}]'

    def __neg__(self):
        return Iv(-self.hi, -self.lo)

    def __add__(self, o):
        o = mk(o)
        return Iv(self.lo + o.lo, self.hi + o.hi)

    __radd__ = __add__

    def __sub__(self, o):
        return self + (-mk(o))

    def __rsub__(self, o):
        return mk(o) + (-self)

    def __mul__(self, o):
        o = mk(o)
        p = (self.lo*o.lo, self.lo*o.hi, self.hi*o.lo, self.hi*o.hi)
        return Iv(min(p), max(p))

    __rmul__ = __mul__

    def inv(self):
        if self.lo <= 0 <= self.hi:
            raise CertificationFailure(f'division by an interval containing 0: {self}')
        return Iv(Q(1)/self.hi, Q(1)/self.lo)

    def __truediv__(self, o):
        return self * mk(o).inv()

    def __rtruediv__(self, o):
        return mk(o) * self.inv()

    def __pow__(self, n):
        if n == 0:
            return Iv(1)
        if n < 0:
            return (self ** (-n)).inv()
        r = Iv(1)
        for _ in range(n):
            r = r * self
        return r

    def sq(self):
        """tighter than self*self when the interval straddles 0"""
        a, b = self.lo, self.hi
        if a >= 0:
            return Iv(a*a, b*b)
        if b <= 0:
            return Iv(b*b, a*a)
        return Iv(0, max(a*a, b*b))

    def contains_zero(self):
        return self.lo <= 0 <= self.hi

    def width(self):
        return self.hi - self.lo

    def mid(self):
        return (self.lo + self.hi)/2


def mk(x):
    return x if isinstance(x, Iv) else Iv(x)


# ----------------------------------------------------------------- proved transcendentals
def iv_sqrt(x):
    """[r, s] with r^2 <= x.lo and s^2 >= x.hi.  Integer-sqrt based, hence exact."""
    if x.lo < 0:
        raise CertificationFailure(f'sqrt of an interval with negative part: {x}')
    import math as _m
    def lower(a):           # largest grid rational r with r^2 <= a
        n = a.numerator * GRID * GRID // a.denominator
        return Q(_m.isqrt(n), GRID)
    def upper(a):
        n = -((-a.numerator) * GRID * GRID // a.denominator)
        r = _m.isqrt(n)
        if r*r < n:
            r += 1
        return Q(r, GRID)
    lo, hi = lower(x.lo), upper(x.hi)
    if lo*lo > x.lo or hi*hi < x.hi:
        raise CertificationFailure('iv_sqrt bracket failed')
    return Iv(lo, hi)


def iv_exp(x, nmax=400):
    r"""exp of an interval.  For a rational t, exp(t) = sum_{k<n} t^k/k! + R_n with
    |R_n| <= |t|^n/n! * 1/(1-|t|/(n+1))  whenever n+1 > |t|  (ratio test, proved).
    Evaluated at both endpoints; exp is increasing, so the hull is the enclosure."""
    def one(t):
        t = Q(t)
        a = abs(t)
        n = 1
        while n <= nmax and not (Q(n + 1) > a and a**n / _fact(n) < Q(1, 1 << 300)):
            n += 1
        if n > nmax:
            raise CertificationFailure(f'iv_exp: no usable truncation for t={float(t)}')
        s = Q(0)
        term = Q(1)
        for k in range(n):
            s += term
            term = term * t / (k + 1)
        rem = a**n / _fact(n) / (1 - a/Q(n + 1))
        return Iv(s - rem, s + rem)
    e_lo, e_hi = one(x.lo), one(x.hi)
    return Iv(min(e_lo.lo, e_hi.lo), max(e_lo.hi, e_hi.hi))


_FACTC = {0: Q(1)}


def _fact(n):
    if n not in _FACTC:
        _FACTC[n] = _fact(n - 1) * n
    return _FACTC[n]


def _sincos_rat(t, nmax=1200):
    r"""sin(t), cos(t) for a RATIONAL t, by Taylor with the proved remainder
    |R_n| <= |t|^n/n!  (valid for sin/cos for every n, since the n-th derivative is
    bounded by 1 and Taylor's theorem with Lagrange remainder applies)."""
    t = Q(t)
    a = abs(t)
    n = 2
    while n <= nmax and a**n / _fact(n) > Q(1, 1 << 300):
        n += 2
    if n > nmax:
        raise CertificationFailure(f'_sincos_rat: no usable truncation for t={float(t)}')
    rem = a**n / _fact(n)
    s = Q(0)
    c = Q(0)
    term = Q(1)                      # t^k/k!
    for k in range(n):
        if k % 4 == 0:
            c += term
        elif k % 4 == 1:
            s += term
        elif k % 4 == 2:
            c -= term
        else:
            s -= term
        term = term * t / (k + 1)
    return Iv(s - rem, s + rem), Iv(c - rem, c + rem)


def iv_sin(x):
    """sin of an interval: enclose via sin at the midpoint plus |x-mid| (|sin'|<=1)."""
    m = x.mid()
    s, _ = _sincos_rat(m)
    r = max(x.hi - m, m - x.lo)
    out = Iv(max(s.lo - r, Q(-1)), min(s.hi + r, Q(1)))
    return out


def iv_cos(x):
    m = x.mid()
    _, c = _sincos_rat(m)
    r = max(x.hi - m, m - x.lo)
    return Iv(max(c.lo - r, Q(-1)), min(c.hi + r, Q(1)))


def iv_sinh(x):
    e = iv_exp(x)
    return (e - e.inv()) * Q(1, 2)


def iv_cosh(x):
    e = iv_exp(x)
    return (e + e.inv()) * Q(1, 2)


def _atan_recip(n, K=60):
    r"""arctan(1/n) for integer n>=2.  Alternating series with decreasing terms, so the
    truncation error is bounded by the first omitted term: proved."""
    s = Q(0)
    for k in range(K):
        s += Q((-1)**k, (2*k + 1) * n**(2*k + 1))
    rem = Q(1, (2*K + 1) * n**(2*K + 1))
    return Iv(s - rem, s + rem)


def _atanh_recip(n, K=120):
    r"""atanh(1/n) = sum_{k>=0} 1/((2k+1) n^{2k+1}), all terms positive.
    Remainder after K terms <= 1/((2K+1) n^{2K+1}) * 1/(1 - 1/n^2): proved by comparison
    with a geometric series."""
    s = Q(0)
    for k in range(K):
        s += Q(1, (2*k + 1) * n**(2*k + 1))
    rem = Q(1, (2*K + 1) * n**(2*K + 1)) / (1 - Q(1, n*n))
    return Iv(s, s + rem)


PI = Iv(16) * _atan_recip(5) - Iv(4) * _atan_recip(239)      # Machin
LOG2 = Iv(2) * _atanh_recip(3)                               # ln2 = 2 atanh(1/3)
LOG3 = LOG2 + Iv(2) * _atanh_recip(5)                        # ln(3/2) = 2 atanh(1/5)


def iv_log(x):
    r"""log of a positive interval, via y = 2^e * m with m in [1,2) and
    log m = 2 atanh((m-1)/(m+1)),  |(m-1)/(m+1)| <= 1/3."""
    if x.lo <= 0:
        raise CertificationFailure(f'log of a non-positive interval: {x}')
    def one(a):
        e = 0
        while a >= 2:
            a /= 2
            e += 1
        while a < 1:
            a *= 2
            e -= 1
        z = (a - 1) / (a + 1)
        s = Q(0)
        K = 200
        for k in range(K):
            s += z**(2*k + 1) / (2*k + 1)
        az = abs(z)
        rem = az**(2*K + 1) / (2*K + 1) / (1 - az*az)
        lm = Iv(2*(s - rem), 2*(s + rem))
        return Iv(e) * LOG2 + lm
    a, b = one(x.lo), one(x.hi)
    return Iv(min(a.lo, b.lo), max(a.hi, b.hi))


LOGPI = iv_log(PI)


# =====================================================================================
#  The witness, the form, and the certificate
# =====================================================================================
L = Q(4, 5)
U = 2 * L                                  # = 8/5, the support of the autocorrelation
GAMMA = Q(14)                              # EXACT chosen ordinate
NDIM = 4

# The frozen EXACT RATIONAL coefficient vector: the 12-decimal printout of the previous
# run's N=4 nodal vector, read as exact rationals (documented choice, see
# CERTIFICATE_THEOREM.md section 2).  NOT normalised: positive scaling cannot change a sign.
XNUM = [99815221908, -354082045264, 928280739889, 54384691304]
XDEN = 10**12
X = [Q(n, XDEN) for n in XNUM]

W = Iv(5, 5) * PI / 8                      # w_k = (2k+1) * W,  W = (5/8) pi
def wk(k):
    return Iv(2*k + 1) * W


def basis_checks():
    """Endpoint / admissibility conditions, read from src/zeta_window.py."""
    out = {}
    # w_k * L = (2k+1) pi/2  =>  cos(w_k L) = 0 exactly, so phi_k(+-L) = 0: the EVEN_D
    # modes are continuous on R when extended by zero.  Verified symbolically, not
    # numerically: w_k L = (2k+1)pi/2 by construction.
    out['phi_vanishes_at_pm_L'] = 'symbolic: w_k L = (2k+1)pi/2 exactly'
    # orthonormality: int phi_j phi_k = delta_jk  (norm2 = L for even_d)
    out['orthonormal'] = 'symbolic: int_{-L}^{L} cos(w_j x)cos(w_k x) dx = L delta_jk'
    # log 5 >= 2L, so n = 5 is outside the support: check e^{8/5} < 5
    e16 = iv_exp(Iv(U))
    if not e16.hi < 5:
        raise CertificationFailure('support cutoff: could not prove log 5 > 2L')
    out['prime_support_cutoff'] = (f'e^(2L) <= {float(e16.hi):.6f} < 5, so log 5 > 2L: '
                                   'prime powers are exactly n = 2,3,4')
    return out


def C_canonical():
    r"""Reduce C(u) = sum_{j,k} x_j x_k C_jk(u) to the canonical form

        C(u) = sum_k [ A_k sin(w_k u) + B_k (U - u) cos(w_k u) ] ,   0 <= u <= U.

    C_jk is src/zeta_window.py::Csym_terms for sector even_d, divided by nrm_j nrm_k = L.
    Every beta occurring is an INTEGER multiple of pi, so sin(beta) = 0 exactly and
    cos(beta) = +-1 exactly; that is what makes this reduction exact.
    """
    A = [Iv(0) for _ in range(NDIM)]
    B = [Iv(0) for _ in range(NDIM)]
    for j in range(NDIM):
        for k in range(NDIM):
            xx = Iv(X[j] * X[k] / L)        # includes the 1/(nrm_j nrm_k) = 1/L
            if j == k:
                # (1/2)(U-u)cos(w_k u)   +   p2-part with p2 = (4k+2)W, cos(p2 L) = -1
                B[k] = B[k] + xx * Q(1, 2)
                p2 = Iv(4*k + 2) * W
                # (1/(2p2)) * (-1) * sin(-w_k u)  +  (-1/(2p2)) * (-1) * sin(w_k u)
                A[k] = A[k] + xx * (Iv(1) / (Iv(2) * p2))      # from sin(-w_k u) term
                A[k] = A[k] + xx * (Iv(1) / (Iv(2) * p2))      # from sin(w_j u) term
            else:
                p1 = Iv(2*(j - k)) * W
                p2 = Iv(2*j + 2*k + 2) * W
                cb1 = Q((-1) ** ((j - k) % 2))                  # cos(p1 L) = (-1)^{j-k}
                cb2 = Q((-1) ** ((j + k + 1) % 2))              # cos(p2 L)
                c1 = Iv(1) / (Iv(2) * p1)
                c2 = Iv(1) / (Iv(2) * p2)
                A[k] = A[k] + xx * (c1 * cb1)                   # + (1/2p1) sin(w_k u)
                A[j] = A[j] - xx * (c1 * cb1)                   # - (1/2p1) sin(w_j u)
                A[k] = A[k] - xx * (c2 * cb2)                   # + (1/2p2) sin(-w_k u)
                A[j] = A[j] - xx * (c2 * cb2)                   # - (1/2p2) sin(w_j u)
    return A, B


def C_at(A, B, t):
    """C(t) for an interval t in [0, U]."""
    tot = Iv(0)
    for k in range(NDIM):
        tot = tot + A[k] * iv_sin(wk(k) * t) + B[k] * (Iv(U) - t) * iv_cos(wk(k) * t)
    return tot


def dC_sup_bound(A, B):
    """An upper bound on sup_{[0,U]} |C'(u)|.
    C' = sum_k [ A_k w_k cos(w_k u) - B_k cos(w_k u) - B_k (U-u) w_k sin(w_k u) ]."""
    tot = Iv(0)
    for k in range(NDIM):
        aw = A[k] * wk(k)
        tot = tot + Iv(max(abs(aw.lo), abs(aw.hi)))
        tot = tot + Iv(max(abs(B[k].lo), abs(B[k].hi)))
        bw = B[k] * wk(k) * Iv(U)
        tot = tot + Iv(max(abs(bw.lo), abs(bw.hi)))
    return Iv(0, tot.hi)


def laplace(A, B, s, E):
    """int_0^U e^{-s u} C(u) du, EXACT closed form.  s rational, E encloses e^{-sU}.
    Uses sin(w_k U) = 0 and cos(w_k U) = -1 exactly (w_k U = (2k+1) pi)."""
    s = Iv(s)
    one_pE = Iv(1) + E
    tot = Iv(0)
    for k in range(NDIM):
        w = wk(k)
        d = s.sq() + w.sq()
        tot = tot + A[k] * (w * one_pE / d)
        tot = tot + B[k] * (s * Iv(U) / d - one_pE * (s.sq() - w.sq()) / (d * d))
    return tot


MEXP = 1 << 17                 # M + 1 = 2^17, so ln(M+1) = 17 log 2 exactly
M = MEXP - 1


def arch(A, B, C0, progress=False):
    r"""Arch(f) = C0 (ln(M+1) - T) - 2 sum_{m<M} I_m + tail,  with

      H_M - euler = ln(M+1) - T_{M+1},   T_{M+1} in [1/(2(M+1)) - 1/(6M^2), 1/(2M)]
    proved from  x^2/2 - x^3/3 <= x - ln(1+x) <= x^2/2  (x = 1/k > 0), and

      |sum_{m>=M} t_m| <= [ (3/4)C0 + (1/2) sup|C'| ] / (M - 3/4)
    proved by integration by parts (C(U) = 0) plus integral bounds; see
    CERTIFICATE_THEOREM.md section 4.
    """
    T = Iv(Q(1, 2*(M + 1)) - Q(1, 6*M*M), Q(1, 2*M))
    acc = Iv(0)
    small = Iv(0, Q(1, 1 << 90))         # e^{-s_m U} <= 2^-90 for m >= 20 (s_m U >= 64.8)
    for m in range(M):
        s = Q(4*m + 1, 2)                # s_m = 2m + 1/2
        E = iv_exp(Iv(-s * U)) if m < 20 else small
        acc = acc + laplace(A, B, s, E)
        if progress and m % 20000 == 0:
            print(f'    arch m={m}/{M}', file=sys.stderr)
    tailc = Iv(Q(3, 4)) * Iv(C0) + Iv(Q(1, 2)) * dC_sup_bound(A, B)
    tail_mag = tailc / Iv(Q(M) - Q(3, 4))
    tail = Iv(-tail_mag.hi, tail_mag.hi)
    return Iv(C0) * (Iv(17) * LOG2 - T) - Iv(2) * acc + tail, tail


SQRT_L = iv_sqrt(Iv(L))
SQRT2 = iv_sqrt(Iv(2))
SQRT3 = iv_sqrt(Iv(3))
COSH_HALF_L = iv_cosh(Iv(L / 2))
SIN_G, COS_G = _sincos_rat(GAMMA * L)         # sin(56/5), cos(56/5): gamma*L = 56/5 exactly


def pole_term():
    """2 (sum_k x_k p_k)^2,  p_k = (-1)^k 2 w_k cosh(L/2)/(w_k^2+1/4)/sqrt(L)."""
    tot = Iv(0)
    for k in range(NDIM):
        w = wk(k)
        p = Iv((-1) ** k) * Iv(2) * w * COSH_HALF_L / (w.sq() + Iv(Q(1, 4))) / SQRT_L
        tot = tot + Iv(X[k]) * p
    return Iv(2) * tot.sq()


def prime_term(A, B):
    """2 sum_{n=2,3,4} Lambda(n) n^{-1/2} C(log n).  Support: log n < 2L only for n=2,3,4."""
    items = [(LOG2, Iv(1) / SQRT2, LOG2),          # n=2: Lambda=log2, 1/sqrt2, t=log2
             (LOG3, Iv(1) / SQRT3, LOG3),          # n=3
             (LOG2, Iv(Q(1, 2)), Iv(2) * LOG2)]    # n=4: Lambda=log2, 1/sqrt4 = 1/2
    tot = Iv(0)
    for lam, inv_sqrt_n, t in items:
        tot = tot + lam * inv_sqrt_n * C_at(A, B, t)
    return Iv(2) * tot


def F_at_gamma():
    """F(14) = (2 cos(56/5)/sqrt(L)) * sum_k x_k (-1)^k w_k/(w_k^2 - 196)."""
    tot = Iv(0)
    for k in range(NDIM):
        w = wk(k)
        tot = tot + Iv(X[k] * Q((-1) ** k)) * w / (w.sq() - Iv(GAMMA * GAMMA))
    return Iv(2) * COS_G / SQRT_L * tot


def Q_zeta(progress=False):
    A, B = C_canonical()
    C0 = sum(x * x for x in X)
    ar, tail = arch(A, B, C0, progress)
    pole = pole_term()
    prime = prime_term(A, B)
    logpi = LOGPI * Iv(C0)
    return dict(Q=pole + ar - logpi - prime, pole=pole, arch=ar, logpi=logpi,
                prime=prime, arch_tail=tail, C0=C0, A=A, B=B,
                dCsup=dC_sup_bound(A, B))


def quartet(delta):
    r"""A_delta(f) and B_delta(f), exact closed forms.

    A_delta = (1/sqrt L) sum_k x_k (1/2)[ Ic(w_k+g) + Ic(w_k-g) ],
      Ic(a) = 2[ a sin(aL) cosh(dL) + d cos(aL) sinh(dL) ]/(a^2+d^2)
    B_delta = (1/sqrt L) sum_k x_k (1/2)[ Is(g+w_k) + Is(g-w_k) ],
      Is(a) = 2[ d sin(aL) cosh(dL) - a cos(aL) sinh(dL) ]/(a^2+d^2)

    With aL = (2k+1)pi/2 +- 56/5 the endpoint trigs are EXACT multiples of
    sin(56/5), cos(56/5):
      a = w_k+g : sin(aL)=(-1)^k cos(56/5),  cos(aL)=-(-1)^k sin(56/5)
      a = w_k-g : sin(aL)=(-1)^k cos(56/5),  cos(aL)=+(-1)^k sin(56/5)
      a = g+w_k : sin(aL)=(-1)^k cos(56/5),  cos(aL)=-(-1)^k sin(56/5)
      a = g-w_k : sin(aL)=-(-1)^k cos(56/5), cos(aL)=+(-1)^k sin(56/5)
    """
    d = mk(delta)
    dL = d * Iv(L)
    ch, sh = iv_cosh(dL), iv_sinh(dL)
    d2 = d.sq()
    Atot = Iv(0)
    Btot = Iv(0)
    for k in range(NDIM):
        w = wk(k)
        sgn = Iv((-1) ** k)
        xk = Iv(X[k])
        # ---- A: a = w_k + g and a = w_k - g
        for a, cosaL in ((w + Iv(GAMMA), -sgn * SIN_G), (w - Iv(GAMMA), sgn * SIN_G)):
            Ic = Iv(2) * (a * (sgn * COS_G) * ch + d * cosaL * sh) / (a.sq() + d2)
            Atot = Atot + xk * Q(1, 2) * Ic
        # ---- B: a = g + w_k and a = g - w_k
        for a, sinaL, cosaL in ((Iv(GAMMA) + w, sgn * COS_G, -sgn * SIN_G),
                                (Iv(GAMMA) - w, -sgn * COS_G, sgn * SIN_G)):
            Is = Iv(2) * (d * sinaL * ch - a * cosaL * sh) / (a.sq() + d2)
            Btot = Btot + xk * Q(1, 2) * Is
    return Atot / SQRT_L, Btot / SQRT_L


def Q_add(Qz, delta):
    """Q_add,delta(f) = Q_zeta(f) + 4[A_delta^2 - B_delta^2].  Returns an enclosure."""
    Ad, Bd = quartet(delta)
    return Qz + Iv(4) * Ad.sq() - Iv(4) * Bd.sq(), Ad, Bd


def Fp_at_gamma():
    r"""F'(14) = -(1/sqrt L) sum_k x_k (1/2)[ J(g+w_k) + J(g-w_k) ],
    J(a) = 2[ -L cos(aL)/a + sin(aL)/a^2 ],  endpoint trigs exact as in quartet()."""
    tot = Iv(0)
    for k in range(NDIM):
        w = wk(k)
        sgn = Iv((-1) ** k)
        for a, sinaL, cosaL in ((Iv(GAMMA) + w, sgn * COS_G, -sgn * SIN_G),
                                (Iv(GAMMA) - w, -sgn * COS_G, sgn * SIN_G)):
            J = Iv(2) * (-Iv(L) * cosaL / a + sinaL / a.sq())
            tot = tot + Iv(X[k]) * Q(1, 2) * J
    return -tot / SQRT_L


def Fpp_at_gamma():
    r"""F''(14) = -(1/sqrt L) sum_k x_k (1/2)[ G(g+w_k) + G(g-w_k) ],
    G(a) = int_{-L}^{L} u^2 cos(au) du = 2[ L^2 sin(aL)/a + 2L cos(aL)/a^2
                                            - 2 sin(aL)/a^3 ]."""
    tot = Iv(0)
    for k in range(NDIM):
        w = wk(k)
        sgn = Iv((-1) ** k)
        for a, sinaL, cosaL in ((Iv(GAMMA) + w, sgn * COS_G, -sgn * SIN_G),
                                (Iv(GAMMA) - w, -sgn * COS_G, sgn * SIN_G)):
            G = Iv(2) * (Iv(L * L) * sinaL / a + Iv(2 * L) * cosaL / a.sq()
                         - Iv(2) * sinaL / (a.sq() * a))
            tot = tot + Iv(X[k]) * Q(1, 2) * G
    return -tot / SQRT_L


def PEr_norm2(gamma_val):
    """||P_E r_gamma||^2 = sum_k <phi_k, r_gamma>^2 = sum_k F_k'(gamma)^2, closed form.
    Used only for the high-frequency remark (correction C); gamma_val must be rational."""
    g = Q(gamma_val)
    sg, cg = _sincos_rat(g * L)
    tot = Iv(0)
    for k in range(NDIM):
        w = wk(k)
        sgn = Iv((-1) ** k)
        s = Iv(0)
        for a, sinaL, cosaL in ((Iv(g) + w, sgn * cg, -sgn * sg),
                                (Iv(g) - w, -sgn * cg, sgn * sg)):
            J = Iv(2) * (-Iv(L) * cosaL / a + sinaL / a.sq())
            s = s + Q(1, 2) * J
        fk = -s / SQRT_L
        tot = tot + fk.sq()
    return tot


# =====================================================================================
#  main
# =====================================================================================
def _d(iv):
    return dict(lo=str(iv.lo), hi=str(iv.hi),
                lo_f=float(iv.lo), hi_f=float(iv.hi), width_f=float(iv.width()))


def main():
    out = dict(model=dict(
        arithmetic='exact Fraction intervals, outward rounding on a 2^-%d dyadic grid' % GRIDBITS,
        transcendentals='pi (Machin), log2/log3 (atanh series), exp/sin/cos/sqrt '
                        '(truncated series with proved remainders); NO libm in the '
                        'certified path',
        L='4/5', gamma='14 (exact)', dim=NDIM,
        x_numerators=XNUM, x_denominator=XDEN,
        arch_truncation_M=M, arch_truncation_note='M+1 = 2^17 so ln(M+1) = 17 log 2'))
    out['basis_checks'] = basis_checks()

    print('[1] Q_zeta(f): certified enclosure  (this is the slow step, ~5 min)',
          file=sys.stderr)
    import os
    key = dict(x=XNUM, xd=XDEN, M=M, grid=GRIDBITS, L=str(L),
               gamma=str(GAMMA), ndim=NDIM, src='even_d cos((2k+1)pi x/2L)/sqrt(L)')
    cache = 'regenerated_qzeta_cache.json'
    r = None
    if os.path.exists(cache):
        cj = json.load(open(cache))
        if cj.get('key') == key:
            A, B = C_canonical()
            r = dict(Q=Iv(Q(cj['Q']['lo']), Q(cj['Q']['hi'])),
                     pole=Iv(Q(cj['pole']['lo']), Q(cj['pole']['hi'])),
                     arch=Iv(Q(cj['arch']['lo']), Q(cj['arch']['hi'])),
                     logpi=Iv(Q(cj['logpi']['lo']), Q(cj['logpi']['hi'])),
                     prime=Iv(Q(cj['prime']['lo']), Q(cj['prime']['hi'])),
                     arch_tail=Iv(Q(cj['arch_tail']['lo']), Q(cj['arch_tail']['hi'])),
                     C0=Q(cj['C0']), A=A, B=B, dCsup=dC_sup_bound(A, B))
            print('    (loaded from qzeta_cache.json; delete it to force a recompute)',
                  file=sys.stderr)
    if r is None:
        r = Q_zeta()
        import platform
        import time as _t
        json.dump(dict(key=key,
                       provenance=dict(
                           produced_by='verify_certificate_consolidated.py (cold path)',
                           python=sys.version, platform=platform.platform(),
                           utc=_t.strftime('%Y-%m-%dT%H:%M:%SZ', _t.gmtime()),
                           note='A hash authenticates file identity, not mathematical '
                                'correctness. The cold computation is the validation.'),
                       **{k: dict(lo=str(r[k].lo), hi=str(r[k].hi))
                          for k in ('Q', 'pole', 'arch', 'logpi', 'prime', 'arch_tail')},
                       C0=str(r['C0'])), open(cache, 'w'), indent=1)
    Qz = r['Q']
    out['Q_zeta'] = _d(Qz)
    out['Q_zeta_parts'] = {k: _d(r[k]) for k in ('pole', 'arch', 'logpi', 'prime')}
    out['arch_tail_allowance'] = _d(r['arch_tail'])
    out['sup_dC_bound'] = float(r['dCsup'].hi)
    out['C0_norm2'] = str(r['C0'])
    if Qz.lo <= 0:
        print('NOTE: the certified Q_zeta enclosure reaches 0 or below.', file=sys.stderr)

    print('[2] node residual and the two functionals', file=sys.stderr)
    Fg = F_at_gamma()
    Fp = Fp_at_gamma()
    Fpp = Fpp_at_gamma()
    a_f = Iv(4) * (Fp.sq() + Fg * Fpp)
    diff = Iv(2) * Fg.sq()                 # Q_add - Q_replacement = +2F(14)^2
    out['F_gamma'] = _d(Fg)
    out['Fp_gamma'] = _d(Fp)
    out['Fpp_gamma'] = _d(Fpp)
    out['a_f_4_Fp2_plus_F_Fpp'] = _d(a_f)
    out['Q_add_minus_Q_replacement'] = _d(diff)

    print('[3] point certificate at delta = 2/5', file=sys.stderr)
    Qa, Ad, Bd = Q_add(Qz, Q(2, 5))
    out['point_delta_2_5'] = dict(delta='2/5', A_delta=_d(Ad), B_delta=_d(Bd),
                                  four_A2=_d(Iv(4)*Ad.sq()), four_B2=_d(Iv(4)*Bd.sq()),
                                  Q_add=_d(Qa), negative=bool(Qa.hi < 0),
                                  margin_eta=str(-Qa.hi), margin_eta_f=float(-Qa.hi))
    print(f'    Q_add,2/5 = {Qa!r}   negative: {Qa.hi < 0}', file=sys.stderr)

    print('[4] interval certificate on delta in [1/4, 49/100] by subdivision',
          file=sys.stderr)
    DLO, DHI = Q(1, 4), Q(49, 100)
    single, Ai, Bi = Q_add(Qz, Iv(DLO, DHI))
    out['interval_single_box_attempt'] = dict(
        Q_add_enclosure=_d(single), negative_everywhere=bool(single.hi < 0),
        note='ONE interval evaluation over the whole delta-box. Valid but far too loose: '
             'delta occurs in several places in the closed forms (d, d^2, cosh(dL), '
             'sinh(dL)), so independent interval occurrences inflate the result. '
             'Recorded as the honest first attempt, not as a result.')
    sub = []
    worst = None
    K = 256
    cells = []
    for j in range(K):
        a = Q(800 + 3*j, 3200)
        b = Q(803 + 3*j, 3200)
        qi, Aj, Bj = Q_add(Qz, Iv(a, b))
        cells.append(dict(j=j, lo=str(a), hi=str(b),
                          Q_upper=str(qi.hi), Q_lower=str(qi.lo),
                          four_A2_upper=str((Iv(4)*Aj.sq()).hi),
                          four_B2_lower=str((Iv(4)*Bj.sq()).lo),
                          negative=bool(qi.hi < 0)))
        if worst is None or qi.hi > worst[1]:
            worst = ((a, b), qi.hi)
        if qi.hi >= 0:
            sub.append(dict(lo=str(a), hi=str(b), Q_upper=str(qi.hi), negative=False))
    # coverage proof: consecutive closed cells ABUT -- hi(I_j) = lo(I_{j+1}) exactly,
    # overlap 0, no gaps
    gaps = [j for j in range(K-1) if Q(803 + 3*j, 3200) < Q(800 + 3*(j+1), 3200)]
    if gaps:
        raise CertificationFailure(f'subdivision has gaps at cells {gaps[:5]}')
    if Q(800, 3200) != DLO or Q(803 + 3*(K-1), 3200) != DHI:
        raise CertificationFailure('subdivision endpoints do not match [1/4, 49/100]')
    json.dump(dict(delta_lo=str(DLO), delta_hi=str(DHI), n_cells=K,
                   cell_formula='I_j = [(800+3j)/3200, (803+3j)/3200], j=0..255',
                   coverage='verified exactly: hi(I_j) = (803+3j)/3200 = '
                            'lo(I_{j+1}) = (800+3(j+1))/3200, so consecutive CLOSED cells '
                            'abut (shared endpoint, zero gap); lo(I_0) = 1/4 and '
                            'hi(I_255) = 49/100, so the union is exactly [1/4, 49/100]',
                   all_negative=bool(not sub), cells=cells),
              open('interval_cells.json', 'w'), indent=1)
    out['interval_1_4__49_100'] = dict(
        delta_lo='1/4', delta_hi='49/100', subdivisions=K,
        all_negative=bool(not sub), failures=sub,
        worst_subinterval=[str(worst[0][0]), str(worst[0][1])],
        worst_Q_upper=str(worst[1]), worst_Q_upper_f=float(worst[1]),
        margin_eta=str(-worst[1]), margin_eta_f=float(-worst[1]),
        method=f'uniform subdivision into {K} closed subintervals; interval arithmetic '
               'on each gives an enclosure valid for EVERY delta in that subinterval, so '
               'the union covers [1/4, 49/100] with no gaps and no grid sampling')
    Qi = Iv(worst[1] - abs(worst[1]), worst[1]) if worst[1] < 0 else Iv(worst[1])
    print(f'    worst subinterval upper bound = {float(worst[1]):.9f} on '
          f'[{float(worst[0][0]):.6f}, {float(worst[0][1]):.6f}]   all negative: '
          f'{not sub}', file=sys.stderr)

    print('[5] corrections A and C', file=sys.stderr)
    Keff_proxy = Fp.sq()                   # = K_eff for an exact node; F(14) ~ 0 here
    out['correction_A_taylor'] = {}
    for dnum, dden in ((1, 10), (3, 20), (166, 1000), (167, 1000), (1, 5), (2, 5)):
        d = Q(dnum, dden)
        R = Iv(Q(16, 3)) * Iv(L)**5 * Iv(d)**4 * iv_exp(Iv(2*L*d))
        ratio = R / (Iv(2) * Keff_proxy * Iv(d*d))
        out['correction_A_taylor'][f'{dnum}/{dden}'] = dict(
            ratio_R_over_2Keff_d2=_d(ratio), holds=bool(ratio.hi <= 1))
    out['correction_C_highfreq'] = {}
    for g in (14, 70, 140, 280):
        pe = PEr_norm2(g)
        out['correction_C_highfreq'][str(g)] = _d(pe)

    with open('cold_certificate.json', 'w') as fh:
        json.dump(out, fh, indent=1)

    # ---- hard gates: the script must fail rather than relax a bound
    gates = []
    gates.append(('Q_zeta enclosure is strictly positive', Qz.lo > 0))
    gates.append(('point delta=2/5 strictly negative', Qa.hi < 0))
    gates.append(('interval [1/4,49/100] strictly negative everywhere',
                  bool(worst[1] < 0)))
    gates.append(('arch tail allowance below the point margin',
                  r['arch_tail'].hi < -Qa.hi))
    print()
    for name, ok in gates:
        print(f'  GATE  {name:52s} {"PASS" if ok else "FAIL"}')
    if not all(ok for _, ok in gates):
        raise CertificationFailure('one or more gates failed; no certificate is claimed')
    print('\nwrote cold_certificate.json, interval_cells.json, regenerated_qzeta_cache.json')
    return out


if __name__ == '__main__':
    main()
