"""Rigorous rational interval arithmetic, Python stdlib only.

Every function returns an Interval [lo,hi] of Fractions that provably contains
the true value.  No floating point is used anywhere in this module.
"""
from fractions import Fraction as F

class IV:
    __slots__ = ("lo", "hi")
    def __init__(self, lo, hi=None):
        if hi is None:
            hi = lo
        self.lo = F(lo); self.hi = F(hi)
        if self.lo > self.hi:
            raise ValueError("empty interval")
    def __repr__(self):
        return "[%s, %s]" % (float(self.lo), float(self.hi))
    def __add__(s, o):
        o = iv(o); return IV(s.lo + o.lo, s.hi + o.hi)
    __radd__ = __add__
    def __neg__(s): return IV(-s.hi, -s.lo)
    def __sub__(s, o): return s + (-iv(o))
    def __rsub__(s, o): return iv(o) + (-s)
    def __mul__(s, o):
        o = iv(o)
        c = [s.lo*o.lo, s.lo*o.hi, s.hi*o.lo, s.hi*o.hi]
        return IV(min(c), max(c))
    __rmul__ = __mul__
    def inv(s):
        if s.lo <= 0 <= s.hi: raise ZeroDivisionError("interval straddles 0")
        return IV(1/s.hi, 1/s.lo)
    def __truediv__(s, o): return s * iv(o).inv()
    def __rtruediv__(s, o): return iv(o) * s.inv()
    def __pow__(s, n):
        assert isinstance(n, int) and n >= 0
        r = IV(1)
        for _ in range(n): r = r * s
        return r
    def contains(s, x): return s.lo <= F(x) <= s.hi
    def pos(s): return s.lo > 0
    def neg(s): return s.hi < 0
    def width(s): return s.hi - s.lo
    def mid(s): return (s.lo + s.hi) / 2

def iv(x):
    return x if isinstance(x, IV) else IV(x)

# ---- elementary functions with rigorous remainders -------------------------
# All series below are evaluated at |x| <= 4; terms shrink monotonically once
# n > |x|, and we bound the tail by a geometric series with ratio 1/2.

def _series(x, sign_alt, odd, N=80):
    """sum_{k} s^k x^{m}/m!  with m = 2k (+1 if odd); s = -1 if sign_alt."""
    x = iv(x)
    B = max(abs(x.lo), abs(x.hi))
    tot = IV(0)
    term_num = 1          # exact rational factorial bookkeeping
    acc = IV(0)
    for k in range(N):
        m = 2*k + (1 if odd else 0)
        fact = 1
        for j in range(2, m+1): fact *= j
        t = (x ** m) * F(1, fact)
        if sign_alt and (k % 2 == 1): t = -t
        acc = acc + t
    # tail bound: next term magnitude * 2  (ratio < 1/2 since (2N)! grows fast)
    m = 2*N + (1 if odd else 0)
    fact = 1
    for j in range(2, m+1): fact *= j
    tail = F(B) ** m / fact * 2
    return IV(acc.lo - tail, acc.hi + tail)

def sin_iv(x): return _series(x, True, True)
def cos_iv(x): return _series(x, True, False)
def sinh_iv(x): return _series(x, False, True)
def cosh_iv(x): return _series(x, False, False)
def tan_iv(x): return sin_iv(x) / cos_iv(x)
def cot_iv(x): return cos_iv(x) / sin_iv(x)

# sqrt by interval Newton / bisection on an exact rational bracket
def sqrt_iv(x, eps=F(1, 10**40)):
    x = iv(x)
    assert x.lo >= 0
    def rt(v):
        if v == 0: return F(0), F(0)
        lo, hi = F(0), F(max(1, v)) + 1
        while hi - lo > eps:
            m = (lo + hi) / 2
            if m*m <= v: lo = m
            else: hi = m
        return lo, hi
    a, _ = rt(x.lo); _, b = rt(x.hi)
    return IV(a, b)

PI = IV(F(31415926535897932384626433832795028841971,  10**40),
        F(31415926535897932384626433832795028841972,  10**40))

def coarsen(x, d=10**25):
    """Widen to an interval with small denominators (keeps rigour, kills blow-up)."""
    x = iv(x)
    lo = F(int(x.lo * d) - 1, d)
    hi = F(int(x.hi * d) + 1, d)
    return IV(lo, hi)
