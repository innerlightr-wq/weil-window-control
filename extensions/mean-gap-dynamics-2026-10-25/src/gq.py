"""Exact Gaussian rationals Q(i), just enough for the toy model."""
from fractions import Fraction as Fr

class GQ:
    __slots__ = ("re", "im")
    def __init__(self, re=0, im=0):
        self.re = Fr(re); self.im = Fr(im)
    def __add__(s, o): o = _c(o); return GQ(s.re+o.re, s.im+o.im)
    def __radd__(s, o): return s + o
    def __neg__(s): return GQ(-s.re, -s.im)
    def __sub__(s, o): return s + (-_c(o))
    def __rsub__(s, o): return _c(o) + (-s)
    def __mul__(s, o):
        o = _c(o); return GQ(s.re*o.re - s.im*o.im, s.re*o.im + s.im*o.re)
    def __rmul__(s, o): return s * o
    def inv(s):
        n = s.re*s.re + s.im*s.im
        if n == 0: raise ZeroDivisionError
        return GQ(s.re/n, -s.im/n)
    def __truediv__(s, o): return s * _c(o).inv()
    def __rtruediv__(s, o): return _c(o) * s.inv()
    def conj(s): return GQ(s.re, -s.im)
    def abs2(s): return s.re*s.re + s.im*s.im
    def __eq__(s, o): o = _c(o); return s.re == o.re and s.im == o.im
    def __repr__(s):
        return f"({s.re}{'+' if s.im>=0 else '-'}{abs(s.im)}i)"

def _c(o):
    return o if isinstance(o, GQ) else GQ(o, 0)

I = GQ(0, 1)
