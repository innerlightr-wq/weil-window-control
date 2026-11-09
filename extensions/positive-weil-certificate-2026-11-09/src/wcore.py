"""Shared primitives: exact rational/integer algebra and Arb ball arithmetic.

Sign-critical paths use python-flint (Arb) balls.  Nothing in this module uses
floating point for a decision; floats appear only in __repr__-style output.
"""
from fractions import Fraction
from flint import arb, ctx

def set_prec(bits):
    ctx.prec = bits
    return bits

# ---------------------------------------------------------------- exact Chebyshev
def cheb_to_poly(coeffs):
    """Expand sum_k coeffs[k] * cos(k*theta) into a polynomial in c = cos(theta).

    Exact: Chebyshev T_k has integer coefficients, computed by the recursion
    T_0 = 1, T_1 = c, T_{k+1} = 2 c T_k - T_{k-1}.  Input coefficients are
    Fractions/ints; output is a list of Fractions, little-endian in c.
    """
    n = len(coeffs)
    T = [[Fraction(1)], [Fraction(0), Fraction(1)]]
    while len(T) < max(n, 2):
        a, b = T[-1], T[-2]
        nxt = [Fraction(0)] * (len(a) + 1)
        for i, v in enumerate(a):
            nxt[i + 1] += 2 * v
        for i, v in enumerate(b):
            nxt[i] -= v
        T.append(nxt)
    out = [Fraction(0)] * n
    for k, ck in enumerate(coeffs):
        for i, v in enumerate(T[k]):
            out[i] += Fraction(ck) * v
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out

def poly_mul(p, q):
    out = [Fraction(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out

def poly_eq(p, q):
    p = list(p) + [Fraction(0)] * max(0, len(q) - len(p))
    q = list(q) + [Fraction(0)] * max(0, len(p) - len(q))
    return all(a == b for a, b in zip(p, q))

# ---------------------------------------------------------------- arithmetic data
def von_mangoldt(n):
    """Lambda(n) as (p, k) with n = p^k, or None if n is not a prime power."""
    m, p = n, 2
    while p * p <= m:
        if m % p == 0:
            q, k = m, 0
            while q % p == 0:
                q //= p
                k += 1
            return (p, k) if q == 1 else None
        p += 1
    return (m, 1) if m >= 2 else None

def prime_powers_below(two_L_num, two_L_den):
    """All n >= 2 with log n < 2L, decided EXACTLY by comparing n^den < e^{num}
    through an Arb ball on e^{num/den}; returns [(n, p, k)] for prime powers."""
    out = []
    n = 2
    while True:
        lg = arb(n).log()
        twoL = arb(two_L_num) / arb(two_L_den)
        if (lg < twoL) is False and (lg > twoL) is False:
            raise RuntimeError("log %d vs 2L undecided at prec %d" % (n, ctx.prec))
        if lg > twoL:
            break
        pk = von_mangoldt(n)
        if pk:
            out.append((n, pk[0], pk[1]))
        n += 1
        if n > 10 ** 6:
            raise RuntimeError("runaway")
    return out

def comb_mass(two_L_num, two_L_den):
    """A_L = sum_{log n < 2L} 2 Lambda(n)/sqrt(n), as an Arb ball, with its term list."""
    terms = prime_powers_below(two_L_num, two_L_den)
    A = arb(0)
    rows = []
    for n, p, k in terms:
        t = arb(2) * arb(p).log() / arb(n).sqrt()
        rows.append((n, p, k, t))
        A += t
    return A, rows

def ball_str(x, d=12):
    return x.str(d, radius=True)
