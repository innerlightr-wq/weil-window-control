"""Brute-force point counts, exact L-polynomial recovery, exact power sums. (T1)"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction
from ffield import GF


def squarefree_over_Fq(f, q):
    """f little-endian integer coeffs. gcd(f, f') == 1 over F_q ?"""
    def norm(a):
        a = [c % q for c in a]
        while len(a) > 1 and a[-1] == 0:
            a.pop()
        return a
    def gcd(a, b):
        a, b = norm(a), norm(b)
        while not (len(b) == 1 and b[0] == 0):
            # a mod b
            inv = pow(b[-1], -1, q)
            A, dB = list(a), len(b) - 1
            for i in range(len(A) - 1, dB - 1, -1):
                c = (A[i] * inv) % q
                if c:
                    for k in range(dB + 1):
                        A[i - dB + k] = (A[i - dB + k] - c * b[k]) % q
            a, b = b, norm(A)
        return norm(a)
    fp = [(i * c) % q for i, c in enumerate(f)][1:]
    g = gcd(f, fp if fp else [0])
    return len(g) == 1 and g[0] != 0


def count_points_hyperelliptic(f, q, n):
    """#C(F_{q^n}) for the smooth model of y^2=f(x), deg f = 2g+1 odd, f squarefree.
    One point at infinity.  Affine: sum_x (1 + chi(f(x)))."""
    assert (len(f) - 1) % 2 == 1, "need odd degree 2g+1"
    F = GF(q, n)
    total = 0
    for x in F.elements():
        total += 1 + F.chi(F.poly_eval(f, x))
    return total + 1


def newton_p_to_A(p, q, g):
    """Newton's identities: power sums p[1..g] of the 2g roots a_j  ->  coefficients
    A[1..g] of P(T)=prod(1-a_j T)=sum A_i T^i.  Then functional equation fills A[g+1..2g].
    Exact: Fractions in, integers out (checked)."""
    # e_k elementary symmetric; A_i = (-1)^i e_i.  Newton: k e_k = sum_{i=1..k} (-1)^{i-1} e_{k-i} p_i
    e = [Fraction(1)] + [Fraction(0)] * g
    for k in range(1, g + 1):
        s = Fraction(0)
        for i in range(1, k + 1):
            s += Fraction(-1) ** (i - 1) * e[k - i] * p[i]
        e[k] = s / k
    A = [Fraction(0)] * (2 * g + 1)
    for i in range(g + 1):
        A[i] = Fraction(-1) ** i * e[i]
    for i in range(g):                       # A_{2g-i} = q^{g-i} A_i
        A[2 * g - i] = q ** (g - i) * A[i]
    out = []
    for a in A:
        assert a.denominator == 1, f"non-integral L-poly coefficient {a}"
        out.append(int(a))
    return out


def A_to_power_sums(A, nmax):
    """Exact integer power sums p_n = sum_j a_j^n from P(T)=sum A_i T^i via
    Newton's identities in the form  p_n = -n*A_n - sum_{i=1}^{n-1} A_i p_{n-i}
    (valid for the roots 1/a_j reciprocal convention; derived from
     -log P(T) = sum_n p_n T^n / n  =>  P'(T)/P(T) = -sum_n p_n T^{n-1}).
    """
    d = len(A) - 1
    p = [0] * (nmax + 1)
    for n in range(1, nmax + 1):
        # coefficient identity from P'(T) = -P(T) * sum_n p_n T^{n-1}
        An = A[n] if n <= d else 0
        s = -n * An
        for i in range(1, n):
            Ai = A[i] if i <= d else 0
            s -= Ai * p[n - i]
        p[n] = s
    return p
