"""Exact inertia of a rational symmetric Toeplitz matrix, over Q. (T1)

Sylvester / Jacobi rule: if all leading principal minors D_1..D_n are nonzero, then
  #negative eigenvalues = #sign changes in the sequence D_0=1, D_1, ..., D_n.
We compute the D_k exactly with fraction-free Gaussian elimination (Bareiss-style via
Fraction), so the inertia -- in particular the SIGN of lambda_min -- is certified
without any floating point.

Rational planted data: t(n) = 2 cos(n phi) (rho^n + rho^-n) is rational whenever
rho and c = cos(phi) are rational, because cos(n phi) = T_n(c) (Chebyshev, integer
coefficients in c).
"""
from fractions import Fraction


def cheb_T(n, c):
    """T_n(c) = cos(n*arccos c), exact for rational c."""
    a, b = Fraction(1), Fraction(c)
    if n == 0:
        return a
    for _ in range(n - 1):
        a, b = b, 2 * Fraction(c) * b - a
    return b


def rational_t(blocks, nmax):
    """blocks: ('online', c) / ('quartet', rho, c) / ('realpair', rho) with rho, c
    rational (c = cos phi). Returns exact Fraction t[0..nmax] and g_eff."""
    t = [Fraction(0)] * (nmax + 1)
    cnt = 0
    for b in blocks:
        if b[0] == 'online':
            c = Fraction(b[1]); cnt += 2
            for n in range(nmax + 1):
                t[n] += 2 * cheb_T(n, c)
        elif b[0] == 'quartet':
            rho, c = Fraction(b[1]), Fraction(b[2]); cnt += 4
            for n in range(nmax + 1):
                t[n] += 2 * cheb_T(n, c) * (rho ** n + rho ** (-n))
        elif b[0] == 'realpair':
            rho = Fraction(b[1]); cnt += 2
            for n in range(nmax + 1):
                t[n] += rho ** n + rho ** (-n)
        else:
            raise ValueError(b)
    return t, cnt // 2


def leading_minors(t, R):
    """Exact leading principal minors D_1..D_{R+1} of Toeplitz(t) via Fraction LU."""
    n = R + 1
    A = [[Fraction(t[abs(i - j)]) for j in range(n)] for i in range(n)]
    D, piv = [], []
    for k in range(n):
        # eliminate using rows 0..k-1 already reduced in place (no row swaps:
        # we need LEADING minors, so pivoting is not allowed)
        for i in range(k):
            if piv[i] == 0:
                return None            # a leading minor vanished: rule inapplicable
            f = A[k][i] / piv[i]
            if f:
                for j in range(i, n):
                    A[k][j] -= f * A[i][j]
        piv.append(A[k][k])
        D.append((D[-1] * A[k][k]) if D else A[k][k])
        if A[k][k] == 0:
            return None
    return D


def exact_inertia(t, R):
    """(n_neg, n_pos) exactly, or None if a leading minor vanishes."""
    D = leading_minors(t, R)
    if D is None:
        return None
    seq = [Fraction(1)] + D
    neg = sum(1 for i in range(len(seq) - 1) if seq[i] * seq[i + 1] < 0)
    return (neg, (R + 1) - neg)


def exact_lambda_min_is_negative(t, R):
    """True / False / None(undecided). lam_min < 0  <=>  n_neg >= 1."""
    r = exact_inertia(t, R)
    return None if r is None else (r[0] >= 1)
