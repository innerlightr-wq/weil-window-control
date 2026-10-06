"""Minimal exact finite field F_{q^n} for odd prime q, by polynomial arithmetic.

Elements are tuples of length n over F_q (little-endian coefficients).
Only what Stage 1 needs: multiplication, powering, iteration, quadratic character.
Exact integer arithmetic throughout -- no floats. (T1)
"""
from itertools import product


def _polymulmod(a, b, modpoly, q, n):
    res = [0] * (2 * n - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    res[i + j] = (res[i + j] + ai * bj) % q
    # reduce by modpoly (monic, degree n, little-endian length n+1)
    for d in range(2 * n - 2, n - 1, -1):
        c = res[d]
        if c:
            res[d] = 0
            for k in range(n):
                res[d - n + k] = (res[d - n + k] - c * modpoly[k]) % q
    return tuple(res[:n])


def _is_irreducible(poly, q, n):
    """poly monic little-endian len n+1. Brute force: no roots in any proper subfield
    is not sufficient in general, so we test by trial division over F_q[x] for all
    monic polys of degree <= n//2. n is tiny here."""
    def polydivmod(A, B):
        A = list(A)
        dB = len(B) - 1
        inv = pow(B[-1], -1, q)
        Q = [0] * max(0, len(A) - dB)
        for i in range(len(A) - 1, dB - 1, -1):
            c = (A[i] * inv) % q
            Q[i - dB] = c
            if c:
                for k in range(dB + 1):
                    A[i - dB + k] = (A[i - dB + k] - c * B[k]) % q
        while len(A) > 1 and A[-1] == 0:
            A.pop()
        return Q, A
    for d in range(1, n // 2 + 1):
        for tail in product(range(q), repeat=d):
            B = list(tail) + [1]
            _, r = polydivmod(poly, B)
            if len(r) == 1 and r[0] == 0:
                return False
    return True


class GF:
    def __init__(self, q, n):
        self.q, self.n = q, n
        self.Q = q ** n
        if n == 1:
            self.modpoly = (0, 1)
        else:
            self.modpoly = None
            for tail in product(range(q), repeat=n):
                poly = list(tail) + [1]
                if _is_irreducible(poly, q, n):
                    self.modpoly = tuple(poly)
                    break
            assert self.modpoly is not None, (q, n)
        self.zero = tuple([0] * n)
        self.one = tuple([1] + [0] * (n - 1))

    def elements(self):
        for coeffs in product(range(self.q), repeat=self.n):
            yield coeffs

    def mul(self, a, b):
        if self.n == 1:
            return ((a[0] * b[0]) % self.q,)
        return _polymulmod(a, b, self.modpoly, self.q, self.n)

    def add(self, a, b):
        return tuple((x + y) % self.q for x, y in zip(a, b))

    def pow(self, a, e):
        r, base = self.one, a
        while e:
            if e & 1:
                r = self.mul(r, base)
            base = self.mul(base, base)
            e >>= 1
        return r

    def from_int(self, m):
        return tuple([m % self.q] + [0] * (self.n - 1))

    def chi(self, a):
        """Quadratic character: 0 if a==0 else +-1."""
        if a == self.zero:
            return 0
        r = self.pow(a, (self.Q - 1) // 2)
        return 1 if r == self.one else -1

    def poly_eval(self, coeffs_int, x):
        """Horner for a polynomial with integer coefficients (little-endian)."""
        acc = self.zero
        for c in reversed(coeffs_int):
            acc = self.add(self.mul(acc, x), self.from_int(c))
        return acc
