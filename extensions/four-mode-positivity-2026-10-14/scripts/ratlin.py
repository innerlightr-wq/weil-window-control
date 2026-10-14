"""Exact rational linear algebra. No floats, no formed inverses."""
from fractions import Fraction as Q


def solve(Amat, bvec):
    """Exact Gaussian elimination with partial pivoting. Returns x with A x = b,
    or None if A is singular."""
    n = len(bvec)
    Am = [[Q(Amat[i][j]) for j in range(n)] + [Q(bvec[i])] for i in range(n)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(Am[r][c]))
        if Am[p][c] == 0:
            return None
        Am[c], Am[p] = Am[p], Am[c]
        pv = Am[c][c]
        for j in range(c, n + 1):
            Am[c][j] /= pv
        for r in range(n):
            if r != c and Am[r][c] != 0:
                f = Am[r][c]
                for j in range(c, n + 1):
                    Am[r][j] -= f * Am[c][j]
    return [Am[i][n] for i in range(n)]


def quadform(Amat, x):
    n = len(x)
    return sum(Q(Amat[i][j]) * Q(x[i]) * Q(x[j]) for i in range(n) for j in range(n))


def dot(u, v):
    return sum(Q(a) * Q(b) for a, b in zip(u, v))


def matvec(Amat, x):
    n = len(x)
    return [sum(Q(Amat[i][j]) * Q(x[j]) for j in range(n)) for i in range(n)]


def charpoly(Amat):
    """Coefficients c[0..n] of det(t I - A) = sum c[i] t^i, exact, by Faddeev-LeVerrier."""
    n = len(Amat)
    A = [[Q(Amat[i][j]) for j in range(n)] for i in range(n)]
    Mk = [[Q(1) if i == j else Q(0) for j in range(n)] for i in range(n)]
    cs = [Q(1)]
    for k in range(1, n + 1):
        AM = [[sum(A[i][t] * Mk[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
        ck = -sum(AM[i][i] for i in range(n)) / k
        cs.append(ck)
        Mk = [[AM[i][j] + (ck if i == j else Q(0)) for j in range(n)] for i in range(n)]
    return list(reversed(cs))        # c[0] + c[1] t + ... + c[n] t^n


def polyval(cs, t):
    t = Q(t); r = Q(0)
    for c in reversed(cs):
        r = r * t + Q(c)
    return r


def sturm_count_neg(cs, lo, hi, depth=80):
    """#sign changes of p on [lo,hi] by bisection on sign changes (p squarefree
    assumed only for counting distinct roots; used only as a diagnostic)."""
    return None


def min_root_bracket(cs, lo, hi, iters=200):
    """If p(lo) and p(hi) have opposite signs, bisect to an exact rational
    bracket of a root. Returns (a,b) with p(a)p(b) < 0, b-a <= (hi-lo)/2^iters."""
    lo, hi = Q(lo), Q(hi)
    fl, fh = polyval(cs, lo), polyval(cs, hi)
    if fl == 0: return (lo, lo)
    if fh == 0: return (hi, hi)
    if (fl > 0) == (fh > 0):
        return None
    for _ in range(iters):
        mid = (lo + hi) / 2
        fm = polyval(cs, mid)
        if fm == 0:
            return (mid, mid)
        if (fm > 0) == (fl > 0):
            lo, fl = mid, fm
        else:
            hi, fh = mid, fm
    return (lo, hi)


def cholesky_pd(Amat):
    """Exact rational LDL^T. Returns the list of pivots d_i if it completes;
    A is positive definite iff every pivot is > 0. None if a zero pivot occurs."""
    n = len(Amat)
    A = [[Q(Amat[i][j]) for j in range(n)] for i in range(n)]
    d = []
    Lm = [[Q(0)]*n for _ in range(n)]
    for j in range(n):
        s = A[j][j] - sum(Lm[j][k]**2 * d[k] for k in range(j))
        if s == 0:
            return None
        d.append(s)
        Lm[j][j] = Q(1)
        for i in range(j+1, n):
            Lm[i][j] = (A[i][j] - sum(Lm[i][k]*Lm[j][k]*d[k] for k in range(j))) / s
    return d


def inertia(Amat):
    """Exact inertia (n_neg, n_zero, n_pos) of a symmetric matrix by LDL^T with
    symmetric (diagonal) pivoting; Sylvester's law of inertia. Returns None only
    if a 2x2 indefinite block is needed and none of the remaining diagonals is
    nonzero (does not occur for the matrices used here)."""
    n = len(Amat)
    A = [[Q(Amat[i][j]) for j in range(n)] for i in range(n)]
    idx = list(range(n))
    neg = zer = pos = 0
    while idx:
        # pick a nonzero diagonal pivot
        pi = None
        for t in idx:
            if A[t][t] != 0:
                pi = t; break
        if pi is None:
            # all remaining diagonals zero; if the whole remaining block is zero
            # the eigenvalues are zero, otherwise a hyperbolic pair contributes
            # one positive and one negative
            nz = any(A[i][j] != 0 for i in idx for j in idx)
            if not nz:
                zer += len(idx); idx = []
                break
            pos += 1; neg += 1
            rest = [t for t in idx]
            if len(rest) <= 2:
                zer += len(rest) - 2; idx = []
                break
            return None
        d = A[pi][pi]
        if d > 0: pos += 1
        else: neg += 1
        idx = [t for t in idx if t != pi]
        for i in idx:
            f = A[i][pi] / d
            if f != 0:
                for j in idx:
                    A[i][j] -= f * A[pi][j]
    return (neg, zer, pos)


def min_eig_bracket(Amat, lo=None, hi=None, iters=120):
    """Bracket the smallest eigenvalue of a symmetric matrix using the exact
    inertia of A - t I (Sylvester). Returns (a, b) with a <= lam_min <= b."""
    n = len(Amat)
    if lo is None:
        lo = min(Q(Amat[i][i]) - sum(abs(Q(Amat[i][j])) for j in range(n) if j != i)
                 for i in range(n)) - 1
    if hi is None:
        hi = max(Q(Amat[i][i]) + sum(abs(Q(Amat[i][j])) for j in range(n) if j != i)
                 for i in range(n)) + 1
    lo, hi = Q(lo), Q(hi)
    def nneg(t):
        B = [[Q(Amat[i][j]) - (t if i == j else 0) for j in range(n)] for i in range(n)]
        r = inertia(B)
        return None if r is None else r[0]
    for _ in range(iters):
        mid = (lo + hi) / 2
        k = nneg(mid)
        if k is None:
            return (lo, hi)
        if k >= 1:
            hi = mid          # at least one eigenvalue < mid
        else:
            lo = mid
    return (lo, hi)


def eig_brackets(Amat, iters=140):
    """All eigenvalues of a symmetric matrix, each as a rational bracket, by
    bisection on the exact inertia:  #{eigenvalues < t} = n_neg(A - t I)
    (Sylvester's law of inertia).  Returns a list of (lo, hi) pairs, sorted
    ascending.  No deflation, no floating point."""
    n = len(Amat)
    A = [[Q(Amat[i][j]) for j in range(n)] for i in range(n)]
    def count_below(t):
        B = [[A[i][j] - (t if i == j else 0) for j in range(n)] for i in range(n)]
        ne, ze, po = inertia(B)
        return ne, ze
    R0 = max(sum(abs(A[i][j]) for j in range(n)) for i in range(n)) + 1   # Gershgorin
    out = []
    for k in range(1, n + 1):          # k-th smallest
        lo, hi = -R0, R0
        for _ in range(iters):
            t = (lo + hi) / 2
            ne, ze = count_below(t)
            if ne >= k: hi = t
            else: lo = t
        out.append((lo, hi))
    return out
