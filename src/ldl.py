r"""LDL^T without pivoting, used as a cheap one-sided indefiniteness certificate.

Why not eigsy:  for a symmetric M, LDL^T without row/column interchanges computes the
LEADING principal minors, so the pivots D_1..D_k are exactly the pivots of the leading
k x k block M_k.  By Sylvester/Jacobi the inertia of M_k is read off the signs of
D_1..D_k.  Hence ONE factorisation of M_{Nmax} gives the inertia of every leading block,
i.e. the whole basis-size sweep -- provided the basis is NESTED, which ours is
(Z_N is the leading N x N block of Z_{Nmax}).

Witness:  if D_j < 0, solve L^T v = e_j.  Then
    v^T M v = (L^T v)^T D (L^T v) = e_j^T D e_j = D_j < 0,
so v is an explicit trial vector with negative Rayleigh quotient, supported on the first
j+1 basis functions.  THAT is what makes a "detected" entry a genuine certificate: a
single trial function with Q(f) < 0 proves the form is indefinite, independently of any
basis-completeness question.
"""
import mpmath as mp


def ldl(M, n=None):
    """Returns (Lmat, D) with M = Lmat * diag(D) * Lmat^T, Lmat unit lower triangular.
    Raises ZeroDivisionError if a pivot vanishes (then the rule does not apply)."""
    n = n or M.rows
    Lm = mp.zeros(n, n)
    D = [mp.mpf(0)] * n
    for j in range(n):
        s = M[j, j]
        for k in range(j):
            s -= Lm[j, k] ** 2 * D[k]
        D[j] = s
        Lm[j, j] = mp.mpf(1)
        if s == 0:
            raise ZeroDivisionError(f"zero pivot at {j}")
        for i in range(j + 1, n):
            t = M[i, j]
            for k in range(j):
                t -= Lm[i, k] * Lm[j, k] * D[k]
            Lm[i, j] = t / s
    return Lm, D


def first_negative_pivot(D, tol=None):
    """Index of the first negative pivot, or None. tol guards against roundoff:
    a pivot is called negative only if it is below -tol * (largest |pivot| so far)."""
    big = mp.mpf(0)
    for j, d in enumerate(D):
        big = max(big, abs(d))
        thr = (tol or mp.mpf(10) ** (-(mp.mp.dps - 8))) * max(big, mp.mpf(1))
        if d < -thr:
            return j
    return None


def witness(Lm, j, n=None):
    """Solve L^T v = e_j by back substitution; v is supported on indices 0..j."""
    n = n or Lm.rows
    v = [mp.mpf(0)] * n
    v[j] = mp.mpf(1)
    for i in range(j - 1, -1, -1):
        s = mp.mpf(0)
        for k in range(i + 1, j + 1):
            s += Lm[k, i] * v[k]
        v[i] = -s
    return v


def rayleigh(M, v, n=None):
    n = n or len(v)
    num = mp.mpf(0)
    for i in range(n):
        if v[i] == 0:
            continue
        for k in range(n):
            if v[k] != 0:
                num += v[i] * M[i, k] * v[k]
    den = sum(x ** 2 for x in v)
    return num / den, num, den
