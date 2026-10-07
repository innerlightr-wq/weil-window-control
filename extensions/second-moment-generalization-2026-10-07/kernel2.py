"""Core library for the reciprocal-pair / second-moment exploration.

Two INDEPENDENT routes to lambda_min, both in Decimal (no mpmath on this box):

  route A  closed form  lambda_pm = r +/- sqrt(p q),  p=<u,u>, q=<w,w>, r=<u,w>
  route B  a general symmetric Jacobi eigensolver that knows nothing about rank 2

Conventions. mu is a finite positive measure given by nodes x_i and weights w_i > 0
(counting measure = all w_i = 1; a quadrature rule approximates Lebesgue). L2(mu) is
represented in the Euclidean orthonormal frame e_i = indicator_i / sqrt(w_i), in which

    (T_a)_ij = sqrt(w_i w_j) * 2 cosh(a (x_i - x_j)).

That frame is used throughout; it is NOT Toeplitz unless nodes are equispaced and
weights constant, and we never call it Toeplitz otherwise.
"""
from decimal import Decimal as D, getcontext

getcontext().prec = 60
TWO, ONE, ZERO = D(2), D(1), D(0)

def dcosh(z):  return (z.exp() + (-z).exp()) / TWO
def dsinh(z):  return (z.exp() - (-z).exp()) / TWO

# ---------------------------------------------------------------- measure helpers
def mass(xs, ws):            return sum(ws, ZERO)
def mean(xs, ws):            return sum((w*x for x, w in zip(xs, ws)), ZERO) / mass(xs, ws)
def centred_moment(xs, ws):
    m = mean(xs, ws)
    return sum((w*(x-m)**2 for x, w in zip(xs, ws)), ZERO)       # V_mu (UNnormalised)

# ---------------------------------------------------------------- route A: closed form
def lambda_pm_closed(xs, ws, a):
    """lambda_pm = r +/- sqrt(p q) with u=e^{ax}, w=e^{-ax} in L2(mu)."""
    p = sum((wt*(TWO*a*x).exp() for x, wt in zip(xs, ws)), ZERO)   # int e^{2ax} dmu
    q = sum((wt*(-TWO*a*x).exp() for x, wt in zip(xs, ws)), ZERO)  # int e^{-2ax} dmu
    r = mass(xs, ws)                                               # <u,w> = int dmu = M
    root = (p*q).sqrt()
    return r + root, r - root

# ---------------------------------------------------------------- route B: Jacobi
def jacobi_eigenvalues(Ain, sweeps=300):
    """Cyclic Jacobi for a real symmetric matrix of Decimals. Returns sorted eigenvalues.

    Independent of any rank assumption. The 2x2 block (p,q) is updated from SAVED
    originals; updating it in place from already-rotated columns corrupts the block and
    silently breaks trace preservation (that bug was caught by the PSD control below).
    """
    n = len(Ain)
    A = [row[:] for row in Ain]
    for _ in range(sweeps):
        off = sum(A[i][j]**2 for i in range(n) for j in range(n) if i != j)
        if off <= D(10) ** (-(getcontext().prec - 8)):
            break
        for p in range(n-1):
            for q in range(p+1, n):
                if A[p][q] == 0:
                    continue
                theta = (A[q][q] - A[p][p]) / (TWO * A[p][q])
                sgn = ONE if theta >= 0 else -ONE
                t = sgn / (abs(theta) + (theta*theta + ONE).sqrt())
                c = ONE / (t*t + ONE).sqrt()
                s = t * c
                app, aqq, apq = A[p][p], A[q][q], A[p][q]          # saved originals
                col_p = [A[k][p] for k in range(n)]
                col_q = [A[k][q] for k in range(n)]
                for k in range(n):
                    if k == p or k == q:
                        continue
                    A[k][p] = A[p][k] = c*col_p[k] - s*col_q[k]
                    A[k][q] = A[q][k] = s*col_p[k] + c*col_q[k]
                A[p][p] = c*c*app - TWO*s*c*apq + s*s*aqq
                A[q][q] = s*s*app + TWO*s*c*apq + c*c*aqq
                A[p][q] = A[q][p] = ZERO
    return sorted(A[i][i] for i in range(n))


def T_matrix(xs, ws, a):
    """(T_a)_ij = sqrt(w_i w_j) 2 cosh(a(x_i - x_j)) in the Euclidean orthonormal frame."""
    n = len(xs)
    return [[(ws[i]*ws[j]).sqrt() * TWO * dcosh(a*(xs[i]-xs[j])) for j in range(n)]
            for i in range(n)]

def lambda_min_jacobi(xs, ws, a, G=None):
    M = T_matrix(xs, ws, a)
    if G is not None:
        M = [[M[i][j] + G[i][j] for j in range(len(M))] for i in range(len(M))]
    return jacobi_eigenvalues(M)[0]
