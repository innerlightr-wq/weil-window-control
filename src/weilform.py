"""Toeplitz Weil form utilities in mpmath (high precision, per project rules)."""
import mpmath as mp


def toeplitz(t, R):
    return mp.matrix(R + 1, R + 1).__class__(
        [[t[abs(i - j)] for j in range(R + 1)] for i in range(R + 1)])


def eigsym(T):
    """Ascending eigenvalues and eigenvectors (columns) of a real symmetric mp.matrix."""
    E, Q = mp.eigsy(T)
    idx = sorted(range(T.rows), key=lambda k: E[k])
    vals = [E[k] for k in idx]
    vecs = [[Q[i, k] for i in range(T.rows)] for k in idx]   # vecs[m] = m-th eigenvector
    return vals, vecs


def parity_of(v, tol=None):
    """Classify v against the flip J: 'even', 'odd', or 'mixed' (+ residual)."""
    n = len(v)
    nrm = mp.sqrt(sum(x ** 2 for x in v))
    w = [x / nrm for x in v]
    rev = list(reversed(w))
    de = mp.sqrt(sum((a - b) ** 2 for a, b in zip(w, rev)))      # ||v - Jv||
    do = mp.sqrt(sum((a + b) ** 2 for a, b in zip(w, rev)))      # ||v + Jv||
    res = min(de, do)
    if tol is None:
        tol = mp.mpf(10) ** (-(mp.mp.dps // 2))
    if de <= tol:
        return 'even', res
    if do <= tol:
        return 'odd', res
    return 'mixed', res


def poly_zeros(c):
    """Zeros of sum_i c_i z^i, as mpmath complex. c little-endian."""
    cc = list(c)
    while len(cc) > 1 and abs(cc[-1]) < mp.mpf(10) ** (-(mp.mp.dps - 5)):
        cc.pop()
    if len(cc) <= 1:
        return []
    coeffs = list(reversed(cc))            # mpmath wants highest-degree first
    return mp.polyroots(coeffs, maxsteps=200, extraprec=200)


def inertia(T, tol=None):
    """(n_neg, n_zero, n_pos) from high-precision eigenvalues, with the scale used."""
    vals, _ = eigsym(T)
    scale = max([abs(v) for v in vals] + [mp.mpf(1)])
    if tol is None:
        tol = scale * mp.mpf(10) ** (-(mp.mp.dps - 10))
    neg = sum(1 for v in vals if v < -tol)
    zer = sum(1 for v in vals if abs(v) <= tol)
    pos = sum(1 for v in vals if v > tol)
    return (neg, zer, pos), vals, tol
