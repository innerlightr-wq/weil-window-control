r"""Numba kernels: digamma trend + the certified cell bound.

scipy cannot be called from numba, so Re psi(1/4 + it/2) is computed here by
recurrence-then-asymptotic: psi(z) = psi(z+m) - sum_{k<m} 1/(z+k), with m chosen so
|z+m| is large enough for the Stirling series to be at full double accuracy.
"""
import numpy as np
from numba import njit, prange

LOG_PI = np.log(np.pi)
_SPLIT = 134217729.0
TWOPI_HI = 6.283185307179586
TWOPI_LO = 2.4492935982947064e-16
SHIFT = 12                     # recurrence shift; |z| >= 12 before the asymptotic

# Stirling: psi(z) ~ log z - 1/(2z) - sum B_2k/(2k z^2k)
_B = np.array([1.0 / 12.0, -1.0 / 120.0, 1.0 / 252.0, -1.0 / 240.0,
               1.0 / 132.0, -691.0 / 32760.0, 1.0 / 12.0])


@njit(cache=True, inline='always')
def _re_psi_quarter(t):
    """Re psi(1/4 + i t/2)."""
    zr = 0.25
    zi = 0.5 * t
    acc = 0.0
    for k in range(SHIFT):
        ar = zr + k
        d = ar * ar + zi * zi
        acc += ar / d                      # Re 1/(z+k)
    zr += SHIFT
    d2 = zr * zr + zi * zi
    out = 0.5 * np.log(d2)                 # Re log z
    out -= 0.5 * zr / d2                   # Re 1/(2z)
    # -(sum) B_2k / (2k z^2k): accumulate powers of 1/z^2 in complex form
    wr = (zr * zr - zi * zi) / (d2 * d2)   # Re 1/z^2
    wi = -(2.0 * zr * zi) / (d2 * d2)      # Im 1/z^2
    pr, pi_ = wr, wi
    for j in range(len(_B)):
        out -= _B[j] * pr
        nr = pr * wr - pi_ * wi
        ni = pr * wi + pi_ * wr
        pr, pi_ = nr, ni
    return out - acc


@njit(cache=True, inline='always')
def _reduced(t, f):
    ph = t * f
    ca = _SPLIT * t; ahi = ca - (ca - t); alo = t - ahi
    cb = _SPLIT * f; bhi = cb - (cb - f); blo = f - bhi
    pe = ((ahi * bhi - ph) + ahi * blo + alo * bhi) + alo * blo
    q = np.rint(ph / TWOPI_HI)
    qh = q * TWOPI_HI
    cq = _SPLIT * q; qhi = cq - (cq - q); qlo = q - qhi
    ct = _SPLIT * TWOPI_HI; thi = ct - (ct - TWOPI_HI); tlo = TWOPI_HI - thi
    qe = ((qhi * thi - qh) + qhi * tlo + qlo * thi) + qlo * tlo
    return (((ph - qh) - qe) - q * TWOPI_LO) + pe


@njit(cache=True, inline='always')
def _P(t, freqs, weights):
    s = 0.0
    for k in range(freqs.shape[0]):
        s += weights[k] * np.cos(_reduced(t, freqs[k]))
    return s


@njit(parallel=True, cache=True)
def psi_many(t, freqs, weights):
    out = np.empty(t.shape[0])
    for i in prange(t.shape[0]):
        out[i] = _re_psi_quarter(t[i]) - LOG_PI - _P(t[i], freqs, weights)
    return out


@njit(parallel=True, cache=True)
def cell_lower_bound(a, h, freqs, weights, A, lipP):
    r"""Rigorous lower bound on Psi_L over each cell [a_i, a_i + h].

    min Psi >= H(a) - min(A, P(mid) + lipP*h/2), using H increasing (notes/monotonicity.md)
    and P <= A unconditionally.
    """
    n = a.shape[0]
    out = np.empty(n)
    half = 0.5 * h
    slack = lipP * half
    for i in prange(n):
        H = _re_psi_quarter(a[i]) - LOG_PI
        pm = _P(a[i] + half, freqs, weights) + slack
        if pm > A:
            pm = A
        out[i] = H - pm
    return out


@njit(parallel=True, cache=True)
def cell_bounds(a, h, freqs, weights, A, lipP):
    r"""Two-sided rigorous bounds on Psi_L over each cell [a_i, a_i+h].

    lower >= H(a) - min(A,  P(mid) + lipP h/2)        (H increasing; P <= A)
    upper <= H(a+h) - max(-A, P(mid) - lipP h/2)

    The upper bound is what makes the scan cheap: a cell certified ENTIRELY BELOW beta
    is accepted whole and never subdivided, so the cost tracks the number of dip EDGES
    rather than the dip measure.
    """
    n = a.shape[0]
    lo = np.empty(n); up = np.empty(n)
    half = 0.5 * h
    slack = lipP * half
    for i in prange(n):
        Ha = _re_psi_quarter(a[i]) - LOG_PI
        Hb = _re_psi_quarter(a[i] + h) - LOG_PI
        pm = _P(a[i] + half, freqs, weights)
        pu = pm + slack
        if pu > A:
            pu = A
        pl = pm - slack
        if pl < -A:
            pl = -A
        lo[i] = Ha - pu
        up[i] = Hb - pl
    return lo, up
