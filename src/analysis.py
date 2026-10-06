"""Shared spectral analysis of a Toeplitz Weil form, with the checks the
project rules demand: precision doubling, float64 contrast, parity, H1/H4 metrics."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
import numpy as np
from weilform import toeplitz, eigsym, parity_of, poly_zeros, inertia

TWO_PI = None


def _angdist(a, b):
    d = abs(a - b) % (2 * mp.pi)
    return min(d, 2 * mp.pi - d)


def lam_min_at(t, R, dps):
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        T = toeplitz([mp.mpf(x) for x in t], R)
        vals, vecs = eigsym(T)
        return vals, vecs
    finally:
        mp.mp.dps = old


def analyse_window(t, R, thetas, dps, dps2=None):
    """One row of diagnostics for the (R+1)x(R+1) Toeplitz form built from t."""
    dps2 = dps2 or 2 * dps
    vals, vecs = lam_min_at(t, R, dps)
    vals2, _ = lam_min_at(t, R, dps2)
    lam, lam2 = vals[0], vals2[0]
    scale = max([abs(v) for v in vals] + [mp.mpf(1)])

    # float64 contrast (the Zhu warning, in this setting)
    Tf = np.array([[float(t[abs(i - j)]) for j in range(R + 1)] for i in range(R + 1)])
    lam_f64 = float(np.linalg.eigvalsh(Tf)[0])

    gap = (vals[1] - vals[0]) if R >= 1 else mp.inf
    # Error-based sign test: estimate the achieved accuracy by precision doubling.
    # The sign of lam_min is only asserted when |lam_min| dominates that estimate.
    err = abs(lam - lam2)
    floor = scale * mp.mpf(10) ** (-(dps2 - 5))      # can't resolve below this anyway
    err = max(err, floor)
    certain = abs(lam2) > 100 * err
    sign = ('neg' if lam2 < 0 else 'pos') if certain else 'zero'
    sign2 = sign
    zero_tol = 100 * err
    simple = bool(R == 0 or gap > zero_tol)

    par, pres = parity_of(vecs[0])
    zs = poly_zeros(vecs[0])
    rad = [abs(z) for z in zs]
    maxoff = max((abs(r - 1) for r in rad), default=mp.mpf(0))
    zang = sorted(mp.arg(z) for z in zs)

    # H4, both directions
    fwd = max((min(_angdist(za, th) for th in thetas) for za in zang), default=None)
    bwd = max((min(_angdist(th, za) for za in zang) for th in thetas), default=None) \
        if zang else None

    (nneg, nzer, npos), _, itol = inertia(toeplitz([mp.mpf(x) for x in t], R))

    return dict(
        R=R, size=R + 1, dps=dps,
        lam_min=mp.nstr(lam, 12), lam_min_sign=sign,
        lam_min_dps2=mp.nstr(lam2, 12), lam_min_err_est=mp.nstr(err, 4),
        sign_certain=bool(certain),
        prec_stable=bool(certain),
        lam_min_float64='%.4e' % lam_f64,
        float64_sign=('neg' if lam_f64 < 0 else 'pos'),
        gap_lam2_lam1=mp.nstr(gap, 8), lam_min_simple=simple,
        parity=par, parity_residual=mp.nstr(pres, 4),
        n_zeros=len(zs),
        H1_max_abs_zero_minus_1=mp.nstr(maxoff, 6),
        H1_applicable=simple,
        H1_holds=(bool(maxoff <= mp.mpf(10) ** (-(dps // 2))) if simple else None),
        H4_zero_to_theta=(mp.nstr(fwd, 6) if fwd is not None else ''),
        H4_theta_to_zero=(mp.nstr(bwd, 6) if bwd is not None else ''),
        inertia=f"({nneg},{nzer},{npos})", inertia_neg=nneg, inertia_zero=nzer,
        inertia_pos=npos, inertia_tol=mp.nstr(itol, 4),
    ), vecs[0], zs
