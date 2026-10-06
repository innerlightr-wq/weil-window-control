r"""Evaluate the zeta geometric-side Weil form at ONE trial vector, without ever
building the N x N matrix.

Why this exists: certifying indefiniteness needs the EXACT on-line contribution.  A
truncated sum over the first K zeros omits positive mass far larger than the tiny negative
Rayleigh quotients near the detection onset, so it cannot certify anything there.  The
geometric side has no truncation error -- but building the matrix is O(N^2) Lerch/digamma
evaluations, which is hopeless at the N needed to resolve high planted heights.

Two reductions make a single-vector evaluation cheap:

 (1) ARCHIMEDEAN AS ONE INTEGRAL.  Summing the m-series against the kernel gives
        Arch(C) = C(0) psi(1/4) + 2 \int_0^inf [C(0) - Ctilde(u)] e^{-u/2}/(1-e^{-2u}) du
     (Ctilde = C on [0,2L], 0 beyond).  The integrand is O(u) at 0 since the kernel is
     ~1/(2u) and C(0)-C(u) = O(u^2).  Verified against the digamma/Lerch closed form to
     1e-30 in both sectors.

 (2) O(N) CORRELATION.  Every phase in C_jk is an integer multiple of pi: for EVEN_D,
     (w_j - w_k)L = (j-k) pi and (w_j + w_k)L = (j+k+1) pi (similarly in the other
     sectors).  Since sin(x + m pi) = (-1)^m sin(x), the whole O(N^2) term list collapses
     to  C_vv(u) = sum_k [ A_k sin(w_k u) + B_k (2L-u) cos(w_k u) ],  only O(N) terms,
     with the A_k, B_k assembled once in O(N^2).
"""
import mpmath as mp
from zeta_window import (Csym_terms, eval_terms, freq, norm2, F_at_half, pole_sign,
                         von_mangoldt_terms, EVEN, EVEN_D, ODD)


def collapse(terms, L, tol=None):
    """Collapse [(c, alpha, beta)] with beta in pi Z into
    (dict w -> A_w  for A sin(wu),  dict w -> B_w  for B (2L-u)cos(wu))."""
    tol = tol or mp.mpf(10) ** (-(mp.mp.dps // 2))
    A, B = {}, {}
    for c, alpha, beta in terms:
        if alpha == 'lin':
            key = mp.nstr(beta, 25)
            B[key] = B.get(key, mp.mpf(0)) + c
            continue
        m = beta / mp.pi
        mr = mp.nint(m)
        assert abs(m - mr) < tol, f"phase {beta} is not an integer multiple of pi"
        s = c * (-1) ** int(mr)
        w = alpha
        if w < 0:                      # sin(-wu) = -sin(wu)
            w, s = -w, -s
        key = mp.nstr(w, 25)
        A[key] = A.get(key, mp.mpf(0)) + s
    return A, B


def Cvv_builder(a, L, N, sector):
    """a = coefficients against the RAW basis f_k (not normalised). Returns C_vv(u)."""
    A, B = {}, {}
    for j in range(N):
        if a[j] == 0:
            continue
        for k in range(N):
            if a[k] == 0:
                continue
            w = a[j] * a[k] if j == k else a[j] * a[k]
            Aj, Bj = collapse(Csym_terms(j, k, L, sector), L)
            for key, val in Aj.items():
                A[key] = A.get(key, mp.mpf(0)) + w * val
            for key, val in Bj.items():
                B[key] = B.get(key, mp.mpf(0)) + w * val
    Aw = [(mp.mpf(k), v) for k, v in A.items()]
    Bw = [(mp.mpf(k), v) for k, v in B.items()]

    def C(u):
        u = abs(mp.mpf(u))
        if u >= 2 * L:
            return mp.mpf(0)
        tot = mp.mpf(0)
        for w, c in Aw:
            tot += c * mp.sin(w * u)
        for w, c in Bw:
            tot += c * (2 * L - u) * mp.cos(w * u)
        return tot
    return C


def arch_integral(C, C0, L):
    ker = lambda u: mp.e ** (-u / 2) / (1 - mp.e ** (-2 * u))
    I1 = mp.quad(lambda u: (C0 - C(u)) * ker(u), [0, L / 4, L, 2 * L])
    I2 = mp.quad(lambda u: C0 * ker(u), [2 * L, 4 * L, 20 * L, mp.inf])
    return C0 * mp.digamma(mp.mpf(1) / 4) + 2 * (I1 + I2)


def Q_geom_at(v, L, dps, sector=EVEN_D):
    """Weil geometric side at f = sum_k v_k * (f_k / ||f_k||).  Returns Q (NOT divided
    by ||v||^2)."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        L = mp.mpf(L)
        N = len(v)
        a = [v[k] / mp.sqrt(norm2(k, L, sector)) for k in range(N)]
        C = Cvv_builder(a, L, N, sector)
        C0 = C(0)
        pole = 2 * pole_sign(sector) * (sum(a[k] * F_at_half(k, L, sector)
                                            for k in range(N))) ** 2
        arch = arch_integral(C, C0, L)
        prime = sum(2 * lam / mp.sqrt(mp.mpf(n)) * C(logn)
                    for n, lam, logn in von_mangoldt_terms(2 * L))
        return pole + arch - mp.log(mp.pi) * C0 - prime, C0
    finally:
        mp.mp.dps = old
