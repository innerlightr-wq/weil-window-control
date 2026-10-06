r"""Interval-arithmetic certification of a negative Rayleigh quotient (mpmath.iv).

WHAT IS AND IS NOT CERTIFIED.  Three distinct levels, kept separate on purpose:

 (1) CERTIFIED (interval arithmetic): that the assembled K-zero + planted form, evaluated
     at the trial vector v, is negative -- i.e. the ARITHMETIC is not the reason for the
     sign.  Enclosures are taken for every zero ordinate and every basis transform.
 (2) BOUNDED, not certified: the contribution of the zeros above gamma_K, which is
     omitted from the assembly and is >= 0, so it can only push the value UP.  We bound it
     with the Riemann-von Mangoldt density and a safety factor; that bound uses a density
     asymptotic, not a proven inequality, so it is T2.
 (3) NOT CERTIFIED AT ALL: anything about zeta.  The whole Stage 3e construction is built
     from zero data with a hand-planted quartet.  It is a diagnostic.

The zero ordinates come from mpmath.zetazero, enclosed at +-10^{-(dps-5)}; mpmath's
zetazero is not itself a certified routine, so (1) is conditional on those ordinates.
"""
import mpmath as mp


def _ivF(k, L, t, iv):
    r"""F_k(t) for the EVEN_D basis in interval arithmetic.

    mpmath's iv.sincpi is broken (it calls ctx.sinpi, which the interval context does not
    define), so we use the exact closed form instead -- which is tighter anyway:
        cos(w_k L) = 0 and sin(w_k L) = (-1)^k, hence
        F_k(t) = (-1)^k * 2 w_k cos(tL) / (w_k^2 - t^2).
    The only singularity is removable, at t = +- w_k.  The caller must not evaluate there;
    `singular_risk` below flags any such near-coincidence.
    """
    w = (2 * k + 1) * iv.pi / (2 * L)
    num = 2 * w * iv.cos(t * L)
    val = num / (w ** 2 - t ** 2)
    return -val if (k % 2) else val


def singular_risk(L, Nuse, gammas_mid, rel=mp.mpf('1e-6')):
    """Flag basis frequencies that sit too close to a zero ordinate for the closed form
    to be evaluated safely in interval arithmetic."""
    L = mp.mpf(L)
    bad = []
    for k in range(Nuse):
        w = (2 * k + 1) * mp.pi / (2 * L)
        for g in gammas_mid:
            if abs(w - g) < rel * max(w, mp.mpf(1)):
                bad.append((k, mp.nstr(w, 10), mp.nstr(g, 10)))
    return bad


def certify_negative(L, Nuse, v, gammas_mid, planted, dps, zero_halfwidth_exp=None):
    """Interval enclosure of  v^T Z v  (Z = K-zero + planted form, EVEN_D, orthonormal).
    Returns (lo, hi, certified_negative)."""
    iv = mp.iv
    old = iv.dps
    iv.dps = dps
    try:
        Lx = iv.mpf(mp.nstr(L, dps))
        hw = iv.mpf(10) ** (-(zero_halfwidth_exp or (dps - 5)))
        nrm = iv.sqrt(Lx)
        acc = iv.mpf(0)
        # on-line zeros: each term is 2 (sum_k v_k F_k(g)/nrm)^2 >= 0
        for gm in gammas_mid:
            c = iv.mpf(mp.nstr(gm, dps))
            g = iv.mpf([c.a - hw.b, c.b + hw.b])
            s = iv.mpf(0)
            for k in range(Nuse):
                if v[k] == 0:
                    continue
                s += iv.mpf(mp.nstr(v[k], dps)) * _ivF(k, Lx, g, iv) / nrm
            acc += 2 * s ** 2
        # planted quartet: 2 Re( Fw_j Fmw_k + Fw_k Fmw_j ) contracted with v  ->  4 Re(S*Sm)
        for gs, de in planted:
            w = iv.mpc(iv.mpf(mp.nstr(gs, dps)), iv.mpf(mp.nstr(de, dps)))
            S = iv.mpc(0, 0)
            Sm = iv.mpc(0, 0)
            for k in range(Nuse):
                if v[k] == 0:
                    continue
                vk = iv.mpf(mp.nstr(v[k], dps))
                S += vk * _ivF(k, Lx, w, iv) / nrm
                Sm += vk * _ivF(k, Lx, -w, iv) / nrm
            acc += 4 * (S * Sm).real
        den = iv.mpf(0)
        for k in range(Nuse):
            if v[k] != 0:
                den += iv.mpf(mp.nstr(v[k], dps)) ** 2
        r = acc / den
        # .a/.b come back as degenerate intervals; convert for reporting
        return mp.mpf(r.a), mp.mpf(r.b), bool(r.b < 0)
    finally:
        iv.dps = old


def tail_bound(L, Nuse, v, gamma_K, safety=2):
    r"""Upper bound on the omitted  2 sum_{gamma > gamma_K} F_v(gamma)^2  (>= 0).

    For EVEN_D,  F_k(t) = (-1)^k 2 w_k cos(tL)/(w_k^2 - t^2), so for t > w_max
        |F_v(t)| <= A / (t^2 - w_max^2),   A = (2/sqrt L) sum_k |v_k| w_k .
    Then with the Riemann-von Mangoldt density dN <= (1/2pi) log(t/2pi) dt (times a
    safety factor, since the explicit error term is not carried here):
        tail <= 2 A^2 * safety * int_{gamma_K}^inf (t^2 - w_max^2)^{-2} log(t/2pi)/(2pi) dt.
    T2: uses a density asymptotic, not a proven inequality.
    """
    L = mp.mpf(L)
    ws = [(2 * k + 1) * mp.pi / (2 * L) for k in range(Nuse)]
    wmax = ws[-1]
    if gamma_K <= wmax:
        return mp.inf
    A = 2 * sum(abs(v[k]) * ws[k] for k in range(Nuse)) / mp.sqrt(L)
    integ = mp.quad(lambda t: mp.log(t / (2 * mp.pi)) / (2 * mp.pi) /
                    (t ** 2 - wmax ** 2) ** 2, [gamma_K, 2 * gamma_K, mp.inf])
    return 2 * A ** 2 * safety * integ
