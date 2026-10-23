r"""The finite-window Weil form for L(s, chi), chi real primitive.

Derivation: notes/dirichlet_window.md.  This module is the Phase-1 implementation and is
deliberately a MINIMAL diff from the parent repository's src/zeta_window.py:

  * the cross-correlation machinery (C_terms, Csym_terms, eval_terms, F, F_at_half, freq,
    norm2) is IMPORTED from zeta_window, not reimplemented -- it depends only on the test
    basis, never on the L-function;
  * the archimedean closed form is generalised by exactly two substitutions (D2 of the
    note), and `arch_closed_a(terms, L, 0)` is asserted equal to the parent's
    `arch_closed(terms, L)` in src/selftest.py;
  * the conductor term is kept as TWO matrices, LogPi and LogQ, so that Gate A's sign-flip
    control can flip the conductor alone.

  M_chi = Arch_a - LogPi + LogQ - Prime_chi            (no pole term: chi is non-principal)

Characters: the five real primitive characters of conductor 3,4,5,7,8 are the Kronecker
symbols of the fundamental discriminants -3, -4, 5, -7, 8.
"""
import os
import sys

import mpmath as mp

_PARENT_SRC = os.path.join(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))))), 'src')
if _PARENT_SRC not in sys.path:
    sys.path.insert(0, _PARENT_SRC)

from zeta_window import (EVEN, EVEN_D, ODD, Csym_terms, C_terms, F, F_at_half,  # noqa: E402
                         eval_terms, freq, norm2, von_mangoldt_terms)

# --------------------------------------------------------------------------- characters

#: fundamental discriminant -> (conductor q, parity a with chi(-1) = (-1)^a, label, sign)
CHARS = {
    -3: dict(q=3, a=1, label='chi_-3', tex=r'\chi_{-3}'),
    -4: dict(q=4, a=1, label='chi_-4', tex=r'\chi_{-4}'),
     5: dict(q=5, a=0, label='chi_5',  tex=r'\chi_{5}'),
    -7: dict(q=7, a=1, label='chi_-7', tex=r'\chi_{-7}'),
     8: dict(q=8, a=0, label='chi_8',  tex=r'\chi_{8}'),
}
DISCRIMINANTS = [-3, -4, 5, -7, 8]


def jacobi(a, n):
    """Jacobi symbol (a/n) for odd n >= 1.  Standard binary-gcd style algorithm."""
    if n <= 0 or n % 2 == 0:
        raise ValueError('jacobi: n must be odd and positive')
    a %= n
    t = 1
    while a != 0:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                t = -t
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            t = -t
        a %= n
    return t if n == 1 else 0


def kronecker(a, n):
    """Kronecker symbol (a/n) for n >= 1."""
    if n <= 0:
        raise ValueError('kronecker: this module only needs n >= 1')
    res = 1
    while n % 2 == 0:
        n //= 2
        if a % 2 == 0:
            return 0
        if a % 8 in (3, 5):
            res = -res
    return res * jacobi(a, n)


def chi(d, n):
    """The real primitive character of fundamental discriminant d, as an int in {-1,0,1}."""
    q = CHARS[d]['q']
    n = int(n) % q
    if n == 0:
        return 0
    return kronecker(d, n)


def char_values(d):
    return [chi(d, r) for r in range(CHARS[d]['q'])]


# ------------------------------------------------------------------- L(s, chi) and zeros

def L_chi(s, d, dps=None):
    """L(s, chi_d) = q^{-s} sum_{r=1}^{q-1} chi(r) zeta(s, r/q)  (Hurwitz zeta).

    Valid for all s != 1 and, since sum_r chi(r) = 0 for non-principal chi, the poles of
    the individual Hurwitz zetas cancel, so the combination is entire."""
    old = mp.mp.dps
    if dps:
        mp.mp.dps = dps
    try:
        q = CHARS[d]['q']
        s = mp.mpmathify(s)
        if s == 1:
            # zeta(s,x) = 1/(s-1) - psi(x) + O(s-1) and sum_r chi(r) = 0, so the poles
            # cancel and the limit is exactly -(1/q) sum_r chi(r) psi(r/q).
            tot = mp.mpf(0)
            for r in range(1, q):
                c = chi(d, r)
                if c:
                    tot -= c * mp.digamma(mp.mpf(r) / q)
            return tot / q
        tot = mp.mpf(0)
        for r in range(1, q):
            c = chi(d, r)
            if c:
                tot += c * mp.zeta(s, mp.mpf(r) / q)
        return tot * mp.power(mp.mpf(q), -s)
    finally:
        mp.mp.dps = old


def theta_chi(t, d):
    r"""arg of the completed gamma factor (q/pi)^{(s+a)/2} Gamma((s+a)/2) at s = 1/2 + it."""
    q, a = CHARS[d]['q'], CHARS[d]['a']
    t = mp.mpf(t)
    return (t / 2) * mp.log(mp.mpf(q) / mp.pi) + mp.im(mp.loggamma(
        mp.mpc(mp.mpf(1) / 4 + mp.mpf(a) / 2, t / 2)))


def Z_chi(t, d):
    r"""Hardy-type Z function: Z(t) = e^{i theta(t)} L(1/2 + it, chi), real for real chi.

    Returns (Re, Im); |Im| is a diagnostic that must stay at the precision floor -- it is a
    direct numerical check of the functional equation and of the root number eps(chi) = 1."""
    v = mp.e ** (mp.mpc(0, theta_chi(t, d))) * L_chi(mp.mpc(mp.mpf(1) / 2, t), d)
    return mp.re(v), mp.im(v)


def Z(t, d):
    return Z_chi(t, d)[0]


def N_count(T, d):
    """Riemann-von Mangoldt for L(s,chi): N(T) ~ (T/pi) log(qT/(2 pi e)).  Used only as a
    COMPLETENESS CHECK on the grid search, never as data."""
    q = CHARS[d]['q']
    T = mp.mpf(T)
    if T <= 0:
        return mp.mpf(0)
    return (T / mp.pi) * mp.log(q * T / (2 * mp.pi * mp.e))


# ------------------------------------------------------------- archimedean, generalised

_ARCH_CACHE_A = {}


def _arch_sf_a(alpha, L, is_lin, a):
    r"""Cached special-function values, generalised to the character parity a.

    (D2) of the note: z = (1/4 + a/2) - i alpha/2, and the truncation prefactor
    e^{-U/2} -> e^{-(1/2+a)U}.  Keyed on (dps, L, alpha, kind, a)."""
    key = (mp.mp.dps, mp.nstr(L, 30), mp.nstr(alpha, 30), is_lin, a)
    hit = _ARCH_CACHE_A.get(key)
    if hit is not None:
        return hit
    U = 2 * L
    w = mp.e ** (-2 * U)
    zr = mp.mpf(1) / 4 + mp.mpf(a) / 2
    if is_lin:
        z0 = zr - 1j * alpha / 2
        val = (mp.re(mp.digamma(z0) + mp.euler),
               mp.lerchphi(w, 2, z0), mp.zeta(2, z0), z0)
    else:
        z = zr - 1j * alpha / 2
        val = (mp.digamma(z) + mp.euler, mp.lerchphi(w, 1, z), z)
    _ARCH_CACHE_A[key] = val
    return val


def arch_closed_a(terms, L, a):
    r"""Archimedean term for Re psi((1/4 + a/2) + it/2).

    Arch = -gamma_E C(0) + sum_{m>=0} [ C(0)/(m+1) - 2 int_0^{2L} C(u) e^{-c_m u} du ],
    c_m = 2m + 1/2 + a, summed in closed form via digamma / Lerch Phi / Hurwitz zeta."""
    U = 2 * L
    pref = mp.e ** (-(mp.mpf(1) / 2 + a) * U)
    C0 = eval_terms(terms, L, 0)
    tot = mp.mpf(0)
    for c, alpha, beta in terms:
        if alpha == 'lin':
            dg, lp2, hz, z0 = _arch_sf_a(beta, L, True, a)
            tot += c * U * dg
            tot -= (c / 2) * mp.re(mp.e ** (1j * beta * U) * pref * lp2 - hz)
        else:
            dgg, lp1, z = _arch_sf_a(alpha, L, False, a)
            A = mp.e ** (1j * (beta + alpha * U))
            tot += c * mp.im(mp.e ** (1j * beta) * dgg)
            tot += c * pref * mp.im(A * lp1)
    return -mp.euler * C0 + tot


# --------------------------------------------------------------------- the window matrix

def von_mangoldt_chi(twoL, d):
    """[(n, Lambda(n) * chi(n), log n)] with log n < 2L; entries with chi(n) = 0 dropped."""
    out = []
    for n, lam, logn in von_mangoldt_terms(twoL):
        c = chi(d, n)
        if c:
            out.append((n, c * lam, logn))
    return out


def build_parts_chi(L, idx, dps, d, sector=EVEN_D):
    """(Arch, LogPi, LogQ, Prime, prim) in the orthonormal basis f_k/sqrt(norm2_k).

    M_chi = Arch - LogPi + LogQ - Prime.  `idx` is an arbitrary index set of basis modes
    (the parent repository's build_parts_idx convention)."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        L = mp.mpf(L)
        a = CHARS[d]['a']
        q = CHARS[d]['q']
        idx = list(idx)
        n = len(idx)
        prim = von_mangoldt_chi(2 * L, d)
        nrm = [mp.sqrt(norm2(k, L, sector)) for k in idx]
        Arch = mp.zeros(n, n); LogPi = mp.zeros(n, n)
        LogQ = mp.zeros(n, n); Prime = mp.zeros(n, n)
        for i in range(n):
            for j in range(i, n):
                tm = Csym_terms(idx[i], idx[j], L, sector)
                den = nrm[i] * nrm[j]
                Arch[i, j] = Arch[j, i] = arch_closed_a(tm, L, a) / den
                if i == j:
                    c0 = eval_terms(tm, L, 0) / den
                    LogPi[i, j] = mp.log(mp.pi) * c0
                    LogQ[i, j] = mp.log(mp.mpf(q)) * c0
                pr = mp.mpf(0)
                for nn, lamchi, logn in prim:
                    pr += 2 * lamchi / mp.sqrt(mp.mpf(nn)) * eval_terms(tm, L, logn)
                Prime[i, j] = Prime[j, i] = pr / den
        return Arch, LogPi, LogQ, Prime, prim
    finally:
        mp.mp.dps = old


def build_matrix_chi(L, idx, dps, d, sector=EVEN_D, flip=None):
    """flip in {None, 'logq', 'arch', 'prime'} negates exactly that term (Gate A control)."""
    Arch, LogPi, LogQ, Prime, _ = build_parts_chi(L, idx, dps, d, sector)
    sA = -1 if flip == 'arch' else 1
    sQ = -1 if flip == 'logq' else 1
    sP = -1 if flip == 'prime' else 1
    return sA * Arch - LogPi + sQ * LogQ - sP * Prime


def Lambda_chi(s, d):
    r"""Completed Lambda(s, chi) = (q/pi)^{(s+a)/2} Gamma((s+a)/2) L(s, chi)."""
    q, a = CHARS[d]['q'], CHARS[d]['a']
    s = mp.mpmathify(s)
    return (mp.power(mp.mpf(q) / mp.pi, (s + a) / 2) * mp.gamma((s + a) / 2)
            * L_chi(s, d))
