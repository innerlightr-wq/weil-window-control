r"""The Davenport--Heilbronn function F(s), its completion, and its Weil window form.

Conventions are fixed once here and used by every other script in this extension.

CHARACTER.  chi is the character mod 5 with

    chi(0..4) = (0, 1, i, -i, -1),

primitive and ODD (chi(-1) = chi(4) = -1), so the gamma factor is Gamma((s+1)/2).

THE FUNCTION.

    F(s) = [ (1 - i kappa) L(s, chi) + (1 + i kappa) L(s, conj chi) ] / 2,

with kappa DERIVED (not imported as a decimal) from the root number: writing
eps_chi = tau(chi)/(i sqrt 5), the completed functional equation for F holds iff
(1 - i kappa) eps_chi = 1 + i kappa, i.e. eps_chi = (1 + i kappa)/(1 - i kappa), i.e.

    kappa = tan( arg(eps_chi) / 2 ).

See notes/phaseA.md.  The closed form is kappa = [sqrt(10 - 2 sqrt 5) - 2]/(sqrt 5 - 1).

DIRICHLET COEFFICIENTS.  a_n depends only on n mod 5:

    (a_1, a_2, a_3, a_4, a_5) = (1, kappa, -kappa, -1, 0),

real and periodic but NOT multiplicative (a_2 a_3 = -kappa^2 while a_6 = a_1 = 1), which
is why F has no Euler product.

COMPLETION.

    Lambda_F(s) = (5/pi)^{(s+1)/2} Gamma((s+1)/2) F(s),      Lambda_F(s) = Lambda_F(1-s).

This is SELF-reciprocal; the individual L(s,chi) equation exchanges chi and conj chi.

WEIL TRANSFORM CONVENTION (matches the parent repository's src/zeta_window.py).
g real, supported in [-L, L];

    G(z) = int g(u) e^{zu} du,   h(v) = int g(u+v) g(u) du,   H(z) = G(z) G(-z),

so H(it) = |G(it)|^2 for real g, and h = g cross-correlated with itself.  The spectral
form is over the zeros of the ENTIRE COMPLETED function,

    Q_F(g) = sum_rho mult(rho) H(rho - 1/2).

The arithmetic form derived in notes/DERIVATION.md (reproduced in DERIVATION.md) is

    Q_F(g) = h(0) log(5/pi)
             + (1/2pi) int |G(it)|^2 Re psi(3/4 + it/2) dt
             - 2 sum_{2 <= n <= e^{2L}} b_n h(log n) / sqrt(n),

with -F'/F(s) = sum_{n>=2} b_n n^{-s}.  There is NO pole term: Lambda_F is entire.
"""
import os
import sys

import mpmath as mp

_HERE = os.path.dirname(os.path.abspath(__file__))
_EXT = os.path.dirname(_HERE)
_ROOT = os.path.dirname(os.path.dirname(_EXT))
_PARENT_SRC = os.path.join(_ROOT, 'src')
_DIRICHLET_SRC = os.path.join(_ROOT, 'extensions', 'dirichlet-windows-2026-10-23', 'src')
for _p in (_PARENT_SRC, _DIRICHLET_SRC):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from zeta_window import (EVEN, EVEN_D, ODD, Csym_terms, eval_terms,  # noqa: E402
                         freq, norm2)
from dirichlet_window import arch_closed_a                            # noqa: E402

MPF_HI_DPS = 220    # precision at which every decimal CONSTANT is parsed


def mpf_hi(x, dps=None):
    r"""Parse a decimal constant at a FIXED high precision, independent of the ambient
    mp.mp.dps.

    mpmath's mp.mpf('1.8') rounds the decimal to the CURRENT working precision and keeps
    that value forever.  A window length built while mp.mp.dps was 15 and then used in a
    computation at 160 digits silently carries a 1e-17 error, and a precision ladder run
    on such a value compares two evaluations of the same wrong number and agrees to 85
    digits while being wrong in the 13th.  Every decimal constant in this extension --
    window lengths, the located zero, perturbation parameters -- goes through this
    function.  (Binary-exact decimals such as 1.0, 2.0, 0.5 are unaffected either way;
    1.2, 1.8, 2.6 are not.)  See RETRACTIONS.md R6.
    """
    old = mp.mp.dps
    mp.mp.dps = dps or MPF_HI_DPS
    try:
        return mp.mpf(x)
    finally:
        mp.mp.dps = old


Q_COND = 5          # conductor
A_PARITY = 1        # chi(-1) = (-1)^a


# ----------------------------------------------------------------- character and kappa

_CHI = {0: None, 1: (1, 0), 2: (0, 1), 3: (0, -1), 4: (-1, 0)}


def chi(n):
    """chi(n) as an mpc; chi(0..4) = (0, 1, i, -i, -1)."""
    v = _CHI[int(n) % 5]
    return mp.mpc(0, 0) if v is None else mp.mpc(v[0], v[1])


def gauss_sum():
    """tau(chi) = sum_{n=1}^{4} chi(n) e^{2 pi i n / 5}."""
    return mp.fsum([chi(n) * mp.e ** (2j * mp.pi * n / 5) for n in range(1, 5)])


def root_number():
    """eps_chi = tau(chi) / (i sqrt q) for odd primitive chi."""
    return gauss_sum() / (1j * mp.sqrt(mp.mpf(Q_COND)))


def kappa():
    """kappa = tan(arg(eps_chi)/2), derived from the Gauss sum.  Cached per precision."""
    key = mp.mp.dps
    hit = _KAPPA_CACHE.get(key)
    if hit is None:
        hit = mp.tan(mp.arg(root_number()) / 2)
        _KAPPA_CACHE[key] = hit
    return hit


_KAPPA_CACHE = {}


def kappa_closed_form():
    """[sqrt(10 - 2 sqrt 5) - 2] / (sqrt 5 - 1).  Checked equal to kappa() in tests."""
    return (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)


def a_coeff(n):
    """a_n = Re chi(n) + kappa Im chi(n); real, period 5."""
    k = kappa()
    r = int(n) % 5
    return (mp.mpf(1) if r == 1 else -mp.mpf(1) if r == 4 else
            k if r == 2 else -k if r == 3 else mp.mpf(0))


def a_vector():
    """(a_1, ..., a_5)."""
    return [a_coeff(r) for r in (1, 2, 3, 4, 5)]


# ------------------------------------------------------------------ evaluation routes

def _hurwitz_combo(s, coeffs):
    r"""5^{-s} sum_{r=1}^{5} coeffs[r-1] zeta(s, r/5), with the s = 1 limit handled.

    sum_r coeffs[r] = 0 is required for entireness; the caller guarantees it."""
    s = mp.mpmathify(s)
    if s == 1:
        # zeta(s, x) = 1/(s-1) - psi(x) + O(s-1); the 1/(s-1) poles cancel.
        return -mp.fsum([coeffs[r - 1] * mp.digamma(mp.mpf(r) / 5)
                         for r in range(1, 6)]) / 5
    return mp.power(mp.mpf(5), -s) * mp.fsum(
        [coeffs[r - 1] * mp.zeta(s, mp.mpf(r) / 5) for r in range(1, 6)])


def L_chi(s, conj=False):
    """L(s, chi) (or L(s, conj chi)) by the Hurwitz-zeta character sum."""
    cs = [mp.conj(chi(r)) if conj else chi(r) for r in range(1, 6)]
    return _hurwitz_combo(s, cs)


def F_route_A(s):
    """Route A: the Dirichlet L-function combination."""
    k = kappa()
    return ((1 - 1j * k) * L_chi(s) + (1 + 1j * k) * L_chi(s, conj=True)) / 2


def F_route_B(s):
    """Route B: the periodic-coefficient Hurwitz expression, coefficients (1,k,-k,-1,0)."""
    return _hurwitz_combo(s, a_vector())


def _hurwitz_EM(s, a, N=None, M=None):
    r"""Hurwitz zeta by Euler--Maclaurin, implemented here and NOT calling mp.zeta.

        zeta(s,a) = sum_{n=0}^{N-1} (n+a)^{-s} + (N+a)^{1-s}/(s-1) + (N+a)^{-s}/2
                    + sum_{k=1}^{M} B_{2k}/(2k)! * (s)_{2k-1} * (N+a)^{-s-2k+1},

    where (s)_{2k-1} = s(s+1)...(s+2k-2).  This is the independent evaluation route: the
    only mpmath primitives it uses are power, bernoulli and factorial."""
    s = mp.mpmathify(s)
    a = mp.mpf(a)
    if s == 1:
        raise ValueError('_hurwitz_EM: pole at s = 1')
    if M is None:
        M = mp.mp.dps // 2 + 6
    if N is None:
        N = max(12, int(mp.mp.dps * 1.6) + int(abs(mp.im(s))) + 10)
    tot = mp.fsum([mp.power(n + a, -s) for n in range(N)])
    x = mp.mpf(N) + a
    tot += mp.power(x, 1 - s) / (s - 1) + mp.power(x, -s) / 2
    poch = s                                    # (s)_1
    for k in range(1, M + 1):
        tot += (mp.bernoulli(2 * k) / mp.factorial(2 * k) * poch
                * mp.power(x, -s - 2 * k + 1))
        poch *= (s + 2 * k - 1) * (s + 2 * k)    # (s)_{2k+1}
    return tot


def F_route_C(s):
    r"""Route C: the periodic-coefficient formula with the Hurwitz zeta supplied by
    `_hurwitz_EM` instead of `mp.zeta`.  Shares no special-function library code with
    routes A and B; the two still share the ALGEBRAIC identity, which is why route D
    exists as well."""
    s = mp.mpmathify(s)
    if s == 1:
        return F_route_B(s)
    av = a_vector()
    return mp.power(mp.mpf(5), -s) * mp.fsum(
        [av[r - 1] * _hurwitz_EM(s, mp.mpf(r) / 5) for r in range(1, 6)])


def F_route_D(s, nmax=20000):
    r"""Route D: direct summation of sum a_n n^{-s}.  Only usable for Re s >= 3, where a
    few thousand terms give ~12 digits; it is a sanity check on the NORMALISATION of the
    Dirichlet coefficients, not a precision route.  Tail bound: nmax^{1-sig}/(sig-1)."""
    s = mp.mpmathify(s)
    sig = mp.re(s)
    if sig < 3:
        raise ValueError('F_route_D needs Re s >= 3')
    tot = mp.mpf(0)
    for n in range(1, int(nmax) + 1):
        c = a_coeff(n)
        if c != 0:
            tot += c * mp.power(mp.mpf(n), -s)
    return tot


F = F_route_B          # the working route


def Lambda_F(s):
    """(5/pi)^{(s+1)/2} Gamma((s+1)/2) F(s)."""
    s = mp.mpmathify(s)
    return (mp.power(mp.mpf(Q_COND) / mp.pi, (s + 1) / 2)
            * mp.gamma((s + 1) / 2) * F(s))


def theta_DH(t):
    r"""arg of the gamma factor at s = 1/2 + it: (t/2) log(5/pi) + Im loggamma(3/4 + it/2)."""
    t = mp.mpf(t)
    return (t / 2) * mp.log(mp.mpf(Q_COND) / mp.pi) + mp.im(
        mp.loggamma(mp.mpc(mp.mpf(3) / 4, t / 2)))


def Z_DH(t):
    r"""Hardy-type Z: e^{i theta(t)} F(1/2 + it).  Real for real t because Lambda_F is
    self-reciprocal and F has real Dirichlet coefficients, so Lambda_F(conj s) =
    conj Lambda_F(s) and Lambda_F(1/2+it) is real.  Returns (Re, Im)."""
    v = mp.e ** (mp.mpc(0, theta_DH(t))) * F(mp.mpc(mp.mpf(1) / 2, t))
    return mp.re(v), mp.im(v)


def Z(t):
    return Z_DH(t)[0]


def N_count(T):
    """Riemann--von Mangoldt MAIN TERM for the ONE-SIDED count #{rho : 0 < Im rho < T}
    of zeros of Lambda_F:  (T/2pi) log(5T/(2 pi e)).

    (The two-sided count over |Im rho| < T is twice this.)  Used ONLY as a completeness
    diagnostic, never as data; the actual counts come from the argument principle."""
    T = mp.mpf(T)
    if T <= 0:
        return mp.mpf(0)
    return (T / (2 * mp.pi)) * mp.log(Q_COND * T / (2 * mp.pi * mp.e))


def dF(s):
    """F'(s), from the Hurwitz representation differentiated term by term:
       F'(s) = -log(5) F(s) + 5^{-s} sum_r a_r d/ds zeta(s, r/5)."""
    s = mp.mpmathify(s)
    av = a_vector()
    d = mp.fsum([av[r - 1] * mp.zeta(s, mp.mpf(r) / 5, 1) for r in range(1, 6)])
    return -mp.log(mp.mpf(5)) * F(s) + mp.power(mp.mpf(5), -s) * d


def dlogLambda(s):
    """Lambda_F'/Lambda_F (s) = (1/2) log(5/pi) + (1/2) psi((s+1)/2) + F'/F (s)."""
    s = mp.mpmathify(s)
    return (mp.log(mp.mpf(Q_COND) / mp.pi) / 2
            + mp.digamma((s + 1) / 2) / 2 + dF(s) / F(s))


# ------------------------------------------------------------ logarithmic-derivative b_n

def b_coeffs(nmax):
    r"""b_n for 2 <= n <= nmax, from  sum_{d|n} b_d a_{n/d} = a_n log n,  b_1 = 0.

    Returns a list indexed 0..nmax with b[0] = b[1] = 0."""
    nmax = int(nmax)
    a = [mp.mpf(0)] * (nmax + 1)
    for n in range(1, nmax + 1):
        a[n] = a_coeff(n)
    b = [mp.mpf(0)] * (nmax + 1)
    for n in range(2, nmax + 1):
        tot = a[n] * mp.log(n)
        d = 2
        while d * d <= n:
            if n % d == 0:
                tot -= b[d] * a[n // d]
                e = n // d
                if e != d and e < n:
                    tot -= b[e] * a[d]
            d += 1
        b[n] = tot
    return b


def b_coeffs_via_inverse(nmax):
    r"""Second, independent route: b = (a log) * a^{-1} under Dirichlet convolution.

    Build the Dirichlet inverse c = a^{-1} (c_1 = 1, c_n = -sum_{d|n, d<n} c_d a_{n/d}),
    then b_n = sum_{d|n} (a_d log d) c_{n/d}."""
    nmax = int(nmax)
    a = [mp.mpf(0)] * (nmax + 1)
    for n in range(1, nmax + 1):
        a[n] = a_coeff(n)
    c = [mp.mpf(0)] * (nmax + 1)
    c[1] = mp.mpf(1)
    for n in range(2, nmax + 1):
        tot = mp.mpf(0)
        for d in _divisors(n):
            if d < n:
                tot -= c[d] * a[n // d]
        c[n] = tot
    b = [mp.mpf(0)] * (nmax + 1)
    for n in range(2, nmax + 1):
        tot = mp.mpf(0)
        for d in _divisors(n):
            if d >= 2:
                tot += a[d] * mp.log(d) * c[n // d]
        b[n] = tot
    return b


def _divisors(n):
    out = []
    d = 1
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            if d != n // d:
                out.append(n // d)
        d += 1
    return sorted(out)


def b_terms(twoL, b=None):
    """[(n, b_n, log n)] for 2 <= n with log n < 2L.  The cutoff is decided in INTEGER
    arithmetic (n <= floor(e^{2L}) then log n < 2L re-checked at full precision), never
    by a float comparison."""
    twoL = mp.mpf(twoL)
    nmax = int(mp.floor(mp.e ** twoL)) + 2
    if b is None or len(b) <= nmax:
        b = b_coeffs(nmax)
    out = []
    for n in range(2, nmax + 1):
        ln = mp.log(n)
        if ln < twoL and b[n] != 0:
            out.append((n, b[n], ln))
    return out


# --------------------------------------------------------------------- the window matrix

def build_parts(L, idx, dps, sector=EVEN_D):
    """(Arch, LogPi, LogQ, Prime, terms) in the orthonormal basis f_k / sqrt(norm2_k).

        M_DH = Arch - LogPi + LogQ - Prime            (no pole term)

    `idx` is an arbitrary index set of basis modes (parent build_parts_idx convention)."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        L = mp.mpf(L)
        idx = list(idx)
        n = len(idx)
        terms = b_terms(2 * L)
        nrm = [mp.sqrt(norm2(k, L, sector)) for k in idx]
        Arch = mp.zeros(n, n); LogPi = mp.zeros(n, n)
        LogQ = mp.zeros(n, n); Prime = mp.zeros(n, n)
        for i in range(n):
            for j in range(i, n):
                tm = Csym_terms(idx[i], idx[j], L, sector)
                den = nrm[i] * nrm[j]
                Arch[i, j] = Arch[j, i] = arch_closed_a(tm, L, A_PARITY) / den
                if i == j:
                    c0 = eval_terms(tm, L, 0) / den
                    LogPi[i, j] = mp.log(mp.pi) * c0
                    LogQ[i, j] = mp.log(mp.mpf(Q_COND)) * c0
                pr = mp.mpf(0)
                for nn, bn, logn in terms:
                    pr += 2 * bn / mp.sqrt(mp.mpf(nn)) * eval_terms(tm, L, logn)
                Prime[i, j] = Prime[j, i] = pr / den
        return Arch, LogPi, LogQ, Prime, terms
    finally:
        mp.mp.dps = old


def build_matrix(L, idx, dps, sector=EVEN_D, drop_n=()):
    """M_DH on the span of `idx`.  `drop_n` deliberately omits those n from the prime sum
    (Gate control B)."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        L = mp.mpf(L)
        idx = list(idx)
        n = len(idx)
        terms = [t for t in b_terms(2 * L) if t[0] not in set(drop_n)]
        nrm = [mp.sqrt(norm2(k, L, sector)) for k in idx]
        M = mp.zeros(n, n)
        logqpi = mp.log(mp.mpf(Q_COND) / mp.pi)
        for i in range(n):
            for j in range(i, n):
                tm = Csym_terms(idx[i], idx[j], L, sector)
                den = nrm[i] * nrm[j]
                v = arch_closed_a(tm, L, A_PARITY) / den
                if i == j:
                    v += logqpi * eval_terms(tm, L, 0) / den
                pr = mp.mpf(0)
                for nn, bn, logn in terms:
                    pr += 2 * bn / mp.sqrt(mp.mpf(nn)) * eval_terms(tm, L, logn)
                v -= pr / den
                M[i, j] = M[j, i] = v
        return M
    finally:
        mp.mp.dps = old


# ----------------------------------------------------------------- transforms of a test g

def G_at(coeffs, idx, L, sector, z):
    r"""G(z) = int g(u) e^{zu} du for g = sum_m coeffs[m] f_{idx[m]}(u)/sqrt(norm2).

    For f = cos(w u) on [-L, L]:  int cos(wu) e^{zu} du = 2 z cosh(zL)... -- we use the
    stable closed forms

        int_{-L}^{L} cos(wu) e^{zu} du = 2 [ z cos(wL) sinh(zL) + w sin(wL) cosh(zL) ]
                                          / (z^2 + w^2)
        int_{-L}^{L} sin(wu) e^{zu} du = 2 [ z sin(wL) cosh(zL) - w cos(wL) sinh(zL) ]
                                          / (z^2 + w^2)

    both of which are derived in DERIVATION.md and unit-tested against quadrature."""
    L = mp.mpf(L)
    z = mp.mpmathify(z)
    tot = mp.mpc(0)
    for m, k in enumerate(idx):
        c = coeffs[m]
        if c == 0:
            continue
        w = freq(k, L, sector)
        den = z ** 2 + w ** 2
        if den == 0:
            # removable singularity at z = +- i w, evaluated in closed form:
            #   int cos(wu) e^{+-iwu} du = L + sin(2wL)/(2w)
            #   int sin(wu) e^{+-iwu} du = +- i [ L - sin(2wL)/(2w) ]
            half = mp.sin(2 * w * L) / (2 * w)
            if sector == ODD:
                val = mp.mpc(0, 1) * (L - half) * (1 if mp.im(z) > 0 else -1)
            else:
                val = mp.mpc(L + half, 0)
            tot += c * val / mp.sqrt(norm2(k, L, sector))
            continue
        if sector == ODD:
            num = 2 * (z * mp.sin(w * L) * mp.cosh(z * L)
                       - w * mp.cos(w * L) * mp.sinh(z * L))
        else:
            num = 2 * (z * mp.cos(w * L) * mp.sinh(z * L)
                       + w * mp.sin(w * L) * mp.cosh(z * L))
        tot += c * num / den / mp.sqrt(norm2(k, L, sector))
    return tot


def H_at(coeffs, idx, L, sector, z):
    """H(z) = G(z) G(-z)."""
    return (G_at(coeffs, idx, L, sector, z)
            * G_at(coeffs, idx, L, sector, -z))


def quartet_contribution(coeffs, idx, L, sector, beta, gamma, mult=1):
    r"""Contribution of the COMPLETE orbit of a non-real off-line zero beta + i gamma:

        O = {beta +- i gamma, 1 - beta +- i gamma},

    i.e. sum over the four shifts rho - 1/2 in {+-delta +- i gamma}, delta = beta - 1/2.
    With A = G(delta + i gamma), B = G(-delta + i gamma) and g real,

        Q_O = 2 [ H(delta + i gamma) + H(-delta + i gamma) ]
            = 4 Re( A conj(B) )
            = 4 ( |P|^2 - |Q|^2 ),
        P = int g cosh(delta u) e^{i gamma u} du,  Q = int g sinh(delta u) e^{i gamma u} du.

    This is NOT |G(delta + i gamma)|^2.  For gamma = 0 the orbit degenerates to the real
    pair {beta, 1 - beta} and the caller must use `real_pair_contribution` instead."""
    d = mp.mpf(beta) - mp.mpf(1) / 2
    g = mp.mpf(gamma)
    z1 = mp.mpc(d, g)
    z2 = mp.mpc(-d, g)
    val = 2 * (H_at(coeffs, idx, L, sector, z1) + H_at(coeffs, idx, L, sector, z2))
    return mult * mp.re(val)


def real_pair_contribution(coeffs, idx, L, sector, beta, mult=1):
    """Contribution of the real off-line pair {beta, 1 - beta} (gamma = 0): two zeros, so
    H(delta) + H(-delta) = 2 H(delta) since H is even."""
    d = mp.mpf(beta) - mp.mpf(1) / 2
    return mult * 2 * mp.re(H_at(coeffs, idx, L, sector, d))


def online_pair_contribution(coeffs, idx, L, sector, gamma, mult=1):
    """Contribution of the on-line conjugate pair {1/2 +- i gamma}: 2 H(i gamma)."""
    g = mp.mpf(gamma)
    return mult * 2 * mp.re(H_at(coeffs, idx, L, sector, mp.mpc(0, g)))


def arch_binet(coeffs, idx, L, sector, dps, x0=None):
    r"""Archimedean term by Gauss's integral for psi -- a representation that shares no
    code and no series with the matrix route's digamma/Lerch closed form.

        psi(z) = int_0^inf [ e^{-x}/x - e^{-zx}/(1 - e^{-x}) ] dx       (Re z > 0)

    At z = 3/4 + it/2 and after (1/2pi) int |G(it)|^2 (.) dt, using
    (1/2pi) int |G(it)|^2 dt = h(0) and (1/2pi) int |G(it)|^2 cos(tx/2) dt = h(x/2):

        Arch = int_0^inf [ h(0) e^{-x}/x - e^{-3x/4} h(x/2) / (1 - e^{-x}) ] dx.

    h(v) vanishes for |v| > 2L, so for x > 4L the integrand is h(0) e^{-x}/x and that
    tail is EXACTLY h(0) E_1(4L).  The integral is one-dimensional and NON-oscillatory,
    unlike the frequency-side form.

    Returns (value, cut_bound) where cut_bound is a rigorous bound on the omitted piece
    over [0, x0]; the quadrature on [x0, 4L] is a library estimate and is reported as
    such by the caller."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        L = mp.mpf(L)
        if x0 is None:
            x0 = mp.mpf(10) ** (-(dps // 3))
        h0 = _h_scalar(coeffs, idx, L, sector, 0)

        def integrand(x):
            return (h0 * mp.e ** (-x) / x
                    - mp.e ** (-3 * x / 4) * _h_scalar(coeffs, idx, L, sector, x / 2)
                    / (1 - mp.e ** (-x)))

        # [0, x0]: the integrand is continuous with value -3 h(0)/4 at 0, so the omitted
        # piece is at most x0 * sup, and we bound sup on a grid plus a safety factor 2.
        sup = max(abs(integrand(x0 * mp.mpf(j) / 8)) for j in range(1, 9))
        cut = 2 * x0 * max(sup, abs(h0))
        # [x0, 4L] by quadrature, subdivided at the one interior kink of h (v = 0 -> x = 0
        # is the endpoint; h is smooth on (0, 2L)), plus the exact tail.
        main = mp.quad(integrand, [x0, mp.mpf(1), 4 * L])
        tail = h0 * mp.e1(4 * L)
        return main + tail, cut
    finally:
        mp.mp.dps = old


def _h_scalar(coeffs, idx, L, sector, v):
    tot = mp.mpf(0)
    for i, ki in enumerate(idx):
        for j, kj in enumerate(idx):
            if coeffs[i] == 0 or coeffs[j] == 0:
                continue
            tm = Csym_terms(ki, kj, L, sector)
            den = mp.sqrt(norm2(ki, L, sector) * norm2(kj, L, sector))
            tot += coeffs[i] * coeffs[j] * eval_terms(tm, L, v) / den
    return tot
