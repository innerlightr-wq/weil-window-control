r"""Suzuki's localized Weil form on the exponential basis, and W_N.

Conventions taken from arXiv:2610.32704v1 and NOT inherited:

    Fourier       uhat(z) = int u(x) e^{izx} dx                           [(2.1)]
    Weil form     Q_W(u,f) = sum_{gamma} m_gamma uhat(gamma) conj(fhat(gamma)),
                  linear in the FIRST variable                           [(2.3)]
    basis         chi_{n,a}(x) = e^{i pi n x/a}/sqrt(2a) on [-a,a]        [(6.1)]
    transform     phi_{n,a}(z) = sqrt(2a) sin(az + n pi)/(az + n pi)      [(6.2)]
    Gram          G_{a,N} = ( Q_W(chi_m, chi_n) )_{m,n=-N..N}             [(6.2)]
    Riesz         G c(z) = phi(zbar)                                      [(6.3)]
    kernel        J_{a,N}(w,z) = phi(z)^T G^{-1} phi(w)                   [(6.4)]
    spectral fn   W_N(a,pi;z) = (z-i) J(-i,z) - (z+i) J(i,z)              [(6.4)]

The geometric side is assembled from the SAME explicit formula this repository already
uses (notes/stage3_assembly.md), specialized to the exponential basis, with every piece
in closed form.  Derivation and the two errors fixed along the way:
notes/normalization.md.  Suzuki cites Yoshida (5.15)/(5.16) for the same object; Yoshida
was NOT obtained, so the assembly here is independent and is reconciled against the
repository's own gated form in Phase 0.

G is computable with no RH input and no zeros of xi.
"""
import mpmath as mp


# --------------------------------------------------------------- basis and transform
def phi(n, a, z):
    """phi_{n,a}(z) = sqrt(2a) sincpi(a z/pi + n); stable at the removable singularity."""
    return mp.sqrt(2 * a) * mp.sincpi(a * z / mp.pi + n)


def chi_corr(m, n, a, v):
    """C_{mn}(v) = int chi_m(x) conj(chi_n(x-v)) dx, exact closed form, 0 for |v| >= 2a."""
    a = mp.mpf(a)
    av = abs(mp.mpf(v))
    if av >= 2 * a:
        return mp.mpc(0)
    L2 = 2 * a - av
    return (mp.e ** (1j * mp.pi * (m + n) * v / (2 * a)) * (L2 / (2 * a))
            * mp.sincpi(mp.mpf(m - n) * L2 / (2 * a)))


def S_corr(m, n, a, v):
    """S(v) = C_{mn}(v) + C_{mn}(-v), REAL, in closed form (notes/normalization.md Sec.3)."""
    a = mp.mpf(a)
    av = abs(mp.mpf(v))
    if av >= 2 * a:
        return mp.mpf(0)
    d = m - n
    if d == 0:
        return 2 * mp.cos(mp.pi * m * v / a) * (2 * a - av) / (2 * a)
    sgn = mp.mpf(-1) ** (d % 2)
    return -(sgn / (mp.pi * d)) * (mp.sin(mp.pi * m * av / a) - mp.sin(mp.pi * n * av / a))


# ------------------------------------------------------------------- arithmetic input
def von_mangoldt(twoa):
    """(n, Lambda(n), log n) for prime powers with log n < 2a.  Cutoff decided in integer
    arithmetic, then re-tested at full precision -- never by a binary float."""
    nmax = int(mp.floor(mp.e ** twoa)) + 2
    out = []
    for n in range(2, nmax + 1):
        m, p, lam = n, 2, None
        while p * p <= m:
            if m % p == 0:
                q = m
                while q % p == 0:
                    q //= p
                lam = mp.log(p) if q == 1 else None
                break
            p += 1
        else:
            lam = mp.log(m)
        if lam is not None and mp.log(n) < twoa:
            out.append((n, lam, mp.log(n)))
    return out


# --------------------------------------------------------------------- the three terms
def pole_term(m, n, a):
    """h(i/2) + h(-i/2) = 2 Re[ phi_m(i/2) phi_n(i/2) ].  See notes/normalization.md Sec.2."""
    h = mp.mpf(1) / 2
    return 2 * mp.re(phi(m, a, 1j * h) * phi(n, a, 1j * h))


def _kern(v):
    """e^{-v/2}/(1-e^{-2v}) = sum_{k>=0} e^{-c_k v}, c_k = 1/2+2k."""
    return mp.e ** (-v / 2) / (-mp.expm1(-2 * v))


def arch_term(m, n, a):
    """The archimedean term in CLOSED FORM -- one analytic 1-D integral, no series tail.

    d != 0:  ((-1)^d/(pi d)) int_0^{2a} [sin(pi m v/a) - sin(pi n v/a)] kern(v) dv
    d == 0:  psi(1/4) + e^{-a} Phi(e^{-4a},1,1/4) - int_0^{2a} [S(v)-2] kern(v) dv

    Both integrands are O(v) against kern's 1/(2v) pole, hence analytic on [0,2a]."""
    a = mp.mpf(a)
    U = 2 * a
    d = m - n
    if d == 0:
        def f(v):
            # S(v) - 2 = 2[cos(th)(1 - v/(2a)) - 1] = -4 sin^2(th/2) - (v/a)cos(th),
            # th = pi m v/a.  The naive difference of two numbers near 2 loses ~|log v|
            # digits: at v = 1e-55 and dps 60 it is wrong by a RELATIVE 2.65e-7, and
            # tanh-sinh quadrature samples v ~ 1e-(2 dps).  The half-angle form has no
            # cancellation.  Measured in REPORT (Retraction 24).
            if v == 0:
                return mp.mpf(0)
            th = mp.pi * m * v / a
            Sm2 = -4 * mp.sin(th / 2) ** 2 - (v / a) * mp.cos(th)
            return Sm2 * _kern(v)
        I = mp.quad(f, [0, U / 8, U / 2, U])
        return (mp.digamma(mp.mpf(1) / 4)
                + mp.e ** (-a) * mp.lerchphi(mp.e ** (-4 * a), 1, mp.mpf(1) / 4) - I)
    sgn = mp.mpf(-1) ** (d % 2)

    def g(v):
        # sin A - sin B = 2 cos((A+B)/2) sin((A-B)/2); exact, and no cancellation
        # (verified: both forms agree to 1.3e-61 at v = 1e-55, dps 60).
        if v == 0:
            return mp.mpf(0)
        return (2 * mp.cos(mp.pi * (m + n) * v / (2 * a))
                * mp.sin(mp.pi * d * v / (2 * a))) * _kern(v)
    I = mp.quad(g, [0, U / 8, U / 2, U])
    return (sgn / (mp.pi * d)) * I


def prime_term(m, n, a, prim):
    """- sum_n Lambda(n) n^{-1/2} S(log n).  Finite sum; S vanishes at |v| = 2a."""
    tot = mp.mpf(0)
    for nn, lam, logn in prim:
        tot += lam / mp.sqrt(mp.mpf(nn)) * S_corr(m, n, a, logn)
    return -tot


def gram_entry(m, n, a, prim):
    """Q_W(chi_m, chi_n).  No RH, no zeros of xi."""
    return (pole_term(m, n, a) + arch_term(m, n, a)
            - mp.log(mp.pi) * mp.re(chi_corr(m, n, a, 0)) + prime_term(m, n, a, prim))


def gram(a, N, dps=60, progress=False):
    """G_{a,N}, (2N+1)x(2N+1), real symmetric.  Reality and the index-reversal symmetry
    G_{-m,-n} = G_{m,n} are DERIVED (notes/normalization.md Sec.3) and re-checked in
    Phase 0 test P0-T4 -- the arithmetic assembly must exhibit them with no reference to
    any zero of xi.

    Only entries with m <= n and m+n >= 0 are computed; the rest follow from symmetry
    and index reversal, which is a 4x saving and makes N = 18 at 180 digits feasible."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        a = mp.mpf(a)
        prim = von_mangoldt(2 * a)
        idx = list(range(-N, N + 1))
        pos = {v: i for i, v in enumerate(idx)}
        G = mp.zeros(len(idx), len(idx))
        done = set()
        for m in idx:
            for n in idx:
                if (m, n) in done:
                    continue
                val = gram_entry(m, n, a, prim)
                for (p, q) in ((m, n), (n, m), (-m, -n), (-n, -m)):
                    G[pos[p], pos[q]] = val
                    done.add((p, q))
            if progress:
                print(f"      row m={m} done", flush=True)
        return G
    finally:
        mp.mp.dps = old


# --------------------------------------------------------------------------- W_N
def W_N(G, a, N, z, brief_convention=False):
    """W_N(a,pi;z) = (z-i) J(-i,z) - (z+i) J(i,z).

    paper (6.4):     J(w,z) = phi(z)^T G^{-1} phi(w)
    brief variant:   J(w,z) = phi(z)^T G^{-1} phi(wbar)   [brief_convention=True]

    They differ by w -> wbar, which swaps the two terms.  Both make -i W_N real on R, so
    reality does not discriminate; Phase 1 test P1-T5 decides empirically."""
    idx = list(range(-N, N + 1))
    pz = mp.matrix([phi(n, a, z) for n in idx])
    wm, wp = (-1j, 1j)
    if brief_convention:
        wm, wp = (1j, -1j)
    vm = mp.lu_solve(G, mp.matrix([phi(n, a, wm) for n in idx]))
    vp = mp.lu_solve(G, mp.matrix([phi(n, a, wp) for n in idx]))
    Jm = sum(pz[i] * vm[i] for i in range(len(idx)))
    Jp = sum(pz[i] * vp[i] for i in range(len(idx)))
    return (z - 1j) * Jm - (z + 1j) * Jp


def W_N_factory(G, a, N, brief_convention=False):
    """Precompute the two Gram solves; return z -> -i W_N(z), real on R."""
    idx = list(range(-N, N + 1))
    wm, wp = (-1j, 1j)
    if brief_convention:
        wm, wp = (1j, -1j)
    vm = mp.lu_solve(G, mp.matrix([phi(n, a, wm) for n in idx]))
    vp = mp.lu_solve(G, mp.matrix([phi(n, a, wp) for n in idx]))

    def f(z):
        pz = [phi(n, a, z) for n in idx]
        Jm = sum(pz[i] * vm[i] for i in range(len(idx)))
        Jp = sum(pz[i] * vp[i] for i in range(len(idx)))
        return mp.re(-1j * ((z - 1j) * Jm - (z + 1j) * Jp))
    return f
