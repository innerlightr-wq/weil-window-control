r"""Stage 3: the zeta window Weil form on [-L, L], both parity sectors.

Derivations in notes/stage3_assembly.md. Everything closed form; mpmath throughout.

  even sector:  phi_k(x) = cos(a_k x), a_k = (2k+1)pi/(2L), k = 0,1,2,...
  odd  sector:  psi_k(x) = sin(b_k x), b_k = k pi / L,       k = 1,2,3,...

Both are continuous on R (they vanish at +-L) and satisfy \int f_j f_k = L delta_jk.
"""
import mpmath as mp

EVEN, ODD, EVEN_D = 'even', 'odd', 'even_d'
# EVEN  : cos(k pi x / L), k = 0,1,2,...   (full Fourier cosine; contains the constant,
#         no boundary condition -- the fast-converging choice)
# EVEN_D: cos((2k+1) pi x / 2L)            (quarter-wave; vanishes at +-L -- kept as an
#         INDEPENDENT basis for cross-validating lambda*)
# ODD   : sin(k pi x / L), k = 1,2,...


def freq(k, L, sector):
    if sector == EVEN_D:
        return (2 * k + 1) * mp.pi / (2 * L)
    if sector == EVEN:
        return k * mp.pi / L                   # k = 0 gives the constant function
    return (k + 1) * mp.pi / L


def norm2(k, L, sector):
    """\\int_{-L}^{L} f_k^2 dx."""
    if sector == EVEN and k == 0:
        return 2 * L
    return L


def F(k, L, t, sector=EVEN):
    """F_k(t) = \\int f_k(x) e^{itx} dx.  Even sectors: real. Odd: i*(returned value).
    The sinc form is valid for any frequency and is numerically stable."""
    w = freq(k, L, sector)
    s = mp.sincpi((w - t) * L / mp.pi)
    u = mp.sincpi((w + t) * L / mp.pi)
    return L * (s - u) if sector == ODD else L * (s + u)


def F_at_half(k, L, sector=EVEN):
    """P_k = F_k(i/2) = \\int f_k(x) e^{-x/2} dx, real in every sector."""
    w = freq(k, L, sector)
    q = w ** 2 + mp.mpf(1) / 4
    if sector == EVEN_D:                       # cos(wL) = 0, sin(wL) = (-1)^k
        return (-1) ** k * 2 * w * mp.cosh(L / 2) / q
    if sector == EVEN:                         # cos(wL) = (-1)^k, sin(wL) = 0
        return (-1) ** k * mp.sinh(L / 2) / q
    kk = k + 1                                 # ODD: sin(wL)=0, cos(wL)=(-1)^k
    return (-1) ** kk * 2 * w * mp.sinh(L / 2) / q


def pole_sign(sector):
    """h(i/2)+h(-i/2) = 2 F_j(i/2) F_k(-i/2) symmetrised = 2 * sign * P_j P_k."""
    return -1 if sector == ODD else 1


def C_terms(j, k, L, sector):
    """C_{jk}(u) = \\int f_j(x) f_k(x-u) dx on u in [0, 2L], as
       (c, alpha, beta)   meaning  c * sin(alpha u + beta)
       (c, 'lin', gamma)  meaning  c * (2L - u) * cos(gamma u)."""
    wj, wk = freq(j, L, sector), freq(k, L, sector)
    p1, p2 = wj - wk, wj + wk
    s2 = -1 if sector == ODD else 1
    out = []
    if p1 == 0:
        out.append((mp.mpf(1) / 2, 'lin', wk))
    else:
        out.append((1 / (2 * p1), wk, p1 * L))
        out.append((-1 / (2 * p1), wj, -p1 * L))
    if p2 == 0:                                # only for the constant mode (j=k=0, EVEN)
        out.append((s2 * mp.mpf(1) / 2, 'lin', wk))
    else:
        out.append((s2 / (2 * p2), -wk, p2 * L))
        out.append((-s2 / (2 * p2), wj, -p2 * L))
    return out


def Csym_terms(j, k, L, sector):
    """Symmetrised (C_jk + C_kj)/2 -- what the quadratic form actually sees."""
    if j == k:
        return C_terms(j, k, L, sector)
    return ([(c / 2, a, b) for c, a, b in C_terms(j, k, L, sector)] +
            [(c / 2, a, b) for c, a, b in C_terms(k, j, L, sector)])


def eval_terms(terms, L, u):
    u = abs(mp.mpf(u))
    if u >= 2 * L:
        return mp.mpf(0)
    tot = mp.mpf(0)
    for c, a, b in terms:
        tot += c * ((2 * L - u) * mp.cos(b * u) if a == 'lin' else mp.sin(a * u + b))
    return tot


def arch_closed(terms, L):
    """Archimedean term: -gamma C(0) + sum_m [C(0)/(m+1) - 2 \\int_0^{2L} C e^{-c_m u} du],
    with the m-sum done analytically via digamma / Lerch / Hurwitz zeta."""
    U = 2 * L
    w = mp.e ** (-2 * U)
    C0 = eval_terms(terms, L, 0)
    tot = mp.mpf(0)
    for c, alpha, beta in terms:
        if alpha == 'lin':
            z0 = mp.mpf(1) / 4 - 1j * beta / 2
            tot += c * U * mp.re(mp.digamma(z0) + mp.euler)
            tot -= (c / 2) * mp.re(mp.e ** (1j * beta * U) * mp.e ** (-U / 2)
                                   * mp.lerchphi(w, 2, z0) - mp.zeta(2, z0))
        else:
            z = mp.mpf(1) / 4 - 1j * alpha / 2
            A = mp.e ** (1j * (beta + alpha * U))
            tot += c * mp.im(mp.e ** (1j * beta) * (mp.digamma(z) + mp.euler))
            tot += c * mp.e ** (-U / 2) * mp.im(A * mp.lerchphi(w, 1, z))
    return -mp.euler * C0 + tot


def von_mangoldt_terms(twoL):
    nmax = int(mp.floor(mp.e ** twoL)) + 2
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
        if lam is not None and mp.log(n) < twoL:
            out.append((n, lam, mp.log(n)))
    return out


def build_parts(L, N, dps, sector=EVEN):
    """Returns (Pole, Arch, LogPi, Prime) matrices in the ORTHONORMAL basis f_k/sqrt(L).
    The Weil form is  M = Pole + Arch - LogPi - Prime.
    Recorder split (Stage 3d):  A := Pole + Arch - LogPi,  2M_rec := Prime."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        L = mp.mpf(L)
        prim = von_mangoldt_terms(2 * L)
        P = [F_at_half(k, L, sector) for k in range(N)]
        nrm = [mp.sqrt(norm2(k, L, sector)) for k in range(N)]
        sg = pole_sign(sector)
        Pole = mp.zeros(N, N); Arch = mp.zeros(N, N)
        LogPi = mp.zeros(N, N); Prime = mp.zeros(N, N)
        for j in range(N):
            for k in range(j, N):
                tm = Csym_terms(j, k, L, sector)
                d = nrm[j] * nrm[k]
                Pole[j, k] = Pole[k, j] = 2 * sg * P[j] * P[k] / d
                Arch[j, k] = Arch[k, j] = arch_closed(tm, L) / d
                if j == k:
                    LogPi[j, k] = mp.log(mp.pi) * eval_terms(tm, L, 0) / d
                pr = mp.mpf(0)
                for n, lam, logn in prim:
                    pr += 2 * lam / mp.sqrt(mp.mpf(n)) * eval_terms(tm, L, logn)
                Prime[j, k] = Prime[k, j] = pr / d
        return Pole, Arch, LogPi, Prime, prim
    finally:
        mp.mp.dps = old


def build_matrix(L, N, dps, sector=EVEN, verbose=False):
    Pole, Arch, LogPi, Prime, prim = build_parts(L, N, dps, sector)
    if verbose:
        print(f"    sector={sector}  prime powers with log n < 2L: "
              f"{[(n, mp.nstr(l, 6)) for n, l, _ in prim]}")
    return Pole + Arch - LogPi - Prime
