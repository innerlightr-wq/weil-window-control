r"""Weil quadratic form Q(f,g) for ARBITRARY compactly supported piecewise-cosine f, g.

The project's zeta_window.py builds Q only in the basis cos(k pi x / L) on [-L, L], which
is NOT nested as L grows.  To test resolution inheritance we need a basis of H_{L'} whose
first block spans H_L, so we need the form on arbitrary pieces.

Same functional as the project (notes/stage3_assembly.md):

    Q(f,g) = 2 sgn P_f P_g + Arch[C] - log(pi) C(0) - sum_n (2 Lambda(n)/sqrt n) C(log n)

with C = symmetrised cross-correlation, P_f = int f(x) e^{-x/2} dx, and

    Arch[C] = -gamma C(0) + sum_{m>=0} [ C(0)/(m+1) - 2 int_0^U C(u) e^{-c_m u} du ],
    c_m = 2m + 1/2.

The m-sum is resummed in closed form here.  Writing C(u) = C(0) + (C(u)-C(0)),

    sum_m [ C(0)/(m+1) - C(0)/(m+1/4) ]        = C(0) (psi(1/4) + gamma)
    sum_m  2 C(0) e^{-c_m U}/c_m               geometric, ratio e^{-2U}
    sum_m  2 int_0^U (C-C(0)) e^{-c_m u} du    = 2 int_0^U (C(u)-C(0)) R(u) du,
                                                 R(u) = e^{-u/2}/(1 - e^{-2u}),

so  Arch[C] = C(0) psi(1/4) + sum_m 2 C(0) e^{-c_m U}/c_m - 2 int_0^U (C-C(0)) R du.
R(u) ~ 1/(2u) at 0 and C(u)-C(0) = O(u), so the integrand is bounded.
"""
import mpmath as mp

# a basis function is a list of pieces; piece = (x0, a, p, q) meaning cos(a(x-x0)) on [p,q]


def _cc_int(a, x0, b, y0, lo, hi):
    """int_lo^hi cos(a(x-x0)) cos(b(x-y0)) dx, exact."""
    tot = mp.mpf(0)
    for s in (1, -1):
        w = a - s * b
        ph = -a * x0 + s * b * y0
        if w == 0:
            tot += mp.cos(ph) * (hi - lo) / 2
        else:
            tot += (mp.sin(w * hi + ph) - mp.sin(w * lo + ph)) / (2 * w)
    return tot


def corr(fp, gp, u):
    """C_fg(u) = int f(x) g(x-u) dx."""
    tot = mp.mpf(0)
    for (x0, a, p, q) in fp:
        for (y0, b, r, s) in gp:
            lo, hi = max(p, r + u), min(q, s + u)
            if hi > lo:
                tot += _cc_int(a, x0, b, y0 + u, lo, hi)
    return tot


def csym(fp, gp, u):
    """Symmetrised (C_fg(u) + C_gf(u))/2; C_gf(u) = C_fg(-u).  Even in u."""
    return (corr(fp, gp, u) + corr(fp, gp, -u)) / 2


def breakpoints(fp, gp, U):
    """u where the overlap endpoints switch -- C is only piecewise smooth there."""
    bs = set()
    for (_, _, p, q) in fp:
        for (_, _, r, s) in gp:
            # C_sym(u) = (C(u) + C(-u))/2, so the kinks of BOTH halves matter: the
            # negative-u kinks land at the negatives of these.  Omitting them leaves
            # undeclared kinks inside the quadrature panels and mp.quad grinds to
            # maxdegree without converging.
            for v in (p - r, p - s, q - r, q - s):
                for w in (v, -v):
                    if 0 < w < U:
                        bs.add(mp.mpf(w))
    return sorted(bs)


def p_half(fp):
    """P_f = int f(x) e^{-x/2} dx."""
    tot = mp.mpf(0)
    for (x0, a, p, q) in fp:
        z = 1j * a - mp.mpf(1) / 2
        tot += mp.re(mp.e ** (-1j * a * x0) * (mp.e ** (z * q) - mp.e ** (z * p)) / z)
    return tot


def l2(fp, gp):
    """<f,g> = C_fg(0)."""
    return corr(fp, gp, mp.mpf(0))


def arch(fp, gp, U, maxdeg=8):
    C0 = csym(fp, gp, mp.mpf(0))
    tot = C0 * mp.digamma(mp.mpf(1) / 4)
    m = 0                                   # geometric tail, ratio e^{-2U}
    while True:
        cm = 2 * m + mp.mpf(1) / 2
        t = 2 * C0 * mp.e ** (-cm * U) / cm
        tot += t
        if m > 3 and abs(t) < mp.mpf(10) ** (-(mp.mp.dps + 5)):
            break
        m += 1
        if m > 10000:
            break

    def integrand(u):
        if u <= 0:
            return mp.mpf(0)
        R = mp.e ** (-u / 2) / (1 - mp.e ** (-2 * u))
        return (csym(fp, gp, u) - C0) * R

    pts = [mp.mpf(0)] + breakpoints(fp, gp, U) + [U]
    # Panels where the integrand is numerically zero must be SKIPPED, not integrated:
    # mp.quad tests a RELATIVE tolerance, which an identically-zero panel can never meet,
    # so it grinds to maxdegree (thousands of nodes) for no information.  Orthogonal
    # annulus modes produce many such panels.
    floor = mp.mpf(10) ** (-(mp.mp.dps + 8))
    val = mp.mpf(0)
    for lo, hi in zip(pts[:-1], pts[1:]):
        if hi <= lo:
            continue
        probe = max(abs(integrand(lo + (hi - lo) * mp.mpf(k) / 6)) for k in range(1, 6))
        if probe < floor:
            continue
        val += mp.quad(integrand, [lo, hi], method='gauss-legendre', maxdegree=maxdeg)
    return tot - 2 * val


def Q_entry(fp, gp, U, prim, maxdeg=8):
    """One entry of the Weil form.  prim = [(n, Lambda(n), log n), ...] with log n < U."""
    C0 = csym(fp, gp, mp.mpf(0))
    pole = 2 * p_half(fp) * p_half(gp)           # even sector: sgn = +1
    pr = mp.mpf(0)
    for n, lam, logn in prim:
        pr += 2 * lam / mp.sqrt(mp.mpf(n)) * csym(fp, gp, logn)
    return pole + arch(fp, gp, U, maxdeg) - mp.log(mp.pi) * C0 - pr


def Q_matrix(basis, U, prim, maxdeg=8):
    n = len(basis)
    M = mp.zeros(n, n)
    for i in range(n):
        for j in range(i, n):
            M[i, j] = M[j, i] = Q_entry(basis[i], basis[j], U, prim, maxdeg)
    return M


# ---------------------------------------------------------------- bases ----
def inner_basis(L, N, kind='even'):
    """EVEN: cos(k pi x / L) on [-L, L].  EVEN_D: cos((2k+1) pi x / 2L)."""
    out = []
    for k in range(N):
        a = (2 * k + 1) * mp.pi / (2 * L) if kind == 'even_d' else k * mp.pi / L
        out.append([(mp.mpf(0), a, -L, L)])
    return out


def annulus_basis(L, Lp, M, kind='even'):
    """Even functions on L <= |x| <= L': cos(b(|x|-L)) on each side."""
    W = Lp - L
    out = []
    for k in range(M):
        b = (2 * k + 1) * mp.pi / (2 * W) if kind == 'even_d' else k * mp.pi / W
        out.append([(L, b, L, Lp), (-L, b, -Lp, -L)])   # cos(b(x-L)) and cos(b(x+L))
    return out


def normalise(basis):
    return [[(x0, a, p, q) for (x0, a, p, q) in f] for f in basis], \
           [mp.sqrt(l2(f, f)) for f in basis]


def scale_matrix(M, nrm):
    n = M.rows
    out = mp.zeros(n, n)
    for i in range(n):
        for j in range(n):
            out[i, j] = M[i, j] / (nrm[i] * nrm[j])
    return out


def fro(M, rows=None, cols=None):
    rows = range(M.rows) if rows is None else rows
    cols = range(M.cols) if cols is None else cols
    return mp.sqrt(sum(M[i, j] ** 2 for i in rows for j in cols))


def sine_annulus_basis(L, Lp, M):
    """sin(k pi (|x|-L)/W), k=1..M -- vanishes at |x| = L and |x| = L' (continuous on R).
    sin(b(x-L)) = cos(b(x - L - pi/2b));  on x<0 the even extension is
    -sin(b(x+L)) = cos(b(x + L + pi/2b))."""
    W = Lp - L
    out = []
    for k in range(1, M + 1):
        b = k * mp.pi / W
        sh = mp.pi / (2 * b)
        out.append([(L + sh, b, L, Lp), (-L - sh, b, -Lp, -L)])
    return out
