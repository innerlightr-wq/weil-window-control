"""Finite-block assembly and certification at L = 4/5, T = 200.

EVERYTHING here is an enclosure.  Any NaN, infinite radius or undecided sign
raises -- nothing is ever skipped (the phase-2 M_sup incident).

------------------------------------------------------------------- OBJECTS

Basis: orthonormal Legendre on [-L,L],  phi_n(u) = sqrt((2n+1)/(2L)) P_n(u/L),
degrees n = 0 .. N-1.  N = 220 therefore means degrees 0..219, i.e. 110 modes
in EACH parity (even n, odd n).

Fourier convention: F(t) = int_{-L}^{L} f(u) e^{i t u} du, int_R|F|^2 = 2pi||f||^2.
Then F_n(t) = sqrt(2L(2n+1)) i^n j_n(Lt).  For real f = sum x_n phi_n,
conj F(t) = F(-t), and with the NORM-PRESERVING diagonal sign change
y_n = (-1)^{floor(n/2)} x_n ... more simply, writing

      v_n(t) := sqrt(2L(2n+1)) j_n(Lt)        (REAL, parity of n),

one has  |F(t)|^2 = ( sum_n y_n v_n(t) )^2  with |y| = |x|, so the quadratic
form is represented by the REAL SYMMETRIC matrix below.  Even and odd n
decouple exactly, because v_n v_m is odd in t when n+m is odd.

The comparison operator (phase 2, DERIVATION sec. 5) is

      R_T(f) = beta ||f||^2 + 2 s F(i/2)^2 + (1/2pi) int_{-T}^{T}(Psi_L-beta)|F|^2

with s = +1 (even sector) and s = -1 (odd sector), so in each sector

      A = beta I + 2 s p p^T + C ,
      C_{nm} = (1/pi) int_0^{T} W(t) v_n(t) v_m(t) dt ,    W = Psi_L - beta ,
      p_n    = sqrt(2L(2n+1)) i_n(L/2)   (up to the same diagonal signs).

WHAT A CERTIFIES.  A is the EXACT compression of R_T to span{phi_n : n < N} --
not of Q, and not of the full infinite comparison.  a := a certified lower bound
for lam_min(A).  The whole-window implication needs a, d and kappa together
(see certify_final).
"""
import os, sys, json, time, math
from flint import arb, acb, arb_mat, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import baseline as B


# --------------------------------------------------------------- strictness
def need_finite(x, what):
    if not x.is_finite():
        raise RuntimeError("NON-FINITE enclosure in %s: %s" % (what, x))
    return x

def need_sign(x, what):
    """Return +1/-1; raise if the sign of the ball is undecided."""
    if x > 0:
        return 1
    if x < 0:
        return -1
    raise RuntimeError("UNDECIDED SIGN in %s: %s" % (what, x.str(10)))


# ------------------------------------------------------- certified GL rule
def gl_rule(k):
    """Certified Gauss-Legendre nodes/weights on (-1,1) from FLINT."""
    nodes, wts = [], []
    for i in range(k):
        x, w = arb.legendre_p_root(k, i, weight=True)
        need_finite(x, "GL node"); need_finite(w, "GL weight")
        nodes.append(x); wts.append(w)
    s = sum(wts, arb(0))
    if not ((s - 2).abs_upper() < arb(10) ** -30):
        raise RuntimeError("GL weights do not sum to 2: %s" % s.str(20))
    return nodes, wts


# -------------------------------------------------------------- the symbol
class Symbol:
    def __init__(self, L_num, L_den, T):
        self.Ln, self.Ld, self.T = L_num, L_den, T
        self.L = arb(L_num) / arb(L_den)
        self.terms = B.comb_terms(L_num, L_den)
        self.A_L = B.A_L(L_num, L_den)
        self.beta = B.H(T) - self.A_L

    def W_real(self, t):
        """Psi_L(t) - beta on the real axis."""
        c = arb(0)
        for n, p in self.terms:
            c += 2 * arb(p).log() / arb(n).sqrt() * (t * arb(n).log()).cos()
        H = acb(arb(1) / 4, t / 2).digamma().real - arb.pi().log()
        return need_finite(H - c - self.beta, "W_real")

    @staticmethod
    def digamma_shift(z, K=24, Kmax=400):
        """psi(z) = psi(z+K) - sum_{k=0}^{K-1} 1/(z+k), evaluated so that the
        digamma call sees an argument far from its poles.

        Arb's acb_digamma returns NaN on balls whose real part has radius
        >~ 0.1 even when the ball lies strictly in the right half plane, because
        its internal argument reduction cannot decide the shift.  Doing the
        recurrence explicitly fixes that: z+K has real part >= K - 1/4, where a
        wide ball is harmless, and the finite sum has no pole because
        Re z > 0 on every box we use.
        """
        if not (z.real > 0):
            raise RuntimeError("digamma_shift needs Re z > 0; got %s" % z.real.str(8))
        while K <= Kmax:
            tail = (z + acb(K)).digamma()
            if tail.is_finite():
                acc = acb(0)
                for k in range(K):
                    acc += 1 / (z + acb(k))
                out = tail - acc
                if out.is_finite():
                    return out
            K *= 2
        raise RuntimeError("NON-FINITE digamma_shift at %s even with K=%d" % (z.str(8), Kmax))

    def W_hol(self, z):
        """HOLOMORPHIC extension: [psi(1/4+iz/2)+psi(1/4-iz/2)]/2 - log pi - P_L - beta.

        Re psi(1/4+iz/2) is NOT holomorphic; this symmetric combination is, and
        agrees with W_real on the real axis.  Nearest poles at z = +- i/2.
        """
        a = acb(arb(1) / 4, 0) + acb(0, 1) * z / 2
        b = acb(arb(1) / 4, 0) - acb(0, 1) * z / 2
        H = (self.digamma_shift(a) + self.digamma_shift(b)) / 2 - acb(arb.pi().log())
        c = acb(0)
        for n, p in self.terms:
            c += acb(2 * arb(p).log() / arb(n).sqrt()) * (z * acb(arb(n).log())).cos()
        out = H - c - acb(self.beta)
        if not out.is_finite():
            raise RuntimeError("NON-FINITE W_hol at %s" % z.str(8))
        return out


def v_real(L, t, N):
    """v_n(t) = sqrt(2L(2n+1)) j_n(Lt), n = 0..N-1, as Arb balls.  t must be > 0."""
    x = L * t
    if not (x > 0):
        raise RuntimeError("v_real requires t > 0 (t=0 is a removable point and is "
                           "never a Gauss-Legendre node of an interior panel)")
    pref = (arb.pi() / (2 * x)).sqrt()
    out = []
    for n in range(N):
        j = pref * acb(x).bessel_j(acb(n) + acb(1) / 2).real
        need_finite(j, "j_%d(%s)" % (n, x.str(6)))
        out.append((2 * L * (2 * n + 1)).sqrt() * j)
    return out

def v_hol_absmax_uniform(L, b_e):
    """UNIFORM bound  |v_n(z)| <= sqrt(2L) exp(L b_e)  for every n and every z
    with |Im z| <= b_e.

    From the Legendre integral representation
        j_n(w) = (1/2)(-i)^n int_{-1}^{1} e^{i w t} P_n(t) dt ,
    valid for all complex w, we get
        |j_n(w)| <= (1/2) e^{|Im w|} int_{-1}^{1}|P_n|
                 <= (1/2) e^{|Im w|} sqrt(2) ||P_n||_2
                  = e^{|Im w|} / sqrt(2n+1) ,
    using ||P_n||_2^2 = 2/(2n+1).  Hence
        |v_n(z)| = sqrt(2L(2n+1)) |j_n(Lz)| <= sqrt(2L) e^{L |Im z|} .

    This REPLACES per-panel ball evaluation of j_n.  Arb's acb_bessel_j on a ball
    of radius ~0.26 at |w| ~ 80 returns an enclosure of width ~1e32 -- useless,
    and not fixable by precision (verified at 900 and 3000 bits).  The bound
    above is cheap, uniform in n, and valid on the whole strip.
    """
    return (2 * L).sqrt() * (L * b_e).exp()


def pole_sigma(n):
    """sigma_n = (-1)^{floor(n/2)}: the phase that converts x-coordinates
    (f = sum x_n phi_n) into the y-coordinates in which the quadratic form is
    sum y_n y_m v_n v_m with v_n = sqrt(2L(2n+1)) j_n(Lt).

    |F|^2 carries a factor (-1)^{(n-m)/2} = sigma_n sigma_m that is dropped when
    the form is written through the real v_n, so EVERY vector expressed in the
    same coordinates -- in particular the pole vector -- must carry sigma too.
    Omitting it is invisible at n = 0 (sigma_0 = +1) and wrong from n = 2 on.
    """
    return -1 if (n // 2) % 2 else 1


def pole_vector(L, N, with_sigma=True):
    """p_n = sqrt(2L(2n+1)) i_n(L/2), optionally in y-coordinates (sigma applied).

    with_sigma=False returns the x-coordinate vector, i.e. the one for which
    F(i/2) = sum_n x_n p_n; that form is used only by the validation test."""
    z = L / 2
    pref = (arb.pi() / (2 * z)).sqrt()
    out = []
    for n in range(N):
        iv = pref * acb(z).bessel_i(acb(n) + acb(1) / 2).real
        need_finite(iv, "i_%d" % n)
        v = (2 * L * (2 * n + 1)).sqrt() * iv
        out.append(pole_sigma(n) * v if with_sigma else v)
    return out


# ------------------------------------------------------ quadrature remainder
def panel_rho(h, b_e):
    """Bernstein parameter for the ellipse of imaginary semi-axis b_e on a panel
    of half-length h:  h (rho - 1/rho)/2 = b_e."""
    r = b_e / h
    return r + (r * r + 1).sqrt()

def panel_error_factor(h, rho, k):
    """Trefethen ATAP Thm 19.3, scaled to a panel of half-length h:
         |I - I_k| <= h * 64 M / (15 (rho^{2k} - 1)) ,
       with M an upper bound for |integrand| on the ellipse."""
    den = 15 * (rho ** (2 * k) - 1)
    if not (den > 0):
        raise RuntimeError("degenerate Bernstein denominator")
    return h * 64 / den


# ------------------------------------------------------------- the assembly
def assemble(L_num=4, L_den=5, T=200, N=220, panels=800, korder=32,
             b_e="0.30", prec_bessel=900, prec_mat=320, chunk=1600,
             progress=None, checkpoint=None):
    """Returns (A_even, A_odd, meta).  Entries are Arb balls ENCLOSING the exact
    matrix, with the quadrature remainder folded in."""
    t_start = time.time()
    ctx.prec = prec_bessel
    S = Symbol(L_num, L_den, T)
    L = S.L
    beta = S.beta
    nodes, wts = gl_rule(korder)
    hw = arb(T) / (2 * panels)                       # panel half-length
    b_e = arb(b_e)
    rho = panel_rho(hw, b_e)
    efac = panel_error_factor(hw, rho, korder)       # multiplies M_panel
    rad_re = hw * (rho + 1 / rho) / 2                # real semi-axis of the ellipse
    VMAX = v_hol_absmax_uniform(L, b_e)              # uniform |v_n| on the strip
    if not (b_e < arb(1) / 2):
        raise RuntimeError("ellipse must stay strictly inside the pole-free strip "
                           "|Im z| < 1/2; got b_e = %s" % b_e.str(10))

    ev = [n for n in range(N) if n % 2 == 0]
    od = [n for n in range(N) if n % 2 == 1]
    ne, no = len(ev), len(od)
    Ce = [[arb(0)] * ne for _ in range(ne)]
    Co = [[arb(0)] * no for _ in range(no)]
    Ee = [arb(0)] * ne                                # error row-vectors (rank-1 bound)
    Eo = [arb(0)] * no
    Wmax_all = arb(0)

    rowsE, rowsO, scalE, scalO = [], [], [], []

    def flush():
        nonlocal rowsE, rowsO, scalE, scalO
        if not rowsE:
            return
        ctx.prec = prec_mat
        V = arb_mat(rowsE)
        DV = arb_mat([[scalE[i] * rowsE[i][j] for j in range(ne)] for i in range(len(rowsE))])
        P = V.transpose() * DV
        for a in range(ne):
            for b in range(ne):
                Ce[a][b] += P[a, b]
        V = arb_mat(rowsO)
        DV = arb_mat([[scalO[i] * rowsO[i][j] for j in range(no)] for i in range(len(rowsO))])
        P = V.transpose() * DV
        for a in range(no):
            for b in range(no):
                Co[a][b] += P[a, b]
        rowsE, rowsO, scalE, scalO = [], [], [], []
        ctx.prec = prec_bessel

    for ip in range(panels):
        c = hw * (2 * ip + 1)                         # panel centre on [0,T]
        # ---- quadrature remainder for this panel
        # BOX containing the Bernstein ellipse: wide in Re, only b_e in Im, so it
        # stays strictly inside the pole-free strip |Im z| < 1/2.
        ball = acb(arb(c.mid(), rad_re.upper()), arb(0, b_e.upper()))
        Wb = S.W_hol(ball)
        Wm = arb(Wb.abs_upper())
        need_finite(Wm, "W_hol sup on panel %d" % ip)
        pe = efac * Wm
        # err_{nm} <= sum_panels pe * vmax^2 <= (sum_p sqrt(pe) vmax)^2, and vmax is
        # uniform in n, so one accumulator suffices per sector.
        for i in range(ne):
            Ee[i] += pe.sqrt() * VMAX
        for i in range(no):
            Eo[i] += pe.sqrt() * VMAX
        if Wm > Wmax_all:
            Wmax_all = Wm
        # ---- nodes
        for q in range(korder):
            t = c + hw * nodes[q]
            W = S.W_real(t)
            vv = v_real(L, t, N)
            rowsE.append([vv[n] for n in ev])
            rowsO.append([vv[n] for n in od])
            s = wts[q] * hw * W / arb.pi()
            scalE.append(s); scalO.append(s)
        if len(rowsE) >= chunk:
            flush()
        if progress and ip % progress == 0:
            print("    panel %d/%d  (%.0fs)" % (ip, panels, time.time() - t_start), flush=True)
    flush()
    ctx.prec = prec_mat

    # ---- fold the quadrature remainder in as entrywise widening
    # err_{nm} <= ( sum_panels pe * vm_n * vm_m ) <= Ee_n * Ee_m  by Cauchy-Schwarz
    # on the panel index; we use the cruder but valid product of the accumulated
    # sqrt(pe)*vm sums.
    p_all = pole_vector(L, N)
    pe_v = [p_all[n] for n in ev]
    po_v = [p_all[n] for n in od]

    def finish(C, E, pv, s, n_):
        M = [[arb(0)] * n_ for _ in range(n_)]
        for a in range(n_):
            for b in range(n_):
                e = (E[a] * E[b]) / arb.pi()
                val = C[a][b] + (2 * s) * pv[a] * pv[b]
                if a == b:
                    val = val + beta
                M[a][b] = arb(val.mid(), (arb(val.rad()) + e).upper())
                need_finite(M[a][b], "A[%d][%d]" % (a, b))
        return M

    Ae = finish(Ce, Ee, pe_v, +1, ne)
    Ao = finish(Co, Eo, po_v, -1, no)
    meta = dict(L="%d/%d" % (L_num, L_den), T=T, N=N, panels=panels, korder=korder,
                b_e=b_e.str(6), rho=rho.str(10), panel_halflength=hw.str(10),
                ellipse_re_semiaxis=rad_re.str(10),
                rho_pow_2k=(rho ** (2 * korder)).str(10),
                prec_bessel=prec_bessel, prec_mat=prec_mat,
                beta=beta.str(25), A_L=S.A_L.str(25),
                W_sup_on_ellipses=Wmax_all.str(10),
                v_uniform_bound=VMAX.str(12),
                quad_error_factor=efac.str(10),
                modes_even=ne, modes_odd=no,
                elapsed_s=round(time.time() - t_start, 1))
    return Ae, Ao, meta


# ------------------------------------------------------------ certification
def ldl_pivots(A, shift):
    """Interval LDL^T of (A - shift*I).  Returns the pivot balls.

    If every returned pivot is certifiably > 0 then A - shift*I is positive
    definite for EVERY matrix in the entrywise enclosure, hence lam_min(A) > shift.
    A pivot whose sign is undecided raises: a failed factorization is NOT a proof
    of a negative eigenvalue.
    """
    n = len(A)
    Lm = [[arb(0)] * n for _ in range(n)]
    D = [arb(0)] * n
    for j in range(n):
        acc = A[j][j] - shift
        for k in range(j):
            acc = acc - Lm[j][k] * Lm[j][k] * D[k]
        need_finite(acc, "pivot %d" % j)
        if need_sign(acc, "pivot %d" % j) < 0:
            raise RuntimeError("NEGATIVE PIVOT at %d: %s -- factorization failed, "
                               "which is NOT a proof of indefiniteness" % (j, acc.str(10)))
        D[j] = acc
        for i in range(j + 1, n):
            t = A[i][j]
            for k in range(j):
                t = t - Lm[i][k] * Lm[j][k] * D[k]
            Lm[i][j] = t / D[j]
            need_finite(Lm[i][j], "L[%d][%d]" % (i, j))
    return D


def certify_block(A, shift, prec):
    ctx.prec = prec
    try:
        D = ldl_pivots(A, arb(shift))
    except RuntimeError as e:
        return False, str(e), None
    return True, None, min(D, key=lambda z: float(z.lower()))


def certify_final(a, d, kappa, mu):
    """Certify lam_min(R_T) >= mu from A >= a I, D >= d I, ||B|| <= kappa.

    Uses the three STRICT interval sign tests  a > mu,  d > mu,
    (a - mu)(d - mu) > kappa^2 -- not the cancellation-prone 2x2 eigenvalue
    formula.  For x = x_N + x_perp,
        <(R - mu) x, x> >= (a-mu)|x_N|^2 - 2 kappa |x_N||x_perp| + (d-mu)|x_perp|^2,
    which is >= 0 exactly when both diagonal terms are >= 0 and the discriminant
    (a-mu)(d-mu) - kappa^2 is >= 0.
    """
    a, d, kappa, mu = arb(a), arb(d), arb(kappa), arb(mu)
    c1 = bool(a > mu)
    c2 = bool(d > mu)
    c3 = bool((a - mu) * (d - mu) > kappa * kappa)
    return dict(a_gt_mu=c1, d_gt_mu=c2, disc_gt_kappa2=c3, ok=(c1 and c2 and c3),
                margin=((a - mu) * (d - mu) - kappa * kappa).str(12))


# -------------------------------------------------- independent cross-check
def quad_entry(S, n, m, a, b, panels, korder):
    """C-entry contribution (1/pi) int_a^b W v_n v_m dt by the SAME panelled GL
    rule used in assemble(), for one entry only."""
    L = S.L
    nodes, wts = gl_rule(korder)
    hw = (arb(b) - arb(a)) / (2 * panels)
    tot = arb(0)
    cn = (2 * L * (2 * n + 1)).sqrt() * (2 * L * (2 * m + 1)).sqrt()
    for ip in range(panels):
        c = arb(a) + hw * (2 * ip + 1)
        for q in range(korder):
            t = c + hw * nodes[q]
            x = L * t
            pref = (arb.pi() / (2 * x)).sqrt()
            jn = pref * acb(x).bessel_j(acb(n) + acb(1) / 2).real
            jm = pref * acb(x).bessel_j(acb(m) + acb(1) / 2).real
            tot += wts[q] * hw * S.W_real(t) * cn * jn * jm
    return tot / arb.pi()


def rigorous_entry(S, n, m, a, b):
    """Same quantity by FLINT's adaptive rigorous integrator, on [a,b] with a > 0.

    The 1/sqrt(w) representation of j_n is singular at w = 0, and the adaptive
    integrator SWALLOWS the resulting exception and returns a wrong value rather
    than failing -- observed during this build.  So the cross-check is run only
    on intervals bounded away from the origin, and the returned value is
    additionally required to be finite.
    """
    if not (arb(a) > 0):
        raise RuntimeError("rigorous_entry needs a > 0 (1/sqrt(w) singularity)")
    L = S.L
    Lc, half, pi_c = acb(L), acb(1) / 2, acb(arb.pi())
    cn = acb((2 * L * (2 * n + 1)).sqrt() * (2 * L * (2 * m + 1)).sqrt())
    nn, mm = acb(n) + half, acb(m) + half

    NAN = acb(arb("nan"), arb("nan"))

    def f(z, analytic):
        # Protocol: return a NON-FINITE value when the box may meet a pole, so the
        # integrator subdivides instead of silently producing a wrong answer.
        # W_hol has poles at Im z = +-1/2; j_n's 1/sqrt(w) form is singular at 0.
        try:
            if not (z.imag.abs_upper() < arb(1) / 2):
                return NAN
            if not (arb((z * Lc).abs_lower()) > 0):
                return NAN
            w = z * Lc
            jj = (pi_c / (2 * w)).sqrt()
            val = S.W_hol(z) * (jj * w.bessel_j(nn)) * (jj * w.bessel_j(mm)) * cn
            return val if val.is_finite() else NAN
        except Exception:
            return NAN

    r = acb.integral(f, arb(a), arb(b))
    if not r.is_finite():
        raise RuntimeError("rigorous integrator returned a non-finite value")
    return r.real / arb.pi()


def validate_conventions(L_num=4, L_den=5, ts=(0.7, 3.5, 11.0), degs=(0, 2, 4), prec=300):
    """Check |F(t)|^2 = (sum_n y_n v_n(t))^2 and F(i/2) = sum_n x_n p_n^{(x)}
    against DIRECT numerical integration of f(u) = sum_n x_n phi_n(u).

    NUMERICAL consistency check: mpmath's quadrature is not certified, so the
    comparison is a relative tolerance (1e-25), not a ball overlap.

    This is the test that catches a wrong sigma phase: it is insensitive at
    degree 0 (sigma_0 = +1) and fails from degree 2 on.
    """
    import mpmath as mp
    ctx.prec = prec
    mp.mp.dps = 40
    L = arb(L_num) / arb(L_den)
    Lf = mp.mpf(L_num) / L_den

    def phi(n, u):
        return mp.sqrt((2 * n + 1) / (2 * Lf)) * mp.legendre(n, u / Lf)

    def f(u):
        return sum(phi(n, u) for n in degs)

    out = []
    for t in ts:
        Fr = mp.quad(lambda u: f(u) * mp.cos(mp.mpf(t) * u), [-Lf, 0, Lf])
        Fi = mp.quad(lambda u: f(u) * mp.sin(mp.mpf(t) * u), [-Lf, 0, Lf])
        direct = arb(mp.nstr(Fr ** 2 + Fi ** 2, 30))
        vv = v_real(L, arb(t), max(degs) + 1)
        model = sum((arb(pole_sigma(n)) * vv[n] for n in degs), arb(0)) ** 2
        rel = ((direct - model) / direct).abs_upper()
        out.append(("|F(%s)|^2" % t, direct, model, bool(arb(rel) < arb(10) ** -25)))
    Fh = mp.quad(lambda u: f(u) * mp.e ** (-u / 2), [-Lf, 0, Lf])
    px = pole_vector(L, max(degs) + 1, with_sigma=False)
    d = arb(mp.nstr(Fh, 30)); m = sum((px[n] for n in degs), arb(0))
    out.append(("F(i/2)", d, m, bool(arb(((d - m) / d).abs_upper()) < arb(10) ** -25)))
    return out


# -------------------------------------------- preconditioned PD verification
def chol_midpoint(B, prec):
    """Approximate upper-triangular Cholesky R of the MIDPOINT of B.

    Midpoints only -- R is a preconditioner, so it needs no enclosure; its
    quality only affects how small ||I - R^{-T} B R^{-1}|| comes out.
    """
    ctx.prec = prec
    n = len(B)
    M = [[arb(B[i][j].mid()) for j in range(n)] for i in range(n)]
    R = [[arb(0)] * n for _ in range(n)]
    for i in range(n):
        s = M[i][i] - sum((R[k][i] * R[k][i] for k in range(i)), arb(0))
        if not (s > 0):
            raise RuntimeError("midpoint Cholesky failed at %d (s=%s): the MIDPOINT "
                               "matrix is not positive definite" % (i, s.str(10)))
        R[i][i] = s.sqrt()
        for j in range(i + 1, n):
            R[i][j] = (M[i][j] - sum((R[k][i] * R[k][j] for k in range(i)), arb(0))) / R[i][i]
    return R


def verify_pd_preconditioned(A, shift, prec):
    """Certify A - shift*I > 0 by Rump's criterion ||I - R^{-T}(A-sI)R^{-1}|| < 1.

    Naive interval LDL^T is unusable here: the entry radii are ~1e-38 and the
    pivot recursion amplifies them by about 6.5 per step, reaching O(1) by step
    ~47 (measured).  Preconditioning replaces that growth by a factor of order
    cond(R) = sqrt(cond(A)), which is ~1e9 here.

    Returns (ok, ||I-G||_infinity enclosure, message).
    """
    n = len(A)
    ctx.prec = prec
    s = arb(shift)
    B = arb_mat([[A[i][j] - (s if i == j else arb(0)) for j in range(n)] for i in range(n)])
    try:
        R = chol_midpoint(A if shift == 0 else
                          [[A[i][j] - (s if i == j else arb(0)) for j in range(n)]
                           for i in range(n)], prec)
    except RuntimeError as e:
        return False, None, str(e)
    ctx.prec = prec
    # Explicit FORWARD SUBSTITUTION for the lower-triangular system R^T X = B.
    # arb_mat.solve does a general ball LU and returned residuals ~1e20 on this
    # well-conditioned triangular system; substitution is exact in structure and
    # amplifies radii only by cond(R).
    Bl = [[B[i, j] for j in range(n)] for i in range(n)]

    def fwd(RHS):
        """Solve R^T X = RHS, R upper triangular so R^T lower triangular."""
        X = [[arb(0)] * n for _ in range(n)]
        for col in range(n):
            for i in range(n):
                acc = RHS[i][col]
                for k in range(i):
                    acc = acc - R[k][i] * X[k][col]
                X[i][col] = acc / R[i][i]
                if not X[i][col].is_finite():
                    raise RuntimeError("non-finite in forward substitution")
        return X

    try:
        X = fwd(Bl)                                   # R^{-T} B
        Xt = [[X[j][i] for j in range(n)] for i in range(n)]
        Y = fwd(Xt)                                   # (R^{-T} B R^{-1})^T
    except RuntimeError as e:
        return False, None, str(e)
    G = arb_mat([[Y[j][i] for j in range(n)] for i in range(n)])
    worst = arb(0)
    for i in range(n):
        row = arb(0)
        for j in range(n):
            d = G[i, j] - (arb(1) if i == j else arb(0))
            if not d.is_finite():
                return False, None, "non-finite residual at (%d,%d)" % (i, j)
            row += arb(d.abs_upper())
        if row > worst:
            worst = row
    return bool(worst < 1), worst, None
