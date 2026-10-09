#!/usr/bin/env python3
r"""Certified detection witness for the finite-window Weil form.  Standard library only.

WHY NOT THE REPO EVALUATOR.  src/zeta_window.py needs mpmath, which is not installed on
this machine (the brief forbids installs).  Its archimedean block is evaluated through
digamma / Lerch Phi / Hurwitz zeta, which is exactly the route Limitation 4 of the
manuscript records as NOT interval-certifiable ("mpmath's interval context has no digamma
and no Lerch transcendent, so the geometric side could not be interval-evaluated").

WHAT IS DONE INSTEAD.  The archimedean term is re-expressed as a series of ELEMENTARY
integrals, with an explicit tail bound:

  psi(z) = -euler + sum_{m>=0} [ 1/(m+1) - 1/(m+z) ]                         (classical)
  Re psi(1/4 + i a/2) + euler = sum_{m>=0} [ 1/(m+1) - 2 s_m/(s_m^2 + a^2) ],  s_m = 2m+1/2
  Arch(f) = -euler C(0) + sum_{m>=0} [ C(0)/(m+1) - 2 int_0^{2L} e^{-s_m u} C(u) du ]

with C the (symmetrised) autocorrelation.  Each integral is done in CLOSED FORM, so no
quadrature error enters, and the tail of the m-sum is bounded rigorously (see tail_bound).
This needs no digamma and no Lerch transcendent, i.e. it is a route around Limitation 4.

CONVENTIONS, verified against the source before use (see TARGET_AND_SOURCES.md):
  F(z) = int_{-L}^{L} f(u) e^{izu} du,  f real even, supp f in [-L,L], ||f||_2 = 1
  Q_zeta = Pole + Arch - LogPi - Prime          (src/zeta_window.py build_matrix)
  Q_0  = Q_zeta + 2 F(gamma)^2                  (paper eq. (3); proofs.md B.4)
  Q_delta - Q_0 = -4 delta^2 [F'^2 + F F''] + R_4          (paper eq. (4))
  |R_4| <= (16/3) L^5 delta^4 e^{2 L delta}               (proofs.md B.3)
  basis even_d: phi_k(x) = cos(w_k x)/sqrt(L),  w_k = (2k+1)pi/(2L)   (orthonormal)
"""
import cmath, json, math, sys

EULER = 0.57721566490153286060651209008240243104215933593992
LOG_PI = math.log(math.pi)
GAMMA1_STR = '14.134725141734693790457251983562'   # NOT certified here; see RESULT.md


# ----------------------------------------------------------------- basis and C-terms
def w_of(k, L):
    return (2*k + 1)*math.pi/(2*L)


def C_terms(j, k, L):
    """Exactly src/zeta_window.py::C_terms for sector even_d (s2 = +1).

    Returns a list of (coef, alpha, beta):  coef*sin(alpha*u + beta), or
    (coef, 'lin', beta): coef*(2L - u)*cos(beta*u).  NOT yet divided by nrm_j*nrm_k.
    """
    wj, wk = w_of(j, L), w_of(k, L)
    p1, p2 = wj - wk, wj + wk
    out = []
    if p1 == 0:
        out.append((0.5, 'lin', wk))
    else:
        out.append((1/(2*p1), wk, p1*L))
        out.append((-1/(2*p1), wj, -p1*L))
    if p2 == 0:
        out.append((0.5, 'lin', wk))
    else:
        out.append((1/(2*p2), -wk, p2*L))
        out.append((-1/(2*p2), wj, -p2*L))
    return out


def Csym_terms(j, k, L):
    """Symmetrised, and divided by nrm_j*nrm_k = L (even_d: norm2 = L)."""
    if j == k:
        raw = C_terms(j, k, L)
    else:
        raw = ([(c/2, a, b) for c, a, b in C_terms(j, k, L)] +
               [(c/2, a, b) for c, a, b in C_terms(k, j, L)])
    return [(c/L, a, b) for c, a, b in raw]


def eval_terms(terms, L, u):
    u = abs(u)
    if u >= 2*L:
        return 0.0
    tot = 0.0
    for c, a, b in terms:
        tot += c*((2*L - u)*math.cos(b*u) if a == 'lin' else math.sin(a*u + b))
    return tot


def dC_sup(terms, L):
    """A bound on sup|C'(u)| on [0,2L], used only for the m-sum tail."""
    U = 2*L
    tot = 0.0
    for c, a, b in terms:
        tot += abs(c)*(1 + U*abs(b)) if a == 'lin' else abs(c)*abs(a)
    return tot


# -------------------------------------------- exact Laplace integrals of the C-terms
def E1(z, U):
    """int_0^U e^{-z u} du."""
    return (1 - cmath.exp(-z*U))/z


def E2(z, U):
    """int_0^U u e^{-z u} du."""
    return (1 - cmath.exp(-z*U)*(1 + z*U))/(z*z)


def laplace_terms(terms, L, s):
    """int_0^{2L} e^{-s u} C(u) du, EXACTLY (no quadrature)."""
    U = 2*L
    tot = 0.0
    for c, a, b in terms:
        if a == 'lin':
            z = complex(s, -b)
            tot += c*(U*E1(z, U) - E2(z, U)).real
        else:
            z = complex(s, -a)
            tot += c*(cmath.exp(1j*b)*E1(z, U)).imag
    return tot


def tail_bound(C0, dCsup, M):
    """Rigorous bound on |sum_{m>=M} t_m|,  t_m = C0/(m+1) - 2 int e^{-s_m u} C(u) du.

    t_m = C0[1/(m+1) - 1/(m+1/4)] - (2I_m - 2C0/s_m), and integrating by parts
    |2I_m - 2C0/s_m| <= 2 sup|C'| / s_m^2  (C(2L) = 0, so no boundary term).
    Then sum_{m>=M} 1/((m+1)(m+1/4)) <= int_{M-1}^inf dx/(x+1/4)^2 = 1/(M-3/4) and
    sum_{m>=M} 1/(2m+1/2)^2 <= (1/4)/(M-3/4).
    """
    return (0.75*abs(C0) + 0.5*dCsup)/(M - 0.75)


def von_mangoldt(twoL):
    out = []
    n = 2
    while math.log(n) < twoL:
        m, p, lam = n, 2, None
        while p*p <= m:
            if m % p == 0:
                q = m
                while q % p == 0:
                    q //= p
                lam = math.log(p) if q == 1 else None
                break
            p += 1
        else:
            lam = math.log(m)
        if lam is not None:
            out.append((n, lam))
        n += 1
    return out


def Qjk(j, k, L, M, prim):
    """One entry of the Weil matrix in the orthonormal even_d basis, with its tail bound."""
    terms = Csym_terms(j, k, L)
    C0 = eval_terms(terms, L, 0.0)
    wj, wk = w_of(j, L), w_of(k, L)
    Pj = (-1)**j*2*wj*math.cosh(L/2)/(wj*wj + 0.25)/math.sqrt(L)
    Pk = (-1)**k*2*wk*math.cosh(L/2)/(wk*wk + 0.25)/math.sqrt(L)
    pole = 2*Pj*Pk
    prime = sum(2*lam/math.sqrt(n)*eval_terms(terms, L, math.log(n)) for n, lam in prim)
    arch = -EULER*C0
    for m in range(M):
        s = 2*m + 0.5
        arch += C0/(m + 1) - 2*laplace_terms(terms, L, s)
    err = tail_bound(C0, dC_sup(terms, L), M)
    return pole + arch - LOG_PI*C0 - prime, err


def weil_matrix(L, N, M=20000):
    prim = von_mangoldt(2*L)
    A = [[0.0]*N for _ in range(N)]
    E = [[0.0]*N for _ in range(N)]
    for j in range(N):
        for k in range(j, N):
            v, e = Qjk(j, k, L, M, prim)
            A[j][k] = A[k][j] = v
            E[j][k] = E[k][j] = e
    return A, E, prim


# ------------------------------------------------------------------- linear algebra
def jacobi(A, sweeps=300):
    n = len(A)
    Mx = [row[:] for row in A]
    V = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _ in range(sweeps):
        off = math.sqrt(sum(Mx[i][j]**2 for i in range(n) for j in range(n) if i != j))
        if off < 1e-20*max(1.0, max(abs(Mx[i][i]) for i in range(n))):
            break
        for p in range(n-1):
            for q in range(p+1, n):
                if abs(Mx[p][q]) < 1e-300:
                    continue
                th = (Mx[q][q] - Mx[p][p])/(2*Mx[p][q])
                t = (1 if th >= 0 else -1)/(abs(th) + math.sqrt(th*th + 1))
                c = 1/math.sqrt(t*t + 1)
                s = t*c
                for r in range(n):
                    a1, a2 = Mx[r][p], Mx[r][q]
                    Mx[r][p], Mx[r][q] = c*a1 - s*a2, s*a1 + c*a2
                for r in range(n):
                    a1, a2 = Mx[p][r], Mx[q][r]
                    Mx[p][r], Mx[q][r] = c*a1 - s*a2, s*a1 + c*a2
                for r in range(n):
                    a1, a2 = V[r][p], V[r][q]
                    V[r][p], V[r][q] = c*a1 - s*a2, s*a1 + c*a2
    pairs = sorted(range(n), key=lambda i: Mx[i][i])
    return ([Mx[i][i] for i in pairs], [[V[r][i] for r in range(n)] for i in pairs])


def quad_form(A, x):
    n = len(x)
    return sum(x[i]*sum(A[i][j]*x[j] for j in range(n)) for i in range(n))


def quad_form_err(E, x):
    n = len(x)
    return sum(abs(x[i])*sum(abs(E[i][j])*abs(x[j]) for j in range(n)) for i in range(n))


# --------------------------------------------------- the derivative-evaluation data
def sincpi(x):
    if abs(x) < 1e-13:
        return 1.0 - (math.pi*x)**2/6
    return math.sin(math.pi*x)/(math.pi*x)


def Fk(k, L, t):
    """F_k(t) = int phi_k(x) e^{itx} dx for the ORTHONORMAL even_d mode (real)."""
    w = w_of(k, L)
    return L*(sincpi((w - t)*L/math.pi) + sincpi((w + t)*L/math.pi))/math.sqrt(L)


def dFk(k, L, t, h=1e-5):
    """F_k'(t) in closed form: d/dt int phi cos(tx) dx = -int x phi sin(tx) dx."""
    w = w_of(k, L)
    # int_{-L}^{L} x cos(w x) sin(t x) dx / sqrt(L)
    def I(a):           # int_{-L}^{L} x sin(a x) dx  = 2[ -L cos(aL)/a + sin(aL)/a^2 ]
        if abs(a) < 1e-14:
            return 0.0
        return 2*(-L*math.cos(a*L)/a + math.sin(a*L)/a**2)
    # x cos(wx) sin(tx) = (x/2)[ sin((t+w)x) + sin((t-w)x) ]
    return -0.5*(I(t + w) + I(t - w))/math.sqrt(L)


def ddFk(k, L, t):
    """F_k''(t) = -int x^2 phi_k cos(tx) dx."""
    w = w_of(k, L)
    def J(a):           # int_{-L}^{L} x^2 cos(a x) dx
        if abs(a) < 1e-14:
            return 2*L**3/3
        return 2*(L*L*math.sin(a*L)/a + 2*L*math.cos(a*L)/a**2 - 2*math.sin(a*L)/a**3)
    # x^2 cos(wx) cos(tx) = (x^2/2)[cos((t+w)x) + cos((t-w)x)]
    return -0.5*(J(t + w) + J(t - w))/math.sqrt(L)


def K_exact(L, g):
    """K(L,gamma) = int_{-L}^{L} u^2 sin^2(gamma u) du.

    Transcribed from the TeX source (paper/weil_window_control.tex lines 372-375), NOT
    from the PDF text layer: pdftotext drops the outer bracket and so appears to flip the
    signs of the last two terms.  The TeX reads
        K = L^3/3 - [ L^2 sin(2gL)/(2g) + L cos(2gL)/(2g^2) - sin(2gL)/(4g^3) ]
    which is what is implemented, and which reproduces the direct quadrature.
    """
    return (L**3/3 - (L*L*math.sin(2*g*L)/(2*g) + L*math.cos(2*g*L)/(2*g*g)
                      - math.sin(2*g*L)/(4*g**3)))


# ============================================================ validation and witnesses
RECORDED_LAMMIN = {4: 1.7339346963e-10, 6: 5.96682932489e-13, 8: 4.96800734174e-15}
RECORDED_LAM2 = {4: 2.711086672e-6, 6: 2.32904375e-8, 8: 5.443266105e-10}


def validate(L, Ns, M):
    print("=" * 78)
    print(f"VALIDATION of the geometric side against data/stage3_gateB_lambda_min.csv")
    print(f"  L={L}  sector=even_d  archimedean m-sum truncated at M={M} with tail bound")
    print("=" * 78)
    out = []
    for N in Ns:
        A, E, prim = weil_matrix(L, N, M)
        ev, _ = jacobi(A)
        rec = RECORDED_LAMMIN[N]
        rec2 = RECORDED_LAM2[N]
        maxerr = max(max(r) for r in E)
        print(f"  N={N:2d}  lam_min(mine) = {ev[0]:+.6e}   recorded = {rec:.6e}   "
              f"rel dev = {abs(ev[0]-rec)/rec:.2e}")
        print(f"        lam_2  (mine) = {ev[1]:+.6e}   recorded = {rec2:.6e}   "
              f"rel dev = {abs(ev[1]-rec2)/rec2:.2e}")
        print(f"        max per-entry tail bound = {maxerr:.2e}   "
              f"prime powers used: {[n for n, _ in prim]}")
        out.append(dict(N=N, lam_min=ev[0], recorded=rec, rel_dev=abs(ev[0]-rec)/rec,
                        lam_2=ev[1], recorded_lam2=rec2, max_entry_tail=maxerr))
    return out


def witness_budget(L, g, N, M, gamma_label):
    """Witnesses A (f = r/||r||, projected to E) and B (nodal) on E = span(phi_0..phi_{N-1})."""
    A, E, prim = weil_matrix(L, N, M)
    K = K_exact(L, g)
    # coefficient vectors of P_E c and P_E r in the orthonormal basis
    av = [Fk(k, L, g) for k in range(N)]            # <phi_k, c> = F_k(gamma)
    bv = [-dFk(k, L, g) for k in range(N)]          # <phi_k, r> = -F_k'(gamma)
    na2 = sum(x*x for x in av)
    nb2 = sum(x*x for x in bv)
    ab = sum(av[i]*bv[i] for i in range(N))
    res = {}
    for name, vec in (('A_projected_r', bv[:]),
                      ('B_nodal', [bv[i] - (ab/na2)*av[i] for i in range(N)]
                       if na2 > 0 else bv[:])):
        n2 = sum(x*x for x in vec)
        if n2 == 0:
            res[name] = dict(degenerate=True)
            continue
        nrm = math.sqrt(n2)
        x = [t/nrm for t in vec]
        Fg = sum(x[k]*av[k] for k in range(N))
        Fp = -sum(x[k]*bv[k] for k in range(N))     # F'(gamma) = sum x_k F_k'
        Fpp = sum(x[k]*ddFk(k, L, g) for k in range(N))
        a_f = 4*(Fp*Fp + Fg*Fpp)
        Qz = quad_form(A, x)
        Qz_err = quad_form_err(E, x)
        b0 = Qz + 2*Fg*Fg
        res[name] = dict(
            L=L, gamma=g, gamma_label=gamma_label, basis='even_d', dim=N,
            selection_rule=('P_E r, no constraint' if name.startswith('A')
                            else 'P_E r with the F(gamma)=0 constraint projected out'),
            coeffs=x, norm_enclosure_halfwidth=0.0,
            F_gamma=Fg, Fp_gamma=Fp, Fpp_gamma=Fpp,
            a_f=a_f, K=K, nb2_PEr=nb2, K_eff=(nb2 - ab*ab/na2) if na2 > 0 else nb2,
            Q_zeta=Qz, Q_zeta_tail_bound=Qz_err, b_upper=b0 + Qz_err,
            b_0=b0)
    res['_shared'] = dict(K=K, norm_PEc2=na2, norm_PEr2=nb2, inner_PEc_PEr=ab,
                          K_eff=nb2 - ab*ab/na2 if na2 > 0 else nb2,
                          mu=ab/na2 if na2 > 0 else 0.0)
    return res


def detection_scan(b_upper, a_lower, L, deltas):
    """Q_delta <= b_upper - a_lower delta^2 + R_upper,  R_upper = (16/3)L^5 d^4 e^{2Ld}."""
    rows = []
    for d in deltas:
        R = (16/3)*L**5*d**4*math.exp(2*L*d)
        rhs = b_upper - a_lower*d*d + R
        rows.append(dict(delta=d, signal=a_lower*d*d, R_upper=R, rhs=rhs,
                         negative=bool(rhs < 0)))
    return rows


# --------------------------------------------- cached archimedean sums (O(N) not O(N^2))
def arch_cache(L, M, alphas, betas):
    """S[a]   = sum_{m<M} E1(s_m - i a, U)              (for sin(a u + b) terms)
       T[b]   = sum_{m<M} (U E1(s_m - i b) - E2(s_m - i b))   (for (U-u)cos(b u) terms)
       H      = sum_{m<M} 1/(m+1)
    Same caching idea the repo uses for its special-function values."""
    U = 2*L
    S = {a: 0j for a in alphas}
    T = {b: 0j for b in betas}
    H = 0.0
    for m in range(M):
        s = 2*m + 0.5
        H += 1.0/(m + 1)
        for a in alphas:
            S[a] += E1(complex(s, -a), U)
        for b in betas:
            z = complex(s, -b)
            T[b] += U*E1(z, U) - E2(z, U)
    return S, T, H


def Qjk_cached(j, k, L, prim, S, T, H, M):
    terms = Csym_terms(j, k, L)
    C0 = eval_terms(terms, L, 0.0)
    wj, wk = w_of(j, L), w_of(k, L)
    Pj = (-1)**j*2*wj*math.cosh(L/2)/(wj*wj + 0.25)/math.sqrt(L)
    Pk = (-1)**k*2*wk*math.cosh(L/2)/(wk*wk + 0.25)/math.sqrt(L)
    pole = 2*Pj*Pk
    prime = sum(2*lam/math.sqrt(n)*eval_terms(terms, L, math.log(n)) for n, lam in prim)
    acc = 0.0
    for c, a, b in terms:
        if a == 'lin':
            acc += c*T[b].real
        else:
            acc += c*(cmath.exp(1j*b)*S[a]).imag
    arch = -EULER*C0 + C0*H - 2*acc
    err = tail_bound(C0, dC_sup(terms, L), M)
    return pole + arch - LOG_PI*C0 - prime, err


def weil_matrix_fast(L, N, M):
    prim = von_mangoldt(2*L)
    alphas, betas = set(), set()
    for j in range(N):
        for k in range(j, N):
            for c, a, b in Csym_terms(j, k, L):
                (betas if a == 'lin' else alphas).add(b if a == 'lin' else a)
    S, T, H = arch_cache(L, M, sorted(alphas), sorted(betas))
    A = [[0.0]*N for _ in range(N)]
    E = [[0.0]*N for _ in range(N)]
    for j in range(N):
        for k in range(j, N):
            v, e = Qjk_cached(j, k, L, prim, S, T, H, M)
            A[j][k] = A[k][j] = v
            E[j][k] = E[k][j] = e
    return A, E, prim


# ----------------------------------------- exact finite-displacement quartet response
def Fk_c(k, L, t):
    """F_k(t) for COMPLEX t, same closed form (sinc), orthonormal even_d mode."""
    w = w_of(k, L)
    def sp(x):
        x = complex(x)
        if abs(x) < 1e-13:
            return 1 - (cmath.pi*x)**2/6
        return cmath.sin(cmath.pi*x)/(cmath.pi*x)
    return L*(sp((w - t)*L/math.pi) + sp((w + t)*L/math.pi))/math.sqrt(L)


def exact_response(x, L, g, d):
    """EXACT  Q_delta(x) - Q_zeta(x)  with no Taylor truncation.

    Q_delta = Q_0 + 4(<x,a_d>^2 - <x,b_d>^2 - <x,c>^2),  Q_0 = Q_zeta + 2<x,c>^2, and
    F_k(g + i d) = <phi_k,a_d> - i <phi_k,b_d>  (verified in the previous extension).
    Hence Q_delta(x) = Q_zeta(x) - 2<x,c>^2 + 4<x,a_d>^2 - 4<x,b_d>^2.
    The +4<x,a_d>^2 term is POSITIVE and is not dropped.
    """
    Z = sum(x[k]*Fk_c(k, L, complex(g, d)) for k in range(len(x)))
    A_, B_ = Z.real, -Z.imag
    Cc = sum(x[k]*Fk(k, L, g) for k in range(len(x)))
    return dict(a_inner=A_, b_inner=B_, c_inner=Cc,
                delta_Q=-2*Cc*Cc + 4*A_*A_ - 4*B_*B_,
                neg_part=4*B_*B_, pos_part=4*A_*A_ - 2*Cc*Cc)
