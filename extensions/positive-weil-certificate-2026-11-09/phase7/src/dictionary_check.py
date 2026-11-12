"""Independent implementation of Andrews-Connes-La Gatta eqs (6),(7),(8)
   [arXiv:2611.13129v1] applied to the repository's Legendre test functions,
   compared against phase5/certs/conventions.json.

   Dictionary under test:  lambda = e^L,  g(u) = f(log u),  QW_lambda(g,g) = Q(f).
"""
from mpmath import mp, mpf, sinh, cosh, exp, log, euler, pi, quad, tanh, legendre, sqrt

mp.dps = 50
L = mpf(4)/5

def autocorr_factory(coeffs):
    """f = sum coeffs[n] * phi_n,  phi_n(x) = sqrt((2n+1)/(2L)) P_n(x/L) on [-L,L]."""
    def f(x):
        if abs(x) > L: return mpf(0)
        return sum(c*sqrt((2*n+1)/(2*L))*legendre(n, x/L) for n, c in coeffs)
    def C(s):
        s = abs(s)
        if s >= 2*L: return mpf(0)
        # int f(t) f(t+s) dt over t in [-L, L-s]
        return quad(lambda t: f(t)*f(t+s), [-L, L-s])
    return f, C

def paper_psi(f, C):
    C0 = C(0)
    # W_{0,2}(F) = 2 * (int R e^{u/2}) * (int R e^{-u/2})
    Ip = quad(lambda u: f(u)*exp(u/2), [-L, L])
    Im = quad(lambda u: f(u)*exp(-u/2), [-L, L])
    W02 = 2*Ip*Im
    # W_R, eq (8), in t = log x
    def integrand(s):
        return (C(s)*exp(s/2) - C0)/sinh(s)
    inner = quad(integrand, [0, 2*L])
    outer = C0*log(tanh(L))          # int_{2L}^inf (-C0)/sinh s ds
    WR = (log(4*pi) + euler)*C0 + inner + outer
    # sum_p W_p, eq (7)
    Wp = mpf(0)
    for p in (2, 3, 5, 7, 11, 13):
        m = 1
        while m*log(p) < 2*L:
            Wp += log(p)*p**(mpf(-m)/2)*2*C(m*log(p))
            m += 1
    return W02 - WR - Wp, dict(W02=W02, WR=WR, Wp=Wp, C0=C0, Ip=Ip, Im=Im)

cases = {
    "phi_0":          [(0, mpf(1))],
    "phi_2":          [(2, mpf(1))],
    "phi_0 - phi_2":  [(0, mpf(1)), (2, mpf(-1))],
}
repo = {   # phase5/certs/conventions.json, ball midpoints
    "phi_0":         "0.0802740549592696112136916593813",
    "phi_2":         "0.389336213205013655883159773666",
    "phi_0 - phi_2": "0.117009437226589934779775252494",
}
repo_pole = {
    "phi_0":         "1.29891290374844724999391529834",
    "phi_2":         "0.0305162242253810235887452910789",
    "phi_0 - phi_2": "1.26839667952306622640517000726",
}

print(f"L = 4/5,  lambda = e^L = {exp(L)},  lambda^2 = {exp(2*L)}")
for name, co in cases.items():
    f, C = autocorr_factory(co)
    val, parts = paper_psi(f, C)
    ref = mpf(repo[name]); pref = mpf(repo_pole[name])
    print(f"\n--- {name} ---")
    print(f"  ||f||^2 (=C(0))            {parts['C0']}")
    print(f"  int f e^{{-u/2}}  (pole mom) {parts['Im']}")
    print(f"     repo pole_moment        {pref}   rel.diff {abs(parts['Im']-pref)/abs(pref)}")
    print(f"  W_02 = {parts['W02']}   (2*poleMom^2 = {2*parts['Im']**2})")
    print(f"  W_R  = {parts['WR']}")
    print(f"  sumWp= {parts['Wp']}")
    print(f"  PAPER  Psi(F) = {val}")
    print(f"  REPO   Q(f)   = {ref}")
    print(f"  abs diff = {abs(val-ref)}   rel diff = {abs(val-ref)/abs(ref)}")
