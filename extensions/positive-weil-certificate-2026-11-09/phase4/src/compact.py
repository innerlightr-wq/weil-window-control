"""Phase 5, part 2: evaluate the ORIGINAL form Q on a frozen vector, in the time
domain, with no long oscillatory frequency integral.

DISCLOSURE: this is a SECOND REPRESENTATION of the same form, not an independent
implementation.  It shares the Arb library, `block.Symbol` for the comb data and
`block.pole_vector`; what is new is the time-domain route (autocorrelation +
compact archimedean identity) replacing the frequency-domain assembly.

--------------------------------------------------------------- THE IDENTITY

With g = f * f~ the autocorrelation, g(u) = int f(t) f(t-u) dt, supported in
|u| <= 2L, and gamma_E Euler's constant (Zhu arXiv:2608.24827v2, Lemma 2.5, in
the repository's normalization):

  Arch(g) = -[gamma_E + log(pi) + log(1 - e^{-4L})] g(0)
            + int_0^{2L} 2[ e^{-2x} g(0) - e^{-x/2} g(x) ] / [1 - e^{-2x} ] dx

and the complete form is

  Q(f) = 2 s F(i/2)^2 + Arch(g) - sum_{n: log n < 2L} (2 Lambda(n)/sqrt n) g(log n) .

Note Arch already carries the -log(pi) g(0) term, so no separate LogPi subtraction.

REMOVABLE ENDPOINT.  Both numerator and denominator vanish at x = 0.  With
N(x) = e^{-2x}g(0) - e^{-x/2}g(x), D(x) = 1 - e^{-2x}: N(0) = 0, D(0) = 0,
D'(0) = 2, and N'(0) = -(3/2) g(0) - g'(0+), so

    lim_{x->0+} 2N/D = -g'(0+) - (3/2) g(0) = c - (3/2) g(0),   c := -g'(0+).

Moreover, for f EVEN, substituting s = t - x,
    g(x) = int_{-L}^{L-x} f(s+x) f(s) ds ,
    g'(x) = -f(L) f(L-x) + int_{-L}^{L-x} f'(s+x) f(s) ds ,
and at x = 0 the integral is (1/2)[f^2]_{-L}^{L} = 0 because f is even.  Hence

    c = -g'(0+) = f(L)^2 ,                                              (EXACT)

which is how the limit is certified here rather than estimated.

QUADRATURE.  On [0,2L] the autocorrelation g is a polynomial of degree at most
2 deg f + 1 -- NOT 2 deg f: for f = 1 one gets g(x) = 2L - x.  What governs the
INNER rule is the degree of the integrand f(s+x)f(s) in s, which is at most
2 deg f, so the inner integral is computed by Gauss-Legendre with deg f + 2
nodes, which is
EXACT -- no inner quadrature error at all.  The outer integrand is analytic in
|Im x| < pi (the poles of 1/(1-e^{-2x}) sit at x = i pi k, k != 0; the k = 0 pole
is cancelled), so panelled Gauss-Legendre is used, with GEOMETRIC panels near 0
so that no ellipse would enclose the origin.

*** WARNING -- INCOMPLETE ERROR BUDGET (RETRACTIONS.md R14) ***
This routine adds the explicit head bound on [0,eta] and NOTHING for the
TRUNCATION error of the outer panelled rule.  An earlier version of this
docstring claimed a Bernstein-ellipse remainder was used; none is computed.
Ball arithmetic encloses nodes, weights and rounding -- not the difference
between a Gauss-Legendre sum and the integral.  Every value this module
produces is therefore a MEASUREMENT, not a certificate, until that bound is
derived and added.  The certified upper bound now comes from the closed-form
route in src/certify.py (no quadrature anywhere); see
phase5/src/replay_upper.py.
"""
import os, sys, json
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "src"))
from flint import arb, acb, ctx
import block as BK


def legendre_vec(z, N):
    """P_0..P_{N-1} at z by the three-term recurrence, in ball arithmetic."""
    out = [arb(1)]
    if N > 1:
        out.append(arb(z))
    for n in range(1, N - 1):
        out.append(((2 * n + 1) * arb(z) * out[n] - n * out[n - 1]) / (n + 1))
    return out


class Poly:
    """f = sum_{n in degs} x_n phi_n, phi_n(u) = sqrt((2n+1)/(2L)) P_n(u/L)."""

    def __init__(self, L, degs, coeffs):
        self.L = L
        self.degs = list(degs)
        self.x = list(coeffs)
        self.maxdeg = max(self.degs)
        self.scale = [(arb(2 * n + 1) / (2 * L)).sqrt() for n in self.degs]

    def at(self, u):
        P = legendre_vec(arb(u) / self.L, self.maxdeg + 1)
        return sum((self.x[i] * self.scale[i] * P[self.degs[i]]
                    for i in range(len(self.degs))), arb(0))

    def norm2(self):
        return sum((c * c for c in self.x), arb(0))


def gl_nodes(k):
    return BK.gl_rule(k)


def autocorr(f, x, k=None):
    """g(x) = int_{x-L}^{L} f(t) f(t-x) dt, EXACT Gauss-Legendre (polynomial)."""
    L = f.L
    if k is None:
        k = f.maxdeg + 2                      # exact for degree <= 2k-1
    nodes, wts = gl_nodes(k)
    a, b = arb(x) - L, L
    hw = (b - a) / 2
    c = (a + b) / 2
    tot = arb(0)
    for q in range(k):
        t = c + hw * nodes[q]
        tot += wts[q] * f.at(t) * f.at(t - arb(x))
    return hw * tot


def outer_integrand(f, g_of, x, g0):
    """Phi(x) = 2[e^{-2x} g(0) - e^{-x/2} g(x)] / (1 - e^{-2x}) for x > 0."""
    x = arb(x)
    num = (-2 * x).exp() * g0 - (-x / 2).exp() * g_of(x)
    den = 1 - (-2 * x).exp()
    if not (den > 0):
        raise RuntimeError("outer denominator not certified positive at x=%s" % x.str(10))
    return 2 * num / den


def arch_compact(f, g_of, g0, eta="1e-40", geo_ratio=4, geo_order=12,
                 bulk_panels=32, bulk_order=32, verbose=False):
    """Certified enclosure of the Arch integral, with the near-0 piece bounded."""
    L = f.L
    twoL = 2 * L
    # --- (i) [0, eta] : crude but valid bound.
    #     |1-e^{-2x}| >= 2x e^{-2x}, and |N(x)| <= (3/2) x g(0) + x sup|g'|,
    #     so |Phi| <= [(3/2) g(0) + sup|g'|] e^{2 eta}.
    eta = arb(eta)
    fL = f.at(L)
    # sup|g'| on [0,eta] <= |f(L)| max|f| + ||f||_2 ||f'||_2 ; Markov: ||f'||<=d^2/L ||f||
    fmax = sum((abs(f.x[i]) * f.scale[i] for i in range(len(f.degs))), arb(0))
    n2 = f.norm2()
    gp = abs(fL) * fmax + n2 * (arb(f.maxdeg) ** 2 / L)
    K = ((arb(3) / 2) * g0 + gp) * (2 * eta).exp()
    head = arb(0, (eta * K).upper())
    # --- (ii) geometric panels from eta up to x0, then uniform to 2L
    x0 = arb("0.001")
    segs = []
    a = eta
    while a < x0:
        b = a * geo_ratio
        if b > x0:
            b = x0
        segs.append((a, b, geo_order))
        a = b
    h = (twoL - x0) / bulk_panels
    for i in range(bulk_panels):
        segs.append((x0 + h * i, x0 + h * (i + 1), bulk_order))
    tot = arb(0)
    for (a, b, order) in segs:
        nodes, wts = gl_nodes(order)
        hw = (b - a) / 2
        c = (a + b) / 2
        s = arb(0)
        for q in range(order):
            s += wts[q] * outer_integrand(f, g_of, c + hw * nodes[q], g0)
        tot += hw * s
    const = -(arb.const_euler() + arb.pi().log() + (1 - (-4 * L).exp()).log()) * g0
    return const + tot + head, dict(segments=len(segs), head_bound=K.str(8),
                                    eta=eta.str(6))


def Q_compact(f, sector_sign, terms, g_of=None, **kw):
    """Q(f) by the time-domain route.  terms = [(n, p)] with log n < 2L."""
    L = f.L
    if g_of is None:
        g_of = lambda x: autocorr(f, x)
    g0 = f.norm2()
    A, meta = arch_compact(f, g_of, g0, **kw)
    # pole: F(i/2) = int f e^{-u/2} du = sum_n x_n p_n   (X-coordinates)
    p = BK.pole_vector(L, f.maxdeg + 1, with_sigma=False)
    Fh = sum((f.x[i] * p[f.degs[i]] for i in range(len(f.degs))), arb(0))
    pole = 2 * sector_sign * Fh * Fh
    prime = arb(0)
    for n, pr in terms:
        prime += 2 * arb(pr).log() / arb(n).sqrt() * g_of(arb(n).log())
    return pole + A - prime, dict(pole=pole.str(20), arch=A.str(20),
                                  prime=prime.str(20), g0=g0.str(20),
                                  F_half=Fh.str(20), **meta)
