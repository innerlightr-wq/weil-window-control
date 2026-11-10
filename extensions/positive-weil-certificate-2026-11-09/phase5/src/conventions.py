"""Item 3: close the original-form validation ABOVE degree 0.

The phi_0 check cannot see the degree-2 pole-phase defect found in phase 3
(sigma_0 = +1).  These checks exercise phi_2 and the signed mixture
phi_0 - phi_2, compute the pole moment DIRECTLY from the physical-space
function without calling block.pole_vector, and evaluate the finite prime
shifts independently of the shared symbol helper.

Basis coordinates are stated explicitly: f = sum_n x_n phi_n with
phi_n(u) = sqrt((2n+1)/(2L)) P_n(u/L), x in the X-COORDINATES (the ones for
which F(i/2) = sum_n x_n int phi_n e^{-u/2}).

EXACTNESS OF THE INNER AUTOCORRELATION.  For f of degree d, the integrand
f(t) f(t-x) has degree 2d in the integration variable, so an m-node Gauss-
Legendre rule is exact when 2m-1 >= 2d, i.e. m >= d + 1/2, i.e. m >= d+1.
`compact.autocorr` uses m = d + 2.  Exactness means no POLYNOMIAL QUADRATURE
TRUNCATION error in exact arithmetic; node, weight, coefficient and rounding
enclosures are still carried as ball radii.  The OUTER archimedean integral
contains exponentials and keeps its own certified Bernstein remainder.
"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "src"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase4", "src"))
sys.path.insert(0, HERE)
from flint import arb, acb
import block as BK
import compact as C
from prec import owned, require_tight

CHECKS = []
def chk(name, ok, extra=""):
    CHECKS.append((name, bool(ok), extra))
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name, ("  -- " + extra) if extra else ""))
    return ok

def pole_moment_direct(f, L, prec):
    """F(i/2) = int_{-L}^{L} f(u) e^{-u/2} du by DIRECT integration of the
    physical-space function.  Does not call block.pole_vector."""
    NAN = acb(arb("nan"), arb("nan"))
    state = {"failed": None}
    def g(z, analytic):
        try:
            u = z
            P = C.legendre_vec(u / acb(L), f.maxdeg + 1) if False else None
            # evaluate f at the complex point by the same recurrence, in acb
            pp = [acb(1)]
            if f.maxdeg >= 1:
                pp.append(u / acb(L))
            for n in range(1, f.maxdeg):
                pp.append(((2 * n + 1) * (u / acb(L)) * pp[n] - n * pp[n - 1]) / (n + 1))
            val = acb(0)
            for i, d in enumerate(f.degs):
                val += acb(f.x[i] * f.scale[i]) * pp[d]
            return val * (-u / 2).exp()
        except BaseException as exc:
            state["failed"] = "%s: %s" % (type(exc).__name__, str(exc)[:60])
            return NAN
    r = acb.integral(g, -L, L)
    if state["failed"] is not None:
        raise RuntimeError("callback failed: %s" % state["failed"])
    if not r.is_finite():
        raise RuntimeError("non-finite pole moment")
    return r.real

def prime_shifts_direct(f, g_of, L, prec):
    """-sum_{log n < 2L} (2 Lambda(n)/sqrt n) g(log n), with the prime powers
    enumerated here rather than taken from block.Symbol."""
    tot = arb(0)
    rows = []
    n = 2
    twoL = 2 * L
    while arb(n).log() < twoL:
        m, p, lam = n, 2, None
        while p * p <= m:
            if m % p == 0:
                q = m
                while q % p == 0:
                    q //= p
                lam = p if q == 1 else None
                break
            p += 1
        else:
            lam = m
        if lam is not None:
            w = 2 * arb(lam).log() / arb(n).sqrt()
            gv = g_of(arb(n).log())
            tot += w * gv
            rows.append((n, lam, w, gv))
        n += 1
    return tot, rows

def run_case(name, degs, xs, L_num=4, L_den=5, prec=900):
    with owned(prec):
        L = arb(L_num) / arb(L_den)
        f = C.Poly(L, degs, [arb(v) for v in xs])
        S = BK.Symbol(L_num, L_den, 200)
        m = f.maxdeg
        nodes_needed = m + 2
        g0 = C.autocorr(f, 0)
        chk("%s: g(0) = ||f||^2" % name,
            bool((g0 - f.norm2()).abs_upper() < arb(10) ** -40),
            "deg=%d, inner GL nodes=%d (exact needs >= %d)" % (m, nodes_needed, m + 1))
        c_id = f.at(L) ** 2
        # c = -g'(0+) by finite difference as an independent read
        h = arb(10) ** -12
        fd = -(C.autocorr(f, h) - g0) / h
        chk("%s: c = -g'(0+) = f(L)^2 (identity vs difference quotient)" % name,
            bool((c_id - fd).abs_upper() < arb(10) ** -6),
            "identity %s vs FD %s" % (c_id.str(10), fd.str(10)))
        pm_direct = pole_moment_direct(f, L, prec)
        p = BK.pole_vector(L, m + 1, with_sigma=False)
        pm_helper = sum((f.x[i] * p[degs[i]] for i in range(len(degs))), arb(0))
        chk("%s: pole moment, direct physical-space integral vs pole_vector" % name,
            bool((pm_direct - pm_helper).abs_upper() < arb(10) ** -25),
            "direct %s" % pm_direct.str(16))
        cache = {}
        def g_of(x):
            k = (prec, arb(x).str(40))
            if k not in cache:
                cache[k] = C.autocorr(f, x)
            return cache[k]
        pr_direct, rows = prime_shifts_direct(f, g_of, L, prec)
        pr_helper = arb(0)
        for n, pp in S.terms:
            pr_helper += 2 * arb(pp).log() / arb(n).sqrt() * g_of(arb(n).log())
        chk("%s: prime shifts, independent enumeration vs block.Symbol" % name,
            bool((pr_direct - pr_helper).abs_upper() < arb(10) ** -30),
            "n = %s" % [r[0] for r in rows])
        Q, meta = C.Q_compact(f, +1, S.terms, g_of=g_of)
        return Q, f.norm2(), pm_direct

def main():
    print("=== higher-degree convention checks ===")
    print("basis: phi_n(u) = sqrt((2n+1)/(2L)) P_n(u/L), X-coordinates\n")
    out, live = {}, {}
    for name, degs, xs in (("phi_0", [0], [1]),
                           ("phi_2", [2], [1]),
                           ("phi_0 - phi_2", [0, 2], [1, -1])):
        print("-- %s" % name)
        Q, n2, pm = run_case(name, degs, xs)
        q = Q / n2
        print("   Q = %s   q = Q/||f||^2 = %s" % (Q.str(18), q.str(18)))
        out[name] = dict(Q=Q.str(30, radius=True), q=q.str(30, radius=True),
                         norm2=n2.str(30, radius=True),
                         pole_moment=pm.str(30, radius=True))
        live[name] = Q
    # The reference constant MUST be built inside a precision scope.  Built at the
    # ambient 53 bits it carries radius ~1e-17 and the comparison is meaningless --
    # the fourth appearance of this defect class in this session.
    with owned(900):
        ref = arb("0.080274054959269611214")
        require_tight(ref, "reference Q[0,0]", 1e-25)
        diff = live["phi_0"] - ref
    chk("phi_0 still matches the repository's exact Q[0,0]",
        bool(diff.abs_upper() < arb(10) ** -19), "diff = %s" % diff.str(8))
    d = os.path.join(HERE, "..", "certs")
    os.makedirs(d, exist_ok=True)
    json.dump(out, open(os.path.join(d, "conventions.json"), "w"), indent=2)
    bad = [n for n, ok, _ in CHECKS if not ok]
    print("\n%d checks, %d failures" % (len(CHECKS), len(bad)))
    return 1 if bad else 0

if __name__ == "__main__":
    sys.exit(main())
