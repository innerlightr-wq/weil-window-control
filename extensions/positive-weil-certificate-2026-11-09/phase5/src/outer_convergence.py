"""MEASURE how much the compact route's outer quadrature moves under refinement.

This is a DIAGNOSTIC, not a bound.  `phase4/src/compact.py::arch_compact` adds no
truncation remainder for its outer panelled Gauss-Legendre rule
(RETRACTIONS.md R14).  Refining the rule and differencing says how large the
omitted term plausibly is; it does not bound it, because the difference of two
approximations is not an enclosure of the error of either.

Result on f_200 (2026-11-10):
    baseline  geo_order=12, bulk 32x32, 94 segments   Arch = -3.14029891294194177220483427962
    refined   geo_order=20, bulk 64x48, 126 segments  Arch = -3.14029891294194177220483427962
    difference = [+/- 2.71e-35]      (i.e. indistinguishable at that level)
    relative to Q ~ 8.3e-17: [+/- 3.26e-19]

So the omitted term is plausibly far below the published radius -- which is what
one expects of a Gauss rule on an analytic integrand, and is still not a proof.

WHY A RIGOROUS BOUND IS NOT CHEAP, recorded so the next attempt starts further on.
A Bernstein-ellipse remainder on a panel needs max|Phi| on that ellipse, hence a
bound on |g(z)| off the real axis.  The naive route -- g is a polynomial of
degree D <= 2 deg f + 1 on [0,2L], so |g| <= max_[0,2L]|g| * R^D on the Bernstein
ellipse of parameter R -- is far too lossy: at D = 877 and a panel of half-width
0.025 it forces panels of width ~1e-3, i.e. ~1600 panels, and the cost explodes.

The route that should work instead is COEFFICIENT-SENSITIVE.  The frozen
minimisers concentrate on LOW Legendre orders (orders 0-10 carry essentially all
the mass), so bounding |f(w)| off the real axis directly from the frozen
coefficients,
    |f(w)| <= sum_n |c_n| sqrt((2n+1)/(2L)) * rho_w^n ,
with rho_w the Bernstein parameter of w with respect to [-L,L], is tight where
R^D is hopeless: the rho_w^n factor only matters at orders where |c_n| has
already collapsed.  That bound, fed through
    |1 - e^{-2z}| >= 1 - e^{-2 Re z}   (valid for Re z > 0)
and the standard Gauss-Legendre ellipse estimate
    |I - G_n| <= (b-a)/2 * 64 M_rho / (15 (rho^2 - 1) rho^(2n-2)) ,
is the derivation to complete.  The denominator's zeros sit at x in i*pi*Z, so
every ellipse must stay inside |Im z| < pi and exclude the origin; the geometric
panels near 0 already do the latter.

Run:
    PYTHONPATH=phase2/src:phase4/src:phase5/src python phase5/src/outer_convergence.py
"""
import os, sys, json, time
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
for p in (os.path.join(ROOT, "phase2", "src"), os.path.join(ROOT, "phase4", "src")):
    sys.path.insert(0, p)
from flint import arb, ctx
import block as BK, compact as C

OUT = os.path.join(ROOT, "phase5", "certs")

def main():
    ctx.prec = 1100
    L = arb(4) / 5
    fm = json.load(open(os.path.join(ROOT, "phase3", "certs", "frozen_minimizer.json")))
    y = [Fraction(s) for s in fm["y"]]
    degs = [2 * i for i in range(len(y))]
    x = [arb(BK.pole_sigma(degs[i]) * v.numerator) / arb(v.denominator)
         for i, v in enumerate(y)]
    f = C.Poly(L, degs, x)
    g0 = f.norm2()
    cache = {}
    def g_of(z):
        k = arb(z).str(40)
        if k not in cache:
            cache[k] = C.autocorr(f, z)
        return cache[k]

    rows = []
    for name, kw in (("baseline", dict(geo_order=12, bulk_panels=32, bulk_order=32)),
                     ("refined",  dict(geo_order=20, bulk_panels=64, bulk_order=48))):
        t = time.time()
        A, meta = C.arch_compact(f, g_of, g0, **kw)
        rows.append((name, A, meta, time.time() - t))
        print("  %-9s Arch = %s   segments=%s  (%.0f s)"
              % (name, A.str(30), meta["segments"], time.time() - t), flush=True)
    d = rows[0][1] - rows[1][1]
    print("\n  difference baseline - refined = %s" % d.str(10))
    print("  relative to Q ~ 8.3e-17       = %s" % (d / arb("8.3e-17")).str(10))
    print("\n  DIAGNOSTIC ONLY: a difference of two approximations is not an")
    print("  enclosure of the error of either.  The bound remains OPEN.")
    rec = dict(status="DIAGNOSTIC -- not a bound",
               baseline=rows[0][1].str(30), refined=rows[1][1].str(30),
               baseline_segments=rows[0][2]["segments"],
               refined_segments=rows[1][2]["segments"],
               difference=d.str(12),
               difference_relative_to_Q=(d / arb("8.3e-17")).str(12),
               elapsed_s=round(sum(r[3] for r in rows), 1),
               open_item="rigorous outer quadrature remainder for compact.arch_compact")
    os.makedirs(OUT, exist_ok=True)
    json.dump(rec, open(os.path.join(OUT, "outer_convergence.json"), "w"), indent=2)
    print("\n-> phase5/certs/outer_convergence.json")
    return 0

if __name__ == "__main__":
    sys.exit(main())
