"""Evaluate the ORIGINAL form Q on the frozen f_200, by the compact time-domain
route, and compare with the certified r_200 and the band lower bound on e_200."""
import os, sys, json, time
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "src"))
sys.path.insert(0, HERE)
from flint import arb, ctx
import block as BK, compact as C

PREC = 1100
FM = os.path.join(HERE, "..", "..", "phase3", "certs", "frozen_minimizer.json")

def load_f200(L):
    fm = json.load(open(FM))
    y = [Fraction(s) for s in fm["y"]]
    degs = [2 * i for i in range(len(y))]
    # y-coordinates -> x-coordinates: x_n = sigma_n y_n, sigma_n = (-1)^{floor(n/2)}
    x = []
    for i, f in enumerate(y):
        s = BK.pole_sigma(degs[i])
        x.append(arb(s * f.numerator) / arb(f.denominator))
    return C.Poly(L, degs, x), fm

def main():
    t0 = time.time()
    ctx.prec = PREC
    L = arb(4) / 5
    S = BK.Symbol(4, 5, 200)
    f, fm = load_f200(L)
    print("f_200: %d modes, degrees 0..%d" % (len(f.degs), f.maxdeg))
    n2 = f.norm2()
    print("  ||f||^2 = %s   (stored %s)" % (n2.str(18), fm["norm_sq"][:24]))
    print("  f(L)    = %s" % f.at(L).str(14))
    print("  c = -g'(0+) = f(L)^2 = %s" % (f.at(L) ** 2).str(14))
    g0 = C.autocorr(f, 0)
    print("  g(0)    = %s   (must equal ||f||^2)" % g0.str(18))
    if not (g0 - n2).abs_upper() < arb(10) ** -40:
        raise RuntimeError("g(0) != ||f||^2")
    print("  evaluating Q by the compact time-domain route ...", flush=True)
    cache = {}
    def g_of(x):
        k = arb(x).str(40)
        if k not in cache:
            cache[k] = C.autocorr(f, x)
        return cache[k]
    Q, meta = C.Q_compact(f, +1, S.terms, g_of=g_of)
    q = Q / n2
    r = arb(fm["r_f"])
    e = q - r
    print()
    print("  Q(f_200)            = %s" % Q.str(20))
    print("  q = Q/||f||^2       = %s   [CERTIFIED]" % q.str(20))
    print("  r_200 (phase 4)     = %s   [CERTIFIED]" % r.str(20))
    print("  e_200 = q - r       = %s   [CERTIFIED]" % e.str(20))
    band = arb(fm["e_f_lower_bound"])
    print("  band lower bound    = %s   (from [200,600])" % band.str(20))
    ok1 = bool(e >= band)
    ok2 = bool(q >= arb("2.561e-17"))
    print("  e >= band lower bound?   %s" % ok1)
    print("  q >= 2.561e-17?          %s" % ok2)
    upper = arb("2.112992612e-17")
    print("  q vs certified upper bound %s : %s" % (upper.str(12),
          "ABOVE, earlier bound retained" if q > upper else "BELOW -- improves it"))
    out = dict(norm_sq=n2.str(25), f_at_L=f.at(L).str(25), c=(f.at(L) ** 2).str(25),
               Q=Q.str(25), q=q.str(25), r_200=r.str(25), e_200=e.str(25),
               band_lower=band.str(25), e_ge_band=ok1, q_ge_2561=ok2,
               improves_upper_bound=bool(q < upper),
               elapsed_s=round(time.time() - t0, 1), **meta)
    d = os.path.join(HERE, "..", "certs")
    os.makedirs(d, exist_ok=True)
    json.dump(out, open(os.path.join(d, "q_f200.json"), "w"), indent=2)
    print("  elapsed %.0fs -> certs/q_f200.json" % (time.time() - t0))
    return 0 if (ok1 and ok2) else 1

if __name__ == "__main__":
    sys.exit(main())
