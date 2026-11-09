"""Certified UPPER bounds on lambda_min(L) from a frozen rational test vector.

lambda_min(L) = inf_{f in H, supp f in [-L,L]} Q(f)/||f||^2  <=  Q(f_x)/||f_x||^2
for ANY single admissible f, so one frozen vector gives a rigorous upper bound.
No completeness, no truncation of the function space, and no appeal to zeros is
involved: this direction is unconditional.

It says NOTHING about positivity.  A lower bound needs the whole window.

Pipeline: mpmath locates a candidate eigenvector (exploratory, uncertified) ->
the vector is FROZEN as exact rationals -> Arb ball arithmetic evaluates
x^T M x / x^T x with the certified closed form of certify.py.
"""
import json, os, sys
from fractions import Fraction
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "src"))
import mpmath as mp
import zeta_window as zw
from flint import arb
from wcore import set_prec, ball_str
from certify import Win

HERE = os.path.dirname(os.path.abspath(__file__))
FROZEN = os.path.join(HERE, "..", "certs", "frozen_vectors.json")

def candidate(Lnum, Lden, N, dps):
    """EXPLORATORY: mpmath eigenvector for the smallest eigenvalue. Uncertified."""
    mp.mp.dps = dps
    M = zw.build_matrix(mp.mpf(Lnum) / Lden, N, dps, sector=zw.EVEN)
    E, V = mp.eigsy(M)
    i = min(range(N), key=lambda k: E[k])
    return [V[r, i] for r in range(N)], E[i]

def freeze(vec, digits=40):
    """Round to exact rationals with denominator 10^digits, normalised so x[0]!=0."""
    scale = 10 ** digits
    big = max(range(len(vec)), key=lambda k: abs(vec[k]))
    v = [x / vec[big] for x in vec]
    return [Fraction(int(mp.nint(x * scale)), scale) for x in v]

def certified_bound(Lnum, Lden, x, prec):
    set_prec(prec)
    W = Win(Lnum, Lden)
    M = W.matrix(len(x))
    return W.quad_form(M, x)

def run(Lnum, Lden, N, dps=120, prec=1200, digits=60):
    label = "%d/%d" % (Lnum, Lden)
    vec, lam = candidate(Lnum, Lden, N, dps)
    x = freeze(vec, digits)
    q = certified_bound(Lnum, Lden, x, prec)
    hi = q.upper()
    print("  L=%-8s N=%-3d  mpmath lambda_min (NUMERICAL) = %s" % (label, N, mp.nstr(lam, 10)))
    print("           frozen rational vector, denominator 10^%d" % digits)
    print("           CERTIFIED enclosure of Q(f_x)/||f_x||^2 = %s" % ball_str(q, 12))
    print("           => lambda_min(%s) <= %s   [CERTIFIED]" % (label, arb(hi).str(10)))
    return dict(L=label, N=N, dps=dps, prec=prec, digits=digits,
                numerical_lambda_min=mp.nstr(lam, 16),
                certified_enclosure=q.str(20, radius=True),
                certified_upper_bound=arb(hi).str(20),
                x=[str(f) for f in x])

def main():
    print("=== CERTIFIED UPPER BOUNDS on lambda_min(L) ===")
    print("Direction: UPPER only. This is NOT a positivity certificate.")
    out = []
    for (a, b, N, dps) in [(4, 5, 16, 80), (4, 5, 24, 100), (4, 5, 28, 110),
                           (81, 100, 16, 80), (81, 100, 24, 100), (81, 100, 28, 110)]:
        out.append(run(a, b, N, dps=dps))
    os.makedirs(os.path.dirname(FROZEN), exist_ok=True)
    json.dump(out, open(FROZEN, "w"), indent=2)
    print("\nfrozen inputs + results -> %s" % os.path.relpath(FROZEN, HERE))
    # comparison against the external reference
    z_hi = arb("2.27e-17")
    best8 = min([arb(r["certified_upper_bound"]) for r in out if r["L"] == "4/5"],
                key=lambda v: float(v.mid()))
    print("\nExternal reference (Zhu v2, abstract): 8.9e-18 <= lambda_min(0.8) <= 2.27e-17")
    print("This build's certified upper bound at L=4/5 : %s" % best8.str(10))
    print("  -> %s the external upper bound"
          % ("BELOW (sharper than)" if best8 < z_hi else "ABOVE (corroborates, not sharper than)"))
    return 0

if __name__ == "__main__":
    sys.exit(main())
