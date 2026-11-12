"""Quadrature-free certified UPPER bound at L~ (lambda_-^2 = 4.999999999999999927).

Reuses phase5/src/upper_highN.py's own two-stage structure and the closed-form
evaluator src/certify.py::Win, with L parameterized.  Stage 1 (mpmath
eigensolver) only LOCATES a trial vector and certifies nothing; stage 2 freezes
it to exact rationals and encloses its Rayleigh quotient in Arb inside an owned
precision scope, publishing the OUTWARD upper endpoint.  No quadrature anywhere.

Consequences, by Prop. 8.17 of arXiv:2611.13129v1 (eps_lambda non-increasing):
  * at L~ itself: a genuine two-sided bracket, with the lower bound of
    certs/certificate_T400_lam5_repaired.json;
  * for EVERY lambda >= lambda_-, including lambda^2 = 5 exactly: a rigorous
    upper bound.  (It does NOT transport the lower bound upward.)
"""
import json, math, os, sys, time
from fractions import Fraction

EXT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
for p in (os.path.join(EXT, "src"), os.path.join(EXT, "phase2", "src"),
          os.path.join(EXT, "phase5", "src"),
          os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))), "src")):
    sys.path.insert(0, p)
from flint import arb
from certify import Win
from prec import owned

INV = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # phase7/
# NOTE: integrated layout -- certs/ src/ logs/ live under this phase7 directory.
# Rewritten from the session path ~/scratch/rh-connes-continuation-2026-11-12 on
# integration; no $HOME path is read.
L_NUM, L_DEN = 80471895621705018, 100000000000000000
N = int(sys.argv[1]) if len(sys.argv) > 1 else 160

def certify_vector(x_fracs, prec):
    with owned(prec):
        W = Win(L_NUM, L_DEN)
        M = W.matrix(len(x_fracs))
        q = W.quad_form(M, x_fracs)
        if not q.is_finite():
            raise RuntimeError("non-finite enclosure at prec=%d" % prec)
        return q, arb(q.upper())

def locate(N, dps):
    """EXPLORATORY ONLY -- floating point, no certified output."""
    import mpmath as mp, zeta_window as zw
    mp.mp.dps = dps
    L = mp.mpf(L_NUM) / L_DEN
    M = zw.build_matrix(L, N, dps, sector=zw.EVEN)
    E, V = mp.eigsy(M)
    i = min(range(N), key=lambda k: E[k])
    return [V[r, i] for r in range(N)], E[i]

def freeze(vec, digits):
    import mpmath as mp
    scale = 10 ** digits
    big = max(range(len(vec)), key=lambda k: abs(vec[k]))
    v = [x / vec[big] for x in vec]
    return [Fraction(int(mp.nint(x * scale)), scale) for x in v]

dps, dig, prec = 100 + N // 8, 20 + N, 1000 + 12 * N
print("=== certified upper bound at L~ = %d/%d, quadrature-free, N = %d ===" % (L_NUM, L_DEN, N))
t0 = time.time()
vec, lam = locate(N, dps)
print("  [exploratory] mpmath lambda_min = %.10e   (%.0f s)" % (float(lam), time.time() - t0))
x = freeze(vec, dig)
q, hi = certify_vector(x, prec)
safe = None
m = float(hi.mid()) * 1e18
for k in range(3, 12):
    cand = "%.*fe-18" % (k, math.ceil(m * 10 ** k) / 10 ** k)
    with owned(prec):
        if bool(arb(cand) >= hi):
            safe = cand; break
if safe is None: raise RuntimeError("no short decimal dominates the endpoint")
print("  [CERTIFIED]  Q/||f||^2            = %s" % q.str(24, radius=True))
print("  [CERTIFIED]  outward upper endpt  = %s" % hi.str(24))
print("  safe decimal (rounds UP)          = %s   dominates: %s"
      % (safe, bool(arb(safe) >= hi)))
lo = arb("7.88500976561711e-15"), arb("8.27343868254787e-18")
print("\n  bracket at lambda_- (even sector):")
print("    [ %s , %s ]" % (lo[1].str(16), hi.str(16)))
print("    width = %s   ratio = %s" % ((hi - lo[1]).str(6), (hi / lo[1]).str(8)))
json.dump(dict(window="L~ = %d/%d" % (L_NUM, L_DEN),
               lambda_sq="4.999999999999999927", N=N,
               basis="even cosine f_k(u)=cos(k pi u/L), k=0..N-1",
               dps_exploratory=dps, digits=dig, prec=prec,
               numerical_lambda_min="%.16e" % float(lam),
               enclosure=q.str(30, radius=True),
               outward_upper_endpoint=hi.str(30),
               safe_decimal=safe,
               route="closed form (digamma / Lerch Phi / Hurwitz zeta) in Arb balls; NO quadrature",
               independent_of="phase4/src/compact.py",
               monotonicity_note=("by Prop 8.17 this upper bound holds for every "
                                  "lambda >= lambda_-, including lambda^2 = 5; it does "
                                  "NOT transport the lower bound upward"),
               x=[str(v) for v in x], elapsed_s=round(time.time() - t0, 1)),
          open(os.path.join(INV, "certs/upper_lam5.json"), "w"), indent=2)
print("\nwrote certs/upper_lam5.json  (%.0f s)" % (time.time() - t0))
