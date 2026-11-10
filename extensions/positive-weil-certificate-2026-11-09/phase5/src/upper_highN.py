"""Tighten the certified upper bound on lambda_Q(4/5) along the QUADRATURE-FREE route.

CONTEXT.  The time-domain evaluator of `phase4/src/compact.py` omits the
truncation error of its outer panelled rule (RETRACTIONS.md R14), so the sharper
upper bound it produced is withdrawn.  Rather than repair that bound, this script
improves the endpoint along a route that has no quadrature to bound at all:
`src/certify.py` evaluates the window Weil form in CLOSED FORM through digamma,
Lerch Phi and Hurwitz zeta, each as an Arb ball.  A better trial function in that
basis gives a better certified endpoint, with nothing new to prove.

The even cosine basis f_k(u) = cos(k pi u / L), k = 0..N-1, is used.  Phase 1 ran
N = 16, 24, 28 and reached 2.1130e-17.  Raising N lowers the bound.

TWO STAGES, kept apart.
  * EXPLORATORY: an mpmath eigensolver locates a candidate vector.  It is
    floating point and contributes NO certified claim.
  * CERTIFIED: the candidate is frozen to exact rationals and its Rayleigh
    quotient is enclosed in Arb inside an owned precision scope.  The OUTWARD
    upper endpoint is the published bound.

The frozen vector is stored, so `verify_upper_highN()` re-certifies from exact
rational data alone and needs neither mpmath nor the eigensolver.

Run:
    PYTHONPATH=src:phase2/src:phase5/src python phase5/src/upper_highN.py [N]
"""
import math, os, sys, json, time
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
for p in (os.path.join(ROOT, "src"), os.path.join(ROOT, "phase2", "src"),
          os.path.join(ROOT, "..", "..", "src"), HERE):
    sys.path.insert(0, p)
from flint import arb, ctx
from certify import Win
from prec import owned, parse, upper as outward_upper

OUT  = os.path.join(ROOT, "phase5", "certs")
CERT = os.path.join(OUT, "upper_highN.json")

def certify_vector(x_fracs, prec):
    """CERTIFIED: enclosure of Q(f_x)/||f_x||^2 and its exact outward endpoint.

    Both are produced inside the owned scope; taking the endpoint afterwards
    would re-round it at the ambient precision (RETRACTIONS.md R13).
    """
    with owned(prec):
        W = Win(4, 5)
        M = W.matrix(len(x_fracs))
        q = W.quad_form(M, x_fracs)
        if not q.is_finite():
            raise RuntimeError("non-finite enclosure at prec=%d" % prec)
        return q, arb(q.upper())

def locate(N, dps):
    """EXPLORATORY ONLY -- floating point, no certified output."""
    import mpmath as mp
    import zeta_window as zw
    mp.mp.dps = dps
    M = zw.build_matrix(mp.mpf(4) / 5, N, dps, sector=zw.EVEN)
    E, V = mp.eigsy(M)
    i = min(range(N), key=lambda k: E[k])
    return [V[r, i] for r in range(N)], E[i]

def freeze(vec, digits):
    import mpmath as mp
    scale = 10 ** digits
    big = max(range(len(vec)), key=lambda k: abs(vec[k]))
    v = [x / vec[big] for x in vec]
    return [Fraction(int(mp.nint(x * scale)), scale) for x in v]

def verify_upper_highN(path=CERT, verbose=True):
    """Re-certify from the stored exact rationals.  No mpmath, no eigensolver."""
    rec = json.load(open(path))
    x = [Fraction(s) for s in rec["x"]]
    prec = rec["prec"]
    q, hi = certify_vector(x, prec)
    stored_hi = outward_upper(parse(rec["outward_upper_endpoint"], prec), prec)
    safe = parse(rec["safe_decimal"], prec)
    ok = (bool(stored_hi >= hi) and bool(safe >= hi)
          and parse(rec["enclosure"], prec).contains(q))
    if verbose:
        print("  re-certified from exact rationals: %s" % q.str(22, radius=True))
        print("  outward endpoint                 : %s" % hi.str(22))
        print("  stored endpoint dominates        : %s" % bool(stored_hi >= hi))
        print("  safe decimal %s dominates : %s" % (rec["safe_decimal"], bool(safe >= hi)))
    return ok, q, hi

def main():
    if "--verify" in sys.argv:
        # Re-certify the stored certificate from its exact rationals.  No
        # eigensolver, no mpmath; this is the entry point the sweep runs.
        print("=== re-certifying phase5/certs/upper_highN.json from exact rationals ===")
        t = time.time()
        ok, q, hi = verify_upper_highN()
        print("  (%.0f s)" % (time.time() - t))
        print("\n%s" % ("HIGH-N UPPER BOUND RE-CERTIFIED" if ok
                         else "RE-CERTIFICATION FAILED"))
        return 0 if ok else 1
    N    = int(sys.argv[1]) if len(sys.argv) > 1 else 128
    # The eigensolver only has to LOCATE a good trial vector; its Rayleigh
    # quotient is certified exactly afterwards, so a modest dps is enough and a
    # poor vector costs a weaker bound, never a wrong one.  lambda_min ~ 1e-17
    # against O(1) entries, so ~17 digits are lost: dps ~ 110 leaves ~90.
    dps  = 100 + N // 8
    dig  = 20 + N
    prec = 1000 + 12 * N
    print("=== certified upper bound, quadrature-free route, N = %d ===" % N)
    t0 = time.time()
    vec, lam = locate(N, dps)
    print("  [exploratory] mpmath lambda_min = %.10e   (%.0f s)" % (float(lam), time.time() - t0))
    x = freeze(vec, dig)
    q, hi = certify_vector(x, prec)
    # Shortest decimal that is CERTAINLY above the endpoint: round the endpoint
    # UP at k significant places and verify the inequality in ball arithmetic.
    safe_dec = None
    m = float(hi.mid()) * 1e17
    for k in range(3, 12):
        cand_str = "%.*fe-17" % (k, math.ceil(m * 10 ** k) / 10 ** k)
        with owned(prec):
            if bool(arb(cand_str) >= hi):
                safe_dec = cand_str
                break
    if safe_dec is None:
        raise RuntimeError("no short decimal dominates the endpoint")
    print("  [CERTIFIED]  Q/||f||^2 = %s" % q.str(22, radius=True))
    print("  [CERTIFIED]  outward upper endpoint = %s" % hi.str(22))
    print("  safe decimal = %s" % safe_dec)
    rec = dict(L="4/5", N=N, basis="even cosine f_k(u)=cos(k pi u/L), k=0..N-1",
               dps_exploratory=dps, digits=dig, prec=prec,
               numerical_lambda_min="%.16e" % float(lam),
               enclosure=q.str(30, radius=True),
               outward_upper_endpoint=hi.str(30),
               safe_decimal=safe_dec,
               route="closed form (digamma / Lerch Phi / Hurwitz zeta) in Arb balls; NO quadrature",
               independent_of="phase4/src/compact.py",
               x=[str(v) for v in x],
               elapsed_s=round(time.time() - t0, 1))
    os.makedirs(OUT, exist_ok=True)
    json.dump(rec, open(CERT, "w"), indent=2)
    print("-> phase5/certs/upper_highN.json   (%.0f s)" % (time.time() - t0))
    ok, _, _ = verify_upper_highN()
    print("\n%s" % ("HIGH-N UPPER BOUND CERTIFIED" if ok else "RE-CERTIFICATION FAILED"))
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
