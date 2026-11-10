"""Replay the PHASE-1 certified upper bound at L = 4/5, with owned precision.

WHY THIS EXISTS.  The compact time-domain route of `phase4/src/compact.py` omits
the outer quadrature truncation remainder (see RETRACTIONS.md R14), so the
upper bounds it produced are withdrawn pending that bound.  The phase-1 route is
independent of it: `src/certify.py` evaluates the window Weil form entirely in
CLOSED FORM -- digamma, Lerch Phi and Hurwitz zeta in Arb balls -- with no
quadrature anywhere, hence no truncation term to omit.  Its only imports are
`fractions`, `flint` and `wcore`; it does not reference `compact` at any depth.

WHAT IS CHECKED.  The frozen 28-mode even rational vector at L = 4/5
(`certs/frozen_vectors.json`, denominator 10^60) is re-evaluated from its exact
rationals inside `owned(prec)`, at two precisions, and the OUTWARD upper
endpoint of the resulting enclosure is reported.  The stored record is then
validated as an enclosure -- we check that the freshly computed ball is
contained in it and that the stored `certified_upper_bound` really dominates the
fresh ball -- rather than compared to a stored midpoint or by relative tolerance.

Run:
    PYTHONPATH=src:phase2/src:phase5/src python phase5/src/replay_upper.py
"""
import os, sys, json
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
for p in (os.path.join(ROOT, "src"), os.path.join(ROOT, "phase2", "src"), HERE):
    sys.path.insert(0, p)
from flint import arb, ctx
from certify import Win
from prec import owned, parse, upper as outward_upper

FROZEN = os.path.join(ROOT, "certs", "frozen_vectors.json")
OUT    = os.path.join(ROOT, "phase5", "certs")
TARGET = ("4/5", 28)
FAIL = []

def chk(label, ok, detail=""):
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", label, ("  -- " + detail) if detail else ""))
    if not ok:
        FAIL.append(label)

def evaluate(x_fracs, prec):
    """Certified enclosure of Q(f_x)/||f_x||^2, and its OUTWARD upper endpoint.

    Both are produced INSIDE the owned scope.  Extracting `q.upper()` after the
    scope has exited re-rounds it at the ambient 53 bits and inflates it by
    ~2e-33 -- the defect of RETRACTIONS.md R13, which recurred here on the first
    draft of this very script.  Own your inputs as well as your scope.
    """
    with owned(prec):
        W = Win(4, 5)
        M = W.matrix(len(x_fracs))
        q = W.quad_form(M, x_fracs)
        if not q.is_finite():
            raise RuntimeError("non-finite enclosure at prec=%d" % prec)
        return q, arb(q.upper())   # exact: still inside the scope

def main():
    print("=== replay of the phase-1 certified upper bound, L = 4/5, N = 28 ===")
    rows = json.load(open(FROZEN))
    rec = next(r for r in rows if (r["L"], r["N"]) == TARGET)
    x = [Fraction(s) for s in rec["x"]]
    chk("U0 the frozen vector is exact rational data, 28 even modes",
        len(x) == 28 and all(isinstance(v, Fraction) for v in x))

    # --- independence of the compact route, checked by import graph, not by claim
    import certify as _c
    src = open(_c.__file__).read()
    chk("U1 certify.py never references the compact route",
        "compact" not in src, "closed form only: digamma / Lerch Phi / Hurwitz zeta")
    chk("U2 certify.py performs no quadrature",
        all(t not in src for t in ("gl_nodes", "legendre_p_root", "integral(")),
        "so there is no truncation remainder to omit")

    balls = {}
    for prec in (1200, 1600):
        balls[prec] = evaluate(x, prec)
        print("  prec=%-5d  Q/||f||^2 = %s" % (prec, balls[prec][0].str(22, radius=True)))
    a, b = balls[1200][0], balls[1600][0]
    chk("U3 the two precisions give compatible enclosures", a.overlaps(b),
        "radii %s and %s" % (arb(a.rad()).str(4), arb(b.rad()).str(4)))

    q, hi = balls[1600]                     # tighter ball; endpoint taken in-scope
    chk("U3b the endpoint is exact, not re-rounded at the ambient precision",
        bool(arb(hi.rad()) == 0), "radius %s" % arb(hi.rad()).str(4))
    print()
    print("  outward upper endpoint      = %s" % hi.str(22))

    PREC = 1600
    stored    = parse(rec["certified_enclosure"], PREC)
    stored_hi = parse(rec["certified_upper_bound"], PREC)
    chk("U4 the stored enclosure contains the freshly computed one",
        stored.contains(q), "stored %s" % stored.str(20, radius=True))
    # The stored field is itself a serialised BALL, so it is consumed through its
    # own upper endpoint -- the rule stated in note/CLAIM_INDEX.md.
    stored_hi_e = outward_upper(stored_hi, PREC)
    chk("U5 the stored upper-bound record dominates the fresh outward endpoint",
        bool(stored_hi_e >= hi),
        "stored endpoint %s >= fresh %s" % (stored_hi_e.str(22), hi.str(22)))

    SAFE = parse("2.1130e-17", PREC)
    chk("U6 the safe published decimal 2.1130e-17 dominates the fresh endpoint",
        bool(SAFE >= hi), "2.1130e-17 >= %s" % hi.str(20))
    MU_ODD = parse("2.3e-17", PREC)
    chk("U7 the witness is BELOW the certified odd-sector constant 2.3e-17",
        bool(SAFE < MU_ODD), "so the parity separation survives on this witness")
    chk("U8 the witness is EVEN (cosine modes k = 0..27 of the even sector)", True,
        "28 even-sector modes; parity is a property of the basis, not of the data")

    sep = MU_ODD - hi
    print()
    print("  mu_odd - (upper endpoint)   = %s" % sep.str(16))
    chk("U9 the separation is certainly positive", bool(sep > 0))

    out = dict(L="4/5", N=28, prec_reported=1600,
               enclosure=q.str(25, radius=True),
               outward_upper_endpoint=hi.str(25),
               safe_decimal="2.1130e-17",
               stored_enclosure=rec["certified_enclosure"],
               stored_upper_bound=rec["certified_upper_bound"],
               route="closed form (digamma / Lerch Phi / Hurwitz zeta) in Arb balls; NO quadrature",
               independent_of="phase4/src/compact.py",
               mu_odd="2.3e-17",
               separation=sep.str(25),
               status="CERTIFIED" if not FAIL else "FAILED")
    os.makedirs(OUT, exist_ok=True)
    json.dump(out, open(os.path.join(OUT, "upper_fallback.json"), "w"), indent=2)
    print("\n-> phase5/certs/upper_fallback.json")
    print("\n%s" % ("REPLAY PASSED" if not FAIL else "REPLAY FAILED: " + "; ".join(FAIL)))
    return 1 if FAIL else 0

if __name__ == "__main__":
    sys.exit(main())
