#!/usr/bin/env python3
"""Targeted tests for precision ownership.  Nonzero exit on failure."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "src"))
sys.path.insert(0, os.path.join(HERE, "..", "phase2", "src"))
from flint import ctx, arb
from prec import owned, owns, PrecCache, require_tight
import constants as K

FAIL = []
def chk(name, ok, extra=""):
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name, ("  -- " + extra) if extra else ""))
    if not ok:
        FAIL.append(name)
    return ok

def main():
    print("=== precision ownership tests ===")

    print("-- P1 scope restores, nests, and survives exceptions")
    ctx.prec = 137
    with owned(500):
        inner1 = ctx.prec
        with owned(900):
            inner2 = ctx.prec
        inner3 = ctx.prec
    chk("P1a nesting: 500 / 900 / 500 and restore to 137",
        inner1 == 500 and inner2 == 900 and inner3 == 500 and ctx.prec == 137)
    try:
        with owned(777):
            raise ValueError("controlled")
    except ValueError:
        pass
    chk("P1b a controlled exception still restores the caller's context", ctx.prec == 137,
        "ctx.prec = %d" % ctx.prec)
    bad = False
    try:
        with owned(1):
            pass
    except ValueError:
        bad = True
    chk("P1c an invalid precision is rejected", bad and ctx.prec == 137)

    print("-- P2 the constants do not depend on ambient precision")
    vals = {}
    for amb in (53, 300, 1200, 2000):
        ctx.prec = amb
        c = K.constants(4, 5, 400, 440, prec=2000)
        vals[amb] = (c["beta"], c["kappa"], c["d_even"])
        if ctx.prec != amb:
            chk("P2 context restored after constants(ambient=%d)" % amb, False)
    ref = vals[2000]
    same = all(all(v[i].overlaps(ref[i]) for i in range(3)) for v in vals.values())
    chk("P2a identical enclosures from ambient 53/300/1200/2000", same,
        "beta = %s" % ref[0].str(14))
    chk("P2b caller context restored in every case",
        all(True for _ in vals))

    print("-- P3 call ORDER does not matter")
    ctx.prec = 64
    a1 = K.constants(4, 5, 200, 220, prec=1600)["kappa"]
    b1 = K.constants(4, 5, 400, 440, prec=1600)["kappa"]
    ctx.prec = 64
    b2 = K.constants(4, 5, 400, 440, prec=1600)["kappa"]
    a2 = K.constants(4, 5, 200, 220, prec=1600)["kappa"]
    chk("P3 reversing the order gives the same enclosures",
        a1.overlaps(a2) and b1.overlaps(b2),
        "kappa(200)=%s kappa(400)=%s" % (a1.str(8), b1.str(8)))

    print("-- P4 the cache is keyed by precision")
    c = PrecCache()
    lo = c.get("k", 64, lambda: arb(2).sqrt())
    hi = c.get("k", 400, lambda: arb(2).sqrt())
    chk("P4a a low-precision entry is NOT served to a high-precision caller",
        float(lo.rad()) > float(hi.rad()) and len(c) == 2,
        "radii %.3e vs %.3e, %d entries" % (float(lo.rad()), float(hi.rad()), len(c)))
    chk("P4b both entries enclose the true value", lo.overlaps(hi))

    print("-- P5 wide or indeterminate enclosures fail closed")
    bad = False
    try:
        require_tight(arb(0, 1e-3), "wide", 1e-20)
    except RuntimeError:
        bad = True
    chk("P5a a too-wide enclosure raises UNRESOLVED", bad)
    bad = False
    try:
        require_tight(arb("nan"), "nan", 1.0)
    except RuntimeError:
        bad = True
    chk("P5b a non-finite enclosure raises UNRESOLVED", bad)
    # the historical failure: S_N at too low a precision must now be refused
    bad = False
    try:
        K.constants(4, 5, 400, 440, prec=200)
    except RuntimeError as e:
        bad = "UNRESOLVED" in str(e)
    chk("P5c constants() at 200 bits is REFUSED, not silently wrong", bad,
        "this is the phase-5 defect, now fail-closed")

    print("-- P6 the decorator form demands an explicit precision")
    @owns("prec")
    def f(x, prec=None):
        return ctx.prec
    bad = False
    try:
        f(1)
    except TypeError:
        bad = True
    chk("P6 calling without prec= is a TypeError", bad and f(1, prec=321) == 321)

    print("\n%s" % ("PRECISION TESTS PASSED" if not FAIL else
                    "PRECISION TESTS FAILED: " + "; ".join(FAIL)))
    return 1 if FAIL else 0

if __name__ == "__main__":
    sys.exit(main())
