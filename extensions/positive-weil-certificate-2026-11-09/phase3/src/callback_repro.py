"""Minimal reproducer for the exception-swallowing failure in acb.integral.

OBSERVED during the phase-2 build: an entry cross-check through
`acb.integral` returned a plausible-looking but WRONG finite value
(-0.000799204100841 instead of -0.7110991221349) while the callback was in fact
raising, and only one of six entries surfaced any error at all.

This file pins the behaviour, locates the layer, and demonstrates the hardening.
"""
import sys, os
from flint import arb, acb, ctx
import flint

def banner(s):
    print("\n-- %s" % s)

def case_raises():
    """Callback raises a Python exception on SOME inputs."""
    calls = {"n": 0, "raised": 0}
    def f(z, analytic):
        calls["n"] += 1
        if float(z.real.mid()) > 0.5:
            calls["raised"] += 1
            raise ValueError("deliberate failure inside the callback")
        return z.sin()
    try:
        r = acb.integral(f, 0, 1)
        return ("returned", r, calls)
    except BaseException as e:
        return ("raised:%s" % type(e).__name__, str(e)[:70], calls)

def case_bad_op():
    """Callback performs an operation that sets an error (sqrt of a ball at 0
    through a 1/sqrt form), exactly as the phase-2 cross-check did."""
    calls = {"n": 0}
    def f(z, analytic):
        calls["n"] += 1
        return (acb(1) / z).sqrt() * z.sin()
    try:
        r = acb.integral(f, 0, 1)
        return ("returned", r, calls)
    except BaseException as e:
        return ("raised:%s" % type(e).__name__, str(e)[:70], calls)

def case_nonfinite():
    """Callback returns a NON-FINITE value -- the DOCUMENTED protocol."""
    NAN = acb(arb("nan"), arb("nan"))
    calls = {"n": 0, "nan": 0}
    def f(z, analytic):
        calls["n"] += 1
        if float(z.real.mid()) > 0.5:
            calls["nan"] += 1
            return NAN
        return z.sin()
    try:
        r = acb.integral(f, 0, 1)
        return ("returned", r, calls)
    except BaseException as e:
        return ("raised:%s" % type(e).__name__, str(e)[:70], calls)

def case_hardened():
    """The hardening used in this repository: catch inside the callback, RECORD
    the failure, return a non-finite enclosure, and REJECT the integration
    afterwards on the recorded flag -- never on the returned value alone."""
    NAN = acb(arb("nan"), arb("nan"))
    state = {"failed": False}
    def f(z, analytic):
        try:
            if float(z.real.mid()) > 0.5:
                raise ValueError("deliberate failure")
            return z.sin()
        except BaseException:
            state["failed"] = True
            return NAN
    r = None
    try:
        r = acb.integral(f, 0, 1)
    except BaseException as e:
        return ("integration raised", str(e)[:50], state)
    if state["failed"]:
        return ("REJECTED by the recorded flag", r, state)
    return ("accepted", r, state)

def main():
    ctx.prec = 200
    print("=== minimal reproducer: exception handling in acb.integral ===")
    print("python-flint %s   FLINT via python-flint bindings" % flint.__version__)
    print("python %s" % sys.version.split()[0])

    banner("A. callback raises a Python exception")
    k, v, c = case_raises()
    print("   outcome: %s" % k)
    print("   value/message: %s" % (v.str(12) if hasattr(v, "str") else v))
    print("   callback invoked %d times, raised %d" % (c["n"], c.get("raised", 0)))

    banner("B. callback performs an invalid operation (1/sqrt at a ball containing 0)")
    k, v, c = case_bad_op()
    print("   outcome: %s" % k)
    print("   value/message: %s" % (v.str(12) if hasattr(v, "str") else v))
    print("   callback invoked %d times" % c["n"])
    print("   true value of int_0^1 sin(z)/sqrt(z) dz is about 0.6205366")

    banner("C. callback returns a NON-FINITE value (documented protocol)")
    k, v, c = case_nonfinite()
    print("   outcome: %s" % k)
    print("   value: %s" % (v.str(12) if hasattr(v, "str") else v))
    print("   finite? %s   (a non-finite result is the SAFE signal)"
          % (v.is_finite() if hasattr(v, "is_finite") else "n/a"))

    banner("D. the hardening actually used here")
    k, v, c = case_hardened()
    print("   outcome: %s" % k)
    print("   recorded failure flag: %s" % c["failed"])

    banner("LAYER ATTRIBUTION")
    print("   FLINT core knows nothing about Python exceptions; it calls a C")
    print("   function pointer and continues on its return value.  The adapter")
    print("   that wraps the Python callable for that pointer is in the")
    print("   python-flint BINDING.  When the Python callback raises, the adapter")
    print("   returns to C with a CPython exception still set; FLINT then keeps")
    print("   integrating with whatever was in the output variable.  The symptom")
    print("   'returned a result with an exception set' is CPython noticing this")
    print("   only once control returns to the interpreter.")
    print("   => BINDING-level defect (callback adapter), not FLINT core, and not")
    print("      this repository's wrapper -- though the repository's wrapper is")
    print("      what must defend against it, which is case D.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
