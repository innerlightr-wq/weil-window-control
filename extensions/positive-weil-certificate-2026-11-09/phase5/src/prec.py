"""Precision ownership.

THE DEFECT THIS REPAIRS.  Arb's working precision is a PROCESS-GLOBAL
(`flint.ctx.prec`).  Routines in this project set it on entry and never restored
it, so a callee silently changed the caller's context.  That bit three times:

  * phase 2, `M_sup`: a covering loop evaluated digamma on wide balls, got NaN,
    and the comparison skipped those cells -- M came out 1831 instead of 8.83.
  * phase 4, `verify_phase3`: `validate_conventions` reset ctx.prec to 300, so
    the complement constants computed afterwards came back as [+/- 7.4e+132].
  * phase 5, `run_block400`: `verify_pd_preconditioned` left ctx.prec at 1400,
    so the EVEN sector's `sector_dk` ran at 900 and gave d = [+/- 2.8e+183]
    while the ODD sector, running after, happened to get a usable value.

In every case the certificate itself was unaffected, but a dependent constant
was silently garbage, and in the first case it was *used*.

THE RULE.  A routine that needs a working precision OWNS it: it takes `prec`
explicitly, enters `with owned(prec):`, computes everything precision-dependent
inside, and the caller's context is restored on normal return AND on exception.
No routine reads the ambient precision as an input.

`flint.ctx.workprec` nests and restores on exception (verified in
`phase5/test_prec.py`); `owned` adds an assertion that it actually did, so a
future binding change cannot reintroduce the leak silently.
"""
import functools
from contextlib import contextmanager
from flint import ctx, arb

@contextmanager
def owned(bits):
    """Enter a precision scope and ASSERT the caller's context is restored."""
    if not isinstance(bits, int) or bits < 2:
        raise ValueError("precision must be an int >= 2, got %r" % (bits,))
    before = ctx.prec
    try:
        with ctx.workprec(bits):
            if ctx.prec != bits:
                raise RuntimeError("workprec did not set the precision")
            yield bits
    finally:
        if ctx.prec != before:
            ctx.prec = before
            raise RuntimeError("precision context was NOT restored (was %d, now %d)"
                               % (before, ctx.prec))

def owns(argname="prec"):
    """Decorator: run the function inside `owned(kwargs[argname])`."""
    def deco(fn):
        @functools.wraps(fn)
        def wrap(*a, **kw):
            if argname not in kw:
                raise TypeError("%s requires an explicit %s=" % (fn.__name__, argname))
            with owned(kw[argname]):
                return fn(*a, **kw)
        return wrap
    return deco


class PrecCache:
    """A cache that is keyed by the working precision as well as the arguments.

    An un-keyed cache is unsound here: a value built at 300 bits must never be
    served to a caller that asked for 2000.  Raising precision on an already
    wide ball does not recover information.
    """
    def __init__(self):
        self._d = {}

    def get(self, key, prec, build):
        k = (prec, key)
        if k not in self._d:
            with owned(prec):
                self._d[k] = build()
        return self._d[k]

    def __len__(self):
        return len(self._d)

    def clear(self):
        self._d.clear()


def require_tight(x, name, max_rad):
    """Fail closed: a value whose enclosure is too wide is UNRESOLVED, not a
    success to be read off its midpoint."""
    if not x.is_finite():
        raise RuntimeError("UNRESOLVED: %s is not finite (%s)" % (name, x))
    if float(x.rad()) > max_rad:
        raise RuntimeError("UNRESOLVED: %s radius %.3e exceeds %.3e -- refine, "
                           "do not read the midpoint" % (name, float(x.rad()), max_rad))
    return x


def parse(s, prec):
    """Parse a SERIALISED Arb ball (or decimal) inside an owned precision scope.

    Reading a stored enclosure back is precision-dependent: at the ambient 53
    bits the midpoint is rounded and the radius inflated, which silently turns a
    tight record into a ball ~1e-33 wide.  That is the defect of
    RETRACTIONS.md R13, and it recurred twice more while repairing it, so
    deserialisation gets its own entry point.  Never call `arb(stored_string)`
    directly on certificate data.
    """
    with owned(prec):
        v = arb(s) if not isinstance(s, arb) else +s
        if not v.is_finite():
            raise RuntimeError("UNRESOLVED: parsed a non-finite enclosure from %r" % (s,))
        return v


def upper(x, prec):
    """The OUTWARD upper endpoint of x, taken inside an owned scope.

    `arb(x.upper())` evaluated after a scope has exited re-rounds the endpoint at
    the ambient precision.  Here the result is exact (radius 0).
    """
    with owned(prec):
        e = arb(x.upper())
        if float(e.rad()) != 0.0:
            raise RuntimeError("endpoint extraction was not exact (radius %.3e)"
                               % float(e.rad()))
        return e


def lower(x, prec):
    """The OUTWARD lower endpoint of x, taken inside an owned scope."""
    with owned(prec):
        e = arb(x.lower())
        if float(e.rad()) != 0.0:
            raise RuntimeError("endpoint extraction was not exact (radius %.3e)"
                               % float(e.rad()))
        return e
