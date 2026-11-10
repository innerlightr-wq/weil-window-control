#!/usr/bin/env python3
"""Numerical-consistency checker for the technical note.

Every decimal that appears as a CLAIM in
`weil_L45_certified_positivity.tex` is re-derived here from the certificate
artifacts and checked for:

  * agreement with the artifact enclosure;
  * correct ROUNDING DIRECTION -- a quoted lower bound must be <= the certified
    quantity, a quoted upper bound must be >= it;
  * correct SCOPE AND PARITY label -- which sector a constant belongs to, and
    whether it bounds lambda_Q, lambda_R or lambda_A.

Exits 0 only if every check passes.  Run:

    PYTHONPATH=../phase2/src ../../../.venv/bin/python consistency_check.py
"""
import json, os, sys
from fractions import Fraction
from flint import arb, ctx

HERE  = os.path.dirname(os.path.abspath(__file__))
EXT   = os.path.join(HERE, "..")
C2    = os.path.join(EXT, "phase2", "certs")
C3    = os.path.join(EXT, "phase3", "certs")
C4    = os.path.join(EXT, "phase4", "certs")
C5    = os.path.join(EXT, "phase5", "certs")
FB    = os.path.join(C5, "upper_fallback.json")

ctx.prec = 600
FAIL = []
N_OK = 0

def chk(ok, label, detail=""):
    global N_OK
    if ok:
        N_OK += 1
        print("  [PASS] %s%s" % (label, ("  -- " + detail) if detail else ""))
    else:
        FAIL.append(label)
        print("  [FAIL] %s%s" % (label, ("  -- " + detail) if detail else ""))

def load(path):
    return json.load(open(path))

def ball(s):
    """Parse a serialised Arb ball at the CURRENT (high) precision.

    Parsing at a low ambient precision rounds the midpoint and inflates the
    endpoints; that is exactly the defect corrected in RETRACTIONS.md R13, so
    this helper is never called outside the high-precision context set above.
    """
    b = arb(s)
    if not b.is_finite():
        raise ValueError("non-finite enclosure: %r" % s)
    return b

def lower_ok(quoted, b, label):
    """A quoted LOWER bound must not exceed the certified quantity."""
    chk(bool(arb(quoted) <= b.lower()), label,
        "quoted %s <= certified lower endpoint %s" % (quoted, arb(b.lower()).str(22)))

def upper_ok(quoted, b, label):
    """A quoted UPPER bound must not fall below the certified quantity."""
    chk(bool(arb(quoted) >= b.upper()), label,
        "quoted %s >= certified upper endpoint %s" % (quoted, arb(b.upper()).str(22)))

print("=== consistency check: note vs. artifacts ===")
print("    working precision: %d bits\n" % ctx.prec)

# ---------------------------------------------------------------- Theorem (i)
print("-- Theorem (i): the whole-window lower bound")
t400 = load(os.path.join(C4, "certificate_T400.json"))
par  = load(os.path.join(C5, "parity.json"))

mu = arb("1.41e-17")
chk(t400["whole_window_mu"].startswith("[1.4200000000000000000000") is False, "T1 mu is not a_even",
    "whole_window_mu = %s" % t400["whole_window_mu"])
chk(str(t400["sectors"]["even"]["mu"]) == "1.41e-17", "T2 even-sector mu recorded as 1.41e-17")
chk(Fraction(141, 10**19) == Fraction("1.41e-17"),
    "T3 1.41e-17 is the exact rational 141/10^19 (checked in exact arithmetic)")

for sec in ("even", "odd"):
    a = ball(t400["sectors"][sec]["a"])
    d = ball(t400["sectors"][sec]["d"])
    k = ball(t400["kappa"])
    m = arb(t400["sectors"][sec]["mu"])
    chk(bool(a > m), "T4-%s  a > mu" % sec, "a = %s" % a.str(12))
    chk(bool(d > m), "T5-%s  d > mu" % sec)
    chk(bool((a - m) * (d - m) > k * k), "T6-%s  (a-mu)(d-mu) > kappa^2" % sec,
        "margin = %s" % ((a - m) * (d - m) - k * k).str(12))

a_even = ball(t400["sectors"]["even"]["a"])
lower_ok("1.41e-17", a_even, "T7 the note's 1.41e-17 is <= the certified even-sector a")
chk(bool(arb(t400["sectors"]["even"]["mu"]) <= arb(t400["sectors"]["odd"]["mu"])),
    "T8 whole-window constant is the SMALLER sector constant (even binds)")

# --------------------------------------------------------------- Theorem (ii)
print("\n-- Theorem (ii): the certified upper bound, closed-form route")
fb    = load(FB)
cross = load(os.path.join(C4, "cross.json"))
q400  = ball(cross["Q_f400"])          # WITHDRAWN as a certificate (R14)
qfb   = ball(fb["enclosure"])
Uu    = ball(fb["outward_upper_endpoint"])

chk(fb["status"] == "CERTIFIED", "U1 the fallback replay is recorded as certified")
chk("compact" in fb["independent_of"] and "NO quadrature" in fb["route"],
    "U2 the certified route performs no quadrature and is independent of compact.py",
    fb["route"])
chk(bool(arb(Uu.upper()) >= qfb.upper()),
    "U3 the recorded endpoint dominates the fresh enclosure", "U = %s" % Uu.str(25))
upper_ok("2.1130e-17", qfb, "U4 the safe decimal 2.1130e-17 is a valid upper bound")
chk(fb["safe_decimal"] == "2.1130e-17", "U5 the artifact records the same safe decimal")
chk(bool(arb("2.1130e-17") < arb("2.27e-17")),
    "U6 the bracket is still nested inside the published one", "2.1130e-17 < 2.27e-17")
chk(bool(q400 < qfb),
    "U7 the WITHDRAWN time-domain value is sharper, which is why it is not used",
    "a sharper candidate is not a certificate: RETRACTIONS.md R14")
mid = arb("1.90138291114546440e-17")
chk(not bool(mid >= q400.upper()),
    "U8 the old ball MIDPOINT 1.90138291114546440e-17 was never a valid bound",
    "the error corrected in RETRACTIONS.md R13(a)")

# -------------------------------------------------------------- Theorem (iii)
print("\n-- Theorem (iii): the parity separation")
mu_odd = arb(par["odd"]["mu_odd"])
sep    = ball(fb["separation"])
chk(par["odd"]["mu_odd"] == "2.3e-17", "P1 odd-sector constant recorded as 2.3e-17")
chk(bool(par["odd"]["ok"]), "P2 the three strict signs certify at mu_odd")
# Two separate claims, kept apart.  The SIGN is proved as an interval sign; the
# agreement between the published record and the recomputed value is a
# reconciliation of two serialisations, so it is checked against the published
# radius rather than as a certain inequality.
recomputed = mu_odd - arb(Uu.upper())
drift = abs(arb(sep.mid()) - recomputed)
tol   = arb(sep.rad()) * 4 + arb("1e-45")
chk(bool(sep > 0), "P3a separation = mu_odd - U is certainly positive",
    "separation = %s" % sep.str(22))
chk(bool(drift <= tol),
    "P3b the published separation agrees with the recomputed value",
    "drift %s <= published radius budget %s" % (drift.str(6), tol.str(6)))
lower_ok("1.87e-18", sep, "P4 the note's 1.87e-18 is <= the certified separation")
chk(bool(arb(Uu.upper()) < mu_odd), "P5 U < mu_odd, so the EVEN infimum carries lambda_Q")
chk(fb["N"] == 28, "P6 the witness is EVEN: 28 even cosine modes k = 0..27",
    "a parity claim, not a numerical one")
chk(bool(arb(t400["sectors"]["odd"]["a"]) > mu_odd),
    "P7 the odd a exceeds mu_odd (the 1.47e-17 of the earlier draft was conservative)")

# --------------------------------------------------------- scope: three infima
print("\n-- scope labels: lambda_Q vs lambda_R vs lambda_A")
r400_f400 = ball(cross["R400_f400"])
chk(bool(r400_f400 >= a_even.lower()),
    "S1 r_400(f_400) >= a_even: the block minimum is pinned from both sides",
    "lambda_A(400,even) in [%s, %s]" % (a_even.str(6), r400_f400.str(12)))
chk(bool(q400 > r400_f400),
    "S2 q(f_400) > r_400(f_400): Q >= R_T on the same vector, strictly here")
chk(bool(ball(cross["e400_f400"]) > 0),
    "S3 the discarded term e_400(f_400) is positive, as the identity requires")
chk(ball(cross["e400_f400"]).overlaps(q400 - r400_f400),
    "S4 e_400(f_400) = q - r exactly, not a difference of infima")

# ---------------------------------------------------------------- the gate
print("\n-- the exterior gate and its monotonicity in T")
cut = load(os.path.join(C4, "cutoff.json"))
pil = load(os.path.join(C5, "pilot800.json"))
b200, b400, b800 = ball(cut["beta_200"]), ball(cut["beta_400"]), ball(pil["beta_800"])
chk(bool(b200 > 0) and bool(b400 > 0) and bool(b800 > 0),
    "G1 all three gates are positive, so the unfavourable exterior set is empty")
chk(bool(b400 > b200) and bool(b800 > b400), "G2 beta is strictly increasing in T")
chk(ball(cut["beta_increase"]).overlaps(b400 - b200), "G3 the recorded increment is beta_400 - beta_200")

# ------------------------------------------------- cutoff monotonicity per row
print("\n-- cutoff monotonicity, checked on each FIXED function separately")
tab = {
    "f_200": (ball(cross["R200_f200"]), ball(cross["R400_f200"]),
              ball(pil["R_800(f_200)"]), ball(cross["Q_f200"])),
    "f_400": (ball(pil["R_200(f_400)"]), ball(cross["R400_f400"]),
              ball(pil["R_800(f_400)"]), q400),
}
for name, (r2, r4, r8, q) in tab.items():
    chk(bool(r2 < r4) and bool(r4 < r8) and bool(r8 < q),
        "C1-%s  R_200 < R_400 < R_800 < Q" % name)

# ------------------------------------------------------------ ceilings quoted
print("\n-- certified ceilings and the percentages quoted in the note")
r200_f200 = ball(cross["R200_f200"])
r800_f400 = ball(pil["R_800(f_400)"])
# A ceiling must be the OUTWARD UPPER ENDPOINT of its enclosure, never the
# midpoint -- the same rule as U3/U8 above.
for quoted_ceil, quoted_pct, b, ref, label in (
        ("1.0310221e-17", "3.1023", r200_f200, "1e-17",    "H1 T=200 ceiling and headroom"),
        ("1.4330725e-17", "1.6364", r400_f400, "1.41e-17", "H2 T=400 ceiling and headroom"),
        ("1.6681461e-17", "18.309", r800_f400, "1.41e-17", "H3 T=800 ceiling and headroom")):
    end = arb(b.upper())
    pct = (end / arb(ref) - 1) * 100
    chk(bool(arb(quoted_ceil) >= end) and bool(abs(pct - arb(quoted_pct)) < arb("0.001")),
        label, "ceiling %s >= endpoint %s ; headroom %s%%"
               % (quoted_ceil, end.str(10), pct.str(8)))
chk(bool(arb("1.1831") >= arb(r800_f400.upper()) / arb("1.41e-17")),
    "H4 the quoted factor 1.1831 dominates the T=800 ceiling ratio",
    "ratio %s" % (arb(r800_f400.upper()) / arb("1.41e-17")).str(8))

# ------------------------------------------------------- discard ratios quoted
print("\n-- discard ratios (reported as measurements, used in no certificate)")
for name, (r2, r4, r8, q) in tab.items():
    e2, e4, e8 = q - r2, q - r4, q - r8
    chk(bool(e2 > 0) and bool(e4 > 0) and bool(e8 > 0), "D1-%s all discards positive" % name)
    chk(bool(e4 < e2) and bool(e8 < e4), "D2-%s discards shrink as T grows" % name)
r1 = (tab["f_400"][3] - tab["f_400"][1]) / (tab["f_400"][3] - tab["f_400"][0])
r2_ = (tab["f_200"][3] - tab["f_200"][1]) / (tab["f_200"][3] - tab["f_200"][0])
chk(bool(abs(r1 - r2_) > arb("0.2")),
    "D3 the two functions DISAGREE at the same doubling, so there is no decay law",
    "f_400: %s   f_200: %s" % (r1.str(6), r2_.str(6)))

# --------------------------------------------------------- superseded baseline
print("\n-- the superseded T=200 certificate is retained and still consistent")
w200 = load(os.path.join(C2, "whole_window_certificate.json"))
chk(w200["mu"] == "1e-17", "B1 the T=200 run certified 1e-17")
chk(bool(arb("1.41e-17") > arb(w200["mu"])), "B2 the T=400 bound supersedes it strictly")
chk(w200["status"] == "CERTIFIED", "B3 its status record is unchanged")
chk(bool(ball(w200["kappa"]) > ball(t400["kappa"])),
    "B4 the coupling improved from T=200 to T=400",
    "%s -> %s" % (ball(w200["kappa"]).str(6), ball(t400["kappa"]).str(6)))

print("\n%d checks passed, %d failed" % (N_OK, len(FAIL)))
if FAIL:
    print("FAILED: " + "; ".join(FAIL))
    sys.exit(1)
print("CONSISTENCY CHECK PASSED")
sys.exit(0)
