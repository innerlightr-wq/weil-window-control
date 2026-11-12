"""Rigorous two-sided bound on eps(lambda_-)/eps(e^{4/5}) from the four certified
   endpoints alone.  No prolate denominator, no parameter transfer, no transport.

   mu_- / U_0  <=  eps_-/eps_0  <=  U_- / mu_0
"""
from fractions import Fraction as F
from flint import arb, ctx
ctx.prec = 400

# the four PACKAGED endpoints (certs/schur_margins.json, certs/upper_lam5.json,
# phase5/certs/upper_highN.json).  Lower bounds are published mu; upper bounds
# are outward endpoints.
mu0 = F("1.43300048827856e-17")                               # lower, L = 4/5
U0  = F("1.67312553407068186807726837736e-17")                # upper, L = 4/5
mum = F("8.27343868254787e-18")                               # lower, lambda_-
Um  = F("9.47476951083049032849709e-18")                      # upper, lambda_-

lo_exact = mum / U0
hi_exact = Um / mu0
print("exact rational arithmetic:")
print("  mu_-/U_0 = %s" % lo_exact)
print("           = %.20f" % float(lo_exact))
print("  U_-/mu_0 = %s" % hi_exact)
print("           = %.20f" % float(hi_exact))

def round_down(x, k):
    from decimal import Decimal, getcontext, ROUND_FLOOR
    getcontext().prec = 40
    d = Decimal(x.numerator) / Decimal(x.denominator)
    return d.quantize(Decimal(1).scaleb(-k), rounding=ROUND_FLOOR)
def round_up(x, k):
    from decimal import Decimal, getcontext, ROUND_CEILING
    getcontext().prec = 40
    d = Decimal(x.numerator) / Decimal(x.denominator)
    return d.quantize(Decimal(1).scaleb(-k), rounding=ROUND_CEILING)

for k in (3, 4, 6):
    lo = round_down(lo_exact, k); hi = round_up(hi_exact, k)
    ok = (F(str(lo)) <= lo_exact) and (F(str(hi)) >= hi_exact)
    print("\n  outward to %d dp:  %s <= eps_-/eps_0 <= %s   direction-safe: %s"
          % (k, lo, hi, ok))
    assert ok
    print("     decrease between %s%% and %s%%"
          % (round_down(F(1) - F(str(hi)), k) * 100, round_up(F(1) - F(str(lo)), k) * 100))

# cross-check in Arb balls, consuming each endpoint through the correct side
lo_b = arb(str(mum.numerator)) / arb(str(mum.denominator)) / (arb(str(U0.numerator)) / arb(str(U0.denominator)))
hi_b = arb(str(Um.numerator)) / arb(str(Um.denominator)) / (arb(str(mu0.numerator)) / arb(str(mu0.denominator)))
print("\nArb cross-check:  [%s , %s]" % (lo_b.str(12), hi_b.str(12)))
print("  3-dp outward [0.494, 0.662] contains it: %s"
      % (bool(arb("0.494") <= lo_b) and bool(arb("0.662") >= hi_b)))

# what the prolate scaling law would predict, and whether this can see it
rung = arb("1.2983880087183625604e-18") / arb("2.2391127080494187067e-18")
print("\nprolate-law prediction delta_4(5)/delta_4(4.953) = %s" % rung.str(10))
print("  inside the certified bracket: %s"
      % (bool(rung > lo_b) and bool(rung < hi_b)))
w = hi_b - lo_b
print("  bracket width = %s ; the law sits %s of the way across"
      % (w.str(6), ((rung - lo_b) / w).str(4)))
print("\n  To TEST the law to +/-5%% the ratio bracket must be narrower than")
print("  %s, i.e. each window's own bracket ratio below about %s"
      % ((rung * arb("0.1")).str(4), (1 + rung * arb("0.1") / 2 / rung).str(5)))
print("  (currently 1.1452 at lambda_- and 1.1676 at L=4/5)")
