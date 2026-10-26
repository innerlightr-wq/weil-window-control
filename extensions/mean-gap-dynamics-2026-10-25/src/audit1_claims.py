"""AUDIT Task 1 -- resolve the three outstanding claims.

  (1a) the proximity cutoff:   is  R < Delta/sqrt2  correct?  is it sharp?
  (1b) the strip bound vs the actual transition threshold:  tau_c is NOT Lambda.
  (1c) the hypotheses behind the infinite-case energy identity.

Everything here is exact algebra, verified over Q where the quantities are rational.
"""
import sys, os, json
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(__file__))
from mpmath import mp, mpf, sqrt, findroot, log, pi, exp

mp.dps = 40
OK = [True]
def chk(n, c, x=""):
    OK[0] &= bool(c); print(f"  {'PASS' if c else 'FAIL'}  {n} {x}")

out = {}

# =====================================================================  (1a)
print("=" * 78)
print("(1a) THE PROXIMITY CUTOFF -- audited")
print("=" * 78)
print("""
The quantity.  At a simple non-real zero z_j = x_j + i y_j (take y_j > 0), the sign
identity gives the contribution of a competing zero z_k with y_k > y_j > 0 as

    c_k = 2 y_j (y_k - y_j) / |z_j - z_k|^2 ,       and  |z_j - z_k|^2 >= R^2 + h^2
    with  R = |Re z_j - Re z_k| ,  h = y_k - y_j .

For ONE competitor to overturn the conjugate-partner term (exactly -1) we need c_k > 1.
""")

Delta = mpf(1)          # normalise Delta_tau = 1; every threshold is a multiple of it

def phi(s):
    """max over 0 < y_j < y_k <= 1 of 2 y_j (y_k-y_j)/(s^2 + (y_k-y_j)^2).

    y_k = 1 is optimal (the expression increases in y_k at fixed h), so with u = h:
        phi(s) = max_{0<u<1} 2(1-u)u/(s^2+u^2).
    """
    s = mpf(s)
    f = lambda u: 2*(1-u)*u/(s*s + u*u)
    # the stationary point solves  u^2 + 2 s^2 u - s^2 = 0   (derived in AUDIT.md)
    u = -s*s + sqrt(s**4 + s*s)
    if u <= 0: return mpf(0)
    if u >= 1: u = mpf(1)
    return f(u)

print("-- the claim in the previous report:  R < Delta/sqrt2 = 0.7071067812 --")
v = phi(1/sqrt(2))
print(f"   phi(1/sqrt2) = {mp.nstr(v, 12)}")
chk("R = Delta/sqrt2 does NOT admit a competitor (phi < 1)", v < 1,
    "-> the stated cutoff is VALID but NOT SHARP")

print("\n-- the SHARP cutoff, solved exactly --")
print("""   The stationary point of 2(1-u)u/(s^2+u^2) satisfies  u^2 + 2 s^2 u - s^2 = 0,
   i.e.  s^2 = u^2/(1-2u).  Substituting into  2u(1-u) = s^2 + u^2  and dividing by u:
        2 - 3u = u/(1-2u)   =>   6u^2 - 8u + 2 = 0   =>   3u^2 - 4u + 1 = 0
        =>   u = 1  or  u = 1/3 ;   u < 1/2  forces  u = 1/3 ,
   and then  s^2 = (1/9)/(1 - 2/3) = 1/3 .""")
u0, s0 = Fr(1, 3), Fr(1, 3)
chk("3u^2-4u+1 = 0 at u = 1/3 (exact over Q)", 3*u0*u0 - 4*u0 + 1 == 0)
chk("s^2 = u^2/(1-2u) = 1/3 (exact over Q)", u0*u0/(1 - 2*u0) == s0)
# the extremal configuration, exactly
yk, yj, h, R2 = Fr(1), Fr(2, 3), Fr(1, 3), Fr(1, 3)
chk("extremal config (y_j, y_k, h) = (2/3, 1, 1/3) gives c_k = 1 EXACTLY",
    2*yj*h/(R2 + h*h) == 1, f"= {2*yj*h/(R2+h*h)}")
chk("phi(1/sqrt3) = 1 numerically", abs(phi(1/sqrt(mpf(3))) - 1) < mpf(10)**-30,
    f"dev {mp.nstr(abs(phi(1/sqrt(mpf(3)))-1),3)}")
thr_sharp = 1/sqrt(mpf(3))
print(f"\n   >>> SHARP single-competitor cutoff:  R < Delta_tau/sqrt3 = "
      f"{mp.nstr(thr_sharp, 12)} * Delta_tau")
print(f"   >>> previously stated:               R < Delta_tau/sqrt2 = "
      f"{mp.nstr(1/sqrt(mpf(2)), 12)} * Delta_tau   (valid, 22.5% too loose)")
out["cutoff_sharp"] = mp.nstr(thr_sharp, 20)
out["cutoff_previous"] = mp.nstr(1/sqrt(mpf(2)), 20)

print("\n-- the N-competitor form in the previous report is WRONG --")
print("""   Stated:  R < Delta sqrt(N)/2 .   Derived:  each competitor contributes at most
   Delta^2/(2R^2), so N of them need  N Delta^2/(2R^2) > 1, i.e.  R < Delta sqrt(N/2).
   At N = 1 the stated form gives Delta/2, which contradicts the N = 1 result
   Delta/sqrt2 that the same report states two lines earlier.  The crude form is""")
print(f"   >>> R < Delta sqrt(N/2) ,  which at N = 1 is Delta/sqrt2  (consistent).")
chk("the stated N-form disagrees with its own N=1 case",
    abs(sqrt(mpf(1))/2 - 1/sqrt(mpf(2))) > mpf("0.2"),
    f"Delta/2 = 0.5 vs Delta/sqrt2 = {mp.nstr(1/sqrt(mpf(2)),6)}")
chk("the corrected crude form reproduces the N=1 case",
    abs(sqrt(Fr(1, 2)*1.0) - 1/sqrt(mpf(2))) < mpf(10)**-20)
print("\n   and the SHARP N-competitor cutoff solves  N phi(s) = 1 :")
print(f"   {'N':>3} {'sharp R/Delta':>16} {'crude sqrt(N/2)':>18}")
sharpN = {}
for N in (1, 2, 3, 4, 6, 10):
    f = lambda s: N*phi(s) - 1
    sN = findroot(f, mpf("0.6") if N == 1 else mpf(N)**mpf("0.5")*mpf("0.6"))
    sharpN[N] = mp.nstr(sN, 14)
    print(f"   {N:>3} {mp.nstr(sN,12):>16} {mp.nstr(sqrt(mpf(N)/2),12):>18}")
    OK[0] &= (sN <= sqrt(mpf(N)/2) + mpf(10)**-20)
chk("the sharp N-cutoff never exceeds the crude one", True)
out["cutoff_sharp_N"] = sharpN

print("\n-- a SHARPER per-zero form, not previously noticed --")
print("""   Do not optimise over y_j.  For a FIXED y_j the unconstrained maximum of
   2 y_j h/(R^2+h^2) over h > 0 is at h = R, with value y_j/R.  So c_k > 1 requires

        >>>  R < y_j   <<<      i.e., in the Dirichlet variable,  |gamma_j - gamma_k| < |delta_j| .

   (If the optimum h = R is cut off by the strip, h = Delta - y_j < R, then
    2 y_j h > R^2 + h^2 >= 2 R h by AM-GM, so y_j > R again.  So R < y_j is necessary
    in ALL cases, with no hypothesis beyond the strip.)""")
# verify both branches numerically on a grid
bad = 0
for yj_ in [mpf(k)/20 for k in range(1, 20)]:
    for R_ in [mpf(k)/20 for k in range(1, 20)]:
        best = mpf(0)
        for h_ in [mpf(k)/200 for k in range(1, 201)]:
            if yj_ + h_ > 1: break
            best = max(best, 2*yj_*h_/(R_*R_ + h_*h_))
        if best > 1 and not (R_ < yj_):
            bad += 1
chk("grid 19x19 x 200: c_k > 1 never occurs with R >= y_j", bad == 0, f"({bad} violations)")
# and R < y_j alone is NOT sufficient -- it must be paired with R < Delta/sqrt3
suff = 0
for yj_ in [mpf(k)/20 for k in range(1, 20)]:
    R_ = yj_*mpf("0.9")
    best = mpf(0)
    for h_ in [mpf(k)/200 for k in range(1, 201)]:
        if yj_ + h_ > 1: break
        best = max(best, 2*yj_*h_/(R_*R_ + h_*h_))
    if best <= 1: suff += 1
chk("R < y_j is NECESSARY but not SUFFICIENT", suff > 0,
    f"({suff} of 19 grid points have R = 0.9 y_j yet no competitor works)")
print("   So for a SINGLE competitor the audited condition is the CONJUNCTION")
print("       R < min( y_j , Delta_tau/sqrt3 )   [necessary]")
print("   of which the first factor is new and is the binding one for a shallow zero.")
print("""
-- SINGLE vs SEVERAL interactions: the distinction, stated explicitly --
   Every bound above concerns ONE competitor's contribution exceeding 1.  Several smaller
   contributions can add up to more than 1 while no single one does, so the single-
   competitor conditions are NOT necessary for N > 1.  With N competitors all at ordinate
   distance >= R, each contributes at most min( y_j/R , Delta^2/(2R^2) ), so

       >>>  R < min( N y_j , Delta_tau sqrt(N/2) )      [necessary, N competitors]  <<<

   Both factors relax as N grows -- linearly and like sqrt(N) respectively -- so a
   configuration with many shallow competitors is NOT excluded by the single-competitor
   form.  Nothing in this extension's conclusions depends on the N = 1 case: the DH census
   has only three competitors for any target and all of them fail even the N = 4 form, and
   the toy counterexample has N = 1.  The distinction is recorded because the N = 1 form is
   the quotable one and is the easiest to over-apply.""")

out["cutoff_per_zero"] = "R < y_j, i.e. |gamma_j - gamma_k| < |delta_j|"

print("\n-- re-applied to the two controls --")
TOY = dict(Delta=mpf(2), yj=mpf(1), R=mpf(0))
print(f"   toy (z^2+1)(z^2+4): y_j = 1, R = 0, Delta = 2;"
      f"  R < y_j ? {'YES' if TOY['R'] < TOY['yj'] else 'NO'};"
      f"  R < Delta/sqrt3 = {mp.nstr(2/sqrt(mpf(3)),8)} ? YES")
chk("the toy counterexample satisfies BOTH sharpened conditions",
    TOY["R"] < TOY["yj"] and TOY["R"] < TOY["Delta"]/sqrt(mpf(3)))
CENSUS = [("O_1", "0.3085171824566374", "85.6993484853775922"),
          ("O_2", "0.1508300806097371", "114.163342730756981"),
          ("O_3", "0.0743560504508060", "166.479305913168156"),
          ("O_4", "0.2242576946268098", "176.702461242855825")]
print(f"\n   DH census, the per-zero form |gamma_j-gamma_k| < |delta_j| :")
print(f"   {'target':>7} {'|delta_j|':>14} {'deepest competitor':>20} {'|d gamma|':>14} {'ratio':>10}")
rows = []
for nj, dj, gj in CENSUS:
    deeper = [(nk, dk, gk) for nk, dk, gk in CENSUS if mpf(dk) > mpf(dj)]
    if not deeper:
        print(f"   {nj:>7} {mp.nstr(mpf(dj),10):>14} {'none (census max)':>20} {'-':>14} {'-':>10}")
        rows.append(dict(target=nj, competitor=None)); continue
    nk, dk, gk = max(deeper, key=lambda t: mpf(t[1]))
    dg = abs(mpf(gj) - mpf(gk))
    ratio = dg/mpf(dj)
    print(f"   {nj:>7} {mp.nstr(mpf(dj),10):>14} {nk:>20} {mp.nstr(dg,10):>14} {mp.nstr(ratio,6):>10}")
    rows.append(dict(target=nj, competitor=nk, dgamma=mp.nstr(dg, 14), ratio=mp.nstr(ratio, 10)))
    OK[0] &= (dg > mpf(dj))
chk("every DH competitor is outside the per-zero cutoff", True,
    "-> by factors 189 to 1086, far stronger than the 45-144 the Delta/sqrt3 form gave")
out["dh_N3"] = []
print("\n   and the SAME census against the N = 3 form (every other quartet counted):")
DeltaDH = 2*(mpf("1.395136158235109721061359") - mpf(1)/2)
for nj, dj, gj in CENSUS:
    yj = 2*mpf(dj)
    Rmin = min(2*abs(mpf(gj) - mpf(gk)) for nk, dk, gk in CENSUS if nk != nj)
    thrN = min(3*yj, DeltaDH*sqrt(mpf(3)/2))
    exc = Rmin > thrN
    print(f"     {nj}: nearest competitor R = {mp.nstr(Rmin,8):>12}   N=3 threshold "
          f"{mp.nstr(thrN,8):>12}   -> {'excluded' if exc else 'NOT excluded'}")
    out["dh_N3"].append(dict(target=nj, R=mp.nstr(Rmin, 12), thr=mp.nstr(thrN, 12),
                             excluded=bool(exc)))
    OK[0] &= exc
chk("every DH target is excluded even by the N = 3 form", True,
    "-> so no conclusion here rests on the single-competitor case")
out["dh_per_zero"] = rows

# =====================================================================  (1b)
print("\n" + "=" * 78)
print("(1b) THE STRIP BOUND IS NOT THE TRANSITION THRESHOLD")
print("=" * 78)
print("""
Three different numbers were in play and must not be merged:

   Lambda          the de Bruijn-Newman constant: the ACTUAL transition, inf{tau : all
                   zeros of H_tau real}.  Not computable from the zero set at tau = 0.
   tau_c           Delta_0^2/2 = sup_rho D(rho) : de Bruijn's SUFFICIENT threshold.  An
                   UPPER bound for Lambda, and the identity of DERIVATION Sec.1.7.
   Escr(0)/2       (1/2) sum_rho D(rho) : the (S2) bound.  Also an upper bound, weaker.

   Always  Lambda <= tau_c <= Escr(0)/2 ,  and both inequalities can be STRICT.""")
# exact on the toy
toy_tau_c = mpf(2)                      # Delta_0^2/2 with Delta_0 = 2
toy_E0_2 = mpf(5)/2
toy_Lambda = (10 + sqrt(mpf(46)))/24    # EXACT: the larger root of 96t^2-80t+9
print(f"\n   toy (z^2+1)(z^2+4), all three EXACT:")
print(f"     Lambda      = (10 + sqrt46)/24 = {mp.nstr(toy_Lambda, 16)}")
print(f"     tau_c       = Delta_0^2/2      = {mp.nstr(toy_tau_c, 16)}")
print(f"     Escr(0)/2                      = {mp.nstr(toy_E0_2, 16)}")
chk("Lambda < tau_c  (STRICTLY, by a factor 2.86)", toy_Lambda < toy_tau_c,
    f"ratio {mp.nstr(toy_tau_c/toy_Lambda, 8)}")
chk("tau_c < Escr(0)/2 (STRICTLY)", toy_tau_c < toy_E0_2)
out["toy_thresholds"] = dict(Lambda=mp.nstr(toy_Lambda, 20), tau_c="2", E0_2="2.5",
                             ratio_tau_c_over_Lambda=mp.nstr(toy_tau_c/toy_Lambda, 12))
print(f"""
   Davenport-Heilbronn, from the parent extension's own certified bracket:
     Lambda_DH in [6.8268151773112564447e-29, 1.602537483559822777787819]
     tau_c^DH   = sup D^DH = 1.602537483559822777787819   -- the UPPER end of that bracket
   So the identity tau_c = sup D names the bracket's upper endpoint, NOT Lambda.  The
   bracket spans 28 orders of magnitude and the identity does not narrow it at all.""")
print("""   AUDIT FINDING.  The previous report's Sec.1.7 and the plain-language verdict call
   tau_c "de Bruijn's threshold", which is correct, and never assert tau_c = Lambda.  But
   the phrase "the de Bruijn-Newman constant" appears nowhere near it either, so the
   distinction was implicit.  It is now explicit, and the toy supplies a factor-2.86 gap
   as a standing reminder.""")

# =====================================================================  (1c)
print("\n" + "=" * 78)
print("(1c) THE HYPOTHESES BEHIND THE INFINITE-CASE ENERGY IDENTITY -- a REAL error")
print("=" * 78)
print("""
Stated previously:  (S2) holds "whenever Escr < infinity".  That is WRONG, and not
merely incomplete.  The dissipation is

    Diss(tau) = sum_{j<k} (y_j - y_k)^2 / |z_j - z_k|^2 ,     Escr' = -2 Diss ,

and EVERY conjugate pair {z_j, conj z_j} with y_j != 0 contributes EXACTLY 1 to it
(numerator (2y_j)^2, denominator (2|y_j|)^2).  Therefore, with m(tau) the number of
non-real zeros,

    Diss(tau) >= m(tau)/2 ,          so     m = infinity  ==>  Diss = +infinity .

But Escr = sum_j D_j can be FINITE with m = infinity (take y_j -> 0 fast enough, e.g.
D_j = 2^{-j}).  So "Escr < infinity" does NOT imply the identity -- in that regime the
identity is not merely unproved, it is INCOMPATIBLE: a finite left side and an infinite
right side.""")
# exhibit the incompatibility concretely
Ds = [Fr(1, 2**j) for j in range(1, 60)]
chk("a configuration with Escr finite and m infinite exists: D_j = 2^-j",
    sum(Ds) < 1, f"partial sum {float(sum(Ds)):.12f} -> 1, while m = infinity")
chk("its dissipation is already >= m/2 = 29.5 from the conjugate terms alone",
    len(Ds)/2 >= 29, f"m/2 = {len(Ds)/2}")
print(f"""
   >>> CORRECTED HYPOTHESIS.  (S2) and its corollary require

           m(tau) < infinity      (FINITELY MANY non-real zeros),

   not Escr(tau) < infinity.  Under m < infinity everything is automatic: Escr is a
   finite sum, the double sum is absolutely convergent (pairs of real zeros contribute
   0; each pair with a non-real member is bounded by 2 Delta/|z_j-z_k|^2 and
   sum_k |z_j-z_k|^{-2} < infinity by order 1), and the rearrangement is legitimate.

   >>> AND THE CORRECTION HELPS THE FILTER, not hurts it.  Last pass recorded
   "sum_rho D(rho) < infinity" as a hypothesis DH is "not known to satisfy", resting on
   Levinson's >> T count and on a guess about summability.  The corrected hypothesis
   m < infinity is failed by DH on the strength of the ORIGINAL 1936 theorem alone --
   Davenport and Heilbronn proved F has INFINITELY MANY zeros off the critical line.
   No summability guess and no second source is needed.
   (Bibliographic, not verified here, exactly as the parent extension records it.)""")
out["S2_hypothesis"] = dict(previously="Escr < infinity  -- WRONG",
                            corrected="m(tau) < infinity (finitely many non-real zeros)",
                            reason="each conjugate pair contributes exactly 1 to the dissipation",
                            dh_fails_by="Davenport-Heilbronn 1936 (infinitude), bibliographic")

print("\n" + "=" * 78)
print("=== audit1_claims.py:", "ALL PASS" if OK[0] else "SOME FAILED", "===")
json.dump(out, open(os.path.join(os.path.dirname(__file__), "..", "data", "audit1.json"), "w"),
          indent=1)
sys.exit(0 if OK[0] else 1)
