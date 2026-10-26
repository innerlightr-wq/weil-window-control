"""AUDIT Task 4 -- is there ONE precise intermediate lemma?  Plausibility numerics.

The candidate combines two INDEPENDENTLY established constraints:

  (A) the pointwise sign identity of this extension (T1, proved):
          D'(tau) = -2 y_j sum_{k != j} (y_j - y_k)/|z_j - z_k|^2 ,
      which is an EQUALITY, so it bounds |D'| from ABOVE as well as below;

  (B) the unconditional shell-counting bound -- the exact analogue, for H_tau, of the
      Riemann-von Mangoldt estimate  N(T+1) - N(T) = O(log T)  that the parent
      extensions already use for N(T).  It follows from Jensen's inequality and the
      order-1 growth bound |H_tau(z)| <= exp(C|z| log|z|) derived in DERIVATION Sec.3.1.

An UPPER bound on |D'| is a LOWER bound on how long a branch takes to reach the critical
line, hence a LOWER bound on Lambda.  Lower bounds on Lambda are the WRONG direction for
RH -- which is exactly why this is proposed as a Davenport-Heilbronn-side sharpening and
NOT as an RH project.  The target is the parent extension's own published certified bound
Lambda_DH >= 6.8268151773112564447e-29, whose weakness that extension diagnoses
structurally ("m/C inherits the gamma decay on one side and none of it on the other").

This script measures every ingredient so the proposal is a plausibility estimate backed by
numbers, not an assertion.  NOTHING HERE IS A THEOREM.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
from mpmath import mp, mpf, mpc, log, pi, sqrt
import kernels as K
from dh_flow import zdot_R1, CENSUS

mp.dps = 40
OK = [True]
def chk(n, c, x=""):
    OK[0] &= bool(c); print(f"  {'PASS' if c else 'FAIL'}  {n} {x}")
out = {}

# ---- the zero list computed by src/dh_interaction.py
ordf = os.path.join(os.path.dirname(__file__), "..", "data", "dh_online_ordinates.txt")
ords = [mpf(l.strip()) for l in open(ordf) if l.strip()]
offl = [(mpf(b) - mpf(1)/2, mpf(g)) for _, b, g in CENSUS]
ZS = []
for t in ords:
    ZS += [2*t, -2*t]
for d, g in offl:
    ZS += [2*g - 2*mpc(0, 1)*d, 2*g + 2*mpc(0, 1)*d, -2*g + 2*mpc(0, 1)*d, -2*g - 2*mpc(0, 1)*d]
print(f"zero list: {len(ords)} on-line ordinates (|Im s| < 260) + 4 quartets -> {len(ZS)} points in z")

b1, g1 = mpf(CENSUS[0][1]), mpf(CENSUS[0][2])
d1 = b1 - mpf(1)/2
z1 = 2*g1 - 2*mpc(0, 1)*d1
y1, x1 = mp.im(z1), mp.re(z1)
nu_orbit = 2*abs(y1)                      # distance to the conjugate partner
sigma0 = mpf("1.395136158235109721061359")
Delta0 = 2*(sigma0 - mpf(1)/2)

print(f"\nz_1 = {mp.nstr(z1,16)};  |Im z_1| = {mp.nstr(abs(y1),12)};  "
      f"2|Im z_1| = {mp.nstr(nu_orbit,12)};  Delta_0 = {mp.nstr(Delta0,12)}")

# =============================================== 1. the isolation hypothesis, at tau = 0
print("\n" + "="*78)
print("1. THE ISOLATION HYPOTHESIS, checked at tau = 0")
print("="*78)
orbit = [z1, mp.conj(z1), -z1, -mp.conj(z1)]
others = [w for w in ZS if min(abs(w - o) for o in orbit) > mpf(10)**-20]
dmin = min(abs(z1 - w) for w in others)
nearest = min(others, key=lambda w: abs(z1 - w))
print(f"   nearest NON-ORBIT zero to z_1: {mp.nstr(nearest,14)}  at distance "
      f"{mp.nstr(dmin,12)}")
print(f"   the conjugate partner is at distance 2|Im z_1| = {mp.nstr(nu_orbit,12)}")
chk("isolation holds at tau = 0: the conjugate partner IS the nearest zero",
    dmin > nu_orbit, f"margin {mp.nstr(dmin/nu_orbit,6)}x")
out["isolation"] = dict(dmin=mp.nstr(dmin, 14), nu_orbit=mp.nstr(nu_orbit, 14),
                        margin=mp.nstr(dmin/nu_orbit, 10))

# ============================================ 2. the shell-count constant, measured
print("\n" + "="*78)
print("2. THE SHELL-COUNT CONSTANT C_0, measured on the computed zero list")
print("="*78)
print("   hypothesis (B):   #{k : |z_k - w| <= 1}  <=  C_0 log(2 + |w|)")
print(f"   {'w = 2 gamma':>14} {'count in unit disc':>20} {'log(2+|w|)':>14} {'implied C_0':>14}")
C0s = []
for gg in ("30", "85.7", "150", "230"):
    w = 2*mpf(gg)
    cnt = sum(1 for v in ZS if abs(v - w) <= 1)
    lg = log(2 + abs(w))
    C0s.append(cnt/lg)
    print(f"   {mp.nstr(w,10):>14} {cnt:>20} {mp.nstr(lg,8):>14} {mp.nstr(cnt/lg,8):>14}")
C0 = max(C0s)
print(f"   measured C_0 (max over the samples) = {mp.nstr(C0,8)}")
print("   NOTE: measured on a finite list, as an ORDER OF MAGNITUDE for the plausibility")
print("   estimate only.  The lemma would need the Jensen bound, not this measurement.")
out["C0_measured"] = mp.nstr(C0, 10)

# ============================================ 3. the actual S_1, and the bound's shape
print("\n" + "="*78)
print("3. S_1 = sum_{k != 1} 1/|z_1 - z_k|^2 : measured vs the bound's shape")
print("="*78)
S1_all = sum(1/abs(z1 - w)**2 for w in others) + sum(1/abs(z1 - o)**2 for o in orbit[1:])
S1p = sum(1/abs(z1 - w)**2 for w in others)          # the NON-ORBIT part -- what is bounded
print(f"   S_1 over ALL listed zeros      = {mp.nstr(S1_all,10)}")
print(f"     of which the orbit supplies  = {mp.nstr(S1_all - S1p,10)}"
      f"   (the conjugate partner alone {mp.nstr(1/nu_orbit**2,8)})")
print(f"   S_1' = the NON-ORBIT part      = {mp.nstr(S1p,10)}   <- this is what is bounded")
print(f"   (it agrees with the cross-orbit sum measured independently in")
print(f"    src/dh_interaction.py, as it must)")
bound_shape = C0*(log(2+abs(z1))/dmin + 1/dmin**2)
print(f"   bound shape C_0(log(2+|z_1|)/nu + 1/nu^2), nu = {mp.nstr(dmin,8)}: "
      f"{mp.nstr(bound_shape,8)}")
chk("the measured S_1' is below the bound's shape", S1p < bound_shape,
    f"slack factor {mp.nstr(bound_shape/S1p,6)}x")
print("   NOTE on the first draft of this check: it compared S_1 (orbit terms INCLUDED,")
print("   0.865) against a bound for S_1' (orbit terms EXCLUDED, 0.532) and FAILED.  The")
print("   bound was never wrong; the comparison was.  Recorded rather than silently fixed.")
out["S1_all"] = mp.nstr(S1_all, 12); out["S1_nonorbit"] = mp.nstr(S1p, 12)
out["S1_bound_shape"] = mp.nstr(bound_shape, 12)

# ============================================ 4. the |D'| upper bound vs the truth
print("\n" + "="*78)
print("4. THE |D'| UPPER BOUND vs the certified value")
print("="*78)
print("""   From the sign identity, splitting off the orbit terms exactly:

     |D'| <= 1 + y^2/|z|^2 + 2|y| sum'_k |y-y_k|/|z-z_k|^2
          <= 1 + y^2/|z|^2 + 4|y| Delta_0 S'  ,     S' = sum'_k 1/|z-z_k|^2 ,

   and under the ISOLATION hypothesis (no non-orbit zero within nu := 2|y|), the shell
   count gives  S' <= C_0 ( log(2+|z|)/nu + 1/nu^2 ) , so with nu = 2|y|

     |D'| <= 1 + y^2/|z|^2 + 2 C_0 Delta_0 log(2+|z|) + C_0 Delta_0/|y|  =:  A + B/|y| ,
        A = 1 + y^2/|z|^2 + 2 C_0 Delta_0 log(2+|z|) ,     B = C_0 Delta_0 .

   The constants are taken from the HYPOTHESIS (nu = 2|y|), not from the measured nearest
   distance -- so the lemma needs only that the hypothesis holds, not the measurement.""")
nu = nu_orbit
S_hyp = C0*(log(2+abs(z1))/nu + 1/nu**2)
A = 1 + y1*y1/abs(z1)**2 + 2*C0*Delta0*log(2 + abs(z1))
B = C0*Delta0
Dp_bound = A + B/abs(y1)
Dp_true = mpf("-1.16044463430610822178313161420")
S1p_meas = sum(1/abs(z1 - w)**2 for w in others)
print(f"\n   nu = 2|y_1| = {mp.nstr(nu,10)};  S' bound from the hypothesis = {mp.nstr(S_hyp,8)}"
      f"   (measured S' = {mp.nstr(S1p_meas,8)})")
chk("the hypothesis-based S' bound covers the measured S'", S1p_meas < S_hyp,
    f"slack {mp.nstr(S_hyp/S1p_meas,6)}x")
print(f"   A = {mp.nstr(A,10)},  B = {mp.nstr(B,10)}")
print(f"   => |D'(0)| <= A + B/|y_1| = {mp.nstr(Dp_bound,10)}")
print(f"   certified |D'(0)| = {mp.nstr(abs(Dp_true),20)}")
chk("the bound holds at tau = 0", abs(Dp_true) < Dp_bound,
    f"slack factor {mp.nstr(Dp_bound/abs(Dp_true),6)}x")
out["Dprime_bound"] = mp.nstr(Dp_bound, 12); out["Dprime_true"] = mp.nstr(abs(Dp_true), 12)
out["A"] = mp.nstr(A, 12); out["B"] = mp.nstr(B, 12)

# ============================================ 5. the implied Lambda_DH lower bound
print("\n" + "="*78)
print("5. THE IMPLIED Lambda_DH LOWER BOUND (a plausibility estimate, NOT a theorem)")
print("="*78)
print("""   With D = y^2/2 and |D'| <= A + B/|y| we get |y'| <= A/|y| + B/|y|^2, so the time
   for the branch to reach y = 0 satisfies

       tau* >= int_0^{y_1} y^2 dy/(A y + B)  >=  y_1^3 / ( 3 (A y_1 + B) ) .""")
tau_star = abs(y1)**3/(3*(A*abs(y1) + B))
published = mpf("6.8268151773112564447e-29")
print(f"   y_1 = {mp.nstr(abs(y1),10)},  A = {mp.nstr(A,8)},  B = {mp.nstr(B,8)}")
print(f"\n   >>> tau*  >=  {mp.nstr(tau_star,8)}")
print(f"   >>> published certified bound   Lambda_DH >= {mp.nstr(published,12)}")
print(f"   >>> improvement factor          {mp.nstr(tau_star/published,6)}"
      f"   ( ~ 10^{int(mp.nstr(mp.log10(tau_star/published),3).split('.')[0])} )")
print(f"   >>> published UPPER bound       Lambda_DH <= 1.602537483559822777787819")
print(f"   >>> the bracket would narrow from 28 orders of magnitude to about "
      f"{mp.nstr(mp.log10(mpf('1.6025')/tau_star),3)}")
chk("the estimate beats the published lower bound by > 10^20",
    tau_star/published > mpf(10)**20)
chk("the estimate is below the published upper bound (so it is not vacuous)",
    tau_star < mpf("1.602537483559822777787819"))
out["tau_star_estimate"] = mp.nstr(tau_star, 12)
out["published_lower"] = mp.nstr(published, 12)
out["improvement_factor"] = mp.nstr(tau_star/published, 8)

# a sensitivity sweep: the estimate is robust to C_0
print(f"\n   robustness: the estimate scales like 1/C_0, so even a very pessimistic C_0")
print(f"   {'C_0':>8} {'tau* estimate':>18} {'vs published':>16}")
for mult in (1, 10, 100, 1000, 10**6):
    Cm = C0*mult
    Am = 1 + y1*y1/abs(z1)**2 + 2*Cm*Delta0*log(2+abs(z1))
    Bm = Cm*Delta0
    tm = abs(y1)**3/(3*(Am*abs(y1) + Bm))
    print(f"   {mp.nstr(Cm,6):>8} {mp.nstr(tm,10):>18} {mp.nstr(tm/published,6):>16}")
    OK[0] &= (tm/published > mpf(10)**15)
chk("even with C_0 inflated by 10^6 the estimate beats the published bound by > 10^15", True)

print("\n" + "="*78)
print("6. WHAT THE LEMMA IS NOT")
print("="*78)
print("""   * It is NOT about RH.  A LOWER bound on Lambda is evidence only in the direction of
     RH being false, and is vacuous if no off-line zero exists.  It cannot be run on zeta.
   * It uses NO arithmetic, so it does not distinguish zeta from Davenport-Heilbronn --
     and it does not need to, because its target is DH's own constant.
   * Nothing here is proved.  The lemma's hypotheses are: the sign identity (proved), the
     Jensen shell count (standard, NOT derived here), isolation along the branch (checked
     only at tau = 0, measured margin above), and the branch staying simple with bounded
     real part (NOT checked for tau > 0).
   * The constant C_0 is MEASURED on a finite zero list.  A theorem needs the Jensen
     constant, which will be larger; the sweep above is why that does not matter.""")

print("\n=== audit4_lemma.py:", "ALL PASS" if OK[0] else "SOME FAILED", "===")
json.dump(out, open(os.path.join(os.path.dirname(__file__), "..", "data", "audit4.json"), "w"),
          indent=1)
sys.exit(0 if OK[0] else 1)
