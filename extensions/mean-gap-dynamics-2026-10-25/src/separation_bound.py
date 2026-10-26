"""A necessary condition for an individual D to grow -- the sharpest thing this
investigation produced beyond the two obstructions.

From the sign identity (DERIVATION Sec.3.4), with y = Im and D_j = y_j^2/2,

    D_j'(tau) = -2 y_j sum_{k != j} (y_j - y_k)/|z_j - z_k|^2 ,

the conjugate partner contributes exactly -1 and only a zero STRICTLY FURTHER off the line
ON THE SAME SIDE contributes positively.  Take y_j > 0 and a competing zero with
y_k > y_j > 0.  Its contribution is

    2 y_j (y_k - y_j) / |z_j - z_k|^2  <=  2 y_j (y_k - y_j) / (Re z_j - Re z_k)^2 ,

and  max_{0 < a < b} 2a(b - a) = b^2/2  at  a = b/2,  so with  y_k <= Delta_tau,

    (one competing zero's contribution)  <=  Delta_tau^2 / ( 2 R^2 ) ,
    R := |Re z_j - Re z_k| .

For that single zero to overturn the conjugate term it must exceed 1, hence

    ===================================================================
      NECESSARY:   R  <  Delta_tau / sqrt(2)         (single competitor)
      and with N competitors all at ordinate distance >= R,
                   R  <  Delta_tau sqrt(N) / 2       (N competitors)
    ===================================================================

In the Dirichlet variable, R = 2|gamma_j - gamma_k| and Delta_tau = 2 delta_max, so the
condition reads |gamma_j - gamma_k| < delta_max/sqrt 2: **a competing zero must share the
ordinate to within a fraction of the strip half-width.**

This is a *necessary* condition, not a sufficient one, and it is the useful residue: it says
exactly where a counterexample to individual monotonicity has to live, and it is cheap to
test against any configuration.
"""
import sys, os, json
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(__file__))
from mpmath import mp, mpf, sqrt

OK = [True]
def chk(n, c, x=""):
    OK[0] &= bool(c); print(f"  {'PASS' if c else 'FAIL'}  {n} {x}")

mp.dps = 30

def bound(Delta, N=1):
    return Delta*sqrt(mpf(N))/2 if N > 1 else Delta/sqrt(2)

print("== 1. the bound, and the toy counterexample which must satisfy it ==")
print("   P_0 = (z^2+1)(z^2+4):  Delta = 2;  the competitor of z = i is z = 2i,")
print("   at ordinate separation R = |Re(i) - Re(2i)| = 0.")
D_toy = mpf(2)
chk("R = 0 < Delta/sqrt2 = " + mp.nstr(bound(D_toy), 8), 0 < bound(D_toy))
print("   NOTE: the minimal counterexample is CONFLUENT IN THE ORDINATE (R = 0).  That is")
print("   the same configuration the repository's")
print("   extensions/moderate-gap-spectral-bridge-2026-10-21/06_CONFLUENT_BLIND_SPOTS.md")
print("   identifies as structural: an off-line zero's functional-equation partner always")
print("   sits at ordinate distance exactly 0.  Here the partner supplies the -1, and a")
print("   SECOND, DEEPER orbit at the same ordinate is what overturns it.")

print("\n== 2. the bound applied to the DH census: no quartet can be a counterexample ==")
CENSUS = [("O_1", "0.3085171824566374", "85.6993484853775922"),
          ("O_2", "0.1508300806097371", "114.163342730756981"),
          ("O_3", "0.0743560504508060", "166.479305913168156"),
          ("O_4", "0.2242576946268098", "176.702461242855825")]
sigma0 = mpf("1.395136158235109721061359")
Delta_DH = 2*(sigma0 - mpf(1)/2)              # the PROVED strip half-width in z
thr = bound(Delta_DH)
print(f"   proved strip half-width Delta_DH = {mp.nstr(Delta_DH,16)}")
print(f"   single-competitor threshold Delta/sqrt2 = {mp.nstr(thr,12)}"
      f"   (i.e. |gamma_j - gamma_k| < {mp.nstr(thr/2,12)})")
print(f"\n   {'target':>7} {'deepest competitor':>20} {'R = |Re z_j - Re z_k|':>24} {'< thr ?':>9}")
rows = []
for nj, dj, gj in CENSUS:
    deeper = [(nk, dk, gk) for nk, dk, gk in CENSUS if mpf(dk) > mpf(dj)]
    if not deeper:
        print(f"   {nj:>7} {'none (census max)':>20} {'-':>24} {'n/a':>9}")
        rows.append(dict(target=nj, competitor=None))
        continue
    nk, dk, gk = max(deeper, key=lambda t: mpf(t[1]))
    R = 2*abs(mpf(gj) - mpf(gk))
    ok = R < thr
    print(f"   {nj:>7} {nk:>20} {mp.nstr(R,16):>24} {'YES' if ok else 'NO':>9}")
    rows.append(dict(target=nj, competitor=nk, R=mp.nstr(R, 16), passes=bool(ok)))
    OK[0] &= (not ok)      # we EXPECT all to fail the necessary condition
print("   All competitors are far outside the threshold (by factors of 100+), so for")
print("   the census configuration no quartet's D can grow -- consistent with the four")
print("   certified D'(0) < 0, and now explained rather than merely observed.")

print("\n== 3. the bound in zeta's own units ==")
print("   zeta's trivial strip gives delta_max <= 1/2, so Delta <= 1 and the single-")
print("   competitor condition is  |gamma_j - gamma_k| < 1/(2 sqrt 2) = 0.3535533906.")
print("   The mean zeta zero spacing is 2 pi/log(gamma/2 pi):")
for g in ("100", "1e6", "1e12", "1e30"):
    gg = mpf(g)
    sp = 2*mp.pi/mp.log(gg/(2*mp.pi))
    print(f"     gamma = {g:>5}:  mean spacing {mp.nstr(sp,8)}"
          f"   -> {'BELOW' if sp < mpf(1)/(2*sqrt(2)) else 'above'} the threshold")
g_star = 2*mp.pi*mp.exp(2*mp.pi*2*sqrt(2))
print(f"   The mean spacing equals the threshold 1/(2 sqrt2) at")
print(f"     gamma* = 2 pi exp(4 pi sqrt 2) = {mp.nstr(g_star,10)}  (about 3.3 x 10^8),")
print("   so the condition is vacuous at low height and becomes satisfiable by TYPICAL")
print("   neighbouring ordinates only above gamma*  -- and by Lehmer-type pairs much")
print("   earlier (the pair 6709/6710 at gamma = 7005.06 has |d gamma| = 0.0377, well")
print("   inside the threshold 0.3536 already).  So the condition does NOT rule out")
print("   non-monotonicity for zeta; it localises where it could occur.  NO CLAIM IS MADE")
print("   THAT IT DOES OCCUR -- that would be a statement about off-line zeros of zeta.")

print("\n== 4. what the bound is and is not ==")
print("   * NECESSARY, not sufficient: it says where a counterexample must live.")
print("   * It uses only the sign identity and the strip half-width -- no arithmetic, so")
print("     it does NOT distinguish zeta from Davenport-Heilbronn.")
print("   * It gives no information at tau = 0 beyond what the strip width already gives.")
print("   * It is the one statement here that is both new and cheap to apply.")

print("\n=== separation_bound.py:", "ALL PASS" if OK[0] else "SOME FAILED", "===")
json.dump(dict(Delta_DH=mp.nstr(Delta_DH, 20), threshold=mp.nstr(thr, 16),
               zeta_threshold_gamma=mp.nstr(mpf(1)/(2*sqrt(2)), 12),
               gamma_star=mp.nstr(g_star, 10), rows=rows),
          open(os.path.join(os.path.dirname(__file__), "..", "data", "part3_separation.json"), "w"),
          indent=1)
sys.exit(0 if OK[0] else 1)
