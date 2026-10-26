"""The zero-interaction representation of D', and its REMAINDER, on the DH control.

The canonical product is set up so that no normalisation is left free:

  H_tau is EVEN in z (Phi real and even) and real on the real axis, so H_tau(z) = h_tau(z^2)
  with h_tau entire.  |H_tau(z)| <= int e^{tau u^2}|Phi(u)| cosh(u |Im z|) du grows like
  exp(C|z| log|z|) because Phi decays double-exponentially, so H_tau has order 1 and
  h_tau has order 1/2.  Order < 1 forces GENUS 0: h_tau(w) = h_tau(0) prod_k (1 - w/w_k),
  absolutely convergent, with NO exponential factor and no normalisation ambiguity.
  (This is the step that would be illegitimate in the z variable, where sum 1/|z_k|
  diverges and genus 1 would force an e^{cz} factor.)

Consequently, with w_k = a_k^2,

  H'/H(z) = 2z sum_k 1/(z^2 - a_k^2),       absolutely convergent off the zeros,
  H''/H'(a_j) = 2 [ 1/(2 a_j) + sum_{k!=j} 2 a_j/(a_j^2 - a_k^2) ]
             = 2 sum_{k != j} 1/(a_j - z_k)   (full symmetric zero list, paired order)

and therefore the SIGN IDENTITY

  D'(tau) = -2 y_j sum_{k != j} (y_j - y_k) / |z_j - z_k|^2 ,     y = Im .

The part of the sum over zeros with y_k = 0 converges ABSOLUTELY (terms y_j^2/|z_j-z_k|^2
with sum 1/|z_j-z_k|^2 < infinity by the order-1 density), so the identity is not a formal
rearrangement.  This script evaluates the partial sums over an explicitly computed finite
zero list and reports the REMAINDER against the exact H''/H' value -- it does not replace
the infinite sum by a finite one.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
from mpmath import mp, mpf, mpc, findroot, log, pi
import kernels as K
from dh_flow import zdot_R1, CENSUS

mp.dps = 30
kap = K.kappa_closed()

def Zline(t):
    return mp.re(K.Lambda_DH(mpf(1)/2 + mpc(0, 1)*t, kap))

def online_zeros(lo, hi, step=mpf("0.05")):
    out = []
    t = mpf(lo); prev = Zline(t)
    while t < hi:
        t2 = t + step
        cur = Zline(t2)
        if (prev < 0) != (cur < 0):
            out.append(findroot(Zline, (t, t2), solver='bisect', tol=mpf(10)**-35))
        t, prev = t2, cur
    return out

def Ncount(T):
    """main term of the one-sided count, (T/2pi) log(5T/(2 pi e))  [parent R1]."""
    return T/(2*pi)*log(5*T/(2*pi*mp.e))

def main():
    HI = mpf(os.environ.get("DH_SCAN_HI", "260"))
    print(f"== on-line zero scan on 0 < Im s < {HI} (step 0.05, bisection to 1e-35) ==")
    zs = online_zeros(mpf(2), HI)
    print(f"   on-line ordinates found: {len(zs)}")
    offl = [(mpf(b) - mpf(1)/2, mpf(g)) for _, b, g in CENSUS]
    print(f"   off-line quartets from the parent census (Im s < 200): {len(offl)}")
    for T in (50, 100, 150, 200, int(HI)):
        n_on = sum(1 for t in zs if t < T)
        n_off = sum(1 for d, g in offl if g < T) * 2       # two zeros per quartet in the upper half
        print(f"   T={T:>4}: on-line {n_on:>3} + 2*off-line {n_off:>3} = {n_on+n_off:>3}"
              f"   vs main term {mp.nstr(Ncount(mpf(T)),8):>10}")

    # the target: the exact cross-orbit part of D'(O_1)
    b1, g1 = mpf(CENSUS[0][1]), mpf(CENSUS[0][2])
    d1 = b1 - mpf(1)/2
    z1 = 2*g1 - 2*mpc(0, 1)*d1
    y1, x1 = mp.im(z1), mp.re(z1)
    mp.dps = 60
    Dp_exact = y1 * mp.im(zdot_R1(z1, K.kappa_closed()))
    self_term = -1 - y1*y1/(x1*x1 + y1*y1)
    cross_exact = Dp_exact - self_term
    mp.dps = 30
    print(f"\n   exact      D'(0)      = {mp.nstr(Dp_exact, 20)}   (from H''/H', no truncation)")
    print(f"   own-orbit  self term  = {mp.nstr(self_term, 20)}   = -1 - y^2/|z|^2, EXACT")
    print(f"   cross-orbit remainder = {mp.nstr(cross_exact, 20)}")

    # build the full zero list (upper and lower half planes, and mirrors)
    def term(zk):
        return (y1 - mp.im(zk)) / abs(z1 - zk)**2

    print(f"\n== partial sums of the interaction identity, and the REMAINDER ==")
    print(f"   {'|Im s| <':>10} {'#zeros used':>12} {'partial cross':>20} {'remainder':>16} {'rem/|D'+chr(39)+'|':>12}")
    rows = []
    for T in (25, 50, 75, 100, 125, 150, 175, 200, 225, int(HI)):
        T = mpf(T)
        if T > HI: continue
        acc = mpf(0)
        used = 0
        for t in zs:
            if t > T: continue
            for zk in (2*t, -2*t):              # z = 2 gamma, real; both halves
                acc += term(zk); used += 1
        for d, g in offl:
            if g > T: continue
            for zk in (2*g - 2*mpc(0,1)*d, 2*g + 2*mpc(0,1)*d,
                       -2*g + 2*mpc(0,1)*d, -2*g - 2*mpc(0,1)*d):
                if abs(zk - z1) < mpf(10)**-20: continue        # skip z1 itself
                if abs(zk - mp.conj(z1)) < mpf(10)**-20: continue   # in the self term
                if abs(zk + z1) < mpf(10)**-20: continue            # in the self term
                if abs(zk + mp.conj(z1)) < mpf(10)**-20: continue   # in the self term
                acc += term(zk); used += 1
        partial = -2*y1*acc
        rem = cross_exact - partial
        rows.append(dict(T=str(T), nzeros=used, partial=mp.nstr(partial, 16),
                         remainder=mp.nstr(rem, 10),
                         rel=mp.nstr(abs(rem/Dp_exact), 6)))
        print(f"   {mp.nstr(T,6):>10} {used:>12} {mp.nstr(partial,14):>20} "
              f"{mp.nstr(rem,8):>16} {mp.nstr(abs(rem/Dp_exact),4):>12}")

    print("\n   The remainder SHRINKS monotonically with the cutoff, which is what an")
    print("   absolutely convergent tail must do; it is reported as a measured residual,")
    print("   NOT as a bound.  No conclusion below rests on truncating the sum: every")
    print("   reported D' comes from H''/H', where the whole spectrum is present.")

    # the zero gap at gamma_1 -- the reason the real-zero drag is weak
    near = sorted(zs, key=lambda t: abs(t - g1))[:2]
    gap = abs(near[0] - near[1])
    dens = log(5*g1/(2*pi))/(2*pi)                 # dN/d gamma
    hz = 2/dens                                     # mean spacing in z = 2 gamma
    # a uniform lattice of spacing hz, averaged over the offset a, gives exactly
    #   < sum_n 1/((n hz + a)^2 + y^2) >_a = pi/(hz |y|)
    S_pred = pi/(hz*abs(y1))
    pred = -2*y1*y1*S_pred
    print(f"\n== why the on-line drag is weak: a zero GAP at gamma_1 ==")
    print(f"   nearest on-line ordinates to gamma_1 = {mp.nstr(g1,10)}: "
          f"{mp.nstr(min(near),12)} and {mp.nstr(max(near),12)}")
    print(f"   gap = {mp.nstr(gap,8)}  vs mean spacing 1/density = {mp.nstr(1/dens,8)}"
          f"   ->  {mp.nstr(gap*dens,5)} x the mean spacing")
    print(f"\n   MEAN-SPACING HEURISTIC, computed rather than guessed.  For a uniform lattice")
    print(f"   of spacing h in z, averaged over the offset, sum 1/|z1-zk|^2 = pi/(h|y|):")
    print(f"     density dN/dgamma = {mp.nstr(dens,8)},  mean spacing in z  h = {mp.nstr(hz,8)}")
    print(f"     predicted  sum 1/|z1-zk|^2 = {mp.nstr(S_pred,8)}")
    print(f"     predicted  contribution to D' = -2 y^2 S = {mp.nstr(pred,8)}")
    print(f"     ACTUAL cross-orbit contribution        = {mp.nstr(cross_exact,8)}")
    print(f"     the heuristic OVERESTIMATES the drag by a factor "
          f"{mp.nstr(abs(pred/cross_exact),6)}")
    print(f"   Note the size: the prediction {mp.nstr(pred,6)} is LARGER in magnitude than the")
    print(f"   exact conjugate-partner term -1, so if it were valid it would roughly double")
    print(f"   the decay rate.  It is not valid, because the quartet sits in a "
          f"{mp.nstr(gap*dens,4)}x zero gap.")
    print(f"   The mean-spacing heuristic is therefore NOT usable as the dissipation lower")
    print(f"   bound that the (S2) route needs.")

    with open(os.path.join(os.path.dirname(__file__), "..", "data", "dh_online_ordinates.txt"), "w") as f:
        for t in zs:
            f.write(mp.nstr(t, 25) + "\n")

    # GATES -- assert the inequalities, do not print a stored label
    ok = True
    print("\n== GATES ==")
    # (i) the count reconciliation of the parent extension, reproduced
    n200 = sum(1 for t in zs if t < 200)
    g1_ = (n200 == 122) and (n200 + 2*4 == 130)
    ok &= g1_
    print(f"   {'PASS' if g1_ else 'FAIL'}  census reconciliation: {n200} on-line below T=200,"
          f" {n200} + 2*4 = {n200+8} (parent contour count 130)")
    # (ii) the remainder decreases monotonically in the cutoff, and is small at the end
    rems = [abs(mpf(r["remainder"])) for r in rows]
    g2_ = all(rems[i] >= rems[i+1] for i in range(len(rems)-1))
    g3_ = rems[-1] < abs(cross_exact)/50
    ok &= g2_ and g3_
    print(f"   {'PASS' if g2_ else 'FAIL'}  the remainder is monotonically decreasing in the cutoff")
    print(f"   {'PASS' if g3_ else 'FAIL'}  the final remainder is below 2% of the cross term"
          f"  ({mp.nstr(rems[-1]/abs(cross_exact)*100, 3)}%)")
    # (iii) the exact decomposition reconciles
    g4_ = abs((self_term + cross_exact) - Dp_exact) < mpf(10)**-25
    ok &= g4_
    print(f"   {'PASS' if g4_ else 'FAIL'}  self + cross = D' exactly"
          f"  (residual {mp.nstr(abs((self_term+cross_exact)-Dp_exact), 3)})")
    # (iv) the mean-spacing heuristic really does overestimate, by a lot
    g5_ = abs(pred/cross_exact) > 4
    ok &= g5_
    print(f"   {'PASS' if g5_ else 'FAIL'}  the mean-spacing heuristic overestimates by > 4x"
          f"  ({mp.nstr(abs(pred/cross_exact), 6)}x)")
    # (v) the zero gap at gamma_1 really is anomalous
    g6_ = gap*dens > 2
    ok &= g6_
    print(f"   {'PASS' if g6_ else 'FAIL'}  the gap at gamma_1 exceeds 2 mean spacings"
          f"  ({mp.nstr(gap*dens, 4)})")
    print("\n=== dh_interaction.py:", "ALL PASS" if ok else "SOME FAILED", "===")

    json.dump(dict(n_online=len(zs), scan_hi=str(HI), Dprime=mp.nstr(Dp_exact, 20),
                   self_term=mp.nstr(self_term, 20), cross=mp.nstr(cross_exact, 20),
                   partials=rows, gap=mp.nstr(gap, 10),
                   gap_in_mean_spacings=mp.nstr(gap*dens, 8),
                   mean_spacing_z=mp.nstr(hz, 10),
                   heuristic_S=mp.nstr(S_pred, 10),
                   heuristic_Dprime_contrib=mp.nstr(pred, 10),
                   heuristic_overestimate_factor=mp.nstr(abs(pred/cross_exact), 8),
                   gates_pass=bool(ok)),
              open(os.path.join(os.path.dirname(__file__), "..", "data", "part3_interaction.json"), "w"),
              indent=1)

    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
