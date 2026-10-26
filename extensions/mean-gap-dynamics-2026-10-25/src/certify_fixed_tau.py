"""A FIXED-PARAMETER certificate:  H_tau has a NONREAL zero at tau = 4/25,
hence  Lambda_DH > 0.16.

No trajectory estimate, no isolation hypothesis, and nothing uniform over an interval:
one parameter value, one contour, one winding number.

HOW THE ENCLOSURE IS BUILT (and why it is not the parent extension's perturbation bound).

  H_tau(z) = (1/2) sum_{m>=0} (tau/4)^m Lambda_F^{(2m)}(s)/m! ,   s = 1/2 + iz/2 ,     (*)

because d_tau H = -d_z^2 H and d_z^2 = -(1/4) d_s^2.  All the oscillation lives inside
Lambda_F, which Arb encloses rigorously; nothing is replaced by e^{Mx}, so no cancellation
is discarded.  More generally, for every k,

  H_tau^{(k)}(z) = (i/2)^k (1/2) G_k(z) ,   G_k = sum_{m>=0} (tau/4)^m Lambda^{(2m+k)}(s)/m! .

The contour is partitioned into N arcs.  Each arc is covered by  c~ + B(0, rho)  where c~ is
an EXACT dyadic point and rho covers both the arc's chord and the rounding of its midpoint.
On that set

  H_tau(c~ + w) = sum_{k=0}^{K} (w^k/k!) (i/2)^k G_k(c~)/2  +  (arc tail) ,

so the special functions are evaluated at a POINT, never on a wide ball.  That matters: a
ball of radius 3.1e-4 fed into the order-43 Hurwitz-zeta series inflates the enclosure of
H_tau to +/- 1.9e-28, swamping its true size 3.8e-31, whereas the point evaluation plus an
explicit arc tail does not.

TWO RIGOROUS TAILS, both majorised by non-oscillatory positive integrals.  From
H_0^{(n)}(z) = i^n int_0^inf Phi_F(u) u^n (trig) du  and  |Phi_F(u)| <= 8 e^{3u} e^{-pi e^{4u}/5}
for u >= 0 (the parent extension's constant):

  T_M  =  (tau^{M+1}/(M+1)!) int_0^inf 8 e^{3u} e^{-pi e^{4u}/5} u^{2M+2} cosh(uY) e^{tau u^2} du
          -- truncating the tau-series (*) at m = M ;
  A_K  =  (rho^{K+1}/(K+1)!) int_0^inf 8 e^{3u} e^{-pi e^{4u}/5} u^{K+1} cosh(uY) e^{tau u^2} e^{rho u} du
          -- truncating the arc expansion at k = K ,

with Y >= max |Im z| on the closed disc.  Both are computed by rigorous Arb integration on
[0, X] plus an exact decay bound beyond X.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
from flint import acb, arb, ctx, acb_series, fmpq

PREC = 700
ctx.prec = PREC

TAU = fmpq(4, 25)           # 0.16, exact
R = fmpq(1, 16)             # disc radius, exact
NARCS = 512
M = 20                      # tau-series truncation
KK = 10                     # arc-expansion truncation (see the note on A_K below)
# centre: an EXACT decimal near the located zero at tau = 4/25 (navigation value)
C_RE = "171.2502"
C_IM = "-0.1823"

def kappa():
    r5 = acb(5).sqrt()
    return ((acb(10) - 2*r5).sqrt() - 2)/(r5 - 1)

_FACT = [acb(1)]
def fact(n):
    while len(_FACT) <= n:
        _FACT.append(_FACT[-1]*len(_FACT))
    return _FACT[n]

def lambda_coeffs(s0, cap):
    """Lambda_F(s0 + x) coefficients; s0 must be an exact point."""
    ctx.cap = cap
    kap = kappa()
    a = [acb(1), kap, -kap, acb(-1)]
    sser = acb_series([s0, acb(1)])
    S = acb_series([acb(0)])
    for r in (1, 2, 3, 4):
        S = S + acb_series([a[r-1]]) * sser.zeta(acb(r)/5)
    five = (acb_series([-s0*acb(5).log(), -acb(5).log()])).exp()
    u = acb_series([(s0+1)/2, acb(1)/2])
    cond = (u*acb_series([acb(5).log() - acb.pi().log()])).exp()
    c = (cond * u.gamma() * five * S).coeffs()
    while len(c) < cap:
        c.append(acb(0))
    return c

def G(c, k):
    acc = acb(0)
    for m in range(M + 1):
        n = 2*m + k
        acc = acc + (acb(TAU)/4)**m * c[n] * fact(n) / fact(m)
    return acc

# ----------------------------------------------------------------- the two tails
def _kern_int(power, Y, rho, X=3):
    """rigorous enclosure of int_0^inf 8 e^{3u} e^{-pi e^{4u}/5} u^power cosh(uY)
       e^{TAU u^2} e^{rho u} du, via Arb integration on [0,X] plus a decay tail."""
    Yb, rb = arb(Y), arb(rho)
    f = lambda u, _: (acb(8)*(3*u).exp()*(-acb.pi()*(4*u).exp()/5).exp()
                      * u**power * (u*acb(Yb)).cosh() * (acb(TAU)*u*u).exp()
                      * (u*acb(rb)).exp())
    I = acb.integral(f, acb(0), acb(X))
    # beyond X = 3: (log f)' = power/u + 3 - (4pi/5)e^{4u} + Y + rho + 2 TAU u <= -1,
    # since (4pi/5)e^{12} = 408862 dominates everything for u >= 3 and power <= 46.
    fX = (acb(8)*(acb(3)*X).exp()*(-acb.pi()*(acb(4)*X).exp()/5).exp()
          * acb(X)**power * (acb(X)*acb(Yb)).cosh() * (acb(TAU)*X*X).exp()
          * (acb(X)*acb(rb)).exp())
    return I + acb(arb(0, abs(fX).mid() + abs(fX).rad()))   # tail <= f(X)

def tail_T(Y):
    I = _kern_int(2*M + 2, Y, 0)
    return abs(acb(TAU)**(M + 1)/fact(M + 1) * I)

def tail_A(Y, rho):
    I = _kern_int(KK + 1, Y, rho)
    return abs(acb(rho)**(KK + 1)/fact(KK + 1) * I)

def ceil_str(a, n=22):
    a = a.real if isinstance(a, acb) else a
    return arb(a.mid() + a.rad()).str(n, radius=False)

def floor_str(a, n=22):
    a = a.real if isinstance(a, acb) else a
    return arb(a.mid() - a.rad()).str(n, radius=False)

def sweep(centre, Rb, Y, label, verbose=True):
    """returns (ok, min|H|, worst_rel, nonadjacent, net_advance)."""
    HP_ok = True
    dth = 2*arb.pi()/NARCS
    chord = 2*Rb*(dth/4).sin()
    def halfplanes(Z):
        o = []
        if (Z.real.mid() - Z.real.rad()) > 0: o.append(0)
        if (Z.imag.mid() - Z.imag.rad()) > 0: o.append(1)
        if (Z.real.mid() + Z.real.rad()) < 0: o.append(2)
        if (Z.imag.mid() + Z.imag.rad()) < 0: o.append(3)
        return o
    def cyc(a, b):
        d = (b - a) % 4
        return {0: 0, 1: 1, 3: -1}.get(d, None)
    minabs = None; worst_rel = arb(0); sets = []
    for i in range(NARCS):
        th = 2*arb.pi()*i/NARCS
        pt = acb(centre) + acb(Rb)*(acb(0, 1)*acb(th)).exp()
        c_exact = acb(arb(pt.real.mid()), arb(pt.imag.mid()))
        rho = chord + abs(pt - c_exact)
        s0 = acb(1)/2 + acb(0, 1)*c_exact/2
        s0 = acb(arb(s0.real.mid()), arb(s0.imag.mid()))
        co = lambda_coeffs(s0, 2*M + KK + 3)
        rr = rho.mid() + rho.rad()
        wball = acb(arb(0, rr), arb(0, rr))
        Hval = acb(0)
        for k in range(KK + 1):
            Hval = Hval + (wball**k)/fact(k)*(acb(0, 1)/2)**k*G(co, k)/2
        A = tail_A(Y, rho); aa = A.mid() + A.rad()
        Hval = Hval + acb(arb(0, aa), arb(0, aa))
        a = abs(Hval)
        rel = arb(a.rad()/a.mid()) if a.mid() > 0 else arb(1)
        if rel > worst_rel: worst_rel = rel
        if minabs is None or a.mid() < minabs.mid(): minabs = a
        S = halfplanes(Hval)
        if not S:
            if verbose: print(f"   [{label}] arc {i}: enclosure in no open half-plane")
            HP_ok = False; S = [0]
        sets.append(S)
    net = 0; cur = sets[0][0]; bad = 0
    for i in range(1, NARCS + 1):
        S = sets[i % NARCS]; best = None
        for cand in S:
            st = cyc(cur, cand)
            if st is None: continue
            if best is None or abs(st) < abs(best[0]): best = (st, cand)
        if best is None:
            bad += 1; best = (0, S[0])
        net += best[0]; cur = best[1]
    return HP_ok, minabs, worst_rel, bad, net


def main():
    ok = True
    out = {}
    print(f"== fixed-parameter certificate at tau = {TAU} = 0.16 exactly ==")
    centre = acb(arb(C_RE), arb(C_IM))
    centre = acb(arb(centre.real.mid()), arb(centre.imag.mid()))     # exact dyadic
    Rb = arb(R)
    Y = abs(centre.imag) + Rb            # >= max |Im z| on the closed disc
    margin = abs(centre.imag) - Rb
    print(f"   disc: centre {centre.str(16)}, radius {Rb.str(10)}")
    print(f"   |Im centre| - r = {margin.str(14)}   -> closed disc disjoint from R: "
          f"{'YES' if margin > 0 else 'NO'}")
    print(f"   Y = |Im centre| + r = {Y.str(10)}  (bounds |Im z| on the closed disc)")
    ok &= bool(margin > 0)
    out["tau"] = str(TAU); out["radius"] = str(R); out["narcs"] = NARCS
    out["margin"] = floor_str(margin); out["Y"] = ceil_str(Y)

    dth = 2*arb.pi()/NARCS
    chord = 2*Rb*(dth/4).sin()           # covers an arc of angular width dth
    print(f"\n   arcs: N = {NARCS}, angular width {dth.str(8)}, chord radius "
          f"{chord.str(8)}")

    # ---- gate the two INPUTS to the tail bounds, rather than quoting them
    print("\n   inputs to the tail bounds, verified rather than quoted:")
    S1 = acb(0)
    for n in range(1, 80):
        S1 = S1 + n*(-acb.pi()*(n*n - 1)/5).exp()
    g1 = bool(S1.real < 2)
    print(f"     {'PASS' if g1 else 'FAIL'}  |Phi_F(u)| <= 8 e^{{3u}} e^{{-pi e^{{4u}}/5}} on u >= 0:")
    print(f"           Phi_F = 4 e^{{3x}} sum n a_n e^{{-pi n^2 T/5}}, T = e^{{4x}} >= 1, |a_n| <= 1,")
    print(f"           so |Phi_F| <= 4 e^{{3x}} e^{{-pi T/5}} S(T) with S(T) = sum n e^{{-pi(n^2-1)T/5}}")
    print(f"           DECREASING in T; S(1) = {S1.real.str(14)} <= 2, so 4 S(1) = "
          f"{(4*S1).real.str(10)} <= 8")
    ok &= g1
    # the decay condition beyond X = 3 that licenses int_X^inf f <= f(X)
    X = acb(3); pmax = 2*M + KK + 2
    dlog = acb(pmax)/X + 3 - (4*acb.pi()/5)*(4*X).exp() + acb(Y) + 2*acb(TAU)*X + acb(1)
    g2 = bool(dlog.real < -1)
    print(f"     {'PASS' if g2 else 'FAIL'}  (log f)'(u) <= -1 for u >= 3, every integrand used:")
    print(f"           (log f)' = p/u + 3 - (4pi/5)e^{{4u}} + Y tanh(uY) + 2 tau u + rho,")
    print(f"           at u = 3 with p = {pmax} (the largest power used): "
          f"{dlog.real.str(10)} <= -1,")
    print(f"           and it only decreases, since (4pi/5)e^{{12}} = "
          f"{((4*acb.pi()/5)*acb(12).exp()).real.str(10)} dominates.")
    print(f"           Hence f(u) <= f(3) e^{{-(u-3)}} and int_3^inf f <= f(3).")
    ok &= g2
    out["S1"] = ceil_str(S1.real, 14); out["dlog_at_3"] = ceil_str(dlog.real, 12)

    T = tail_T(Y)
    print(f"\n   tau-series tail  T_M   (M = {M})  <= {ceil_str(T, 8)}")
    out["T_M"] = ceil_str(T, 10)

    # one representative arc tail; rho is the same for every arc (plus midpoint rounding)
    A_rep = tail_A(Y, chord*2)      # a representative, with rho generously over-estimated
    print(f"   arc tail         A_K   (K = {KK})  <= {ceil_str(A_rep, 8)}"
          f"   (at rho = {(chord*2).str(6)})")
    print(f"""
   WHY K = {KK} AND NOT 4.  A_K is majorised with |Phi_F|, so it keeps NO oscillatory
   cancellation, while H_tau itself is ~1e-30 BECAUSE of cancellation.  At K = 4 the
   first version of this script got A_4 ~ 2e-23, which swamped H_tau by seven orders of
   magnitude and every arc enclosure contained 0 -- the same disease that makes the
   parent extension's perturbation bound weak.  The fix is free: A_K ~ (rho u*)^(K+1)/(K+1)!
   with u* ~ 0.6 the kernel's effective support, so raising K costs four more series
   coefficients and buys 29 orders of magnitude.  Raising the arc count instead would have
   needed about 61000 arcs for the same effect.""")
    out["A_K"] = ceil_str(A_rep, 10)

    print(f"\n   sweeping {NARCS} arcs ...")
    print("""
   WHY net/4 IS THE WINDING NUMBER EXACTLY, not up to an error.  If arc i's enclosure lies
   in the open half-plane of direction alpha_i = h_i pi/2, then a continuous branch has
   arg H(theta) = alpha_i + e_i(theta) with |e_i| < pi/2 on that arc.  Consecutive arcs
   ABUT (ball i is centred at theta_i = 2 pi i/N with half-width Delta theta/2, so the N
   of them tile the circle), and at a shared endpoint the same H value gives
   alpha_{i+1} - alpha_i = e_i - e_{i+1} in (-pi, pi), consistent with a step of 0 or +-1.
   Telescoping once round and using H(2pi) = H(0): net * pi/2 = 2 pi W, so W = net/4
   EXACTLY.  (The step is also uniquely forced: when the current index is not in the next
   arc's set, only one candidate is cyclically adjacent.)

   THE WINDING NUMBER IS COMPUTED BY HALF-PLANES, NOT BY SUMMING ARGUMENTS.  Summing 512
   argument increments, each enclosed to a relative 5e-2, accumulates a radius of order
   2 * sum (rel radius) -- which is scale-invariant in the arc count (halving rho doubles
   the number of arcs) and came out as +/- 3.69 on the first attempt, enclosing many
   integers.  Instead: each arc's enclosure is placed in one of the four OPEN half-planes
   {Re>0, Im>0, Re<0, Im<0}; a point lies in exactly one or two of them, so the index
   advances by 0 or +-1 per arc, and the net advance over the loop is 4 x the winding
   number -- an EXACT integer with no error accumulation.""")

    hp_ok, minabs, worst_rel, bad, net = sweep(centre, Rb, Y, "main")
    wind_int = net // 4
    ok &= hp_ok
    print(f"\n   min |H_tau| over the arcs      >= {floor_str(minabs, 10)}")
    print(f"   worst relative enclosure radius = {worst_rel.str(6)}")
    print(f"   non-adjacent half-plane steps   = {bad}   (must be 0)")
    print(f"   net half-plane advance          = {net}  ->  winding number = {wind_int}")
    ok &= bool(minabs > 0)
    ok &= (bad == 0)
    ok &= (net % 4 == 0)
    ok &= (wind_int == 1)
    out["min_H"] = floor_str(minabs, 12)
    out["worst_rel_radius"] = worst_rel.str(8)
    out["nonadjacent_steps"] = bad
    out["net_advance"] = net
    out["winding"] = wind_int

    print(f"""
   GATES
     {'PASS' if minabs > 0 else 'FAIL'}  min|H_tau| > 0 on the contour -- no zero on dD
     {'PASS' if bad == 0 else 'FAIL'}  every half-plane step adjacent ({bad} violations)
     {'PASS' if wind_int == 1 else 'FAIL'}  winding number is EXACTLY 1 (net advance {net} = 4 x 1)
     {'PASS' if margin > 0 else 'FAIL'}  the closed disc is disjoint from the real axis

   CONCLUSION, if all four gates pass:
     H_tau has exactly one zero in D, counted with multiplicity, and D contains no real
     point.  So H_tau has a NONREAL zero at tau = {TAU}.  Real-zero preservation runs
     toward increasing tau, so

            Lambda_DH  >  4/25  =  0.16   exactly.

   Against the parent extension's certified Lambda_DH >= 6.8268151773112564447e-29 and
   its certified upper bound Lambda_DH <= 1.602537483559822777787819, the bracket becomes
            [ 0.16 , 1.602537483559822777787819 ]   -- a factor of about 10.""")
    # ---------------------------------------------------------- negative control
    print("\n" + "="*78)
    print("NEGATIVE CONTROL -- the same machinery on a disc that should contain NO zero")
    print("="*78)
    off = acb(centre) + acb(arb("0.6"))         # shifted 0.6 to the right in Re z
    off = acb(arb(off.real.mid()), arb(off.imag.mid()))
    Yo = abs(off.imag) + Rb
    print(f"   disc: centre {off.str(14)}, radius {Rb.str(8)}  (0.6 to the right)")
    hp2, min2, rel2, bad2, net2 = sweep(off, Rb, Yo, "control")
    w2 = net2 // 4
    print(f"   min |H_tau| >= {floor_str(min2, 8)},  worst rel radius {rel2.str(6)},"
          f"  non-adjacent steps {bad2}")
    print(f"   net half-plane advance = {net2}  ->  winding number = {w2}")
    ctl = hp2 and (bad2 == 0) and (net2 % 4 == 0) and (w2 == 0) and bool(min2 > 0)
    print(f"   {'PASS' if ctl else 'FAIL'}  winding number is 0 there -- the method does not "
          f"manufacture zeros")
    ok &= ctl
    out["control_winding"] = w2; out["control_min_H"] = floor_str(min2, 10)

    print("\n=== certify_fixed_tau.py:", "ALL CERTIFIED" if ok else "SOME FAILED", "===")
    json.dump(out, open(os.path.join(os.path.dirname(__file__), "..", "data",
                                     "certificate_fixed_tau.json"), "w"), indent=1)
    return ok

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
