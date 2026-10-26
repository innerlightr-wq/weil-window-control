"""INDEPENDENT cross-check of the representation the fixed-parameter certificate rests on.

The certificate evaluates

    H_tau(z) = (1/2) sum_{m>=0} (tau/4)^m Lambda_F^{(2m)}(s)/m! ,     s = 1/2 + i z/2 ,   (*)

derived from d_tau H = -d_z^2 H together with d_z^2 = -(1/4) d_s^2.  A wrong constant in (*)
-- (tau/4)^m for (tau/2)^m, a sign, a factor of 2 -- would NOT be caught by:

  * the internal M-consistency of the series (M = 12 vs M = 20 agree for any constant);
  * the stage-0 dictionary audit, which checks d_tau H = -d_z^2 H through the QUADRATURE
    route only, and H_0 = Lambda_F/2 at tau = 0 only;
  * the winding-number gates, which are self-consistent whatever H is.

So (*) is checked here at tau != 0 against the DEFINING integral

    H_tau(z) = int_0^infty e^{tau u^2} Phi_F(u) cos(z u) du ,

evaluated by oscillation-resolving mpmath quadrature at high precision -- a route that
shares no code with the series and calls no special function beyond exp and cos.  Note that
|H_tau| ~ 1e-30 while the integrand is O(1), so the quadrature must resolve ~30 digits of
cancellation; the working precision is set accordingly and a precision-doubling check is
reported.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
from flint import acb, arb, ctx, acb_series
from mpmath import mp, mpf, mpc
import kernels as K

PREC = 500
M = 20

def kappa():
    r5 = acb(5).sqrt()
    return ((acb(10) - 2*r5).sqrt() - 2)/(r5 - 1)

_F = [acb(1)]
def fact(n):
    while len(_F) <= n:
        _F.append(_F[-1]*len(_F))
    return _F[n]

def H_series(z_re, z_im, tau, M=M):
    """(*) evaluated with Arb; z given as exact decimal strings."""
    ctx.prec = PREC
    z = acb(arb(z_re), arb(z_im))
    s0 = acb(1)/2 + acb(0, 1)*z/2
    s0 = acb(arb(s0.real.mid()), arb(s0.imag.mid()))
    ctx.cap = 2*M + 2
    kap = kappa()
    a = [acb(1), kap, -kap, acb(-1)]
    sser = acb_series([s0, acb(1)])
    S = acb_series([acb(0)])
    for r in (1, 2, 3, 4):
        S = S + acb_series([a[r-1]]) * sser.zeta(acb(r)/5)
    five = (acb_series([-s0*acb(5).log(), -acb(5).log()])).exp()
    u = acb_series([(s0+1)/2, acb(1)/2])
    cond = (u*acb_series([acb(5).log() - acb.pi().log()])).exp()
    c = (cond*u.gamma()*five*S).coeffs()
    while len(c) < 2*M + 2:
        c.append(acb(0))
    acc = acb(0)
    for m in range(M + 1):
        acc = acc + (acb(tau)/4)**m * c[2*m] * fact(2*m) / fact(m)
    return acc/2

def H_quad(z_re, z_im, tau, dps):
    """the defining integral, oscillation-resolving nodes, no special functions."""
    mp.dps = dps
    z = mpc(mpf(z_re), mpf(z_im))
    return K.H_tau(z, mpf(tau), "dh")

def main():
    ok = True
    rows = []
    PTS = [
        ("171.2502", "-0.1823", "4/25", "0.16"),      # the certificate's disc centre
        ("171.2502", "-0.1823", "1/10", "0.1"),
        ("171.3100", "-0.2500", "4/25", "0.16"),      # another point on/near the disc
        ("20.0",     "-0.4000", "4/25", "0.16"),      # a low-height point, |H| is O(1)-ish
        ("5.0",      "0.3000",  "1/5",  "0.2"),       # small z, large tau
    ]
    print("== (*) vs the defining integral, at tau != 0 ==")
    print(f"   Arb series M = {M}; quadrature at dps 90 and 140 (precision doubling)\n")
    for zr, zi, tfrac, tdec in PTS:
        hs = H_series(zr, zi, tdec)
        a = H_quad(zr, zi, tdec, 90)
        b = H_quad(zr, zi, tdec, 140)
        mp.dps = 60
        sm = mpc(mpf(arb(hs.real.mid()).str(50, radius=False)),
                 mpf(arb(hs.imag.mid()).str(50, radius=False)))
        dbl = abs(a - b)/max(abs(b), mpf(10)**-120)
        rel = abs(sm - b)/max(abs(b), mpf(10)**-120)
        good = (rel < mpf(10)**-25) and (dbl < mpf(10)**-25)
        ok &= good
        print(f"   z = {zr} {zi} i,  tau = {tfrac}")
        print(f"     series (*)  {mp.nstr(sm, 20)}")
        print(f"     integral    {mp.nstr(b, 20)}")
        print(f"     |H|         {mp.nstr(abs(b), 8)}")
        print(f"     relative deviation series vs integral : {mp.nstr(rel, 4)}")
        print(f"     quadrature precision doubling 90->140 : {mp.nstr(dbl, 4)}"
              f"    {'PASS' if good else 'FAIL'}")
        rows.append(dict(z_re=zr, z_im=zi, tau=tfrac, absH=mp.nstr(abs(b), 10),
                         rel_series_vs_integral=mp.nstr(rel, 6),
                         rel_precision_doubling=mp.nstr(dbl, 6), pass_=bool(good)))

    # a WRONG-CONSTANT control: the same series with (tau/2)^m must DISAGREE
    print("\n== negative control: the same series with (tau/2)^m instead of (tau/4)^m ==")
    ctx.prec = PREC
    z = acb(arb("171.2502"), arb("-0.1823"))
    s0 = acb(1)/2 + acb(0, 1)*z/2
    s0 = acb(arb(s0.real.mid()), arb(s0.imag.mid()))
    ctx.cap = 2*M + 2
    kap = kappa(); a4 = [acb(1), kap, -kap, acb(-1)]
    sser = acb_series([s0, acb(1)]); S = acb_series([acb(0)])
    for r in (1, 2, 3, 4):
        S = S + acb_series([a4[r-1]])*sser.zeta(acb(r)/5)
    five = (acb_series([-s0*acb(5).log(), -acb(5).log()])).exp()
    uu = acb_series([(s0+1)/2, acb(1)/2])
    cond = (uu*acb_series([acb(5).log() - acb.pi().log()])).exp()
    c = (cond*uu.gamma()*five*S).coeffs()
    while len(c) < 2*M + 2: c.append(acb(0))
    bad = acb(0)
    for m in range(M + 1):
        bad = bad + (acb("0.16")/2)**m * c[2*m]*fact(2*m)/fact(m)
    bad = bad/2
    mp.dps = 60
    badm = mpc(mpf(arb(bad.real.mid()).str(50, radius=False)),
               mpf(arb(bad.imag.mid()).str(50, radius=False)))
    ref = H_quad("171.2502", "-0.1823", "0.16", 140)
    relbad = abs(badm - ref)/abs(ref)
    print(f"   (tau/2)^m version : {mp.nstr(badm, 14)}")
    print(f"   the integral      : {mp.nstr(ref, 14)}")
    print(f"   relative deviation: {mp.nstr(relbad, 6)}   "
          f"{'PASS (it disagrees, as it must)' if relbad > mpf('0.01') else 'FAIL'}")
    ok &= bool(relbad > mpf("0.01"))

    print("\n=== check_taylor_vs_integral.py:", "ALL PASS" if ok else "SOME FAILED", "===")
    json.dump(dict(M=M, points=rows, wrong_constant_rel_dev=mp.nstr(relbad, 8)),
              open(os.path.join(os.path.dirname(__file__), "..", "data",
                                "taylor_vs_integral.json"), "w"), indent=1)
    return ok

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
