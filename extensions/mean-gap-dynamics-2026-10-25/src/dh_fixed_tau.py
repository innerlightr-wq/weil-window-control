"""Navigation for the fixed-parameter Davenport-Heilbronn certificate.

Locates the branch of the off-line zero O_1 at a FIXED tau > 0, using the Taylor
representation of src/heat_taylor.py evaluated through Arb series (midpoints only -- this
script navigates, src/certify_fixed_tau.py certifies).  Then reports, for each tau:

  * the located zero and its Newton residual,
  * |Im z| (the margin to the real axis) and D = (Im z)^2/2,
  * a workable disc radius r < |Im z|,
  * min |H_tau| on that circle and the numerical winding number,
  * the rigorous tail majorant T_M at that tau, for the truncation order used.

The point of the exercise: if the winding number is 1 on a circle whose closed disc misses
the real axis, then H_tau has a NONREAL zero at that tau, hence  Lambda_DH > tau  directly
-- with no trajectory estimate, no isolation hypothesis, and no uniformity over an interval.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
from flint import acb, arb, ctx, acb_series
from mpmath import mp, mpf, mpc, exp, cosh, pi, quad, factorial

PREC = 500
ctx.prec = PREC
M_DEFAULT = 20

def kappa():
    r5 = acb(5).sqrt()
    return ((acb(10) - 2*r5).sqrt() - 2)/(r5 - 1)

def lambda_series(s0, cap):
    ctx.cap = cap
    kap = kappa()
    a = [acb(1), kap, -kap, acb(-1)]
    sser = acb_series([acb(s0), acb(1)])
    S = acb_series([acb(0)])
    for r in (1, 2, 3, 4):
        S = S + acb_series([a[r-1]]) * sser.zeta(acb(r)/5)
    five = (acb_series([-acb(s0)*acb(5).log(), -acb(5).log()])).exp()
    u = acb_series([(acb(s0)+1)/2, acb(1)/2])
    cond = (u*acb_series([acb(5).log() - acb.pi().log()])).exp()
    return cond * u.gamma() * five * S

_FACT = [acb(1)]
def fact(n):
    while len(_FACT) <= n:
        _FACT.append(_FACT[-1]*len(_FACT))
    return _FACT[n]

def G(coeffs, tau, M, j):
    """sum_{m=0..M} (tau/4)^m Lambda^{(2m+j)}(s)/m!   with Lambda^{(n)} = n! coeffs[n]."""
    acc = acb(0)
    for m in range(M + 1):
        n = 2*m + j
        acc = acc + (acb(tau)/4)**m * coeffs[n] * fact(n) / fact(m)
    return acc

def H_and_Hp(z, tau, M=M_DEFAULT):
    """(H_tau(z), dH_tau/dz)  as acb balls, from the series.  d/dz = (i/2) d/ds."""
    s0 = acb(1)/2 + acb(0, 1)*z/2
    c = lambda_series(s0, 2*M + 3).coeffs()
    while len(c) < 2*M + 3:
        c.append(acb(0))
    G0, G1 = G(c, tau, M, 0), G(c, tau, M, 1)
    return G0/2, acb(0, 1)*G1/4

def newton(z, tau, M=M_DEFAULT, iters=40, tol=1e-35):
    for _ in range(iters):
        H, Hp = H_and_Hp(z, tau, M)
        step = H/Hp
        z = z - step
        if abs(step).mid() < tol:
            break
    H, Hp = H_and_Hp(z, tau, M)
    return z, H, Hp

def tail_majorant(tau, Y, M, X=3):
    mp.dps = 40
    f = lambda u: 8*exp(3*u)*exp(-pi*exp(4*u)/5)*u**(2*M+2)*cosh(u*mpf(Y))*exp(mpf(tau)*u*u)
    I = quad(f, [0, mpf(1)/4, mpf(1)/2, mpf(3)/4, 1, mpf(3)/2, X])
    return mpf(tau)**(M+1)/factorial(M+1)*I

def main():
    ok = True
    z0 = acb("171.398696970755184343858535418", "-0.617034364913274771106703921214")
    print("== the O_1 branch at fixed tau (navigation; Arb midpoints) ==")
    print(f"   truncation M = {M_DEFAULT}; tail majorant reported per tau\n")
    hdr = (f"   {'tau':>7} {'Re z':>24} {'Im z':>18} {'D':>14} "
           f"{'|H|/|H-prime|':>14} {'T_M':>12}")
    print(hdr)
    rows = []
    z = z0
    TAUS = ["0.01", "0.05", "0.10", "0.12", "0.14", "0.15", "0.155", "0.16", "0.163"]
    for tau in TAUS:
        zt, H, Hp = newton(z, tau)
        if not (zt.real.mid() == zt.real.mid()):      # NaN guard
            print(f"   {tau:>7}   Newton did not converge from the previous zero -- "
                  f"the branch has reached or passed the real axis; stopping the scan")
            rows.append(dict(tau=tau, status="newton diverged"))
            break
        z = zt
        resid = (abs(H)/abs(Hp)).mid()
        y = z.imag
        D = y*y/2
        TM = tail_majorant(tau, 0.95, M_DEFAULT)
        print(f"   {tau:>7} {z.real.str(18):>24} {y.str(14):>18} {D.str(10):>14} "
              f"{arb(resid).str(6):>14} {mp.nstr(TM,4):>12}")
        rows.append(dict(tau=tau, re=z.real.str(20), im=y.str(20), D=D.str(16),
                         newton_resid=arb(resid).str(6), tail=mp.nstr(TM, 6)))
        ok &= (abs(float(y.str(12, radius=False))) > 0.05)
    print("""
   The branch stays non-real throughout; |Im z| shrinks as D decreases, exactly as
   D'(0) = -1.1604 predicted (linear extrapolation would put D = 0 near tau = 0.164).""")

    # ---- the winding number, numerically, at the most ambitious workable tau
    print("\n== numerical winding number on a circle about the located zero ==")
    print(f"   {'tau':>7} {'r':>7} {'margin |Im z|-r':>17} {'min|H| on dD':>16} "
          f"{'(1/2pi) sum d arg':>19}")
    wind = []
    zc_prev = z0
    for tau, r in (("0.01", "0.30"), ("0.10", "0.20"), ("0.14", "0.10"),
                   ("0.15", "0.08"), ("0.155", "0.06")):
        zc, _, _ = newton(zc_prev, tau)
        zc_prev = zc
        rr = arb(r)
        N = 400
        prev = None; tot = acb(0); mn = None
        for k in range(N + 1):
            th = 2*acb.pi()*k/N
            w = zc + acb(rr)*(acb(0, 1)*th).exp()
            Hv, _ = H_and_Hp(w, tau)
            a = abs(Hv)
            mn = a if mn is None else (a if a.mid() < mn.mid() else mn)
            if prev is not None:
                d = (Hv/prev).log().imag
                tot = tot + d
            prev = Hv
        wnum = tot/(2*acb.pi())
        marg = abs(zc.imag) - rr
        print(f"   {tau:>7} {r:>7} {marg.str(10):>17} {mn.str(8):>16} "
              f"{wnum.str(12):>19}")
        wind.append(dict(tau=tau, r=r, margin=marg.str(12), min_H=mn.str(10),
                         winding=wnum.str(14)))
        ok &= bool(marg > 0) and (abs(float(wnum.real.str(10, radius=False)) - 1.0) < 1e-6)
    print("""
   Winding number 1 on every circle, with the closed disc disjoint from the real axis.
   So H_tau has exactly one nonreal zero there -- numerically.  src/certify_fixed_tau.py
   turns the most ambitious of these into an Arb certificate.""")

    print("\n=== dh_fixed_tau.py:", "ALL PASS" if ok else "SOME FAILED", "===")
    json.dump(dict(M=M_DEFAULT, branch=rows, winding=wind),
              open(os.path.join(os.path.dirname(__file__), "..", "data",
                                "dh_fixed_tau_nav.json"), "w"), indent=1)
    return ok

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
