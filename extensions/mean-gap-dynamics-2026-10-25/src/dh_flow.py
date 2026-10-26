"""Control 2 -- the Davenport-Heilbronn family, which has certified off-line zeros at tau = 0.

D'(0) is computed from  z' = H''(z)/H'(z)  by ANALYTIC differentiation, so the
infinite zero sum is never truncated: the whole spectrum is present implicitly in
H_0 = Lambda_F(1/2+iz/2)/2.  (That is why this route is used instead of the
zero-interaction sum, whose tail would have to be bounded separately.)

Two independent derivative routes are compared:
  (R1)  the chain rule in s:  z' = (i/2) Lambda''(s)/Lambda'(s),  s = 1/2 + iz/2,
        with Lambda', Lambda'' assembled from analytic Hurwitz-zeta derivatives;
  (R2)  the quadrature route: H, H', H'' = int e^{tau u^2} Phi_F(u) u^k trig(zu) du,
        differentiating cos(zu) under the integral sign.
R2 shares no special-function code with R1 (it uses only exp/cos on the kernel).
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
from mpmath import mp, mpf, mpc, log, pi, gamma, psi, zeta, exp
import kernels as K

# the parent extension's census (T2 numerical inputs, NOT certified zero locations)
CENSUS = [
    ("O_1", "0.808517182456637386", "85.6993484853775922"),
    ("O_2", "0.650830080609737082", "114.163342730756981"),
    ("O_3", "0.574356050450805991", "166.479305913168156"),
    ("O_4", "0.724257694626809780", "176.702461242855825"),
]

def S_derivs(s, kap):
    """S(s) = sum_{r=1}^{4} a_r zeta(s, r/5) and its first two s-derivatives."""
    a = K.dh_coeffs(kap)
    out = []
    for n in (0, 1, 2):
        out.append(sum(a[r-1] * mp.zeta(s, mpf(r)/5, n) for r in range(1, 5)))
    return out

def Lambda_log_derivs(s, kap=None):
    """(Lambda, Lambda'/Lambda, Lambda''/Lambda) at s, analytically.

        Lambda(s) = (5/pi)^{(s+1)/2} Gamma((s+1)/2) F(s),   F(s) = 5^{-s} S(s)
        Lambda'/Lambda      = (1/2)log(5/pi) + (1/2)psi((s+1)/2) - log 5 + S'/S
        (Lambda'/Lambda)'   = (1/4)psi'((s+1)/2) + S''/S - (S'/S)^2
        Lambda''/Lambda     = (Lambda'/Lambda)^2 + (Lambda'/Lambda)'
    """
    if kap is None:
        kap = K.kappa_closed()
    S0, S1, S2 = S_derivs(s, kap)
    L = (mpf(5)/pi)**((s+1)/2) * gamma((s+1)/2) * mpf(5)**(-s) * S0
    l1 = log(mpf(5)/pi)/2 + psi(0, (s+1)/2)/2 - log(5) + S1/S0
    l1p = psi(1, (s+1)/2)/4 + S2/S0 - (S1/S0)**2
    l2 = l1*l1 + l1p
    return L, l1, l2

def zdot_R1(z, kap=None):
    """z'(0) = (i/2) Lambda''(s)/Lambda'(s),  s = 1/2 + i z/2."""
    s = mpf(1)/2 + mpc(0, 1)*z/2
    L, l1, l2 = Lambda_log_derivs(s, kap)
    return mpc(0, 1)/2 * (l2/l1)          # Lambda''/Lambda' = (L''/L)/(L'/L)

def main():
    ok = True
    rows = []
    mp.dps = 60
    kap = K.kappa_closed()

    print("== DH: D'(0) at the four census orbits, two independent routes ==")
    print("   D = (Im z)^2/2 = 2 delta^2 = A - M_H;   D'(0) = Im(z) * Im(z'(0))")
    print()
    hdr = f"{'orbit':>5} {'delta':>14} {'D':>14} {'Dprime (R1)':>20} {'orbit self-term':>18} {'cross-orbit rest':>18}"
    print(hdr)
    for name, bs, gs in CENSUS:
        b, g = mpf(bs), mpf(gs)
        d = b - mpf(1)/2
        z = 2*g - 2*mpc(0, 1)*d
        y, x = mp.im(z), mp.re(z)
        D = y*y/2
        zd = zdot_R1(z, kap)
        Dp = y * mp.im(zd)
        # the exact own-orbit contribution:  -1 (conjugate partner) - y^2/|z|^2 (mirror -z)
        selfterm = -1 - y*y/(x*x + y*y)
        rows.append(dict(orbit=name, beta=bs, gamma=gs, delta=mp.nstr(d, 18),
                         D=mp.nstr(D, 18), Dprime=mp.nstr(Dp, 18),
                         self_term=mp.nstr(selfterm, 18),
                         cross=mp.nstr(Dp - selfterm, 18)))
        print(f"{name:>5} {mp.nstr(d,10):>14} {mp.nstr(D,10):>14} {mp.nstr(Dp,14):>20} "
              f"{mp.nstr(selfterm,12):>18} {mp.nstr(Dp-selfterm,12):>18}")
        ok &= (Dp < 0)

    # ---- precision doubling
    print("\n== precision doubling on O_1 (dps 60 -> 120) ==")
    vals = []
    for dps in (60, 120):
        mp.dps = dps
        kap2 = K.kappa_closed()
        b, g = mpf(CENSUS[0][1]), mpf(CENSUS[0][2])
        z = 2*g - 2*mpc(0, 1)*(b - mpf(1)/2)
        vals.append(mp.im(z) * mp.im(zdot_R1(z, kap2)))
    mp.dps = 60
    dev = abs(vals[0] - vals[1])
    print(f"   dps  60: D'(0) = {mp.nstr(vals[0], 30)}")
    print(f"   dps 120: D'(0) = {mp.nstr(vals[1], 30)}")
    print(f"   |difference| = {mp.nstr(dev, 4)}")
    ok &= (dev < mpf(10)**-45)

    # ---- route R2: independent quadrature
    print("\n== R2: independent quadrature route on O_1 (no Hurwitz zeta, no psi) ==")
    mp.dps = 90
    b, g = mpf(CENSUS[0][1]), mpf(CENSUS[0][2])
    z = 2*g - 2*mpc(0, 1)*(b - mpf(1)/2)
    H0, H1, H2 = K.H_derivs(z, 0, "dh", 2)
    zd2 = H2/H1
    Dp2 = mp.im(z) * mp.im(zd2)
    mp.dps = 60
    zd1 = zdot_R1(z, K.kappa_closed())
    rel = abs(zd2 - zd1)/abs(zd1)
    print(f"   R1  z'(0) = {mp.nstr(zd1, 24)}")
    print(f"   R2  z'(0) = {mp.nstr(zd2, 24)}")
    print(f"   relative deviation R1 vs R2 = {mp.nstr(rel, 4)}")
    print(f"   R2  D'(0) = {mp.nstr(Dp2, 24)}")
    ok &= (rel < mpf(10)**-20)

    # ---- extremality
    print("\n== extremality: which orbits does the SIGN LEMMA actually apply to? ==")
    ds = [(n, abs(mpf(bs) - mpf(1)/2)) for n, bs, gs in CENSUS]
    ds_sorted = sorted(ds, key=lambda t: -t[1])
    print("   |delta| over the census, descending: " +
          ", ".join(f"{n} {mp.nstr(v,8)}" for n, v in ds_sorted))
    print(f"   O_1 is the census maximum, so WITHIN THE CENSUS every other zero")
    print(f"   contributes <= 0 to D'(O_1) and the lemma gives D'(O_1) <= -1.")
    print(f"   BUT the census covers only 0 < Im s < 200, while the proved strip is")
    print(f"   |delta| <= sigma_0 - 1/2 = 0.895136158235109721061359 -- almost three times")
    print(f"   |delta_1| = 0.3085.  Global extremality of O_1 is therefore NOT established,")
    print(f"   and no finite census can establish it.  For O_2, O_3, O_4 the hypothesis is")
    print(f"   FALSE already inside the census (O_1 lies further off the line), so for them")
    print(f"   the negative sign observed above is a NUMERICAL FACT, not a forced one.")

    print("\n=== dh_flow.py:", "ALL PASS" if ok else "SOME FAILED", "===")
    json.dump(dict(rows=rows, dps_doubling=[mp.nstr(v, 30) for v in vals],
                   R1_vs_R2_rel=mp.nstr(rel, 6)),
              open(os.path.join(os.path.dirname(__file__), "..", "data", "part4_dh.json"), "w"),
              indent=1)
    return ok

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
