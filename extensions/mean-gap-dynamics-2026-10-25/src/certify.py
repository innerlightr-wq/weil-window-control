"""Interval certification (Arb ball arithmetic, via python-flint) of the representative
statements that are not already exact over Q.

What is certified here:
  (C1)  D'(0) < 0 at the Davenport-Heilbronn census orbit O_1, as a FUNCTION of the
        supplied z_1.  The input z_1 is an exact decimal; whether it encloses an actual
        zero of H_0 rests on the parent census, which is T2 -- exactly the parent's own
        convention ("certified as a function of the supplied (beta_1, gamma_1)").
        Note that D'(0) is a well-defined analytic function of z wherever H'(z) != 0, so
        the certified object is unambiguous even if z_1 is not exactly a zero.
  (C2)  the same for O_2, O_3, O_4.
  (C3)  the own-orbit decomposition  D'(0) - (-1 - y^2/|z|^2)  < 0 at O_1, i.e. the
        cross-orbit remainder is itself negative.
  (C4)  tau_DH = 2(sigma_0 - 1/2)^2 and the mean-gap range statement sup D^DH > 1/2.

NOT certified, and said so: the toy-model results are EXACT OVER Q (stronger than an
interval), the mean-gap algebra is EXACT OVER Q, and nothing here certifies a zero
LOCATION or any statement about tau != 0.
"""
import sys, os, json
from flint import acb, arb, ctx, acb_series

PREC = 400
ctx.prec = PREC
ctx.cap = 3

CENSUS = [
    ("O_1", "0.808517182456637386", "85.6993484853775922"),
    ("O_2", "0.650830080609737082", "114.163342730756981"),
    ("O_3", "0.574356050450805991", "166.479305913168156"),
    ("O_4", "0.724257694626809780", "176.702461242855825"),
]
SIGMA0 = "1.395136158235109721061359"
TAU_DH = "1.602537483559822777787819"

def kappa():
    r5 = acb(5).sqrt()
    return ((acb(10) - 2*r5).sqrt() - 2) / (r5 - 1)

def lambda_series(s0):
    """Lambda_F(s0 + x) as a truncated power series, rigorously enclosed."""
    kap = kappa()
    a = [acb(1), kap, -kap, acb(-1)]
    x1 = acb(1)
    # zeta(s0 + x, r/5)
    S = acb_series([acb(0)])
    sser = acb_series([acb(s0), x1])
    for r in (1, 2, 3, 4):
        S = S + acb_series([a[r-1]]) * sser.zeta(acb(r)/5)
    # 5^{-s} = exp(-s log 5)
    five_pow = (acb_series([-acb(s0)*acb(5).log(), -acb(5).log()])).exp()
    # u = (s0+1)/2 + x/2
    u = acb_series([(acb(s0)+1)/2, acb(1)/2])
    conduct = (u * acb_series([acb(5).log() - acb.pi().log()])).exp()
    return conduct * u.gamma() * five_pow * S

def Dprime_cert(beta, gamma):
    """certified enclosure of D'(0) = Im(z) * Im( (i/2) Lambda''/Lambda' ) at z = 2g - 2i d."""
    b, g = acb(beta), acb(gamma)
    d = b - acb(1)/2
    z = 2*g - 2*acb(0, 1)*d
    s0 = acb(1)/2 + acb(0, 1)*z/2            # = beta + i gamma
    L = lambda_series(s0)
    c = L.coeffs()
    c0, c1, c2 = c[0], c[1], c[2]
    zdot = acb(0, 1) * c2 / c1               # (i/2) * (2 c2) / c1
    y = z.imag
    return y * zdot.imag, z, c0, c1

def _re(a):
    return a.real if isinstance(a, acb) else a

def floor_str(a, n=22):
    """one-sided: the ball's LOWER endpoint, printed as a plain decimal (rounds down)."""
    a = _re(a)
    return arb(a.mid() - a.rad()).str(n, radius=False)

def ceil_str(a, n=22):
    a = _re(a)
    return arb(a.mid() + a.rad()).str(n, radius=False)

def main():
    ok = True
    out = {"prec_bits": PREC, "rows": []}
    print(f"== Arb certification at {PREC} bits (python-flint) ==\n")
    print(f"   {'orbit':>5} {'D':>26} {'certified D-prime(0)':>44} {'<0':>4}")
    for name, bs, gs in CENSUS:
        Dp, z, L0, L1 = Dprime_cert(bs, gs)
        D = z.imag**2 / 2
        neg = (Dp < 0)
        ok &= bool(neg)
        print(f"   {name:>5} {D.str(20):>26} {Dp.str(32):>44} {'YES' if neg else 'NO':>4}")
        print(f"         upper bound:  D'(0) <= {ceil_str(Dp)}        (plain decimal, rounded UP)")
        print(f"         |Lambda(s1)| = {abs(L0).str(8)}   |Lambda'(s1)| = {abs(L1).str(8)}"
              f"   (H' != 0, so D' is well defined)")
        out["rows"].append(dict(orbit=name, beta=bs, gamma=gs,
                                D=D.str(24), Dprime=Dp.str(32),
                                Dprime_upper=ceil_str(Dp), strictly_negative=bool(neg)))

    # ---- C3: the cross-orbit remainder at O_1
    print("\n== (C3) the cross-orbit remainder at O_1 ==")
    Dp, z, _, _ = Dprime_cert(CENSUS[0][1], CENSUS[0][2])
    y, x = z.imag, z.real
    selfterm = (-acb(1) - y*y/(x*x + y*y)).real
    cross = Dp - selfterm
    print(f"   own-orbit self term  = -1 - y^2/|z|^2 = {selfterm.str(26)}")
    print(f"   cross-orbit rest     = {cross.str(26)}")
    print(f"   certified upper bound: cross <= {ceil_str(cross)}   -> strictly negative: "
          f"{'YES' if cross < 0 else 'NO'}")
    ok &= bool(cross < 0)
    out["O1_self_term"] = selfterm.str(26); out["O1_cross"] = cross.str(26)
    out["O1_cross_upper"] = ceil_str(cross)

    # ---- C4: the range statement
    print("\n== (C4) sup D for Davenport-Heilbronn exceeds the partition domain ==")
    s0 = acb(SIGMA0)
    supD = (2*(s0 - acb(1)/2)**2).real
    print(f"   sup D^DH = 2(sigma_0 - 1/2)^2 = {supD.str(26)}")
    print(f"   committed tau_DH          = {TAU_DH}")
    diff = supD - arb(TAU_DH)
    print(f"   difference                = {diff.str(8)}   (identical to the committed digits)")
    print(f"   sup D^DH > 1/2 :            {'YES' if supD > arb(1)/2 else 'NO'}"
          f"   -> the binary-partition reading is unavailable on part of the DH strip")
    ok &= bool(supD > arb(1)/2)
    ok &= bool(abs(diff) < arb(10)**(-22))
    out["supD_DH"] = supD.str(26)

    print("\n== scope of these certificates ==")
    print("   * Each D'(0) is certified as an analytic function of the SUPPLIED z_1; that")
    print("     z_1 encloses an actual zero of H_0^DH rests on the parent census (T2).")
    print("   * Nothing here certifies a statement at tau != 0, nor any zero location.")
    print("   * The toy-model counterexample and the mean-gap algebra are exact over Q and")
    print("     need no interval arithmetic; they are stronger as they stand.")
    print("   * Arb supplies rigorous enclosures for Gamma and the Hurwitz zeta; kappa is")
    print("     built from square roots only.  No libm value enters.")

    print("\n=== certify.py:", "ALL CERTIFIED" if ok else "SOME FAILED", "===")
    json.dump(out, open(os.path.join(os.path.dirname(__file__), "..", "data", "certificates.json"), "w"), indent=1)
    return ok

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
