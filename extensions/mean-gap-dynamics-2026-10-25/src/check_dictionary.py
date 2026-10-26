"""Stage 0 -- verify the normalisation dictionary BEFORE any dynamics is computed.

Checks, in order:
  A  H_0^zeta(z) = xi(1/2 + iz/2)/8                       (the zeta dictionary)
  B  d_tau H = -d_z^2 H for both kernels                  (the PDE and its sign)
  C  kappa from the Gauss sum = kappa closed form = the parent's committed value
  D  H_0^DH(z) = Lambda_F(1/2 + iz/2)/2                   (the DH dictionary)
  E  Lambda_F(s) = Lambda_F(1-s)                          (the functional equation)
  F  z = 2 gamma - 2 i delta, so delta = -Im z / 2, D = (Im z)^2/2
  G  the committed off-line zero really is a zero of H_0^DH
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
from mpmath import mp, mpf, mpc, exp, cos, pi, sqrt, quad
import kernels as K

OUT = []
def rec(tag, what, val, tol, extra=""):
    ok = (val is not None) and (abs(val) <= tol)
    OUT.append(dict(tag=tag, check=what, residual=mp.nstr(mp.mpf(val), 8), tol=mp.nstr(mp.mpf(tol),3), pass_=bool(ok), note=extra))
    print(f"[{tag}] {'PASS' if ok else 'FAIL'}  {what}: residual {mp.nstr(mp.mpf(val),8)}  (tol {mp.nstr(mp.mpf(tol),3)}) {extra}")
    return ok

# committed values of the parent extensions, quoted as CHECK TARGETS
COMMITTED = dict(
    kappa="0.2840790438404122960282918323931261690911",
    sigma0="1.395136158235109721061359",
    beta1="0.8085171824566373855533519606068441278507",
    gamma1="85.69934848537759217192926770894172903799",
    Re_z1="171.398696970755184343858535418",
    Im_z1="-0.617034364913274771106703921214",
    tau_DH="1.602537483559822777787819",
)

allok = True
mp.dps = 90

# ---- A  zeta dictionary
print("\n== A. zeta dictionary  H_0(z) = xi(1/2+iz/2)/8 ==")
for zv in ["0", "1", "5", "14.134725", "28.2694545"]:
    z = mp.mpf(zv)
    lhs = K.H_tau(z, 0, "zeta")
    s = mp.mpf(1)/2 + mpc(0,1)*z/2
    rhs = K.xi(s)/8
    r = abs(lhs-rhs)/max(abs(rhs), mp.mpf(10)**(-30))
    allok &= rec("A", f"z={zv}", r, mp.mpf(10)**(-25), "relative")

# ---- B  the PDE, both kernels, by exact differentiation in z and finite difference in tau
print("\n== B. d_tau H = - d_z^2 H ==")
h = mp.mpf(10)**(-30)
for kern, zs in (("zeta", ["3.0","28.2694545"]), ("dh", ["2.0","171.4"])):
    for zv in zs:
        for tau in (mp.mpf(0), mp.mpf("0.2")):
            z = mp.mpf(zv)
            dtau = (K.H_tau(z, tau+h, kern) - K.H_tau(z, tau-h, kern))/(2*h)
            _,_,H2 = K.H_derivs(z, tau, kern, 2)
            scale = max(abs(dtau), abs(H2), mp.mpf(10)**(-40))
            allok &= rec("B", f"{kern} z={zv} tau={mp.nstr(tau,3)}", (dtau + H2)/scale,
                         mp.mpf(10)**(-18), "relative; FD step 1e-30, dps 90")

# ---- C  kappa
print("\n== C. kappa: Gauss sum vs closed form vs the committed value ==")
mp.dps = 60
kg, kc = K.kappa_gauss(), K.kappa_closed()
allok &= rec("C", "Gauss sum vs closed form", kg-kc, mp.mpf(10)**(-45))
allok &= rec("C", "closed form vs committed", kc - mp.mpf(COMMITTED["kappa"]), mp.mpf(10)**(-38))

# ---- D/E  DH dictionary and functional equation
print("\n== D. DH dictionary  H_0(z) = Lambda_F(1/2+iz/2)/2 ==")
mp.dps = 90
for zv in ["0", "2", "10"]:
    z = mp.mpf(zv)
    lhs = K.H_tau(z, 0, "dh")
    rhs = K.Lambda_DH(mp.mpf(1)/2 + mpc(0,1)*z/2)/2
    allok &= rec("D", f"z={zv}", abs(lhs-rhs)/max(abs(rhs),mp.mpf(10)**(-30)), mp.mpf(10)**(-20), "relative")

print("\n== E. Lambda_F(s) = Lambda_F(1-s) ==")
mp.dps = 60
for sv in [mpc("0.3","4.7"), mpc("1.1","12.3"), mpc("0.81","85.7")]:
    L1, L2 = K.Lambda_DH(sv), K.Lambda_DH(1-sv)
    allok &= rec("E", f"s={mp.nstr(sv,6)}", abs(L1-L2)/max(abs(L1),mp.mpf(10)**(-60)), mp.mpf(10)**(-45), "relative")

# ---- F  the coordinate map
print("\n== F. coordinate map  z = 2 gamma - 2 i delta ==")
b1, g1 = mp.mpf(COMMITTED["beta1"]), mp.mpf(COMMITTED["gamma1"])
d1 = b1 - mp.mpf(1)/2
z1 = 2*g1 - 2*mpc(0,1)*d1
allok &= rec("F", "Re z1 vs committed", mp.re(z1) - mp.mpf(COMMITTED["Re_z1"]), mp.mpf(10)**(-25))
allok &= rec("F", "Im z1 vs committed", mp.im(z1) - mp.mpf(COMMITTED["Im_z1"]), mp.mpf(10)**(-25))
allok &= rec("F", "delta = -Im z/2", K.delta_of_z(z1) - d1, mp.mpf(10)**(-50))
allok &= rec("F", "D = 2 delta^2 = (Im z)^2/2", K.D_of_z(z1) - 2*d1**2, mp.mpf(10)**(-50))
s0 = mp.mpf(COMMITTED["sigma0"])
allok &= rec("F", "tau_DH = Delta^2/2 = 2(sigma0-1/2)^2 = D_max",
             2*(s0-mp.mpf(1)/2)**2 - mp.mpf(COMMITTED["tau_DH"]), mp.mpf(10)**(-22),
             "-> the de Bruijn threshold IS sup D")

# ---- G  the committed zero is a zero
print("\n== G. H_0^DH(z1) = 0 ==")
mp.dps = 60
s1 = mp.mpf(1)/2 + mpc(0,1)*z1/2      # = beta1 + i gamma1
allok &= rec("G", "s1 == beta1 + i gamma1", abs(s1 - (b1 + mpc(0,1)*g1)), mp.mpf(10)**(-50))
LF = K.Lambda_DH(s1)
FF = K.F_DH(s1)
allok &= rec("G", "|F(rho1)| (scale-free)", abs(FF), mp.mpf(10)**(-30), "absolute; F is O(1) nearby")

print("\n=== OVERALL:", "ALL PASS" if allok else "SOME FAILED", "===")
json.dump(OUT, open(os.path.join(os.path.dirname(__file__), "..", "data", "stage0_dictionary.json"), "w"), indent=1)
sys.exit(0 if allok else 1)
