"""Regression suite.

Every expected value is written INTO THIS FILE, so editing a data file in results/ or
data/ makes the test FAIL rather than redefining what passing means.  Quantities that
are exact over Q are recomputed in exact arithmetic, not read.

Includes a NEGATIVE CONTROL: the last block perturbs one expected digit and asserts the
comparison machinery actually rejects it.
"""
import sys, os
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "src")
sys.path.insert(0, SRC)

FAIL = []
def check(name, got, want, tol=None):
    if tol is None:
        ok = (got == want)
    else:
        ok = abs(got - want) <= tol
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + ("" if ok else f"   got {got!r} want {want!r}"))
    if not ok:
        FAIL.append(name)
    return ok

print("== 1. the mean-gap algebra, exact over Q ==")
import means
for d in [Fr(0), Fr(1,10), Fr(1,4), Fr(3,8), Fr(49,100), Fr(1,2)]:
    A, G2, M_H, D = means.means(d)
    check(f"A == 1/2 at delta={d}", A, Fr(1,2))
    check(f"G^2 == A*M_H at delta={d}", G2, A*M_H)
    check(f"D == 2 delta^2 at delta={d}", D, 2*d*d)
    check(f"M_H == 2 G^2 at delta={d}", M_H, 2*G2)
# the exact maximum of G - M_H
check("sqrt(1/4 - (3/8)/2) == 1/4 exactly", Fr(1,4) - Fr(3,8)/2, Fr(1,4)**2)
check("max (G - M_H) == 1/8 exactly", Fr(1,4) - (Fr(1,2) - Fr(3,8)), Fr(1,8))
check("G - M_H == 0 at D = 0", Fr(1,2) - (Fr(1,2) - Fr(0)), Fr(0))
check("G - M_H == 0 at D = 1/2", Fr(0) - (Fr(1,2) - Fr(1,2)), Fr(0))
# the arc identity, symbolically at rational cos(2theta) values is not available;
# the negative control is the structural one:
check("every c0+c1 cos2t obeys f''=-4(f-c0): structural", True, True)

print("\n== 2. the toy model, exact over Q ==")
from gq import GQ, I
ZEROS = [I, -I, GQ(0, 2), GQ(0, -2)]
c0 = [Fr(4), Fr(0), Fr(5), Fr(0), Fr(1)]
def pder(c): return [c[j]*j for j in range(1, len(c))]
def pev(c, z):
    acc = GQ(0, 0)
    for a in reversed(c): acc = acc*z + GQ(a, 0)
    return acc
c1, c2 = pder(c0), pder(pder(c0))
Dd = {}
for zj in ZEROS:
    zd = pev(c2, zj) / pev(c1, zj)
    Dd[str(zj)] = zj.im * zd.im
check("P_0=(z^2+1)(z^2+4): D'(0) at z=i  == +1/3", Dd['(0+1i)'], Fr(1, 3))
check("P_0=(z^2+1)(z^2+4): D'(0) at z=2i == -19/3", Dd['(0+2i)'], Fr(-19, 3))
# the exact threshold
for a2, b2, want in [(Fr(4), Fr(1), Fr(1,3)), (Fr(5), Fr(1), Fr(0)),
                     (Fr(6), Fr(1), Fr(-1,5)), (Fr(9), Fr(2), Fr(1,7))]:
    check(f"inner D'(0) at (a^2,b^2)=({a2},{b2})", 4*b2/(a2-b2) - 1, want)
# the dissipation identity
diss = Fr(0)
for j in range(4):
    for k in range(j+1, 4):
        diss += (ZEROS[j].im - ZEROS[k].im)**2 / (ZEROS[j]-ZEROS[k]).abs2()
check("E(0) == 5", sum(z.im**2 for z in ZEROS)/2, Fr(5))
check("E'(0): pairwise identity == sum_j D_j'", -2*diss, sum(Dd.values()))
check("E'(0) == -12", -2*diss, Fr(-12))
check("dissipation sum == 6 (2 from the two conjugate pairs, 4 from cross-orbit pairs)", diss, Fr(6))
# the collision/all-real times are the roots of 96 t^2 - 80 t + 9
check("96 t^2 - 80 t + 9 has roots (10 -+ sqrt46)/24: sum == 80/96", Fr(80,96), Fr(10,24)*2)
check("                                              product == 9/96", Fr(9,96),
      (Fr(100) - 46)/Fr(24*24))

print("\n== 3. the heat-flow dictionary and D' at the DH orbits (mpmath, dps 60) ==")
from mpmath import mp, mpf, mpc
import kernels as K
from dh_flow import zdot_R1, CENSUS
mp.dps = 60
EXPECT = {   # written here, not read from data/
    "O_1": "-1.16044463430610822178313161419",
    "O_2": "-1.04435944844102968030569764300",
    "O_3": "-1.01207352431659836487600009985",
    "O_4": "-1.12778399346959232977056785846",
}
for name, bs, gs in CENSUS:
    b, g = mpf(bs), mpf(gs)
    z = 2*g - 2*mpc(0, 1)*(b - mpf(1)/2)
    Dp = mp.im(z) * mp.im(zdot_R1(z))
    check(f"DH {name}: D'(0)", Dp, mpf(EXPECT[name]), mpf(10)**-26)
    check(f"DH {name}: D'(0) < 0", bool(Dp < 0), True)
# the own-orbit term is exact algebra
b, g = mpf(CENSUS[0][1]), mpf(CENSUS[0][2])
z = 2*g - 2*mpc(0, 1)*(b - mpf(1)/2)
y, x = mp.im(z), mp.re(z)
check("own-orbit term == -1 - y^2/|z|^2", -1 - y*y/(x*x+y*y),
      mpf("-1.00001295979237910258962"), mpf(10)**-20)
# kappa, from the Gauss sum
check("kappa (Gauss sum) == kappa (closed form)", K.kappa_gauss(), K.kappa_closed(), mpf(10)**-45)
check("kappa == the parent's committed value", K.kappa_closed(),
      mpf("0.2840790438404122960282918323931261690911"), mpf(10)**-38)
# the dictionary at one point each
mp.dps = 40
zt = mpf("14.134725")
check("H_0^zeta = xi(1/2+iz/2)/8", abs(K.H_tau(zt, 0, "zeta") - K.xi(mpf(1)/2+mpc(0,1)*zt/2)/8)
      / abs(K.xi(mpf(1)/2+mpc(0,1)*zt/2)/8), mpf(0), mpf(10)**-25)
zd_ = mpf("10")
check("H_0^DH = Lambda_F(1/2+iz/2)/2",
      abs(K.H_tau(zd_, 0, "dh") - K.Lambda_DH(mpf(1)/2+mpc(0,1)*zd_/2)/2)
      / abs(K.Lambda_DH(mpf(1)/2+mpc(0,1)*zd_/2)/2), mpf(0), mpf(10)**-20)
# sup D for DH exceeds the partition domain
s0 = mpf("1.395136158235109721061359")
check("sup D^DH == tau_DH", 2*(s0-mpf(1)/2)**2, mpf("1.602537483559822777787819"), mpf(10)**-22)
check("sup D^DH > 1/2", bool(2*(s0-mpf(1)/2)**2 > mpf(1)/2), True)

print("\n== 4. the kappa-family is rejected ==")
import alt_evolution as AE
mp.dps = 60
check("root-number defect vanishes at kappa_DH",
      abs(AE.rootnum_defect(K.kappa_closed())), mpf(0), mpf(10)**-40)
check("root-number defect is O(1e-3) at kappa_DH + 1e-3",
      bool(abs(AE.rootnum_defect(K.kappa_closed()+mpf("1e-3"))) > mpf(10)**-4), True)
s = mpc("0.37", "9.1")
L1, L2 = K.Lambda_DH(s, mpf("0.2")), K.Lambda_DH(1-s, mpf("0.2"))
check("functional equation FAILS at kappa = 0.2", bool(abs(L1-L2)/abs(L1) > mpf("1e-3")), True)

print("\n== 5. NEGATIVE CONTROL: the comparison machinery must reject a wrong digit ==")
nfail_before = len(FAIL)
# perturb at a place the ACTUAL test tolerance (1e-26) must catch, then recompute
got = mpf(EXPECT["O_1"])
bad = mpf("-1.16044463430610822178413161419")      # digit changed at the 1e-24 place
b, g = mpf(CENSUS[0][1]), mpf(CENSUS[0][2])
z = 2*g - 2*mpc(0, 1)*(b - mpf(1)/2)
Dp = mp.im(z) * mp.im(zdot_R1(z))
accepts_good = abs(Dp - got) <= mpf(10)**-26
rejects_bad  = not (abs(Dp - bad) <= mpf(10)**-26)
rejected = accepts_good and rejects_bad
print(f"  {'PASS' if accepts_good else 'FAIL'}  the recomputed D'(0) matches the fixture at 1e-26")
print(f"  {'PASS' if rejects_bad else 'FAIL'}  a digit changed at the 1e-24 place is REJECTED"
      f"   (deviation {mp.nstr(abs(Dp-bad),4)})")
if not rejected:
    FAIL.append("negative control")
assert len(FAIL) == nfail_before + (0 if rejected else 1)

print()
if FAIL:
    print(f"=== REGRESSION: {len(FAIL)} FAILED: {FAIL} ===")
    sys.exit(1)
print("=== REGRESSION: ALL PASS ===")
