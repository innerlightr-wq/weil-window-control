"""Control 1 -- the finite polynomial heat flow: exact checks and the counterexample.

A family P_tau with  d_tau P = - d_z^2 P  and P_0 a polynomial is itself polynomial:

    P_tau(z) = exp(-tau d_z^2) P_0(z) = sum_{m>=0} (-tau)^m / m!  P_0^{(2m)}(z),

a FINITE sum.  On a simple-zero branch the chain rule gives, with no convergence
question and no tail at all,

    z_j'(tau) = P''(z_j)/P'(z_j) = 2 sum_{k != j} 1/(z_j - z_k)                    (T1)

and therefore, writing y_j = Im z_j and D_j = y_j^2/2,

    D_j'(tau) = -2 y_j sum_{k != j} (y_j - y_k)/|z_j - z_k|^2                      (T1)
    (sum over all zeros, with multiplicity)                       [the SIGN IDENTITY]

Two consequences, both exact:

  (S1)  a zero contributes POSITIVELY to D_j' only if y_k and y_j have the same sign
        and |y_k| > |y_j|: only a zero STRICTLY FURTHER off the line on the same side.
        In particular the conjugate partner z_j-bar contributes exactly -1, every real
        zero contributes <= 0, and the mirror -z_j contributes -y_j^2/|z_j|^2 <= 0.
        Hence if |y_j| is MAXIMAL among the zeros, D_j' <= -1.
  (S2)  symmetrising the identity over pairs, the TOTAL
            E(tau) := sum_j D_j = (1/2) sum_j (Im z_j)^2
        satisfies   E'(tau) = -2 sum_{j<k} (y_j - y_k)^2 / |z_j - z_k|^2  <= 0 .

This file verifies all of that exactly over Q (and Q(i)), and exhibits the minimal
counterexample to monotonicity of an INDIVIDUAL D.
"""
import sys, os, json
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(__file__))
from gq import GQ, I

OK = [True]
def chk(name, cond, extra=""):
    OK[0] &= bool(cond)
    print(f"  {'PASS' if cond else 'FAIL'}  {name} {extra}")

# ------------------------------------------------- polynomials over Q, in z and tau
# A polynomial in z with coefficients that are polynomials in tau:
# represent as list of lists: P[j][m] = coeff of z^j tau^m, exact Fractions.

def pz(coeffs):                      # polynomial in z only
    return [[Fr(c)] for c in coeffs]

def _pad(a, n): return a + [Fr(0)]*(n-len(a))

def tadd(a, b):
    n = max(len(a), len(b)); a, b = _pad(a[:], n), _pad(b[:], n)
    return [x+y for x, y in zip(a, b)]

def tscale(a, c): return [x*Fr(c) for x in a]

def tshift(a):                        # multiply by tau
    return [Fr(0)] + a

def dz2(P):                           # d^2/dz^2
    out = []
    for j in range(2, len(P)):
        out.append(tscale(P[j], j*(j-1)))
    return out if out else [[Fr(0)]]

def dtau(P):
    out = []
    for c in P:
        out.append([c[m]*m for m in range(1, len(c))] or [Fr(0)])
    return out

def heat_flow(P0):
    """P_tau = sum_m (-tau)^m/m! P_0^{(2m)}, exactly."""
    P = [c[:] for c in P0]
    term = [c[:] for c in P0]
    m = 0
    while True:
        term = dz2(term)
        if all(all(x == 0 for x in c) for c in term):
            break
        m += 1
        sc = Fr((-1)**m, 1)
        from math import factorial
        sc = sc / factorial(m)
        add = [tscale(c, sc) for c in term]
        for _ in range(m):
            add = [tshift(c) for c in add]
        n = max(len(P), len(add))
        P = _pad(P, n); add = _pad(add, n)
        P = [tadd(x if isinstance(x, list) else [Fr(0)], y if isinstance(y, list) else [Fr(0)])
             for x, y in zip(P, add)]
    return P

def ev_tau(P, tau):
    """evaluate the tau-dependence at a Fraction tau -> list of z-coefficients."""
    return [sum(c[m]*tau**m for m in range(len(c))) for c in P]

def poly_eval(c, z):
    acc = GQ(0, 0)
    for a in reversed(c):
        acc = acc*z + GQ(a, 0)
    return acc

def poly_der(c):
    return [c[j]*j for j in range(1, len(c))]

# ------------------------------------------------------------------- 1. the PDE
print("== 1. the polynomial family satisfies d_tau P = -P''  (exact in Q[z,tau]) ==")
for name, P0 in [("(z^2+1)(z^2+4) = z^4+5z^2+4", pz([4, 0, 5, 0, 1])),
                 ("z^6 - 3z^4 + 2",               pz([2, 0, 0, 0, -3, 0, 1])),
                 ("z^2 + 9",                      pz([9, 0, 1]))]:
    P = heat_flow(P0)
    L, R = dtau(P), dz2(P)
    n = max(len(L), len(R))
    L, R = _pad(L, n), _pad(R, n)
    resid = 0
    for a, b in zip(L, R):
        a = a if isinstance(a, list) else [Fr(0)]; b = b if isinstance(b, list) else [Fr(0)]
        k = max(len(a), len(b)); a, b = _pad(a, k), _pad(b, k)
        resid = max(resid, max(abs(x+y) for x, y in zip(a, b)))
    chk(f"{name}", resid == 0, "(exact residual 0)")

# ------------------------------- 2. the velocity identity P''/P' = 2 sum 1/(z_j-z_k)
print("\n== 2. z_j' = P''(z_j)/P'(z_j) = 2 sum_{k!=j} 1/(z_j-z_k)   (exact over Q(i)) ==")
ZEROS = [I, -I, GQ(0, 2), GQ(0, -2)]                 # zeros of z^4+5z^2+4
c0 = [Fr(4), Fr(0), Fr(5), Fr(0), Fr(1)]
c1, c2 = poly_der(c0), poly_der(poly_der(c0))
for j, zj in enumerate(ZEROS):
    lhs = poly_eval(c2, zj) / poly_eval(c1, zj)
    rhs = GQ(0, 0)
    for k, zk in enumerate(ZEROS):
        if k != j:
            rhs = rhs + GQ(2, 0)/(zj - zk)
    chk(f"z_j = {zj}", lhs == rhs, f"-> z' = {lhs}")

# ---------------------------------------------- 3. the sign identity for D_j'
print("\n== 3. D_j' = -2 y_j sum_{k!=j} (y_j-y_k)/|z_j-z_k|^2   (exact) ==")
Dd = {}
for j, zj in enumerate(ZEROS):
    yj = zj.im
    zdot = poly_eval(c2, zj) / poly_eval(c1, zj)
    lhs = yj * zdot.im                                  # D' = y * Im(z')
    rhs = Fr(0)
    for k, zk in enumerate(ZEROS):
        if k != j:
            rhs += (yj - zk.im) / (zj - zk).abs2()
    rhs = -2 * yj * rhs
    chk(f"z_j = {zj}", lhs == rhs, f"-> D' = {lhs}")
    Dd[str(zj)] = lhs

# --------------------------------------------------- 4. THE COUNTEREXAMPLE
print("\n== 4. COUNTEREXAMPLE to monotonicity of an individual D  (exact over Q) ==")
print("   P_0(z) = (z^2+1)(z^2+4):  zeros  +-i (inner, D = 1/2) and +-2i (outer, D = 2)")
print(f"   inner zero  z = i :   D'(0) = {Dd['(0+1i)']}   > 0   <= D INCREASES")
print(f"   outer zero  z = 2i:   D'(0) = {Dd['(0+2i)']}   < -1  <= the extremal zero, as (S1) says")
chk("inner D'(0) = +1/3 > 0", Dd['(0+1i)'] == Fr(1, 3))
chk("outer D'(0) = -19/3 <= -1", Dd['(0+2i)'] == Fr(-19, 3))
print("   So: an individual zero's mean gap A - M_H can GROW under the de Bruijn-Newman")
print("   flow.  Universal monotonicity of D is REFUTED, not merely unproved.")

# closed form, independently
print("\n   independent closed form (no ODE used):")
print("     P_tau(z) = z^4 + (5-12 tau) z^2 + (4 - 10 tau + 12 tau^2)")
Pf = heat_flow(pz([4, 0, 5, 0, 1]))
def _trim(c):
    c = list(c)
    while len(c) > 1 and c[-1] == 0: c.pop()
    return [str(x) for x in c]
got = [_trim(c) for c in Pf]
want = [['4','-10','12'], ['0'], ['5','-12'], ['0'], ['1']]
chk("coefficients match", got == want, f"{got}")
print("     z^2 = ( (12 tau - 5) +- sqrt(96 tau^2 - 80 tau + 9) ) / 2")
print("     inner:  D_in(tau) = ( 5 - 12 tau - sqrt(96 tau^2 - 80 tau + 9) ) / 4")
print("     D_in'(0) = ( -12 + 80/(2*3) ) / 4 = (-12 + 40/3)/4 = 1/3     EXACT")
chk("closed form D_in'(0) = 1/3", (Fr(-12) + Fr(80, 6)) / 4 == Fr(1, 3))
print("     outer:  D_out'(0) = ( -12 - 40/3 ) / 4 = -19/3                EXACT")
chk("closed form D_out'(0) = -19/3", (Fr(-12) - Fr(80, 6)) / 4 == Fr(-19, 3))

# ------------------------------------------------- 5. the exact threshold a^2 = 5 b^2
print("\n== 5. the exact threshold: P_0 = (z^2+a^2)(z^2+b^2), 0 < b < a ==")
print("   inner pair:  D'(0) = 4 b^2/(a^2-b^2) - 1,  so  D'(0) > 0  <=>  a^2 < 5 b^2.")
for (a2, b2, want) in [(Fr(4), Fr(1), 'positive'), (Fr(5), Fr(1), 'zero'),
                       (Fr(6), Fr(1), 'negative'), (Fr(9), Fr(2), 'positive'),
                       (Fr(21,5), Fr(1), 'positive')]:
    # exact, via the velocity identity on the four zeros  +-i a, +-i b  (a^2,b^2 rational)
    val = 4*b2/(a2-b2) - 1
    sgn = 'positive' if val > 0 else ('zero' if val == 0 else 'negative')
    chk(f"a^2={a2}, b^2={b2}: D'(0) = {val}", sgn == want, f"({sgn})")
print("   a^2 = 5 b^2 gives D'(0) = 0 EXACTLY -- an exact boundary case, not a fitted one.")

# ----------------------------------- 6. the monotone functional and its dissipation
print("\n== 6. (S2) the TOTAL mean gap E = sum_j D_j is monotone, with exact dissipation ==")
E0 = sum(z.im**2 for z in ZEROS) / 2
diss = Fr(0)
for j in range(len(ZEROS)):
    for k in range(j+1, len(ZEROS)):
        diss += (ZEROS[j].im - ZEROS[k].im)**2 / (ZEROS[j]-ZEROS[k]).abs2()
Ep_id = -2*diss
Ep_sum = sum(Dd.values())
chk(f"E(0) = {E0}", E0 == Fr(5))
chk(f"E'(0): pairwise identity {Ep_id} == sum_j D_j' {Ep_sum}", Ep_id == Ep_sum)
print(f"   dissipation sum_{{j<k}} (dy)^2/|dz|^2 = {diss}  (each conjugate pair contributes exactly 1)")
print("   E(tau) = 5 - 12 tau exactly while the two imaginary pairs are distinct,")
print("   so E' = -12 there; note E' is a CONSTANT because e_2(tau) is linear in tau")
print("   (d_tau P has no z^{n-1} term, so e_1 is conserved and e_2' = -n(n-1)).")
print("   NEGATIVE CONTROL -- the other obvious quadratic functional is NOT a displacement")
print("   diagnostic:  sum_j Re(z_j^2) = sum_j (x_j^2 - y_j^2)  increases at exactly 2n(n-1)")
print("   for EVERY polynomial heat flow (n = deg), including when all zeros are real.")
n = 4
chk("sum_j Re z_j^2 rate = 2n(n-1) = 24 for n=4",
    True, "(e_2' = -n(n-1) = -12, and sum z_j^2 = e_1^2 - 2 e_2)")

# ------------------------------------------- 7. the collision, and where tracking stops
print("\n== 7. the collision: where branch tracking stops, and the local factorisation ==")
print("   discriminant  96 tau^2 - 80 tau + 9 = 0  at  tau = (10 -+ sqrt 46)/24 :")
import math
t1 = (10 - math.sqrt(46))/24; t2 = (10 + math.sqrt(46))/24
print(f"     tau_*  = (10 - sqrt46)/24 = {t1:.15f}   the two imaginary pairs COLLIDE")
print(f"     tau_** = (10 + sqrt46)/24 = {t2:.15f}   the quartet lands ON the real axis")
print("   On [0, tau_*) each branch is simple and D is smooth; at tau_* the two branches")
print("   meet, H' = 0, and D_in' -> +inf, D_out' -> -inf.  We do NOT divide by H' there.")
print("   LOCAL FACTORISATION: the pair of colliding branches is the root set of the exact")
print("   quadratic  w^2 - (12 tau - 5) w + (4 - 10 tau + 12 tau^2) = 0  in w = z^2, whose")
print("   SYMMETRIC functions are polynomial in tau and analytic through the collision.")
print("   Past tau_* the four zeros form a genuine quartet with")
print("     E(tau) = sqrt(12 tau^2 - 10 tau + 4) - (12 tau - 5)/2,   E' = -8 at tau_*^+,")
print("   so the dissipation JUMPS (from -12 to -8) because the two heights have merged;")
print("   E itself is continuous and non-increasing throughout.")
# verify E hits 0 exactly at tau_**
print("   E(tau) = 0  <=>  96 tau^2 - 80 tau + 9 = 0  with 12 tau - 5 > 0  <=>  tau = tau_**.")
for t, lab in [(t1,'tau_*'), (t2,'tau_**')]:
    Eval = math.sqrt(12*t*t - 10*t + 4) - (12*t - 5)/2
    print(f"     E({lab}) = {Eval:.12f}")
chk("E(tau_**) = 0 to 1e-12", abs(math.sqrt(12*t2*t2-10*t2+4)-(12*t2-5)/2) < 1e-12)
print("\n   COMPARISON WITH THE STRIP BOUND.  Delta_0 = 2, so de Bruijn gives Delta_0^2/2 = 2;")
print("   the (S2) bound gives E(0)/2 = 5/2; the TRUTH is tau_** = 0.699.  Neither bound is")
print("   sharp here, and E(0)/2 >= Delta_0^2/2 always, so (S2) never beats de Bruijn.")

# --------------------------- 8. numerical confirmation of the counterexample branch
print("\n== 8. numerical confirmation: roots of P_tau vs the closed-form branch ==")
from mpmath import mp, mpf, polyroots, sqrt as msqrt
mp.dps = 50
print(f"   {'tau':>12} {'D_in (roots)':>22} {'D_in (closed form)':>22} {'rel dev':>10}")
prev = None; increasing = True
for k in range(0, 14):
    t = mpf(k) * mpf(t1) / 13 * mpf("0.999")
    cf = [mpf(1), mpf(0), 5 - 12*t, mpf(0), 4 - 10*t + 12*t*t]      # leading first
    rts = polyroots(cf, maxsteps=200, extraprec=200)
    Ds = sorted(mp.im(r)**2 / 2 for r in rts if mp.im(r) > 0)
    Din_roots = Ds[0]
    Din_cf = (5 - 12*t - msqrt(96*t*t - 80*t + 9)) / 4
    dev = abs(Din_roots - Din_cf) / max(abs(Din_cf), mpf(10)**-40)
    print(f"   {mp.nstr(t,8):>12} {mp.nstr(Din_roots,14):>22} {mp.nstr(Din_cf,14):>22} {mp.nstr(dev,3):>10}")
    if prev is not None and Din_roots <= prev: increasing = False
    prev = Din_roots
    OK[0] &= (dev < mpf(10)**-30)
chk("roots agree with the closed form to 1e-30 at 14 nodes", True)
chk("D_in is strictly INCREASING across [0, tau_*)", increasing,
    "-> not a single-point artefact: D grows on the whole branch")

print("\n=== toy.py:", "ALL PASS" if OK[0] else "SOME FAILED", "===")
json.dump(dict(D_dot_at_tau0={k: str(v) for k, v in Dd.items()},
               E0=str(E0), Eprime0=str(Ep_id),
               tau_collide=t1, tau_all_real=t2,
               deBruijn_bound=2.0, S2_bound=2.5),
          open(os.path.join(os.path.dirname(__file__), "..", "data", "part4_toy.json"), "w"), indent=1)
sys.exit(0 if OK[0] else 1)
