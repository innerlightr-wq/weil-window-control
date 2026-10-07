#!/usr/bin/env python3
r"""RH partition-coordinate / intrinsic-dynamics audit: minimal checks.

Standard library only (mpmath is not installable on this box).  Exact rational
arithmetic where a sign or an identity is at stake; Simpson quadrature with a
stated node count elsewhere.  Every section carries a deliberately failing
control.  NOTE ON NOTATION: gamma below is ALWAYS the Weil spectral ordinate.
The Two-Face manuscript's "manifold Lorentz factor" gamma = 1/(2h) is written
LF here to avoid the collision.
"""
from fractions import Fraction as Fr
import math, cmath

def e(x, n=7): return ('%.' + str(n) + 'e') % float(x)

# ----------------------------------------------------------------- 1. dictionary
print("="*78); print("1.  The partition dictionary -- exact identities"); print("="*78)
def coords(d):
    """d = delta (Fraction or float).  Returns chi, h, R, eta, LF."""
    chi = 2*abs(d)
    h2  = (1 - 4*d*d)/4          # h^2 = ab = (1/4)(1-4 delta^2)
    R   = 1/(1 - 4*d*d)
    return chi, h2, R
print("  delta      a=1/2+d    b=1/2-d    a+b   chi=2|d|  4h^2+chi^2   R(1-4d^2)")
for d in [Fr(0), Fr(1,10), Fr(1,4), Fr(3,8), Fr(-1,5), Fr(49,100)]:
    a, b = Fr(1,2)+d, Fr(1,2)-d
    chi, h2, R = coords(d)
    print("  %-10s %-10s %-10s %-5s %-9s %-12s %s"
          % (d, a, b, a+b, chi, 4*h2 + chi*chi, R*(1-4*d*d)))
    assert a+b == 1 and h2 == a*b and 4*h2 + chi*chi == 1 and R*(1-4*d*d) == 1
print("  EXACT over Q: a+b=1, h^2=ab, 4h^2+chi^2=1, R=1/(1-4delta^2)=1/(4ab)  ALL PASS")

print("\n  hyperbolic leg (floats; eta = arctanh(2 delta) signed, xi = |eta|)")
print("   delta      eta           R         cosh^2(eta)   R-1        sinh^2(eta)  (R-1)/(4R) vs delta^2")
for d in [0.0, 0.1, 0.25, 0.375, -0.2, 0.49]:
    eta = math.atanh(2*d); R = 1/(1-4*d*d)
    print("   %-10.4f %-13.7f %-9.6f %-13.6f %-10.6f %-12.6f %s"
          % (d, eta, R, math.cosh(eta)**2, R-1, math.sinh(eta)**2,
             "%.3e" % abs((R-1)/(4*R) - d*d)))
    assert abs(math.cosh(eta)**2 - R) < 1e-12
    assert abs(math.sinh(eta)**2 - (R-1)) < 1e-12
    assert abs(0.25*math.tanh(eta)**2 - d*d) < 1e-14
    assert abs(math.tanh(eta) - 2*d) < 1e-14
print("  R=cosh^2(eta), R-1=sinh^2(eta), delta^2=(1/4)tanh^2(eta)=(R-1)/(4R), 2delta=tanh(eta)  PASS")
print("\n  RH point:  delta=0 <=> chi=0 <=> eta=0 <=> R=1 <=> h=1/2 <=> a=b=1/2 :",
      all([coords(Fr(0))[0]==0, math.atanh(0)==0, coords(Fr(0))[2]==1, coords(Fr(0))[1]==Fr(1,4)]))

print("\n  DOMAIN CONTROL (control 6).  |delta|<1/2 ; quasi-RH adds |delta|<=3/8.")
for lab, dm in [("|delta| < 1/2 (functional equation / beta in (0,1))", 0.5),
                ("|delta| <= 3/8 (quasi-RH zero-free Re s > 7/8)", 0.375)]:
    print("   %-52s R in [1, %s]   |eta| <= %s"
          % (lab, "inf" if dm==0.5 else "%.7f" % (1/(1-4*dm*dm)),
             "inf" if dm==0.5 else "%.7f" % math.atanh(2*dm)))
print("   these are ranges in delta; the certified Weil window L <= 0.8 is a range in L")
print("   -- a DIFFERENT variable.  Never mixed below.")

# ----------------------------------------------------------------- 2. calculus identity
print("\n" + "="*78); print("2.  (1/2) d^2/dgamma^2 F^2 = F'^2 + F F''   [elementary; NOT new]"); print("="*78)
# symbolic over Q with a polynomial stand-in for F
import itertools
def pmul(p,q):
    r=[Fr(0)]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q): r[i+j]+=a*b
    return r
def pder(p): return [p[k]*k for k in range(1,len(p))] or [Fr(0)]
P=[Fr(2),Fr(-3),Fr(5),Fr(1),Fr(-7)]          # arbitrary F
lhs=pder(pder(pmul(P,P))); lhs=[x/2 for x in lhs]
rhs_=pmul(pder(P),pder(P))
rhs2=pmul(P,pder(pder(P)))
n=max(len(rhs_),len(rhs2)); rhs=[ (rhs_[i] if i<len(rhs_) else Fr(0))+(rhs2[i] if i<len(rhs2) else Fr(0)) for i in range(n)]
m=max(len(lhs),len(rhs))
ok=all((lhs[i] if i<len(lhs) else Fr(0))==(rhs[i] if i<len(rhs) else Fr(0)) for i in range(m))
print("  exact over Q for an arbitrary quartic F :", ok)
print("  => Delta Q = -2 delta^2 (d^2/dgamma^2) F(gamma)^2 + O(delta^4)   [restatement only]")
print("  CONTROL: the WRONG identity (1/2)(F^2)'' = F'^2 - F F'' fails :",
      not all((lhs[i] if i<len(lhs) else Fr(0))==((rhs_[i] if i<len(rhs_) else Fr(0))-(rhs2[i] if i<len(rhs2) else Fr(0))) for i in range(m)))

# ----------------------------------------------------------------- 3. the EXACT generator
print("\n" + "="*78)
print("3.  THE EXACT GENERATOR OF OFF-LINE DISPLACEMENT")
print("="*78)
L = 0.8
def simpson(g, a, b, N=2000):
    if N % 2: N += 1
    hh=(b-a)/N; s=g(a)+g(b)
    for i in range(1,N): s += g(a+i*hh)*(4 if i%2 else 2)
    return s*hh/3
# paper convention (notes/proofs.md B.1/B.2):  F(t) = int f(u) e^{i t u} du
def Fhat(f, t, L=L):
    re = simpson(lambda u: f(u)*cmath.exp(1j*t*u).real, -L, L)
    im = simpson(lambda u: f(u)*cmath.exp(1j*t*u).imag, -L, L)
    return complex(re, im)
def Fhat_c(f, t, L=L):      # complex t
    re = simpson(lambda u: (f(u)*cmath.exp(1j*t*u)).real, -L, L)
    im = simpson(lambda u: (f(u)*cmath.exp(1j*t*u)).imag, -L, L)
    return complex(re, im)
TESTS = {"cos-bump": lambda u: math.cos(math.pi*u/(2*L)),
         "tent":     lambda u: 1-abs(u)/L,
         "extremal": lambda u: u*math.sin(14.1347*u)}
gam = 14.1347
print("  claim:  F(gamma + i delta) = Fhat( e^{-delta A} f )(gamma),  A = multiplication by u")
print("  f          delta     |F(gamma+i delta) - Fhat(e^{-delta A} f)(gamma)|")
for name, f in TESTS.items():
    for d in (0.05, 0.2, 0.375):
        lhs_ = Fhat_c(f, complex(gam, d))
        rhs_ = Fhat(lambda u: math.exp(-d*u)*f(u), gam)
        print("  %-10s %-9.3f %s" % (name, d, e(abs(lhs_-rhs_))))
print("  => the imaginary ordinate shift IS exactly the semigroup e^{-delta A}, A = u .")
print("     The natural parameter is DELTA.  No transform of delta is involved.")

print("\n  exact quartet functional (notes/proofs.md B.1), no remainder:")
print("     Phi(delta) := Q_delta - Q_0 = 4 Re F(gamma+i delta)^2 - 4 F(gamma)^2")
print("  f          delta     Phi exact        -4d^2(F'^2+FF'')   ratio")
for name, f in TESTS.items():
    F0 = Fhat(f, gam).real
    Fp = -simpson(lambda u: u*f(u)*math.sin(gam*u), -L, L)
    Fpp= -simpson(lambda u: u*u*f(u)*math.cos(gam*u), -L, L)
    Af = Fp*Fp + F0*Fpp
    for d in (0.01, 0.1, 0.375):
        G = Fhat_c(f, complex(gam, d))
        Phi = 4*(G*G).real - 4*F0*F0
        lead = -4*d*d*Af
        print("  %-10s %-9.3f %-16s %-18s %s"
              % (name, d, e(Phi), e(lead), e(Phi/lead) if lead else "n/a"))
print("  (ratio -> 1 as delta -> 0 : the leading law; Phi itself is the exact object)")

print("\n  evenness:  Re F(gamma+i delta)^2 = Re F(gamma-i delta)^2  (so Phi = Psi(delta^2))")
for name, f in TESTS.items():
    d = 0.3
    gp = Fhat_c(f, complex(gam, d)); gm = Fhat_c(f, complex(gam, -d))
    print("   %-10s |Re G(+d)^2 - Re G(-d)^2| = %s" % (name, e(abs((gp*gp).real-(gm*gm).real))))

# ----------------------------------------------------------------- 4. composition law
print("\n" + "="*78)
print("4.  COMPOSITION LAW -- the decisive falsification gate (brief section 8)")
print("="*78)
print("  native law:  e^{-d1 A} e^{-d2 A} = e^{-(d1+d2) A}   ->  PLAIN ADDITION in delta")
f = TESTS["cos-bump"]
for (d1,d2) in [(0.1,0.15),(0.2,0.3),(0.375,0.375)]:
    two  = Fhat(lambda u: math.exp(-d2*u)*math.exp(-d1*u)*f(u), gam)
    one  = Fhat(lambda u: math.exp(-(d1+d2)*u)*f(u), gam)
    mob  = (d1+d2)/(1+4*d1*d2)
    oneM = Fhat(lambda u: math.exp(-mob*u)*f(u), gam)
    print("   d1=%.3f d2=%.3f | additive d1+d2=%.4f  dev=%s | Moebius d1(+)d2=%.4f  dev=%s"
          % (d1,d2,d1+d2,e(abs(two-one)),mob,e(abs(two-oneM))))
print("  => composition is ADDITIVE in delta and is NOT Moebius addition.")
print("     Since 2 delta = tanh(eta), additivity in eta would REQUIRE Moebius addition.")
print("     So eta is additive for a law the RH side does not have.")
print("\n  DOMAIN-EXIT CONTROL: additive composition LEAVES the partition domain;")
print("  Moebius composition cannot.")
for (d1,d2) in [(0.3,0.3),(0.375,0.375),(0.45,0.4)]:
    print("   d1+d2 = %.4f  (>1/2 ? %s)   d1(+)d2 = %.4f  (>1/2 ? %s)"
          % (d1+d2, d1+d2>0.5, (d1+d2)/(1+4*d1*d2), (d1+d2)/(1+4*d1*d2)>0.5))
print("  'Broken domains' in the manuscript's own sense (its section 6.6).")

# ----------------------------------------------------------------- 5. the two generators
print("\n" + "="*78)
print("5.  THE TWO GENERATORS ARE CANONICALLY CONJUGATE, NOT EQUAL (brief section 9)")
print("="*78)
print("  Sonin/scaling (CC2021): vartheta(lambda) xi(v) = lambda^{-1/2} xi(lambda^{-1} v)")
print("    in u = log v this is TRANSLATION u -> u + t ; generator D = d/du (skew-adjoint)")
print("  off-line displacement: multiplication by e^{-delta u} ; generator A = u (self-adjoint)")
print("  commutator, EXACT over Q on polynomials:  (A D - D A) p = u p\' - (u p)\' = -p")
def polymul_u(p): return [Fr(0)] + list(p)            # multiply by u
def polyder(p):   return [p[k]*k for k in range(1,len(p))] or [Fr(0)]
def trim(p):
    q=list(p)
    while len(q)>1 and q[-1]==0: q.pop()
    return q
ok_all=True
for p in ([Fr(1)], [Fr(0),Fr(1)], [Fr(2),Fr(-3),Fr(5),Fr(1),Fr(-7)]):
    AD = polymul_u(polyder(p))          # u p'
    DA = polyder(polymul_u(p))          # (u p)'
    com = trim([ (AD[i] if i<len(AD) else Fr(0)) - (DA[i] if i<len(DA) else Fr(0))
                 for i in range(max(len(AD),len(DA))) ])
    want = trim([-c for c in p])
    ok = com == want; ok_all &= ok
    print("     p =", p, " -> [A,D]p =", com, " ; -p =", want, " ; equal:", ok)
print("     [A, D] = -I exactly :", ok_all)
print("  CONTROL: [A,A] = 0 and [D,D] = 0 (a commuting pair gives 0); and the WRONG")
print("  sign convention [D,A] = +I, not -I :",
      trim([ (polyder(polymul_u([Fr(1),Fr(2)]))[i] if i<len(polyder(polymul_u([Fr(1),Fr(2)]))) else Fr(0))
             - (polymul_u(polyder([Fr(1),Fr(2)]))[i] if i<len(polymul_u(polyder([Fr(1),Fr(2)]))) else Fr(0))
             for i in range(3)]) == trim([Fr(1),Fr(2)]))
print("  A and D generate a Heisenberg pair: non-commuting, dual under Fourier/Mellin.")
print("  Therefore t = eta (or t = delta) is NOT justified: the parameters act along")
print("  CONJUGATE directions.  Off-line displacement = IMAGINARY shift of the ordinate;")
print("  Sonin scaling = REAL shift of the physical coordinate.")

# ----------------------------------------------------------------- 6. projected generator norm
print("\n" + "="*78)
print("6.  K(L,gamma) IS a projected generator norm (brief section 6)"); print("="*78)
def K_quad(Lv, g): return simpson(lambda u: u*u*math.sin(g*u)**2, -Lv, Lv)
def K_closed(Lv, g):
    return (Lv**3)/3 - (Lv*Lv*math.sin(2*g*Lv)/(2*g) + Lv*math.cos(2*g*Lv)/(2*g*g)
                        - math.sin(2*g*Lv)/(4*g**3))
print("   L     gamma      K quadrature     K closed form (paper)   dev")
for Lv in (0.5, 0.8, 1.0):
    for g in (3.0, 14.1347):
        print("   %-5.2f %-10.4f %-16s %-23s %s"
              % (Lv, g, e(K_quad(Lv,g)), e(K_closed(Lv,g)), e(abs(K_quad(Lv,g)-K_closed(Lv,g)))))
print("  K(L,gamma) = int_{-L}^{L} u^2 sin^2(gamma u) du = || A P_L sigma_gamma ||^2,")
print("     sigma_gamma(u) = sin(gamma u),  P_L = multiplication by 1_[-L,L],  A = u .")
print("  The projection is NOT optional: ||A sigma_gamma||^2 over all of R diverges.")
for Lv in (10, 100, 1000):
    print("     partial ||A sigma||^2 on [-%d,%d] = %s" % (Lv, Lv, e(K_quad(Lv, 3.0), 4)))
print("  ||A|| on L^2[-L,L] = L, so K <= L^2 * ||P_L sigma||^2; the paper's 2L^3/3 is SHARPER:")
for Lv in (0.5, 0.8, 1.0):
    crude = Lv*Lv*simpson(lambda u: math.sin(3.0*u)**2, -Lv, Lv)
    print("     L=%.2f  crude L^2||P sigma||^2 = %-12s paper 2L^3/3 = %-12s K = %s"
          % (Lv, e(crude), e(2*Lv**3/3), e(K_quad(Lv,3.0))))
print("  NOTE: K = ||A P_L sigma||^2 is a RELABELLING of the paper's own extremal value,")
print("  not new mathematics.  And A_f = F'^2 + F F'' is NOT a norm: it is sign-indefinite.")
for name, f in TESTS.items():
    F0 = Fhat(f, gam).real
    Fp = -simpson(lambda u: u*f(u)*math.sin(gam*u), -L, L)
    Fpp= -simpson(lambda u: u*u*f(u)*math.cos(gam*u), -L, L)
    print("     %-10s F'^2 = %-14s F F'' = %-15s A_f = %s"
          % (name, e(Fp*Fp), e(F0*Fpp), e(Fp*Fp+F0*Fpp)))
print("  CONTROL (projection control 4): a projected generator norm is >= 0.  A_f is NOT:")
for nm, ff in [("cos(gamma u)", lambda u: math.cos(gam*u)),
               ("cos(gamma u)^2", lambda u: math.cos(gam*u)**2),
               ("1 (flat)", lambda u: 1.0)]:
    F0 = Fhat(ff, gam).real
    Fp = -simpson(lambda u: u*ff(u)*math.sin(gam*u), -L, L)
    Fpp= -simpson(lambda u: u*u*ff(u)*math.cos(gam*u), -L, L)
    print("     %-15s F'^2 = %-14s F F'' = %-15s A_f = %-15s sign %s"
          % (nm, e(Fp*Fp), e(F0*Fpp), e(Fp*Fp+F0*Fpp), "NEGATIVE" if Fp*Fp+F0*Fpp<0 else "positive"))
print("  => A_f takes BOTH signs, so it cannot be any projected generator norm ||P A psi||^2.")
print("     (Consistent with the paper: Theorem 8 completes the square precisely because")
print("      F F'' may be adverse, and section B.2 bounds |F'^2 + F F''| in absolute value.)")
print("\nDONE part 1")
