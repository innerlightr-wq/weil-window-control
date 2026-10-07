#!/usr/bin/env python3
r"""Part 2: the mandatory controls (brief section 11) and the coordinate-response tests."""
from fractions import Fraction as Fr
import math, cmath

def e(x, n=7): return ('%.' + str(n) + 'e') % float(x)
L = 0.8; gam = 14.1347
def simpson(g,a,b,N=2000):
    if N%2: N+=1
    h=(b-a)/N; s=g(a)+g(b)
    for i in range(1,N): s+=g(a+i*h)*(4 if i%2 else 2)
    return s*h/3
def Phi(f, d, Lv=L, g=gam):      # exact quartet functional  4 Re F(g+id)^2 - 4F(g)^2
    def G(t):
        re=simpson(lambda u:(f(u)*cmath.exp(1j*t*u)).real,-Lv,Lv)
        im=simpson(lambda u:(f(u)*cmath.exp(1j*t*u)).imag,-Lv,Lv)
        return complex(re,im)
    Gd=G(complex(g,d)); F0=G(g).real
    return 4*(Gd*Gd).real - 4*F0*F0
TESTS={"cos-bump":lambda u: math.cos(math.pi*u/(2*L)),
       "extremal":lambda u: u*math.sin(gam*u),
       "cos(gam u)":lambda u: math.cos(gam*u)}

print("="*78); print("7.  REMAINDER CONTROL (control 5): same physical displacement"); print("="*78)
cL = 8*L**3*(1/3+1/math.sqrt(5))
def B_delta(d): return cL*d*d + (16/3)*L**5*d**4*math.exp(2*L*d)
print("  paper budget  B_L(delta) = c(L) delta^2 + (16/3) L^5 delta^4 e^{2L delta},")
print("  c(L) = 8L^3(1/3 + 1/sqrt5) = %s at L = %.1f" % (e(cL), L))
print("\n  delta    B(delta)        B via eta EXACT  B via R EXACT    B eta TRUNC(eta^2)  B R TRUNC(R-1)")
for d in (0.05, 0.1, 0.25, 0.375):
    eta = math.atanh(2*d); R = 1/(1-4*d*d)
    b   = B_delta(d)
    b_eta_exact = B_delta(0.5*math.tanh(eta))              # exact rewriting
    b_R_exact   = B_delta(0.5*math.sqrt(1-1/R))            # exact rewriting
    b_eta_trunc = cL*eta*eta/4 + (16/3)*L**5*(eta/2)**4*math.exp(L*eta)   # tanh -> eta
    b_R_trunc   = cL*(R-1)/4   + (16/3)*L**5*((R-1)/4)**2*math.exp(L*math.sqrt(R-1))
    print("  %-8.4f %-15s %-16s %-16s %-19s %s"
          % (d, e(b), e(b_eta_exact), e(b_R_exact), e(b_eta_trunc), e(b_R_trunc)))
print("\n  EXACT rewritings agree to machine precision -> NO gain and NO loss:")
print("   max dev over the table = %s"
      % e(max(abs(B_delta(d)-B_delta(0.5*math.tanh(math.atanh(2*d)))) for d in (0.05,0.1,0.25,0.375))))
print("  TRUNCATED forms are strictly LARGER (worse) at the same physical delta:")
for d in (0.25, 0.375):
    eta=math.atanh(2*d); R=1/(1-4*d*d)
    bt=cL*eta*eta/4+(16/3)*L**5*(eta/2)**4*math.exp(L*eta)
    bR=cL*(R-1)/4+(16/3)*L**5*((R-1)/4)**2*math.exp(L*math.sqrt(R-1))
    print("   delta=%.4f : eta-trunc / exact = %.5f    R-trunc / exact = %.5f"
          % (d, bt/B_delta(d), bR/B_delta(d)))
print("  because tanh^2(eta) <= eta^2 and (R-1)/R <= R-1 on the domain.")
print("  VERDICT on remainder: UNCHANGED in difficulty if rewritten exactly,")
print("  STRICTLY WORSE if truncated in eta or R.  No coordinate gives a simpler remainder.")

print("\n" + "="*78); print("8.  MONOTONICITY TRANSFERS; CONVEXITY DOES NOT"); print("="*78)
print("  on 0 <= delta < 1/2 the four coordinates are STRICTLY MONOTONE in each other:")
print("   delta^2    chi^2=4delta^2  eta=arctanh(2d)  R=1/(1-4d^2)   all increasing?")
prev=None; mono=True
for d in [i*0.05 for i in range(10)]:
    row=(d*d, 4*d*d, math.atanh(2*d), 1/(1-4*d*d))
    if prev: mono &= all(row[k]>prev[k] for k in range(4))
    prev=row
    print("   %-10.6f %-15.6f %-16.6f %-14.6f" % row)
print("   strictly increasing in lockstep :", mono)
print("  => ANY monotonicity statement in one coordinate is EQUIVALENT to it in the")
print("     others.  No monotonicity theorem can be gained by the coordinate change.")
print("     (This is the Two-Face manuscript's own point, section 7.5: strict monotonic")
print("      transformations give identical variance-explained under monotonic fits.)")
print("\n  CONVEXITY is the one property NOT preserved by a nonlinear monotone change.")
print("  Explicit witness: g(x) = x is linear (hence convex AND concave) in x = delta^2,")
print("  but as a function of R it is x = (R-1)/(4R), which is CONCAVE in R:")
for R in (1.0,1.2,1.5,2.0,3.0,5.0):
    print("     R=%-5.2f  (R-1)/(4R) = %-10.6f  second derivative d^2/dR^2 = %-12.6f (<0)"
          % (R,(R-1)/(4*R), -1/(2*R**3)))
print("  and as a function of eta it is tanh^2(eta)/4, CONVEX near 0 then CONCAVE:")
def t2(x): return math.tanh(x)**2/4
for et in (0.1,0.4,0.6,0.8,1.2):
    hh=1e-4; d2=(t2(et+hh)-2*t2(et)+t2(et-hh))/hh**2
    print("     eta=%-5.2f  tanh^2/4 = %-10.6f  d^2/deta^2 = %-12.6f  %s"
          % (et,t2(et),d2,"convex" if d2>0 else "CONCAVE"))
print("  inflection of tanh^2(eta)/4 at tanh^2 = 1/3, i.e. eta = arctanh(1/sqrt3) = %.7f,"
      % math.atanh(1/math.sqrt(3)))
print("  i.e. delta = %.7f -- INSIDE the physical domain and inside |delta| <= 3/8."
      % (math.tanh(math.atanh(1/math.sqrt(3)))/2))
print("  So a convexity claim is coordinate-DEPENDENT and therefore not a theorem about")
print("  the Weil form.  No convexity payoff is available this way.")
print("\n  exact Phi: is it monotone in |delta| on the physical range?  (f-dependent)")
for nm,f in TESTS.items():
    vals=[Phi(f,d) for d in (0.0,0.05,0.1,0.2,0.3,0.375)]
    signs=["+" if v>0 else ("0" if v==0 else "-") for v in vals]
    dec=all(vals[i+1]<=vals[i] for i in range(len(vals)-1))
    inc=all(vals[i+1]>=vals[i] for i in range(len(vals)-1))
    print("   %-11s Phi = %s" % (nm, " ".join(e(v,3) for v in vals)))
    print("   %-11s signs %s  monotone decreasing: %s  increasing: %s"
          % ("", "".join(signs), dec, inc))
print("  Phi has NO universal sign and NO universal monotonicity across f -- consistent")
print("  with A_f being sign-indefinite (checks.py section 6).")

print("\n" + "="*78); print("9.  FAKE-GENERATOR CONTROL (control 2)"); print("="*78)
print("  Any even analytic function of a parameter is a power series in that parameter")
print("  squared.  So an even Taylor series is NOT evidence for a scaling generator.")
def fake(x): return math.cos(x) + math.exp(-x*x)        # even; nothing to do with scaling
co=[2.0, -1.0-1.0, 1/24.0+0.5]      # cos: 1 - x^2/2 + x^4/24 ; exp(-x^2): 1 - x^2 + x^4/2
print("   witness g(x) = cos(x) + exp(-x^2):  g(x) = %.4f %+.4f x^2 %+.4f x^4 + ..."
      % (co[0],co[1],co[2]))
for x in (0.1,0.3,0.5):
    series=co[0]+co[1]*x*x+co[2]*x**4
    print("     x=%.2f  g = %-12s even series = %-12s dev = %s"
          % (x, e(fake(x)), e(series), e(abs(fake(x)-series))))
print("   g is even and analytic, so g(x) = sum a_k x^{2k}, and one can ALWAYS write")
print("   g(x) 'formally' as 2I + x^2 A^2 + ... by DEFINING A^2 := a_1 I, A^4 := 12 a_2 I, ...")
print("   but those definitions need not come from one operator: consistency requires")
print("   a_2 = a_1^2/12 for a single A with A^4 = (A^2)^2.  Check the witness:")
print("     a_1 = %.6f  a_1^2/12 = %.6f  a_2 = %.6f  consistent? %s"
      % (co[1], co[1]**2/12, co[2], abs(co[1]**2/12-co[2])<1e-12))
print("   -> NOT consistent, so no single generator reproduces it.  Taylor matching alone")
print("      proves nothing.  CONTRAST: the real quartet case passes this test exactly,")
print("      because F(gamma+i delta) = Fhat(e^{-delta A} f)(gamma) holds at ALL orders")
print("      (checks.py section 3) -- an identity, not a coefficient match.")

print("\n" + "="*78); print("10. COMPOSITION CONTROL (control 3)"); print("="*78)
print("  A system where eta = arctanh(2 delta) is valid algebra but has NO native")
print("  additive dynamics: a single Bernoulli trial with p = 1/2 + delta.")
print("  There is no binary operation on two independent trials that composes their")
print("  p's by Moebius addition; the natural operations are:")
for (p1,p2) in [(0.6,0.7),(0.75,0.55)]:
    d1,d2=p1-0.5,p2-0.5
    print("   p1=%.2f p2=%.2f | product (AND) p1 p2 = %.4f -> delta = %+.4f"
          % (p1,p2,p1*p2,p1*p2-0.5))
    print("                  | mixture (1/2,1/2)  = %.4f -> delta = %+.4f  (= (d1+d2)/2)"
          % ((p1+p2)/2,(p1+p2)/2-0.5))
    print("                  | Moebius d1(+)d2    = %+.4f  <- matches NEITHER"
          % ((d1+d2)/(1+4*d1*d2)))
print("  So the hyperbolic coordinate is exact as algebra and empty as dynamics here.")
print("  This is the manuscript's own 'exact translation, no shared dynamics' case.")

print("\n" + "="*78); print("11. COORDINATE-INVARIANCE CONTROL (control 1)"); print("="*78)
print("  sign of Phi must not change under delta -> eta -> R at the same physical point")
bad=0
for nm,f in TESTS.items():
    for d in (0.05,0.2,0.375):
        eta=math.atanh(2*d); R=1/(1-4*d*d)
        a=Phi(f,d); b=Phi(f,0.5*math.tanh(eta)); c=Phi(f,0.5*math.sqrt(1-1/R))
        same=(a>0)==(b>0)==(c>0) and abs(a-b)<1e-12*max(1,abs(a)) and abs(a-c)<1e-10*max(1,abs(a))
        bad += 0 if same else 1
        print("   %-11s delta=%-7.4f Phi_delta=%-13s Phi_eta=%-13s Phi_R=%-13s same=%s"
              % (nm,d,e(a,4),e(b,4),e(c,4),same))
print("   all sign statements invariant :", bad==0)

print("\n" + "="*78); print("12. ANALYTICITY: delta is the best coordinate"); print("="*78)
print("  Phi is EVEN in delta (checks.py section 3), so Phi = Psi(delta^2) with Psi ENTIRE")
print("  (f compactly supported => Fhat entire of exponential type L).")
print("   coordinate   substitution for delta^2      singularities of the composition")
print("   delta        delta^2                       NONE (entire)")
print("   chi          chi^2/4                       NONE (entire)")
print("   eta          tanh^2(eta)/4                 POLES at eta = +- i pi/2 (+ i k pi)")
print("                                              -> finite radius pi/2 = %.7f" % (math.pi/2))
print("   R            (R-1)/(4R)                    pole at R = 0 only; ANALYTIC at R = 1")
print("  NOTE, against the obvious error: the R-substitution is RATIONAL, and because Psi")
print("  is a function of delta^2 the square root never appears.  There is NO branch point")
print("  at R = 1.  (delta = (1/2)sqrt(1-1/R) alone does have one; Phi does not.)")
print("  CONTROL: an ODD function of delta WOULD acquire a branch point in R.")
print("   witness delta itself:  delta(R) = (1/2)sqrt(1-1/R), two-valued around R = 1:")
for R in (1.0001,1.01,1.1):
    print("     R=%-8.4f  +branch = %+.6f   -branch = %+.6f" % (R, 0.5*math.sqrt(1-1/R), -0.5*math.sqrt(1-1/R)))
print("  Phi is even, so it is single-valued -- the branch point is invisible to it.")
print("\nDONE part 2")
