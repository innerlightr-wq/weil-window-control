"""Section 6: why 'symmetric spectral problem' alone is not enough."""
import sys, math; sys.path.insert(0,'.')
from kernel2 import *
getcontext().prec=40
print("="*84); print("6(a)  spectral symmetry is NOT A(a)=A(-a); ordered min is not an analytic branch"); print("="*84)
print("  A(a) = diag(a, -a).  spec A(a) = {a,-a} = spec A(-a) as a SET, but A(a) != A(-a).")
print(f"  {'a':>8} {'eigs':>22} {'ordered min':>14}   analytic branches a and -a are smooth;")
for a in (-0.3,-0.1,0,0.1,0.3):
    print(f"  {a:>8} {str((a,-a)):>22} {min(a,-a):>14}   their MINIMUM is -|a|: not differentiable at 0.")
print("  => the ordered minimum is only Lipschitz; a quadratic law needs the branch, not min.")

print()
print("="*84); print("6(b)  even analytic family whose minimum INCREASES quadratically"); print("="*84)
print("  A(a) = diag(1+a^2, 2):  A(a)=A(-a), lambda_min = 1 + a^2 > lambda_min(0). Coefficient +1.")

print()
print("="*84); print("6(c)  even analytic, vanishing quadratic term, higher-order leading response"); print("="*84)
print("  Trivial: A(a) = diag(-a^4, 1) -> lambda_min = -a^4.")
print("  Non-trivial, inside our own family: section 4 case (C), where P_E x = 0.")
print("  Measured there: (lam-lam0)/a^6 -> -0.8 exactly, i.e. NO a^2 and NO a^4 term.")
print("  Reason (computed): the first-order a^4 term of T_a contributes +(1/2)|<x^2,v>|^2 a^4")
print("  and second-order perturbation theory in a^2 Q contributes -(1/2)|<x^2,v>|^2 a^4;")
print("  they cancel, leaving a^6.")

print()
print("="*84); print("6(d)  strictly positive baseline: small perturbations keep positivity"); print("="*84)
print("  Section 4 case (D): lambda_min(A_0) = 5 > 0, response -20a^2, so lambda_min(A_a) > 0")
print("  for |a| < sqrt(5/20) = 0.5. A negative eigenvalue is CREATED only once the")
print("  perturbation exceeds the baseline gap. 'Expansion about a positive baseline' and")
print("  'creation of a negative eigenvalue' are different statements.")

print()
print("="*84); print("6(e)  same support length, different mass and variance"); print("="*84)
def coeff(xs,ws): return 2*centred_moment(xs,ws)
cases=[("Lebesgue [-1,1], N=200",      *( (lambda L,N: ([-L+2*L/N*(D(i)+D('0.5')) for i in range(N)],[2*L/N]*N))(D(1),200) )),
       ("same support, mass x3",       *( (lambda L,N: ([-L+2*L/N*(D(i)+D('0.5')) for i in range(N)],[3*2*L/N]*N))(D(1),200) )),
       ("two atoms at +-1, mass 2",    [D(-1),D(1)], [ONE,ONE]),
       ("two atoms at +-1, mass 20",   [D(-1),D(1)], [D(10),D(10)])]
for name,xs,ws in cases:
    print(f"  {name:<30} M={float(mass(xs,ws)):>7.3f}  V_mu={float(centred_moment(xs,ws)):>9.5f}  "
          f"coeff={float(coeff(xs,ws)):>9.5f}")
print("  Support length alone fixes neither the coefficient nor its L-power.")

print()
print("="*84); print("6(f)  reparametrisation changes the apparent exponents"); print("="*84)
print("  a = c*delta  =>  coefficient multiplied by c^2.  In the note a = delta*log q, so the")
print("  function-field coefficient in delta is (log q)^2 times the one in a.")
print("  a = L*delta  =>  -2V_mu a^2 = -2V_mu L^2 delta^2, which ADDS 2 to the apparent window")
print("  power. Exponents are only comparable after the physical parameters are fixed.")
for L in (1,2,4):
    V=2*L**3/3
    print(f"    L={L}: coeff in a = {2*V:.4f} (~L^3);  coeff in delta with a=L*delta = {2*V*L*L:.4f} (~L^5)")
