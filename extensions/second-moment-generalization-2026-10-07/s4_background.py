"""Section 4: A_a = G + T_a with a fixed Hermitian background G >= 0.

Claim under test (derived, then checked):
  T_a = T_0 + a^2 Q + O(a^4),  T_0 = 2|1><1|,  Q = kernel (x-y)^2.
  For v in E = ker(A_0) we have <1,v> = 0 and G v = 0, so
      <v, Q v> = <v,x^2><1,v> + <v,1><x^2,v> - 2|<x,v>|^2 = -2|<x,v>|^2,
  i.e.  P_E Q P_E = -2 |P_E x><P_E x| ,  hence
      lambda_min(A_a) = -2 a^2 ||P_E x||^2 + O(a^4)   when E != 0 and P_E x != 0.
"""
import sys; sys.path.insert(0,'.')
from kernel2 import *
from fractions import Fraction as Fr
getcontext().prec = 60

# ---- exact linear algebra over Q for the projections (no floating point here)
def gram_schmidt_Q(vs):
    out=[]
    for v in vs:
        w=list(v)
        for u in out:
            c=sum(a*b for a,b in zip(w,u))/sum(b*b for b in u)
            w=[a-c*b for a,b in zip(w,u)]
        if any(c!=0 for c in w): out.append(w)
    return out
def proj_perp_Q(x, basis):
    """exact projection of x onto the orthogonal complement of span(basis)"""
    w=list(x)
    for u in gram_schmidt_Q(basis):
        c=sum(a*b for a,b in zip(w,u))/sum(b*b for b in u)
        w=[a-c*b for a,b in zip(w,u)]
    return w
def norm2_Q(v): return sum(c*c for c in v)

def run(label, xs_f, G_rows_f, extra_span):
    """xs_f: Fractions nodes (weights 1 => Euclidean frame = coordinate frame).
       extra_span: list of exact vectors spanning (range of G), so ker G = that^perp."""
    n=len(xs_f)
    one=[Fr(1)]*n
    xv =list(xs_f)
    E_basis_killers = [one]+extra_span            # E = {1} ^perp  intersect  ker G
    PEx = proj_perp_Q(xv, E_basis_killers)
    pred = 2*norm2_Q(PEx)                          # |coefficient|
    dimE = n - len(gram_schmidt_Q(E_basis_killers))
    xs=[D(str(float(t))) if t.denominator>10**12 else D(t.numerator)/D(t.denominator) for t in xs_f]
    ws=[ONE]*n
    G=[[ (D(G_rows_f[i][j].numerator)/D(G_rows_f[i][j].denominator)) for j in range(n)] for i in range(n)]
    lam0 = lambda_min_jacobi(xs,ws,ZERO,G)
    print(f"\n  {label}")
    print(f"    dim E = {dimE},  ||P_E x||^2 = {norm2_Q(PEx)} , predicted |coeff| = {pred}")
    # NOTE: always print Decimals in scientific notation. Slicing str(Decimal) drops the
    # exponent and makes 1e-25 look like a number of order 1 -- that bug cost a wrong
    # reading of case (C) during this exploration.
    print(f"    lambda_min(A_0) = {float(lam0):.6e}")
    for aa in ('1e-2','1e-3','1e-4'):
        av=D(aa); lam=lambda_min_jacobi(xs,ws,av,G)
        q2=(lam-lam0)/av**2; q4=(lam-lam0)/av**4; q6=(lam-lam0)/av**6
        print(f"    a={aa:<5} lam_min={float(lam):>14.6e}  d/a^2={float(q2):>14.6e}  "
              f"d/a^4={float(q4):>12.4e}  d/a^6={float(q6):>12.4e}   predicted d/a^2 = {-float(pred):.6g}")

xs_f=[Fr(0),Fr(1),Fr(3),Fr(4)]              # 4 nodes, counting measure, mean 2
n=4
Z=[[Fr(0)]*n for _ in range(n)]
def outer(v,c=Fr(1)): return [[c*v[i]*v[j] for j in range(n)] for i in range(n)]

print("="*84); print("4.  background G >= 0 :  coefficient = 2 ||P_E x||^2  (projected second moment)"); print("="*84)
# (A) G = 0  -> E = 1^perp, ||P_E x||^2 = V_mu
run("(A) G = 0        -> full centred moment", xs_f, Z, [])
m=Fr(sum(xs_f),n); print(f"      check: V_mu = sum (x-mean)^2 = {sum((t-m)**2 for t in xs_f)}  (should equal ||P_E x||^2)")

# (B) G = 3 |v><v| with v generic -> E = {1,v}^perp, projected moment, strictly smaller
v=[Fr(1),Fr(0),Fr(-1),Fr(2)]
run("(B) G = 3|v><v|, v generic -> PROJECTED moment", xs_f, outer(v,Fr(3)), [v])

# (C) G = 3 |xc><xc| with xc = centred x  -> P_E x = 0, quadratic term DIES
xc=[t-m for t in xs_f]
run("(C) G = 3|x_c><x_c|      -> P_E x = 0, quadratic term vanishes", xs_f, outer(xc,Fr(3)), [xc])

# (D) G = 5 (I - |1><1|/M) -> E = {0}, A_0 strictly positive
G_D=[[ (Fr(5) if i==j else Fr(0)) - Fr(5,n) for j in range(n)] for i in range(n)]
run("(D) G = 5 P_{1^perp}     -> E = {0}, A_0 > 0", xs_f, G_D, [[Fr(1)]*n])
