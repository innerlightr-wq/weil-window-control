"""Section 2: the finite-rank baseline, proved then calibrated."""
import sys; sys.path.insert(0, '.')
from kernel2 import *
from fractions import Fraction as Fr

getcontext().prec = 80
def dd(x): return D(str(x))

print("="*80); print("2.1  rank<=2 closed form vs an independent Jacobi eigensolve"); print("="*80)
cases = [
 ("counting {0..4}",        [dd(i) for i in range(5)],        [ONE]*5),
 ("irregular nodes+weights",[dd(v) for v in (0,'0.3','1.7',5,'5.1')], [dd(v) for v in ('0.4',2,'1.25','0.1',3)]),
 ("translated interval",    [dd(v) for v in (10,'10.5',11,12)],[ONE]*4),
 ("one point (mass 3)",     [dd(7)],                          [D(3)]),
 ("two equal points",       [dd(-1), dd(1)],                  [ONE,ONE]),
]
a = D('0.001')
for name, xs, ws in cases:
    lp, lm = lambda_pm_closed(xs, ws, a)
    jac = jacobi_eigenvalues(T_matrix(xs, ws, a))
    agree_min = abs(lm - jac[0]) < D(10)**-50
    agree_max = abs(lp - jac[-1]) < D(10)**-50
    nz = sum(1 for e in jac if abs(e) > D(10)**-45)
    print(f"  {name:<26} rank={nz}  lam_- closed={str(lm)[:16]:<16} jacobi={str(jac[0])[:16]:<16} "
          f"{'OK' if agree_min and agree_max else 'MISMATCH'}")

print()
print("="*80); print("2.2  lambda_-(a) = -2 V_mu a^2 - (2/3) V_4 a^4 + O(a^6)"); print("="*80)
def V4(xs, ws):
    m = mean(xs, ws); return sum((w*(x-m)**4 for x, w in zip(xs, ws)), ZERO)
for name, xs, ws in cases:
    V, V4_ = centred_moment(xs, ws), V4(xs, ws)
    print(f"\n  {name}:  M={str(mass(xs,ws))[:10]}  mean={str(mean(xs,ws))[:10]}  V_mu={str(V)[:12]}")
    for aa in ('1e-2','1e-3','1e-4'):
        av = D(aa); lm = lambda_pm_closed(xs, ws, av)[1]
        pred2 = -2*V*av**2
        pred4 = pred2 - (D(2)/3)*V4_*av**4
        print(f"    a={aa:<5} lam_-={str(lm)[:22]:<22} -2V a^2={str(pred2)[:22]:<22} "
              f"resid/a^4={str((lm-pred2)/av**4)[:14]:<14} pred -2V4/3={str(-(D(2)/3)*V4_)[:12]}")

print()
print("="*80); print("2.3  calibration (a): counting measure on {0..R}  ->  binom(R+2,3)"); print("="*80)
for R in (1,2,4,8,16,32):
    xs=[dd(i) for i in range(R+1)]; ws=[ONE]*(R+1)
    V = centred_moment(xs, ws)
    exact_V = Fr(R*(R+1)*(R+2),12)
    coeff = 2*V                                     # |coefficient| = 2 V_mu
    binom = Fr((R+2)*(R+1)*R,6)
    print(f"  R={R:<3} V_mu={str(V)[:14]:<14} exact R(R+1)(R+2)/12={str(exact_V)!s:<12} "
          f"2V_mu={str(coeff)[:12]:<12} binom(R+2,3)={binom}  match={abs(coeff-D(binom.numerator)/D(binom.denominator))<D(10)**-40}")
print("  => coefficient is binom(R+2,3) = R(R+1)(R+2)/6, i.e. the /6 form. NOT /3.")

print()
print("="*80); print("2.4  calibration (b): Lebesgue on [-L,L] -> exact 2L - sinh(2aL)/a"); print("="*80)
def lam_lebesgue_closed(L, a):
    if a == 0: return ZERO
    return 2*L - dsinh(2*a*L)/a
def lam_lebesgue_from_formula(L, a):
    """M - sqrt(pq) with p=q=sinh(2aL)/a for Lebesgue on [-L,L]."""
    p = dsinh(2*a*L)/a
    return 2*L - (p*p).sqrt()
for L in ('0.8','1.6','3.0'):
    Lv = D(L)
    for aa in ('1e-2','1e-3'):
        av = D(aa)
        e1, e2 = lam_lebesgue_closed(Lv,av), lam_lebesgue_from_formula(Lv,av)
        pred = -(D(4)/3)*Lv**3*av**2 - (D(4)/15)*Lv**5*av**4
        print(f"  L={L:<4} a={aa:<5} 2L-sinh(2aL)/a={str(e1)[:22]:<22} M-sqrt(pq)={str(e2)[:22]:<22} "
              f"-4L^3a^2/3-4L^5a^4/15={str(pred)[:22]:<22} agree={abs(e1-e2)<D(10)**-50}")
print("  continuous extension at a=0: sinh(2aL)/a -> 2L, so lambda_-(0)=0.")
print(f"  V_mu for Lebesgue[-L,L] = 2L^3/3, so -2V_mu = -4L^3/3  (matches the a^2 term).")

print()
print("="*80); print("2.5  when is lambda_- < 0 ?  (Cauchy-Schwarz, sharp)"); print("="*80)
print("  sqrt(<u,u><w,w>) >= |<u,w>| = M with equality iff u || w in L2(mu),")
print("  i.e. iff e^{2ax} is mu-a.e. constant, i.e. a = 0 or mu is a single atom.")
print("  => lambda_- < 0  strictly  <=>  a != 0 AND supp(mu) has at least two points.")
xs,ws=[dd(7)],[D(3)]
print(f"  one-point check: lambda_pm = {[str(v)[:8] for v in lambda_pm_closed(xs,ws,D('0.5'))]}  (rank 1, no negative eigenvalue)")

print()
print("="*80); print("2.6  explicit remainder bound:  lambda_-(a) >= -2 V_mu a^2 cosh(2 a R_mu)"); print("="*80)
print("  R_mu = max |x - mean| on supp mu  (<= support diameter).")
for name, xs, ws in cases:
    if len(xs) == 1: continue
    m = mean(xs,ws); R = max(abs(x-m) for x in xs); V = centred_moment(xs,ws)
    ok = True
    for aa in ('1e-1','1e-2','0.5'):
        av=D(aa); lm = lambda_pm_closed(xs,ws,av)[1]; bnd = -2*V*av**2*dcosh(2*av*R)
        ok &= (lm >= bnd)
    print(f"  {name:<26} R_mu={str(R)[:8]:<8} bound holds at a=0.1,0.01,0.5: {ok}")
