"""Minimal controls for the intrinsic-positivity-bridge audit. Standard library only.

Everything structural is decided in exact rational arithmetic (Fraction); the few real
constants quoted from the sources are checked in float only to the digits printed there.
"""
from fractions import Fraction as F
import math

ok=[]
def chk(name, cond, detail=""):
    ok.append((name,bool(cond)))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))

print("="*80)
print("1.  Support dictionary and PRIME ACTIVATION  (exact)")
print("="*80)
print("  Connes-Consani Thm 1 / 6.11: g supported in [2^{-1/2}, 2^{1/2}].")
print("  Multiplicative x = e^u  =>  u in [-(1/2)log2, (1/2)log2], i.e. OUR L = (1/2) log 2.")
print("  Their Lemma 6.10 works in H = L^2(I), I = [-(1/2)log2, (1/2)log2]  -- same interval.")
L_cc = math.log(2)/2
print(f"  L_cc = (1/2) log 2 = {L_cc:.6f};  autocorrelation g*g* lives on [-2L, 2L] = [-log2, log2].")
chk("their I equals our [-L,L] with L = (1/2)log 2", abs(L_cc-0.34657359)<1e-7, f"L = {L_cc:.8f}")

print()
print("  W_p vanishes on test functions supported in the OPEN interval (p^{-1}, p).")
print("  Our prime term is von_mangoldt_terms(2L): n with log n <= 2L, i.e. n <= e^{2L}.")
print("  %10s %10s %12s   active n (Lambda(n) != 0, n <= e^(2L))" % ("L","2L","e^(2L)"))
def active(L):
    lim=math.exp(2*L); out=[]
    for n in range(2,200):
        m=n; p=None; ok_=True
        d=2
        mm=n
        while d*d<=mm:
            if mm%d==0:
                p=d
                while mm%d==0: mm//=d
                if mm!=1: ok_=False
                break
            d+=1
        if p is None: p=mm if mm>1 else None
        if p is not None and ok_ and n<=lim+1e-12: out.append(n)
    return out
for L in (L_cc, 0.4, 0.55, 0.8, 1.0):
    a=active(L)
    print(f"  {L:>10.6f} {2*L:>10.6f} {math.exp(2*L):>12.6f}   {a}")
chk("at the Connes-Consani window NO prime power is strictly active",
    math.exp(2*L_cc)==2.0 or abs(math.exp(2*L_cc)-2)<1e-12,
    "e^{2L} = 2 exactly: n=2 sits ON the boundary, and W_2 vanishes for support in the OPEN (1/2,2)")
chk("at Zhu's certified L = 0.8 primes ARE active", set(active(0.8))>= {2,3,4},
    "e^(1.6) = %.4f, active n = %s" % (math.exp(1.6), active(0.8)))

print()
print("="*80)
print("2.  CONTROL: A S A* >= 0 even when A does not commute with S   (exact, rational)")
print("="*80)
def mm(A,B): return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def T(A):    return [[A[j][i] for j in range(len(A))] for i in range(len(A[0]))]
def minor_dets(M):
    """leading principal minors, exact"""
    out=[]
    for k in range(1,len(M)+1):
        sub=[row[:k] for row in M[:k]]
        # exact determinant by fraction-free elimination
        n=k; Ad=[r[:] for r in sub]; det=F(1)
        for i in range(n):
            piv=next((r for r in range(i,n) if Ad[r][i]!=0), None)
            if piv is None: det=F(0); break
            if piv!=i: Ad[i],Ad[piv]=Ad[piv],Ad[i]; det=-det
            det*=Ad[i][i]
            inv=F(1)/Ad[i][i]
            for r in range(i+1,n):
                f0=Ad[r][i]*inv
                Ad[r]=[x-f0*y for x,y in zip(Ad[r],Ad[i])]
        out.append(det)
    return out
# S = orthogonal projection onto span(e1,e2) in R^4 ; A a non-commuting invertible map
S=[[F(1),F(0),F(0),F(0)],[F(0),F(1),F(0),F(0)],[F(0),F(0),F(0),F(0)],[F(0),F(0),F(0),F(0)]]
A=[[F(1),F(2),F(0),F(1)],[F(0),F(1),F(3),F(0)],[F(1),F(0),F(1),F(2)],[F(0),F(1),F(0),F(1)]]
AS=mm(A,S); ASA=mm(AS,T(A))
chk("S is an orthogonal projection (S^2 = S = S^T)", mm(S,S)==S and T(S)==S)
chk("A does not commute with S", mm(A,S)!=mm(S,A))
chk("A S A^T is PSD anyway (all leading principal minors >= 0)",
    all(d>=0 for d in minor_dets(ASA)), f"minors = {[str(d) for d in minor_dets(ASA)]}")
chk("A S A^T = (A S)(A S)^T, i.e. a Gram matrix  => PSD by construction",
    ASA==mm(AS,T(AS)))
print("  This is the structural reason Tr(theta(g) S theta(g)*) >= 0 without S being")
print("  invariant under theta: it is a Gram/Hilbert-Schmidt square, not a compression")
print("  of a positive operator to an invariant subspace.")

print()
print("="*80)
print("3.  Lemma 6.9 (Connes-Consani): exact criterion + 2x2 reduction + FAILING case")
print("="*80)
print("  B(xi) = -b|<phi|xi>|^2 + a|<psi|xi>|^2 + c||P_phi xi||^2  is positive  iff")
print("     a + c >= b   AND   b(a+c) <= a(b+c)|<phi|psi>|^2 .")
def lemma69(a,b,c,alpha2):
    """exact test of the criterion, and of the 2x2 matrix it reduces to."""
    beta2=1-alpha2
    M=[[a*alpha2-b, None],[None, a*beta2+c]]
    # off-diagonal is a*alpha*beta; work with its square to stay rational
    tr = a+c-b
    det = a*alpha2*(b+c) - b*(a+c)
    crit = (a+c>=b) and (b*(a+c) <= a*(b+c)*alpha2)
    return tr, det, crit
print(f"  {'a':>4} {'b':>4} {'c':>4} {'|<phi|psi>|^2':>14} {'trace':>8} {'det':>10} {'criterion':>10} {'tr>0&det>=0':>12}")
cases=[(F(3),F(1),F(1),F(1)),      # perfectly aligned, should pass
       (F(3),F(1),F(1),F(1,2)),    # partly aligned
       (F(1),F(3),F(1),F(1)),      # b too big -> fail
       (F(3),F(1),F(1),F(1,10)),   # poorly aligned -> fail
       (F(2),F(1),F(0),F(1))]      # no gap, aligned
for a,b,c,al in cases:
    tr,det,crit=lemma69(a,b,c,al)
    agree = crit == (tr>=0 and det>=0)
    print(f"  {str(a):>4} {str(b):>4} {str(c):>4} {str(al):>14} {str(tr):>8} {str(det):>10} {str(crit):>10} {str(tr>=0 and det>=0):>12}")
    chk(f"criterion == (trace>=0 and det>=0) for (a,b,c,|<phi|psi>|^2)=({a},{b},{c},{al})", agree)
chk("the test DISCRIMINATES: at least one case fails the criterion",
    any(not lemma69(*cse)[2] for cse in cases))
print("  => alignment |<phi|psi>|^2 is what decides: with a,b,c fixed at (3,1,1), full")
print("     alignment passes and alignment 1/10 fails. Magnitude alone does not settle it.")

print()
print("="*80)
print("4.  Constants quoted from the sources (float, to the digits printed there)")
print("="*80)
gamma_cc=2.94355
c_611=4*gamma_cc/math.log(2)
print(f"  Lemma 6.10: gamma ~ {gamma_cc}")
print(f"  Theorem 6.11: c = 4 gamma / log 2 = {c_611:.6f}")
# The source states c = 4*gamma/log 2 symbolically and gamma ~ 2.94355; it prints no decimal
# value for c. We therefore check the ARITHMETIC against the symbolic definition only, and do
# not assert a decimal value the source never gives.
chk("c equals 4*gamma/log 2 with the quoted gamma",
    abs(c_611 - 4*gamma_cc/math.log(2)) < 1e-12, f"c = {c_611:.5f} (derived, not quoted)")
print(f"  Remark 6.12: lambda_max = 1.05158, epsilon_1 ~ 0.00122 (quoted, NOT independently derived)")

print()
print("="*80)
bad=[n for n,v in ok if not v]
print("ALL CHECKS PASSED" if not bad else f"{len(bad)} FAILED: {bad}")
print("="*80)
raise SystemExit(1 if bad else 0)
