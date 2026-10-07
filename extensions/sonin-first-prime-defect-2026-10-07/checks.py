#!/usr/bin/env python3
"""Minimal checks for the one-prime Sonin discrepancy test.

Standard library only (mpmath is not installable on this box).  Exact rational
arithmetic (fractions.Fraction) wherever a sign or an identity is at stake;
Decimal only for printing transcendental constants.  Every section carries at
least one deliberately failing control, so that a silent no-op is visible.
"""
from fractions import Fraction as Fr
from decimal import Decimal, getcontext
import math, random, cmath

getcontext().prec = 50
def e(x, n=6):            # never slice str(Decimal): it drops the exponent
    return ('%.' + str(n) + 'e') % float(x)

# ---------------------------------------------------------------- linear algebra over Q
def mm(A, B):
    n, k, m = len(A), len(B), len(B[0])
    return [[sum(A[i][t] * B[t][j] for t in range(k)) for j in range(m)] for i in range(n)]
def mt(A):       return [list(r) for r in zip(*A)]
def msub(A, B):  return [[A[i][j] - B[i][j] for j in range(len(A[0]))] for i in range(len(A))]
def eye(n):      return [[Fr(int(i == j)) for j in range(n)] for i in range(n)]
def iszero(A, tol=0):
    return all(abs(x) <= tol for r in A for x in r)
def inv(A):
    n = len(A); Ai = [list(A[i]) + eye(n)[i] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if Ai[r][c] != 0)
        Ai[c], Ai[p] = Ai[p], Ai[c]
        pv = Ai[c][c]; Ai[c] = [x / pv for x in Ai[c]]
        for r in range(n):
            if r != c and Ai[r][c] != 0:
                f = Ai[r][c]; Ai[r] = [a - f * b for a, b in zip(Ai[r], Ai[c])]
    return [row[n:] for row in Ai]
def trace(A):    return sum(A[i][i] for i in range(len(A)))

print("=" * 78)
print("1.  Pi_S = J P K^-1 P J*  IS the orthogonal projection onto J(M)  [exact over Q]")
print("=" * 78)
# M = span(e1,e2) inside Q^3 ; P = diag(1,1,0) ; J bounded invertible, NOT unitary,
# and deliberately NOT commuting with P.
P = [[Fr(1),Fr(0),Fr(0)],[Fr(0),Fr(1),Fr(0)],[Fr(0),Fr(0),Fr(0)]]
J = [[Fr(2),Fr(1),Fr(0)],[Fr(0),Fr(3),Fr(1)],[Fr(1),Fr(0),Fr(2)]]
Jt = mt(J)
print("  det J =", (lambda A: A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])
                           - A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
                           + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]))(J), " (invertible)")
print("  J unitary?            J^T J == I :", iszero(msub(mm(Jt,J), eye(3))))
print("  [J,P] == 0 ?                     :", iszero(msub(mm(J,P), mm(P,J))))
H  = mm(Jt, J)                                   # H = J* J
K2 = [[H[i][j] for j in (0,1)] for i in (0,1)]   # K = (P H P) restricted to M
K2i = inv(K2)
Kinv = [[Fr(0)]*3 for _ in range(3)]             # extended by zero on M-perp
for i in (0,1):
    for j in (0,1): Kinv[i][j] = K2i[i][j]
Pi = mm(mm(mm(J,P), mm(Kinv,P)), Jt)
print("  Pi self-adjoint   Pi == Pi^T     :", iszero(msub(Pi, mt(Pi))))
print("  Pi idempotent     Pi^2 == Pi     :", iszero(msub(mm(Pi,Pi), Pi)))
# range: Pi fixes J(M) pointwise and kills (J(M))^perp
JM = [[J[i][c] for c in (0,1)] for i in range(3)]
print("  Pi J|_M == J|_M  (fixes J(M))    :", iszero(msub(mm(Pi,JM), JM)))
n_ = [Fr(0)]*3                                   # normal of J(M) = span(Je1,Je2)
a, b = [J[i][0] for i in range(3)], [J[i][1] for i in range(3)]
n_ = [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
print("  Pi n == 0  (kills J(M)^perp)     :",
      all(sum(Pi[i][j]*n_[j] for j in range(3)) == 0 for i in range(3)))
print("  rank check  tr Pi == dim M == 2  :", trace(Pi) == 2, " tr =", trace(Pi))

print("\n  CONTROL 2 of the brief: the two objects the formula is NOT")
A1 = mm(mm(J,P), inv(J))      # J P J^-1
A2 = mm(mm(J,P), Jt)          # J P J*
print("   J P J^-1 : idempotent =", iszero(msub(mm(A1,A1),A1)),
      " self-adjoint =", iszero(msub(A1, mt(A1))))
print("   J P J*   : idempotent =", iszero(msub(mm(A2,A2),A2)),
      " self-adjoint =", iszero(msub(A2, mt(A2))))
print("   J P J^-1 == Pi ?", iszero(msub(A1,Pi)), "   J P J* == Pi ?", iszero(msub(A2,Pi)))
print("   tr(J P J^-1) == tr P  (similarity invariance) :", trace(A1) == trace(P),
      " ->", trace(A1), "==", trace(P))
print("   J P J* >= 0 but NOT a projection: tr =", trace(A2))
print("  CONTROL 5 (outside the hypotheses): J singular -> K not invertible")
Jb = [[Fr(1),Fr(1),Fr(0)],[Fr(1),Fr(1),Fr(0)],[Fr(0),Fr(0),Fr(1)]]   # rank 2, kills e1-e2
Hb = mm(mt(Jb), Jb); Kb = [[Hb[i][j] for j in (0,1)] for i in (0,1)]
detKb = Kb[0][0]*Kb[1][1] - Kb[0][1]*Kb[1][0]
print("   det K =", detKb, "-> K^-1 does not exist, formula inapplicable:", detKb == 0)

print("\n" + "=" * 78)
print("2.  Two-projection structure of  G = Phat P - Q   (Route A defect)")
print("=" * 78)
def rand_proj(n, k, rng):
    V = [[rng.gauss(0,1) for _ in range(k)] for _ in range(n)]    # Gram-Schmidt
    cols = []
    for j in range(k):
        v = [V[i][j] for i in range(n)]
        for c in cols:
            d = sum(v[i]*c[i] for i in range(n)); v = [v[i]-d*c[i] for i in range(n)]
        nr = math.sqrt(sum(t*t for t in v)); cols.append([t/nr for t in v])
    return [[sum(c[i]*c[j] for c in cols) for j in range(n)] for i in range(n)]
def fmul(A,B):
    n,k,m=len(A),len(B),len(B[0])
    return [[sum(A[i][t]*B[t][j] for t in range(k)) for j in range(m)] for i in range(n)]
def ftr(A): return sum(A[i][i] for i in range(len(A)))
def sym_eig(A, iters=400):           # Jacobi, updating from SAVED originals
    n=len(A); B=[list(r) for r in A]
    for _ in range(iters):
        off=0.0; p=q=0
        for i in range(n):
            for j in range(i+1,n):
                if abs(B[i][j])>off: off=abs(B[i][j]); p,q=i,j
        if off<1e-14: break
        app,aqq,apq=B[p][p],B[q][q],B[p][q]
        th=0.5*math.atan2(2*apq, aqq-app); c,s=math.cos(th),math.sin(th)
        rowp=list(B[p]); rowq=list(B[q])
        for k in range(n):
            B[p][k]=c*rowp[k]-s*rowq[k]; B[q][k]=s*rowp[k]+c*rowq[k]
        colp=[B[i][p] for i in range(n)]; colq=[B[i][q] for i in range(n)]
        for k in range(n):
            B[k][p]=c*colp[k]-s*colq[k]; B[k][q]=s*colp[k]+c*colq[k]
    return sorted(B[i][i] for i in range(n))
# PSD control for the eigensolver (this is the bug that once returned -8.16)
g=[1.0,2.0,-1.0,0.5]; G1=[[g[i]*g[j] for j in range(4)] for i in range(4)]
ev=sym_eig(G1); print("  eigensolver PSD control (rank-1 g g^T): min ev =", e(min(ev)),
      " trace preserved =", abs(sum(ev)-ftr(G1))<1e-10)
rng=random.Random(20261007)
print("  n  k1 k2 | dim(ran P n ran Ph) | tr(PhP)-trQ = sum c^2 | ||PhP-Q||_1 = sum c")
for (n,k1,k2) in [(6,3,3),(6,4,3),(8,5,4),(8,4,4)]:
    Pp=rand_proj(n,k1,rng); Ph=rand_proj(n,k2,rng)
    R=fmul(Ph,Pp)
    # B = ran P cap ran Ph  <=>  eigenvalue 1 of P Ph P ; generic part: eigenvalues in (0,1)
    PPhP=fmul(fmul(Pp,Ph),Pp); lam=sym_eig(PPhP)
    ones=[t for t in lam if t>1-1e-9]; gen=[t for t in lam if 1e-9<t<1-1e-9]
    sum_c2=sum(gen); sum_c=sum(math.sqrt(t) for t in gen)
    # Q = spectral projection of P Ph P at eigenvalue 1 ; here generically dim 0
    trQ=len(ones)
    print("  %d  %d  %d  |        %d            |  %s   (trR-trQ=%s) |  %s"
          %(n,k1,k2,len(ones),e(sum_c2),e(ftr(R)-trQ),e(sum_c)))
print("  COMMUTING CONTROL (must give G == 0):")
Pp=[[1.0 if i==j and i<3 else 0.0 for j in range(6)] for i in range(6)]
Ph=[[1.0 if i==j and i in (1,2,4) else 0.0 for j in range(6)] for i in range(6)]
R=fmul(Ph,Pp); lam=sym_eig(fmul(fmul(Pp,Ph),Pp))
gen=[t for t in lam if 1e-9<t<1-1e-9]
print("   [P,Ph]==0 :", max(abs(fmul(Pp,Ph)[i][j]-fmul(Ph,Pp)[i][j])
                            for i in range(6) for j in range(6))<1e-12,
      "  generic part empty :", len(gen)==0, "  tr(PhP)=", e(ftr(R)),
      " = dim(ran cap ran) =", sum(1 for t in lam if t>1-1e-9))

print("\n" + "=" * 78)
print("3.  The theta_S multiplier at p = 2 and its exact constants")
print("=" * 78)
l2 = Decimal(2).ln()
def Mfun(s):   # |m_2(s)|^2 with m_2(s) = 1 - 2^(-1/2 - i s)
    m = 1 - 2**(-0.5) * cmath.exp(-1j*s*math.log(2))
    return abs(m)**2
def Mclosed(s): return 1.5 - math.sqrt(2)*math.cos(s*math.log(2))
worst = max(abs(Mfun(s)-Mclosed(s)) for s in [i*0.37 for i in range(400)])
print("  |m_2(s)|^2 == 3/2 - sqrt(2) cos(s log 2)   max abs dev over 400 pts:", e(worst))
mmin = 1.5-math.sqrt(2); mmax = 1.5+math.sqrt(2)
print("  min M = 3/2-sqrt2 =", e(mmin), "   max M = 3/2+sqrt2 =", e(mmax))
print("  mean of M over a period (cos averages to 0)  = 3/2 =",
      e(sum(Mclosed(i*2*math.pi/math.log(2)/2000) for i in range(2000))/2000))
print("  ambient (worst-case) two-sided factor  max/min =", e(mmax/mmin))
print("  amplitude condition number sqrt(max/min)        =", e(math.sqrt(mmax/mmin)))
print("  worst-case  ||K^-1|| <= 1/(3/2-sqrt2) =", e(1/mmin),
      " ; mean-field 1/(3/2) =", e(2/3.0), " ; GAP FACTOR =", e((1/mmin)/(2/3.0)))
print("  zeros of M on R? m_2(s)=0 needs 2^(-1/2)=1 :", False,
      " -> M > 0 everywhere; min attained only on s in (2pi/log2)Z, i.e. spacing",
      e(2*math.pi/math.log(2)))
print("\n  SIGN / OMISSION CONTROL (must change the numbers):")
print("   wrong sign  3/2 + sqrt2 cos : min =", e(min(1.5+math.sqrt(2)*math.cos(i*0.01) for i in range(1000))),
      " (same min by symmetry, but the operator differs: see section 4)")
print("   omitting the oscillation (M := 3/2) : max/min = 1.0 != ", e(mmax/mmin))
print("   using 2^-1 instead of 2^-1/2 : max/min =", e((1+0.5)**2/(1-0.5)**2),
      " (wrong: that is |1-2^-1-is|^2 with the WRONG normalisation)")

print("\n" + "=" * 78)
print("4.  M as a three-term difference operator;  the L = (log2)/2 threshold")
print("=" * 78)
print("  M(s) = 3/2 - 2^(-1/2) (2^(is) + 2^(-is))      [exact trig identity]")
chk = max(abs(Mclosed(s) - (1.5 - 2**-0.5*(cmath.exp(1j*s*math.log(2))
                                          + cmath.exp(-1j*s*math.log(2))).real))
          for s in [i*0.41 for i in range(400)])
print("  max abs dev:", e(chk))
print("  mult by 2^(-is)  <->  scaling lambda -> lambda/2   (CCM2024 proof of Prop 4.6(ii))")
print("  so M F  has the three shifted pieces  f, f_2, f_(1/2)  at  u = 0, +log2, -log2")
log2 = math.log(2)
print("   L          supp f        supp f_2              disjoint from supp f ?")
for L in [0.20, 0.30, 0.3465735, 0.35, 0.40, 0.5493061]:
    dis = (L + L) < log2 - 1e-12
    print("   %-10.7f [%7.4f,%7.4f] [%7.4f,%7.4f]   %s"
          %(L, -L, L, log2-L, log2+L, "YES (K is scalar)" if dis else "no (K genuinely non-scalar)"))
print("  threshold  2L = log 2  <=>  L = (log2)/2 =", e(log2/2, 7))

print("\n" + "=" * 78)
print("5.  Prime activation on the calibration window  (log2)/2 < L < (log3)/2")
print("=" * 78)
def vonmangoldt(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return (p, math.log(p)) if m == 1 else None
    return None
print("   L        e^(2L)   active n (log n < 2L, Lambda(n)!=0)")
for L in [0.3465735, 0.34657359, 0.36, 0.45, 0.5493061, 0.55, 0.8, 1.0]:
    act = [n for n in range(2, 40) if vonmangoldt(n) and math.log(n) < 2*L]
    print("   %-9.8f %-8.5f %s" % (L, math.exp(2*L), act))
print("  exact boundaries: (log2)/2 =", e(log2/2,8), "  (log3)/2 =", e(math.log(3)/2,8))
print("  BOUNDARY CONTROL: at L = (log2)/2 exactly, log 2 < 2L is FALSE ->",
      math.log(2) < 2*(log2/2), " : prime block EMPTY, test vacuous (window must be open)")
print("  Zhu's certified window L <= 0.8 contains n =", [n for n in range(2,40)
      if vonmangoldt(n) and math.log(n) < 1.6])

print("\n" + "=" * 78)
print("6.  SEPARATION: |L_2(1/2+is)|^-2  versus  -(L_2'/L_2)(1/2+is)")
print("=" * 78)
print("  |m_2(s)|^2 = |L_2(1/2+is)|^-2 has EXACTLY the shifts {0, +-log2}:")
print("    coefficient at k=0 : 3/2 ;  at k=+-1 : -2^(-1/2) =", e(-2**-0.5), " ; at |k|>=2 : 0")
print("  -(L_2'/L_2)(1/2+is) = -log2 * sum_{k>=1} 2^(-k/2) 2^(-i k s)  -> ALL k>=1:")
for k in range(1, 7):
    print("    k=%d  weight log2 * 2^(-k/2) = %s" % (k, e(math.log(2)*2**(-k/2.0))))
def dlogM(s, hh=1e-6):      # d/ds log M(s)
    return (math.log(Mclosed(s+hh)) - math.log(Mclosed(s-hh)))/(2*hh)
def series(s, N=4000):      # 2 sum_k log2 2^(-k/2) sin(k s log2)
    return 2*sum(math.log(2)*2**(-k/2.0)*math.sin(k*s*math.log(2)) for k in range(1, N))
dev = max(abs(dlogM(s)-series(s)) for s in [0.13+i*0.29 for i in range(300)])
print("  d/ds log M(s) == 2 sum_{k>=1} log2 * 2^(-k/2) sin(k s log2) : max dev =", e(dev))
print("  -> the full W_2 coefficient sequence log2*2^(-k/2) lives in M only through")
print("     a LOGARITHMIC DERIVATIVE, not linearly.  M itself carries 3 modes; W_2 needs all.")
print("\n  CONTROL: a 3-mode trig polynomial cannot equal the geometric series.")
best = None
for a0 in [i/20.0 for i in range(-60, 61)]:
    for a1 in [i/20.0 for i in range(-60, 61)]:
        err = max(abs((a0 + 2*a1*math.cos(s*math.log(2))) - series(s)) for s in
                  [0.0, 0.7, 1.4, 2.1, 2.8, 3.5, 4.2])
        if best is None or err < best[0]: best = (err, a0, a1)
print("   best L^inf fit of a0 + 2 a1 cos(s log2) to the W_2 series on 7 nodes:")
print("   err =", e(best[0]), " at a0 =", best[1], " a1 =", best[2],
      "  (series is ODD, the polynomial is EVEN -> no fit)")
print("   oddness check: series(-s) == -series(s) :",
      max(abs(series(-s)+series(s)) for s in [0.3,1.1,2.7])<1e-12,
      " ; M(-s) == M(s) :", max(abs(Mclosed(-s)-Mclosed(s)) for s in [0.3,1.1,2.7])<1e-15)

print("\n" + "=" * 78)
print("7.  Targeted versus ambient bound on the SAME objects")
print("=" * 78)
print("  ambient   : (m/Mx) P_inf <= P_S <= (Mx/m) P_inf  with Mx/m =", e(mmax/mmin))
print("  targeted  : K = (3/2)P - sqrt2 * (P C P),  C = (vartheta(2)+vartheta(2)^-1)/2, ||C|| = 1")
print("              so  K >= (3/2 - sqrt2 ||PCP||) P  and the gain is governed by ||PCP||")
for t in [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]:
    lo = 1.5 - math.sqrt(2)*t
    print("     ||PCP|| = %.2f  ->  K >= %s P  ->  ||K^-1|| <= %s   (ambient %s)"
          % (t, e(lo), e(1/lo) if lo > 0 else "inf", e(1/mmin)))
print("  NOTE: ||PCP|| <= 1 always; ||PCP|| < 1 is NOT established by any checked source.")
print("  CONTROL: at ||PCP|| = 1 the targeted bound DEGENERATES to the ambient one:",
      abs(1/(1.5-math.sqrt(2)*1.0) - 1/mmin) < 1e-12)
print("\nDONE")
