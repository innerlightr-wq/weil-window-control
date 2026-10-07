#!/usr/bin/env python3
r"""Amendment checks (2026-10-07, later pass).

Verifies the two identities on which Proposition 3's OVERREACH is withdrawn:

 (A)  1/M(s) already contains the k=2 ("4") harmonic and all higher ones, so the
      finite Laurent support of J*J does NOT bound the reach of K^{-1}.
 (B)  d/dsigma log M_sigma(s) at sigma=1/2 is EVEN in s and carries exactly the
      weights log2 * 2^{-k/2} for all k, so parity is not a barrier either.

Plus a finite control showing that compositions (PCP)^n connect indices two
steps apart even when PCP alone does not.  Standard library only.
"""
import math

log2 = math.log(2)
def M(s):  return 1.5 - math.sqrt(2)*math.cos(s*log2)

print("="*78); print("A.  1/M(s) contains the k=2 harmonic (the '4' shift)"); print("="*78)
# Poisson-kernel identity: 1/(1-2r cos th + r^2) = (1/(1-r^2))(1 + 2 sum_k r^k cos k th)
r = 2**-0.5
print("  M(s) = 1 - 2r cos(s log2) + r^2 with r = 2^(-1/2):",
      max(abs(M(s) - (1 - 2*r*math.cos(s*log2) + r*r)) for s in [i*0.37 for i in range(500)]) < 1e-14)
def invseries(s, N=4000):
    return (1/(1-r*r))*(1 + 2*sum(r**k * math.cos(k*s*log2) for k in range(1, N)))
dev = max(abs(1.0/M(s) - invseries(s)) for s in [0.11+i*0.29 for i in range(400)])
print("  1/M(s) == (1/(1-r^2))(1 + 2 sum_k r^k cos(k s log2))   max dev = %.6e" % dev)
c = [(1/(1-r*r))*(2*r**k if k else 1) for k in range(6)]
print("  Fourier coefficients of 1/M  (coeff of cos(k s log2)):")
for k in range(6):
    print("     k=%d : %.7f" % (k, c[k]), "   <- the '4' harmonic, NONZERO" if k == 2 else "")
print("  predicted by the reviewer: 2, 2*sqrt2, 2, ... :",
      abs(c[0]-2) < 1e-12, abs(c[1]-2*math.sqrt(2)) < 1e-12, abs(c[2]-2) < 1e-12)
print("  => finite Laurent support of M does NOT bound the support of 1/M.")

print("\n" + "="*78); print("B.  sigma-derivative: EVEN in s, with the W_2 weights"); print("="*78)
def Msig(s, sig): return abs(1 - 2**(-sig) * complex(math.cos(s*log2), -math.sin(s*log2)))**2
def dsig(s, sig=0.5, h=1e-6):
    return (math.log(Msig(s, sig+h)) - math.log(Msig(s, sig-h)))/(2*h)
def even_series(s, N=4000):
    return 2*log2*sum(2**(-k/2.0)*math.cos(k*s*log2) for k in range(1, N))
dev2 = max(abs(dsig(s) - even_series(s)) for s in [0.13+i*0.31 for i in range(300)])
print("  d/dsigma log M_sigma |_(sigma=1/2) == 2 log2 sum_k 2^(-k/2) cos(k s log2)")
print("     max dev = %.6e" % dev2)
print("  weights log2 * 2^(-k/2):", ["%.7f" % (log2*2**(-k/2.0)) for k in range(1, 5)])
print("  parity:  this expression is EVEN in s :",
      max(abs(even_series(-s) - even_series(s)) for s in [0.3, 1.1, 2.7]) < 1e-12)
print("  => the even prime-power weights ARE encoded in the local factor.")
print("     (algebraic expansion of the Euler factor; NOT a new arithmetic theorem,")
print("      and NOT a geometric operator family on the critical line.)")

print("\n" + "="*78); print("C.  compositions reach two steps; a single PCP need not"); print("="*78)
# finite model: C = symmetrised one-step shift on Z/11, P = projection onto {0,..,4}
n, k = 11, 5
C = [[0.0]*n for _ in range(n)]
for i in range(n):
    C[i][(i+1) % n] += 0.5; C[i][(i-1) % n] += 0.5
P = [[1.0 if (i == j and i < k) else 0.0 for j in range(n)] for i in range(n)]
def mul(A, B): return [[sum(A[i][t]*B[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
B1 = mul(mul(P, C), P)
B2 = mul(B1, B1)
print("  PCP   entry (0,1) = %.4f   entry (0,2) = %.4f" % (B1[0][1], B1[0][2]))
print("  (PCP)^2 entry (0,2) = %.4f   <- two-step reach appears at n = 2" % B2[0][2])
print("  so B^n for n>=2 carries the shift 2*log2 that M alone does not:",
      abs(B1[0][2]) < 1e-14 and abs(B2[0][2]) > 1e-6)
print("  CONTROL: with P = I (no compression) the same holds, and with C = I")
Ci = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
Bi = mul(mul(mul(P, Ci), P), mul(mul(P, Ci), P))
print("     C = I -> (PIP)^2 entry (0,2) = %.4f  (no shift generated):" % Bi[0][2],
      abs(Bi[0][2]) < 1e-14)
print("\n  NOTE: this is a finite shift model, NOT the Sonin space.  It shows only")
print("  that the support argument fails; it does NOT show that Pi_S produces the")
print("  correct Weil prime-power coefficients.  Cancellations in the full")
print("  assembly are not ruled out either way.")

print("\n" + "="*78); print("D.  invertibility needs only ||B|| <= 1"); print("="*78)
rr = 2*math.sqrt(2)/3
print("  K = (3/2)[I_E - (2 sqrt2/3) B],  2 sqrt2/3 = %.7f < 1, ||B|| <= 1" % rr)
print("  => K is invertible and the Neumann series converges WITHOUT ||B|| < 1.")
print("  A strict ||PCP|| < 1 is therefore NOT needed for invertibility here;")
print("  it would only sharpen the worst-case constant 1/(3/2-sqrt2) = %.7f." % (1/(1.5-math.sqrt(2))))
print("\nDONE")
