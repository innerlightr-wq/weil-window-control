"""Exact prolate rungs delta_n(lambda) = 1 - |chi_n(lambda)| by Bouwkamp's method,
   following Appendix A of arXiv:2611.13129v1 (eqs 74-78).

   Validation ladder (paper's App. A): (i) stability under K -> 2K; (ii) gamma -> 0
   reproduces Legendre; (iv) eq (78) against the paper's printed chi_0, chi_4 at
   lambda^2 = 5.
"""
import sys
from mpmath import mp, mpf, matrix, eigsy, sqrt, log, exp, pi, fabs, mpmathify

def pbar0(k):
    """normalized Legendre P-bar_k(0) for even k, eq (78)."""
    num = mpf(1); den = mpf(1)
    for j in range(1, k, 2): num *= j          # (k-1)!!
    for j in range(2, k+1, 2): den *= j        # k!!
    return sqrt(mpf(2*k+1)/2) * (-1)**(k//2) * num/den

def prolate_even(gamma, K):
    """Bouwkamp tridiagonal (75) on even Legendre indices k=0,2,...,K.
       Returns (eigenvalues mu, eigenvector matrix columns d)."""
    ks = list(range(0, K+1, 2))
    m = len(ks)
    A = matrix(m, m)
    g2 = gamma**2
    for i, k in enumerate(ks):
        k = mpf(k)
        A[i, i] = k*(k+1) + g2*(2*k*k + 2*k - 1)/((2*k-1)*(2*k+3))
        if i+1 < m:
            off = g2*(k+1)*(k+2)/((2*k+3)*sqrt((2*k+1)*(2*k+5)))
            A[i, i+1] = off; A[i+1, i] = off
    return ks, eigsy(A)

def chis(lam, K, nmax=10):
    """chi_n for even n = 0,2,4,..., via eq (78)."""
    gamma = 2*pi*lam**2
    ks, (mu, Dm) = prolate_even(gamma, K)
    out = {}
    pb = [pbar0(k) for k in ks]
    for idx in range(0, nmax//2 + 1):
        d = [Dm[i, idx] for i in range(len(ks))]
        phi0 = sum(dk*p for dk, p in zip(d, pb))
        out[2*idx] = lam*sqrt(mpf(2))*d[0]/phi0
    return out

def report(lam2, K, dps):
    mp.dps = dps
    lam = sqrt(mpf(lam2)) if not isinstance(lam2, str) else mpmathify(lam2)
    c = chis(lam, K)
    print(f"  lambda^2 = {lam**2}   gamma = 2 pi lambda^2 = {2*pi*lam**2}   K={K} dps={dps}")
    for n in (0, 2, 4, 6, 8):
        d = 1 - fabs(c[n])
        print(f"    chi_{n} = {mp.nstr(c[n], 32)}   delta_{n} = {mp.nstr(d, 12)}"
              f"   -log delta_{n} = {mp.nstr(-log(d), 10)}")
    return c

print("=== (ii) gamma -> 0 reproduces Legendre: chi_n(0) should be ... (lam->0) ===")
mp.dps = 30
ks, (mu, Dm) = prolate_even(mpf('1e-20'), 40)
print("    mu_n at gamma~0:", [mp.nstr(mu[i], 8) for i in range(4)], " expected n(n+1):", [0, 6, 20, 42])

print("\n=== (iv) VALIDATION at lambda^2 = 5 against paper App. A printed values ===")
PAPER_CHI0 = "0.99999999999999999999999998989"
PAPER_CHI4 = "0.99999999999999999870161199128"
for K in (140, 280):
    c = report(5, K, 80)
    mp.dps = 80
    print(f"    vs paper chi_0: diff = {mp.nstr(c[0]-mpmathify(PAPER_CHI0), 6)}")
    print(f"    vs paper chi_4: diff = {mp.nstr(c[4]-mpmathify(PAPER_CHI4), 6)}")
    print(f"    delta_4(lambda^2=5) = {mp.nstr(1-c[4], 20)}")

print("\n=== TARGET: lambda = e^{4/5}  (our certified window L = 4/5) ===")
for K, dps in ((140, 80), (280, 80), (280, 110)):
    mp.dps = dps
    lam = exp(mpf(4)/5)
    c = chis(lam, K)
    print(f"  K={K} dps={dps}  lambda^2={mp.nstr(lam**2,20)}")
    for n in (0, 2, 4, 6):
        d = 1-fabs(c[n])
        print(f"    delta_{n} = {mp.nstr(d, 20)}    -log delta_{n} = {mp.nstr(-log(d), 12)}")
