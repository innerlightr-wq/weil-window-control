#!/usr/bin/env python3
r"""Independent cross-check of Q_zeta(x) for the N=4 nodal witness.

Route G (used for the certificate): geometric side with the archimedean term as an
  m-sum of EXACT elementary integrals  -> certify_witness.weil_matrix_fast
Route S (independent):               geometric side with the archimedean term as the
  spectral integral (1/2pi) int |F_x(r)|^2 Re psi(1/4 + i r/2) dr, using an independently
  implemented digamma (asymptotic + recurrence), validated on psi(1), psi(1/4), psi(1/2).

The pole, log pi and prime pieces are shared; only the archimedean piece differs, and it is
the piece Limitation 4 of the manuscript says cannot be interval-evaluated with mpmath.
"""
import cmath, math, sys
sys.path.insert(0, '.')
import certify_witness as cw

EULER = cw.EULER


def digamma(z):
    s = 0j
    z = complex(z)
    while z.real < 15:
        s -= 1.0/z
        z += 1
    B = [1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6]
    r = cmath.log(z) - 1/(2*z)
    zz = z*z
    p = zz
    for n, b in enumerate(B, start=1):
        r -= b/(2*n*p)
        p *= zz
    return s + r


def simp(fn, lo, hi, n):
    h = (hi-lo)/n
    t = fn(lo) + fn(hi)
    for i in range(1, n):
        t += (4 if i % 2 else 2)*fn(lo+i*h)
    return t*h/3


L, g, N, M = 0.8, 14.0, 4, 2000000
av = [cw.Fk(k, L, g) for k in range(N)]
bv = [-cw.dFk(k, L, g) for k in range(N)]
na2 = sum(t*t for t in av)
ab = sum(av[i]*bv[i] for i in range(N))
vec = [bv[i] - (ab/na2)*av[i] for i in range(N)]
nrm = math.sqrt(sum(t*t for t in vec))
x = [t/nrm for t in vec]
print("N=4 nodal witness coefficients (orthonormal even_d basis), L=0.8, gamma=14 exactly:")
for k, t in enumerate(x):
    print(f"   x_{k} = {t:+.15f}")
print(f"   sum x_k^2 = {sum(t*t for t in x):.15f}")
print(f"   F(gamma) = sum x_k F_k(gamma) = {sum(x[k]*av[k] for k in range(N)):+.3e}"
      f"   (zero to rounding; charged, not assumed)")

Fx = lambda r: sum(x[k]*cw.Fk(k, L, r) for k in range(N))
print(f"\n(1/2pi) int |F_x|^2 dr over |r|<=6000 = "
      f"{simp(lambda r: Fx(r)**2, -6000, 6000, 3000000)/(2*math.pi):.10f}   (must be 1)")

# shared pieces
pole = 2*sum(x[k]*(-1)**k*2*cw.w_of(k, L)*math.cosh(L/2)
             / (cw.w_of(k, L)**2 + 0.25)/math.sqrt(L) for k in range(N))**2
prim = cw.von_mangoldt(2*L)
terms = {}
for j in range(N):
    for k in range(N):
        terms[(j, k)] = cw.Csym_terms(j, k, L)
def Cx(u):
    return sum(x[j]*x[k]*cw.eval_terms(terms[(j, k)], L, u)
               for j in range(N) for k in range(N))
prime = sum(2*lam/math.sqrt(n)*Cx(math.log(n)) for n, lam in prim)
print(f"\nshared pieces:  pole = {pole:.10f}   prime = {prime:.10f}   "
      f"logpi*C(0) = {cw.LOG_PI*Cx(0.0):.10f}   C(0) = {Cx(0.0):.12f}")

# route S
for R, n in ((2000, 2000000), (8000, 6000000)):
    archS = simp(lambda r: Fx(r)**2*digamma(0.25+0.5j*r).real, -R, R, n)/(2*math.pi)
    print(f"  route S arch (|r|<={R}) = {archS:+.10f}    "
          f"Q_zeta = {pole + archS - cw.LOG_PI*Cx(0.0) - prime:+.10f}")

# route G
A, E, _ = cw.weil_matrix_fast(L, N, M)
QzG = cw.quad_form(A, x)
QeG = cw.quad_form_err(E, x)
print(f"  route G Q_zeta (M={M}) = {QzG:+.10f} +- {QeG:.2e}")
