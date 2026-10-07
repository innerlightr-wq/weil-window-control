"""Validation for the Retraction-18 revision. Pure standard library.

The project's own regression scripts require mpmath (requirements.txt pins 1.3.0), which is
not installable in this environment (no pip, no uv, no system package). This script therefore
re-checks, independently and in double precision, exactly the three things the revision
touches -- the perturbation expansion, K(L,gamma), and the constants -- plus the logical
inference that was withdrawn.

Conventions follow paper/weil_window_control.tex 263-264 and src/zeta_window.py:
  f real, even, supp f in [-L,L], ||f||_2 = 1;  F(xi) = int f(u) e^{i xi u} du.
"""
import math, cmath

# ---------------------------------------------------------------- quadrature
def simpson(g, a, b, n=200001):
    if n % 2 == 0: n += 1
    h = (b - a) / (n - 1); s = 0.0
    for i in range(n):
        w = 1 if i in (0, n-1) else (4 if i % 2 else 2)
        s += w * g(a + i*h)
    return s * h / 3

def csimpson(g, a, b, n=200001):
    if n % 2 == 0: n += 1
    h = (b - a) / (n - 1); s = 0j
    for i in range(n):
        w = 1 if i in (0, n-1) else (4 if i % 2 else 2)
        s += w * g(a + i*h)
    return s * h / 3

# ---------------------------------------------------------------- admissible test functions
def normalise(f, L):
    nrm = math.sqrt(simpson(lambda u: f(u)**2, -L, L))
    return lambda u: f(u) / nrm

def F(f, L, xi):      return simpson(lambda u: f(u)*math.cos(xi*u), -L, L)
def Fp(f, L, xi):     return -simpson(lambda u: u*f(u)*math.sin(xi*u), -L, L)
def Fpp(f, L, xi):    return -simpson(lambda u: u*u*f(u)*math.cos(xi*u), -L, L)
def G(f, L, g, d):    return csimpson(lambda u: f(u)*cmath.exp(1j*(g+1j*d)*u), -L, L)

def c_of_L(L):   return 8*L**3*(1.0/3 + 1/math.sqrt(5))
def B(L, d):     return c_of_L(L)*d*d + (16.0/3)*L**5*d**4*math.exp(2*L*d)

FAIL = []
def check(name, cond, detail=""):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not cond: FAIL.append(name)

print("=" * 78)
print("1. Perturbation expansion (paper eq. pert) and the remainder bound (notes B.3)")
print("=" * 78)
print("   Q_delta - Q_0 = 4 Re G(delta)^2 - 4 F(gamma)^2")
print("                 = -4 delta^2 (F'^2 + F F'') + R_4,   |R_4| <= (16/3) L^5 d^4 e^{2Ld}")
for L, g, name in [(0.8, 14.1347, "u*sin(gamma u)"), (0.8, 14.1347, "cos(pi u/2L)"), (1.6, 14.1347, "u*sin(gamma u)")]:
    raw = (lambda u: u*math.sin(g*u)) if name.startswith("u*") else (lambda u: math.cos(math.pi*u/(2*L)))
    f = normalise(raw, L)
    f0, f1, f2 = F(f,L,g), Fp(f,L,g), Fpp(f,L,g)
    print(f"\n   L={L}, gamma={g}, f proportional to {name}")
    print(f"   {'delta':>8} {'actual drop':>16} {'-4d^2(F^2+FF\'\')':>18} {'|R_4|':>12} {'bound':>12} {'|drop|<=B_L':>12}")
    for d in (1e-1, 1e-2, 1e-3):
        actual = 4*(G(f,L,g,d)**2).real - 4*f0*f0
        quad   = -4*d*d*(f1*f1 + f0*f2)
        R4     = actual - quad
        bnd    = (16.0/3)*L**5*d**4*math.exp(2*L*d)
        ok     = abs(R4) <= bnd and abs(actual) <= B(L,d)
        print(f"   {d:>8.0e} {actual:>16.9e} {quad:>18.9e} {abs(R4):>12.3e} {bnd:>12.3e} {str(abs(actual)<=B(L,d)):>12}")
        check(f"remainder bound, L={L}, {name}, d={d:g}", abs(R4) <= bnd)
        check(f"Thm 7 per-f inequality Q_d >= Q_0 - B_L(d), L={L}, {name}, d={d:g}", abs(actual) <= B(L,d))

print()
print("=" * 78)
print("2. Cauchy-Schwarz constants (notes B.2) and c(L)")
print("=" * 78)
for L in (0.8, 1.6):
    f = normalise(lambda u: u*math.sin(14.1347*u), L)
    g = 14.1347
    check(f"|F|^2 <= 2L           (L={L})", F(f,L,g)**2   <= 2*L + 1e-9,      f"{F(f,L,g)**2:.6g} <= {2*L:.6g}")
    check(f"|F'|^2 <= 2L^3/3      (L={L})", Fp(f,L,g)**2  <= 2*L**3/3 + 1e-9, f"{Fp(f,L,g)**2:.6g} <= {2*L**3/3:.6g}")
    check(f"|F''|^2 <= 2L^5/5     (L={L})", Fpp(f,L,g)**2 <= 2*L**5/5 + 1e-9, f"{Fpp(f,L,g)**2:.6g} <= {2*L**5/5:.6g}")
    bound = 8*L**3*(1.0/3 + 1/math.sqrt(5))
    lhs = 4*abs(Fp(f,L,g)**2 + F(f,L,g)*Fpp(f,L,g))
    check(f"4|F'^2+FF''| <= c(L)  (L={L})", lhs <= bound + 1e-9, f"{lhs:.6g} <= {bound:.6g}")
check("c(0.8) ~ 6.244*0.8^3", abs(c_of_L(0.8) - 6.244*0.8**3) < 2e-3, f"{c_of_L(0.8):.6f}")

print()
print("=" * 78)
print("3. K(L,gamma): definition, closed form AS BRACKETED IN THE TEX, extremal function")
print("=" * 78)
def K_def(L, g):  return simpson(lambda u: u*u*math.sin(g*u)**2, -L, L)
def K_tex(L, g):
    return L**3/3 - ( L**2*math.sin(2*g*L)/(2*g) + L*math.cos(2*g*L)/(2*g**2)
                      - math.sin(2*g*L)/(4*g**3) )
g = 14.1347
for L in (0.8, 1.0, 1.3, 1.6, 2.0):
    check(f"K_tex == integral (L={L})", abs(K_tex(L,g)-K_def(L,g)) < 1e-9*max(1,K_def(L,g)),
          f"{K_tex(L,g):.12f} vs {K_def(L,g):.12f}")
for L, printed in [(0.8,0.7419),(1.0,1.3427),(1.3,3.1158),(1.6,5.1130),(2.0,10.652)]:
    check(f"4K matches printed table (L={L})", abs(4*K_tex(L,g)-printed) < 5e-4, f"{4*K_tex(L,g):.4f} vs {printed}")
for L in (0.8, 1.6):
    f = normalise(lambda u: u*math.sin(g*u), L)
    check(f"f prop. u sin(gamma u) attains |F'|^2 = K (L={L})",
          abs(Fp(f,L,g)**2 - K_def(L,g)) < 1e-8*max(1,K_def(L,g)),
          f"{Fp(f,L,g)**2:.10f} vs K={K_def(L,g):.10f}")
    check(f"K <= 2L^3/3 (Cor. 9, L={L})", K_def(L,g) <= 2*L**3/3 + 1e-12)

print()
print("=" * 78)
print("4. LOGICAL REGRESSION CHECK for the withdrawn inference (Retraction 18)")
print("=" * 78)
print("   'X >= b - D' permits replacing b by 0 only given a separate b >= 0.")
cases = [(-1.0, 0.0, -1.0, "negative scalar baseline b=-1, D=0, X=-1"),
         ( 8.9e-18, 1e-20, 8.8e-18, "certified-positive baseline (the L<=0.8 case)")]
for b, D, X, label in cases:
    premise   = X >= b - D - 1e-18
    concl     = X >= -D - 1e-18
    valid_sub = b >= 0
    print(f"\n   {label}")
    print(f"     premise  X >= b - D : {premise}")
    print(f"     b >= 0 available   : {valid_sub}")
    print(f"     substituted X >= -D: {concl}")
    check("substitution sound exactly when b >= 0", (concl if valid_sub else not concl) or (premise and valid_sub and concl))
print("\n   NOTE: the b=-1 row is a statement about SCALARS. It exhibits the failure of the")
print("   discarded inference; it is NOT an example of a negative actual Weil form.")

print()
print("=" * 78)
print(("ALL CHECKS PASSED" if not FAIL else f"{len(FAIL)} CHECK(S) FAILED: {FAIL}"))
print("=" * 78)
raise SystemExit(1 if FAIL else 0)
