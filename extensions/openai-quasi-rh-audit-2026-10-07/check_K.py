"""Verification of K(L,gamma) and eq. (5) of paper/weil_window_control.tex.

The manuscript (tex 308-316) defines

    K(L,gamma) = sup{ |F'(gamma)|^2 : f real even, supp f in [-L,L], ||f||_2 = 1 }
               = int_{-L}^{L} u^2 sin^2(gamma u) du      [attained at f prop. u sin(gamma u)]

with the closed form, brackets AS IN THE TEX:

    K = L^3/3 - [ L^2 sin(2gL)/(2g) + L cos(2gL)/(2g^2) - sin(2gL)/(4g^3) ]

VERDICT: CORRECT. Four independent checks below.

NOTE (recorded deliberately). An earlier pass of this audit reported two sign errors here.
That was wrong, and the cause is worth keeping: `pdftotext` renders the tall \left[ \right]
delimiters as stray '"' and '#' characters on their own line, so the bracket disappears from
the extracted text and the last two signs appear flipped. The variant K_nobracket below is
that misreading; it is kept only as the negative control that the ODE test discriminates.
Formulas in this audit were re-read from the TeX source after this.
"""
import math

def K_quad(L, g, n=400001):
    """Simpson on [0,L], doubled (integrand even)."""
    h = L / (n - 1); s = 0.0
    for i in range(n):
        u = i * h
        w = 1 if i in (0, n-1) else (4 if i % 2 else 2)
        s += w * u*u * math.sin(g*u)**2
    return 2 * s * h / 3

def K_tex(L, g):
    """Closed form exactly as the TeX reads (bracketed)."""
    return L**3/3 - ( L**2*math.sin(2*g*L)/(2*g)
                      + L*math.cos(2*g*L)/(2*g**2)
                      - math.sin(2*g*L)/(4*g**3) )

def K_nobracket(L, g):
    """NEGATIVE CONTROL: the bracket-dropped misreading. Should fail every test."""
    return L**3/3 - L**2*math.sin(2*g*L)/(2*g) \
                  + L*math.cos(2*g*L)/(2*g**2) \
                  - math.sin(2*g*L)/(4*g**3)

def d_dL(f, L, g, h=1e-6):
    return (f(L+h, g) - f(L-h, g)) / (2*h)

g = 14.1347                      # planted ordinate used in the manuscript's tables
LS = (0.6, 0.8, 1.0, 1.3, 1.6, 2.0)

print("CHECK 1  dK/dL must equal 2 L^2 sin^2(gL)  [decisive]")
print(f"{'L':>5} {'target':>16} {'K_tex':>16} {'K_nobracket':>16}")
for L in LS:
    print(f"{L:>5} {2*L*L*math.sin(g*L)**2:>16.10f} {d_dL(K_tex,L,g):>16.10f} {d_dL(K_nobracket,L,g):>16.10f}")

print("\nCHECK 2  value against quadrature")
print(f"{'L':>5} {'quadrature':>16} {'|K_tex-quad|':>14} {'|K_nobr-quad|':>15}")
for L in (0.8, 1.0, 1.6):
    q = K_quad(L, g)
    print(f"{L:>5} {q:>16.10f} {abs(K_tex(L,g)-q):>14.2e} {abs(K_nobracket(L,g)-q):>15.2e}")

print("\nCHECK 3  the manuscript's printed table of 4K(L,gamma_1)  (tex 350)")
printed = {0.8:0.7419, 1.0:1.3427, 1.3:3.1158, 1.6:5.1130, 2.0:10.652}
for L, v in printed.items():
    print(f"  L={L:>4}  printed {v:>8.4f}   4*K_tex {4*K_tex(L,g):>8.4f}   4*K_nobracket {4*K_nobracket(L,g):>8.4f}")

print("\nCHECK 4  the manuscript's deviation table |K/(L^3/3)-1|  (tex 361)")
print("  printed: 17.3%, 8.68%, 0.699%, 6.37%, 6.38%, 0.136%")
print("  K_tex:  " + ", ".join(f"{abs(K_tex(L,g)/(L**3/3)-1)*100:.3f}%" for L in LS))
print("  K_nobr: " + ", ".join(f"{abs(K_nobracket(L,g)/(L**3/3)-1)*100:.3f}%" for L in LS))

print("\nCHECK 5  Corollary 9: K <= 2L^3/3 for every gamma, and K -> L^3/3 as gamma -> inf")
for L in (0.8, 2.0):
    print(f"  L={L}:  K={K_tex(L,g):.6f}  2L^3/3={2*L**3/3:.6f}  L^3/3={L**3/3:.6f}"
          f"   K at gamma=1e6: {K_tex(L,1e6):.6f}")
