"""Cost ledger for the only quantity the 7/8 theorem supplies to the window form.

OpenAI Thm 1.1 (pinned adc7f124) => every nontrivial zero has Re in [1/8, 7/8],
i.e. |delta_rho| <= 3/8 where delta = Re rho - 1/2.  In the window Weil form the
depth delta enters ONLY through the quartet contribution

    Q_delta(f) - Q_0(f) = 4 Re G(delta)^2 - 4 F(gamma)^2,   G(delta) = F(gamma + i delta),

with, by Cauchy-Schwarz on [-L,L] for ||f||_2 = 1,
    |G(delta)|^2 <= sinh(2 L delta)/delta ,      F(gamma)^2 <= 2L .

Two valid lower bounds for the drop:
  (P)  paper Thm 7 :  -( c(L) delta^2 + (16/3) L^5 delta^4 e^{2 L delta} ),  c(L)=8L^3(1/3+1/sqrt5)
  (D)  direct      :  -( 4 sinh(2 L delta)/delta + 8L )
(P) is a rigorous series-tail bound, valid for all delta -- not an asymptotic
expansion -- but it degrades as L delta grows. (D) does not vanish as delta->0
but is flat in delta. The useful bound is the max (least negative).
"""
import math

def cL(L):  return 8*L**3*(1.0/3 + 1/math.sqrt(5))
def P(L,d): return cL(L)*d*d + (16.0/3)*L**5*d**4*math.exp(2*L*d)
def D(L,d): return 4*math.sinh(2*L*d)/d + 8*L

print("=== 1. Is the paper's Thm-7 bound still informative at delta = 3/8 (the 7/8 cap)? ===")
print(f"{'L':>5} {'c(L)d^2':>12} {'(16/3)L^5 d^4 e^2Ld':>22} {'P total':>12} {'D total':>12} {'better':>8}")
d = 3.0/8
for L in (0.5,0.8,1.0,1.2,1.5,1.7,2.0,2.5,3.0,4.0,6.0,8.0):
    main, rem = cL(L)*d*d, (16.0/3)*L**5*d**4*math.exp(2*L*d)
    p, dd = P(L,d), D(L,d)
    print(f"{L:>5} {main:>12.4f} {rem:>22.4f} {p:>12.4f} {dd:>12.4f} {('P' if p<dd else 'D'):>8}")

print()
print("=== 2. Where does the delta^4 remainder overtake the delta^2 main term, at delta=3/8? ===")
lo, hi = 0.1, 10.0
for _ in range(200):
    mid = (lo+hi)/2
    if (16.0/3)*mid**5*d**4*math.exp(2*mid*d) < cL(mid)*d*d: lo = mid
    else: hi = mid
print(f"  remainder = main term at L = {lo:.4f}  (delta = 3/8); beyond this the")
print(f"  Thm-7 bound is dominated by its own O(delta^4) allowance.")

print()
print("=== 3. Where does the direct bound (D) beat the paper bound (P), at delta=3/8? ===")
lo, hi = 0.1, 20.0
for _ in range(300):
    mid=(lo+hi)/2
    if P(mid,d) < D(mid,d): lo=mid
    else: hi=mid
print(f"  crossover at L = {lo:.4f}: for L > {lo:.4f} the elementary bound (D) is sharper.")

print()
print("=== 4. What the 7/8 theorem actually buys: delta<=3/8 instead of the trivial delta<=1/2 ===")
print(f"{'L':>5} {'D at d=1/2 (trivial)':>22} {'D at d=3/8 (7/8 thm)':>22} {'gain factor':>12}")
for L in (1.0,2.0,3.0,4.0,6.0,8.0,12.0):
    a, b = D(L,0.5), D(L,0.375)
    print(f"{L:>5} {a:>22.3f} {b:>22.3f} {a/b:>12.3f}")
print()
print("  asymptotically the gain factor ~ (0.5/0.375) e^{(1-0.75)L} = 1.333 e^{L/4}:")
for L in (8.0,12.0,16.0):
    print(f"    L={L:>4}  measured {D(L,0.5)/D(L,0.375):>10.3f}   predicted {1.3333*math.exp(L/4):>10.3f}")
