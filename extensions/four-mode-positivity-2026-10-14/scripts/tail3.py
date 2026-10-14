"""Third-order archimedean tail  O(1/M^3)  and a sharper gamma_E interval T.

Route O needs the error allowance BELOW lambda_min(G) ~ 1e-10.  The O(1/M^2)
tail of tail_sharp.py gives 1.6e-9 per unit ||x||^2, and the published gamma_E
interval T has width 3.9e-11 per diagonal entry.  Both are sharpened here.

-------------------------------------------------------------------- 1. T
T := sum_{k>M} [ 1/k - log(1+1/k) ].  For 0 < x <= 1,
   x - log(1+x) = sum_{n>=2} (-1)^n x^n / n
is alternating with decreasing terms, so
   x^2/2 - x^3/3  <=  x - log(1+x)  <=  x^2/2 - x^3/3 + x^4/4 .
Hence with P_j := sum_{k>M} k^{-j},
   P_2/2 - P_3/3  <=  T  <=  P_2/2 - P_3/3 + P_4/4 ,
and, by midpoint/integral comparison (all proved below),
   P_2 in [ 1/(M+1/2) - 1/(9 M^3) , 1/(M+1/2) ] ,
   P_3 in [ 1/(2(M+1)^2) , 1/(2(M+1/2)^2) ] ,
   P_4 in [ 1/(3(M+1)^3) , 1/(3 M^3) ] .
Proofs.  1/k^2 <= int_{k-1/2}^{k+1/2} x^-2 dx = 1/(k^2-1/4), so
P_2 <= int_{M+1/2}^inf = 1/(M+1/2); and the defect is
sum (1/4)/(k^2(k^2-1/4)) <= (1/3) P_4 <= 1/(9M^3) since k^2-1/4 >= (3/4)k^2.
x^-3 is convex so the midpoint rule underestimates: 1/k^3 <=
int_{k-1/2}^{k+1/2} x^-3 dx, giving P_3 <= 1/(2(M+1/2)^2); and 1/k^3 >=
int_k^{k+1} x^-3 dx gives P_3 >= 1/(2(M+1)^2).  Same for P_4 with x^-4.
Width of T is then O(1/M^3) ~ 3e-16 instead of 2/(3M^2) ~ 3.9e-11.

------------------------------------------------------- 2. the O(1/M^3) tail
Three integrations by parts of I_m = int_0^U e^{-s u} C(u) du, using C(U) = 0:
   I_m = C(0)/s + C'(0)/s^2 + C''(0)/s^3
         - e^{-sU}[ C'(U)/s^2 + C''(U)/s^3 ]  +  (1/s^3) int_0^U e^{-su} C''' du .
With 2C(0)/s = C(0)/(m+1/4) and s = s_m = 2m+1/2,
   t_m = -(3/4)C(0)/((m+1)(m+1/4)) - 2C'(0)/s^2 - 2C''(0)/s^3 + e_m ,
   |e_m| <= 2 sup|C'''| / s^4 + 2 e^{-sU}( |C'(U)|/s^2 + |C''(U)|/s^3 ) .
For this basis BOTH endpoint terms vanish identically:
   sin(w_k U) = 0  =>  C''(U) = 0 ,     cos(w_k U) = -1  =>  C'(U) = -C'(0) = 0 ,
and C'(0) = 0 because C is even and C^1.  Moreover
   C''(0) = -||f'||^2 = - sum_k w_k^2 x_k^2        (an EXACT diagonal form,
because int phi_j' phi_k' = w_k^2 delta_jk for this basis).
So
   sum_{m>=M} t_m = -(3/4) C(0) S_1 + 2 ( sum_k w_k^2 x_k^2 ) S_3 + R' ,
   S_1 = sum_{m>=M} 1/((m+1)(m+1/4)) ,   S_3 = sum_{m>=M} s_m^{-3} ,
   |R'| <= 2 sup|C'''| sum_{m>=M} s_m^{-4} ,
with  S_3 = (1/8) sum_{m>=M}(m+1/4)^{-3} in
      (1/8)[ 1/(2(M+1/4)^2) , (M+1/4)^{-3} + 1/(2(M+1/4)^2) ]  and
      sum_{m>=M} s_m^{-4} = (1/16) sum (m+1/4)^{-4}
                         <= (1/16)[ (M+1/4)^{-4} + 1/(3(M+1/4)^3) ] ,
both by the same integral comparison.  Finally
   sup|C'''| <= sum_k [ |A_k| w_k^3 + |B_k| ( 3 w_k^2 + U w_k^3 ) ] .
"""
import sys, os
from fractions import Fraction as Q
sys.path.insert(0, os.path.dirname(__file__))
import weil_core as C
import tail_sharp as TS
V = C.V
Iv = C.Iv
U, NDIM, M = C.U, C.NDIM, C.M
wk = C.wk


def P_enclosures():
    P2 = Iv(Q(1) / Q(2*M + 1, 2) - Q(1, 9*M**3), Q(1) / Q(2*M + 1, 2))
    P3 = Iv(Q(1, 2) / Q(M + 1)**2, Q(1, 2) / Q(2*M + 1, 2)**2)
    P4 = Iv(Q(1, 3) / Q(M + 1)**3, Q(1, 3) / Q(M)**3)
    return P2, P3, P4


def T_sharp():
    P2, P3, P4 = P_enclosures()
    lo = (P2 * Iv(Q(1, 2)) - P3 * Iv(Q(1, 3)))
    hi = lo + P4 * Iv(Q(1, 4))
    return Iv(lo.lo, hi.hi)


def S1_sharp():
    """S_1 = sum_{m>=M} g(m), g(x) = 1/((x+1)(x+1/4)) = (4/3)[1/(x+1/4) - 1/(x+1)],
    to O(1/M^2) -- an order better than the monotone comparison of tail_sharp.py.

    g is positive, decreasing and CONVEX on [M-1/2, inf) (g'' = (8/3)[(x+1/4)^-3
    - (x+1)^-3] > 0), so by Hermite-Hadamard the midpoint rule underestimates:
        g(m) <= int_{m-1/2}^{m+1/2} g ,
    whence  S_1 <= int_{M-1/2}^inf g = (4/3) log( (M+1/2)/(M-1/4) ) .
    The defect obeys  0 <= int_{m-1/2}^{m+1/2} g - g(m) = g''(xi)/24
        <= (8/3)/(24 (m-1/4)^3) = 1/(9 (m-1/4)^3) ,
    so  S_1 >= (4/3) log((M+1/2)/(M-1/4)) - D,
        D := (1/9) sum_{m>=M} (m-1/4)^{-3} <= (1/9)[ (M-1/4)^{-3} + 1/(2(M-5/4)^2) ] .
    """
    top = Iv(Q(4, 3)) * V.iv_log(Iv(Q(2*M + 1, 2)) / Iv(Q(4*M - 1, 4)))
    D = Q(1, 9) * (Q(1) / Q(4*M - 1, 4)**3 + Q(1, 2) / Q(4*M - 5, 4)**2)
    return Iv(top.lo - D, top.hi)


def S3_enclosure():
    a = Q(4*M + 1, 4)
    return Iv(Q(1, 8) * (Q(1, 2) / a**2), Q(1, 8) * (Q(1) / a**3 + Q(1, 2) / a**2))


def S4_bound():
    a = Q(4*M + 1, 4)
    return Iv(0, Q(1, 16) * (Q(1) / a**4 + Q(1, 3) / a**3))


def dddC_sup(A, B):
    tot = Iv(0)
    for k in range(NDIM):
        w = wk(k)
        aa = Iv(max(abs(A[k].lo), abs(A[k].hi)))
        bb = Iv(max(abs(B[k].lo), abs(B[k].hi)))
        tot = tot + aa * w * w.sq() + bb * (Iv(3) * w.sq() + Iv(U) * w * w.sq())
    return tot


def dC_at_U(A, B):
    return TS.dC_at(A, B, Iv(U))


def ddC_at(A, B, u_iv):
    """C''(u) = sum_k [ -A_k w_k^2 sin + B_k ( 2 w_k sin - (U-u) w_k^2 cos ) ]."""
    tot = Iv(0)
    for k in range(NDIM):
        w = wk(k)
        su = V.iv_sin(w * u_iv); cu = V.iv_cos(w * u_iv)
        tot = tot + (-A[k]) * w.sq() * su + B[k] * (Iv(2) * w * su - (Iv(U) - u_iv) * w.sq() * cu)
    return tot


def arch_core3(x):
    """The x-quadratic part of Arch with the SHARP T: C0*(17 log2 - T_sharp) - 2*sum."""
    C._set_X(x)
    A, B = V.C_canonical()
    C0 = sum(Q(t) * Q(t) for t in x)
    SA, SB = C.arch_sums()
    acc = Iv(0)
    for k in range(NDIM):
        acc = acc + A[k] * SA[k] + B[k] * SB[k]
    return Iv(C0) * (Iv(17) * C.LOG2 - T_sharp()) - Iv(2) * acc, A, B, C0


def Qcore3(x):
    ar, A, B, C0 = arch_core3(x)
    pole = V.pole_term()
    prime = V.prime_term(A, B)
    logpi = C.LOGPI * Iv(C0)
    return dict(core=pole + ar - logpi - prime, A=A, B=B, C0=C0,
                pole=pole, arch_core=ar, logpi=logpi, prime=prime)


def tail3(A, B, C0, x):
    """Enclosure of sum_{m>=M} t_m to O(1/M^3)."""
    S1 = S1_sharp(); S3 = S3_enclosure()
    w2x2 = Iv(0)
    for k in range(NDIM):
        w = wk(k)
        w2x2 = w2x2 + w.sq() * Iv(Q(x[k]) * Q(x[k]))
    main = -(Iv(Q(3, 4)) * Iv(C0) * S1) + Iv(2) * w2x2 * S3
    Rm = Iv(2) * dddC_sup(A, B) * S4_bound()
    Rmag = Iv(0, Rm.hi)
    return main + Iv(-Rmag.hi, Rmag.hi), dict(S1=S1, S3=S3, Rmag=Rmag,
                                              dddC=dddC_sup(A, B), w2x2=w2x2)


def Q_zeta3(x):
    d = Qcore3(x)
    t, info = tail3(d['A'], d['B'], d['C0'], x)
    d['tail3'] = t; d['tail3_info'] = info
    d['Q'] = d['core'] + t
    return d


def Q_add3(x, delta):
    d = Q_zeta3(x)
    Ad, Bd = C.quartet(x, delta)
    return d['Q'] + Iv(4) * Ad.sq() - Iv(4) * Bd.sq(), Ad, Bd, d
