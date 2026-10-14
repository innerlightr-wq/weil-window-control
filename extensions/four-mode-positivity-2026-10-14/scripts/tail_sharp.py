"""A SHARPER archimedean tail: the leading 1/M term is computed instead of bounded.

Published bound (CERTIFICATE_THEOREM.md section 5.2, reused unchanged elsewhere):
    | sum_{m>=M} t_m |  <=  [ (3/4)C0 + (1/2) sup|C'| ] / (M - 3/4)      = O(1/M)
with  t_m = C0/(m+1) - 2 I_m ,  I_m = int_0^U e^{-s_m u} C(u) du ,  s_m = 2m+1/2.

NEW (derived and proved here).  Two integrations by parts, using C(U) = 0:
    I_m = C0/s_m + J_m/s_m ,          J_m = int_0^U e^{-s_m u} C'(u) du ,
    J_m = C'(0)/s_m - e^{-s_m U} C'(U)/s_m + (1/s_m) int_0^U e^{-s_m u} C''(u) du .
Since 2 C0/s_m = C0/(m+1/4),
    t_m = -(3/4) C0 / ((m+1)(m+1/4))  -  2 C'(0)/s_m^2  +  r_m ,
    |r_m| <= 2 ( sup|C''| + |C'(U)| e^{-s_m U} ) / s_m^3 .
Hence
    sum_{m>=M} t_m = -(3/4) C0 S1(M) - 2 C'(0) S2(M) + R ,
    S1(M) = sum_{m>=M} 1/((m+1)(m+1/4)) ,   S2(M) = sum_{m>=M} s_m^{-2} ,
    |R|  <=  sup|C''| / ( 4 (M-3/4)^2 )  +  (tiny e^{-s_M U} part) .

Both S1 and S2 are enclosed to O(1/M^2) by monotone integral comparison:
  the summand of S1 is positive and decreasing, so
    (4/3) log((M+1)/(M+1/4))  <=  S1(M)  <=  (4/3) log(M/(M-3/4)) ,
  using int_a^inf dx/((x+1)(x+1/4)) = (4/3) log((a+1)/(a+1/4));
  likewise  (1/4)/(M+1/4) <= S2(M) <= (1/4)/(M-3/4), from
    sum_{m>=M}(m+1/4)^{-2} in [ 1/(M+1/4), 1/(M-3/4) ].
  For |R|: sum_{m>=M} s_m^{-3} = (1/8) sum (m+1/4)^{-3} <= (1/8)/(2(M-3/4)^2).

So the allowance drops from O(1/M) to O(1/M^2) at NO extra truncation cost.
Everything is enclosed outward; the only new inputs are C'(0), C'(U), sup|C''|,
all exact closed forms in (A,B), and one rational logarithm.
"""
import sys, os
from fractions import Fraction as Q
sys.path.insert(0, os.path.dirname(__file__))
import weil_core as C
V = C.V
Iv, mk = C.Iv, C.mk
U, NDIM, M = C.U, C.NDIM, C.M
wk = C.wk


def dC_at(A, B, u_iv):
    """C'(u) for an interval u in [0,U]:
    C' = sum_k [ A_k w_k cos(w_k u) - B_k cos(w_k u) - B_k (U-u) w_k sin(w_k u) ]."""
    tot = Iv(0)
    for k in range(NDIM):
        w = wk(k)
        cu = V.iv_cos(w * u_iv)
        su = V.iv_sin(w * u_iv)
        tot = tot + A[k] * w * cu - B[k] * cu - B[k] * (Iv(U) - u_iv) * w * su
    return tot


def ddC_sup(A, B):
    """sup_{[0,U]} |C''| <= sum_k [ |A_k| w_k^2 + |B_k| (2 w_k + U w_k^2) ], where
    C'' = sum_k [ -A_k w_k^2 sin(w_k u) + B_k( 2 w_k sin(w_k u) - (U-u) w_k^2 cos(w_k u) ) ]."""
    tot = Iv(0)
    for k in range(NDIM):
        w = wk(k)
        aa = Iv(max(abs(A[k].lo), abs(A[k].hi)))
        bb = Iv(max(abs(B[k].lo), abs(B[k].hi)))
        tot = tot + aa * w.sq() + bb * (Iv(2) * w + Iv(U) * w.sq())
    return tot


def S1_enclosure():
    """sum_{m>=M} 1/((m+1)(m+1/4)), enclosed to O(1/M^2)."""
    lo = Iv(Q(4, 3)) * V.iv_log(Iv(Q(M + 1)) / Iv(Q(4*M + 1, 4)))
    hi = Iv(Q(4, 3)) * V.iv_log(Iv(Q(M)) / Iv(Q(4*M - 3, 4)))
    return Iv(lo.lo, hi.hi)


def S2_enclosure():
    """sum_{m>=M} (2m+1/2)^{-2} = (1/4) sum (m+1/4)^{-2}, enclosed to O(1/M^2)."""
    return Iv(Q(1, 4) / Q(4*M + 1, 4), Q(1, 4) / Q(4*M - 3, 4))


def tail_sharp(A, B, C0):
    """Enclosure of sum_{m>=M} t_m, accurate to O(1/M^2)."""
    dC0 = dC_at(A, B, Iv(0))
    dCU = dC_at(A, B, Iv(U))
    S1, S2 = S1_enclosure(), S2_enclosure()
    main = -(Iv(Q(3, 4)) * Iv(C0) * S1) - Iv(2) * dC0 * S2
    Rmag = ddC_sup(A, B) / (Iv(4) * Iv(Q(4*M - 3, 4)).sq())
    # the e^{-s_M U} piece: s_M U >= 2*131070*1.6, utterly negligible, but bounded
    ex = Iv(0, Q(1, 1 << 200))
    Rmag = Rmag + Iv(2) * Iv(max(abs(dCU.lo), abs(dCU.hi))) * ex
    return main + Iv(-Rmag.hi, Rmag.hi), dict(dC0=dC0, dCU=dCU, S1=S1, S2=S2,
                                              Rmag=Rmag, ddCsup=ddC_sup(A, B))


def Q_zeta_sharp(x):
    """Q_zeta(f_x) with the sharper tail. Same Qcore, new tail only."""
    d = C.Qcore(x)
    t, info = tail_sharp(d['A'], d['B'], d['C0'])
    d['tail_sharp'] = t
    d['tail_info'] = info
    d['Q'] = d['core'] + t
    return d


def Q_add_sharp(x, delta):
    d = Q_zeta_sharp(x)
    Ad, Bd = C.quartet(x, delta)
    return d['Q'] + Iv(4) * Ad.sq() - Iv(4) * Bd.sq(), Ad, Bd, d
