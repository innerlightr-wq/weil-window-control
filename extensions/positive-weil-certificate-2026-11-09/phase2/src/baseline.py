"""BASELINE: the finite-to-infinite comparison at L = 4/5, in the repository
normalization, derived here rather than imported.

------------------------------------------------------------------ THE CHAIN

Normalization (verified against src/zeta_window.py in phase 1):

    F(t) = int_{-L}^{L} f(u) e^{itu} du ,    int_R |F|^2 = 2 pi ||f||_2^2 ,
    Q(f) = 2 s F(i/2)^2 + (1/2pi) int_R |F(t)|^2 Psi_L(t) dt ,
    Psi_L(t) = Re psi(1/4 + it/2) - log pi - P_L(t) ,
    P_L(t)   = sum_{log n < 2L} (2 Lambda(n)/sqrt n) cos(t log n) ,
    s = +1 in the even sectors, s = -1 in the odd sector (zeta_window.pole_sign).

STEP 1 (exterior gate).  H(t) := Re psi(1/4 + it/2) - log pi is strictly
increasing on t > 0 (phase 1, G1-G3), and P_L <= A_L.  So with

    beta(T) := H(T) - A_L ,        Psi_L(t) >= beta(T)  for all |t| >= T .

STEP 2 (split).  For any T > 0, since int_{|t|>=T}(Psi_L - beta)|F|^2 >= 0,

    Q(f) >= R_T(f) := beta ||f||_2^2 + 2 s F(i/2)^2
                      + (1/2pi) int_{-T}^{T} (Psi_L(t) - beta) |F(t)|^2 dt .

STEP 3 (basis).  Legendre, because its Fourier transforms are spherical Bessel
functions and these are super-exponentially small past the Shannon number:

    phi_n(u) = sqrt((2n+1)/(2L)) P_n(u/L) ,   F_n(t) = i^n sqrt(2L(2n+1)) j_n(tL) ,
    sum_{n>=0} (2n+1) j_n(x)^2 = 1   (exactly; Bessel closure relation).

STEP 4 (complement).  Let V_N = span{phi_n : n < N} and V_N^perp its complement.
For f in V_N^perp with ||f|| = 1, Cauchy-Schwarz in n gives
|F(t)|^2 <= 2L sum_{n>=N}(2n+1) j_n(tL)^2, so

    g_N := sup_{f in V_N^perp, ||f||=1} (1/2pi) int_{-T}^{T} |F|^2 dt
         <= (2 L T / pi) * S_N(TL) ,     S_N(x) := sum_{n>=N} (2n+1) j_n(x)^2 ,

using that S_N is increasing in x on 0 <= x < N.  (This is the Cauchy-Schwarz
surrogate for the (N+1)-st eigenvalue of the Slepian concentration operator; it
is lossy but rigorous, and super-exponentially small once N > e T L / 2.)

STEP 5 (block bound).  Write M := sup_{|t|<=T} |Psi_L(t) - beta|.  For
x = x_N + x_perp with a = ||x_N||, b = ||x_perp||:

    <R x_perp, x_perp>  >= beta b^2 - M g_N b^2            (the pole term is >= 0
                                                            when s = +1)
    |cross|             <= 2 ( 2 |p_N| |p_perp| + M sqrt(g_N) ) a b

so with  lam_D := beta - M g_N  and  kappa := 2 |p_N| |p_perp| + M sqrt(g_N),

    lam_min(R_T) >= ( lam_A + lam_D - sqrt( (lam_A - lam_D)^2 + 4 kappa^2 ) ) / 2 ,

where lam_A := lam_min of the N x N compression of R_T.  THIS is the only step
that needs the finite matrix.  Everything above is certified here.

STATUS.  Steps 1-5 are PROVED and the constants beta, g_N, lam_D, kappa are
CERTIFIED below.  lam_A is NOT certified in this build: see REPORT phase 2.
The chain therefore yields a CERTIFIED CONDITIONAL statement, not a theorem.
"""
import sys, os, json, math
from flint import arb, acb, ctx

HERE = os.path.dirname(os.path.abspath(__file__))

def set_prec(b):
    ctx.prec = b

def sph_j(n, x):
    """j_n(x) = sqrt(pi/(2x)) J_{n+1/2}(x) as an Arb ball."""
    return (arb.pi() / (2 * x)).sqrt() * acb(x).bessel_j(acb(n) + acb(1) / 2).real

def comb_terms(L_num, L_den):
    """[(n, p, k)] with log n < 2L, by certified comparison."""
    out, n = [], 2
    twoL = 2 * arb(L_num) / arb(L_den)
    while True:
        lg = arb(n).log()
        if lg > twoL:
            break
        m, p = n, 2
        lam = None
        while p * p <= m:
            if m % p == 0:
                q = m
                while q % p == 0:
                    q //= p
                lam = p if q == 1 else None
                break
            p += 1
        else:
            lam = m
        if lam is not None:
            out.append((n, lam))
        n += 1
    return out

def A_L(L_num, L_den):
    return sum((2 * arb(p).log() / arb(n).sqrt() for n, p in comb_terms(L_num, L_den)), arb(0))

def H(t):
    """H(t) = Re psi(1/4 + it/2) - log pi, strictly increasing on t > 0."""
    return acb(arb(1) / 4, arb(t) / 2).digamma().real - arb.pi().log()

def Psi(t, L_num, L_den):
    c = arb(0)
    for n, p in comb_terms(L_num, L_den):
        c += 2 * arb(p).log() / arb(n).sqrt() * (arb(t) * arb(n).log()).cos()
    return H(t) - c

def S_tail(N, x, top=None):
    """Upper enclosure for S_N(x) = sum_{n>=N}(2n+1) j_n(x)^2.

    Terms are summed explicitly to `top`; beyond `top` the ratio test is used:
    for n > x the sequence (2n+1)j_n(x)^2 is decreasing by a factor < (x/(2n))^2
    per step (from j_{n+1}(x) <= x j_n(x)/(2n+3)), so the remainder is bounded by
    a geometric series.  Everything is an Arb ball.
    """
    if top is None:
        top = N + 160
    tot = arb(0)
    last = None
    for n in range(N, top + 1):
        v = sph_j(n, x)
        t = (2 * n + 1) * v * v
        tot += t
        last = t
    # geometric remainder: ratio bound r = (x/(2*top+3))^2 * (2top+3)/(2top+1) < 1
    r = (x / (2 * top + 3)) ** 2 * arb(2 * top + 3) / arb(2 * top + 1)
    if not (r < arb(1) / 2):
        raise RuntimeError("tail ratio not small enough; raise `top`")
    tot += last * r / (1 - r)
    return tot

def M_sup(L_num, L_den, beta, T):
    """Rigorous bound for sup_{|t|<=T} |Psi_L(t) - beta|, WITHOUT any covering.

    H(t) = Re psi(1/4+it/2) - log pi is strictly increasing on t > 0 (G1-G3) and
    even, so on |t| <= T it lies in [H(0), H(T)].  And |P_L| <= A_L.  Hence

        H(0) - A_L  <=  Psi_L(t)  <=  H(T) + A_L      for |t| <= T,

    and M = max( |H(0) - A_L - beta| , |H(T) + A_L - beta| ).

    An earlier version covered [0,T] with 3000 balls and evaluated digamma on
    each; that was both lossy and unsound -- digamma on a wide ball returned nan
    on some cells and the comparison silently skipped them, inflating M to 1831
    instead of 8.84.  Monotonicity removes the need for any covering.
    """
    A = A_L(L_num, L_den)
    lo = H(0) - A - beta
    hi = H(T) + A - beta
    m = lo.abs_upper()
    m2 = hi.abs_upper()
    return arb(m if arb(m) > arb(m2) else m2)

def pole_vec_tail(N, L, top=None):
    """|p_perp|^2 = sum_{n>=N} p_n^2 with p_n = sqrt(2L(2n+1)) * i_n(L/2) * (-1)^n,
    i_n the modified spherical Bessel function.  Returns an upper enclosure."""
    if top is None:
        top = N + 80
    tot = arb(0)
    z = L / 2
    for n in range(N, top + 1):
        # i_n(z) = sqrt(pi/(2z)) I_{n+1/2}(z)
        iv = (arb.pi() / (2 * z)).sqrt() * acb(z).bessel_i(acb(n) + acb(1) / 2).real
        tot += 2 * L * (2 * n + 1) * iv * iv
    return tot

def report(L_num=4, L_den=5, T=200, Ns=(180, 200, 220, 250, 300), prec=1600):
    set_prec(prec)
    L = arb(L_num) / arb(L_den)
    A = A_L(L_num, L_den)
    beta = H(T) - A
    print("L = %d/%d   T = %d" % (L_num, L_den, T))
    print("  A_L    = %s" % A.str(14))
    print("  beta   = H(T) - A_L = %s     (> 0 required)" % beta.str(14))
    assert beta > 0
    Msup = M_sup(L_num, L_den, beta, T)
    print("  M      = sup_{|t|<=T} |Psi_L - beta| <= %s   (monotonicity of H, no covering)"
          % Msup.str(12))
    pnorm = (2 * L.sinh()).sqrt()
    print("  |p|    = sqrt(2 sinh L) = %s   (exact: sum_n p_n^2 = int_{-L}^L e^{-u} du)"
          % pnorm.str(12))
    x = arb(T) * L
    print("  T*L    = %s   (Shannon number 2TL/pi = %s)"
          % (x.str(8), (2 * x / arb.pi()).str(8)))
    rows = []
    for N in Ns:
        S = S_tail(N, x)
        g = (2 * L * arb(T) / arb.pi()) * S
        lamD = beta - Msup * g
        pt = pole_vec_tail(N, L).sqrt()
        kappa = 2 * pnorm * pt + Msup * g.sqrt()          # |p_N| <= |p| = sqrt(2 sinh L)
        print("  N=%-4d S_N=%-15s g_N=%-15s |p_perp|=%-13s lam_D=%-13s kappa=%s"
              % (N, S.str(7), g.str(7), pt.str(7), lamD.str(9), kappa.str(7)))
        # required lam_A for the Schur bound to give a positive lower bound:
        # need kappa^2 <= lam_A * lam_D  =>  lam_A >= kappa^2 / lam_D
        need = kappa * kappa / lamD
        print("        -> requires lam_A >= kappa^2/lam_D = %s" % need.str(8))
        rows.append(dict(N=N, S=S.str(20), g=g.str(20), lamD=lamD.str(20),
                         kappa=kappa.str(20), lamA_needed=need.str(20)))
    return dict(L="%d/%d" % (L_num, L_den), T=T, A_L=A.str(20), beta=beta.str(20),
                M=Msup.str(20), rows=rows)

if __name__ == "__main__":
    out = report()
    d = os.path.join(HERE, "..", "certs")
    os.makedirs(d, exist_ok=True)
    json.dump(out, open(os.path.join(d, "baseline_constants.json"), "w"), indent=2)
    print("\n-> certs/baseline_constants.json")
