"""Phase 3 audit: beta monotonicity, the three spectral quantities, and the
Cauchy-Schwarz provenance of the uniform transform bound.

------------------------------------------------- 1. BETA MONOTONICITY (PROVED)

For fixed T the complete comparison form is

    R_{T,beta}(f) = pole(f) + beta ||f||_2^2
                    + (1/2pi) int_{|t|<=T} (Psi_L(t) - beta)|F(t)|^2 dt .

By Parseval, (1/2pi) int_R |F|^2 = ||f||_2^2, so

    beta ||f||^2 - (beta/2pi) int_{|t|<=T} |F|^2 = (beta/2pi) int_{|t|>T} |F|^2 ,

and therefore the EXACT equivalent expression

    R_{T,beta}(f) = pole(f) + (1/2pi) int_{|t|<=T} Psi_L |F|^2
                            + (beta/2pi) int_{|t|>T} |F|^2 .

Hence for beta_2 >= beta_1 and every f,

    R_{T,beta_2}(f) - R_{T,beta_1}(f)
        = ((beta_2 - beta_1)/2pi) int_{|t|>T} |F|^2 dt   >=  0 .          (*)

So R is MONOTONE NON-DECREASING in beta, and so are its spectral infimum and the
least eigenvalue of its compression to any fixed subspace.  **Lowering beta
cannot increase lambda_A.**  Looking only at `Psi_L - beta` inside the integral
and concluding that a smaller beta makes the interior "less negative" misses the
`beta ||f||^2` term, which moves the other way and dominates by exactly (*).

FINITE-MATRIX COUNTERPART.  In the orthonormal basis of one sector,

    A(beta)_{nm} = beta delta_{nm} + 2 s p_n p_m
                   + (1/2pi) int_{|t|<=T} (Psi_L - beta) v_n v_m dt ,

so the beta-dependence is  beta delta_{nm} - beta (K_T)_{nm}  with

    (K_T)_{nm} := (1/2pi) int_{-T}^{T} v_n(t) v_m(t) dt ,

giving exactly

    A(beta_2) - A(beta_1) = (beta_2 - beta_1) (I - K_T) .                  (**)

K_T is the compression of the Fourier energy-concentration operator:
<K_T y, y> = (1/2pi) int_{-T}^{T} |F|^2 dt with F = sum y_n v_n.  Therefore

    0 <= <K_T y, y> <= (1/2pi) int_R |F|^2 = ||f||_2^2 = ||y||^2 ,

i.e. 0 <= K_T <= I, so (**) is positive semidefinite for beta_2 >= beta_1.

CONSEQUENCE FOR THIS PROJECT.  The "next step" proposed at the end of the phase-3
run -- lower beta to buy a larger certified a -- is REFUTED by (**).  Any
improvement seen after lowering beta would have to come from enclosure quality or
another changed parameter, not from an operator improvement.

NOTE ON SUFFICIENCY.  beta = H(T) - A_L is a SUFFICIENT exterior bound
(Psi_L >= beta on |t| >= T).  It is not automatically the exact infimum of Psi_L
on the exterior, and nothing here claims it is optimal.

---------------------------------------- 2. THREE SPECTRAL QUANTITIES (distinct)

    lambda_Q = inf { Q(f)/||f||^2 }        over the admissible whole window
    lambda_R = inf { R_{T,beta}(f)/||f||^2 } over the same window
    lambda_A = min eigenvalue of the finite compression of R in ONE sector

Because Q >= R pointwise (the discarded exterior integrand is >= 0 by the gate),
and because an infimum over a subspace is >= the infimum over the whole space,

    lambda_R <= lambda_Q        and        lambda_R <= lambda_A .

There is NO inequality between lambda_A and lambda_Q in general.  In particular
the numerical even-block value 1.0310220e-17 estimates **lambda_A**, and is NOT
an upper bound for lambda_Q.  A certified upper bound for lambda_Q needs a
variational evaluation of Q itself on some admissible f.

-------------------------------- 6. THE UNIFORM TRANSFORM BOUND IS CAUCHY-SCHWARZ

For any phi supported in [-L,L] with ||phi||_2 = 1 and any complex z,

    | int phi(u) e^{i z u} du |
        <= ||phi||_2 ( int_{-L}^{L} |e^{izu}|^2 du )^{1/2}
         = ( int_{-L}^{L} e^{-2 u Im z} du )^{1/2}
         = ( sinh(2 L Im z) / Im z )^{1/2}
        <= sqrt(2L) e^{L |Im z|} .

That is the whole proof.  The bound used in the assembly is this, specialised to
the orthonormal Legendre basis; it is an effective implementation choice, NOT a
new inequality, and the earlier phrasing ("from the Legendre integral
representation") should not be read as a novelty claim.
"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "src"))
from flint import arb, acb, ctx
import block as BK

CHECKS = []
def chk(name, ok, extra=""):
    CHECKS.append((name, bool(ok), extra))
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name, ("  -- " + extra) if extra else ""))
    return ok

def K_T(S, N, T, panels, korder):
    """(K_T)_{nm} = (1/2pi) int_{-T}^{T} v_n v_m dt = (1/pi) int_0^T v_n v_m dt,
    on the EVEN indices (same parity decoupling as the main assembly)."""
    L = S.L
    nodes, wts = BK.gl_rule(korder)
    hw = arb(T) / (2 * panels)
    ev = [n for n in range(N) if n % 2 == 0]
    ne = len(ev)
    K = [[arb(0)] * ne for _ in range(ne)]
    for ip in range(panels):
        c = hw * (2 * ip + 1)
        for q in range(korder):
            t = c + hw * nodes[q]
            vv = BK.v_real(L, t, N)
            w = wts[q] * hw / arb.pi()
            for a in range(ne):
                for b in range(a, ne):
                    K[a][b] = K[b][a] = K[a][b] + w * vv[ev[a]] * vv[ev[b]]
    return K

def A_of_beta(S, beta, N, T, panels, korder):
    """Even-sector A(beta) at a GIVEN beta (overriding S.beta), small N."""
    L = S.L
    nodes, wts = BK.gl_rule(korder)
    hw = arb(T) / (2 * panels)
    ev = [n for n in range(N) if n % 2 == 0]
    ne = len(ev)
    C = [[arb(0)] * ne for _ in range(ne)]
    for ip in range(panels):
        c = hw * (2 * ip + 1)
        for q in range(korder):
            t = c + hw * nodes[q]
            W = S.W_real_at(t, beta)
            vv = BK.v_real(L, t, N)
            w = wts[q] * hw * W / arb.pi()
            for a in range(ne):
                for b in range(a, ne):
                    C[a][b] = C[b][a] = C[a][b] + w * vv[ev[a]] * vv[ev[b]]
    p = BK.pole_vector(L, N)
    pe = [p[n] for n in ev]
    return [[C[a][b] + 2 * pe[a] * pe[b] + (beta if a == b else arb(0))
             for b in range(ne)] for a in range(ne)]


def main():
    ctx.prec = 700
    T, N, PAN, KO = 200, 16, 200, 32
    S = BK.Symbol(4, 5, T)
    L = S.L
    b1 = arb("0.45")
    b2 = arb("0.55")
    print("=== PHASE 3 AUDIT ===")
    print("-- M1 the finite-matrix beta identity  A(b2)-A(b1) = (b2-b1)(I-K_T)")
    A1 = A_of_beta(S, b1, N, T, PAN, KO)
    A2 = A_of_beta(S, b2, N, T, PAN, KO)
    K = K_T(S, N, T, PAN, KO)
    n = len(K)
    worst = arb(0)
    for a in range(n):
        for b in range(n):
            lhs = A2[a][b] - A1[a][b]
            rhs = (b2 - b1) * ((arb(1) if a == b else arb(0)) - K[a][b])
            d = (lhs - rhs).abs_upper()
            if arb(d) > worst:
                worst = arb(d)
    chk("M1 identity holds entrywise (N=%d, %d panels, order %d)" % (N, PAN, KO),
        bool(worst < arb(10) ** -40), "max |LHS-RHS| = %s" % worst.str(8))

    print("-- M2  0 <= K_T <= I")
    # diagonal check plus a Rayleigh sweep on frozen rational vectors
    diag_ok = all(bool(K[a][a] >= 0) and bool(K[a][a] <= 1) for a in range(n))
    chk("M2a every diagonal entry of K_T lies in [0,1]", diag_ok,
        "min %s max %s" % (min(K[a][a] for a in range(n)).str(8),
                           max(K[a][a] for a in range(n)).str(8)))
    import random
    random.seed(7)
    okA = okB = True
    worst_lo, worst_hi = arb(1), arb(0)
    for trial in range(12):
        y = [arb(random.randint(-9, 9)) for _ in range(n)]
        nn = sum((v * v for v in y), arb(0))
        if not (nn > 0):
            continue
        q = sum((y[a] * K[a][b] * y[b] for a in range(n) for b in range(n)), arb(0)) / nn
        okA &= bool(q >= 0)
        okB &= bool(q <= 1)
        if q < worst_lo:
            worst_lo = q
        if q > worst_hi:
            worst_hi = q
    chk("M2b Rayleigh quotients of K_T lie in [0,1] on 12 frozen integer vectors",
        okA and okB, "range [%s, %s]" % (worst_lo.str(8), worst_hi.str(8)))

    print("-- M3 monotonicity consequence: A(b2) - A(b1) is PSD")
    D = [[A2[a][b] - A1[a][b] for b in range(n)] for a in range(n)]
    ok, res, err = BK.verify_pd_preconditioned(D, arb(0), 700)
    chk("M3 A(0.55) - A(0.45) is positive definite (so lam_min is non-decreasing in beta)",
        ok, "||I-G||_inf = %s" % ("n/a" if res is None else res.str(8)))

    print("-- M4 the retracted optimisation: lowering beta cannot raise lam_A")
    okd, resd, _ = BK.verify_pd_preconditioned(
        [[A1[a][b] - A2[a][b] for b in range(n)] for a in range(n)], arb(0), 700)
    chk("M4 the REVERSE difference A(0.45)-A(0.55) is NOT positive definite", not okd,
        "as required by (**)")

    print("-- M5 Cauchy-Schwarz provenance of the uniform transform bound")
    b_e = arb("0.30")
    cs = (2 * L).sqrt() * (L * b_e).exp()
    used = BK.v_hol_absmax_uniform(L, b_e)
    chk("M5 the bound used equals sqrt(2L) exp(L b_e) from Cauchy-Schwarz",
        bool((cs - used).abs_upper() < arb(10) ** -40), used.str(12))

    out = dict(T=T, N=N, panels=PAN, korder=KO,
               beta_1=b1.str(12), beta_2=b2.str(12),
               identity_max_residual=worst.str(20),
               K_diag_min=min(K[a][a] for a in range(n)).str(20),
               K_diag_max=max(K[a][a] for a in range(n)).str(20),
               uniform_bound=used.str(20))
    d = os.path.join(HERE, "..", "certs")
    os.makedirs(d, exist_ok=True)
    json.dump(out, open(os.path.join(d, "beta_monotonicity.json"), "w"), indent=2)
    bad = [x for x, ok, _ in CHECKS if not ok]
    print("\nAUDIT: %d checks, %d failures" % (len(CHECKS), len(bad)))
    return 1 if bad else 0

if __name__ == "__main__":
    sys.exit(main())
