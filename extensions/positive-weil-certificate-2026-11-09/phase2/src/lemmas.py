"""Lemmas A and B: local-energy control for window-supported functions.

STATUS: standard estimates / derived tools.  NOT novel arithmetic results.
They are elementary consequences of Cauchy-Schwarz, Parseval and a Bernstein-type
derivative bound, and are recorded here because they CORRECT an overbroad claim
made in the first build (see ../RETRACTIONS.md R6).

Fourier convention, fixed once and used everywhere:

    F(t) = int_{-L}^{L} f(u) e^{i t u} du ,      int_R |F|^2 dt = 2 pi ||f||_2^2 .

--------------------------------------------------------------------- LEMMA A
For measurable E with meas(E) < infinity and any f in L^2[-L,L],

    int_E |F|^2  <=  min(1, L meas(E)/pi) * int_R |F|^2 .

Proof.  By Cauchy-Schwarz, |F(t)|^2 <= (int_{-L}^L 1)(int_{-L}^L |f|^2)
= 2L ||f||_2^2 for every t.  Hence int_E |F|^2 <= meas(E) * 2L ||f||_2^2.  By
Parseval ||f||_2^2 = (1/2pi) int_R |F|^2, so the right side is
(L meas(E)/pi) int_R |F|^2.  The bound by 1 is trivial.  []

--------------------------------------------------------------------- LEMMA B
Let h > 0 and 0 <= eps <= 1 satisfy  meas(E ∩ I) <= eps h  for every cell I of
ONE fixed partition of R into intervals of length h.  Then

    int_E |F|^2  <=  q * int_R |F|^2 ,      q = min(1, eps (1 + 2 L h)) .

Proof.  Put G = |F|^2.  Since |u| <= L on supp f, both f and u f(u) lie in
L^2[-L,L]; so F and F' = FT[i u f] lie in L^2(R), F is absolutely continuous
(indeed entire of exponential type L), and G is C^1 with G' = 2 Re(F' conj F).
On a cell I of length h, for the point t* where G attains its max on the closed
cell, G(t*) <= G(t) + int_I |G'| for every t in I; integrating in t over I,

    h sup_I G  <=  int_I G + h int_I |G'| .

Therefore int_{E ∩ I} G <= meas(E ∩ I) sup_I G <= eps h sup_I G
<= eps ( int_I G + h int_I |G'| ).  Summing over the cells,

    int_E G  <=  eps ( int_R G + h int_R |G'| ) .

Finally int_R |G'| <= 2 int_R |F'||F| <= 2 ||F'||_2 ||F||_2, and by Parseval
||F'||_2^2 = 2 pi || u f ||_2^2 <= 2 pi L^2 ||f||_2^2 = L^2 ||F||_2^2, so
||F'||_2 <= L ||F||_2 and int_R |G'| <= 2 L int_R G.  Hence
int_E G <= eps (1 + 2 L h) int_R G.  The bound by 1 is trivial.  []

REMARK (why this is the correction).  The first build asserted that "no bound on
meas(E) can ever control int_E |F|^2".  That is FALSE: Lemma A does exactly that
from finite measure alone, and Lemma B does it from a UNIFORM LOCAL occupation
bound.  What is true is the weaker statement that *global asymptotic density*
supplies no uniform control, because density says nothing about where the mass
sits at any finite scale.  The hypothesis "for every real a" is stronger than
needed for the proof (one fixed partition suffices) but is the robust form, and
is the shape of the uniform-in-translate windowed bounds in Curry's 3x+1 note.
"""
import math, sys
from flint import arb, acb, ctx

def q_lemma_A(L, measE):
    """min(1, L*meas(E)/pi) as an Arb ball, rounded OUTWARD (upper)."""
    L, measE = arb(L), arb(measE)
    v = L * measE / arb.pi()
    return arb(1) if (v > 1) else v

def q_lemma_B(L, h, eps):
    """min(1, eps*(1+2*L*h)) as an Arb ball, rounded OUTWARD (upper)."""
    L, h, eps = arb(L), arb(h), arb(eps)
    v = eps * (1 + 2 * L * h)
    return arb(1) if (v > 1) else v

# ------------------------------------------------------------------ numerical tests
def _F(f_coeffs, L, t):
    """F(t) for f(u) = sum_k c_k cos(k pi u / L) on [-L,L]  (real, even)."""
    tot = 0.0
    for k, c in enumerate(f_coeffs):
        w = k * math.pi / L
        if abs(abs(t) - w) < 1e-12:
            tot += c * L
        else:
            a = math.sin((w - t) * L) / ((w - t) * L) if abs(w - t) > 1e-300 else 1.0
            b = math.sin((w + t) * L) / ((w + t) * L) if abs(w + t) > 1e-300 else 1.0
            tot += c * L * (a + b)
    return tot

def _norm2(f_coeffs, L):
    s = 0.0
    for k, c in enumerate(f_coeffs):
        s += c * c * (2 * L if k == 0 else L)
    return s

def test_lemmas(L=0.8, trials=6, M=240000, Tmax=400.0):
    """Check A and B numerically on explicit f and explicit sets E.  NUMERICAL."""
    import random
    random.seed(20261109)
    out = []
    for tr in range(trials):
        n = random.randint(1, 6)
        c = [random.uniform(-1, 1) for _ in range(n)]
        nrm2 = _norm2(c, L)
        total = 2 * math.pi * nrm2                     # int_R |F|^2
        # --- E = a union of random intervals inside [-Tmax,Tmax]
        h = random.choice([0.5, 1.0, 2.0])
        eps = random.choice([0.05, 0.1, 0.25])
        # E occupies EXACTLY eps*h of every cell, placed at a random offset in the
        # cell -- so the uniform local occupation hypothesis of Lemma B holds with
        # this eps, and q_B is genuinely < 1.
        off = random.uniform(0.0, 1.0 - eps)
        dt = 2 * Tmax / M
        intE = 0.0
        measE = 0.0
        for i in range(M):
            t = -Tmax + (i + 0.5) * dt
            frac_in_cell = ((t + Tmax) % h) / h
            if off <= frac_in_cell < off + eps:
                intE += _F(c, L, t) ** 2 * dt
                measE += dt
        qA = float(q_lemma_A(L, measE).mid())
        qB = float(q_lemma_B(L, h, eps).mid())
        okA = intE <= qA * total * (1 + 1e-9)
        okB = intE <= qB * total * (1 + 1e-9)
        out.append(dict(trial=tr, h=h, eps=eps, measE=measE, ratio=intE / total,
                        qA=qA, qB=qB, okA=okA, okB=okB))
    return out

def main():
    ctx.prec = 200
    print("=== LEMMAS A and B (standard estimates; correction to RETRACTIONS R6) ===")
    print("convention: F(t) = int_{-L}^{L} f e^{itu} du,  int_R|F|^2 = 2 pi ||f||^2\n")
    L = 0.8
    print("Lemma A, L = 0.8:")
    for m in (0.01, 0.1, 1.0, 4.0):
        print("   meas(E) = %-5s -> q_A = %s" % (m, q_lemma_A(L, m).str(12)))
    print("Lemma B, L = 0.8:")
    for h, e in ((0.1, 0.01), (1.0, 0.01), (1.0, 0.1), (10.0, 0.001)):
        print("   h = %-5s eps = %-6s -> q_B = %s   (eps*(1+2Lh) = %s)"
              % (h, e, q_lemma_B(L, h, e).str(12), (arb(e) * (1 + 2 * arb(L) * arb(h))).str(12)))
    print("\nnumerical consistency tests (NUMERICAL, not proof):")
    bad = 0
    for r in test_lemmas(L):
        flag = "ok" if (r["okA"] and r["okB"]) else "VIOLATION"
        if flag != "ok":
            bad += 1
        print("   trial %d  h=%.1f eps=%.2f meas(E)=%6.1f  actual=%.6f  q_A=%.4f q_B=%.4f  "
              "B informative:%s  %s"
              % (r["trial"], r["h"], r["eps"], r["measE"], r["ratio"], r["qA"], r["qB"],
                 "yes" if r["qB"] < 1 else "no", flag))
    print("\n%d violations" % bad)
    return 1 if bad else 0

if __name__ == "__main__":
    sys.exit(main())
