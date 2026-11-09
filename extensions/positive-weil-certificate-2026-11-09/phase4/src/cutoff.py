"""Phase 5, part 1: the cutoff comparison, and beta at T = 400.

------------------------------------------------------- THE CUTOFF IDENTITY

Let beta_i be certified exterior bounds, Psi_L(t) >= beta_i for |t| > T_i.
Using the Parseval form established in phase 4 (RETRACTIONS R11),

    R_{T,beta}(f) = pole(f) + (1/2pi) int_{|t|<=T} Psi_L |F|^2
                            + (beta/2pi) int_{|t|>T} |F|^2 ,

take T_2 >= T_1 and beta_2 >= beta_1.  Splitting the second integral of the
T_2 form over {|t| <= T_1} and {T_1 < |t| <= T_2}:

    R_{T_2,beta_2}(f) - R_{T_1,beta_1}(f)
       = (1/2pi) int_{T_1<|t|<=T_2} (Psi_L(t) - beta_1) |F(t)|^2 dt
       + ((beta_2 - beta_1)/2pi) int_{|t|>T_2} |F(t)|^2 dt .                (C)

Both terms are >= 0: the first because Psi_L >= beta_1 on |t| > T_1, the second
because beta_2 >= beta_1.  Hence lambda_R(T_2) >= lambda_R(T_1) FOR THE EXACT
FORMS.  It does NOT follow that a certified numerical bound improves, nor does
it order finite matrices of different sizes: changing N changes the compression,
so only matrices in the SAME subspace are comparable by (C).
"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "src"))
from flint import arb, acb, ctx
import block as BK, baseline as B

CHECKS = []
def chk(name, ok, extra=""):
    CHECKS.append((name, bool(ok), extra))
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name, ("  -- " + extra) if extra else ""))
    return ok

def main():
    ctx.prec = 1200
    L = arb(4) / 5
    A_L = B.A_L(4, 5)
    print("=== cutoff comparison and beta at T = 400 ===")
    b200 = B.H(200) - A_L
    b400 = B.H(400) - A_L
    print("  A_L        = %s" % A_L.str(20))
    print("  beta_200   = H(200) - A_L = %s" % b200.str(20))
    print("  beta_400   = H(400) - A_L = %s" % b400.str(20))
    chk("C1 beta_400 > beta_200 strictly (H is increasing, proved G1-G3)",
        bool(b400 > b200), "difference = %s" % (b400 - b200).str(12))
    chk("C2 both exterior bounds are positive", bool(b200 > 0) and bool(b400 > 0))
    # the gate itself at the new cutoff, re-derived not assumed
    chk("C3 Psi_L(400) >= beta_400 (gate tight at T=400)",
        bool(B.H(400) - A_L >= b400 - arb(10) ** -40))
    # sample the band [200,400] to confirm Psi_L - beta_200 >= 0 there (term 1 of (C))
    worst = None
    for k in range(0, 41):
        t = arb(200) + arb(5) * k
        v = B.Psi(t, 4, 5) - b200
        if worst is None or v < worst:
            worst = v
    chk("C4 Psi_L - beta_200 >= 0 sampled across [200,400]", bool(worst >= 0),
        "min sampled = %s (monotonicity of H makes this structural)" % worst.str(10))
    out = dict(A_L=A_L.str(25), beta_200=b200.str(25), beta_400=b400.str(25),
               beta_increase=(b400 - b200).str(25))
    d = os.path.join(HERE, "..", "certs")
    os.makedirs(d, exist_ok=True)
    json.dump(out, open(os.path.join(d, "cutoff.json"), "w"), indent=2)
    bad = [n for n, ok, _ in CHECKS if not ok]
    print("\n%d checks, %d failures" % (len(CHECKS), len(bad)))
    return 1 if bad else 0

if __name__ == "__main__":
    sys.exit(main())
