"""LOCAL-OCCUPATION PILOT: what the energy allowance can actually buy at L=4/5.

The penalized comparison the brief asks to consume is

    Q(f) >= (beta - b q) ||f||^2 + 2 s F(i/2)^2
            + (1/2pi) int_{|t|<T} (Psi_L - beta) |F|^2 dt ,

where E is the unfavourable exterior set {|t| >= T : Psi_L < beta},
b >= sup_E (beta - Psi_L), and q bounds int_E |F|^2 / int_R |F|^2 (Lemma A or B).
NOTE the interior subtraction stays beta, not beta - b q.

The allowance is only useful if b*q can be ABSORBED, i.e. if

    b q  <  lam_A - kappa^2/lam_D   ~  the finite-block margin,

which at L = 4/5 is of order 1e-17 (the window minimum itself).  This module
computes the q that would be required and compares it with what Lemmas A and B
can deliver.
"""
import os, sys, json
from flint import arb, acb, ctx
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import baseline as B
from lemmas import q_lemma_A, q_lemma_B

def exterior_set_band(L_num, L_den, T, Told_hi=10**6):
    """E = {|t| >= T : Psi_L(t) < beta} is contained in a FINITE band.

    With beta := H(T) - A_L the gate is tight at T and H is increasing, so
    Psi_L(t) >= H(t) - A_L >= beta for every t >= T: E is EMPTY.
    If instead a SMALLER cutoff T' < T is used with the same beta, then E is
    contained in [T', T_old] where T_old is the least t with H(t) - A_L >= beta,
    i.e. T_old = T.  This function returns that band for a trial cutoff.
    """
    A = B.A_L(L_num, L_den)
    beta = B.H(T) - A
    return beta, arb(T)

def main():
    ctx.prec = 600
    L_num, L_den, T = 4, 5, 200
    L = arb(L_num) / arb(L_den)
    A = B.A_L(L_num, L_den)
    beta = B.H(T) - A
    print("=== LOCAL-OCCUPATION PILOT, L = 4/5 ===")
    print("beta = H(T) - A_L = %s  at T = %d" % (beta.str(12), T))
    print()
    print("-- 1. At the gate cutoff T = 200 the unfavourable set is EMPTY.")
    print("   Psi_L(t) >= H(t) - A_L >= H(T) - A_L = beta for every t >= T,")
    print("   by monotonicity of H.  So E = {} and b*q = 0: there is nothing")
    print("   for an occupation bound to buy at the baseline cutoff.")
    print()
    print("-- 2. The allowance only matters for a SMALLER cutoff T' < T.")
    print("   Then E is contained in the finite band [T', T], of measure <= T - T'.")
    print("   The finite-block margin that must absorb b*q is of order lam_A ~ 1e-17")
    print("   (the window minimum itself).  Required q:")
    b = B.M_sup(L_num, L_den, beta, T)          # b <= sup|Psi - beta| <= M
    print("   b <= M = %s" % b.str(10))
    for margin in ("1e-17", "9e-18", "1e-20"):
        need = arb(margin) / b
        print("     margin %-7s  ->  need q <= %s" % (margin, need.str(8)))
    print()
    print("-- 3. What Lemmas A and B can deliver on a band of measure meas(E):")
    print("   Lemma A: q = L*meas(E)/pi.  To reach q = 1e-18 needs")
    for qt in ("1e-18",):
        m = arb(qt) * arb.pi() / L
        print("     meas(E) <= %s   -- i.e. the band must be shorter than %s" % (m.str(8), m.str(8)))
    print("   Lemma B: q = eps*(1+2*L*h).  Even with h -> 0 this needs eps <= 1e-18,")
    print("     i.e. the unfavourable set must occupy a 1e-18 fraction of EVERY cell.")
    print()
    print("-- 4. VERDICT")
    print("   The elementary local-energy lemmas are ~18 orders of magnitude too weak")
    print("   to be consumed at this margin.  They are not useless in general -- they")
    print("   are the correct tools and they do control int_E|F|^2 -- but the quantity")
    print("   they must beat here is the window minimum itself, ~1e-17, while b ~ 8.8.")
    print("   A Logvinenko-Sereda / Kovrijkine thick-set bound has the same shape")
    print("   (q bounded away from 0 by a power of the thickness) and cannot reach")
    print("   1e-18 either on a set of positive density.")
    print()
    print("   So: NO occupation bound was consumed in a certificate, and none can be")
    print("   at L = 4/5 with this margin.  The pilot's honest outcome is that the")
    print("   exterior is already handled EXACTLY by the gate (E empty), and the")
    print("   allowance route would only trade a smaller T for a penalty that the")
    print("   margin cannot pay.")
    out = dict(L="4/5", T=T, beta=beta.str(20), b=b.str(20),
               E_at_gate="empty", q_needed_for_1e17_margin=(arb("1e-17") / b).str(20),
               verdict="not consumable at this margin")
    d = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "certs")
    json.dump(out, open(os.path.join(d, "pilot.json"), "w"), indent=2)
    return 0

if __name__ == "__main__":
    sys.exit(main())
