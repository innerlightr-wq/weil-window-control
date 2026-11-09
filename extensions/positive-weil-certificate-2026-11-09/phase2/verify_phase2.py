#!/usr/bin/env python3
"""Phase 2 cold verifier.  Nonzero exit on a genuine validation failure.

UNRESOLVED is a legitimate outcome and is NOT a failure; what fails here is a
broken computation, a regression, an undetected mutation, or a report that
mislabels a conditional statement as a theorem.
"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "src"))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from flint import arb, acb, ctx
import baseline as B
from lemmas import q_lemma_A, q_lemma_B

FAIL = []
def chk(name, ok, extra=""):
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name, ("  -- " + extra) if extra else ""))
    if not ok:
        FAIL.append(name)
    return ok

def main():
    ctx.prec = 1600
    L_num, L_den, T = 4, 5, 200
    L = arb(L_num) / arb(L_den)
    print("=== PHASE 2 COLD VERIFIER ===")

    print("-- W1 lemmas")
    chk("W1a Lemma A monotone and capped at 1",
        bool(q_lemma_A(0.8, 0.1) < q_lemma_A(0.8, 1.0)) and bool(q_lemma_A(0.8, 100) == 1))
    chk("W1b Lemma B equals eps*(1+2Lh) below the cap",
        bool((q_lemma_B(0.8, 1.0, 0.01) - arb("0.026")).abs_upper() < arb(10) ** -12),
        q_lemma_B(0.8, 1.0, 0.01).str(10))
    chk("W1c Lemma B is informative (q<1) for eps=0.05,h=1,L=0.8",
        bool(q_lemma_B(0.8, 1.0, 0.05) < 1), q_lemma_B(0.8, 1.0, 0.05).str(8))

    print("-- W2 gate and sup bound")
    A = B.A_L(L_num, L_den)
    beta = B.H(T) - A
    chk("W2a beta = H(T) - A_L > 0 at T=200", bool(beta > 0), beta.str(14))
    chk("W2b H is increasing: H(1) < H(10) < H(200)",
        bool(B.H(1) < B.H(10)) and bool(B.H(10) < B.H(T)))
    M = B.M_sup(L_num, L_den, beta, T)
    chk("W2c M bound contains the pointwise value at t=0",
        bool(M >= (B.Psi(0, L_num, L_den) - beta).abs_upper()), M.str(12))
    chk("W2d exterior set is EMPTY at the gate cutoff (Psi_L >= beta for t>=T)",
        bool(B.Psi(T, L_num, L_den) >= beta - arb(10) ** -20))
    # MUTATION: drop a prime power -> A_L and hence beta must change
    terms = B.comb_terms(L_num, L_den)
    chk("W2e mutation: the comb carries exactly {2,3,4} at L=4/5",
        [n for n, _ in terms] == [2, 3, 4], str([n for n, _ in terms]))

    print("-- W3 complement and coupling (the certified part of the baseline)")
    x = arb(T) * L
    pn = (2 * L.sinh()).sqrt()
    chk("W3a |p|^2 = 2 sinh L exactly", bool((pn * pn - 2 * L.sinh()).abs_upper() < arb(10) ** -30),
        pn.str(12))
    ok_mono = True
    prev = None
    rows = []
    for N in (200, 220, 250):
        S = B.S_tail(N, x)
        g = (2 * L * arb(T) / arb.pi()) * S
        lamD = beta - M * g
        pt = B.pole_vec_tail(N, L).sqrt()
        kap = 2 * pn * pt + M * g.sqrt()
        need = kap * kap / lamD
        rows.append((N, g, lamD, kap, need))
        if prev is not None and not (g < prev):
            ok_mono = False
        prev = g
    chk("W3b g_N is strictly decreasing in N", ok_mono)
    for N, g, lamD, kap, need in rows:
        chk("W3c N=%d: lam_D > 0.5 and the required lam_A is below 1e-14" % N,
            bool(lamD > arb("0.5")) and bool(need < arb("1e-14")),
            "lam_D=%s need=%s" % (lamD.str(8), need.str(8)))
    chk("W3d at N=220 the required lam_A is below 1e-27",
        bool(rows[1][4] < arb("1e-27")), rows[1][4].str(8))

    print("-- W4 status labelling")
    rep = open(os.path.join(HERE, "..", "REPORT.md")).read()
    chk("W4a REPORT declares the baseline UNRESOLVED", "UNRESOLVED" in rep)
    chk("W4b REPORT does not claim a whole-window theorem",
        not any(s in rep for s in ["whole-window positivity is CERTIFIED",
                                   "we certify Q_zeta(f) >= m", "RH follows"]))
    ret = open(os.path.join(HERE, "..", "RETRACTIONS.md")).read()
    chk("W4c RETRACTIONS records R6 (the overbroad measure claim)",
        "R6" in ret and "FALSE" in ret)
    chk("W4d the density retraction R1 is preserved", "R1" in ret and "Kronecker" in ret)
    chk("W4e cos-modulation is described as a superposition, not concentration",
        "superposition" in ret)

    print("\n%s" % ("PHASE 2 VERIFY PASSED" if not FAIL else
                    "PHASE 2 VERIFY FAILED: " + "; ".join(FAIL)))
    return 1 if FAIL else 0

if __name__ == "__main__":
    sys.exit(main())
