#!/usr/bin/env python3
"""Phase 5 cold verifier.  Nonzero exit on failure.  Fail-closed throughout."""
import os, sys, json
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "src"))
sys.path.insert(0, os.path.join(HERE, "..", "phase2", "src"))
from flint import arb, ctx
import block as BK, baseline as B, compact as C

FAIL = []
def chk(name, ok, extra=""):
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name, ("  -- " + extra) if extra else ""))
    if not ok:
        FAIL.append(name)
    return ok

CERTS = os.path.join(HERE, "certs")
PREC = 2000

def main():
    ctx.prec = PREC
    L = arb(4) / 5
    print("=== PHASE 5 COLD VERIFIER (T = 400) ===")

    print("-- T1 cutoff ordering")
    A_L = B.A_L(4, 5)
    b200, b400 = B.H(200) - A_L, B.H(400) - A_L
    chk("T1a beta_400 > beta_200 > 0", bool(b400 > b200) and bool(b200 > 0),
        "%s > %s" % (b400.str(14), b200.str(14)))
    chk("T1b Psi_L(400) >= beta_400 (gate at the new cutoff)",
        bool(B.Psi(arb(400), 4, 5) - b400 >= -arb(10) ** -40))

    print("-- T2 compact evaluator validated against the repository closed form")
    S = BK.Symbol(4, 5, 200)
    f0 = C.Poly(L, [0], [arb(1)])
    q0, _ = C.Q_compact(f0, +1, S.terms)
    ref = arb("0.080274054959269611214")
    chk("T2a Q_compact(phi_0) matches the repo's exact Q[0,0] to 1e-20",
        bool((q0 - ref).abs_upper() < arb(10) ** -20), q0.str(22))
    chk("T2b c = -g'(0+) = f(L)^2 exactly on phi_0",
        bool((f0.at(L) ** 2 - arb(1) / (2 * L)).abs_upper() < arb(10) ** -30),
        (f0.at(L) ** 2).str(14))
    g0 = C.autocorr(f0, 0)
    chk("T2c g(0) = ||f||^2", bool((g0 - f0.norm2()).abs_upper() < arb(10) ** -30))

    print("-- T3 the T=400 complement and coupling, recomputed")
    T, N = 400, 440
    beta = b400
    M = B.M_sup(4, 5, beta, T)
    x = arb(T) * L
    St = B.S_tail(N, x, top=N + 200)
    g = (2 * L * arb(T) / arb.pi()) * St
    pt2 = B.pole_vec_tail(N, L)
    pn = (2 * L.sinh()).sqrt()
    kappa = 2 * pn * pt2.sqrt() + M * g.sqrt()
    d_e, d_o = beta - M * g, beta - M * g - 2 * pt2
    chk("T3a S_N(TL) is tight and tiny", float(St.rad()) < 1e-60 and bool(St < arb("1e-60")),
        St.str(8))
    chk("T3b d > 0 in both sectors", bool(d_e > 0) and bool(d_o > 0), d_e.str(14))
    chk("T3c kappa^2/d is below 1e-50", bool((kappa * kappa) / d_e < arb("1e-50")),
        ((kappa * kappa) / d_e).str(8))

    print("-- T4 the T=400 blocks are PD at the certified shifts")
    cert = json.load(open(os.path.join(CERTS, "certificate_T400.json")))
    for sec, d in (("even", d_e), ("odd", d_o)):
        a = arb(cert["sectors"][sec]["a"])
        Mt = [[arb(v) for v in row] for row in
              json.load(open(os.path.join(CERTS, "block400_%s.json" % sec)))]
        ok, res, err = BK.verify_pd_preconditioned(Mt, a, 1400)
        chk("T4 %s: A - %s I is PD" % (sec, a.str(8)), ok,
            "||I-G||_inf = %s" % ("n/a" if res is None else res.str(8)))
        mu = arb(cert["sectors"][sec]["mu"])
        r = BK.certify_final(a, d, kappa, mu)
        chk("T4 %s: a>mu, d>mu, (a-mu)(d-mu)>kappa^2 at mu=%s" % (sec, mu.str(8)), r["ok"],
            "margin %s" % r["margin"])

    print("-- T5 the cross table and the new bracket")
    cr = json.load(open(os.path.join(CERTS, "cross.json")))
    r200 = arb(cr["R200_f200"]); r400 = arb(cr["R400_f200"])
    q200 = arb(cr["Q_f200"]); q400 = arb(cr["Q_f400"])
    r400f400 = arb(cr["R400_f400"])
    chk("T5a R_400(f_200) > R_200(f_200) (cutoff identity, same vector)",
        bool(r400 > r200), "%s > %s" % (r400.str(12), r200.str(12)))
    chk("T5b e_400(f_200) <= e_200(f_200)", bool((q200 - r400) <= (q200 - r200)),
        "%s <= %s" % ((q200 - r400).str(12), (q200 - r200).str(12)))
    chk("T5c lambda_A(400,even) <= R_400(f_400)", bool(arb(cert["sectors"]["even"]["a"]) < r400f400),
        "a=%s <= %s" % (cert["sectors"]["even"]["a"][:12], r400f400.str(12)))
    LO = arb(cert["whole_window_mu"])
    chk("T5d the new upper bound q(f_400) improves 2.112992612e-17",
        bool(q400 < arb("2.112992612e-17")), q400.str(20))
    chk("T5e the bracket is nonempty: mu <= q(f_400)", bool(LO < q400),
        "%s <= %s" % (LO.str(10), q400.str(12)))
    chk("T5f the new bracket is narrower than the old",
        bool((q400 - LO) < (arb("2.112992612e-17") - arb("1e-17"))),
        "width %s vs %s" % ((q400 - LO).str(8),
                            (arb("2.112992612e-17") - arb("1e-17")).str(8)))

    print("\n%s" % ("PHASE 5 VERIFY PASSED" if not FAIL else
                    "PHASE 5 VERIFY FAILED: " + "; ".join(FAIL)))
    return 1 if FAIL else 0

if __name__ == "__main__":
    sys.exit(main())
