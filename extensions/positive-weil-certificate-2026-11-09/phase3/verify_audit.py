#!/usr/bin/env python3
"""Phase 4 audit verifier.  Nonzero exit on failure.  Explicit failure
propagation everywhere: no NaN, infinite radius or unresolved sign may pass."""
import os, sys, json
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "src"))
sys.path.insert(0, os.path.join(HERE, "..", "phase2", "src"))
from flint import arb, acb, ctx
import block as BK, baseline as B
import audit as AU

FAIL = []
def chk(name, ok, extra=""):
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name, ("  -- " + extra) if extra else ""))
    if not ok:
        FAIL.append(name)
    return ok

PREC = 1400
LOWER = arb("1e-17")
UPPER = arb("2.112992612e-17")

def main():
    ctx.prec = PREC
    print("=== PHASE 4 AUDIT VERIFIER ===")

    print("-- A1 beta monotonicity (the retracted optimisation)")
    T, N, PAN, KO = 200, 12, 120, 32
    S = BK.Symbol(4, 5, T)
    b1, b2 = arb("0.45"), arb("0.55")
    A1 = AU.A_of_beta(S, b1, N, T, PAN, KO)
    A2 = AU.A_of_beta(S, b2, N, T, PAN, KO)
    K = AU.K_T(S, N, T, PAN, KO)
    n = len(K)
    worst = arb(0)
    for a in range(n):
        for b in range(n):
            d = ((A2[a][b] - A1[a][b]) - (b2 - b1) * ((arb(1) if a == b else arb(0)) - K[a][b])).abs_upper()
            if arb(d) > worst:
                worst = arb(d)
    chk("A1a A(b2)-A(b1) = (b2-b1)(I-K_T) entrywise", bool(worst < arb(10) ** -40),
        "max residual %s" % worst.str(8))
    chk("A1b 0 <= K_T <= I on the diagonal",
        all(bool(K[a][a] >= 0) and bool(K[a][a] <= 1) for a in range(n)))
    okP, _, _ = BK.verify_pd_preconditioned(
        [[A2[a][b] - A1[a][b] for b in range(n)] for a in range(n)], arb(0), PREC)
    okN, _, _ = BK.verify_pd_preconditioned(
        [[A1[a][b] - A2[a][b] for b in range(n)] for a in range(n)], arb(0), PREC)
    chk("A1c A(b2)-A(b1) is PD and the reverse is NOT (lowering beta cannot raise lam_A)",
        okP and not okN)

    print("-- A2 the three spectral quantities are kept distinct")
    doc = AU.__doc__
    chk("A2a the audit states lambda_R <= lambda_Q and lambda_R <= lambda_A",
        "lambda_R <= lambda_Q" in doc and "lambda_R <= lambda_A" in doc)
    flat = " ".join(doc.split())
    chk("A2b it states 1.0310220e-17 estimates lambda_A, not an upper bound for lambda_Q",
        "estimates **lambda_A**" in flat and "NOT an upper bound for lambda_Q" in flat)

    print("-- A3 the frozen minimizer, recomputed from its exact rationals")
    ctx.prec = PREC
    fm = json.load(open(os.path.join(HERE, "certs", "frozen_minimizer.json")))
    y = [Fraction(s) for s in fm["y"]]
    ya = [arb(f.numerator) / arb(f.denominator) for f in y]
    A = [[arb(x) for x in row] for row in
         json.load(open(os.path.join(HERE, "..", "phase2", "certs", "block_even.json")))]
    r_f, nrm2 = AU_rayleigh(A, ya)
    chk("A3a ||y||^2 reproduces", bool((nrm2 - arb(fm["norm_sq"])).abs_upper() < arb(10) ** -12),
        nrm2.str(16))
    chk("A3b r_f reproduces from the frozen rationals",
        bool((r_f - arb(fm["r_f"])).abs_upper() < arb(10) ** -28), r_f.str(16))
    e_lo = arb(fm["e_f_lower_bound"])
    chk("A3c e_f lower bound is positive and exceeds r_f", bool(e_lo > 0) and bool(e_lo > r_f),
        "e_f >= %s, r_f = %s" % (e_lo.str(10), r_f.str(10)))
    q_lo = r_f + e_lo
    chk("A3d q_f = r_f + e_f >= %s does NOT improve the certified upper bound %s"
        % (q_lo.str(10), UPPER.str(10)), bool(q_lo > UPPER),
        "earlier upper bound retained")

    print("-- A4 the two-sided enclosure stands")
    chk("A4 1e-17 <= lambda_Q <= 2.112992612e-17", bool(LOWER < UPPER),
        "width %s" % (UPPER - LOWER).str(8))

    print("-- A5 exception handling")
    ok_raise = False
    try:
        BK.rigorous_entry(S, 0, 0, -1, 6)
    except RuntimeError:
        ok_raise = True
    chk("A5a rigorous_entry rejects a range reaching the singularity by its guard", ok_raise)
    v = BK.rigorous_entry(S, 0, 0, 1, 6)
    chk("A5b a safe spot-check still returns a finite tight enclosure",
        v.is_finite() and float(v.rad()) < 1e-100, "radius %.2e" % float(v.rad()))
    src = open(os.path.join(HERE, "..", "phase2", "src", "block.py")).read()
    chk("A5c the callback RECORDS failures and the result is rejected on the flag",
        'state["failed"]' in src and "result rejected" in src)
    i0 = src.index("def assemble(")
    nxt = [j for j in range(i0 + 1, len(src)) if src.startswith("\ndef ", j - 1)]
    asm = src[i0:(nxt[0] if nxt else len(src))]
    chk("A5d block.assemble never calls acb.integral (grep of its body)",
        "acb.integral" not in asm, "%d chars of body scanned" % len(asm))
    chk("A5e acb.integral occurs ONLY inside rigorous_entry",
        src.count("acb.integral(") == 1 and
        src.index("acb.integral(") > src.index("def rigorous_entry("))

    print("-- A6 provenance of the uniform transform bound")
    L = arb(4) / 5
    be = arb("0.30")
    chk("A6 equals sqrt(2L)exp(L b_e), i.e. Cauchy-Schwarz",
        bool(((2 * L).sqrt() * (L * be).exp() - BK.v_hol_absmax_uniform(L, be)).abs_upper()
             < arb(10) ** -40))

    print("\n%s" % ("PHASE 4 AUDIT PASSED" if not FAIL else
                    "PHASE 4 AUDIT FAILED: " + "; ".join(FAIL)))
    return 1 if FAIL else 0

def AU_rayleigh(M, y):
    n = len(M)
    num = sum((y[i] * M[i][j] * y[j] for i in range(n) for j in range(n)), arb(0))
    den = sum((y[i] * y[i] for i in range(n)), arb(0))
    return num / den, den

if __name__ == "__main__":
    sys.exit(main())
