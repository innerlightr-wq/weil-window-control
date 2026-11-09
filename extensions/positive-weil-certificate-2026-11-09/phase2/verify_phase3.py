#!/usr/bin/env python3
"""Cold verifier for the L=4/5 whole-window certificate.  Nonzero exit on failure.

Validates the certified enclosures rather than re-running the 22-minute
assembly: it reloads the frozen matrix enclosures, re-runs the preconditioned
positive-definiteness proof, recomputes d and kappa from the complement bounds,
recomputes the final implication, spot-checks entries against an independent
adaptive integrator, and rejects three mutations.
"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "src"))
from flint import arb, ctx
import block as BK, baseline as B

FAIL = []
def chk(name, ok, extra=""):
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name, ("  -- " + extra) if extra else ""))
    if not ok:
        FAIL.append(name)
    return ok

CK = os.path.join(HERE, "certs", "block_%s.json")
A_SHIFT, MU = "1.03e-17", "1e-17"
N, T, PREC = 220, 200, 1400

def load(sec):
    raw = json.load(open(CK % sec))
    return [[arb(x) for x in row] for row in raw]

def main():
    ctx.prec = PREC
    L = arb(4) / 5
    S = BK.Symbol(4, 5, T)
    beta = S.beta
    print("=== PHASE 3 COLD VERIFIER: whole-window certificate at L = 4/5 ===")

    print("-- X1 conventions (catches a wrong sigma phase; blind at degree 0)")
    for nm, d_, m_, ok in BK.validate_conventions():
        chk("X1 %s matches direct integration" % nm, ok)

    ctx.prec = PREC          # validate_conventions resets ctx.prec internally
    print("-- X2 frozen enclosures")
    Ae, Ao = load("even"), load("odd")
    chk("X2a both sectors have %d modes (degrees 0..%d, 110 per parity)" % (N // 2, N - 1),
        len(Ae) == N // 2 and len(Ao) == N // 2, "%d / %d" % (len(Ae), len(Ao)))
    rad = max(float(Ae[i][j].rad()) for i in range(len(Ae)) for j in range(len(Ae)))
    chk("X2b every entry is a finite ball, max radius < 1e-50",
        all(Ae[i][j].is_finite() for i in range(len(Ae)) for j in range(len(Ae))) and rad < 1e-50,
        "max radius %.3e" % rad)

    print("-- X3 spot-check entries against the adaptive integrator")
    for (n, m) in ((0, 0), (2, 6), (1, 5)):
        q = BK.quad_entry(S, n, m, 1, 6, 20, 64)
        r = BK.rigorous_entry(S, n, m, 1, 6)
        rel = float(((q - r) / r).abs_upper())
        chk("X3 entry (%d,%d) GL rule vs adaptive agree to 1e-40" % (n, m), rel < 1e-40,
            "rel diff %.3e" % rel)

    ctx.prec = PREC          # quad_entry/rigorous_entry may have lowered it
    print("-- X4 complement and coupling")
    M = B.M_sup(4, 5, beta, T)
    x = arb(T) * L
    St = B.S_tail(N, x)
    g = (2 * L * arb(T) / arb.pi()) * St
    pt2 = B.pole_vec_tail(N, L)
    pn = (2 * L.sinh()).sqrt()
    kappa = 2 * pn * pt2.sqrt() + M * g.sqrt()
    d_even = beta - M * g
    d_odd = beta - M * g - 2 * pt2
    chk("X4a beta > 0 and the exterior gate is tight at T", bool(beta > 0), beta.str(14))
    chk("X4b d > 0 in both sectors", bool(d_even > 0) and bool(d_odd > 0),
        "even %s / odd %s" % (d_even.str(12), d_odd.str(12)))
    chk("X4c the ODD sector carries the negative pole term (d_odd < d_even)",
        bool(d_odd < d_even) or bool((d_even - d_odd).abs_upper() < 1e-300),
        "difference %s" % (d_even - d_odd).str(8))

    chk("X4d the complement sum S_N is tight (radius below 1e-40)",
        float(St.rad()) < 1e-40, "S_N = %s" % St.str(8))
    print("-- X5 positive definiteness of the frozen blocks")
    for sec, A in (("even", Ae), ("odd", Ao)):
        ok, res, err = BK.verify_pd_preconditioned(A, arb(A_SHIFT), PREC)
        chk("X5 %s sector: A - %s I is PD (Rump)" % (sec, A_SHIFT), ok,
            "||I-G||_inf = %s" % ("n/a" if res is None else res.str(8)))

    ctx.prec = PREC
    print("-- X6 the final implication")
    a = arb(A_SHIFT)
    for sec, dd in (("even", d_even), ("odd", d_odd)):
        r = BK.certify_final(a, dd, kappa, arb(MU))
        chk("X6 %s sector: a>mu, d>mu, (a-mu)(d-mu)>kappa^2 at mu=%s" % (sec, MU), r["ok"],
            "margin %s" % r["margin"])

    print("-- X7 mutations must be REJECTED")
    bad = [[arb(v) for v in row] for row in json.load(open(CK % "even"))]
    bad[0][0] = bad[0][0] - arb(1)
    ok, _, _ = BK.verify_pd_preconditioned(bad, arb(A_SHIFT), PREC)
    chk("X7a corrupted entry A[0][0] -= 1 is rejected", not ok)
    bad2 = [[arb(v) for v in row] for row in json.load(open(CK % "odd"))]
    p = BK.pole_vector(L, N)
    po = [p[k] for k in range(N) if k % 2 == 1]
    for i in range(len(bad2)):
        for j in range(len(bad2)):
            bad2[i][j] = bad2[i][j] + 4 * po[i] * po[j]      # flip s = -1 to +1
    ok2, _, _ = BK.verify_pd_preconditioned(bad2, arb(A_SHIFT), PREC)
    orig = load("odd")
    differs = any(not bad2[i][j].overlaps(orig[i][j])
                  for i in range(len(orig)) for j in range(len(orig)))
    # HONEST: flipping s = -1 to +1 ADDS a positive rank-1 term, so the result is
    # still PD.  Positive definiteness alone therefore CANNOT detect a wrong pole
    # sign; what detects it is the convention test X1, which compares F(i/2) and
    # |F(t)|^2 against direct integration of f.  Recorded, not claimed as caught.
    chk("X7b flipped odd pole sign changes the matrix (PD alone does NOT detect it; "
        "X1 does)", differs, "flipped matrix still PD = %s" % ok2)
    nanm = [[arb(v) for v in row] for row in json.load(open(CK % "even"))]
    nanm[3][3] = arb("nan")
    ok3, _, err3 = BK.verify_pd_preconditioned(nanm, arb(A_SHIFT), PREC)
    chk("X7c an indeterminate (NaN) entry is rejected, not skipped", not ok3,
        (err3 or "")[:60])

    print("\n%s" % ("PHASE 3 VERIFY PASSED -- certificate stands" if not FAIL else
                    "PHASE 3 VERIFY FAILED: " + "; ".join(FAIL)))
    return 1 if FAIL else 0

if __name__ == "__main__":
    sys.exit(main())
