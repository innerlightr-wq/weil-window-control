#!/usr/bin/env python3
"""Cold verifier.  Exits NONZERO on a genuine validation failure.

An inconclusive research outcome is NOT a failure: whole-window status is
reported as UNRESOLVED and that is a PASS.  What fails this verifier is a broken
computation, a regression against a frozen input, an undetected mutation, or a
report that confuses finite-subspace positivity with whole-window positivity.

Run:  PYTHONPATH=src:../../src ../../.venv/bin/python verify.py
"""
import json, os, sys
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "src"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))

import mpmath as mp
from flint import arb
import zeta_window as zw
from wcore import set_prec, comb_mass, prime_powers_below, ball_str
from certify import Win

FAIL = []
def chk(name, ok, extra=""):
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name, ("  -- " + extra) if extra else ""))
    if not ok:
        FAIL.append(name)
    return ok

def main():
    set_prec(400)
    print("=== COLD VERIFIER: positive-weil-certificate-2026-11-09 ===")

    # ---------------------------------------------------------------- V1 exact algebra
    print("-- V1 exact algebra")
    from wcore import cheb_to_poly, poly_mul, poly_eq
    P = cheb_to_poly([3, 4, 1])
    tgt = [2 * x for x in poly_mul([Fraction(1), Fraction(1)], [Fraction(1), Fraction(1)])]
    chk("V1a 3+4cos+cos2 == 2(1+cos)^2 exactly", poly_eq(P, tgt))
    bad = cheb_to_poly([3, 4, 2])          # MUTATION: wrong coefficient
    chk("V1b mutation detected: 3+4cos+2cos2 is NOT 2(1+cos)^2", not poly_eq(bad, tgt))

    # ---------------------------------------------------------- V2 prime activation
    print("-- V2 prime activation at L_5 = (log 5)/2")
    L5 = arb(5).log() / 2
    chk("V2a L_5 strictly between 4/5 and 81/100", bool(L5 > arb(4) / 5) and bool(L5 < arb(81) / 100),
        L5.str(14))
    p8 = [n for n, _, _ in prime_powers_below(8, 5)]
    p81 = [n for n, _, _ in prime_powers_below(162, 100)]
    chk("V2b L=4/5 activates exactly {2,3,4}", p8 == [2, 3, 4], str(p8))
    chk("V2c L=81/100 activates exactly {2,3,4,5}", p81 == [2, 3, 4, 5], str(p81))
    A8, _ = comb_mass(8, 5); A81, _ = comb_mass(162, 100)
    exact = 2 * arb(5).log() / arb(5).sqrt()
    chk("V2d A_L jump == 2 log 5 / sqrt 5 exactly",
        bool((A81 - A8 - exact).abs_upper() < arb(10) ** -30), ball_str(A81 - A8, 12))
    chk("V2e A_L(4/5) reproduces the external 2.9420",
        bool((A8 - arb("2.9420")).abs_upper() < arb("0.0001")), ball_str(A8, 10))

    # --------------------------------------------- V3 two backends agree on the form
    print("-- V3 two implementations of the Weil form agree (Arb balls vs mpmath)")
    W = Win(4, 5); N = 6
    A = W.matrix(N)
    mp.mp.dps = 60
    B = zw.build_matrix(mp.mpf(4) / 5, N, 60, sector=zw.EVEN)
    worst = max(abs(float(A[i][j].mid()) - float(B[i, j])) for i in range(N) for j in range(N))
    chk("V3a Arb and mpmath agree to < 1e-30 on every entry (N=6, L=4/5)", worst < 1e-30,
        "max diff %.3g" % worst)
    chk("V3b every Arb entry is a finite ball", all(A[i][j].is_finite() for i in range(N) for j in range(N)))
    # MUTATION: drop the prime term and require the form to change
    import copy
    pp = prime_powers_below(8, 5)
    chk("V3c mutation detected: removing a prime power changes A_L",
        bool((comb_mass(8, 5)[0] - comb_mass(6, 5)[0]).abs_lower() > arb("0.3")),
        "dropping n=4 shifts A_L by %s" % ball_str(comb_mass(8, 5)[0] - comb_mass(6, 5)[0], 8))

    # ------------------------------------------------------------ V4 both parity sectors
    print("-- V4 both parity sectors assembled")
    mp.mp.dps = 50
    for sec, nm in ((zw.EVEN, "even"), (zw.EVEN_D, "even_d"), (zw.ODD, "odd")):
        M = zw.build_matrix(mp.mpf(4) / 5, 5, 50, sector=sec)
        E = mp.eigsy(M, eigvals_only=True)
        chk("V4 %s sector assembles with finite spectrum" % nm,
            all(mp.isfinite(e) for e in E), "lambda_min = %s" % mp.nstr(min(E), 8))
    chk("V4 pole sign is -1 in the odd sector and +1 in the even sectors",
        zw.pole_sign(zw.ODD) == -1 and zw.pole_sign(zw.EVEN) == 1 and zw.pole_sign(zw.EVEN_D) == 1)

    # ---------------------------------------------- V5 frozen certified upper bounds
    print("-- V5 frozen certified upper bounds recomputed from the frozen vectors")
    fp = os.path.join(HERE, "certs", "frozen_vectors.json")
    if not os.path.exists(fp):
        chk("V5 frozen_vectors.json present", False, "missing -- run src/upper_bound.py")
    else:
        rows = json.load(open(fp))
        for r in rows:
            a, b = (int(v) for v in r["L"].split("/"))
            x = [Fraction(s) for s in r["x"]]
            set_prec(r["prec"])
            Wk = Win(a, b)
            q = Wk.quad_form(Wk.matrix(len(x)), x)
            ref = arb(r["certified_upper_bound"])
            ok = bool(q.upper() <= ref * arb("1.000001")) and bool(q.is_finite())
            chk("V5 L=%s N=%d: frozen vector reproduces its certified bound"
                % (r["L"], r["N"]), ok, "recomputed %s vs stored %s"
                % (arb(q.upper()).str(8), ref.str(8)))
            chk("V5 L=%s N=%d: the bound is an UPPER bound and is labelled so"
                % (r["L"], r["N"]), "certified_upper_bound" in r and "lower" not in json.dumps(r).lower())

    # ------------------------------------- V6 status cannot be confused by reporting
    print("-- V6 reporting cannot confuse finite with whole-window")
    rep = open(os.path.join(HERE, "REPORT.md")).read() if os.path.exists(os.path.join(HERE, "REPORT.md")) else ""
    hyp = open(os.path.join(HERE, "HYPOTHESES.md")).read()
    chk("V6a HYPOTHESES marks the reduction's matrix step NOT executed",
        "NOT EXECUTED HERE" in hyp or "NOT REPRODUCED" in hyp)
    chk("V6b whole-window status is declared UNRESOLVED somewhere",
        "UNRESOLVED" in hyp and ("UNRESOLVED" in rep or rep == ""))
    bad_claim = any(s in rep for s in ["we prove Q_zeta(f) >= m_L", "whole-window positivity is established",
                                       "RH follows"])
    chk("V6c no whole-window lower-bound claim appears in REPORT.md", not bad_claim)

    print("\n%s" % ("COLD VERIFY PASSED" if not FAIL else
                    "COLD VERIFY FAILED: " + "; ".join(FAIL)))
    return 1 if FAIL else 0

if __name__ == "__main__":
    sys.exit(main())
