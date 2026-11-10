"""Items 4 and 5: the missing scalar comparison R_200(f_400), and a FIXED-FUNCTION
T=800 pilot.  No T=800 matrix is assembled and none is authorised."""
import os, sys, json, time
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "src"))
sys.path.insert(0, HERE)
from flint import arb
import baseline as B
import scalar as SC, constants as K
from prec import owned

C3 = os.path.join(HERE, "..", "..", "phase3", "certs")
C4 = os.path.join(HERE, "..", "..", "phase4", "certs")
OUT = os.path.join(HERE, "..", "certs")
PREC = 1600

def main():
    t0 = time.time()
    res = {}
    with owned(PREC):
        A_L = B.A_L(4, 5)
        b200, b400, b800 = B.H(200) - A_L, B.H(400) - A_L, B.H(800) - A_L
    print("beta_200 = %s" % b200.str(18))
    print("beta_400 = %s" % b400.str(18))
    print("beta_800 = %s" % b800.str(18))
    print("beta_800 > beta_400 : %s   (increase %s)"
          % (bool(b800 > b400), (b800 - b400).str(14)))
    res["beta_200"], res["beta_400"], res["beta_800"] = (b200.str(25), b400.str(25), b800.str(25))
    res["beta_800_gt_400"] = bool(b800 > b400)

    f200 = [Fraction(s) for s in json.load(open(os.path.join(C3, "frozen_minimizer.json")))["y"]]
    d200 = [2 * i for i in range(len(f200))]
    f400 = [Fraction(s) for s in json.load(open(os.path.join(C4, "cross.json")))["f400_y"]]
    d400 = [2 * i for i in range(len(f400))]
    q200 = arb(json.load(open(os.path.join(C4, "q_f200.json")))["q"])
    q400 = arb(json.load(open(os.path.join(C4, "cross.json")))["Q_f400"])
    print("\nQ(f_200)/||.||^2 = %s   Q(f_400)/||.||^2 = %s" % (q200.str(20), q400.str(20)))

    jobs = [("R_200(f_400)", f400, d400, 200, b200, 200),
            ("R_800(f_400)", f400, d400, 800, b800, 800),
            ("R_800(f_200)", f200, d200, 800, b800, 800)]
    for name, y, degs, T, beta, panels in jobs:
        t1 = time.time()
        print("\n-- %s  (T=%d, %d panels, order 64, %d modes) ..." % (name, T, panels, len(y)),
              flush=True)
        R, meta = SC.R_scalar(4, 5, y, degs, T, beta, panels, 64, PREC)
        print("   %s = %s   [CERTIFIED]  (quad err <= %s, %.0fs)"
              % (name, R.str(20), meta["quad_err"], time.time() - t1))
        res[name] = R.str(25)
        res[name + "_quad_err"] = meta["quad_err"]

    with owned(PREC):
        r400_f400 = arb(json.load(open(os.path.join(C4, "cross.json")))["R400_f400"])
        r400_f200 = arb(json.load(open(os.path.join(C4, "cross.json")))["R400_f200"])
        e400_f400 = q400 - r400_f400
        e800_f400 = q400 - arb(res["R_800(f_400)"])
        e400_f200 = q200 - r400_f200
        e800_f200 = q200 - arb(res["R_800(f_200)"])
        print("\n-- fixed-function discarded contributions")
        for nm, e4, e8 in (("f_400", e400_f400, e800_f400), ("f_200", e400_f200, e800_f200)):
            ratio = e8 / e4
            print("   %s : e_400 = %-24s e_800 = %-24s  e_800/e_400 = %s"
                  % (nm, e4.str(14), e8.str(14), ratio.str(12)))
            res["e400_" + nm] = e4.str(25)
            res["e800_" + nm] = e8.str(25)
            res["ratio_" + nm] = ratio.str(25)
        print("\n   One ratio per function is NOT a decay law in T, and comparing")
        print("   different optimised vectors would not establish a rate either.")
        print("   R_800 on a fixed vector is a variational UPPER bound for lambda_R800;")
        print("   it gives no whole-window lower bound and does not predict f_800.")
    res["elapsed_s"] = round(time.time() - t0, 1)
    os.makedirs(OUT, exist_ok=True)
    json.dump(res, open(os.path.join(OUT, "pilot800.json"), "w"), indent=2)
    print("\n-> certs/pilot800.json   (%.0fs)" % (time.time() - t0))
    return 0

if __name__ == "__main__":
    sys.exit(main())
