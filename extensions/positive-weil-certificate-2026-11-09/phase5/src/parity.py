"""Replay the T=400 certificate through the repaired precision path, and extract
the certified parity separation.

WHY THE CACHED MATRIX ENCLOSURES REMAIN VALID.  The defect was that routines
READ the ambient precision.  `block.assemble` does not: it sets its own
`prec_bessel` on entry and `prec_mat` around each flush, so every arithmetic
operation that produced a matrix entry ran at a precision chosen by the
assembly, not inherited from a caller.  (It LEAKED that precision outward, which
is what bit the constants computed afterwards, but it never read one in.)  The
entries were then serialised as Arb balls with explicit radii and are re-read by
`arb(str)`, which is exact and precision-independent.  So the matrices need no
reassembly; only the DEPENDENT constants -- beta, M, S_N, g, kappa, d -- were
affected, and those are recomputed here at an owned precision.
"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "src"))
sys.path.insert(0, HERE)
from flint import arb, ctx
import block as BK
import constants as K
from prec import owned

CERT4 = os.path.join(HERE, "..", "..", "phase4", "certs")
OUT = os.path.join(HERE, "..", "certs")
PREC = 2000
LDL = 1400

def load(p):
    return [[arb(x) for x in row] for row in json.load(open(p))]

def main():
    print("=== T=400 replay through the repaired path, and parity separation ===")
    c = K.constants(4, 5, 400, 440, prec=PREC)
    beta, M, kappa = c["beta"], c["M"], c["kappa"]
    print("  beta   = %s" % beta.str(20))
    print("  M      = %s" % M.str(16))
    print("  g_N    = %s" % c["g"].str(12))
    print("  kappa  = %s   kappa^2 = %s" % (kappa.str(12), (kappa * kappa).str(12)))
    print("  d_even = %s" % c["d_even"].str(20))
    print("  d_odd  = %s" % c["d_odd"].str(20))

    res = {}
    with owned(LDL):
        for sec, dkey, shifts in (("even", "d_even", ("1.42e-17",)),
                                  ("odd", "d_odd", ("2.4e-17",))):
            A = load(os.path.join(CERT4, "block400_%s.json" % sec))
            a = arb(shifts[0])
            ok, resid, err = BK.verify_pd_preconditioned(A, a, LDL)
            print("\n  %s: A - %s I  PD = %s   ||I-G||_inf = %s"
                  % (sec, a.str(8), ok, "n/a" if resid is None else resid.str(8)))
            if not ok:
                raise RuntimeError("replay failed in the %s sector: %s" % (sec, err))
            res[sec] = (a, c[dkey], resid)

    print("\n-- the even sector (binding) retains mu = 1.41e-17")
    ae, de, _ = res["even"]
    re_ = BK.certify_final(ae, de, kappa, arb("1.41e-17"))
    print("     a=%s d=%s mu=1.41e-17 ok=%s margin=%s"
          % (ae.str(8), de.str(14), re_["ok"], re_["margin"]))

    print("-- the ODD sector, re-evaluated against mu_odd = 2.3e-17")
    ao, do, _ = res["odd"]
    ro = BK.certify_final(ao, do, kappa, arb("2.3e-17"))
    print("     a>mu:   %s" % ro["a_gt_mu"])
    print("     d>mu:   %s" % ro["d_gt_mu"])
    print("     (a-mu)(d-mu) > kappa^2: %s" % ro["disc_gt_kappa2"])
    print("     margin = %s" % ro["margin"])
    if not ro["ok"]:
        raise RuntimeError("mu_odd = 2.3e-17 does NOT certify")

    print("\n-- the separation")
    cross = json.load(open(os.path.join(CERT4, "cross.json")))
    qU = arb(cross["Q_f400"])
    U = arb(qU.upper())                       # exact upper endpoint of the enclosure
    mu_odd = arb("2.3e-17")
    sep = mu_odd - U
    print("     U = upper endpoint of q(f_400) = %s" % U.str(20))
    print("     mu_odd                          = %s" % mu_odd.str(20))
    print("     separation mu_odd - U           = %s" % sep.str(16))
    ok_sep = bool(sep > 0)
    print("     certified U < mu_odd : %s" % ok_sep)
    print()
    print("     f_400 is EVEN, so   inf_even Q/||f||^2 <= U < mu_odd <= inf_odd Q/||f||^2 .")
    print("     With the orthogonal parity decomposition this identifies the EVEN")
    print("     infimum as the whole-window infimum.  It does NOT prove simplicity,")
    print("     uniqueness, or that the infimum is attained.")
    print("     Parity separation itself is not new (it is in Zhu); this is an")
    print("     implication of THIS pipeline's own certified data.")

    out = dict(prec=PREC, beta=beta.str(25), M=M.str(25), kappa=kappa.str(25),
               d_even=c["d_even"].str(25), d_odd=c["d_odd"].str(25),
               even=dict(a=ae.str(25), mu="1.41e-17", margin=re_["margin"], ok=re_["ok"]),
               odd=dict(a=ao.str(25), mu_odd="2.3e-17", margin=ro["margin"], ok=ro["ok"]),
               U_upper_endpoint=U.str(25), separation=sep.str(25),
               separation_certified=ok_sep,
               whole_window_lower="1.41e-17 (even sector binds)")
    os.makedirs(OUT, exist_ok=True)
    json.dump(out, open(os.path.join(OUT, "parity.json"), "w"), indent=2)
    print("\n-> certs/parity.json")
    return 0 if (re_["ok"] and ro["ok"] and ok_sep) else 1

if __name__ == "__main__":
    sys.exit(main())
