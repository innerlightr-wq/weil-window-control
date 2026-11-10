"""Regenerate phase4/certs/certificate_T400.json from first principles.

The file was originally written ad hoc, which left it without a committed
producer even though `phase4/verify_t400.py` never trusts it: that verifier
recomputes beta, M, S_N, g, kappa and d from `baseline`, and reads only the
proposed shift `a` and the proposed constant `mu` from the file before re-running
the preconditioned PD check and the three strict Schur signs itself.

This script closes the provenance gap.  It rebuilds every derived field through
the OWNED-precision path of `phase5/src/constants.py`, re-runs
`block.certify_final` for both sectors, and refuses to write anything if a
sector fails.  The two proposed shifts are the inputs, exactly as they were:

    a_even = 1.42e-17   (certified PD, ||I-G||_inf = 1.4973903e-14)
    a_odd  = 2.4e-17    (certified PD, ||I-G||_inf = 4.3535838e-22)

Run:
    PYTHONPATH=phase2/src:phase4/src:phase5/src python phase4/src/emit_certificate_T400.py
"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "src"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase5", "src"))
from flint import arb
import block as BK
import constants as K
from prec import owned

CERTS = os.path.join(HERE, "..", "certs")
L_NUM, L_DEN, T, N = 4, 5, 400, 440
PREC = 2000
SHIFTS = {"even": "1.42e-17", "odd": "2.4e-17"}
MUS    = {"even": "1.41e-17", "odd": "1.47e-17"}

def main():
    c = K.constants(L_NUM, L_DEN, T, N, prec=PREC)
    # OWN THE INPUTS, not only the scope: arb("1.42e-17") built at the ambient
    # precision carries a radius of ~3e-33 and serialises a needlessly wide
    # record.  See REPORT.md section 11.1.
    with owned(PREC):
        shifts = {k: arb(v) for k, v in SHIFTS.items()}
        mus    = {k: arb(v) for k, v in MUS.items()}
        return _emit(c, shifts, mus)

def _emit(c, shifts, mus):
    out = {
        "beta":   c["beta"].str(25),
        "M":      c["M"].str(25),
        "g":      c["g"].str(25),
        "kappa":  c["kappa"].str(25),
        "d_even": c["d_even"].str(25),
        "d_odd":  c["d_odd"].str(25),
        "T": T, "N": N, "sectors": {},
    }
    ok_all = True
    for sec in ("even", "odd"):
        a  = shifts[sec]
        d  = c["d_%s" % sec]
        mu = mus[sec]
        r  = BK.certify_final(a, d, c["kappa"], mu)
        ok_all = ok_all and bool(r["ok"])
        print("  %-4s a=%s d=%s mu=%s  ok=%s  margin=%s"
              % (sec, a.str(10), d.str(16), mu.str(6), r["ok"], r["margin"]))
        out["sectors"][sec] = {
            "a": a.str(25), "d": d.str(25), "mu": MUS[sec],
            "margin": r["margin"],
            "kappa2_over_d": ((c["kappa"] * c["kappa"]) / d).str(25),
        }
    if not ok_all:
        print("REFUSING to write: a sector did not certify")
        return 1
    lo = mus["even"] if bool(mus["even"] < mus["odd"]) else mus["odd"]
    out["whole_window_mu"] = lo.str(25)
    path = os.path.join(CERTS, "certificate_T400.json")
    json.dump(out, open(path, "w"), indent=2)
    print("-> %s" % os.path.relpath(path, os.path.join(HERE, "..", "..")))
    return 0

if __name__ == "__main__":
    sys.exit(main())
