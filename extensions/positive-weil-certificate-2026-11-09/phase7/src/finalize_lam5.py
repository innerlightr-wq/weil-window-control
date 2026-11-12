"""Re-certify the lambda_- run from its CHECKPOINTED blocks, using the
   repository's repaired owned-precision path (phase5/src/constants.py), which
   exists precisely to prevent the leak that corrupted work/run_block_lam5.py's
   even-sector d and kappa (REPORT.md section 11.1, the same defect class as R13).

   Two-sided bracket search, so a bad start guess cannot silently certify nothing.
"""
import json, os, time
from flint import arb, ctx
import block as BK
import constants as C

INV = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # phase7/
# NOTE: integrated layout -- certs/ src/ logs/ live under this phase7 directory.
# Rewritten from the session path ~/scratch/rh-connes-continuation-2026-11-12 on
# integration; no $HOME path is read.
L_NUM, L_DEN, T, N = 80471895621705018, 100000000000000000, 400, 440
PREC_LDL, PREC_C = 1400, 2000

K = C.constants(L_NUM, L_DEN, T, N, prec=PREC_C)
print("OWNED-PRECISION CONSTANTS at L~ (prec=%d):" % PREC_C)
for k in ("L", "A_L", "beta", "M", "S", "g", "p_perp2", "kappa", "d_even", "d_odd"):
    print("   %-9s = %s" % (k, K[k].str(22)))
print("   lambda^2  = %s" % (2 * K["L"]).exp().str(22))
print("   kappa^2/d_even = %s" % ((K["kappa"] ** 2) / K["d_even"]).str(10))

def pd(A, a):
    ok, res, err = BK.verify_pd_preconditioned(A, a, PREC_LDL)
    return ok, res

def bracket(A, start):
    """two-sided: find (pass, fail) with pass < fail"""
    a = arb(start)
    ok, res = pd(A, a)
    print("      start a=%-24s PD=%s" % (a.str(10), ok), flush=True)
    if ok:
        lo, hi = a, a * 2
        while True:
            ok2, _ = pd(A, hi)
            print("      up    a=%-24s PD=%s" % (hi.str(10), ok2), flush=True)
            if not ok2: return lo, hi
            lo, hi = hi, hi * 2
    else:
        hi = a
        lo = a / 2
        for _ in range(40):
            ok2, _ = pd(A, lo)
            print("      down  a=%-24s PD=%s" % (lo.str(10), ok2), flush=True)
            if ok2: return lo, hi
            hi, lo = lo, lo / 2
        return None, None

out = {"window": "L~ = %d/%d" % (L_NUM, L_DEN), "L": K["L"].str(25),
       "lambda_sq": (2 * K["L"]).exp().str(25), "T": T, "N": N,
       "source": "certs/block400_lam5_{even,odd}.json (checkpointed by run_block_lam5.py)",
       "constants_path": "phase5/src/constants.py constants(prec=2000), owned",
       "beta": K["beta"].str(25), "M": K["M"].str(25), "kappa": K["kappa"].str(25),
       "d_even": K["d_even"].str(25), "d_odd": K["d_odd"].str(25), "sectors": {}}

for sec, start in (("even", "8.27e-18"), ("odd", "4e-15")):
    # R13(b): a serialized ball re-parsed OUTSIDE an owned scope carries the
    # AMBIENT precision's radius.  constants() exits owned(prec) and restores
    # ctx.prec to 53 bits, so the parse must set it explicitly first.
    ctx.prec = PREC_C
    A = [[arb(v) for v in row] for row in
         json.load(open(os.path.join(INV, "certs/block400_lam5_%s.json" % sec)))]
    d = K["d_even"] if sec == "even" else K["d_odd"]
    print("\n===== %s sector (d=%s) =====" % (sec.upper(), d.str(14)), flush=True)
    lo, hi = bracket(A, start)
    if lo is None:
        print("   NO BRACKET"); out["sectors"][sec] = None; continue
    for _ in range(18):
        mid = (lo + hi) / 2
        ok, res = pd(A, mid)
        print("      bisect a=%-24s PD=%-5s ||I-G||=%s"
              % (mid.str(12), ok, "n/a" if res is None else res.str(6)), flush=True)
        if ok: lo = mid
        else:  hi = mid
    ok, res = pd(A, lo); assert ok
    a = lo
    mlo, mhi = arb(0), a
    for _ in range(300):
        m = (mlo + mhi) / 2
        if BK.certify_final(a, d, K["kappa"], m)["ok"]: mlo = m
        else: mhi = m
    mu = arb(mlo.lower())
    print("   >>> CERTIFIED a = %s   ||I-G|| = %s" % (a.str(16), res.str(8)))
    print("   >>> mu (LOWER endpoint) = %s" % mu.str(20))
    out["sectors"][sec] = dict(a=a.str(25), rump_residual=res.str(12),
                               mu_lower_endpoint=mu.str(25))

json.dump(out, open(os.path.join(INV, "certs/certificate_T400_lam5_repaired.json"), "w"), indent=2)
print("\nwrote certs/certificate_T400_lam5_repaired.json")
