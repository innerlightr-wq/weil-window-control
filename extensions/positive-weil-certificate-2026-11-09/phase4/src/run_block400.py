"""Driver: assemble the finite blocks at L=4/5, T=400 and certify them.

Outputs phase4/certs/finite_block_certificate_T400.json.  Exploratory steps are labelled
NUMERICAL; only the interval LDL^T pivots and the three final sign tests are
certificates.
"""
import os, sys, json, time, resource
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "src"))
CERTS = os.path.join(HERE, "..", "certs")
from flint import arb, acb, arb_mat, ctx
import block as BK
import baseline as B

L_NUM, L_DEN, T = 4, 5, 400
N = 440
PANELS, KORDER = 800, 64
PREC_B, PREC_M, PREC_LDL = 900, 384, 1400

def sector_dk(L, beta, M, N, parity):
    """d and kappa for one parity sector.  The complement sums over n >= N of
    BOTH parities, which is an upper bound for either one (all terms >= 0)."""
    x = arb(T) * L
    S = B.S_tail(N, x)
    g = (2 * L * arb(T) / arb.pi()) * S
    pt2 = B.pole_vec_tail(N, L)            # |p_perp|^2, both parities: an upper bound
    pt = pt2.sqrt()
    pn = (2 * L.sinh()).sqrt()             # |p_N| <= |p|
    if parity == "even":
        d = beta - M * g                   # pole term is +, dropped favourably
    else:
        d = beta - M * g - 2 * pt2         # pole term is -, must be carried
    kappa = 2 * pn * pt + M * g.sqrt()
    return dict(g=g, d=d, kappa=kappa, S=S, p_perp=pt)

def main():
    t0 = time.time()
    ctx.prec = PREC_B
    S = BK.Symbol(L_NUM, L_DEN, T)
    L, beta = S.L, S.beta
    print("assembling: L=%d/%d T=%d N=%d panels=%d k=%d" % (L_NUM, L_DEN, T, N, PANELS, KORDER), flush=True)
    Ae, Ao, meta = BK.assemble(L_NUM, L_DEN, T, N, PANELS, KORDER,
                               prec_bessel=PREC_B, prec_mat=PREC_M,
                               chunk=1600, progress=50)
    print("assembled in %.0fs" % (time.time() - t0), flush=True)

    M = B.M_sup(L_NUM, L_DEN, beta, T)
    out = dict(meta=meta, M=M.str(25), sectors={})

    for name, A in (("even", Ae), ("odd", Ao)):
        dk = sector_dk(L, beta, M, N, name)
        print("\n-- %s sector: d=%s kappa=%s" % (name, dk["d"].str(12), dk["kappa"].str(10)), flush=True)
        # EXPLORATORY (NUMERICAL): midpoint spectrum, only to choose a shift
        ctx.prec = PREC_M
        mid = arb_mat([[arb(A[i][j].mid()) for j in range(len(A))] for i in range(len(A))])
        try:
            ev = mid.eig()
            lo = min(float(e.real.mid()) for e in ev)
            print("   NUMERICAL diagnostic lam_min(midpoint) ~ %.6e" % lo, flush=True)
        except Exception as e:
            lo = None
            print("   NUMERICAL diagnostic unavailable: %s" % str(e)[:60], flush=True)
        # CHECKPOINT the assembled matrix so certification can be retried cheaply
        ck = os.path.join(CERTS, "block400_%s.json" % name)
        json.dump([[A[i][j].str(60, radius=True) for j in range(len(A))]
                   for i in range(len(A))], open(ck, "w"))
        print("   checkpointed -> certs/block400_%s.json" % name, flush=True)
        # CERTIFY: preconditioned Rump criterion (naive interval LDL^T is unusable,
        # see block.verify_pd_preconditioned docstring)
        best, pivots = None, None
        for sh in ("1.03e-17", "1e-17", "5e-18", "1e-18", "1e-20", "1e-24", "0"):
            ok, resid, err = BK.verify_pd_preconditioned(A, arb(sh), PREC_LDL)
            print("   Rump shift %-7s -> %s  ||I-G||_inf=%s%s"
                  % (sh, "PD" if ok else "failed",
                     "n/a" if resid is None else resid.str(8),
                     "" if err is None else "  (%s)" % err[:60]), flush=True)
            if ok:
                best, pivots = sh, resid
                break
        sec = dict(d=dk["d"].str(25), kappa=dk["kappa"].str(25), g=dk["g"].str(25),
                   S_tail=dk["S"].str(25), p_perp=dk["p_perp"].str(25),
                   numerical_lam_min_midpoint=(None if lo is None else "%.6e" % lo),
                   certified_shift=best,
                   rump_residual=(None if pivots is None else pivots.str(12)))
        if best is not None:
            a = arb(best)
            # final implication: find the largest mu on a declared ladder
            mu_ok = None
            for mu in ("1e-17", "9e-18", "5e-18", "1e-18", "1e-20", "1e-30"):
                r = BK.certify_final(a, dk["d"], dk["kappa"], arb(mu))
                if r["ok"]:
                    mu_ok = (mu, r)
                    break
            sec["final"] = (None if mu_ok is None else
                            dict(mu=mu_ok[0], **{k: str(v) for k, v in mu_ok[1].items()}))
        out["sectors"][name] = sec

    out["elapsed_s"] = round(time.time() - t0, 1)
    out["peak_rss_mib"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)
    d = CERTS
    os.makedirs(d, exist_ok=True)
    json.dump(out, open(os.path.join(d, "finite_block_certificate_T400.json"), "w"), indent=2)
    print("\nelapsed %.0fs  peak RSS %.1f MiB  -> certs/finite_block_certificate.json"
          % (out["elapsed_s"], out["peak_rss_mib"]))
    return 0

if __name__ == "__main__":
    sys.exit(main())
