"""Certified whole-window lower bound at the window of arXiv:2611.13129v1 Table 1/2,
   row lambda^2 = 5, i.e. L* = (log 5)/2 = 0.80471895621705018730...

L is taken as the rational L~ = 80471895621705018 / 100000000000000000, which is
STRICTLY BELOW L*, so that (a) log 5 > 2 L~ keeps n = 5 out of the prime comb and
A_L~ = A_{4/5} = 2.941973525..., hence T_gate = 119.0869 and T = 400 is admissible
exactly as at L = 4/5; and (b) |L~ - L*| < 1.3e-17.

Writes only into the investigation directory.  No repository artifact is touched.
"""
import os, sys, json, time, resource
sys.path.insert(0, "phase2/src"); sys.path.insert(0, "phase4/src")
from flint import arb, arb_mat, ctx
import block as BK, baseline as B

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "certs")
L_NUM, L_DEN = 80471895621705018, 100000000000000000
T, N = 400, 440
PANELS, KORDER = 800, 64
PREC_B, PREC_M, PREC_LDL = 900, 384, 1400

def sector_dk(L, beta, M, N, parity):
    St = B.S_tail(N, arb(T) * L, top=N + 200)
    g = (2 * L * arb(T) / arb.pi()) * St
    pt2 = B.pole_vec_tail(N, L)
    pn = (2 * L.sinh()).sqrt()
    d = beta - M * g if parity == "even" else beta - M * g - 2 * pt2
    return dict(g=g, d=d, kappa=2 * pn * pt2.sqrt() + M * g.sqrt(), S=St, p_perp=pt2.sqrt())

def push(A, start):
    """largest PASSing Rump shift: double up, then 16 bisections."""
    best, trial = None, arb(start)
    for _ in range(30):
        ok, res, err = BK.verify_pd_preconditioned(A, trial, PREC_LDL)
        print("      a=%-24s PD=%-5s ||I-G||=%s" % (trial.str(12), ok,
              "n/a" if res is None else res.str(6)), flush=True)
        if ok: best, trial = trial, trial * 2
        else: break
    if best is None: return None, None
    lo, hi = best, trial
    for _ in range(16):
        mid = (lo + hi) / 2
        ok, res, err = BK.verify_pd_preconditioned(A, mid, PREC_LDL)
        print("      bisect a=%-24s PD=%-5s ||I-G||=%s" % (mid.str(12), ok,
              "n/a" if res is None else res.str(6)), flush=True)
        if ok: lo, keep = mid, res
        else:  hi = mid
    ok, res, _ = BK.verify_pd_preconditioned(A, lo, PREC_LDL)
    assert ok
    return lo, res

t0 = time.time()
ctx.prec = PREC_B
S = BK.Symbol(L_NUM, L_DEN, T)
L, beta = S.L, S.beta
print("L~ = %s" % L.str(22))
print("lambda^2 = e^{2L~} = %s   (target 5)" % (2 * L).exp().str(22))
print("comb terms (n,p) with log n < 2L~: %s" % [(n, p) for n, p in S.terms])
print("A_L~ = %s   (A_{4/5} = 2.941973525223620455503909)" % S.A_L.str(25))
print("beta_400 = %s" % beta.str(25))
assert all(n != 5 for n, _ in S.terms), "n=5 leaked into the comb -- gate invalid at T=400"
assert bool(B.Psi(arb(400), L_NUM, L_DEN) - beta >= -arb(10) ** -40), "gate fails at T=400"
print("GATE OK: Psi_L(400) >= beta_400, n=5 excluded\n", flush=True)

Ae, Ao, meta = BK.assemble(L_NUM, L_DEN, T, N, PANELS, KORDER,
                           prec_bessel=PREC_B, prec_mat=PREC_M, chunk=1600, progress=50)
print("assembled in %.0fs" % (time.time() - t0), flush=True)
M = B.M_sup(L_NUM, L_DEN, beta, T)
out = dict(L_num=L_NUM, L_den=L_DEN, L=L.str(25), lambda_sq=(2 * L).exp().str(25),
           A_L=S.A_L.str(25), beta=beta.str(25), M=M.str(25), T=T, N=N, meta=meta, sectors={})
START = {"even": "8e-18", "odd": "8e-15"}
for name, A in (("even", Ae), ("odd", Ao)):
    dk = sector_dk(L, beta, M, N, name)
    print("\n-- %s: d=%s kappa=%s" % (name, dk["d"].str(12), dk["kappa"].str(10)), flush=True)
    json.dump([[A[i][j].str(60, radius=True) for j in range(len(A))] for i in range(len(A))],
              open(os.path.join(OUT, "block400_lam5_%s.json" % name), "w"))
    print("   checkpointed block400_lam5_%s.json" % name, flush=True)
    a, res = push(A, START[name])
    rec = dict(d=dk["d"].str(25), kappa=dk["kappa"].str(25), g=dk["g"].str(25),
               S_tail=dk["S"].str(25), p_perp=dk["p_perp"].str(25))
    if a is None:
        rec["certified_shift"] = None
    else:
        lo, hi = arb(0), a
        for _ in range(200):
            mid = (lo + hi) / 2
            if BK.certify_final(a, dk["d"], dk["kappa"], mid)["ok"]: lo = mid
            else: hi = mid
        rec.update(certified_shift=a.str(25), rump_residual=res.str(12),
                   mu_lower_endpoint=arb(lo.lower()).str(25))
        print("   >>> certified a = %s   mu >= %s" % (a.str(14), arb(lo.lower()).str(20)))
    out["sectors"][name] = rec
out["elapsed_s"] = round(time.time() - t0, 1)
out["peak_rss_mib"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)
json.dump(out, open(os.path.join(OUT, "certificate_T400_lam5.json"), "w"), indent=2)
print("\nelapsed %.0fs peak RSS %.1f MiB -> certs/certificate_T400_lam5.json"
      % (out["elapsed_s"], out["peak_rss_mib"]))
