"""How far can the FROZEN T=400 blocks be pushed?  Bisect the Rump shift a.

Nothing is reassembled: the enclosures in phase4/certs/block400_{even,odd}.json
are consumed as given.  A PASS at shift a certifies A - aI > 0, hence
lam_min(compression of R_400 to V_N) >= a for that sector.
"""
import json, os, time
from flint import arb, ctx
import block as BK, baseline as B

CERTS = "phase4/certs"
PREC_LDL = 1400

def load(sec):
    return [[arb(v) for v in row] for row in
            json.load(open(os.path.join(CERTS, "block400_%s.json" % sec)))]

def chain_mu(a, d, kappa):
    lo, hi = arb(0), a
    for _ in range(100):
        mid = (lo + hi) / 2
        if BK.certify_final(a, d, kappa, mid)["ok"]: lo = mid
        else:                                       hi = mid
    return lo

ctx.prec = 2000
L = arb(4) / 5
T, N = 400, 440
A_L = B.A_L(4, 5)
beta = B.H(T) - A_L
M = B.M_sup(4, 5, beta, T)
St = B.S_tail(N, arb(T) * L, top=N + 200)
g = (2 * L * arb(T) / arb.pi()) * St
pt2 = B.pole_vec_tail(N, L)
pn = (2 * L.sinh()).sqrt()
kappa = 2 * pn * pt2.sqrt() + M * g.sqrt()
D = {"even": beta - M * g, "odd": beta - M * g - 2 * pt2}
print("beta=%s  M=%s  kappa=%s  kappa^2/d=%s"
      % (beta.str(14), M.str(10), kappa.str(8), ((kappa*kappa)/D['even']).str(8)))

cert = json.load(open(os.path.join(CERTS, "certificate_T400.json")))
for sec in ("odd", "even"):
    Mt = load(sec)
    a0 = arb(cert["sectors"][sec]["a"])
    print("\n===== %s sector: committed a = %s =====" % (sec.upper(), a0.str(8)), flush=True)
    best, trial = None, a0
    for k in range(40):
        t0 = time.time()
        ok, res, err = BK.verify_pd_preconditioned(Mt, trial, PREC_LDL)
        print("   a = %-24s PD=%-5s ||I-G||=%-14s %.1fs %s"
              % (trial.str(12), ok, "n/a" if res is None else res.str(6),
                 time.time()-t0, "" if ok else "(" + (err or "")[:50] + ")"), flush=True)
        if ok:
            best = trial; trial = trial * 2
        else:
            break
    if best is None:
        print("   committed shift itself FAILED -- stop"); continue
    lo, hi = best, trial
    for k in range(14):
        mid = (lo + hi) / 2
        ok, res, err = BK.verify_pd_preconditioned(Mt, mid, PREC_LDL)
        print("   bisect a = %-24s PD=%-5s ||I-G||=%s"
              % (mid.str(12), ok, "n/a" if res is None else res.str(6)), flush=True)
        if ok: lo = mid
        else:  hi = mid
    mu = chain_mu(lo, D[sec], kappa)
    print("   >>> certified lam_min(block,%s) >= %s   (committed %s; factor %s)"
          % (sec, lo.str(14), a0.str(8), (lo / a0).str(6)))
    print("   >>> sector mu through the Schur chain:  %s" % mu.str(14))
