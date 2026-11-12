"""CERTIFIED conditioned coercivity at L = 4/5 -- the rigorous counterpart of
   Fact 5.3 of arXiv:2611.13129v1, with NO prolate construction and NO E-map.

Mechanism.  For any vectors v_1..v_K and any c >= 0,

        A + c * sum_i v_i v_i^T  >=  a I        (certified by Rump)
   ==>  <A f, f> >= a ||f||^2    for every f with <v_i, f> = 0 for all i,

because the bump contributes exactly 0 on that subspace.  The v_i are frozen
EXACT RATIONALS, so the conditioning subspace is explicit; their provenance
(approximate eigenvectors of the midpoint matrix) needs no rigour at all --
only the interval test does.  By Courant-Fischer the same test also certifies
lambda_{K+1}(A) >= a.

Then the unchanged complement/coupling constants of DERIVATION.md section 5
(d = beta - M g_N, kappa) carry it to the whole window.
"""
import json, os, sys, time
from fractions import Fraction
from flint import arb, arb_mat, ctx
import block as BK, baseline as B

INV = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # phase7/
# NOTE: integrated layout -- certs/ src/ logs/ live under this phase7 directory.
# Rewritten from the session path ~/scratch/rh-connes-continuation-2026-11-12 on
# integration; no $HOME path is read.
CERTS = "phase4/certs"
PREC_LDL = 1400
DEN = 10 ** 34          # freezing denominator for the conditioning vectors
C_BUMP = arb(4)

ctx.prec = 2000
L = arb(4) / 5
T, N = 400, 440
A_L = B.A_L(4, 5)
beta = B.H(T) - A_L
M = B.M_sup(4, 5, beta, T)
St = B.S_tail(N, arb(T) * L, top=N + 200)
g = (2 * L * arb(T) / arb.pi()) * St
pt2 = B.pole_vec_tail(N, L)
kappa = 2 * (2 * L.sinh()).sqrt() * pt2.sqrt() + M * g.sqrt()
d_odd = beta - M * g - 2 * pt2   # odd sector: the pole term is carried, not dropped
print("d_odd = %s   kappa = %s   kappa^2/d = %s"
      % (d_odd.str(16), kappa.str(8), ((kappa * kappa) / d_odd).str(8)))

A = [[arb(v) for v in row] for row in
     json.load(open(os.path.join(CERTS, "block400_odd.json")))]
n = len(A)
low = json.load(open(os.path.join(INV, "src/lowvecs_odd.json")))
print("Ritz values from the numerical selection step: %s" % low["ritz"])

def freeze(vec):
    """exact rationals num/DEN, from the full decimal strings (no float step)"""
    return [int(Fraction(str(x).strip()) * DEN) for x in vec]

def bumped(A, vs, c):
    """A + c * sum v v^T  with v exact rational (num/DEN)"""
    ctx.prec = PREC_LDL
    out = [[A[i][j] for j in range(n)] for i in range(n)]
    for num in vs:
        nz = [(i, num[i]) for i in range(n) if num[i] != 0]
        sq = sum(p * p for _, p in nz)
        if sq == 0: continue
        scale = c / arb(sq)                     # normalise v to a unit vector
        for i, pi in nz:
            for j, pj in nz:
                out[i][j] = out[i][j] + scale * arb(pi * pj)
    return out

vs_all = [freeze(v) for v in low["vectors"]]
results = {}
for K in (3, 4):
    vs = vs_all[:K]
    Ab = bumped(A, vs, C_BUMP)
    print("\n===== conditioning on K = %d frozen rational directions =====" % K, flush=True)
    # bisect the largest certifiable shift in [0, 2]
    lo, hi, res_lo = arb(0), arb(2), None
    for k in range(20):
        mid = (lo + hi) / 2
        t0 = time.time()
        ok, res, err = BK.verify_pd_preconditioned(Ab, mid, PREC_LDL)
        print("   a = %-20s PD=%-5s ||I-G||=%-14s %.0fs" %
              (mid.str(10), ok, "n/a" if res is None else res.str(6), time.time() - t0),
              flush=True)
        if ok: lo, res_lo = mid, res
        else:  hi = mid
    if res_lo is None:
        print("   nothing certified"); results[K] = None; continue
    ok, res_lo, _ = BK.verify_pd_preconditioned(Ab, lo, PREC_LDL)
    assert ok
    a = lo
    mlo, mhi = arb(0), a
    for _ in range(200):
        m = (mlo + mhi) / 2
        if BK.certify_final(a, d_odd, kappa, m)["ok"]: mlo = m
        else: mhi = m
    mu = arb(mlo.lower())
    print("   >>> CERTIFIED  a = %s   ||I-G|| = %s" % (a.str(14), res_lo.str(8)))
    print("   >>> whole-window conditioned constant  mu = %s" % mu.str(16))
    print("   >>> hence also  lambda_%d(A_odd) >= %s  (Courant-Fischer)" % (K + 1, a.str(10)))
    results[K] = dict(K=K, a=a.str(25), rump=res_lo.str(12), mu=mu.str(25),
                      c_bump=C_BUMP.str(4), denominator="10^34")

json.dump(dict(window="L = 4/5", source="phase4/certs/block400_odd.json @ 7515396",
               d_odd=d_odd.str(25), kappa=kappa.str(25),
               ritz_numerical=low["ritz"], results={str(k): v for k, v in results.items()},
               vectors_frozen_denominator="10^34",
               vectors=[[str(x) for x in v] for v in vs_all]),
          open(os.path.join(INV, "certs/conditioned_coercivity_L45_odd.json"), "w"), indent=1)
print("\nwrote certs/conditioned_coercivity_L45_odd.json")
