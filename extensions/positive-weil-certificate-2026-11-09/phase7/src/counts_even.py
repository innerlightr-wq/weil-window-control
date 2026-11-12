"""CERTIFIED: at LEAST four even eigenvalues of Q below 1e-2 at L = 4/5.

Mechanism.  Win(4,5).matrix(N) is the Weil form on the first N even cosine modes
in an ORTHONORMAL basis, in closed form (digamma / Lerch Phi / Hurwitz zeta), as
Arb balls -- no quadrature.  It is a COMPRESSION of Q, so by min-max

        lambda_j( Win.matrix(N) )  >=  eps_lambda(j)    for every j,

i.e. its eigenvalues are UPPER bounds for the Weil eigenvalues.  Then for any
explicit 4-dimensional W = span{w_1..w_4},

        eps_lambda(4)  <=  lambda_4(M)  <=  lambda_max( M|_W ) ,

and lambda_max(M|_W) < c is certified by the 4x4 test  c*G - B  positive
definite, where G = W^T W and B = W^T M W.  The w_i are frozen exact rationals;
their provenance (numerical eigenvectors of the midpoint) needs no rigour.

Combined with the conditioned certificate eps_lambda(5) >= 0.6175537109 > 1e-2
(certs/conditioned_coercivity_L45.json, which gives AT MOST four), this yields
EXACTLY four even eigenvalues below 1e-2.
"""
import json, os, sys, time
from fractions import Fraction
EXT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
for p in (os.path.join(EXT, "src"), os.path.join(EXT, "phase2", "src"),
          os.path.join(EXT, "phase5", "src"),
          os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))), "src")):
    sys.path.insert(0, p)
from flint import arb, ctx
from certify import Win
from prec import owned
import block as BK

INV = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # phase7/
# NOTE: integrated layout -- certs/ src/ logs/ live under this phase7 directory.
# Rewritten from the session path ~/scratch/rh-connes-continuation-2026-11-12 on
# integration; no $HOME path is read.
N = int(sys.argv[1]) if len(sys.argv) > 1 else 160
PREC = 1000 + 12 * N
THRESH = arb("1e-2")
K = 4
DEN = 10 ** 30

t0 = time.time()
print("building Win(4,5).matrix(%d) in Arb (closed form, no quadrature)..." % N, flush=True)
with owned(PREC):
    M = Win(4, 5).matrix(N)
print("  built in %.0fs" % (time.time() - t0), flush=True)
rad = max((M[i][j].rad() for i in range(N) for j in range(N)),
          key=lambda r: float(r.mid().str(30, radius=False) or 0) if r != 0 else 0.0)
print("  max entry radius = %s   (prints as 0.0e0 in float: it underflows)" % rad.str(6))

print("locating the %d lowest eigenvectors (NUMERICAL, selects W only)..." % K, flush=True)
import mpmath as mp
mp.mp.dps = 120
Mid = mp.matrix(N, N)
for i in range(N):
    for j in range(N):
        Mid[i, j] = mp.mpf(M[i][j].mid().str(90, radius=False))
E, V = mp.eigsy(Mid)
idx = sorted(range(N), key=lambda k: E[k])[:K]
print("  lowest %d midpoint eigenvalues: %s"
      % (K, ", ".join(mp.nstr(E[k], 8) for k in idx)))
print("  5th (for context): %s" % mp.nstr(E[sorted(range(N), key=lambda k: E[k])[K]], 8))

W = [[Fraction(int(mp.nint(V[r, k] * DEN)), DEN) for r in range(N)] for k in idx]
ctx.prec = PREC
Wa = [[arb(c.numerator) / arb(c.denominator) for c in w] for w in W]
G = [[sum(Wa[p][i] * Wa[q][i] for i in range(N)) for q in range(K)] for p in range(K)]
MW = [[sum(M[i][j] * Wa[q][j] for j in range(N)) for i in range(N)] for q in range(K)]
B = [[sum(Wa[p][i] * MW[q][i] for i in range(N)) for q in range(K)] for p in range(K)]
print("\n4x4 pencil (B, G):")
for p in range(K):
    print("   B[%d] = %s" % (p, [B[p][q].str(8) for q in range(K)]))
T = [[THRESH * G[p][q] - B[p][q] for q in range(K)] for p in range(K)]
ok, resid, err = BK.verify_pd_preconditioned(T, arb(0), PREC)
print("\n  c*G - B positive definite at c = 1e-2 : %s   ||I-Gm||=%s %s"
      % (ok, "n/a" if resid is None else resid.str(8), "" if ok else (err or "")[:70]))
assert ok, "the 4-dim upper certificate FAILED"
# how far below the threshold can we push c?
lo, hi = arb(0), THRESH
for _ in range(40):
    mid = (lo + hi) / 2
    T2 = [[mid * G[p][q] - B[p][q] for q in range(K)] for p in range(K)]
    o, _, _ = BK.verify_pd_preconditioned(T2, arb(0), PREC)
    if o: hi = mid
    else: lo = mid
print("  smallest certified c with c*G - B > 0 : %s" % hi.str(12))
print("\n  => eps_lambda(4) <= lambda_4(M) <= lambda_max(M|_W) <= %s  < 1e-2" % hi.str(12))
print("  => AT LEAST 4 even eigenvalues of Q below 1e-2 at L = 4/5")
print("  with eps_lambda(5) >= 0.6175537109 (conditioned certificate: AT MOST 4)")
print("  => EXACTLY 4.")
json.dump(dict(window="L = 4/5", N=N, prec=PREC, K=K,
               mechanism="Win.matrix(N) is a compression of Q in an orthonormal basis, so "
                         "lambda_j(M) >= eps_lambda(j); lambda_max(M|_W) bounds lambda_4(M)",
               max_entry_radius=rad.str(8),
               midpoint_lowest=[mp.nstr(E[k], 20) for k in idx],
               certified_lambda4_upper=hi.str(25),
               threshold="1e-2", conclusion="at least 4 even eigenvalues below 1e-2",
               combined_with="conditioned_coercivity_L45.json K=4 -> at most 4",
               W_denominator="10^30",
               W=[[str(c) for c in w] for w in W],
               elapsed_s=round(time.time() - t0, 1)),
          open(os.path.join(INV, "certs/counts_even_L45.json"), "w"), indent=1)
print("\nwrote certs/counts_even_L45.json  (%.0fs)" % (time.time() - t0))
