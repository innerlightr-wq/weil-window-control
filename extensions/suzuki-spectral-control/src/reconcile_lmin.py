r"""Reconcile the two lambda_min(G) values reported at X = 3.

Phase 0 reports  lambda_min = 1.042312874e-7   at a = 0.5493061443, N = 6.
Phase 1/2/3 report lambda_min = 6.9320137e-8   at the same a, N = 18.

These are NOT in conflict and neither is a slip: they are the same quantity on two
different basis sizes.  G_{a,N} is the Weil form's Gram on the L^2-orthonormal
{chi_{n,a} : |n| <= N}, a (2N+1)-dimensional subspace of the admissible class, so
lambda_min(G_{a,N}) is a Rayleigh quotient minimized over that subspace.  Enlarging
the subspace can only LOWER the minimum:

    lambda_min(G_{a,N})  is non-increasing in N,  and  >= lambda*(a)  for every N.

N = 6 gives a 13x13 matrix, N = 18 a 37x37 one, so 6.93e-8 < 1.04e-7 is required,
not merely permitted.  This script prints the whole ladder so the two published
numbers can be read off the same table, and so the monotonicity that Phase 3's
"variational upper bound" language depends on is MEASURED rather than asserted.

That the limit is exactly lambda*(a) is Connes-Consani Corollary 2.4
(arXiv:2106.01715 p.13), which is RH-free -- see notes/prop61_dependency_chain.md.
The ladder below is the finite-N shadow of that corollary.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
import mpmath as mp
import suzuki as S
from suzuki import fp

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
NS = (2, 4, 6, 8, 10, 12, 14, 16, 18)
PHASE0_N, PHASE1_N = 6, 18


def main():
    mp.mp.dps = 140
    a = mp.log(mp.mpf(3)) / 2
    print("=" * 78)
    print("RECONCILIATION -- the two lambda_min(G) values at X = 3")
    print("=" * 78)
    print(f"   a = (1/2) log 3 = {mp.nstr(a,15)},  dps 140, identical assembly throughout.")
    print(f"   Phase 0 quotes N = {PHASE0_N}; Phase 1/2/3 quote N = {PHASE1_N}.\n")
    print(f"   {'N':>4} {'dim':>5} {'lambda_min(G_(a,N))':>24} {'ratio to prev':>14}  note")
    rows, prev, ok = [], None, True
    for N in NS:
        G = S.gram(a, N, dps=140)
        lm = min(mp.eigsy(G, eigvals_only=True))
        ratio = lm / prev if prev is not None else None
        if prev is not None and not (lm < prev):
            ok = False
        note = ("<- Phase 0 quotes this" if N == PHASE0_N else
                "<- Phase 1/2/3 quote this" if N == PHASE1_N else "")
        print(f"   {N:>4} {2*N+1:>5} {mp.nstr(lm,12):>24} "
              f"{(mp.nstr(ratio,6) if ratio else '--'):>14}  {note}")
        rows.append(dict(N=N, dim=2 * N + 1, lmin=fp(lm),
                         ratio_to_prev=fp(ratio) if ratio is not None else None))
        prev = lm
    l0 = mp.mpf([r for r in rows if r["N"] == PHASE0_N][0]["lmin"])
    l1 = mp.mpf([r for r in rows if r["N"] == PHASE1_N][0]["lmin"])
    print(f"\n   {'PASS' if ok else 'FAIL'}  strictly decreasing in N at every step "
          f"(required by the variational principle)")
    print(f"   N={PHASE0_N} -> N={PHASE1_N}:  {mp.nstr(l0,10)} -> {mp.nstr(l1,10)}"
          f"   (ratio {mp.nstr(l1/l0,6)})")
    print("""
   Both numbers are VARIATIONAL UPPER BOUNDS on lambda*(0.5493); neither is a lower
   bound, and neither certifies positivity.  Positivity at this window comes from
   Zhu, not from here.""")
    json.dump(dict(a=fp(a), dps=140, rows=rows, monotone=ok,
                   phase0_N=PHASE0_N, phase1_N=PHASE1_N,
                   lmin_phase0=fp(l0), lmin_phase1=fp(l1), ratio=fp(l1 / l0)),
              open(os.path.join(DATA, "lmin_ladder.json"), "w"), indent=1)
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
