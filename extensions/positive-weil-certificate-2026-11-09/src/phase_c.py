"""Phase C: L = 81/100, the next prime activation, and the measured deficit.

Two independent things are computed.

(1) THE METHOD'S COST, exactly.  The one-stroke reduction needs T# above the
    root of its gate function, and that root is driven by the comb mass
    A_L = sum_{log n < 2L} 2 Lambda(n)/sqrt n.  A_L is a STEP function of L with
    a jump of exactly 2 log p / sqrt(p^?) at each L = (log n)/2, so the required
    T# ~ 2 pi e^{A_L} jumps by a MULTIPLICATIVE factor there.  Crossing
    L_5 = (log 5)/2 multiplies it by exactly exp(2 log 5 / sqrt 5).

(2) AN UPPER BOUND on the true window minimum, from the repository's exact
    Weil-form assembly.  For any finite subspace V_N,

        lambda_min(L) = inf_{f in H} Q(f)/||f||^2  <=  inf_{f in V_N} Q(f)/||f||^2
                      = lambda_min(M_N) ,

    so the finite Gram matrix gives an UPPER bound on the whole-window constant,
    never a lower one.  A NEGATIVE lambda_min(M_N) would REFUTE whole-window
    positivity at that L; a positive one proves nothing about the complement.
    This is the (ii) gap of the brief and the code labels it as such.
"""
import os, sys, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "src"))
import mpmath as mp
import zeta_window as zw
from flint import arb
from wcore import set_prec, comb_mass, ball_str

CHECKS = []
def chk(name, ok, extra=""):
    CHECKS.append((name, bool(ok), extra))
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name, ("  -- " + extra) if extra else ""))
    return ok

def activation_table():
    """L at which each prime power switches on, exactly: L_n = (log n)/2."""
    print("\n-- prime activation thresholds L_n = (log n)/2")
    rows = []
    for n in (2, 3, 4, 5, 7, 8, 9):
        Ln = arb(n).log() / 2
        rows.append((n, Ln))
        print("     n=%-2d  L_n = %s" % (n, Ln.str(14)))
    L5 = arb(5).log() / 2
    chk("C1 L_5 = (log 5)/2 lies strictly between 4/5 and 81/100",
        bool(L5 > arb(4) / 5) and bool(L5 < arb(81) / 100),
        "4/5 = 0.8 < %s < 0.81" % L5.str(12))
    return rows

def cost_jump():
    A8, _ = comb_mass(8, 5)
    A81, _ = comb_mass(162, 100)
    jump = A81 - A8
    exact = 2 * arb(5).log() / arb(5).sqrt()
    chk("C2 A_L jump across L_5 equals exactly 2 log 5 / sqrt 5",
        bool((jump - exact).abs_upper() < arb(10) ** -30), ball_str(jump, 12))
    fac = exact.exp()
    print("     multiplicative cost factor on T# : exp(2 log 5/sqrt 5) = %s" % fac.str(14))
    # Zhu used T# = 200 at L = 4/5 with gate root 120.0824; keep the same safety ratio.
    ratio = arb(200) / arb("120.082403914464154")
    T81 = ratio * arb("503.388460910693812")
    print("     same safety ratio %s  =>  T# (L=81/100) ~ %s" % (ratio.str(8), T81.str(8)))
    # N ~ e L T# / 2 (paper's super-exponential cutoff for the Bessel/Legendre entries)
    N8 = arb(1).exp() * arb(4) / 5 * arb(200) / 2
    N81 = arb(1).exp() * arb(81) / 100 * T81 / 2
    print("     mode count   N ~ e L T#/2 :  L=4/5 -> %s      L=81/100 -> %s"
          % (N8.str(6), N81.str(6)))
    cost = (N81 / N8) ** 2 * (T81 / arb(200))
    print("     assembly cost ~ N^2 T#    :  factor %s" % cost.str(6))
    chk("C3 cost factor crossing L_5 exceeds 50x", bool(cost > 50), cost.str(6))
    return fac, T81, N81, cost

def gram_upper_bounds(Ls, Nlist, dps=60):
    """lambda_min of the exact Weil Gram matrix on the first N EVEN modes."""
    out = {}
    for name, L in Ls:
        out[name] = []
        for N in Nlist:
            t0 = time.time()
            mp.mp.dps = dps
            M = zw.build_matrix(mp.mpf(L[0]) / L[1], N, dps, sector=zw.EVEN)
            vals, _ = __import__("weilform").eigsym(M) if False else (None, None)
            E = mp.eigsy(M, eigvals_only=True)
            lam = min(E)
            dt = time.time() - t0
            out[name].append((N, mp.nstr(lam, 12), dt))
            print("     L=%-8s N=%-3d  lambda_min(M_N) = %-24s  (%.1fs)"
                  % (name, N, mp.nstr(lam, 12), dt))
    return out

def main():
    set_prec(400)
    print("=== PHASE C: L = 81/100, prime activation, measured deficit ===")
    activation_table()
    print("\n-- the method's cost across the activation")
    fac, T81, N81, cost = cost_jump()

    print("\n-- UPPER bounds on lambda_min(L) from the exact finite Gram matrix")
    print("   (inf over a subspace >= inf over the space, so these bound lambda_min ABOVE)")
    res = gram_upper_bounds([("4/5", (4, 5)), ("81/100", (81, 100))], [4, 8, 12, 16], dps=60)

    # the decisive question: does any finite block go NEGATIVE?
    neg = {k: [r for r in v if mp.mpf(r[1]) < 0] for k, v in res.items()}
    for k, v in neg.items():
        chk("C4 L=%s: no finite block found with lambda_min < 0 (no refutation)" % k,
            len(v) == 0, "%d negative blocks of %d tested" % (len(v), len(res[k])))
    chk("C5 finite lambda_min is an UPPER bound only; whole-window status UNRESOLVED",
        True, "the (ii) infinite-dimensional complement is NOT discharged here")
    out = {"activation_L5": str(arb(5).log() / 2),
           "A_L_jump_exact": "2 log 5 / sqrt 5",
           "T_sharp_cost_factor": fac.str(16),
           "gram_upper_bounds": {k: [(n, s) for n, s, _ in v] for k, v in res.items()}}
    d = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
    os.makedirs(d, exist_ok=True)
    json.dump(out, open(os.path.join(d, "phase_c.json"), "w"), indent=2)
    bad = [n for n, ok, _ in CHECKS if not ok]
    print("\nPHASE C: %d checks, %d failures" % (len(CHECKS), len(bad)))
    return 1 if bad else 0

if __name__ == "__main__":
    sys.exit(main())
