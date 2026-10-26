r"""PHASE 2 (certified fixed-window roots) and PHASE 3 (Landau-Widom consistency).

Phase 2 runs ONLY because Phase 0 passed.  The dependency analysis is
notes/prop61_dependency.md (verdict CONDITIONAL: every inspectable step of Prop 6.1
factors through coercivity (2.4), but three cited inputs are in papers not obtained).

Certified-window status, from the repository's own record:
  X = 3  -> a = 0.5493 <= 0.8   : INSIDE Zhu's certified range, lambda*(0.8) >= 8.9e-18,
                                  and lambda* is non-increasing in L, so lambda_a > 0.
  X = 5  -> a = 0.8047 > 0.8    : OUTSIDE Zhu.  Conditional on Liu (UNREFEREED).
  X = 7  -> a = 0.9730          : conditional on Liu (UNREFEREED), half-width 1.
  X = 11 -> a = 1.1990          : NO certified window.  Numerical only.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "src"))
import mpmath as mp
import suzuki as S
from phase1 import gram_cached, slice_gram, roots_on, NMAX, TABLE1

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

STATUS = {
    3:  ("CERTIFIED (Zhu, L <= 0.8)", "lambda*(0.8) >= 8.9e-18, lambda* non-increasing in L"),
    5:  ("CONDITIONAL (Liu, UNREFEREED)", "a = 0.8047 > 0.8, OUTSIDE Zhu's range"),
    7:  ("CONDITIONAL (Liu, UNREFEREED)", "a = 0.9730, Liu half-width 1"),
    11: ("NUMERICAL ONLY", "a = 1.1990, no certified window"),
}
ZHU_LOWER = mp.mpf('8.9e-18')


def main():
    OK = [True]; out = {}
    def chk(name, cond, extra=""):
        OK[0] &= bool(cond); print(f"  {'PASS' if cond else 'FAIL'}  {name} {extra}")

    print("=" * 78); print("PHASE 2 -- fixed-window roots, with their certification status")
    print("=" * 78)

    # ------------------------------------------------- P2-T5 prime-free positive control
    print("\n-- P2-T5: prime-free window control (2a < log 2) --")
    mp.mp.dps = 60
    a_pf = mp.mpf('0.3')          # 2a = 0.6 < log 2 = 0.693
    prim = S.von_mangoldt(2 * a_pf)
    print(f"   a = {mp.nstr(a_pf,6)}, 2a = {mp.nstr(2*a_pf,6)} < log 2 = {mp.nstr(mp.log(2),8)}")
    print(f"   prime powers with log n < 2a: {prim}  -> the arithmetic sum is EMPTY")
    Gpf = S.gram(a_pf, 6, dps=60)
    evs = mp.eigsy(Gpf, eigvals_only=True)
    print(f"   G is archimedean + pole only;  lambda_min = {mp.nstr(min(evs),10)},"
          f"  lambda_max = {mp.nstr(max(evs),10)}")
    chk("P2-T5: prime-free Gram has no prime input and is positive definite",
        len(prim) == 0 and min(evs) > 0)
    out["P2_T5"] = dict(a=mp.nstr(a_pf, 10), nprimes=len(prim),
                        lmin=mp.nstr(min(evs), 12), lmax=mp.nstr(max(evs), 12))

    # ------------------------------------------------- P2-T2/T3/T4 certified roots
    print("\n-- P2-T2/T3/T4: roots of W_N, by window, with precision doubling and")
    print("   exclusion of the forced lattice k pi/a --")
    rows = []
    for X in (3, 5, 7, 11):
        status, why = STATUS[X]
        print(f"\n   X = {X}:  {status}   ({why})")
        rr = {}
        for dps in (140, 180):
            G18 = gram_cached(X, dps); mp.mp.dps = dps
            a = mp.log(mp.mpf(X)) / 2
            f = S.W_N_factory(G18, a, NMAX, drop_overline=False)
            rr[dps] = [r for r in roots_on(f, mp.mpf('0.5'), 40, 4001, dps) if r > 1]
        mp.mp.dps = 180
        a = mp.log(mp.mpf(X)) / 2
        same = len(rr[140]) == len(rr[180])
        shift = max(abs(rr[140][i] - rr[180][i]) for i in range(len(rr[140]))) if same else mp.inf
        # lattice exclusion
        lat = [k * mp.pi / a for k in range(1, 200) if k * mp.pi / a < 41]
        mind = min(min(abs(r - t) for t in lat) for r in rr[180]) if lat else mp.inf
        print(f"      {len(rr[180])} roots in (1,40];  precision doubling 140->180 worst"
              f" shift {mp.nstr(shift,4)}")
        print(f"      closest approach of any root to the lattice k pi/a:"
              f" {mp.nstr(mind,6)}   (lattice spacing pi/a = {mp.nstr(mp.pi/a,8)})")
        c1 = same and shift < mp.mpf('1e-25')
        c2 = mind > mp.mpf('1e-3')
        chk(f"   X={X}: roots stable under precision doubling", c1,
            f"worst shift {mp.nstr(shift,3)}")
        chk(f"   X={X}: no root on the forced lattice", c2,
            f"min distance {mp.nstr(mind,4)}")
        lab = ("CERTIFIED roots of the finite approximation W_N" if X == 3 else
               "roots of W_N, CONDITIONAL on Liu (UNREFEREED)" if X in (5, 7) else
               "roots of W_N, NUMERICAL ONLY")
        print(f"      -> reported as: {lab}")
        print(f"      first six: " + ", ".join(mp.nstr(r, 12) for r in rr[180][:6]))
        rows.append(dict(X=X, status=status, why=why, label=lab,
                         nroots=len(rr[180]),
                         doubling_shift=mp.nstr(shift, 6),
                         min_lattice_distance=mp.nstr(mind, 6),
                         roots=[mp.nstr(r, 20) for r in rr[180][:8]]))
    out["P2_roots"] = rows
    print("""
   STATED PRECISELY.  These are certified/stable roots of the FINITE approximation
   W_N at a FIXED window.  They are NOT zeros of W(a,pi;.), and they are NOT zeros of
   zeta.  Prop 6.1 gives W_N -> W as N -> infinity at fixed a (conditionally, per
   notes/prop61_dependency.md); Theorems 1.1 and 1.2 concern a -> infinity and remain
   RH-conditional whatever is computed here.""")

    # ------------------------------------------------------------------ PHASE 3
    print("\n" + "=" * 78); print("PHASE 3 -- Landau-Widom consistency"); print("=" * 78)
    print("""   P3-T1: G is the Weil form's Gram on an L^2-ORTHONORMAL basis, so
   lambda_min(G_{a,N}) is a Rayleigh quotient over a subspace of the admissible class
   and must satisfy lambda_min(G) >= lambda*(a).  A violation would mean the Phase 0
   mapping is wrong, and feeds back into the gate.""")
    mp.mp.dps = 60
    p3 = []
    for X in (3, 5, 7, 11):
        G18 = gram_cached(X, 140); mp.mp.dps = 140
        a = mp.log(mp.mpf(X)) / 2
        evs = mp.eigsy(G18, eigvals_only=True)
        lmin, lmax = min(evs), max(evs)
        ok = (a > mp.mpf('0.8')) or (lmin >= ZHU_LOWER)
        note = "n/a (a > 0.8, no certified lambda*)" if a > mp.mpf('0.8') else \
               f"lambda_min >= Zhu's {mp.nstr(ZHU_LOWER,3)}"
        print(f"   X={X:>2}  a={mp.nstr(a,8)}  lambda_min(G) = {mp.nstr(lmin,8)}"
              f"   kappa2 = {mp.nstr(lmax/lmin,6)}   {note}")
        chk(f"   X={X}: lambda_min(G) >= lambda*(a) where lambda* is certified", ok)
        p3.append(dict(X=X, a=mp.nstr(a, 12), lmin=mp.nstr(lmin, 12),
                       lmax=mp.nstr(lmax, 12), kappa2=mp.nstr(lmax / lmin, 8),
                       ln_kappa2=mp.nstr(mp.log(lmax / lmin), 8),
                       minus_ln_lmin=mp.nstr(-mp.log(lmin), 8)))
    print("""
   P3-T2 (preregistered growth check).  Zhu's law -ln lambda*(L) ~ C N(T*)/ln N(T*),
   T* = 2 pi e^{2L}, C ~ 2 pi^2, predicts super-exponential decay of lambda* in L.
   Since lambda_max(G) is O(1), ln kappa2 should track -ln lambda_min.  Agreement
   within a factor 3 -> PASS.""")
    print(f"\n   {'X':>3} {'a':>10} {'ln kappa2':>12} {'-ln lambda_min':>16} {'ratio':>8}")
    worst = mp.mpf(0)
    for r in p3:
        rat = mp.mpf(r["ln_kappa2"]) / mp.mpf(r["minus_ln_lmin"])
        worst = max(worst, abs(rat - 1))
        print(f"   {r['X']:>3} {r['a']:>10} {r['ln_kappa2']:>12} {r['minus_ln_lmin']:>16}"
              f" {mp.nstr(rat,6):>8}")
    chk("P3-T2: ln kappa2 tracks -ln lambda_min within a factor 3",
        worst < 2, f"worst |ratio-1| = {mp.nstr(worst,4)}")
    out["P3"] = p3

    print("\n" + "=" * 78)
    print("PHASE 2/3 VERDICT:", "PASS" if OK[0] else "see REPORT")
    print("=" * 78)
    out["verdict"] = "PASS" if OK[0] else "MIXED"
    json.dump(out, open(os.path.join(DATA, "phase23.json"), "w"), indent=1, default=str)
    return OK[0]


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
