r"""PHASE 2 (certified fixed-window roots) and PHASE 3 (Landau-Widom consistency).

Phase 2 runs ONLY because Phase 0 passed.  The dependency analysis is in
notes/prop61_dependency.md and notes/prop61_dependency_chain.md.  Verdict: PASS, scoped
to X = 3.  Every inspectable step of Prop 6.1 factors through coercivity (2.4), and the
three cited inputs that were previously in papers not obtained -- [2, Lemma 2.2] and
[10, Sec 3.2] (the two form cores) and [10, Sec 3] (the form-domain characterization) --
have since been located and read from page renders and are each RH-FREE.  So at a window
where lambda_a > 0 is CERTIFIED, Prop 6.1's conclusions hold without assuming RH.

Not erased by that upgrade: Prop 2.1's third de Branges axiom does need Gamma subset R,
it simply is not invoked by Prop 6.1's chain; and nothing here touches the a -> infinity
Theorems 1.1/1.2, which need Q_W >= 0 on all of C_c^inf(R), i.e. RH itself.

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


from suzuki import fp

STATUS = {
    3:  ("CERTIFIED (Zhu, L <= 0.8)", "lambda*(0.8) >= 8.9e-18, lambda* non-increasing in L"),
    5:  ("CONDITIONAL (Liu, UNREFEREED)", "a = 0.8047 > 0.8, OUTSIDE Zhu's range"),
    7:  ("CONDITIONAL (Liu, UNREFEREED)", "a = 0.9730, Liu half-width 1"),
    11: ("NUMERICAL ONLY", "a = 1.1990, no certified window"),
}
# Literature constants.  These are CITED, not computed: Zhu publishes two significant
# figures and the repo's Gate B four/five.  They are kept as their literal decimal
# strings in the JSON.  Running them through fp() would pad them to the working
# precision and print e.g. 8.900000000000000405...e-18 -- binary noise from the dps the
# mpf happened to be built at, displayed as though Zhu had certified 140 digits.
# Full precision is for quantities we computed; cited numbers keep the digits they came
# with.  (Companion to Retraction 27, same principle, opposite direction.)
ZHU_LOWER_LIT = '8.9e-18'          # Zhu, certified LOWER bound on lambda*(0.8)
GATEB_LIT     = '2.2702e-17'       # repo Gate B, N=24 even_d basis, variational UPPER bound
GATEB2_LIT    = '2.3357e-17'       # repo Gate B, N=24 second independent basis
ZHU_LOWER = mp.mpf(ZHU_LOWER_LIT)


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
    out["P2_T5"] = dict(a=fp(a_pf), nprimes=len(prim),
                        lmin=fp(min(evs)), lmax=fp(max(evs)))

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
                         doubling_shift=fp(shift),
                         min_lattice_distance=fp(mind),
                         roots=[fp(r) for r in rr[180]]))
    out["P2_roots"] = rows
    print("""
   STATED PRECISELY.  These are certified/stable roots of the FINITE approximation
   W_N at a FIXED window.  They are NOT zeros of W(a,pi;.), and they are NOT zeros of
   zeta.  Prop 6.1 gives W_N -> W as N -> infinity at fixed a (conditionally, per
   notes/prop61_dependency.md); Theorems 1.1 and 1.2 concern a -> infinity and remain
   RH-conditional whatever is computed here.""")

    # ------------------------------------------------------------------ PHASE 3
    print("\n" + "=" * 78); print("PHASE 3 -- matched-window comparison against the repo")
    print("=" * 78)
    print("""   P3-T1: G is the Weil form's Gram on an L^2-ORTHONORMAL basis, so
   lambda_min(G_{a,N}) is a Rayleigh quotient over an N-dimensional SUBSPACE of the
   admissible class.  Restricting to a subspace RAISES the minimum, so

        lambda_min(G_{a,N})  >=  lambda*(a)        for every N,

   i.e. lambda_min(G) is a VARIATIONAL UPPER BOUND on lambda*(a), decreasing to it as
   N -> infinity.  That the limit is exactly lambda*(a) is Connes-Consani Corollary 2.4
   (arXiv:2106.01715 p.13, read from the page render), which is RH-free -- see
   notes/prop61_dependency_chain.md.  Upper and lower bounds are never interchanged
   below.""")
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
              f"   lambda_max = {mp.nstr(lmax,6)}   {note}")
        chk(f"   X={X}: lambda_min(G) >= lambda*(a) where lambda* is certified", ok)
        p3.append(dict(X=X, a=fp(a), lmin=fp(lmin), lmax=fp(lmax),
                       kappa2=fp(lmax / lmin),
                       ln_kappa2=fp(mp.log(lmax / lmin)),
                       minus_ln_lmin=fp(-mp.log(lmin))))

    # ---------------------------------------------------- P3-T2 (demoted to a sanity check)
    print("""
   P3-T2 -- SANITY CHECK ONLY, NOT A RESULT.  PREREG framed the agreement of
   ln kappa2 with -ln lambda_min as a test of Zhu's growth law.  It is not one.
   kappa2 = lambda_max/lambda_min by definition, so

        ln kappa2 = ln lambda_max - ln lambda_min,   and lambda_max ~ 3 here,

   which makes the ratio  ln kappa2 / (-ln lambda_min)  equal to
   1 + ln(lambda_max)/(-ln lambda_min) = 1 + O(1.1/40).  It is forced to be near 1 by
   arithmetic alone, for ANY matrix with a bounded top eigenvalue and a tiny bottom one.
   It tests that lambda_max is O(1).  It carries no information about the growth law,
   and it is NOT reported as the Phase 3 result.  See REPORT Retraction 28.""")
    print(f"\n   {'X':>3} {'a':>10} {'lambda_max':>12} {'ln kappa2':>12}"
          f" {'-ln lambda_min':>16} {'ratio':>9}")
    worst = mp.mpf(0)
    for r in p3:
        rat = mp.mpf(r["ln_kappa2"]) / mp.mpf(r["minus_ln_lmin"])
        worst = max(worst, abs(rat - 1))
        print(f"   {r['X']:>3} {mp.nstr(mp.mpf(r['a']),8):>10}"
              f" {mp.nstr(mp.mpf(r['lmax']),6):>12}"
              f" {mp.nstr(mp.mpf(r['ln_kappa2']),6):>12}"
              f" {mp.nstr(mp.mpf(r['minus_ln_lmin']),6):>16} {mp.nstr(rat,6):>9}")
    chk("P3-T2 (sanity): lambda_max is O(1), so ln kappa2 tracks -ln lambda_min",
        worst < 2, f"worst |ratio-1| = {mp.nstr(worst,4)}")
    out["P3_sanity_worst_ratio_dev"] = fp(worst)

    # ---------------------------------------------------- P3-T3  THE PHASE 3 RESULT
    print("\n" + "-" * 78)
    print("   P3-T3 -- THE PHASE 3 RESULT: matched-window comparison")
    print("-" * 78)
    mp.mp.dps = 140
    a5 = mp.log(mp.mpf(5)) / 2
    lmin5 = mp.mpf([r for r in p3 if r["X"] == 5][0]["lmin"])
    GATEB = mp.mpf(GATEB_LIT)
    GATEB2 = mp.mpf(GATEB2_LIT)
    print(f"""
   Three independently computed numbers in the same decade, at nearly the same window:

     this extension   lambda_min(G_(a,18)) = {mp.nstr(lmin5,8)}   at a = {mp.nstr(a5,10)}
                      closed-form Gram, Suzuki's exponential basis, dps 140
                      -> a VARIATIONAL UPPER BOUND on lambda*({mp.nstr(a5,6)})

     repo Gate B      lambda*(0.8) <= {mp.nstr(GATEB,8)}        at L = 0.8
                      (second basis: {mp.nstr(GATEB2,8)}), N = 24
                      -> also a VARIATIONAL UPPER BOUND, different code, different basis

     Zhu (certified)  lambda*(0.8) >= {mp.nstr(ZHU_LOWER,3)}             at L = 0.8
                      -> a CERTIFIED LOWER BOUND""")
    offset = a5 - mp.mpf('0.8')
    ratio  = lmin5 / GATEB
    print(f"""
   THE WINDOW OFFSET IS REAL AND IS NOT SMOOTHED OVER.
     a = {mp.nstr(a5,12)} vs L = 0.8 differ by {mp.nstr(offset,4)} ({mp.nstr(100*offset/mp.mpf('0.8'),3)}%).
     These are DIFFERENT windows, so this is a consistency comparison, not a
     reproduction.  lambda* is non-increasing in the window width, hence
         lambda*({mp.nstr(a5,6)}) <= lambda*(0.8),
     and the observed ordering {mp.nstr(lmin5,6)} < {mp.nstr(GATEB,6)} (ratio
     {mp.nstr(ratio,4)}) is in the direction monotonicity predicts.  That is a weak
     sign check, not a derivation: both quantities are upper bounds at different N,
     so their ordering is not forced.

   WHAT THIS DOES *NOT* ESTABLISH -- the Retraction-18 failure mode, named in PREREG:
     Zhu's LOWER bound is at L = 0.8 and a = {mp.nstr(a5,6)} > 0.8.  Monotonicity runs the
     WRONG WAY to carry a lower bound outward: lambda*({mp.nstr(a5,6)}) <= lambda*(0.8), so
     {mp.nstr(ZHU_LOWER,3)} is NOT a lower bound at a = {mp.nstr(a5,6)}.  Nothing here certifies
     lambda_a > 0 at X = 5.  The Phase 2 PASS is at X = 3 (a = 0.5493 < 0.8), where the
     direction does work.  X = 5 remains CONDITIONAL on Liu (UNREFEREED).

   WHAT IT DOES ESTABLISH:
     a third independent implementation -- different basis, different quadrature,
     different codebase from the repo's Gate B -- lands within a factor
     {mp.nstr(1/ratio,4)} of the repo's upper bound at a window 0.6% away, and inside
     Zhu's certified interval [{mp.nstr(ZHU_LOWER,3)}, 2.27e-17].  It is a cross-implementation
     consistency check on the assembly, at the tier of a numerical (T2) agreement.""")
    consistent = lmin5 < GATEB and lmin5 > 0
    chk("P3-T3: matched-window comparison consistent with monotonicity and Gate B",
        consistent, f"ratio lambda_min(G)/GateB = {mp.nstr(ratio,6)}")
    out["P3_matched_window"] = dict(
        a=fp(a5), L_repo="0.8", offset=fp(offset),
        lmin_extension=fp(lmin5),
        gateB_repo=GATEB_LIT, gateB_repo_basis2=GATEB2_LIT, zhu_lower=ZHU_LOWER_LIT,
        cited_constants_note="gateB_repo, gateB_repo_basis2 and zhu_lower are CITED "
                             "literature values kept at their published number of "
                             "significant figures; every other number in this file is "
                             "computed here and serialized at full working precision",
        ratio_to_gateB=fp(ratio),
        bound_direction="all three are bounds, never interchanged: extension and Gate B "
                        "are variational UPPER bounds, Zhu is a certified LOWER bound",
        lower_bound_transfers="NO -- a = 0.8047 > 0.8 and lambda* is non-increasing, so "
                              "Zhu's lower bound at 0.8 does not carry to 0.8047")
    out["P3"] = p3

    print("\n" + "=" * 78)
    print("PHASE 2/3 VERDICT:", "PASS" if OK[0] else "see REPORT")
    print("=" * 78)
    out["verdict"] = "PASS" if OK[0] else "MIXED"
    json.dump(out, open(os.path.join(DATA, "phase23.json"), "w"), indent=1, default=str)
    return OK[0]


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
