"""AUDIT Tasks 2 and 3 -- map the available constraints, and apply the DH control.

Task 2.  For each constraint: what it restricts, what freedom remains, whether it adds
         information beyond the others.  The three means count as ONE diagnostic.
Task 3.  Which constraints Davenport-Heilbronn ALSO satisfies.  Any combination drawn
         only from that set cannot, by itself, exclude off-line zeros -- because DH
         satisfies the combination and has certified off-line zeros and a certified
         NEGATIVE finite-window Weil form.

Computable entries are computed; the rest are sourced to the file that establishes them.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
from mpmath import mp, mpf, sqrt, log

mp.dps = 30
OK = [True]
def chk(n, c, x=""):
    OK[0] &= bool(c); print(f"  {'PASS' if c else 'FAIL'}  {n} {x}")

# key = (id, name, source, restricts, freedom, dh_status, block)
# dh_status: SAT = DH satisfies it; FAIL = DH provably fails it; SAT-N = numerically
K = [
 ("K0", "the mean-gap diagnostic (ONE scalar D = 2 delta^2)",
  "this ext., exact/Q",
  "nothing -- A = 1/2 IS the functional equation; G^2 = A M_H",
  "all of delta",
  "SAT", "diagnostic"),
 ("K1", "functional equation / reflection symmetry",
  "classical, T1",
  "zeros in orbits {rho, 1-rho, conj rho, 1-conj rho}",
  "delta free in (-1/2, 1/2)",
  "SAT", "foundational"),
 ("K2", "bounded zero strip (abscissa)",
  "classical / parent REPORT Sec.3, T1",
  "sup D <= 1/2 (zeta), 1.6025374836 (DH)",
  "the whole interior of the strip",
  "SAT", "strip"),
 ("K3", "de Bruijn forward bound  Lambda <= Delta_0^2/2 = sup D(0)",
  "de Bruijn 1950 + this ext. Sec.1.7, T1",
  "the parameter at which reality is GUARANTEED",
  "everything at tau = 0",
  "SAT", "flow"),
 ("K4", "extremal sign rule  D' <= -1 at the deepest zero",
  "this ext. Sec.3.5 (S1c), T1",
  "the motion of the deepest zero",
  "every non-extremal zero (refuted exactly)",
  "SAT", "flow"),
 ("K5", "proximity cutoff  |gamma_j - gamma_k| < |delta_j|  and  R < Delta/sqrt3",
  "this ext. Sec.3.5b, sharpened in this audit, T1",
  "WHICH zeros can have a growing D",
  "configurations with a nearby deeper zero",
  "SAT", "flow"),
 ("K6", "total-energy dissipation  Escr' = -2 sum (dy)^2/|dz|^2,  Lambda <= Escr(0)/2",
  "this ext. Sec.3.6 (S2), T1 conditional on m < infinity",
  "the aggregate displacement",
  "-- (weaker than K3, proved)",
  "FAIL", "flow"),
 ("K7", "delta^2-L^3 sensitivity law  lambda_min(Q_delta) >= lambda*(L) - c(L) delta^2 - O(delta^4)",
  "parent paper Thm 7 / notes/proofs.md B.4, T1 (RELATIVE, Retraction 18)",
  "how much an off-line zero can COST the window form (a LOWER bound)",
  "the sign of the baseline lambda*(L) outside L <= 0.8",
  "SAT", "window"),
 ("K8", "certified finite-window positivity for zeta, L <= 0.8: lambda*(0.8) >= 8.9e-18",
  "Zhu, cited T1",
  "the window form is positive there",
  "nothing about zero locations directly",
  "SAT-N", "window"),
 ("K9", "certified four-mode detection: Q_add,delta < 0 on [9/100, 49/100]; certified NON-detection on [0, 0.05986]",
  "parent Sec.6.1 + four-mode-positivity ext., T1 (T2-C)",
  "an ADDED (or replacement) quartet is / is not seen by one 4-mode window",
  "delta in (0.05986, 0.09); every other window and basis",
  "SAT", "window"),
 ("K10", "moderate-gap budget: depth enters only via kappa; trace delta-blind; off-line pair has ordinate distance exactly 0; log-additive",
  "parent moderate-gap ext., T1, verdict NO USABLE TRANSFER",
  "the ARTZ budget's delta-sensitivity",
  "--",
  "SAT", "detectability"),
 ("K11", "confluent cluster bound; NG-2, NG-2', NG-3",
  "parent 06_CONFLUENT_BLIND_SPOTS.md, audited GREEN / class VI",
  "detectability kappa(Gamma) under competing ordinates",
  "--",
  "SAT", "detectability"),
 ("K12", "Euler product / multiplicativity of b_n",
  "definitional for zeta; DH has none (parent Sec.4: 15 non-prime-powers below 50)",
  "the support and structure of the arithmetic side",
  "--",
  "FAIL", "distinguishing"),
 ("K13", "finitely many off-line zeros  (m < infinity)",
  "HYPOTHESIS: unknown for zeta; FALSE for DH (Davenport-Heilbronn 1936, bibliographic)",
  "enables K6",
  "--",
  "FAIL", "distinguishing"),
]

print("=" * 100)
print("TASK 2 -- THE CONSTRAINT MAP")
print("=" * 100)
for kid, name, src, restr, free, dh, blk in K:
    print(f"\n[{kid}]  {name}")
    print(f"      source     : {src}")
    print(f"      restricts  : {restr}")
    print(f"      freedom    : {free}")
    print(f"      block      : {blk}")
    print(f"      DH         : {dh}")

print("\n" + "=" * 100)
print("TASK 2b -- WHICH ADD INFORMATION BEYOND THE OTHERS")
print("=" * 100)
print("""
DIAGNOSTIC BLOCK.  K0 collapses to the single scalar D (exact over Q, this extension
Sec.1.1).  The three means are ONE diagnostic, as the brief says.  It adds NOTHING to
anything else: it is a relabelling of delta^2.  Its only non-vacuous consequence is the
dictionary identity tau_c = sup D, which names K3's threshold without narrowing it.

FLOW BLOCK {K3, K4, K5, K6}.  These are NOT four facts.
   * K4 ==> K3 by the envelope argument: integrating D' <= -1 at the deepest zero gives
     Delta_tau^2 <= Delta_0^2 - 2 tau, which IS K3.  So K3 is a consequence of K4.
   * K6 ==> a bound Lambda <= Escr(0)/2 which this extension PROVES is never better than
     K3's (Escr(0)/2 >= Delta_0^2/2 always).  So K6 adds no strength, and its hypothesis
     m < infinity is K13.
   * K5 is NOT implied by K4.  K4 is an ENVELOPE statement about one zero (the deepest);
     K5 is a PER-ZERO statement about every zero.  K4 says nothing about a zero that is
     not extremal; K5 does.
   >>> The flow block reduces to TWO independent facts: K4 (= de Bruijn) and K5 (new).

WINDOW BLOCK {K7, K8, K9}.  Each is independent of the flow block (they are about the
quadratic form, not about where zeros are).  Do they combine with each other?
   * K7 + K8 give  lambda_min(Q_delta) >= 8.9e-18 - c(L) delta^2  for L <= 0.8: a LOWER
     bound.  Excluding an off-line zero needs an UPPER bound.  Wrong direction.
   * K8 + K9: the certificate's added-quartet form and the count-preserving replacement
     differ by exactly 2 F(gamma)^2 >= 0, so Q_rep,delta <= Q_add,delta and the
     certificate covers the count-preserving convention too.  THIS COMBINATION IS REAL --
     AND IT IS ALREADY IN THE REPOSITORY (checked below).  It still does not reach RH,
     because concluding anything about zeta needs Q_rest >= 0, which is Theorem 8's
     hypothesis and is exactly what is unavailable.
   >>> The window block has no usable NEW internal combination.

DETECTABILITY BLOCK {K10, K11}.  Both are structural facts about the detector, both
already carry the verdicts NO USABLE TRANSFER and class VI in the parent extensions.
K11's NG-2' (margin exactly 0 at multiplicity >= 2) and K5 describe the SAME
configuration from two sides -- an off-line zero whose competitor shares its ordinate --
which is the one place the flow block and the detectability block touch.
   >>> Noted, and it is a coincidence of configuration, not a new estimate.

DISTINGUISHING BLOCK {K12, K13}.  The ONLY two items in the whole inventory that
separate zeta from Davenport-Heilbronn.""")

print("=" * 100)
print("TASK 2c -- the one real combination, and the check that it is NOT new")
print("=" * 100)
F_gamma = mpf("4.129358070779164e-13")
committed = mpf("3.410319615341804e-25")
margin = mpf("0.014510765749822711")
print(f"""   Q_add,delta  = Q_zeta + 4 Re F(gamma+i delta)^2          (quartet ADDED)
   Q_rep,delta  = Q_zeta - 2 F(gamma)^2 + 4 Re F(gamma+i delta)^2   (count-preserving)
   difference   = 2 F(gamma)^2 >= 0   ==>   Q_rep,delta <= Q_add,delta .

   So the certified  Q_add,delta < -1/500  also gives  Q_rep,delta < -1/500 .  More
   generally, REMOVING any on-line pairs only lowers the form (each contributes
   +2F(gamma)^2 >= 0), so the certificate covers every multiset obtained from zeta's by
   adding the quartet and deleting on-line pairs.""")
print(f"\n   frozen witness (L = 4/5, gamma = 14, dim 4): F(14) = {mp.nstr(F_gamma,12)}")
print(f"   2 F(14)^2                      = {mp.nstr(2*F_gamma**2, 16)}")
print(f"   parent's committed field\n     cold_certificate.json : 'Q_add_minus_Q_replacement' = {mp.nstr(committed,16)}")
chk("2F(14)^2 equals the parent's committed Q_add - Q_replacement",
    abs(2*F_gamma**2 - committed)/committed < mpf(10)**-14,
    f"ratio {mp.nstr(2*F_gamma**2/committed, 16)}")
print(f"   and it is {mp.nstr(2*F_gamma**2/margin, 4)} of the certified margin "
      f"{mp.nstr(margin,8)} -- utterly negligible.")
print("""
   >>> AUDIT FINDING: NOT NEW.  The parent extension already computes this quantity and
   already states the comparison, in
     extensions/planted-witness-closure-2026-10-08/CERTIFICATE_THEOREM.md  (its
       "Comparison with the count-preserving replacement convention") and
     extensions/planted-witness-consolidation-2026-10-08/VERIFICATION_REPORT.md Sec.168.
   I reached it independently in this audit and it is pre-existing.  Recorded as such.
   The brief's "whether it adds information beyond the others" is answered NO for the one
   window-block combination that works at all.""")

print("\n" + "=" * 100)
print("TASK 3 -- THE DAVENPORT-HEILBRONN FILTER")
print("=" * 100)
sat = [k for k in K if k[5].startswith("SAT")]
fail = [k for k in K if k[5] == "FAIL"]
print(f"\n   DH SATISFIES ({len(sat)} of {len(K)}):  " + ", ".join(k[0] for k in sat))
print(f"   DH FAILS     ({len(fail)} of {len(K)}):  " + ", ".join(k[0] for k in fail))
print(f"""
   K8 is marked SAT-N, not SAT: DH's own window form is numerically positive at small L
   (lambda_min = 4.03e-5 at L = 1.0, parent level-2 scan) but no certificate exists for
   it, so "small-window positivity" is satisfied by DH as far as anyone has measured and
   is certainly not a separator.

   AND DH SIMULTANEOUSLY HAS:
     * a certified NEGATIVE finite-window Weil form,
       Q_F(g_*)/||g_*||^2 = -0.679071346441480590986377586093..., Arb radius 1.4e-116;
     * four reconciled off-line quartets below T = 200, the deepest at delta = 0.3085;
     * and (this extension, certified) all four of them with D'(0) < 0.

   >>> THEREFORE: no combination drawn from {{K0,K1,K2,K3,K4,K5,K7,K8,K9,K10,K11}} can,
   >>> by itself, exclude off-line zeros.  DH satisfies every one of them and has
   >>> off-line zeros.  This is the parent source-flow audit's FILTER, applied to the
   >>> full inventory rather than to one theorem's hypotheses.

   Any future argument must use K12 (the Euler product) or K13 (finiteness of the
   off-line zero count) ESSENTIALLY.  The repository's status on each:
     K12: parent Phase H showed that splitting the arithmetic side by the support of b_n
          does NOT separate "has off-line zeros" from "has no Euler product"
          (lambda_min of the truncation is negative at 17 of 18 configurations, including
          windows where the true form is positive).  The route is closed; the property
          is not.
     K13: newly identified in this audit as the correct hypothesis for K6 (Task 1c), and
          DH fails it on the strength of the 1936 theorem alone.  But K6's conclusion is
          PROVED weaker than K3's, so K13 + the flow buys nothing.  K13 is a genuine
          separator attached to a useless conclusion.""")
chk("the DH-satisfied set contains the entire flow block except K6",
    all(k[5].startswith("SAT") for k in K if k[6] == "flow" and k[0] != "K6"))
# K6 is failed by DH only BECAUSE its hypothesis is K13; the primitive separators are two
prim = [k for k in K if k[6] == "distinguishing"]
chk("exactly two PRIMITIVE separators (the distinguishing block)", len(prim) == 2,
    f"({', '.join(k[0] for k in prim)}); K6 is failed only via K13, so it is derived")
print(f"""
   The three FAIL rows are K6, K12, K13 -- but K6 fails only BECAUSE its hypothesis is
   K13, so the PRIMITIVE separators number two: K12 and K13.  Equivalently: of the
   fourteen constraints, eleven are DH-satisfied, two are primitive separators, and one
   (K6) is a derived consequence of a separator and is weaker than something DH satisfies.""")

print("\n" + "=" * 100)
print("=== audit2_constraints.py:", "ALL PASS" if OK[0] else "SOME FAILED", "===")
json.dump(dict(constraints=[dict(id=k[0], name=k[1], source=k[2], restricts=k[3],
                                 freedom=k[4], dh=k[5], block=k[6]) for k in K],
               dh_satisfies=[k[0] for k in sat], dh_fails=[k[0] for k in fail],
               flow_block_reduces_to=["K4", "K5"],
               window_combination_preexisting=True,
               two_F_gamma_sq=mp.nstr(2*F_gamma**2, 18),
               parent_committed=mp.nstr(committed, 18)),
          open(os.path.join(os.path.dirname(__file__), "..", "data", "audit2.json"), "w"), indent=1)
sys.exit(0 if OK[0] else 1)
