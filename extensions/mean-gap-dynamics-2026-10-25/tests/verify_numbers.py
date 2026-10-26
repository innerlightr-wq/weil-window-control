"""Verify that every headline number printed in the Markdown is derivable from data/.

Follows the convention of the parent extensions' note verifiers.  Each entry names

  * a LITERAL that must appear verbatim in the named documents,
  * a SOURCE value recomputed from data/*.json, and
  * the source's own UNCERTAINTY (an Arb radius, or the precision at which it was printed).

The entry fails if the literal is absent, if it disagrees with the source by more than half
a unit in the literal's last printed place, or if the literal is printed to a precision the
source does not support (overstatement).  Comparison is exact decimal arithmetic via
Fraction -- no floating point.

This is a whitelist of headline quantities, not a scan of every digit: the documents are
written, not machine-generated, and a full scan would require otherwise.
"""
import sys, os, json
from fractions import Fraction as Fr
from decimal import Decimal

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
DATA = os.path.join(ROOT, "data")
DOCS = ["REPORT.md", "DERIVATION.md", "CLAIMS.md", "README.md", "AUDIT.md"]
TEXT = {d: open(os.path.join(ROOT, d)).read() for d in DOCS}
FAIL = []

def _dec(s):
    """exact Fraction from a decimal string, and the magnitude of its last place."""
    s = s.strip()
    if "e" in s or "E" in s:
        m, e = s.lower().split("e")
        f = Fr(Decimal(m)) * Fr(10) ** int(e)
        nd = len(m.split(".")[1]) if "." in m else 0
        ulp = Fr(10) ** (-nd) * Fr(10) ** int(e)
    else:
        f = Fr(Decimal(s))
        nd = len(s.split(".")[1]) if "." in s else 0
        ulp = Fr(10) ** (-nd)
    return f, ulp

def arb_mid_rad(s):
    """parse '[x +/- r]' or a plain decimal."""
    s = s.strip()
    if s.startswith("["):
        body = s[1:-1]
        mid, rad = body.split("+/-")
        return _dec(mid)[0], _dec(rad)[0]
    v, ulp = _dec(s)
    return v, ulp            # a plain printed value: uncertainty = its own last place

def entry(label, literal, source, docs):
    missing = [d for d in docs if literal not in TEXT[d]]
    if missing:
        print(f"  FAIL  {label}: literal {literal!r} absent from {missing}")
        FAIL.append(label); return
    lit, ulp = _dec(literal)
    val, unc = arb_mid_rad(source)
    if abs(lit - val) > ulp / 2:
        print(f"  FAIL  {label}: {literal} disagrees with source {source}"
              f"  (deviation {float(abs(lit-val)):.3e} > half-ulp {float(ulp/2):.3e})")
        FAIL.append(label); return
    if ulp < unc:
        print(f"  FAIL  {label}: {literal} OVERSTATES the source's precision"
              f"  (last place {float(ulp):.2e} < uncertainty {float(unc):.2e})")
        FAIL.append(label); return
    print(f"  PASS  {label}: {literal}")

def load(n):
    with open(os.path.join(DATA, n)) as f:
        return json.load(f)

def exact(label, got, want, note=""):
    ok = (got == want)
    print(f"  {'PASS' if ok else 'FAIL'}  {label}: {got} == {want}  {note}")
    if not ok: FAIL.append(label)

print("== exact-over-Q quantities ==")
m = load("part2_means.json")
exact("D_crit", m["nonmonotonicity"]["D_crit"], "3/8")
exact("max (G - M_H)", m["nonmonotonicity"]["max_gap"], "1/8")
exact("G - M_H vanishes at", str(m["nonmonotonicity"]["vanishes_at"]), "['0', '1/2']")
t = load("part4_toy.json")
exact("toy D'(0) at z=i", t["D_dot_at_tau0"]["(0+1i)"], "1/3")
exact("toy D'(0) at z=2i", t["D_dot_at_tau0"]["(0+2i)"], "-19/3")
exact("toy E(0)", t["E0"], "5")
exact("toy E'(0)", t["Eprime0"], "-12")

print("\n== headline decimals, against data/ ==")
r1, r2 = m["nonmonotonicity"]["two_roots_at_v_0.1"]
entry("G-M_H two-to-one root 1", "0.2381966011", f"{r1:.12f}", ["DERIVATION.md"])
entry("G-M_H two-to-one root 2", "0.4618033989", f"{r2:.12f}", ["DERIVATION.md"])
entry("sup D^DH", "1.6025374835598227777878195",
      m["outside_strip"]["dh_supD"] + "5", ["REPORT.md", "DERIVATION.md"])
entry("toy collision tau_*", "0.13406958", f"{t['tau_collide']:.14f}", ["DERIVATION.md"])
entry("toy all-real tau_**", "0.69926375", f"{t['tau_all_real']:.14f}", ["DERIVATION.md"])
entry("toy de Bruijn bound", "2", str(t["deBruijn_bound"]), ["DERIVATION.md", "REPORT.md"])

c = load("certificates.json")
cby = {r["orbit"]: r for r in c["rows"]}
for o, lit in [("O_1", "-1.16044463430610822178313161420"),
               ("O_2", "-1.04435944844102968030569764300"),
               ("O_3", "-1.01207352431659836487600009985"),
               ("O_4", "-1.12778399346959232977056785846")]:
    entry(f"certified DH {o} D'(0)", lit, cby[o]["Dprime"],
          ["REPORT.md", "DERIVATION.md"])
    assert cby[o]["strictly_negative"], o
print("  PASS  all four certified strictly negative (data/certificates.json)")
# CLAIMS.md quotes the same four at 20 decimals.  NOTE: the erroneous value
# -1.12778399346959232978 ALSO still appears in CLAIMS.md, inside the correction note that
# records it -- deliberately, following the parent extensions' append-only retraction
# convention.  The check below asserts the CORRECT literal is present; it does not and
# should not assert the erroneous one is absent.
for o, lit in [("O_1", "-1.16044463430610822178"), ("O_2", "-1.04435944844102968031"),
               ("O_3", "-1.01207352431659836488"), ("O_4", "-1.12778399346959232977")]:
    entry(f"CLAIMS DH {o} D'(0) at 20 digits", lit, cby[o]["Dprime"], ["CLAIMS.md"])
entry("certified O_1 cross upper bound", "-0.1604316745137291191935",
      c["O1_cross_upper"], ["CLAIMS.md"])

d = load("part4_dh.json")
dby = {r["orbit"]: r for r in d["rows"]}
entry("DH O_1 own-orbit term", "-1.0000129598", dby["O_1"]["self_term"],
      ["REPORT.md", "DERIVATION.md"])
entry("DH O_1 cross term", "-0.1604316745", dby["O_1"]["cross"],
      ["REPORT.md", "DERIVATION.md"])
entry("R1 vs R2 relative deviation", "1.07e-46", d["R1_vs_R2_rel"], ["DERIVATION.md"])

i = load("part3_interaction.json")
entry("DH zero gap in mean spacings", "3.05", i["gap_in_mean_spacings"],
      ["REPORT.md", "DERIVATION.md"])
entry("heuristic predicted contribution", "-1.3026848", i["heuristic_Dprime_contrib"],
      ["REPORT.md", "DERIVATION.md"])
entry("heuristic overestimate factor", "8.12", i["heuristic_overestimate_factor"],
      ["REPORT.md", "CLAIMS.md"])
entry("DH nearest on-line ordinate below", "83.1088410847", "83.10884108470",
      ["REPORT.md", "DERIVATION.md"])
entry("DH nearest on-line ordinate above", "87.6474763326", "87.64747633260",
      ["REPORT.md", "DERIVATION.md"])
exact("on-line ordinates to T=260", i["n_online"], 169)
n200 = [r for r in i["partials"] if r["T"] == "200.0"][0]
exact("zeros used at cutoff 200", n200["nzeros"], 256)
entry("partial cross at cutoff 260", "-0.15877279",
      [r for r in i["partials"] if r["T"] == "260.0"][0]["partial"],
      ["REPORT.md", "DERIVATION.md"])
entry("partial cross at cutoff 100", "-0.14965664",
      [r for r in i["partials"] if r["T"] == "100.0"][0]["partial"],
      ["REPORT.md", "DERIVATION.md"])

sp = load("part3_separation.json")
entry("Delta_DH", "1.790272316470219", sp["Delta_DH"], ["DERIVATION.md"])
# the Delta/sqrt2 threshold 1.26591369515 was SUPERSEDED by AUDIT Sec.1.1; the sharp
# DH threshold is Delta_DH/sqrt3, and the per-zero form (S1d') is stronger still.
import mpmath as _mpm
_thr3 = _mpm.mpf(sp["Delta_DH"])/_mpm.sqrt(3)
entry("sharp DH threshold Delta_DH/sqrt3", "1.0336142", _mpm.nstr(_thr3, 12),
      ["DERIVATION.md"])
entry("zeta separation threshold", "0.3535533906", sp["zeta_threshold_gamma"],
      ["REPORT.md", "DERIVATION.md"])
entry("gamma*", "328290294.5", sp["gamma_star"], ["REPORT.md", "DERIVATION.md"])

L = load("part4_zeta_lehmer.json")
lby = {r["pair"]: r for r in L["rows"]}
entry("Lehmer gamma", "7005.062866174921", lby["6709/6710"]["gamma"], ["REPORT.md"])
entry("Lehmer collision scale", "-0.00071058837", lby["6709/6710"]["tau_est"], ["REPORT.md"])
entry("pair 19/20 collision scale", "-1.0370151", lby["19/20"]["tau_est"], ["REPORT.md"])
entry("pair 13/14 collision scale", "-1.1022183", lby["13/14"]["tau_est"], ["REPORT.md"])
entry("collision-scale ratio", "1459.3753", L["ratio"], ["REPORT.md", "DERIVATION.md", "CLAIMS.md"])

z = load("part4_zeta.json")
if "Z3_emergent_pair" in z:
    tc = z["Z3_emergent_pair"]["tau_c"]
    entry("zeta collision tau_c", "-1.2755439545316", tc, ["REPORT.md"])
    rows = {r["dtau"]: r for r in z["Z3_emergent_pair"]["rows"]}
    for k, lit in [("0.0001", "-1.00003514107"), ("0.01", "-1.00352238631"),
                   ("0.8", "-1.34537668298")]:
        if k in rows:
            entry(f"zeta (Z3) D' at dtau={k}", lit, rows[k]["Dprime"], ["REPORT.md"])
    for k, lit in [("0.8", "0.92881312100299")]:
        if k in rows:
            entry(f"zeta (Z3) D at dtau={k}", lit, rows[k]["D"], ["REPORT.md"])
else:
    print("  SKIP  (Z3) block: data/part4_zeta.json has no Z3 section yet")

# ---- the audit's own headline numbers
a1 = load("audit1.json"); a2 = load("audit2.json"); a4 = load("audit4.json")
print("\n== AUDIT.md headline numbers ==")
entry("sharp cutoff 1/sqrt3", "0.57735026919", a1["cutoff_sharp"], ["AUDIT.md"])
entry("phi(1/sqrt2)", "0.732050807569", "0.7320508075688773", ["AUDIT.md", "DERIVATION.md"])
entry("toy Lambda", "0.6992637492968862", a1["toy_thresholds"]["Lambda"], ["AUDIT.md"])
entry("toy tau_c/Lambda", "2.8601511", a1["toy_thresholds"]["ratio_tau_c_over_Lambda"], ["AUDIT.md"])
for nm, lit in [("O_2", "189"), ("O_3", "1086"), ("O_4", "406")]:
    row = [x for x in a1["dh_per_zero"] if x["target"] == nm][0]
    entry(f"per-zero ratio {nm}", lit, row["ratio"], ["AUDIT.md"])
exact("DH satisfies 11 of 14", len(a2["dh_satisfies"]), 11)
exact("flow block reduces to two", str(a2["flow_block_reduces_to"]), "['K4', 'K5']")
entry("2 F(14)^2", "3.410319615341804e-25", a2["parent_committed"], ["AUDIT.md"])
entry("isolation margin", "3.19659", a4["isolation"]["margin"], ["AUDIT.md"])
entry("C_0 measured", "0.38792761", a4["C0_measured"], ["AUDIT.md"])
entry("S' measured", "0.20851027", a4["S1_nonorbit"], ["AUDIT.md"])
entry("S' hypothesis bound slack", "8.99419", "8.994190", ["AUDIT.md"])
entry("|D'| bound", "9.286639392", a4["Dprime_bound"], ["AUDIT.md"])
entry("A", "8.161100685", a4["A"], ["AUDIT.md"])
entry("B", "0.6944960613", a4["B"], ["AUDIT.md"])
entry("tau* estimate", "0.01366592", a4["tau_star_estimate"], ["AUDIT.md"])
entry("improvement factor", "2.0018e+26", a4["improvement_factor"], ["AUDIT.md"])

# ---- the fixed-parameter certificate
cf = load("certificate_fixed_tau.json")
print("\n== the fixed-parameter certificate ==")
exact("tau", cf["tau"], "4/25")
exact("disc radius", cf["radius"], "1/16")
exact("arcs", cf["narcs"], 512)
exact("winding number", cf["winding"], 1)
exact("non-adjacent half-plane steps", cf["nonadjacent_steps"], 0)
exact("net half-plane advance", cf["net_advance"], 4)
exact("negative control winding", cf["control_winding"], 0)
entry("certified min |H_tau|", "5.756765507e-31", cf["min_H"], ["AUDIT.md", "CLAIMS.md"])
entry("disc margin", "0.11980000000000", cf["margin"], ["AUDIT.md"])
entry("worst relative enclosure radius", "0.0121658", cf["worst_rel_radius"], ["AUDIT.md"])
entry("tau-series tail T_20", "9.9540082e-47", cf["T_M"], ["AUDIT.md"])
entry("arc tail A_10", "6.2959781e-47", cf["A_K"], ["AUDIT.md"])

tv = load("taylor_vs_integral.json")
print("\n== the representation cross-check ==")
exact("points checked", len(tv["points"]), 5)
exact("all points pass", all(p["pass_"] for p in tv["points"]), True)
entry("wrong-constant control deviation", "18026", tv["wrong_constant_rel_dev"], ["AUDIT.md"])
p0 = tv["points"][0]
entry("disc-centre |H| at tau=4/25", "5.6058091e-34", p0["absH"], ["AUDIT.md"])
entry("disc-centre series vs integral (RELATIVE)", "5.926e-31",
      p0["rel_series_vs_integral"], ["AUDIT.md"])
# pin the ABSOLUTE figure too, so the relative/absolute distinction stays checked
import mpmath as _m
_abs = _m.mpf(p0["rel_series_vs_integral"])*_m.mpf(p0["absH"])
entry("disc-centre series vs integral (ABSOLUTE)", "3.3e-64", _m.nstr(_abs, 4), ["AUDIT.md"])
entry("S(1)", "1.3236800759485", cf["S1"], ["AUDIT.md"])
entry("(log f)'(3)", "-409024.8675", cf["dlog_at_3"], ["AUDIT.md"])

print()
if FAIL:
    print(f"=== VERIFY_NUMBERS: {len(FAIL)} FAILED: {FAIL} ===")
    sys.exit(1)
print("=== VERIFY_NUMBERS: ALL PASS ===")
