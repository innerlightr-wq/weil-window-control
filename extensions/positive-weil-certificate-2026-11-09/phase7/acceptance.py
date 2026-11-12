#!/usr/bin/env python3
"""RELEASE GATE for the two-window Weil positivity package.

Verifies the ALREADY-INTEGRATED candidate.  It applies and reverts nothing, so
it cannot silently test the pre-integration baseline, and it can be run twice.

Every path is resolved from this file's own location.  No $HOME path, no
external cache, no network.  Run it from an extracted archive and it must still
pass, provided the documented dependencies are installed.

Stages
  A  the repository's own verifier suite, plus the exact-rational closed-form
     re-certification of the L=4/5 upper bound
  B2 Rump replay of BOTH improved unconditioned L=4/5 shifts.  phase4's
     verify_t400.py tests the COMMITTED 1.42e-17 / 2.4e-17 shifts, and
     certify_final ASSUMES A >= a I rather than checking it, so stage A does
     not cover these two matrix bounds -- the release's lower endpoints rest on
     them.
  B  all eight strict Schur checks, each from its own stored a, d, kappa, mu,
     with rejection required at mu = a
  C  the four conditioned rank-K bumps, rebuilt from the stored rational vectors
  D  both lambda_- blocks, from their packaged checkpoints
  E  the lambda_- upper quotient, from its stored exact rationals
  F  directed-decimal checks and the two-window ratio, from packaged endpoints
  G  the four-dimensional even count certificate, with its Gram matrix and
     subspace interpretation checked explicitly
  H  provenance and analytic-remainder fields the certified chain depends on

Exit status is nonzero on any failure.
"""
import json, os, subprocess, sys, time
from fractions import Fraction

PH7  = os.path.dirname(os.path.abspath(__file__))
EXT  = os.path.dirname(PH7)
REPO = os.path.dirname(os.path.dirname(EXT))
CERTS = os.path.join(PH7, "certs")
for p in (os.path.join(EXT, "phase2", "src"), os.path.join(EXT, "phase4", "src"),
          os.path.join(EXT, "phase5", "src"), os.path.join(EXT, "src"),
          os.path.join(REPO, "src")):
    sys.path.insert(0, p)

FAIL, NCHK, T0 = [], [0], time.time()
def chk(name, ok, extra=""):
    NCHK[0] += 1
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name, ("  -- " + extra) if extra else ""),
          flush=True)
    if not ok: FAIL.append(name)
    return ok
def head(t):
    print("\n" + "=" * 74 + "\n" + t + "\n" + "=" * 74, flush=True)

def jload(*parts):
    return json.load(open(os.path.join(*parts)))

# ---------------------------------------------------------------- A
head("A.  repository verifier suite (frozen-artifact replay, not a cold rebuild)")
PP = os.pathsep.join([os.path.join(EXT, "phase2", "src"), os.path.join(EXT, "phase4", "src"),
                      os.path.join(EXT, "phase5", "src"), os.path.join(EXT, "src"),
                      os.path.join(REPO, "src")])
SUITE = [("verify.py", ["verify.py"]),
         ("phase2/verify_phase2.py", ["phase2/verify_phase2.py"]),
         ("phase2/verify_phase3.py", ["phase2/verify_phase3.py"]),
         ("phase3/verify_audit.py", ["phase3/verify_audit.py"]),
         ("phase4/verify_t400.py", ["phase4/verify_t400.py"]),
         ("phase5/test_prec.py", ["phase5/test_prec.py"]),
         ("note/consistency_check.py", ["note/consistency_check.py"]),
         ("phase5/src/upper_highN.py --verify", ["phase5/src/upper_highN.py", "--verify"])]
npass = 0
for name, args in SUITE:
    t = time.time()
    r = subprocess.run([sys.executable] + args, cwd=EXT, capture_output=True, text=True,
                       env=dict(os.environ, PYTHONPATH=PP))
    npass += r.stdout.count("[PASS]")
    chk("%-34s exit 0" % name, r.returncode == 0,
        "%.0fs%s" % (time.time() - t, "" if r.returncode == 0 else "  " + r.stderr[-200:]))
chk("repository suite [PASS] total", npass >= 176, "%d" % npass)

from flint import arb, ctx
import block as BK
import constants as C
from prec import owned, parse, upper as outward_upper

PREC = 2000
K45 = C.constants(4, 5, 400, 440, prec=PREC)                       # owned precision
LN, LD = 80471895621705018, 100000000000000000
KLM = C.constants(LN, LD, 400, 440, prec=PREC)
SM = jload(CERTS, "schur_margins.json")

def load_block(path):
    ctx.prec = PREC            # R13(b): name the precision BEFORE parsing balls
    return [[arb(v) for v in row] for row in json.load(open(path))]

# ---------------------------------------------------------------- B2
head("B2. the two IMPROVED unconditioned L=4/5 shifts (stage A does NOT cover these)")
for sec, key in (("even", "L45_even_uncond"), ("odd", "L45_odd_uncond")):
    A = load_block(os.path.join(EXT, "phase4", "certs", "block400_%s.json" % sec))
    a = arb(SM[key]["a"])
    ok, resid, err = BK.verify_pd_preconditioned(A, a, 1400)
    chk("L=4/5 %-4s: A - %s I is PD" % (sec, SM[key]["a"]), ok,
        "||I-G||=%s" % ("n/a" if resid is None else resid.str(8)))
    comm = arb(jload(EXT, "phase4", "certs", "certificate_T400.json")["sectors"][sec]["a"])
    chk("L=4/5 %-4s: improved shift > committed shift" % sec, bool(a > comm),
        "%s > %s" % (a.str(12), comm.str(8)))

# ---------------------------------------------------------------- B
head("B.  eight strict Schur checks from their own stored a, d, kappa, mu")
for key, rec in SM.items():
    ctx.prec = PREC
    a, d, k = arb(rec["a"]), arb(rec["d"]), arb(rec["kappa"])
    q = Fraction(rec["mu_quotable"])
    mu = arb(q.numerator) / arb(q.denominator)
    r, at = BK.certify_final(a, d, k, mu), BK.certify_final(a, d, k, a)
    chk("%-16s mu certifies, mu=a rejected, mu<a" % key,
        r["ok"] and not at["ok"] and bool(mu < a), "margin %s" % r["margin"])

# ---------------------------------------------------------------- C
head("C.  four conditioned rank-K bumps, from the stored rational vectors")
for fn, blk in (("conditioned_coercivity_L45.json", "even"),
                ("conditioned_coercivity_L45_odd.json", "odd")):
    rec = jload(CERTS, fn)
    DEN = 10 ** 34
    chk("%s records denominator 10^34" % fn, rec["vectors_frozen_denominator"] == "10^34")
    vs = [[int(x) for x in v] for v in rec["vectors"]]
    A = load_block(os.path.join(EXT, "phase4", "certs", "block400_%s.json" % blk))
    n = len(A)
    for Kstr, res in sorted(rec["results"].items()):
        Kk = int(Kstr)
        ctx.prec = 1400
        B = [row[:] for row in A]
        for num in vs[:Kk]:
            nz = [(i, num[i]) for i in range(n) if num[i] != 0]
            scale = arb(res["c_bump"]) / arb(sum(p * p for _, p in nz))
            for i, pi in nz:
                for j, pj in nz:
                    B[i][j] = B[i][j] + scale * arb(pi * pj)
        ok, resid, err = BK.verify_pd_preconditioned(B, arb(res["a"]), 1400)
        chk("%-4s K=%d: A + c*sum v_i v_i^T - a I is PD" % (blk, Kk), ok,
            "||I-G||=%s" % ("n/a" if resid is None else resid.str(8)))

# ---------------------------------------------------------------- D
head("D.  both lambda_- blocks, from their packaged checkpoints")
CL = jload(CERTS, "certificate_T400_lam5_repaired.json")
chk("lambda_- window is the rational L~ = %d/%d" % (LN, LD),
    CL["window"] == "L~ = %d/%d" % (LN, LD), CL["window"])
for sec in ("even", "odd"):
    A = load_block(os.path.join(CERTS, "block400_lam5_%s.json" % sec))
    ok, resid, err = BK.verify_pd_preconditioned(A, arb(CL["sectors"][sec]["a"]), 1400)
    chk("lambda_- %-4s block PD at the stored a" % sec, ok,
        "||I-G||=%s" % ("n/a" if resid is None else resid.str(8)))

# ---------------------------------------------------------------- E
head("E.  the lambda_- upper quotient, from its stored exact rationals")
from certify import Win
UP = jload(CERTS, "upper_lam5.json")
x = [Fraction(s) for s in UP["x"]]
with owned(UP["prec"]):
    q = Win(LN, LD).quad_form(Win(LN, LD).matrix(len(x)), x)
    hi = arb(q.upper())
chk("re-certified enclosure is inside the stored one",
    parse(UP["enclosure"], UP["prec"]).contains(q), q.str(24, radius=True))
chk("stored outward endpoint dominates",
    bool(outward_upper(parse(UP["outward_upper_endpoint"], UP["prec"]), UP["prec"]) >= hi))
chk("safe decimal %s dominates (rounds UP)" % UP["safe_decimal"],
    bool(parse(UP["safe_decimal"], UP["prec"]) >= hi))
chk("route is quadrature-free and independent of compact.py",
    "NO quadrature" in UP["route"] and UP["independent_of"] == "phase4/src/compact.py")

# ---------------------------------------------------------------- F
head("F.  directed-decimal checks and the two-window ratio")
UH = jload(EXT, "phase5", "certs", "upper_highN.json")
ctx.prec = 600
U0_full = outward_upper(parse(UH["outward_upper_endpoint"], 600), 600)
Um_full = outward_upper(parse(UP["outward_upper_endpoint"], 600), 600)
SHORT_UP = [("1.674e-17", U0_full), (UP["safe_decimal"], Um_full), ("0.003119", None)]
for s, tgt in SHORT_UP[:2]:
    chk("short UPPER %-12s dominates its endpoint" % s, bool(arb(s) >= tgt))
chk("1.6731e-17 is correctly REJECTED as an upper bound", not bool(arb("1.6731e-17") >= U0_full))
mu0 = Fraction(SM["L45_even_uncond"]["mu_quotable"])
mum = Fraction(SM["lam5_even"]["mu_quotable"])
for s, m in (("1.4330e-17", mu0), ("8.2734e-18", mum)):
    chk("short LOWER %-12s is dominated by its mu" % s, Fraction(s) <= m)
chk("0.6176 is correctly REJECTED as a short lower bound",
    not (Fraction("0.6176") <= Fraction(SM["L45_even_K4"]["mu_quotable"])))
chk("0.6175 is a valid short lower bound",
    Fraction("0.6175") <= Fraction(SM["L45_even_K4"]["mu_quotable"]))
U0, Um = Fraction(UH["safe_decimal"]), Fraction(UP["safe_decimal"])
lo, hi_r = mum / U0, Um / mu0
chk("0.494 <= mu_-/U_0", Fraction("0.494") <= lo, "%.12f" % float(lo))
chk("U_-/mu_0 <= 0.662", hi_r <= Fraction("0.662"), "%.12f" % float(hi_r))
chk("ratio bracket lies strictly in (0,1)", 0 < lo and hi_r < 1)

# ---------------------------------------------------------------- G
head("G.  the four-dimensional even count certificate (Gram + subspace)")
CT = jload(CERTS, "counts_even_L45.json")
Nc, precc = CT["N"], CT["prec"]
with owned(precc):
    Mc = Win(4, 5).matrix(Nc)
ctx.prec = precc
Wv = [[arb(Fraction(c).numerator) / arb(Fraction(c).denominator) for c in w] for w in CT["W"]]
Kc = len(Wv)
chk("W has %d columns of length %d" % (Kc, Nc), Kc == 4 and all(len(w) == Nc for w in Wv))
Gm = [[sum(Wv[p][i] * Wv[r][i] for i in range(Nc)) for r in range(Kc)] for p in range(Kc)]
ok, _, _ = BK.verify_pd_preconditioned(Gm, arb(0), precc)
chk("the Gram matrix W^T W is positive definite (W has rank 4)", ok)
MWv = [[sum(Mc[i][j] * Wv[r][j] for j in range(Nc)) for i in range(Nc)] for r in range(Kc)]
Bm = [[sum(Wv[p][i] * MWv[r][i] for i in range(Nc)) for r in range(Kc)] for p in range(Kc)]
c_up = arb(CT["certified_lambda4_upper"])
T4 = [[c_up * Gm[p][r] - Bm[p][r] for r in range(Kc)] for p in range(Kc)]
okc, _, _ = BK.verify_pd_preconditioned(T4, arb(0), precc)
chk("c*(W^T W) - W^T M W is PD, so lambda_max(M|_W) <= c", okc, "c = %s" % c_up.str(16))
chk("c is below the 1e-2 threshold", bool(c_up < arb("1e-2")), c_up.str(14))
chk("short UPPER 0.003119 dominates c", bool(arb("0.003119") >= c_up))
chk("0.0031183989 is correctly REJECTED as an upper bound for c",
    not bool(arb("0.0031183989") >= c_up))
mu5 = Fraction(SM["L45_even_K4"]["mu_quotable"])
chk("at most 4 (mu_5 >= %s > 1e-2) and at least 4 (c < 1e-2) => EXACTLY 4"
    % SM["L45_even_K4"]["mu_quotable"],
    mu5 > Fraction("1/100") and bool(c_up < arb("1e-2")))
chk("no odd count is claimed: certificate is even-sector only",
    CT["window"] == "L = 4/5" and "even" in CT["conclusion"])

# ---------------------------------------------------------------- H
head("H.  provenance and analytic-remainder fields the chain depends on")
fb = jload(EXT, "phase4", "certs", "finite_block_certificate_T400.json")["meta"]
NEED = ("panels", "korder", "panel_halflength", "b_e", "rho", "ellipse_re_semiaxis",
        "W_sup_on_ellipses", "rho_pow_2k", "quad_error_factor", "v_uniform_bound")
chk("T=400 assembly records every Bernstein-remainder parameter",
    all(k in fb for k in NEED), "panels=%s korder=%s" % (fb.get("panels"), fb.get("korder")))
ctx.prec = PREC
be, rho, qef = arb(fb["b_e"]), arb(fb["rho"]), arb(fb["quad_error_factor"])
chk("b_e = %s is strictly inside the pole-free strip |Im z| < 1/2" % be.str(6),
    bool(be > 0) and bool(be < arb(1) / 2))
chk("rho = %s > 1, so the Bernstein series converges" % rho.str(8), bool(rho > 1))
chk("the folded quadrature remainder factor is below 1e-50",
    bool(qef < arb("1e-50")) and bool(qef > 0), qef.str(10))
ctx.prec = PREC
chk("S_N tail is tight and below 1e-60 (complement estimate)",
    bool(K45["S"] < arb("1e-60")) and float(K45["S"].rad()) < 1e-70, K45["S"].str(8))
chk("kappa^2/d is below 1e-50 at both windows",
    bool((K45["kappa"] ** 2) / K45["d_even"] < arb("1e-50")) and
    bool((KLM["kappa"] ** 2) / KLM["d_even"] < arb("1e-50")),
    ((K45["kappa"] ** 2) / K45["d_even"]).str(8))
chk("A_L at L~ equals A_L at 4/5 (n=5 excluded, so T_gate is unchanged)",
    bool((K45["A_L"] - KLM["A_L"]).abs_upper() < arb("1e-20")), KLM["A_L"].str(18))
chk("the withdrawn time-domain route is excluded from every release claim",
    all("compact" not in json.dumps(jload(CERTS, f)) or
        jload(CERTS, f).get("independent_of") == "phase4/src/compact.py"
        for f in ("upper_lam5.json", "counts_even_L45.json", "schur_margins.json")))
chk("constants() was used under owned precision (require_tight passed)",
    all(float(K45[k].rad()) < 1e-20 for k in ("beta", "M", "kappa", "d_even", "d_odd")))

print("\n" + "=" * 74)
print("ACCEPTANCE: %s   (%d checks, %d failed, %.0fs)"
      % ("PASS" if not FAIL else "FAIL", NCHK[0], len(FAIL), time.time() - T0))
if FAIL:
    print("failed:")
    for f in FAIL: print("   -", f)
print("Evidence class: FROZEN-ARTIFACT REPLAY.  Every inequality above is recomputed")
print("in ball arithmetic from the packaged enclosures and exact rationals.  A COLD")
print("REBUILD of the three matrix assemblies is NOT performed here -- measured at")
print("1343 s (T=200, L=4/5), 3532 s (T=400, L=4/5) and 3678 s (T=400, L~) -- and is")
print("regenerated by phase2/src/run_block.py, phase4/src/run_block400.py and")
print("phase7/src/run_block_lam5.py respectively.")
sys.exit(1 if FAIL else 0)
