"""Repair and audit the Schur margin for every certificate produced today.

Defect found 2026-11-12 after review: the bisection drove mu to within
kappa^2/d (~1.4e-58) of a, then serialized it with .str(25), whose truncation
radius (~3e-46) is 1e12 times LARGER than the margin.  The published mu is
therefore a ball straddling a, and certify_final's strict `a > mu` is False on
replay.  The Schur step WAS performed at full precision; the artifact could not
be replayed.

Repair: publish mu as a short EXACT decimal strictly below a, with the margin
shown.  For each certificate this prints
  - a,
  - the true supremum mu* (located by bisection) and a - mu* against kappa^2/d,
  - the quotable mu, its margin, and certify_final's three signs,
  - that certify_final FAILS at mu = a.
"""
import json, os
from fractions import Fraction
from flint import arb, ctx
import block as BK
import constants as C

INV = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # phase7/
# NOTE: integrated layout -- certs/ src/ logs/ live under this phase7 directory.
# Rewritten from the session path ~/scratch/rh-connes-continuation-2026-11-12 on
# integration; no $HOME path is read.
PREC = 2000
ctx.prec = PREC
K45 = C.constants(4, 5, 400, 440, prec=PREC)
KLM = C.constants(80471895621705018, 100000000000000000, 400, 440, prec=PREC)

def q_decimal(q):
    """exact decimal string for a Fraction with denominator a power of 10"""
    from decimal import Decimal, getcontext
    getcontext().prec = 60
    return str(Decimal(q.numerator) / Decimal(q.denominator))

def sup_mu(a, d, kappa):
    lo, hi = arb(0), arb(a)
    for _ in range(400):
        m = (lo + hi) / 2
        if BK.certify_final(a, d, kappa, m)["ok"]: lo = m
        else: hi = m
    return lo

def quotable(a_str, d, kappa, rel=Fraction(1, 10**12), sig=15):
    """mu = a*(1 - rel) rounded DOWN to `sig` significant decimal digits.

    Uniform rule: the published mu gives up `rel` = 1e-12 of a relatively,
    which is ~1e30 times the true Schur margin kappa^2/d (~1e-41 of a), so it
    certifies with enormous room while staying a short exact decimal that a
    reader can quote and a replay can parse exactly.
    """
    F = Fraction(a_str) * (1 - rel)
    e = 0
    t = F
    while t >= 1: t /= 10; e += 1
    while t < Fraction(1, 10): t *= 10; e -= 1       # F = t * 10^e, 0.1 <= t < 1
    scale = 10 ** (sig - e)
    q = Fraction(int(F * scale), scale)              # truncate downward
    assert q < Fraction(a_str) and q > 0
    cand = arb(q.numerator) / arb(q.denominator)
    r = BK.certify_final(arb(a_str), d, kappa, cand)
    return q, cand, r, sig

OUT = {}
def show(key, tag, a_str, d, kappa):
    a = arb(a_str)
    mus = sup_mu(a, d, kappa)
    q, mu, r, sig = quotable(a_str, d, kappa)
    at = BK.certify_final(a, d, kappa, a)
    print("\n--- %s" % tag)
    print("    a                    = %s" % a.str(26))
    print("    kappa^2/d            = %s" % ((kappa * kappa) / d).str(8))
    print("    sup mu*  (bisected)  = %s" % mus.str(26))
    print("    a - mu*              = %s   <- the real Schur margin" % (a - mus).str(6))
    print("    QUOTABLE mu          = %s" % q_decimal(q))
    print("    a - mu (quotable)    = %s" % (a - mu).str(6))
    print("    a>mu %s  d>mu %s  (a-mu)(d-mu)>kappa^2 %s  -> ok=%s"
          % (r["a_gt_mu"], r["d_gt_mu"], r["disc_gt_kappa2"], r["ok"]))
    print("    certify_final margin = %s" % r["margin"])
    print("    at mu = a:  ok=%s   margin %s   <- must be False" % (at["ok"], at["margin"]))
    assert r["ok"] and not at["ok"], tag
    OUT[key] = dict(label=tag, a=a_str, d=d.str(25), kappa=kappa.str(25),
                    kappa2_over_d=((kappa*kappa)/d).str(12),
                    sup_mu=mus.str(30), a_minus_sup_mu=(a - mus).str(8),
                    mu_quotable=q_decimal(q), mu_quotable_exact="%d/%d" % (q.numerator, q.denominator),
                    a_minus_mu=(a - mu).str(8), certify_final_margin=r["margin"],
                    signs=dict(a_gt_mu=r["a_gt_mu"], d_gt_mu=r["d_gt_mu"],
                               disc_gt_kappa2=r["disc_gt_kappa2"], ok=r["ok"]))

print("=" * 78); print("L = 4/5"); print("=" * 78)
show("L45_even_uncond", "L=4/5 even, unconditioned", "1.43300048828e-17", K45["d_even"], K45["kappa"])
show("L45_odd_uncond",  "L=4/5 odd,  unconditioned", "1.3857e-14",        K45["d_odd"],  K45["kappa"])
ce = json.load(open(os.path.join(INV, "certs/conditioned_coercivity_L45.json")))
co = json.load(open(os.path.join(INV, "certs/conditioned_coercivity_L45_odd.json")))
show("L45_even_K4", "L=4/5 even, 4 conditions", "0.61755371093750", K45["d_even"], K45["kappa"])
show("L45_even_K5", "L=4/5 even, 5 conditions", "0.85349273681641", K45["d_even"], K45["kappa"])
show("L45_odd_K3",  "L=4/5 odd,  3 conditions", "0.10680580139160", K45["d_odd"],  K45["kappa"])
show("L45_odd_K4",  "L=4/5 odd,  4 conditions", "0.78413009643555", K45["d_odd"],  K45["kappa"])

print("\n" + "=" * 78); print("lambda_-"); print("=" * 78)
show("lam5_even", "lambda_- even", "8.27343868255615234375e-18", KLM["d_even"], KLM["kappa"])
show("lam5_odd",  "lambda_- odd",  "7.885009765625e-15",         KLM["d_odd"],  KLM["kappa"])

json.dump(OUT, open(os.path.join(INV, "certs/schur_margins.json"), "w"), indent=1)
print("\nALL 8 MARGIN CHECKS PASS.  wrote certs/schur_margins.json")
