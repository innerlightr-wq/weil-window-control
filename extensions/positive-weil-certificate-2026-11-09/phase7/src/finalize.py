"""Re-certify the two sharpened shifts ONCE MORE and emit outward-rounded
   endpoints plus the Hypothesis-3.4 ratios.  Rounding direction is explicit:
   a lower bound is consumed through .lower(), an upper bound through .upper().
"""
import json, os
from flint import arb, ctx
import block as BK, baseline as B

CERTS = "phase4/certs"
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
D = {"even": beta - M * g, "odd": beta - M * g - 2 * pt2}

# shifts that PASSED in the bisection, re-stated as exact decimals
SHIFT = {"odd": arb("1.3857e-14"), "even": arb("1.43300048828e-17")}
out = {}
for sec in ("even", "odd"):
    Mt = [[arb(v) for v in row] for row in
          json.load(open(os.path.join(CERTS, "block400_%s.json" % sec)))]
    a = SHIFT[sec]
    ok, res, err = BK.verify_pd_preconditioned(Mt, a, 1400)
    assert ok, (sec, err)
    lo, hi = arb(0), a
    for _ in range(200):
        mid = (lo + hi) / 2
        if BK.certify_final(a, D[sec], kappa, mid)["ok"]: lo = mid
        else: hi = mid
    out[sec] = dict(a=a, mu=lo, rump=res)
    print("%-5s  PD at a=%s  ||I-G||=%s  -> mu >= %s (lower endpt %s)"
          % (sec, a.str(14), res.str(6), lo.str(20), arb(lo.lower()).str(20)))

mu_even = arb(out["even"]["mu"].lower())
mu_odd  = arb(out["odd"]["mu"].lower())
U = arb("1.67312553407068186807726837736e-17")      # phase5/certs/upper_highN.json
U_up = arb(U.upper())
print("\n=== whole-window L=4/5, T=400, N=440 ===")
print("  mu_even (LOWER endpoint) = %s" % mu_even.str(20))
print("  mu_odd  (LOWER endpoint) = %s" % mu_odd.str(20))
print("  whole-window mu = min    = %s" % min(mu_even, mu_odd).str(20))
print("  certified upper U        = %s" % U_up.str(20))
print("  new bracket width        = %s   (old 2.64e-18)" % (U_up - mu_even).str(8))
print("  improvement over 1.41e-17: factor %s" % (mu_even / arb("1.41e-17")).str(8))

# ---- Hypothesis 3.4 ratios, with the exact prolate rungs (work/prolate_rungs.py)
d4 = arb("2.2391127080494187067e-18")
d6 = arb("3.4175764349345855425e-15")
print("\n=== Hypothesis 3.4 (arXiv:2611.13129v1 eq (32)) at lambda = e^{4/5} ===")
print("  delta_4 = %s   delta_6 = %s   (lambda^2 = %s)"
      % (d4.str(14), d6.str(14), (2*L).exp().str(14)))
print("  EVEN  eps_lambda(1)/delta_4  in  [ %s , %s ]"
      % ((mu_even / d4).str(8), (U_up / d4).str(8)))
print("  ODD   eps^odd_lambda(1)/delta_6  >=  %s" % (mu_odd / d6).str(8))
print("  (paper Table 2/4 at lambda^2=5, N=960: even 7.06, odd 4.84, both UPPER bounds)")
print("  Zhu Thm 6.2 odd bound 8.206e-15 -> ratio >= %s ; ours is %s x Zhu"
      % ((arb("8.206e-15") / d6).str(6), (mu_odd / arb("8.206e-15")).str(6)))

json.dump({
  "window": "L = 4/5 (lambda = e^{4/5}, lambda^2 = 4.95303242439511480365)",
  "source_blocks": "phase4/certs/block400_{even,odd}.json @ 7515396 (frozen, not reassembled)",
  "mu_even_lower_endpoint": mu_even.str(25),
  "mu_odd_lower_endpoint": mu_odd.str(25),
  "whole_window_mu_lower_endpoint": min(mu_even, mu_odd).str(25),
  "rump_residual_even": out["even"]["rump"].str(8),
  "rump_residual_odd": out["odd"]["rump"].str(8),
  "upper_U_outward": U_up.str(30),
  "delta_4_highprec": d4.str(20),
  "delta_6_highprec": d6.str(20),
  "ratio_even_interval": [ (mu_even/d4).str(10), (U_up/d4).str(10) ],
  "ratio_odd_lower": (mu_odd/d6).str(10),
}, open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "certs", "sharpened_T400.json"), "w"), indent=2)
print("\nwrote work/sharpened_T400.json")
