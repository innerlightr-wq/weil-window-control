"""Addendum to Control 3: can the tracked low-height pair be the EXTREMAL non-real zero?

The forced-sign lemma (DERIVATION Sec.3.5, S1c) applies only to a zero of maximal distance
from the critical line.  For the emergent pair of (Z3) that hypothesis is not just unverified
-- it is almost certainly FALSE, and the reason is cheap to state: the gap between adjacent
zeta zeros shrinks like 2 pi / log(gamma/2 pi), and a close pair collides under the BACKWARD
flow at roughly tau ~ -(gap in z)^2 / 8 (the exact two-zero model).  Lehmer-type pairs at
large height therefore collide at |tau| orders of magnitude smaller, so by the time the
low-height pair has left the axis at tau ~ -1.28 there are already many non-real pairs,
some of them deeper.

This script makes that quantitative with actual zeta zeros.
"""
import sys, os, json
from mpmath import mp, mpf, zetazero

mp.dps = 30
PAIRS = [(19, 20), (13, 14), (9, 10), (6709, 6710), (6708, 6709)]

def main():
    print("== adjacent zeta zero pairs and the two-zero collision estimate ==")
    print("   exact two-zero model: P = z^2 - (g/2)^2 under d_tau P = -P'' collides at")
    print("   tau = -(g)^2/8 with g the gap in z = 2 gamma.  Used only as a SCALE.")
    print()
    print(f"   {'pair':>14} {'gamma_n':>22} {'gap in gamma':>16} {'gap in z':>12} {'tau_c ~':>14}")
    rows = []
    for a, b in PAIRS:
        ga, gb = mp.im(zetazero(a)), mp.im(zetazero(b))
        gap = gb - ga
        gz = 2*gap
        tc = -gz**2/8
        rows.append(dict(pair=f"{a}/{b}", gamma=mp.nstr(ga, 16), gap=mp.nstr(gap, 12),
                         gap_z=mp.nstr(gz, 12), tau_est=mp.nstr(tc, 10)))
        print(f"   {f'{a}/{b}':>14} {mp.nstr(ga,16):>22} {mp.nstr(gap,10):>16} "
              f"{mp.nstr(gz,10):>12} {mp.nstr(tc,8):>14}")
    print()
    print("   The Lehmer pair 6709/6710 collides at an estimated tau ~ -7e-4, against")
    print("   tau ~ -1.0 for the closest pair among the first twenty zeros: a ratio of")
    print("   more than 1000.  Under the backward flow, zeros therefore leave the axis in")
    print("   order of increasing gap, i.e. essentially at the TOP of the computed range")
    print("   first, and the low-height pair of (Z3) is one of the LAST to leave.")
    print()
    print("   CONSEQUENCE, and it is a correction to the natural reading of (Z3):")
    print("   at tau ~ -1.28 the tracked pair is certainly NOT the extremal non-real zero,")
    print("   so (S1c) does NOT apply to it.  Any observed D' <= -1 there is a numerical")
    print("   fact CONSISTENT with the lemma, not an instance of it.  Combined with the")
    print("   Davenport-Heilbronn side -- where O_1 is extremal only within a finite census")
    print("   and O_2..O_4 are provably not extremal -- the position is uniform:")
    print()
    print("     IN NEITHER CONTROL IS THE ONE SIGN-FORCING HYPOTHESIS VERIFIABLE.")
    print()
    print("   This is not a defect of the computation.  Extremality is a statement about")
    print("   ALL zeros, so no finite computation can establish it; and that is exactly why")
    print("   the forced-sign statement integrates only to the global strip bound.")
    # GATES -- assert the inequalities, do not print a stored label
    ok = True
    by = {r["pair"]: r for r in rows}
    t_lehmer = abs(mpf(by["6709/6710"]["tau_est"]))
    t_low    = abs(mpf(by["19/20"]["tau_est"]))
    g1 = mpf(by["6709/6710"]["gap"]); g2 = mpf(by["19/20"]["gap"])
    checks = [
        ("the Lehmer pair's gap is smaller than the low pair's", g1 < g2),
        ("its collision scale is smaller by a factor > 100", t_low/t_lehmer > 100),
        ("so the low-height pair is NOT the first to leave the axis", t_lehmer < t_low),
        ("the low pair's scale brackets the located collision -1.2755",
         t_low < mpf("1.2755439545316") < 2*t_low),
    ]
    print("   GATES")
    for name, c in checks:
        ok &= bool(c); print(f"     {'PASS' if c else 'FAIL'}  {name}")
    print(f"     ratio of collision scales = {mp.nstr(t_low/t_lehmer, 8)}")
    json.dump(dict(rows=rows, ratio=mp.nstr(t_low/t_lehmer, 10),
                   gates={n: bool(c) for n, c in checks}),
              open(os.path.join(os.path.dirname(__file__), "..",
              "data", "part4_zeta_lehmer.json"), "w"), indent=1)
    print("\n=== zeta_lehmer.py:", "ALL PASS" if ok else "SOME FAILED", "===")
    return ok

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
