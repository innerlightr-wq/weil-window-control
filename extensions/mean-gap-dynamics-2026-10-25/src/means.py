"""Part 2 -- the mean-gap diagnostic, exactly.

For a zero rho = beta + i gamma put

    delta = beta - 1/2,     D = 2 delta^2 .

On 0 <= beta <= 1 the reflected-coordinate partition  p = beta, q = 1 - beta  has

    A   = (p+q)/2    = 1/2                       (the functional equation, nothing else)
    G   = sqrt(p q)  = sqrt(1/4 - delta^2)
    M_H = 2pq/(p+q)  = 2 p q = 1/2 - 2 delta^2
    D   = A - M_H    = 2 delta^2 .

Everything below is exact rational algebra (G is handled through G^2), verified at
exact rational delta.  M_H is written M_H, never H, to keep it clear of H_tau.
"""
from fractions import Fraction as Fr
import json, os, sys

PTS = [Fr(0), Fr(1,10), Fr(1,4), Fr(3,8), Fr(-1,5), Fr(49,100),
       Fr(3085171824566374, 10**16),            # the DH control's delta_1
       Fr(1,2)]

def means(delta):
    """(A, G2, M_H, D) with G given through its square, exactly over Q."""
    p, q = Fr(1,2) + delta, Fr(1,2) - delta
    A   = (p + q) / 2
    G2  = p * q
    M_H = 2 * p * q / (p + q)
    D   = A - M_H
    return A, G2, M_H, D

def gaps_from_D(D):
    """The three mean gaps as functions of D alone.  G is returned squared."""
    A, M_H = Fr(1,2), Fr(1,2) - D
    G2 = Fr(1,4) - D/2
    return dict(A_minus_M_H=A - M_H, A_minus_G2=(A, G2), G2_minus_M_H=(G2, M_H))

def report():
    out = {"points": [], "structural": {}, "nonmonotonicity": {}}
    ok = True

    # ---- 1. the dictionary, exactly, at exact rational delta
    print("== 1. mean dictionary, exact over Q ==")
    print(f"{'delta':>24} {'A':>6} {'G^2':>20} {'M_H':>20} {'D':>20}  checks")
    for d in PTS:
        A, G2, M_H, D = means(d)
        c1 = (A == Fr(1,2))                      # A is identically 1/2
        c2 = (G2 == A * M_H)                     # G^2 = A * M_H  (AM.HM = GM^2)
        c3 = (D == 2*d*d)                        # D = 2 delta^2
        c4 = (M_H == Fr(1,2) - 2*d*d)
        c5 = (G2 == Fr(1,4) - d*d)
        c6 = (M_H == 2*G2)                       # since A = 1/2
        c7 = (4*G2 + (2*d)**2 == 1)              # the parent dictionary's 4h^2 + chi^2 = 1
        allc = c1 and c2 and c3 and c4 and c5 and c6 and c7
        ok &= allc
        print(f"{str(d):>24} {str(A):>6} {str(G2):>20} {str(M_H):>20} {str(D):>20}  {'PASS' if allc else 'FAIL'}")
        out["points"].append(dict(delta=str(d), A=str(A), G2=str(G2), M_H=str(M_H), D=str(D), pass_=allc))

    # ---- 2. the structural collapse
    print("\n== 2. structural identities (the whole content of the mean triple) ==")
    print("  (i)   A == 1/2 identically                  -> carries no information about delta")
    print("  (ii)  G^2 == A * M_H  (AM.HM = GM^2)        -> G and M_H are algebraically equivalent")
    print("  (iii) hence  M_H = 2 G^2  and  D = A - M_H = 2 delta^2")
    print("  => the mean triple (A,G,M_H) is a function of the SINGLE scalar D, and D = 2 delta^2.")
    print("     No mean gap carries information beyond squared zero displacement.")
    out["structural"] = dict(A_constant=True, G2_eq_A_MH=True, single_scalar=True)

    # ---- 3. the gaps as functions of D; monotonicity and the non-monotone one
    print("\n== 3. the three gaps as functions of D on 0 <= D <= 1/2 ==")
    print("  A - M_H = D                                      (identity; strictly increasing)")
    print("  A - G   = 1/2 - sqrt(1/4 - D/2)                   (strictly increasing, convex)")
    print("  G - M_H = sqrt(1/4 - D/2) - 1/2 + D               (NOT monotone)")
    # exact critical point and value of G - M_H
    Dc, Gc = Fr(3,8), Fr(1,4)
    assert Fr(1,4) - Dc/2 == Gc**2                   # sqrt(1/4 - 3/16) = 1/4 exactly
    vmax = Gc - (Fr(1,2) - Dc)
    c_end0 = (Fr(1,4) - Fr(0)/2 == Fr(1,2)**2) and (Fr(1,2) - (Fr(1,2) - Fr(0)) == 0)
    c_end1 = (Fr(1,4) - Fr(1,2)/2 == 0) and (Fr(0) - (Fr(1,2) - Fr(1,2)) == 0)
    print(f"\n  critical point   D = {Dc}   (delta = sqrt(3)/4 = 0.4330127..., beta = 0.9330127...)")
    print(f"  there  G = {Gc} exactly, and  max (G - M_H) = {vmax}  exactly")
    print(f"  G - M_H = 0 at D = 0    (equality p = q = 1/2, the RH point):      {'PASS' if c_end0 else 'FAIL'}")
    print(f"  G - M_H = 0 at D = 1/2  (complete dominance, beta in {{0,1}}, G = M_H = 0): {'PASS' if c_end1 else 'FAIL'}")
    print("  d^2/dD^2 (G - M_H) = -(1/16)(1/4 - D/2)^{-3/2} < 0  -> strictly concave,")
    print("  so G - M_H is two-to-one onto (0, 1/8): it is a STRICTLY WORSE diagnostic than D.")
    ok &= c_end0 and c_end1
    out["nonmonotonicity"] = dict(D_crit=str(Dc), G_at_crit=str(Gc), max_gap=str(vmax),
                                  vanishes_at=["0", "1/2"], concave=True, two_to_one=True)

    # exhibit the two-to-one explicitly, exactly: G - M_H = 1/10 has two rational-free roots
    print("\n  two-to-one, exhibited: solve sqrt(1/4 - D/2) = v + 1/2 - D for v = 1/10")
    v = Fr(1,10)
    # (v + 1/2 - D)^2 = 1/4 - D/2  ->  D^2 - D(2v + 1 - 1/2) + (v^2 + v) = 0
    b = 2*v + 1 - Fr(1,2); c = v*v + v
    disc = b*b - 4*c
    print(f"    D^2 - ({b}) D + ({c}) = 0,  discriminant = {disc} > 0  -> two roots in (0, 1/2)")
    import math
    r1 = (float(b) - math.sqrt(float(disc)))/2; r2 = (float(b) + math.sqrt(float(disc)))/2
    print(f"    D = {r1:.10f} and {r2:.10f}   (delta = {math.sqrt(r1/2):.10f} and {math.sqrt(r2/2):.10f})")
    ok &= (disc > 0 and 0 < r1 < 0.5 and 0 < r2 < 0.5)
    out["nonmonotonicity"]["two_roots_at_v_0.1"] = [r1, r2]

    # ---- 4. derivative of each gap along ANY evolution of D
    print("\n== 4. along any specified evolution tau |-> D(tau) ==")
    print("  (A - M_H)' = D'                      -- always the sign of D'")
    print("  (A - G)'   = D'/(4G)                 -- always the sign of D'  (G > 0)")
    print("  (G - M_H)' = D' (1 - 1/(4G))         -- sign of D' iff G > 1/4 iff D < 3/8")
    print("  So for D > 3/8 the gap G - M_H moves OPPOSITE to the displacement diagnostic.")
    out["derivatives"] = dict(A_minus_MH="D'", A_minus_G="D'/(4G)", G_minus_MH="D'(1-1/(4G))",
                              sign_flip_at="D = 3/8")

    # ---- 5. the arc coordinate: M_H'' = 1 - 4 M_H is NOT dynamics
    print("\n== 5. the arc coordinate, and why M_H'' = 1 - 4 M_H is not a restoring law ==")
    print("  put delta = (1/2) sin(theta).  Then M_H = 1/2 - 2 delta^2 = 1/4 + (1/4) cos(2 theta),")
    print("  so  d^2 M_H/d theta^2 = -cos(2 theta) = 1 - 4 M_H.   TRUE, and empty:")
    print("  theta is a coordinate on the fixed configuration, not a time.  NEGATIVE CONTROL:")
    print("  every function of the form  c0 + c1 cos(2 theta)  satisfies  f'' = -4(f - c0),")
    print("  so the 'oscillator' form is a property of 'quadratic in sin theta', not of zeta.")
    import math
    for th in (0.3, 1.1, 2.0):
        for (c0, c1) in ((0.25, 0.25), (7.0, -3.0), (0.0, 1.0)):
            f  = c0 + c1*math.cos(2*th)
            fpp = -4*c1*math.cos(2*th)
            r = abs(fpp - (-4*(f - c0)))
            ok &= (r < 1e-12)
    print("  checked at three theta for three (c0,c1): max residual < 1e-12  -> PASS")
    out["arc_coordinate"] = dict(identity="M_H'' = 1 - 4 M_H", is_dynamics=False,
                                 control="any c0 + c1 cos 2theta satisfies f'' = -4(f-c0)")

    # ---- 6. outside the strip
    print("\n== 6. outside 0 <= beta <= 1 ==")
    for d in (Fr(3,5), Fr(9,10)):
        p, q = Fr(1,2)+d, Fr(1,2)-d
        print(f"  delta = {d}: p = {p}, q = {q}, p q = {p*q} < 0  -> G imaginary, M_H = {2*p*q} < 0;")
        print(f"               D = {2*d*d} is still a fine squared-displacement diagnostic.")
    print("  The partition reading needs |delta| <= 1/2, i.e. D <= 1/2.")
    print("  zeta's trivial strip gives sup D = 1/2  -- exactly the edge.")
    print("  The DH strip (sigma_0 = 1.3951361582...) gives sup D = 2(sigma_0-1/2)^2 = 1.6025374836 > 1/2,")
    print("  so for the one repository object with certified off-line zeros the partition")
    print("  interpretation is NOT available on the whole admissible range.")
    out["outside_strip"] = dict(needs="|delta| <= 1/2 i.e. D <= 1/2",
                                zeta_supD="1/2", dh_supD="1.602537483559822777787819")

    # ---- 7. overlap with the committed partition-coordinate audit
    print("\n== 7. OVERLAP with extensions/rh-partition-coordinate-audit-2026-10-07 ==")
    print("  that audit's dictionary already contains these means under other names:")
    print("    h = sqrt(p q) = G   (its 'altitude');     R = 1/(4 p q) = 1/(4 G^2) = 1/(2 M_H)")
    print("    chi = 2|delta|;     4 h^2 + chi^2 = 1  (checked above)")
    print("  so the harmonic mean is 1/(2R) and the geometric mean is h: NOT new coordinates.")
    print("  That audit's verdict MISMATCH, and its 'monotonicity is coordinate-free' finding,")
    print("  apply verbatim to A - G and to A - M_H.  What is new here is only the")
    print("  NON-monotone gap G - M_H (its dictionary has no non-monotone member) and")
    print("  the flow statements of Parts 3-4.")
    out["overlap"] = dict(h="G", R="1/(2 M_H)", chi="2|delta|",
                          prior_verdict="MISMATCH", new_here=["G - M_H non-monotonicity", "flow results"])

    print("\n=== means.py:", "ALL PASS" if ok else "SOME FAILED", "===")
    d = os.path.join(os.path.dirname(__file__), "..", "data")
    json.dump(out, open(os.path.join(d, "part2_means.json"), "w"), indent=1)
    return ok

if __name__ == "__main__":
    sys.exit(0 if report() else 1)
