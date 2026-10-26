r"""PHASE 4 -- function-field and planted-spectrum control.

Construction fixed in notes/ff_construction.md BEFORE computing:

   phi_i(w) = w^i on the unit circle        <->   phi_{n,a}(z)
   T_R[i,k] = t(|i-k|)                      <->   G_{a,N}
   K_R(w,z) = phi(z)^T conj(T_R^{-1} phi(w)) <->  J_{a,N}(w,z)      (same formula)
   B_R(z;beta) = Phi_R(z) - beta Phi*_R(z)  <->   W(a,theta;z)
   zeros on the unit circle                 <->   zeros real

Phi_R is the monic Szego polynomial (T_{R-1} a = (t(R),...,t(1))^T), Phi*_R its
reversal; B_R is para-orthogonal, and its zeros are the Gauss-quadrature nodes of the
Frobenius measure -- the characteristic polynomial of the rank-one unitary perturbations,
as W(a,theta;.) is the spectral function of the rank-one self-adjoint ones.

The repository's own planting module is reused (src/planted.t_from_blocks).
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "src"))
import mpmath as mp
from planted import t_from_blocks
from suzuki import fp
from weilform import toeplitz

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")


def szego(t, R):
    """Monic Szego polynomial Phi_R and its reversal Phi*_R, as coefficient lists
    (ascending powers).  T_{R-1} a = (t(R), ..., t(1))^T, Phi_R = z^R - sum a_i z^i."""
    T = toeplitz(t, R - 1)                      # (R x R), entries t(|i-j|)
    rhs = mp.matrix([t[R - i] for i in range(R)])
    a = mp.lu_solve(T, rhs)
    coeff = [-a[i] for i in range(R)] + [mp.mpf(1)]
    rev = [mp.conj(c) for c in coeff[::-1]]     # real moments => plain reversal
    return coeff, rev


def popuc_zeros(t, R, beta):
    """Zeros of B_R(z;beta) = Phi_R(z) - beta Phi*_R(z)."""
    c, cs = szego(t, R)
    b = [c[i] - beta * cs[i] for i in range(R + 1)]
    while len(b) > 1 and b[-1] == 0:
        b.pop()
    return mp.polyroots(b[::-1], maxsteps=200, extraprec=200)


def lmin(t, R):
    T = toeplitz(t, R)
    return min(mp.eigsy(T, eigvals_only=True))


def main():
    OK = [True]; out = {}
    def chk(n, c, x=""):
        OK[0] &= bool(c); print(f"  {'PASS' if c else 'FAIL'}  {n} {x}")

    mp.mp.dps = 60
    print("=" * 78); print("PHASE 4 -- function-field / planted-spectrum control")
    print("=" * 78)

    # --------------------------------------------------------------- 4a + 4b
    print("\n-- 4a/4b: genuine spectra (RH here is Weil's theorem) --")
    print("""   Blocks: three ONLINE pairs at angles 0.7, 1.9, 2.8 (so 6 betas, d = 6).
   T_R is nonsingular for R < d; the repository records ker T_R of dim R+1-d for
   R >= d, so the construction runs on R <= 5.""")
    ANG = [mp.mpf('0.7'), mp.mpf('1.9'), mp.mpf('2.8')]
    blocks = [('online', A) for A in ANG]
    t, _geff = t_from_blocks(blocks, 40)
    rows4b = []
    for R in (2, 3, 4, 5, 6):
        lm = lmin(t, R) if R < 6 else mp.nan   # T_6 is singular (ker dim R+1-d = 1)
        zs = popuc_zeros(t, R, mp.mpf(-1))
        offc = max(abs(abs(z) - 1) for z in zs)
        angs = sorted(mp.arg(z) for z in zs)
        # does the zero set interlace / approach the Frobenius angles?
        # mpmath, not float64.  This distance carries the 4b headline ("recovers the
        # Frobenius angles EXACTLY at R = 2g"), and that claim must not rest on a
        # float64 cast rounding a small nonzero distance to 0.0.  See REPORT Retraction 29.
        tgt = sorted([mp.mpf(A) for A in ANG] + [-mp.mpf(A) for A in ANG])
        near = [min(abs(an - s) for s in tgt) for an in angs]
        print(f"   R={R}: lambda_min(T_R) = {mp.nstr(lm,8):>14}, "
              f"max||z|-1| = {mp.nstr(offc,4):>10}, "
              f"max distance to a Frobenius angle = {mp.nstr(max(near),4)}")
        if R < 6:
            chk(f"   R={R}: T_R positive definite", lm > 0)
        chk(f"   R={R}: all POPUC zeros ON the unit circle", offc < mp.mpf('1e-40'),
            f"max||z|-1| = {mp.nstr(offc,3)}")
        rows4b.append(dict(R=R, lmin=fp(lm), max_off_circle=fp(offc),
                           worst_angle_distance=fp(max(near)),
                           zeros_arg=[fp(an) for an in angs]))
    exact6 = [r for r in rows4b if r["R"] == 6]
    # NOT "== 0".  The float64 version of this line printed 0.000e+00 and the claim went
    # out as "EXACTLY"; in mpmath the distance is 7.78e-62 -- the precision floor, not
    # zero.  The tolerance below matches the on-circle test.  See REPORT Retraction 29.
    chk("4b: at R = 2g = 6 the zeros agree with the Frobenius angles to 1e-61",
        bool(exact6) and mp.mpf(exact6[0]["worst_angle_distance"]) < mp.mpf('1e-40'),
        f"max distance = {mp.nstr(mp.mpf(exact6[0]['worst_angle_distance']),4)}"
        if exact6 else "")
    print("""
   The zeros lie on the circle to the precision floor at EVERY R, as the classical
   theorem requires for a positive definite moment matrix.  At R < 2g they do NOT
   reproduce the Frobenius angles -- they are the GAUSS-QUADRATURE NODES of the
   Frobenius measure, which interlace its support.  At R = 2g = 6, the number of atoms,
   the Gauss rule becomes EXACT and the zeros agree with the six angles TO 1e-61
   (max distance 7.779e-62), against 0.34 at R < 2g -- sixty orders of magnitude.
   NOT "exactly", NOT 0.000e+00: that figure was a float64 cast rounding 7.78e-62
   down.  7.78e-62 is the resolution of the computation, not a property of the
   spectrum.  See REPORT Retraction 29.
   (T_5, the 6x6 matrix B_6 is built from, is the LARGEST nonsingular one: ker T_R has
   dimension R+1-d with d = 6.)  So the analogue of Suzuki's W_N recovers the Frobenius
   spectrum to 1e-61, at the one R where it can.""")
    out["P4b"] = rows4b

    # ------------------------------------------------------------------- 4c
    print("\n-- 4c: planted off-circle spectra; two thresholds and the rate --")
    print("""   Blocks: two ONLINE pairs (0.7, 1.9) + one QUARTET at angle 2.8 with
   rho = e^delta.  8 betas, d = 8, so T_R is nonsingular for R <= 7.""")
    print("""   WHICH Toeplitz matrix governs -- an off-by-one, caught and recorded.
   B_R has degree R and is built from T_{R-1} (size R x R) through the Szego normal
   equations.  The first run of this phase compared the zero-departure threshold against
   lambda_min(T_R) (size R+1) and got 28.7% and 187% disagreement.  Against T_{R-1} the
   agreement at R = 3 and 4 is below the bisection resolution.  Both are reported
   below.  See Retraction 26, and Retraction 30 on why that is a BOUND not an identity.""")
    rows4c = []
    for R in (3, 4, 5):
        def tt(d):
            bl = [('online', ANG[0]), ('online', ANG[1]), ('quartet', mp.e ** d, ANG[2])]
            return t_from_blocks(bl, 40)[0]

        # RESOLUTION, not exactness.  Every call uses the SAME bracket and the SAME
        # halving sequence and returns (lo+hi)/2, so two different predicates that flip
        # inside the same final subinterval return BIT-IDENTICAL values.  When that
        # happens it means the two thresholds agree to within BISECT_RES -- it does NOT
        # mean they are equal, and the printed "0.0" relative difference is an artifact
        # of the shared bisection path.  See REPORT Retraction 30.
        BISECT_RES = mp.mpf('1e-35')

        def thresh(f, pred, lo='1e-6', hi='3'):
            lo, hi = mp.mpf(lo), mp.mpf(hi)
            for _ in range(200):
                m = (lo + hi) / 2
                if hi - lo < BISECT_RES:
                    break
                if pred(f(m)):
                    lo = m
                else:
                    hi = m
            return (lo + hi) / 2

        TOL = mp.mpf('1e-30')
        dR = thresh(lambda d: lmin(tt(d), R), lambda v: v > 0)
        dR1 = thresh(lambda d: lmin(tt(d), R - 1), lambda v: v > 0)
        dz = thresh(lambda d: max(abs(abs(z) - 1) for z in popuc_zeros(tt(d), R, mp.mpf(-1))),
                    lambda v: v < TOL)
        rel_R = abs(dz - dR) / dR
        rel_R1 = abs(dz - dR1) / dR1
        # a measured 0 means "below the bisection resolution", so report that bound
        rel_R1_bound = max(rel_R1, BISECT_RES / dR1)
        identical = (dz == dR1)
        print(f"\n   R={R}:")
        print(f"      delta* from lambda_min(T_R)     = {mp.nstr(dR,16)}   rel to zeros: "
              f"{mp.nstr(rel_R,4)}")
        print(f"      delta* from lambda_min(T_(R-1)) = {mp.nstr(dR1,16)}   rel to zeros: "
              f"{mp.nstr(rel_R1,4)}")
        print(f"      delta* where zeros leave circle = {mp.nstr(dz,16)}")
        if identical:
            print(f"      the two bisections returned the SAME value; that bounds the"
                  f" relative gap by {mp.nstr(rel_R1_bound,3)} (= bisection resolution"
                  f" {mp.nstr(BISECT_RES,2)} / delta*), it does NOT show equality")
        agree = rel_R1 < mp.mpf('1e-20')
        chk(f"   R={R}: zero-departure threshold agrees with the lambda_min(T_(R-1))"
            f" threshold to within the bisection resolution", agree,
            f"rel <= {mp.nstr(rel_R1_bound,3)}" + ("" if agree else
            "  <-- the CONVERSE of the classical implication FAILS here"))
        # the rate, measured past the CORRECT threshold
        rates = []
        for k in (1, 2, 5, 10, 20):
            dd = dz * (1 + mp.mpf(k) / 100)
            e = max(abs(abs(z) - 1) for z in popuc_zeros(tt(dd), R, mp.mpf(-1)))
            rates.append((f"+{k}%", fp(e)))
        print(f"      excursion past delta*: " + ", ".join(f"{a} {b}" for a, b in rates))
        try:
            xs = [mp.log(dz * mp.mpf(k) / 100) for k in (1, 2, 5, 10, 20)]
            ys = [mp.log(mp.mpf(b)) for a, b in rates]
            p = (ys[-1] - ys[0]) / (xs[-1] - xs[0])
            pm = (ys[2] - ys[0]) / (xs[2] - xs[0])
            print(f"      fitted exponent max||z|-1| ~ (delta-delta*)^p :  p = {mp.nstr(p,6)}"
                  f"   (on the first three points: {mp.nstr(pm,6)})")
        except Exception:
            p = mp.nan
            print(f"      exponent not fittable (excursion at or below the floor)")
        rows4c.append(dict(R=R, delta_TR=fp(dR), delta_TR1=fp(dR1),
                           delta_zeros=fp(dz), rel_to_TR=fp(rel_R),
                           rel_to_TR1=fp(rel_R1), rel_to_TR1_bound=fp(rel_R1_bound),
                           returned_identical=bool(identical), agree=bool(agree),
                           rates=rates, exponent=fp(p)))
    out["P4c"] = rows4c
    print("""
   READ THIS WITH notes/ff_construction.md Sec.4 AND its correction in REPORT.
   The classical theorem is ONE-DIRECTIONAL: moment matrix positive definite ==> POPUC
   zeros on the circle.  notes/ff_construction.md asserted an EQUIVALENCE; that is an
   overstatement (Retraction 26).  The data above shows the converse failing at R = 5,
   where the zeros stay on the circle well past the loss of positivity.  So:
     - at R = 3 and 4 the two thresholds agree to within the bisection resolution,
       i.e. a relative gap below ~1e-35.  The two bisections return bit-identical
       values, but they share a bracket and a halving sequence, so that is a
       RESOLUTION BOUND and not a proof of equality (Retraction 30);
     - at R = 5 they do NOT, by a factor of about 2.4.
   Either way the preregistered READING -- "the delta^2-L^3 law governs a second
   pipeline" -- is not supported: where the thresholds coincide, the coincidence is the
   classical implication, not the delta^2-L^3 law.""")

    print("\n" + "=" * 78)
    print("PHASE 4 VERDICT:", "see REPORT" if not OK[0] else "PASS (with the 4c amendment)")
    print("=" * 78)
    out["verdict"] = "PASS" if OK[0] else "MIXED"
    json.dump(out, open(os.path.join(DATA, "phase4.json"), "w"), indent=1, default=str)
    return OK[0]


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
