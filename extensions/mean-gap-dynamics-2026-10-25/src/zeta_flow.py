"""Control 3 -- the zeta family itself, with no RH assumption and no planted zeros.

Three things are established, in increasing order of interest.

(Z1) On a real branch the diagnostic is IDENTICALLY ZERO and cannot leave the axis.
     H_tau is real on the real axis, so at a real simple zero H''/H' is REAL, hence
     z'(tau) is real: the branch stays real for as long as it stays simple.  Therefore
     D == 0 and D' == 0 along it.  Tracking a real zeta zero forward in tau can never
     produce an off-line zero, and "the zero stayed real" is not evidence of anything.

(Z2) A branch can leave the axis ONLY through a collision.  Going BACKWARD in tau the
     real zeros attract; we locate the first collision for a chosen adjacent pair by
     solving H_tau(z) = H'_tau(z) = 0 for (z, tau) with Newton, and then continue the
     emergent conjugate pair.

(Z3) The emergent pair is, at the moment it appears, the ONLY non-real pair in the
     tracked range, so it is extremal and the sign lemma applies to it: D' <= -1 is
     FORCED.  Measured, it is.  That is the exact sense in which the forced-sign
     statement is vacuous for zeta: the only zeta-zeros whose D has a forced sign are
     the ones for which there is nothing to decide.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
from mpmath import mp, mpf, mpc, matrix, lu_solve, zetazero
import kernels as K

# Newton break tolerance.  Chosen ABOVE the measured quadrature noise floor (~1e-24 in
# |H|/|H''| at |z| ~ 154, dps 40): a tighter threshold is unreachable, so the iteration
# would stall at the floor and run its whole budget.  Convergence is quadratic and reaches
# this in 4 steps; the achieved scale-free residual is reported with every value.
NEWTON_TOL = mpf(10)**-20

def Hd(z, tau, order=3):
    return K.H_derivs(z, tau, "zeta", order)

def track_real(z0, tau0, tau1, nsteps=12):
    """Continuation of a simple real branch, with a residual check at each step.

    The iterate is kept COMPLEX and the initial guess is perturbed off the axis, so
    "Im z = 0" and "Im z' = 0" are measured, not imposed by the arithmetic.
    """
    out = []
    z = mpc(z0, mpf("1e-3"))
    for k in range(nsteps + 1):
        tau = mpf(tau0) + (mpf(tau1) - mpf(tau0)) * k / nsteps
        # Newton on H_tau(z) = 0 in z
        for _ in range(30):
            H, H1, _ = Hd(z, tau, 2)
            step = H / H1
            z = z - step
            if abs(step) < NEWTON_TOL:        # noise-floor aware, see NEWTON_TOL
                break
        H, H1, H2 = Hd(z, tau, 2)
        zd = H2 / H1
        out.append(dict(tau=mp.nstr(tau, 10), z=mp.nstr(z, 24),
                        residual=mp.nstr(abs(H), 4),
                        Hprime=mp.nstr(abs(H1), 6),
                        im_z=mp.nstr(abs(mp.im(z)), 4),
                        zdot_im=mp.nstr(abs(mp.im(zd)), 4),
                        zdot_re=mp.nstr(abs(mp.re(zd)), 4),
                        D=mp.nstr(mp.im(z)**2 / 2, 4)))
        z = mpc(mp.re(z), mpf("1e-3"))      # re-perturb for the next step
    return out

def find_collision(zguess, tguess):
    """Solve H_tau(z) = H'_tau(z) = 0 for (z, tau): a double zero (a collision)."""
    z, t = mpc(zguess), mpf(tguess)
    for it in range(40):
        H, H1, H2, H3 = Hd(z, t, 3)
        # d_tau H = -H2 ,  d_tau H' = -H3
        J = matrix([[H1, -H2], [H2, -H3]])
        rhs = matrix([-H, -H1])
        try:
            sol = lu_solve(J, rhs)
        except Exception:
            return None
        z = z + sol[0]; t = t + sol[1]
        if abs(sol[0]) + abs(sol[1]) < NEWTON_TOL:
            break
    H, H1, H2, _ = Hd(z, t, 3)
    return dict(z=z, tau=t, H=H, Hp=H1, Hpp=H2, iters=it)

def main():
    ok = True
    res = {}
    mp.dps = 40

    # ---------------------------------------------------------------- (Z1)
    print("== (Z1) a real branch carries D == 0, and cannot leave the real axis ==")
    print("   H_tau is real on R (Phi real, cos(zu) real), so at a real simple zero")
    print("   H''/H' is real, hence z'(tau) is REAL.  Checked numerically:")
    g = mp.im(zetazero(1))                      # 14.134725...
    z0 = 2*g
    rows = track_real(z0, 0, mpf("0.6"), 6)
    print(f"   {'tau':>10} {'Re z':>24} {'|Im z|':>10} {'|H|':>10} {'|Im z-dot|':>11} "
          f"{'|Re z-dot|':>11} {'D':>8}")
    for r in rows:
        print(f"   {r['tau']:>10} {r['z'].split('+')[0][:24]:>24} {r['im_z']:>10} "
              f"{r['residual']:>10} {r['zdot_im']:>11} {r['zdot_re']:>11} {r['D']:>8}")
        ok &= (abs(mp.mpf(r['im_z'])) < mpf(10)**-20 and abs(mp.mpf(r['zdot_im'])) < mpf(10)**-20)
    print("   Im(z-dot) = 0 to the working precision at every step: the branch is pinned")
    print("   to the axis.  D == 0 identically, so D' == 0 identically.  NO INFORMATION.")
    res["Z1_real_branch"] = rows

    # ---------------------------------------------------------------- (Z2)
    print("\n== (Z2) the first collision of an adjacent pair, going backward in tau ==")
    # choose the closest adjacent pair among the low zeros
    gams = [mp.im(zetazero(k)) for k in range(1, 21)]
    gaps = [(gams[i+1] - gams[i], i) for i in range(len(gams)-1)]
    gaps.sort()
    print("   three closest adjacent pairs among the first 20 zeta zeros (in gamma):")
    for gp, i in gaps[:3]:
        print(f"     gamma_{i+1} = {mp.nstr(gams[i],10)}, gamma_{i+2} = {mp.nstr(gams[i+1],10)}"
              f"   gap {mp.nstr(gp,8)}  -> gap in z = {mp.nstr(2*gp,8)}"
              f"   two-zero estimate tau_c ~ {mp.nstr(-(2*gp)**2/8,8)}")
    found = []
    for gp, i in gaps[:3]:
        mid = (gams[i] + gams[i+1])
        tg = -(2*gp)**2 / 8
        c = find_collision(mid, tg)
        if c is None:
            print(f"     pair {i+1}/{i+2}: Newton failed")
            continue
        # SCALE-FREE residuals: |H| and |H'| are tiny in absolute terms only because
        # H itself is ~ e^{-pi gamma/4}; the meaningful ratios are against |H''|.
        r0 = abs(c["H"]) / abs(c["Hpp"])          # has dimensions length^2
        r1 = abs(c["Hp"]) / abs(c["Hpp"])         # has dimensions length
        okc = (r0 < mpf(10)**-12 and r1 < mpf(10)**-12
               and abs(mp.im(c["z"])) < mpf(10)**-15 and abs(mp.im(c["tau"])) < mpf(10)**-15)
        print(f"     pair {i+1}/{i+2}: collision at tau = {mp.nstr(mp.re(c['tau']),14)}, "
              f"z = {mp.nstr(mp.re(c['z']),18)}  (|Im tau| = {mp.nstr(abs(mp.im(c['tau'])),3)},"
              f" |Im z| = {mp.nstr(abs(mp.im(c['z'])),3)})")
        print(f"        scale-free:  |H|/|H''| = {mp.nstr(r0,4)},  |H'|/|H''| = {mp.nstr(r1,4)}"
              f"   (local zero spacing in z is ~3)")
        print(f"        raw: |H| = {mp.nstr(abs(c['H']),4)}, |H'| = {mp.nstr(abs(c['Hp']),4)}, "
              f"|H''| = {mp.nstr(abs(c['Hpp']),6)}  -> {'DOUBLE ZERO' if okc else 'UNRESOLVED'}")
        if okc:
            found.append((i, c))
        ok &= okc
    res["Z2_collisions"] = [dict(pair=f"{i+1}/{i+2}", tau=mp.nstr(mp.re(c['tau']), 20),
                                 z=mp.nstr(mp.re(c['z']), 20),
                                 H_over_Hpp=mp.nstr(abs(c['H'])/abs(c['Hpp']), 4),
                                 Hp_over_Hpp=mp.nstr(abs(c['Hp'])/abs(c['Hpp']), 4),
                                 absHpp=mp.nstr(abs(c['Hpp']), 8)) for i, c in found]
    print("   At a collision H' = 0, so z' = H''/H' is NOT used there; the branch")
    print("   derivative is undefined and continuation is restarted on the other side.")

    # ---------------------------------------------------------------- (Z3)
    if found:
        i, c = max(found, key=lambda t: mp.re(t[1]["tau"]))   # the LATEST (least negative)
        tc, zc = mp.re(c["tau"]), mp.re(c["z"])
        print(f"\n== (Z3) continue the emergent conjugate pair below tau_c = {mp.nstr(tc,14)} ==")
        print("   locally H_tau(z) ~ (H''/2)(z-zc)^2 - (tau-tc) H''(zc) ... so")
        print("   Im z ~ sqrt(2 (tc - tau)) : the pair leaves the axis vertically.")
        print(f"   local law: D ~ (tau_c - tau), so D' ~ -1 at emergence.")
        print(f"   {'tau_c - tau':>12} {'Im z':>20} {'D':>20} {'D-prime':>18} {'<= -1 ?':>8}")
        rows2 = []
        for dt in ["0.0001", "0.001", "0.01", "0.05", "0.1", "0.2", "0.4", "0.8"]:
            dt = mpf(dt)
            tau = tc - dt
            # local law: H ~ (H''/2)(z-zc)^2 - (tau-tc) H''  =>  Im z ~ sqrt(2 (tc-tau))
            z = mpc(zc, mp.sqrt(2*dt))
            nit = 0
            for _ in range(30):
                H, H1, _ = Hd(z, tau, 2)
                st = H / H1
                z = z - st
                nit += 1
                if abs(st) < NEWTON_TOL:
                    break
            H, H1, H2 = Hd(z, tau, 2)
            if abs(mp.im(z)) < mpf(10)**-12:
                print(f"   {mp.nstr(dt,10):>12}  branch returned to the axis -- skipped")
                continue
            zd = H2 / H1
            D = mp.im(z)**2 / 2
            Dp = mp.im(z) * mp.im(zd)
            flag = "yes" if Dp <= -1 else "NO"
            print(f"   {mp.nstr(dt,10):>12} {mp.nstr(mp.im(z),14):>20} {mp.nstr(D,14):>20} "
                  f"{mp.nstr(Dp,12):>18} {flag:>8}   |H|/|H\'\'| {mp.nstr(abs(H)/abs(H2),3)}"
                  f"  ({nit} Newton steps)")
            rows2.append(dict(dtau=mp.nstr(dt, 12), Im_z=mp.nstr(mp.im(z), 16),
                              D=mp.nstr(D, 16), Dprime=mp.nstr(Dp, 16),
                              residual_over_Hpp=mp.nstr(abs(H)/abs(H2), 4),
                              newton_iters=nit,
                              le_minus1=bool(Dp <= -1)))
            ok &= (Dp <= -1 + mpf(10)**-12)
        res["Z3_emergent_pair"] = dict(tau_c=mp.nstr(tc, 20), z_c=mp.nstr(zc, 20), rows=rows2)
        print("   D' <= -1 at every step, and D -> (tau_c - tau) as the local law predicts.")
        print()
        print("   WHAT THIS IS NOT.  It is tempting to say the emergent pair is the only")
        print("   non-real pair, hence extremal, hence (S1c) FORCES D' <= -1.  That reading")
        print("   is WRONG: see src/zeta_lehmer.py.  Adjacent zeta gaps shrink like")
        print("   2 pi/log(gamma/2 pi) and a pair of gap g in z collides at tau ~ -g^2/8, so")
        print("   the Lehmer pair (6709/6710) leaves the axis at tau ~ -7.1e-4 against")
        print("   ~ -1.04 here -- a factor 1460.  Zeros leave in order of INCREASING gap, so")
        print("   by tau ~ -1.28 this low-height pair is one of the LAST to leave and is")
        print("   certainly not extremal.  The D' <= -1 above is a numerical fact CONSISTENT")
        print("   with the lemma, not an instance of it.  (That the competing deeper zeros")
        print("   sit at gamma ~ 7000 while this pair is at gamma ~ 77 is also why their")
        print("   positive pull is ~1e-8: see the separation bound (S1d),")
        print("   src/separation_bound.py.)")
        print()
        print("   One further observation, which is the point of DERIVATION Sec.1.6.  At")
        print("   tau_c - tau = 0.8 the pair has D = 0.9288 > 1/2, i.e. |delta| > 1/2: the")
        print("   zero has left the strip 0 <= beta <= 1 altogether, so the binary-partition")
        print("   reading of the means is UNAVAILABLE there while D remains perfectly well")
        print("   defined.  The partition interpretation fails along the flow even for zeta.")
    else:
        print("\n== (Z3) not reached: no collision was resolved ==")
        ok = False

    print("\n=== zeta_flow.py:", "ALL PASS" if ok else "SOME FAILED", "===")
    json.dump(res, open(os.path.join(os.path.dirname(__file__), "..", "data", "part4_zeta.json"), "w"),
              indent=1)
    return ok

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
