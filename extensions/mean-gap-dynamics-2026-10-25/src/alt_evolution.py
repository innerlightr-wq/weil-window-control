"""Part 5 -- the one alternative evolution considered, and why it is rejected.

Candidate: the kappa-family.  The Davenport-Heilbronn function is

    F_kappa(s) = [ (1 - i kappa) L(s, chi) + (1 + i kappa) L(s, chi-bar) ] / 2
               = sum_n a_n(kappa) n^{-s},   a = (1, kappa, -kappa, -1, 0)  periodic mod 5,

and kappa is an ARITHMETIC parameter -- tan(arg(eps_chi)/2), derived from the Gauss sum.
It is the only family in this repository whose dependence on the Dirichlet coefficients
can be tracked explicitly while a parameter moves, so it is the one candidate that meets
the brief's stated preference.  This script tests the single question that decides it:

    does  Lambda_{F_kappa}(s) = Lambda_{F_kappa}(1-s)  hold for kappa != kappa_DH ?

If it does not, the family has no functional equation off the DH point, hence no explicit
formula, hence no spectral side, hence its "zeros" are not the zeros of anything carrying
a finite-window Weil form -- and the evolution is not usable.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
from mpmath import mp, mpf, mpc, exp, pi, sqrt
import kernels as K

mp.dps = 50
kDH = K.kappa_closed()

def rootnum_defect(kap):
    """(1 - i kappa) eps_chi - (1 + i kappa);  zero exactly at kappa = kappa_DH."""
    chi = {1: mpc(1), 2: mpc(0, 1), 3: mpc(0, -1), 4: mpc(-1)}
    tau = sum(chi[a]*exp(2j*pi*a/5) for a in range(1, 5))
    eps = tau/(mpc(0, 1)*sqrt(5))
    return (1 - mpc(0, 1)*kap)*eps - (1 + mpc(0, 1)*kap)

def main():
    ok = True
    print("== Part 5: the kappa-family ==")
    print(f"   kappa_DH = {mp.nstr(kDH, 30)}")
    print()
    print(f"   {'kappa':>22} {'|root-number defect|':>22} {'rel |Lam(s)-Lam(1-s)|':>24}")
    rows = []
    s = mpc("0.37", "9.1")
    for kv in [kDH, kDH + mpf("1e-3"), mpf("0.2"), mpf("0.3"), mpf(0), mpf(1)]:
        kap = mpf(kv)
        d = abs(rootnum_defect(kap))
        L1, L2 = K.Lambda_DH(s, kap), K.Lambda_DH(1 - s, kap)
        rel = abs(L1 - L2)/max(abs(L1), mpf(10)**-40)
        rows.append(dict(kappa=mp.nstr(kap, 24), defect=mp.nstr(d, 6), fe_rel=mp.nstr(rel, 6)))
        print(f"   {mp.nstr(kap,18):>22} {mp.nstr(d,6):>22} {mp.nstr(rel,6):>24}")
    # the defect vanishes only at kappa_DH
    ok &= abs(rootnum_defect(kDH)) < mpf(10)**-40
    ok &= abs(rootnum_defect(kDH + mpf("1e-3"))) > mpf(10)**-5
    L1, L2 = K.Lambda_DH(s, mpf("0.2")), K.Lambda_DH(1 - s, mpf("0.2"))
    ok &= abs(L1 - L2)/abs(L1) > mpf("1e-3")

    print("\n   VERDICT on the kappa-family: REJECTED, on a checkable criterion.")
    print("   * what varies:      kappa, a genuine arithmetic parameter (Gauss-sum derived)")
    print("   * what is preserved: the Dirichlet series, its periodicity mod 5, its")
    print("                        abscissa of convergence -- and NOTHING ELSE")
    print("   * functional equation: holds ONLY at kappa = kappa_DH.  Moving kappa by 1e-3")
    print("                        already breaks Lambda(s) = Lambda(1-s) at the 1e-3 level,")
    print("                        and at kappa = 0.2 by O(1).  So off the DH point there is")
    print("                        no completed function, no explicit formula, no spectral")
    print("                        side, and no finite-window Weil form at all.")
    print("   * does zeta occur at a specified parameter?  NO.  zeta is not in the family at")
    print("                        any kappa (conductor 1 vs 5, and zeta has a pole).")
    print("   * are the zeros zeros of a deformed function?  They are zeros of an actual")
    print("                        Dirichlet series, but of one with no functional equation,")
    print("                        so they are not the spectrum of any Weil form.")
    print("   * what inference could return to zeta?  NONE.  There is no parameter path from")
    print("                        the family to zeta, so no limit argument is available.")
    print()
    print("   Three of the brief's five criteria fail outright.  The branch is CLOSED.")
    print()
    print("   No second alternative is pursued.  The repository's own prior findings already")
    print("   cover the obvious ones: extensions/source-flow-detector-audit-2026-10-24 Sec.7")
    print("   records that no appropriate prime-side representation of the deformed family")
    print("   has been established at tau != 0, and its Phase-H experiment (parent")
    print("   extensions/davenport-heilbronn-2026-10-23 Sec.15) records that splitting the")
    print("   arithmetic side by the support of b_n does not separate 'has off-line zeros'")
    print("   from 'has no Euler product'.  Nothing here reopens either.")

    print("\n=== alt_evolution.py:", "ALL PASS" if ok else "SOME FAILED", "===")
    json.dump(dict(kappa_DH=mp.nstr(kDH, 40), rows=rows, verdict="REJECTED"),
              open(os.path.join(os.path.dirname(__file__), "..", "data", "part5_alt.json"), "w"), indent=1)
    return ok

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
