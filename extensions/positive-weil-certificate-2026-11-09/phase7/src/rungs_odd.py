"""Odd-sector rungs delta_{4j+2} at lambda = e^{4/5}, against this repository's
   frozen T=400 odd block.  Hypothesis 3.4 pairs eps^odd(j) with delta_{4j+2}."""
from mpmath import mp, mpf, sqrt, log, exp, fabs
import sys; sys.path.insert(0, "work")
from prolate_rungs import chis
mp.dps = 90
for tag, lam in (("lambda^2 = 5", sqrt(mpf(5))), ("lambda = e^{4/5}", exp(mpf(4)/5))):
    c = chis(lam, 160, nmax=20)
    d = {n: 1 - fabs(c[n]) for n in range(0, 18, 2)}
    print("\n=== %s ===" % tag)
    for n in (6, 10, 14):
        print("   delta_%-2d = %-24s -log = %s" % (n, mp.nstr(d[n], 14), mp.nstr(-log(d[n]), 11)))
    print("   log(delta_10/delta_6) = %s" % mp.nstr(log(d[10]/d[6]), 10))
    if tag.startswith("lambda = e"):
        blk = [mpf('1.38574029573e-14'), mpf('2.42537508217e-9'), mpf('6.88937205273e-5')]
        print("\n   --- frozen T=400 ODD block, lowest 3 ---")
        for j, v in enumerate(blk, start=1):
            print("     lam_%d(A_odd) = %-16s  lam_%d/delta_%d = %s"
                  % (j, mp.nstr(v, 7), j, 4*j+2, mp.nstr(v/d[4*j+2], 7)))
        print("     log(lam_2/lam_1) = %s  vs exact log(delta_10/delta_6) = %s  shortfall %s"
              % (mp.nstr(log(blk[1]/blk[0]), 8), mp.nstr(log(d[10]/d[6]), 8),
                 mp.nstr(log(d[10]/d[6]) - log(blk[1]/blk[0]), 6)))
        print("     paper, lambda^2=5: odd j=1 ratio 4.84 (range 2.92-4.84);")
        print("       log(eps^odd(1)/eps(1)) shortfall vs log(delta_6/delta_4): 0.29-0.39")
