"""Exact rungs delta_{4j} at lambda = e^{4/5} and at lambda^2 = 5, j = 1..4,
   for the eq.(33) test of arXiv:2611.13129v1 against this repository's ladder."""
from mpmath import mp, mpf, sqrt, log, exp, fabs
import sys; sys.path.insert(0, "work")
from prolate_rungs import chis
mp.dps = 90
for tag, lam in (("lambda^2 = 5", sqrt(mpf(5))), ("lambda = e^{4/5}", exp(mpf(4)/5))):
    c = chis(lam, 160, nmax=20)
    print("\n=== %s   (lambda^2 = %s) ===" % (tag, mp.nstr(lam**2, 18)))
    d = {n: 1 - fabs(c[n]) for n in range(0, 18, 2)}
    for n in (4, 6, 8, 12, 16):
        print("   delta_%-2d = %-26s  -log = %s" % (n, mp.nstr(d[n], 14), mp.nstr(-log(d[n]), 11)))
    print("   log(delta_8/delta_4)  = %s" % mp.nstr(log(d[8]/d[4]), 10))
    print("   log(delta_12/delta_4) = %s" % mp.nstr(log(d[12]/d[4]), 10))
    print("   log(delta_16/delta_4) = %s" % mp.nstr(log(d[16]/d[4]), 10))
    print("   log(delta_6/delta_4)  = %s   (odd-vs-even test)" % mp.nstr(log(d[6]/d[4]), 10))
    print("   8 log lambda + 8.2425 = %s   (Fuchs, for comparison)" % mp.nstr(8*log(lam)+mpf('8.2425'), 10))
    if tag.startswith("lambda = e"):
        lam_block = [mpf('1.433072e-17'), mpf('7.227961e-12'), mpf('5.598341e-7'), mpf('0.002863372')]
        print("\n   --- this repository's frozen T=400 even block, lowest 4 ---")
        for j, v in enumerate(lam_block, start=1):
            print("     lam_%d(A) = %-14s   lam_%d/delta_%d = %s"
                  % (j, mp.nstr(v, 7), j, 4*j, mp.nstr(v/d[4*j], 7)))
        print("     log(lam_2/lam_1) = %s   vs exact log(delta_8/delta_4) = %s   shortfall %s"
              % (mp.nstr(log(lam_block[1]/lam_block[0]), 8), mp.nstr(log(d[8]/d[4]), 8),
                 mp.nstr(log(d[8]/d[4]) - log(lam_block[1]/lam_block[0]), 6)))
        print("     (paper at lambda^2=5: log ratio 13.179, exact 13.950, shortfall 0.771;")
        print("      paper's reported shortfall range over all lambda: 0.53 - 0.83)")
