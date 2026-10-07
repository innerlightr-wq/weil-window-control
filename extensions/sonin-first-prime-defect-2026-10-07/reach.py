#!/usr/bin/env python3
r"""Section 8: the REACH of the theta_S multiplier.

prod_{p in S\{inf}} |L_p(1/2+is)|^{-2} = prod_p (1-p^{-1/2-is})(1-p^{-1/2+is})
is a Laurent polynomial of degree <= 1 in each p^{is}.  So the dilates of f it can
generate are exactly  d/d'  with d,d' coprime SQUAREFREE divisors of rad(S).
The Weil prime block needs the atoms h(n^{+-1}) for every prime power n = p^k with
log n < 2L.

*** AMENDED 2026-10-07 (after review).  The COMPUTATIONS below are unchanged and
*** correct: they concern the BARE multiplier J*J only.  The inference that higher
*** prime powers are therefore "out of reach" for the full projection, and the
*** resulting ceiling L < log 2, are RETRACTED: Pi_S contains K^{-1}, and
*** 1/M(s) = 2 + 2sqrt2 cos(s log2) + 2 cos(2 s log2) + ... has INFINITE Laurent
*** support.  See amend_checks.py and FIRST_PRIME_DEFECT.md section 3.2.  The words
*** "unreachable", "out of reach" and "CEILING" in the output below are retained as
*** a record of what was run and must be read as RETRACTED.
"""
import math, itertools
from fractions import Fraction as Fr

def vm(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return (p, math.log(p)) if m == 1 else None
    return None

print("="*78); print("8.  Reach of the multiplier: squarefree shifts only"); print("="*78)
print("  *** AMENDED 2026-10-07: the expansion below is correct for the BARE multiplier")
print("  *** J*J.  The inference to the FULL projection Pi_S (which contains K^-1) is")
print("  *** RETRACTED -- 1/M has infinite Laurent support.  See amend_checks.py.")

# 8a. exact expansion of prod |1-p^{-1/2-is}|^2 as a Laurent polynomial in p^{is}
def expand(primes):
    """returns {(e_p,...): coefficient} for prod_p (1 - p^-1/2 x_p)(1 - p^-1/2 /x_p)."""
    terms = {tuple(0 for _ in primes): Fr(1)}
    for i, p in enumerate(primes):
        new = {}
        # (1 - p^-1/2 x)(1 - p^-1/2 x^-1) = (1 + 1/p) - p^-1/2 x - p^-1/2 x^-1
        # keep p^-1/2 symbolic by tracking the exponent of p^-1/2 separately:
        # coefficients are (rational part, half-power of p).  Use floats only to print.
        for key, c in terms.items():
            for sh, lab in ((0, 'const'), (1, 'up'), (-1, 'dn')):
                k = list(key); k[i] = sh; k = tuple(k)
                if lab == 'const':   co = c * (1 + Fr(1, p))
                else:                co = -c                      # times p^{-1/2}
                new[k] = new.get(k, Fr(0)) + co
        terms = new
    return terms

for primes in [(2,), (2, 3), (2, 3, 5)]:
    t = expand(primes)
    shifts = sorted(t.keys())
    print("\n  S\\{inf} = %s   -> %d monomials, exponents in {-1,0,1} per prime"
          % (list(primes), len(shifts)))
    print("     max |exponent| per prime:",
          [max(abs(k[i]) for k in shifts) for i in range(len(primes))])
    dil = sorted({Fr(1,1) * Fr(int(__import__('math').prod(p**max(k[i],0) for i,p in enumerate(primes))),
                               int(__import__('math').prod(p**max(-k[i],0) for i,p in enumerate(primes))))
                  for k in shifts})
    print("     dilates generated (d/d'):", [str(x) for x in dil])
    nonsf = [str(x) for x in dil if x.numerator > 1 and
             any(x.numerator % (q*q) == 0 for q in primes)]
    print("     any non-squarefree numerator?", nonsf if nonsf else "NONE")

# 8b. which active prime powers are reachable, as a function of L
print("\n  active prime powers vs reachability (reachable <=> squarefree)")
print("   L        2L       active n (log n < 2L)      unreachable (k>=2)")
for L in [0.36, 0.45, 0.5493061, 0.60, 0.6931472, 0.70, 0.75, 0.80, 1.00]:
    act = [n for n in range(2, 60) if vm(n) and math.log(n) < 2*L]
    bad = [n for n in act if any(n % (q*q) == 0 for q in range(2, n+1))]
    print("   %-8.6f %-8.6f %-24s %s" % (L, 2*L, act, bad if bad else "none"))
print("\n  first non-squarefree prime power is n = 4 ; it switches on at 2L > log 4")
print("   CEILING  L < log 2 =", '%.7f' % math.log(2),
      "  (vs the archimedean CC2021 window L = (log2)/2 =", '%.7f' % (math.log(2)/2), ")")
print("   Zhu's certified window L <= 0.8 has 0.8 >", '%.7f' % math.log(2),
      "-> n = 4 IS active there and IS out of reach:",
      4 in [n for n in range(2,60) if vm(n) and math.log(n) < 1.6])

# 8c. the k=1 atom weight that W_2 needs on the calibration window
print("\n  on the calibration window (log2)/2 < L < (log3)/2 the prime-2 block has")
print("   ONLY the k=1 atoms, weight log2 * 2^(-1/2) = %.7f" % (math.log(2)*2**-0.5))
print("   (k=2 would need 2L > log 4 = %.7f, but 2L < log 3 = %.7f)"
      % (math.log(4), math.log(3)))
print("   -> on THIS window the shift supports COINCIDE; the separation below is")
print("      a statement about weights and sign, not about support.")

# 8d. Neumann radius for K^-1 around the mean-field value
print("\n" + "="*78); print("9.  Neumann expansion of K^-1 around (3/2)P"); print("="*78)
r = 2*math.sqrt(2)/3
print("  K = (3/2)P - sqrt2 * (PCP) = (3/2)[ P - (2 sqrt2 / 3) PCP ]")
print("  ratio  (2 sqrt2)/3 = %.7f  < 1 ;  ||PCP|| <= ||C|| = 1" % r)
print("  so the Neumann series K^-1 = (2/3) sum_n ((2 sqrt2/3) PCP)^n CONVERGES")
print("  UNCONDITIONALLY in operator norm, for every S-independent bound ||PCP|| <= 1.")
print("  worst-case tail after n terms (||PCP||=1):")
for n in (0, 1, 2, 4, 8, 16, 32, 64):
    print("     n = %2d   remaining factor r^(n+1)/(1-r) = %.6e" % (n, r**(n+1)/(1-r)))
print("  CONTROL: a multiplier with a real zero would break this.")
for a in (2**-0.5, 0.9, 1.0):
    lo = 1 - a; hi = 1 + a; rr = 2*a/(1+a*a)
    print("     |1 - a e^{-i t}|^2 : min = %.6f  ratio 2a/(1+a^2) = %.6f  %s"
          % (lo*lo, rr, "ok" if rr < 1 else "DIVERGES (a=1 gives a zero at t=0)"))
print("\n" + "="*78)
print("RETRACTION NOTICE (2026-10-07)")
print("="*78)
print("  Section 8b/8c above printed 'unreachable', 'out of reach' and 'CEILING'.")
print("  Those WORDS ARE RETRACTED.  What is proved is only that the bare multiplier")
print("  J*J = prod_p |L_p(1/2+is)|^-2 has Laurent degree <= 1 in each p^(is).")
print("  Pi_S = J P K^-1 P J* contains the INVERSE, and")
print("     1/M(s) = 2 + 2sqrt2 cos(s log2) + 2 cos(2 s log2) + ...")
print("  carries the k=2 ('4') harmonic with coefficient 2 and every higher one.")
print("  So no general prime-power obstruction is proved, and equally Pi_S is NOT")
print("  shown to produce the correct Weil coefficients.  Run amend_checks.py.")
print("\nDONE")
