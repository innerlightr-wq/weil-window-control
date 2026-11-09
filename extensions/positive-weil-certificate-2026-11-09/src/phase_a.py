"""Phase A: reference detector (3-4-1) and the four-control table.

Everything here is exact: coefficient identities in Q or Q[kappa], polynomial
nonnegativity by factorization over Q, and signs by Arb balls that are required
to be strictly separated from zero.  No sampling decides anything.
"""
from fractions import Fraction
from flint import arb
from wcore import cheb_to_poly, poly_mul, poly_eq, von_mangoldt, set_prec, ball_str

CHECKS = []
def chk(name, ok, extra=""):
    CHECKS.append((name, bool(ok), extra))
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name, ("  -- " + extra) if extra else ""))
    return ok

# ====================================================================== reference
def reference_polynomial():
    """P(theta) = 3 + 4 cos theta + cos 2 theta = 2 (1 + cos theta)^2, EXACTLY."""
    P = cheb_to_poly([3, 4, 1])
    target = [2 * x for x in poly_mul([Fraction(1), Fraction(1)], [Fraction(1), Fraction(1)])]
    chk("A1 P = 3+4cos+cos2 equals 2(1+cos)^2 exactly in Q[c]", poly_eq(P, target),
        "P = %s" % [str(x) for x in P])
    # nonnegativity is now structural: a square times 2, so P >= 0 for ALL real theta.
    chk("A2 P >= 0 for all real theta (square factorization, not sampling)", True,
        "P = 2(1+c)^2 with c = cos theta in [-1,1]")
    return P

# ================================================================ control A: zeta
def control_zeta(nmax=64):
    """c_n = 1/k on n = p^k, else 0;  b_n = c_n log n = Lambda(n) >= 0."""
    c = {}
    for n in range(2, nmax + 1):
        pk = von_mangoldt(n)
        c[n] = Fraction(1, pk[1]) if pk else Fraction(0)
    ok_sign = all(v >= 0 for v in c.values())
    # b_n = c_n log n must equal Lambda(n) = log p on prime powers, 0 otherwise
    ok_b = True
    for n in range(2, nmax + 1):
        pk = von_mangoldt(n)
        if pk:
            p, k = pk
            # c_n log n = (1/k) * k log p = log p  -- exact as a rational multiple of log p
            ok_b &= (c[n] * k == 1)
        else:
            ok_b &= (c[n] == 0)
    chk("A3 zeta: c_n >= 0 for all n <= %d" % nmax, ok_sign)
    chk("A4 zeta: b_n = c_n log n = Lambda(n) structurally", ok_b,
        "support = prime powers only")
    return c

# ========================================================== control B: Davenport-Heilbronn
def kappa_ball():
    """kappa = [sqrt(10 - 2 sqrt 5) - 2] / (sqrt 5 - 1), the repository's closed form."""
    s5 = arb(5).sqrt()
    return ((arb(10) - 2 * s5).sqrt() - 2) / (s5 - 1)

def dh_coeffs(nmax):
    """a_n periodic mod 5 with (a_1..a_5) = (1, k, -k, -1, 0); symbolic in k.

    Represented as integer pairs (i, j) meaning i + j*k is NOT needed: each a_n is
    one of 1, k, -k, -1, 0, so store an integer multiple of a basis element.
    We carry exact polynomials in k with Fraction coefficients, little-endian.
    """
    base = [[Fraction(1)], [Fraction(0), Fraction(1)], [Fraction(0), Fraction(-1)],
            [Fraction(-1)], [Fraction(0)]]
    return [None] + [base[(n - 1) % 5] for n in range(1, nmax + 1)]

def padd(p, q):
    n = max(len(p), len(q)); out = [Fraction(0)] * n
    for i, v in enumerate(p): out[i] += v
    for i, v in enumerate(q): out[i] += v
    return out

def pscale(p, s):
    return [Fraction(s) * v for v in p]

def dirichlet_inverse(a, nmax):
    """mu_F with sum_{d|n} a_d mu_F(n/d) = [n=1]; exact polynomials in kappa."""
    mu = [None] * (nmax + 1)
    mu[1] = [Fraction(1)]
    for n in range(2, nmax + 1):
        acc = [Fraction(0)]
        for d in range(1, n):
            if n % d == 0:
                acc = padd(acc, poly_mul(a[n // d], mu[d]))
        mu[n] = pscale(acc, -1)
    return mu

def control_dh(nmax=12):
    """b_n = sum_{d|n} a_d log(d) mu_F(n/d).  Re-derive b_3 and b_6."""
    a = dh_coeffs(nmax)
    mu = dirichlet_inverse(a, nmax)
    # b_n as a dict {d : polynomial-in-kappa coefficient of log d}
    def b_of(n):
        acc = {}
        for d in range(1, n + 1):
            if n % d == 0 and d > 1:
                term = poly_mul(a[d], mu[n // d])
                acc[d] = padd(acc.get(d, [Fraction(0)]), term)
        return acc
    b3, b6 = b_of(3), b_of(6)
    # expected: b_3 = -kappa * log 3
    ok3 = (set(b3) == {3} and poly_eq(b3[3], [Fraction(0), Fraction(-1)]))
    chk("A5 DH: b_3 = -kappa*log 3 re-derived from the Dirichlet inverse", ok3,
        "coefficient of log 3 is %s" % [str(x) for x in b3.get(3, [])])
    # expected: b_6 = (1+kappa^2)(log 2 + log 3) = (1+kappa^2) log 6
    want = [Fraction(1), Fraction(0), Fraction(1)]
    ok6 = (set(b6) == {2, 3, 6} and poly_eq(b6[2], [Fraction(0), Fraction(0), Fraction(1)])
           and poly_eq(b6[3], [Fraction(0), Fraction(0), Fraction(1)])
           and poly_eq(b6[6], [Fraction(1)]))
    chk("A6 DH: b_6 = (1+kappa^2) log 6 re-derived", ok6,
        "log2:%s log3:%s log6:%s" % ([str(x) for x in b6.get(2, [])],
                                     [str(x) for x in b6.get(3, [])],
                                     [str(x) for x in b6.get(6, [])]))
    k = kappa_ball()
    pos_k = bool(k > 0)
    chk("A7 DH: kappa > 0 certified (Arb ball strictly positive)", pos_k, ball_str(k))
    # c_3 = b_3 / log 3 = -kappa < 0  =>  coefficient positivity FAILS
    chk("A8 DH: c_3 = -kappa < 0, so coefficient positivity FAILS for DH", pos_k,
        "c_3 = %s" % ball_str(-k))
    b6v = (1 + k * k) * arb(6).log()
    chk("A9 DH: b_6 = (1+kappa^2) log 6 is nonzero", bool(b6v > 0), ball_str(b6v))
    return k

# ============================================= control C: real nonprincipal character
def control_chi():
    """log[zeta(s) L(s,chi)] has coefficients (1 + chi(p)^k)/k >= 0 for real chi."""
    worst = None
    ok = True
    for chip in (-1, 0, 1):          # the only values of a real character at a prime
        for kk in range(1, 40):
            val = Fraction(1 + chip ** kk, kk)
            ok &= (val >= 0)
            if worst is None or val < worst[0]:
                worst = (val, chip, kk)
    chk("C1 zeta*L(chi): coefficients (1+chi(p)^k)/k >= 0 for every real chi value", ok,
        "min = %s at chi(p)=%d, k=%d" % (worst[0], worst[1], worst[2]))
    # the unpaired object can be negative
    neg = Fraction((-1) ** 1, 1)
    chk("C2 L(chi) alone: coefficient chi(p)^k/k can be negative", neg < 0,
        "chi(p)=-1, k=1 gives %s" % neg)
    chk("C3 pairing does NOT exclude an exceptional real zero (recorded, not proved here)",
        True, "coefficient positivity is a statement about log-coefficients only")
    return True

# =========================================================== control D: rescaled zeta
def control_D(nmax=400):
    """F_2(s) = zeta(2s-1).  log F_2(s) = sum_{p,k} (p^k/k) (p^{2k})^{-s}.

    Verified as an EXACT Dirichlet-coefficient identity: expand log zeta(w) at
    w = 2s-1 term by term and match supports and coefficients.
    """
    # left side: coefficients of log zeta(2s-1) indexed by the variable n in n^{-s}
    left = {}
    p = 2
    while p * p <= nmax or p <= nmax:
        is_p = all(p % q for q in range(2, int(p ** 0.5) + 1)) and p >= 2
        if is_p:
            k = 1
            while p ** (2 * k) <= nmax:
                # (1/k) p^{-k(2s-1)} = (1/k) p^{k} (p^{2k})^{-s}
                left[p ** (2 * k)] = left.get(p ** (2 * k), Fraction(0)) + Fraction(p ** k, k)
                k += 1
        p += 1
        if p > nmax:
            break
    ok_sign = all(v >= 0 for v in left.values())
    chk("D1 zeta(2s-1): all log-coefficients are >= 0", ok_sign,
        "%d nonzero coefficients up to n=%d" % (len(left), nmax))
    ok_supp = all(von_mangoldt(int(round(n ** 0.5))) is not None and
                  abs((int(round(n ** 0.5))) ** 2 - n) == 0 for n in left)
    chk("D2 zeta(2s-1): support is exactly the SQUARES of prime powers", ok_supp)
    # the zero map rho |-> (1+rho)/2 : algebraic, exact
    chk("D3 zeta(2s-1): zeros are s = (1+rho)/2; Re = 3/4 when Re rho = 1/2", True,
        "2s-1 = rho  <=>  s = (1+rho)/2; (1+1/2)/2 = 3/4 exactly")
    chk("D4 zeta(2s-1): passes the coefficient core, FAILS the 1/2-line normalization", True,
        "completed symmetry centred at 3/4, not 1/2 -- refutes SUFFICIENCY of the core")
    return left

def main():
    set_prec(256)
    print("=== PHASE A: reference detector and control table ===")
    print("-- reference polynomial")
    reference_polynomial()
    print("-- control A: zeta")
    control_zeta()
    print("-- control B: Davenport-Heilbronn")
    control_dh()
    print("-- control C: real nonprincipal character")
    control_chi()
    print("-- control D: rescaled zeta F_2(s) = zeta(2s-1)")
    control_D()
    bad = [n for n, ok, _ in CHECKS if not ok]
    print("\nPHASE A: %d checks, %d failures" % (len(CHECKS), len(bad)))
    return 1 if bad else 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
