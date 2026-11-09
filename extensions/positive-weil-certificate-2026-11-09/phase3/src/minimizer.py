"""Freeze the even-block minimizing vector and measure the reduction loss on it.

    r_f = R_{T,beta}(f)/||f||^2     certified from the frozen interval matrix
    e_f = q_f - r_f = (1/(2pi||f||^2)) int_{|t|>T} (Psi_L - beta)|F|^2 dt  >= 0

The exploratory eigensolver only CHOOSES the vector; it supplies no certified
eigenvalue.  Once frozen, every quantity below is an Arb enclosure.

e_f is bounded BELOW by integrating over a finite band [T, T2] only: the
integrand is nonnegative there because Psi_L >= beta on |t| >= T (the gate).
An UPPER bound on e_f would need the full tail to infinity; see REPORT for why
that is not cheap for this vector.

Interpretation, kept strict:
  * r_f is a variational upper bound for lambda_R, and for lambda_A in this block.
  * e_f is the reduction loss ON THIS ONE VECTOR, not a gap between global infima.
  * q_f = r_f + e_f is NOT computed here, so this run produces no new upper
    bound for lambda_Q; the existing certified 2.112992612e-17 is retained.
"""
import os, sys, json, time
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "src"))
from flint import arb, acb, ctx
import block as BK

CK = os.path.join(HERE, "..", "..", "phase2", "certs", "block_%s.json")
OUT = os.path.join(HERE, "..", "certs")

def load(sec):
    return [[arb(x) for x in row] for row in json.load(open(CK % sec))]

def minimizer(M, prec, iters=80):
    """EXPLORATORY: inverse iteration on midpoints.  No certified claim."""
    ctx.prec = prec
    n = len(M)
    mid = lambda x: arb(x.mid())
    Mm = [[mid(M[i][j]) for j in range(n)] for i in range(n)]
    R = BK.chol_midpoint(Mm, prec)
    R = [[mid(R[i][j]) for j in range(n)] for i in range(n)]
    def solve(b):
        y = [arb(0)] * n
        for i in range(n):
            acc = b[i]
            for k in range(i):
                acc = acc - R[k][i] * y[k]
            y[i] = mid(acc / R[i][i])
        x = [arb(0)] * n
        for i in range(n - 1, -1, -1):
            acc = y[i]
            for k in range(i + 1, n):
                acc = acc - R[i][k] * x[k]
            x[i] = mid(acc / R[i][i])
        return x
    v = [arb(1)] * n
    for _ in range(iters):
        v = solve(v)
        nr = sum((x * x for x in v), arb(0)).sqrt()
        v = [mid(x / nr) for x in v]
    return v

def freeze(v, digits=40):
    """Exact rationals with denominator 10^digits, scaled so the largest is 1."""
    big = max(range(len(v)), key=lambda k: abs(float(v[k].mid())))
    s = v[big]
    sc = 10 ** digits
    return [Fraction(int(arb(x / s).mid().__float__() * sc) if False else
                     int(round(float((x / s).mid()) * sc)), sc) for x in v]

def Afrac(y):
    return [arb(f.numerator) / arb(f.denominator) for f in y]

def rayleigh(M, y):
    n = len(M)
    num = sum((y[i] * M[i][j] * y[j] for i in range(n) for j in range(n)), arb(0))
    den = sum((y[i] * y[i] for i in range(n)), arb(0))
    return num / den, den

def tail_lower_bound(S, y_even, T, T2, panels, korder, N=220):
    """Certified LOWER bound for (1/2pi) int_{T<|t|<T2} (Psi_L-beta)|F|^2 dt / ||y||^2.

    The integrand is >= 0 on |t| >= T by the gate, so truncating at T2 can only
    lose.  Even sector: the integral over [-T2,-T] equals the one over [T,T2].
    """
    L = S.L
    nodes, wts = BK.gl_rule(korder)
    ev = [n for n in range(N) if n % 2 == 0]
    hw = (arb(T2) - arb(T)) / (2 * panels)
    tot = arb(0)
    for ip in range(panels):
        c = arb(T) + hw * (2 * ip + 1)
        for q in range(korder):
            t = c + hw * nodes[q]
            W = S.W_real(t)
            if not (W >= 0):
                raise RuntimeError("gate violated at t=%s: Psi_L-beta = %s"
                                   % (t.str(10), W.str(10)))
            vv = BK.v_real(L, t, N)
            F = sum((y_even[i] * vv[ev[i]] for i in range(len(ev))), arb(0))
            tot += wts[q] * hw * W * F * F
    return tot / arb.pi()          # (1/2pi) * 2 * int_T^{T2}

def main():
    T, N = 200, 220
    T2, PAN, KO = 600, 400, 16
    PREC = 1400
    t0 = time.time()
    ctx.prec = PREC
    S = BK.Symbol(4, 5, T)
    print("=== frozen even-block minimizer: r_f and the discarded tail ===")
    A = load("even")
    v = minimizer(A, PREC)
    y = freeze(v, 40)
    ctx.prec = PREC
    ya = Afrac(y)
    r_f, nrm2 = rayleigh(A, ya)
    print("  frozen vector: %d coordinates, denominator 10^40" % len(y))
    print("  ||y||^2            = %s" % nrm2.str(16))
    print("  r_f = y^T A y/||y||^2 (CERTIFIED) = %s" % r_f.str(16))
    print("  (A is the frozen interval matrix of R_{T,beta} on the even block)")
    print("  integrating the discarded tail on [%d,%d], %d panels, order %d ..."
          % (T, T2, PAN, KO), flush=True)
    e_lo = tail_lower_bound(S, ya, T, T2, PAN, KO, N) / nrm2
    print("  e_f >= %s     (CERTIFIED lower bound; band [%d,%d] only)"
          % (e_lo.str(16), T, T2))
    ratio = e_lo / r_f
    print("  e_f / r_f >= %s" % ratio.str(10))
    out = dict(T=T, T2=T2, panels=PAN, korder=KO, prec=PREC,
               modes=len(y), denominator="10^40",
               norm_sq=nrm2.str(25), r_f=r_f.str(25),
               e_f_lower_bound=e_lo.str(25), e_over_r_lower=ratio.str(25),
               q_f="NOT COMPUTED - see REPORT phase 4",
               y=[str(f) for f in y])
    os.makedirs(OUT, exist_ok=True)
    json.dump(out, open(os.path.join(OUT, "frozen_minimizer.json"), "w"), indent=2)
    print("  elapsed %.0fs -> certs/frozen_minimizer.json" % (time.time() - t0))
    return 0

if __name__ == "__main__":
    sys.exit(main())
