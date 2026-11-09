"""Phase 5, part 3: the 2x2 table.  R_200 and R_400 evaluated on f_200 and f_400,
with each function's original Q value from the compact time-domain route."""
import os, sys, json, time
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "src"))
sys.path.insert(0, HERE)
from flint import arb, ctx
import block as BK, compact as C
from eval_f200 import load_f200

PREC = 1100
C2 = os.path.join(HERE, "..", "..", "phase2", "certs")
C4 = os.path.join(HERE, "..", "certs")

def load_block(path):
    return [[arb(x) for x in row] for row in json.load(open(path))]

def rayleigh(M, y):
    n = len(M)
    num = sum((y[i] * M[i][j] * y[j] for i in range(n) for j in range(n)), arb(0))
    den = sum((y[i] * y[i] for i in range(n)), arb(0))
    return num / den

def minimizer_of(M, prec, iters=80):
    import math
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
        nr = sum((z * z for z in v), arb(0)).sqrt()
        v = [mid(z / nr) for z in v]
    return v

def main():
    t0 = time.time()
    ctx.prec = PREC
    L = arb(4) / 5
    S = BK.Symbol(4, 5, 200)
    out = {}

    f200, fm = load_f200(L)
    y200 = [arb(Fraction(s).numerator) / arb(Fraction(s).denominator) for s in fm["y"]]
    A200 = load_block(os.path.join(C2, "block_even.json"))
    A400p = os.path.join(C4, "block400_even.json")
    have400 = os.path.exists(A400p)
    print("=== 2x2 table: R_200, R_400 on f_200, f_400 ===")
    r200_f200 = rayleigh(A200, y200)
    print("  R_200(f_200) = %s  [CERTIFIED]" % r200_f200.str(20))
    out["R200_f200"] = r200_f200.str(25)
    q200 = arb(json.load(open(os.path.join(C4, "q_f200.json")))["q"])
    print("  Q(f_200)     = %s  [CERTIFIED, compact route]" % q200.str(20))
    out["Q_f200"] = q200.str(25)

    if not have400:
        print("  (T=400 block not yet available: %s)" % A400p)
        json.dump(out, open(os.path.join(C4, "cross.json"), "w"), indent=2)
        return 0

    A400 = load_block(A400p)
    n4 = len(A400)
    y200_pad = y200 + [arb(0)] * (n4 - len(y200))
    r400_f200 = rayleigh(A400, y200_pad)
    print("  R_400(f_200) = %s  [CERTIFIED]" % r400_f200.str(20))
    out["R400_f200"] = r400_f200.str(25)
    e400_f200 = q200 - r400_f200
    print("  e_400(f_200) = Q - R_400 = %s  [CERTIFIED]" % e400_f200.str(20))
    out["e400_f200"] = e400_f200.str(25)

    # f_400: freeze the new block's minimizer
    v = minimizer_of(A400, PREC)
    ctx.prec = PREC
    big = max(range(n4), key=lambda k: abs(float(v[k].mid())))
    sc = 10 ** 40
    y400 = [Fraction(int(round(float((z / v[big]).mid()) * sc)), sc) for z in v]
    ya = [arb(f.numerator) / arb(f.denominator) for f in y400]
    r400_f400 = rayleigh(A400, ya)
    print("  R_400(f_400) = %s  [CERTIFIED]" % r400_f400.str(20))
    out["R400_f400"] = r400_f400.str(25)
    out["f400_y"] = [str(f) for f in y400]

    degs = [2 * i for i in range(n4)]
    x400 = [arb(BK.pole_sigma(degs[i]) * y400[i].numerator) / arb(y400[i].denominator)
            for i in range(n4)]
    f400 = C.Poly(L, degs, x400)
    print("  evaluating Q(f_400) by the compact route (deg <= %d) ..." % f400.maxdeg, flush=True)
    cache = {}
    def g_of(x):
        k = arb(x).str(40)
        if k not in cache:
            cache[k] = C.autocorr(f400, x)
        return cache[k]
    Q4, meta = C.Q_compact(f400, +1, S.terms, g_of=g_of)
    q400 = Q4 / f400.norm2()
    print("  Q(f_400)     = %s  [CERTIFIED, compact route]" % q400.str(20))
    out["Q_f400"] = q400.str(25)
    out["e400_f400"] = (q400 - r400_f400).str(25)
    print("  e_400(f_400) = %s  [CERTIFIED]" % (q400 - r400_f400).str(20))
    out["elapsed_s"] = round(time.time() - t0, 1)
    json.dump(out, open(os.path.join(C4, "cross.json"), "w"), indent=2)
    print("  elapsed %.0fs -> certs/cross.json" % (time.time() - t0))
    return 0

if __name__ == "__main__":
    sys.exit(main())
