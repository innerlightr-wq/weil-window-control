"""Phase B: the one-stroke reduction's SCALAR GATE, reproduced and sharpened.

Imported result (Zhu, arXiv:2608.24827v2, Thm 1.1): for real even f with
supp f in [-L,L],

    Q(f) >= 2F(i/2)^2 + (1/pi) int_0^{T#} (Psi_L - beta*) |F|^2 dt + beta* ||f||_2^2 ,

valid whenever Psi_L(t) >= beta* for every t >= T#, where

    Psi_L(t) = Re psi(1/4 + i t/2) - log pi - P_L(t),
    P_L(t)   = sum_{log n < 2L} (2 Lambda(n)/sqrt n) cos(t log n),
    A_L      = sum_{log n < 2L} 2 Lambda(n)/sqrt n      (Lemma 3.2: sup_t P_L = A_L).

The paper takes beta* = log(T#/2pi) - 1/T# - A_L, from a lower bound on the digamma.

PROVED HERE (G1-G3), replacing that lower bound by an exact monotonicity fact:

  G1. Re psi(1/4 + i t/2) = -gamma + sum_{n>=0} [ 1/(n+1) - (n+1/4)/((n+1/4)^2 + t^2/4) ].
  G2. Each summand is strictly increasing in t > 0, since
         d/dt [ -(n+1/4)/((n+1/4)^2 + t^2/4) ] = (n+1/4)(t/2)/((n+1/4)^2+t^2/4)^2 > 0.
      Hence t |-> Re psi(1/4 + i t/2) is strictly increasing on t > 0.
  G3. Therefore, with  beta_sharp(T#) := Re psi(1/4 + i T#/2) - log pi - A_L,
         Psi_L(t) >= beta_sharp(T#)   for every t >= T# ,
      because the digamma part is >= its value at T# and P_L(t) <= A_L.

beta_sharp >= beta_paper always (the paper's expression is a lower bound for the
same quantity), so the gate is reached no later.  Define the EXACT GATE

      T_gate(L) := the unique t > 0 with  Re psi(1/4 + i t/2) = log pi + A_L ,

the least admissible T# for which beta > 0 at all.  It is located here by
certified bisection on Arb balls.
"""
from flint import arb, acb, ctx
from wcore import set_prec, comb_mass, ball_str

CHECKS = []
def chk(name, ok, extra=""):
    CHECKS.append((name, bool(ok), extra))
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name, ("  -- " + extra) if extra else ""))
    return ok

def re_digamma_half(t):
    """Re psi(1/4 + i t/2) as an Arb ball (t an arb ball or number)."""
    z = acb(arb(1) / 4, arb(t) / 2)
    return z.digamma().real

def beta_sharp(T, A):
    return re_digamma_half(T) - arb.pi().log() - A

def beta_paper(T, A):
    T = arb(T)
    return (T / (2 * arb.pi())).log() - 1 / T - A

def T1_paper(A):
    return 2 * arb.pi() * A.exp()

def root_of(fn, A, lo=arb(1), hi=arb(10) ** 6, iters=200):
    """Certified bisection for the unique sign change of fn(., A) on [lo,hi].

    Both beta_sharp and beta_paper are strictly increasing in T (G2 for the first;
    d/dT[log(T/2pi) - 1/T] = 1/T + 1/T^2 > 0 for the second), so the bracket is valid.
    """
    lo, hi = arb(lo), arb(hi)
    if not (fn(lo, A) < 0):
        raise RuntimeError("lower bracket is not below the root")
    if not (fn(hi, A) > 0):
        raise RuntimeError("upper bracket is not above the root")
    for _ in range(iters):
        mid = (lo + hi) / 2
        v = fn(mid, A)
        if v > 0:
            hi = mid
        elif v < 0:
            lo = mid
        else:
            break
        if (hi - lo) < arb(10) ** (-40):
            break
    return lo, hi

def gate_threshold(A, lo=arb(1), hi=arb(10) ** 6, iters=200):
    """Certified bisection for the unique root of Re psi(1/4+it/2) = log pi + A.

    Uses G2 (strict monotonicity), so the sign of beta_sharp brackets the root.
    Returns a ball containing T_gate.
    """
    lo, hi = arb(lo), arb(hi)
    if not (beta_sharp(lo, A) < 0):
        raise RuntimeError("lower bracket is not below the gate")
    if not (beta_sharp(hi, A) > 0):
        raise RuntimeError("upper bracket is not above the gate")
    for _ in range(iters):
        mid = (lo + hi) / 2
        b = beta_sharp(mid, A)
        if b > 0:
            hi = mid
        elif b < 0:
            lo = mid
        else:
            break
        if (hi - lo) < arb(10) ** (-40):
            break
    return lo, hi

def report_for(L_num, L_den, Tsharp_list, label):
    A, rows = comb_mass(2 * L_num, L_den)
    print("\n--- L = %d/%d   (2L = %d/%d) ---" % (L_num, L_den, 2 * L_num, L_den))
    print("    prime powers with log n < 2L: %s" % [n for n, _, _, _ in rows])
    print("    A_L  = %s" % ball_str(A))
    T1 = T1_paper(A)
    print("    T_1  = 2 pi exp(A_L) = %s   (paper's sufficient threshold)" % ball_str(T1))
    lo, hi = root_of(beta_sharp, A)
    plo, phi = root_of(beta_paper, A)
    print("    T_gate   (root of beta_sharp, certified) = %s" % hi.str(18))
    print("    T_paper  (root of beta_paper, certified) = %s" % phi.str(18))
    print("    T_1 = 2 pi e^{A_L} is BELOW both: it is the root of log(T/2pi) - A_L,")
    print("        which drops the -1/T term, so it is not itself an admissible T#.")
    chk("%s: T_1 < T_gate < T_paper (sharp gate admits a smaller T#)" % label,
        bool(T1 < lo) and bool(hi < plo),
        "T_1=%s  T_gate=%s  T_paper=%s" % (T1.str(12), hi.str(12), phi.str(12)))
    chk("%s: beta_paper(T_1) = -1/T_1 < 0, so T_1 itself is inadmissible" % label,
        bool(beta_paper(T1, A) < 0), "beta_paper(T_1) = %s" % ball_str(beta_paper(T1, A), 6))
    for T in Tsharp_list:
        bs, bp = beta_sharp(T, A), beta_paper(T, A)
        print("    T# = %-6s  beta_paper = %-28s beta_sharp = %s"
              % (T, ball_str(bp, 10), ball_str(bs, 10)))
        if T > 0:
            chk("%s: beta_sharp >= beta_paper at T#=%s" % (label, T), bool(bs >= bp),
                "gain = %s" % ball_str(bs - bp, 6))
    return A, T1, (lo, hi), rows

def main():
    set_prec(400)
    print("=== PHASE B: the one-stroke reduction, scalar gate ===")
    print("Imported: Zhu arXiv:2608.24827v2 Thm 1.1 (reduction) and Lemma 3.2 (sup P_L = A_L).")
    print("Proved here: G1-G3 monotonicity, giving beta_sharp >= beta_paper.")

    A8, T1_8, g8, rows8 = report_for(4, 5, [200], "L=4/5")
    # paper's reported values, as a reproduction check
    chk("L=4/5: A_L reproduces the paper's 2.9420", bool((A8 - arb("2.9420")).abs_lower() < arb("0.0001")),
        ball_str(A8, 8))
    bp = beta_paper(200, A8)
    chk("L=4/5: beta_paper(T#=200) reproduces the paper's 0.5135",
        bool((bp - arb("0.5135")).abs_lower() < arb("0.0002")), ball_str(bp, 8))
    chk("L=4/5: paper's T# = 200 exceeds T_1", bool(arb(200) > T1_8),
        "T_1 = %s" % T1_8.str(10))

    A81, T1_81, g81, rows81 = report_for(81, 100, [200, 600], "L=81/100")
    chk("L=81/100: prime 5 activates (n=5 present, absent at L=4/5)",
        (5 in [n for n, _, _, _ in rows81]) and (5 not in [n for n, _, _, _ in rows8]))
    jump = (A81 - A8)
    exact_jump = 2 * arb(5).log() / arb(5).sqrt()
    chk("L=81/100: A_L jump equals exactly 2 log 5 / sqrt 5", bool((jump - exact_jump).abs_upper() < arb(10) ** -30),
        "jump = %s" % ball_str(jump, 12))
    ratio = T1_81 / T1_8
    chk("L=81/100: T_1 multiplies by exp(2 log 5/sqrt 5) at the activation",
        bool((ratio - exact_jump.exp()).abs_upper() < arb(10) ** -20),
        "T_1 ratio = %s" % ball_str(ratio, 12))
    print("\n    T_1(81/100)/T_1(4/5) = %s   (exact factor exp(2 log 5 / sqrt 5))"
          % ball_str(ratio, 12))
    bad = [n for n, ok, _ in CHECKS if not ok]
    print("\nPHASE B: %d checks, %d failures" % (len(CHECKS), len(bad)))
    return 1 if bad else 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
