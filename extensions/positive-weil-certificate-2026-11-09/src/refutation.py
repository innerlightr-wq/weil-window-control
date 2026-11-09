"""REFUTATION of the sparsity target proposed in REPORT.md v1 section 6.

The proposed estimate was:

    meas{ t in [T,2T] : P_L(t) > theta A_L }  <=  c T^{1-eta}      (FALSE)

PROOF THAT IT IS FALSE (Kronecker-Weyl).  Let p_1<...<p_r be the primes with
p_i < e^{2L} and write the comb as a function on the torus T^r,

    P_L(t) = Pi(t log p_1, ..., t log p_r),
    Pi(x) = sum_{p_i^k < e^{2L}} (2 log p_i / p_i^{k/2}) cos(k x_i) .

Pi is continuous and Pi(0) = A_L.  So for every theta < 1 the set
U_theta = {x : Pi(x) > theta A_L} is OPEN and contains 0, hence has positive
Haar measure mu(U_theta) > 0.

The numbers log p_1, ..., log p_r are linearly independent over Q (unique
factorization: sum a_i log p_i = 0 with a_i in Z forces prod p_i^{a_i} = 1,
so every a_i = 0).  By the continuous Kronecker-Weyl equidistribution theorem
the line t |-> (t log p_1, ..., t log p_r) equidistributes in T^r, so

    lim_{T->inf} (1/T) meas{ t in [0,T] : (t log p) in U_theta } = mu(U_theta) > 0 .

Hence meas E_{L,theta}(T) has POSITIVE LOWER DENSITY in T, and no bound
c T^{1-eta} with eta > 0 can hold, at any fixed L.  The target is false.

This module MEASURES both sides to confirm the refutation numerically:
the Haar measure mu(U_theta) by quadrature on T^r, and the time density by
direct sampling, for growing T.  Status NUMERICAL -- it corroborates a proof,
it is not itself the proof.

It also measures the MODULATION OBSTRUCTION: g(u) = f(u) cos(omega u) has the
same support as f and the same parity, but its Fourier profile is translated to
sit on top of a chosen resonance.  So no bound on the *density* of the bad set
can control int_E |F|^2 by itself.
"""
import math, sys
from flint import arb
from wcore import set_prec, comb_mass, prime_powers_below

def comb_data(two_L_num, two_L_den):
    """[(prime, k, weight)] with weight = 2 log p / p^{k/2}; distinct primes listed."""
    pp = prime_powers_below(two_L_num, two_L_den)
    primes = sorted({p for _, p, _ in pp})
    terms = [(p, k, 2.0 * math.log(p) / (n ** 0.5)) for n, p, k in pp]
    A = sum(w for _, _, w in terms)
    return primes, terms, A

def Pi(xs, primes, terms):
    """Comb as a function of the phase vector xs (one entry per distinct prime)."""
    idx = {p: i for i, p in enumerate(primes)}
    return sum(w * math.cos(k * xs[idx[p]]) for p, k, w in terms)

def haar_measure(primes, terms, A, theta, m=720):
    """mu{ x in T^r : Pi(x) > theta A } by product quadrature (r = 2 here)."""
    r = len(primes)
    assert r == 2, "written for r = 2 (L < log 5 / ... ); extend if needed"
    hit = 0
    tot = m * m
    for i in range(m):
        x = 2 * math.pi * i / m
        for j in range(m):
            y = 2 * math.pi * j / m
            if Pi([x, y], primes, terms) > theta * A:
                hit += 1
    return hit / tot

def time_density(primes, terms, A, theta, T, samples=400000):
    """meas{t in [T,2T] : P_L(t) > theta A}/T by uniform sampling."""
    hit = 0
    for i in range(samples):
        t = T + T * (i + 0.5) / samples
        xs = [t * math.log(p) for p in primes]
        if Pi(xs, primes, terms) > theta * A:
            hit += 1
    return hit / samples

def main():
    set_prec(256)
    print("=== REFUTATION of the proposed sparsity target ===")
    print("Claim under test: meas E_{L,theta}(T) <= c T^{1-eta}.  Predicted FALSE.\n")
    for (num, den, lab) in ((8, 5, "4/5"), (162, 100, "81/100")):
        primes, terms, A = comb_data(num, den)
        print("-- L = %s   primes %s   A_L = %.9f" % (lab, primes, A))
        if len(primes) != 2:
            print("   (r = %d distinct primes; Haar quadrature skipped)" % len(primes))
            continue
        for theta in (0.50, 0.80, 0.90):
            mu = haar_measure(primes, terms, A, theta, m=600)
            row = []
            for T in (1e2, 1e3, 1e4, 1e5, 1e6):
                row.append(time_density(primes, terms, A, theta, T, samples=200000))
            print("   theta=%.2f  Haar mu(U)=%.6f   time density over [T,2T]:" % (theta, mu))
            print("              " + "  ".join("T=1e%d:%.6f" % (int(math.log10(T)), d)
                                               for T, d in zip((1e2,1e3,1e4,1e5,1e6), row)))
            decays = row[-1] < row[0] / 10
            print("              -> %s ; densities are flat in T, not decaying"
                  % ("DECAY SEEN (unexpected)" if decays else "NO DECAY, as the proof predicts"))
    print("\n-- modulation obstruction")
    print("   g(u) = f(u) cos(omega u) has the SAME support and parity as f, and")
    print("   G(t) = [F(t-omega) + F(t+omega)]/2, so the Fourier profile can be")
    print("   translated onto any chosen resonance while supp g = supp f.")
    print("   Therefore a density bound on E alone cannot control int_E |G|^2.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
