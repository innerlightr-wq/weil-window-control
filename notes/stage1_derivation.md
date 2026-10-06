# Stage 1 derivation (T1 = exact / proved)

## Setup
`C/F_q` smooth projective curve, genus `g`, `N_n = #C(F_{q^n})`.
Weil zeta function:

    Z_C(T) = exp( sum_{n>=1} N_n T^n / n ) = P(T) / ((1-T)(1-qT)),
    P in Z[T], deg P = 2g, P(0) = 1,  P(T) = prod_{j=1}^{2g} (1 - a_j T).

## The point-count identity
Take `log` of both sides:

    sum_n N_n T^n/n = log P(T) - log(1-T) - log(1-qT)
                    = -sum_j sum_n a_j^n T^n/n + sum_n T^n/n + sum_n q^n T^n/n.

Matching coefficients of `T^n/n`:

    (1)   N_n = 1 + q^n - p_n ,     p_n := sum_j a_j^n .

Weil's theorem (RH for curves, a THEOREM here -- this is why Stage 1 is a control):
`|a_j| = sqrt(q)`. Write `a_j = sqrt(q) e^{i th_j}`. `P` has real coefficients, so the
multiset `{a_j}` is closed under conjugation, hence `{th_j}` is closed under negation.
Therefore

    p_n = q^{n/2} sum_j e^{i n th_j} = q^{n/2} sum_j cos(n th_j)   (imaginary parts cancel).

Define the window data

    (2)   t(n) := sum_j cos(n th_j) = p_n q^{-n/2}
               = q^{n/2} + q^{-n/2} - N_n q^{-n/2}      (n >= 1),
          t(0) = 2g.

This matches the identity in the project brief. VERIFIED numerically in
`src/stage1_curves.py` (`--verify`) against brute-force point counts.

## Weil positivity is automatic (T1)
For real `c = (c_0..c_R)` put `chat(th) = sum_i c_i e^{i i th}` and
`T_R[i][j] = t(|i-j|)`. Because `{th_j}` is negation-symmetric,
`t(i-j) = sum_m e^{i (i-j) th_m}`, so

    (3)   Q(c) = c^T T_R c = sum_m (sum_i c_i e^{i i th_m})(sum_j c_j e^{-i j th_m})
                           = sum_{m=1}^{2g} |chat(th_m)|^2  >= 0.

Equivalently `T_R = sum_m w_m w_m^H` with `(w_m)_i = e^{i i th_m}`, so

    (4)   rank T_R <= min(R+1, #distinct th_m) <= min(R+1, 2g),   T_R PSD.

So the geometric-side Toeplitz form is PSD *because* `|a_j| = sqrt(q)`.
Planting `|a| != sqrt(q)` (Stage 2) is exactly what can break (3).

## Why H1 should be free of RH (Caratheodory-Fejer / Pisarenko)
`T_R - lam_min I` is still Toeplitz (diagonal shift), is PSD by construction, and is
singular. If `lam_min` is simple then it has rank exactly `R` at size `R+1`.
Caratheodory-Fejer: a PSD Toeplitz matrix of size `n+1` with rank `n` is the moment
matrix of a *unique* measure with `n` atoms on the circle; its kernel vector's
polynomial has all `n` zeros simple and on `|z| = 1`. This is Pisarenko harmonic
decomposition. Nothing about `|a_j| = sqrt(q)` enters. H1 is therefore predicted TRUE
even for planted off-circle data -- the "witness" is free.

Note also `J T_R J = T_R^T = T_R` for the flip `J`, so `T_R` commutes with `J`;
a simple eigenvalue forces its eigenvector to be exactly even or exactly odd.

## H3 (only N_1..N_g are needed)
Functional equation of `Z_C`: `a_j -> q/a_j` permutes the roots, giving for
`P(T) = sum_i A_i T^i`:

    (5)   A_{2g-i} = q^{g-i} A_i ,   A_0 = 1.

Newton's identities convert `p_1..p_g` <-> `A_1..A_g`, and `p_n = 1 + q^n - N_n`.
So `N_1..N_g` determine all of `P`, hence all `t(n)`. Tested exactly in code.
