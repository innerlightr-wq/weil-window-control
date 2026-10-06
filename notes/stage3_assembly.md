# Stage 3 assembly (derivation; verify each step numerically)

## The form
f real, even, supported in [-L, L];  F(t) = \int f(x) e^{itx} dx.  Zhu's geometric side:

    Q(f) = 2 F(i/2)^2 + (1/2pi) \int_R |F(t)|^2 Psi_L(t) dt
    Psi_L(t) = Re psi(1/4 + it/2) - log pi - sum_{log n < 2L} 2 Lambda(n) n^{-1/2} cos(t log n)

## Basis (closed-form transforms)
    phi_k(x) = cos(a_k x) 1_{|x|<=L},   a_k = (2k+1) pi / (2L),   k = 0,1,2,...

Even, and phi_k(+-L) = cos((2k+1)pi/2) = 0, so phi_k is CONTINUOUS on R (needed below).
Orthogonal:  \int_{-L}^{L} phi_j phi_k = L delta_{jk}.  Orthonormal basis: phi_k / sqrt(L).

    F_k(t) = L [ sinc((a_k - t)L) + sinc((a_k + t)L) ]        (sinc y = sin y / y)
           = (-1)^k  2 a_k cos(tL) / (a_k^2 - t^2)
    F_k(i/2) = (-1)^k 2 a_k cosh(L/2) / (a_k^2 + 1/4)

(The sinc form is the numerically stable one: the second has removable 0/0 at t = +-a_k.)

## Everything moves to u-space by Parseval
With F_j F_k = (phi_j * phi_k)^ , Fourier inversion gives, for G := phi_j * phi_k
(supported in [-2L, 2L], even, C^1 because the phi are continuous):

    (1/2pi) \int F_j F_k dt            = G(0) = L delta_{jk}
    (1/2pi) \int F_j F_k cos(tu) dt    = G(u)

so the prime term is a FINITE SUM and the -log pi term is diagonal. No oscillatory
quadrature anywhere.

### Closed form for G
For u in [0, 2L]  (G even, so this is all of it):
    G_{jk}(u) = (1/2) [ I(a_j - a_k,  a_k u) + I(a_j + a_k, -a_k u) ],
    I(p, q)   = \int_{u-L}^{L} cos(p x + q) dx
              = [ sin(pL + q) - sin(p(u-L) + q) ] / p      (p != 0)
              = (2L - u) cos(q)                            (p = 0)
with a_j - a_k = (j-k) pi / L  and  a_j + a_k = (j+k+1) pi / L.

### Archimedean term
From psi(z) = -gamma + sum_{m>=0} [ 1/(m+1) - 1/(m+z) ] at z = 1/4 + it/2:

    Re psi(1/4 + it/2) = -gamma + sum_{m>=0} [ 1/(m+1) - 1/(c_m + it) - 1/(c_m - it) ],
    c_m = 2m + 1/2.

Since 1/(c+it) + 1/(c-it) has inverse Fourier transform e^{-c|u|}, the archimedean
distribution is  K(u) = (-gamma + sum 1/(m+1)) delta(u) - sum_m e^{-c_m |u|}, and pairing
with G (even) gives the CONVERGENT series

    Arch_{jk} = -gamma G(0) + sum_{m>=0} [ G(0)/(m+1) - 2 \int_0^{2L} G(u) e^{-c_m u} du ].

Each bracket is O(1/m^2) with a full asymptotic expansion in 1/m (G is C^1 with
G'(0) = 0), so Richardson/Euler-Maclaurin acceleration applies.

## Assembled matrix (orthonormal basis phi_k/sqrt L)
    Qtilde_{jk} = (1/L) [ 2 F_j(i/2) F_k(i/2) + Arch_{jk} - log(pi) L delta_{jk}
                          - sum_{n: log n < 2L} 2 Lambda(n) n^{-1/2} G_{jk}(log n) ]

and lambda_min(Qtilde) = min over ||f||_{L^2} = 1 of Q(f).

## Independent gate (do this BEFORE trusting any eigenvalue)
Under RH the explicit formula says  Q(f) = sum_rho |F(gamma_rho)|^2.  Compute the right
side directly from mpmath.zetazero for many zeros and compare to the assembled left side
for random f in the span. This validates the assembly without reference to Zhu.

## Both parity sectors

The pole term is really `h(i/2) + h(-i/2) = 2 F(i/2) F(-i/2)`.  For real f,
`F(-i/2) = +F(i/2)` if f is EVEN and `-F(i/2)` if f is ODD, so

    pole term = +2 F_j(i/2) F_k(i/2)   (even sector, the brief's formula)
              = -2 P_j P_k             (odd sector),  P_k := F_k(i/2) real.

Odd basis (continuous, vanishing at +-L):  psi_k(x) = sin(b_k x), b_k = k pi / L, k>=1.

    F_k(t) = i L [ sinc((b_k - t)L) - sinc((b_k + t)L) ]        (purely imaginary on R)
    P_k    = F_k(i/2) = 2 b_k (-1)^k sinh(L/2) / (b_k^2 + 1/4)
    \int psi_j psi_k = L delta_{jk}

Parseval, stated once for both sectors with the CROSS-CORRELATION
`C_{jk}(u) = \int f_j(x) f_k(x-u) dx`  (which equals the convolution when both are even):

    F_j(t) conj(F_k(t)) = Chat_{jk}(t)
    (1/2pi) \int F_j conj(F_k) dt           = C_{jk}(0) = L delta_{jk}
    (1/2pi) \int F_j conj(F_k) cos(tu) dt   = (C_{jk}(u) + C_{jk}(-u))/2
    archimedean term                        = pairing of K with the symmetrised C.

For the odd basis, for u in [0, 2L]:
    C_{jk}(u) = (1/2) [ I(b_j - b_k, b_k u) - I(b_j + b_k, -b_k u) ],   same I as before.
Note C is NOT even in u when j != k (C_{jk}(-u) = C_{kj}(u)); the matrix uses the
symmetrisation, which is what the form sees.

NOTE on the odd sector and RH.  Weil positivity concerns the full space; the odd sector
is part of it.  The CvS hypothesis is about a simple lowest eigenvalue with an EVEN
eigenfunction, so the even/odd gap `lambda_min^odd - lambda_min^even` is the quantity
that decides whether the global ground state is even.
