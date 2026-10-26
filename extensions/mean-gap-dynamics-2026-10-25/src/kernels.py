"""Kernels and the coordinate dictionary for the de Bruijn--Newman heat family.

Normalisation is the repository's (extensions/source-flow-detector-audit-2026-10-24,
Polymath15 convention), reproduced here from scratch and re-verified, not inherited:

    H_tau(z) = int_0^infty e^{tau u^2} Phi(u) cos(z u) du,     d_tau H = - d_z^2 H
    s        = 1/2 + i z / 2,        so   z = -2i (s - 1/2)
    for s = beta + i gamma :   Re z = 2 gamma ,  Im z = -2 (beta - 1/2) = -2 delta

Hence the displacement diagnostic of this investigation,

    delta = beta - 1/2 = -Im(z)/2 ,     D = 2 delta^2 = (Im z)^2 / 2 .

zeta:   Phi(u)   = sum_{n>=1} (2 pi^2 n^4 e^{9u} - 3 pi n^2 e^{5u}) exp(-pi n^2 e^{4u})
        H_0(z)   = xi(1/2 + i z/2) / 8
DH:     Phi_F(x) = 4 e^{3x} sum_{n>=1} n a_n exp(-pi n^2 e^{4x} / 5)
        H_0(z)   = Lambda_F(1/2 + i z/2) / 2

Nothing in this module is imported from elsewhere in the repository.  Agreement with
the committed values of the parent extensions is treated as a *check*, never as input.
"""
from mpmath import mp, mpf, mpc, exp, cos, pi, sqrt, gamma, zeta, quad, log

# ---------------------------------------------------------------- zeta kernel

def Phi_zeta(u, nmax=None):
    """The de Bruijn--Newman kernel for zeta (Polymath15 Phi)."""
    u = mp.mpf(u) if not isinstance(u, mp.mpf) else u
    if nmax is None:
        # exp(-pi n^2 e^{4u}) below the working epsilon
        e4 = exp(4 * u)
        nmax = max(8, int(mp.sqrt(mp.mpf(mp.dps) * 2.4 / float(e4)) + 4))
    tot = mp.mpf(0)
    e4u, e9u, e5u = exp(4 * u), exp(9 * u), exp(5 * u)
    for n in range(1, nmax + 1):
        n2 = mp.mpf(n) ** 2
        term = (2 * pi ** 2 * n2 ** 2 * e9u - 3 * pi * n2 * e5u) * exp(-pi * n2 * e4u)
        tot += term
        if n > 4 and abs(term) < mp.eps * (abs(tot) + 1) / 100:
            break
    return tot


def xi(s):
    """The completed Riemann zeta, xi(s) = (1/2) s (s-1) pi^{-s/2} Gamma(s/2) zeta(s)."""
    return s * (s - 1) / 2 * pi ** (-s / 2) * gamma(s / 2) * zeta(s)


# ------------------------------------------------------------------ DH kernel

def kappa_closed():
    """kappa = [sqrt(10 - 2 sqrt 5) - 2] / (sqrt 5 - 1)."""
    r5 = sqrt(5)
    return (sqrt(10 - 2 * r5) - 2) / (r5 - 1)


def kappa_gauss():
    """kappa derived from the Gauss sum of the odd primitive character mod 5.

    chi(1,2,3,4) = (1, i, -i, -1) is the odd primitive character mod 5
    (chi(2) = i since 2 is a primitive root mod 5 and chi has order 4).
    tau(chi) = sum_a chi(a) e(a/5);  eps_chi = tau(chi) / (i sqrt 5)  for odd chi;
    kappa = tan( arg(eps_chi) / 2 ).
    """
    chi = {1: mpc(1), 2: mpc(0, 1), 3: mpc(0, -1), 4: mpc(-1)}
    tau = sum(chi[a] * exp(2j * pi * a / 5) for a in range(1, 5))
    eps = tau / (mpc(0, 1) * sqrt(5))
    return mp.tan(mp.arg(eps) / 2)


def dh_coeffs(kap):
    """a_n for n = 1..5, periodic mod 5:  (1, kappa, -kappa, -1, 0)."""
    return [mp.mpf(1), kap, -kap, mp.mpf(-1), mp.mpf(0)]


def Phi_DH(x, kap=None, nmax=None):
    x = mp.mpf(x) if not isinstance(x, mp.mpf) else x
    if kap is None:
        kap = kappa_closed()
    a = dh_coeffs(kap)
    e4 = exp(4 * x)
    if nmax is None:
        nmax = max(16, int(mp.sqrt(mp.mpf(mp.dps) * 3.7 * 5 / float(e4)) + 6))
    tot = mp.mpf(0)
    for n in range(1, nmax + 1):
        an = a[(n - 1) % 5]
        if an == 0:
            continue
        tot += n * an * exp(-pi * mp.mpf(n) ** 2 * e4 / 5)
    return 4 * exp(3 * x) * tot


def F_DH(s, kap=None):
    """F(s) = sum a_n n^{-s} = 5^{-s} sum_{r=1}^{4} a_r zeta(s, r/5)."""
    if kap is None:
        kap = kappa_closed()
    a = dh_coeffs(kap)
    return 5 ** (-s) * sum(a[r - 1] * zeta(s, mp.mpf(r) / 5) for r in range(1, 5))


def Lambda_DH(s, kap=None):
    """Lambda_F(s) = (5/pi)^{(s+1)/2} Gamma((s+1)/2) F(s), entire, = Lambda_F(1-s)."""
    return (mp.mpf(5) / pi) ** ((s + 1) / 2) * gamma((s + 1) / 2) * F_DH(s, kap)


# ------------------------------------------------------------- the heat family

def _nodes(z, xmax):
    """Quadrature nodes resolving cos(zu): about two nodes per oscillation.

    Without this the integrand at z ~ 171 oscillates ~170 times per unit and
    mp.quad on a coarse split silently returns only ~10 correct digits -- which
    is how the first version of the PDE check at z = 171.4 appeared to FAIL.
    """
    xmax = mp.mpf(xmax)
    n = max(8, int(float(abs(z)) * float(xmax) / 2.0) + 1)
    return [xmax * k / n for k in range(n + 1)]


def H_tau(z, tau, kernel="zeta", xmax=3):
    """H_tau(z) = int_0^infty e^{tau u^2} Phi(u) cos(z u) du  by adaptive quadrature.

    The kernel decays double-exponentially, so the omitted tail beyond xmax is
    astronomically small; the quadrature error inside is a *library estimate*,
    not a bound.  (Same limitation as the parent audit.)
    """
    Phi = Phi_zeta if kernel == "zeta" else Phi_DH
    f = lambda u: exp(tau * u ** 2) * Phi(u) * cos(z * u)
    return quad(f, _nodes(z, xmax))


def H_derivs(z, tau, kernel="zeta", order=2, xmax=3):
    """(H, H\', H\'\') in z, by differentiating cos(zu) under the integral sign.

    d_z^k cos(zu) = u^k * {cos, -sin, -cos, sin}[k mod 4], so each derivative is a
    quadrature of the same shape -- no finite differences anywhere.
    """
    Phi = Phi_zeta if kernel == "zeta" else Phi_DH
    nd = _nodes(z, xmax)
    out = []
    for k in range(order + 1):
        trig = [cos, lambda w: -mp.sin(w), lambda w: -cos(w), mp.sin][k % 4]
        f = lambda u, k=k, trig=trig: exp(tau * u ** 2) * Phi(u) * u ** k * trig(z * u)
        out.append(quad(f, nd))
    return out


# ------------------------------------------------------------- the diagnostic

def delta_of_z(z):
    """delta = beta - 1/2 = -Im(z)/2."""
    return -mp.im(z) / 2


def D_of_z(z):
    """D = 2 delta^2 = (Im z)^2 / 2 = A - M_H."""
    return mp.im(z) ** 2 / 2


def zdot(z, tau, kernel="zeta", xmax=3):
    """z'(tau) = H''(z)/H'(z) on a simple-zero branch."""
    _, H1, H2 = H_derivs(z, tau, kernel, 2, xmax)
    return H2 / H1


def Ddot(z, tau, kernel="zeta", xmax=3):
    """D'(tau) = Im(z) * Im( H''(z)/H'(z) )."""
    return mp.im(z) * mp.im(zdot(z, tau, kernel, xmax))
