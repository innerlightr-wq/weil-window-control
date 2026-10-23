r"""Planting blocks: the matrix of the contribution of a planted zero configuration.

The explicit formula's zero side is  sum_rho h(z_rho)  with  h(z) = F(z) F(-z)  and
rho = 1/2 + i z (so a critical-line zero has real z, and h(t) = |F(t)|^2 there).  h is even
in z, which is why ordinates come in +- pairs.

Three configurations, all as matrices in the orthonormal basis f_k/||f_k||:

  quartet_block(gamma, delta)  the four zeros 1/2 +- delta +- i gamma, i.e. z = +-gamma +- i delta
                               contributing  4 Re F(gamma + i delta)^2 = 4 (A^2 - B^2)
  pair_block(gamma)            the on-line pair z = +-gamma, contributing 2 F(gamma)^2
  real_block(delta)            the REAL pair 1/2 +- delta, i.e. z = +- i delta,
                               contributing  2 F(i delta) F(-i delta)

The three conventions of notes/proofs.md B.1 are then
  added quartet : + quartet_block
  moved pair    : - pair_block + quartet_block          (zero count preserved)
  moved double  : - 2 pair_block + quartet_block

real_block is the gamma = 0 degeneration and is the Siegel configuration (Phase 5).  Its
delta = 1/2 member is EXACTLY the parent repository's pole-term matrix, which is already
validated end-to-end by the parent's zeta Gate A -- that identity is the regression test in
src/selftest_planting.py and it is what makes the Phase-5 sector sign a checked fact rather
than a bookkeeping assertion.
"""
import os
import sys

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
_PARENT_SRC = os.path.join(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))))), 'src')
if _PARENT_SRC not in sys.path:
    sys.path.insert(0, _PARENT_SRC)

from zeta_window import EVEN, EVEN_D, ODD, F, freq, norm2                # noqa: E402


def F_at_imag(k, L, delta, sector):
    r"""F_k(i delta) = \int_{-L}^{L} f_k(x) e^{-delta x} dx, in closed form.

    even basis: 2 [ w sin(wL) cosh(dL) + d cos(wL) sinh(dL) ] / (w^2 + d^2)   (even in d)
    odd  basis: -2 [ d sin(wL) cosh(dL) - w cos(wL) sinh(dL) ] / (w^2 + d^2)  (odd in d,
                since sin(b_k L) = sin(k pi) = 0 exactly)

    The same expression is used for +-delta, so the parity of F(i delta) in delta is a
    computed consequence here, not an assumption."""
    L, delta = mp.mpf(L), mp.mpf(delta)
    w = freq(k, L, sector)
    den = w ** 2 + delta ** 2
    if sector == ODD:
        return -2 * (delta * mp.sin(w * L) * mp.cosh(delta * L)
                     - w * mp.cos(w * L) * mp.sinh(delta * L)) / den
    return 2 * (w * mp.sin(w * L) * mp.cosh(delta * L)
                + delta * mp.cos(w * L) * mp.sinh(delta * L)) / den


def _nrm(L, idx, sector):
    return [mp.sqrt(norm2(k, mp.mpf(L), sector)) for k in idx]


def real_block(L, idx, delta, sector=EVEN_D):
    r"""2 F(i delta) F(-i delta), symmetrised: the REAL zero pair 1/2 +- delta (gamma = 0).

    Returns (block, Aplus) with Aplus_k = F_k(i delta)/||f_k||, so the caller can report
    the rank-one direction and its norm."""
    L = mp.mpf(L)
    idx = list(idx)
    n = len(idx)
    nrm = _nrm(L, idx, sector)
    Ap = [F_at_imag(k, L, delta, sector) / nrm[i] for i, k in enumerate(idx)]
    Am = [F_at_imag(k, L, -delta, sector) / nrm[i] for i, k in enumerate(idx)]
    B = mp.zeros(n, n)
    for i in range(n):
        for j in range(i, n):
            B[i, j] = B[j, i] = Ap[i] * Am[j] + Ap[j] * Am[i]
    return B, Ap


def h_sign(sector):
    r"""h(z) = F(z) F(-z) computed from the values zeta_window.F returns.

    For the ODD basis zeta_window.F returns \tilde F = F/i, so
    h(z) = F(z)F(-z) = i^2 \tilde F(z) \tilde F(-z) = -\tilde F(z)\tilde F(-z).
    Terms built from F at REAL or purely imaginary z where \tilde F has a definite parity
    absorb this sign automatically; a term built at a GENERAL complex z does not, and must
    carry it explicitly.  (The parent repository's stage3e_certify.move_block is only ever
    called on the even sectors, where h_sign = +1, so this factor is invisible there; it is
    not a correction to any published parent result.)"""
    return -1 if sector == ODD else 1


def pair_block(L, idx, gamma, sector=EVEN_D):
    """2 F_j(gamma) F_k(gamma): the on-line ordinate pair +- gamma."""
    L = mp.mpf(L)
    idx = list(idx)
    n = len(idx)
    nrm = _nrm(L, idx, sector)
    Fr = [F(k, L, mp.mpf(gamma), sector) / nrm[i] for i, k in enumerate(idx)]
    B = mp.zeros(n, n)
    for i in range(n):
        for j in range(i, n):
            B[i, j] = B[j, i] = 2 * Fr[i] * Fr[j]
    return B


def quartet_block(L, idx, gamma, delta, sector=EVEN_D):
    """4 Re F(gamma + i delta)^2, symmetrised: the off-line quartet 1/2 +- delta +- i gamma."""
    L = mp.mpf(L)
    idx = list(idx)
    n = len(idx)
    nrm = _nrm(L, idx, sector)
    w = mp.mpc(mp.mpf(gamma), mp.mpf(delta))
    Fw = [F(k, L, w, sector) / nrm[i] for i, k in enumerate(idx)]
    Fm = [F(k, L, -w, sector) / nrm[i] for i, k in enumerate(idx)]
    sg = h_sign(sector)
    B = mp.zeros(n, n)
    for i in range(n):
        for j in range(i, n):
            B[i, j] = B[j, i] = 2 * sg * mp.re(Fw[i] * Fm[j] + Fw[j] * Fm[i])
    return B


def K_closed(L, g):
    r"""K(L,gamma) = int_{-L}^{L} u^2 sin^2(gamma u) du, the exact extremal of Theorem 8.
    Character-independent.  K(L, 0) = 0."""
    L, g = mp.mpf(L), mp.mpf(g)
    if g == 0:
        return mp.mpf(0)
    return (L ** 3 / 3
            - (L ** 2 * mp.sin(2 * g * L) / (2 * g)
               + L * mp.cos(2 * g * L) / (2 * g ** 2)
               - mp.sin(2 * g * L) / (4 * g ** 3)))
