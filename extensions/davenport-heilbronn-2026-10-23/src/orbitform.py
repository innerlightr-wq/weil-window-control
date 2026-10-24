r"""The exact rank-<=4 matrix of one zero orbit, in the window basis.

For an off-line quartet O = {beta +- i gamma, 1-beta +- i gamma} (DERIVATION.md (4.1)),

    Q_O(g) = 4 ( |P|^2 - |S|^2 ),
    P = int g cosh(delta u) e^{i gamma u} du,   S = int g sinh(delta u) e^{i gamma u} du,

and with A = G(delta + i gamma), B = G(-delta + i gamma) we have P = (A+B)/2,
S = (A-B)/2.  Writing p_m, s_m for the values of P, S on the normalised basis function
f_{idx[m]}/sqrt(norm2), the form is EXACTLY

    Q_O = 4 ( (Re p . c)^2 + (Im p . c)^2 ) - 4 ( (Re s . c)^2 + (Im s . c)^2 ),

a signed rank-<= 4 matrix

    M_O = 4 (pr pr^T + pi pi^T) - 4 (sr sr^T + si si^T).

In a pure-parity sector two of the four vectors vanish identically and the rank is 2.
This is what Phase E subtracts and Phase F removes; it is never approximated.
"""
import mpmath as mp

import dh


def orbit_vectors(idx, L, sector, beta, gamma):
    """(pr, pi, sr, si) as lists of mpf, one entry per basis index."""
    d = mp.mpf(beta) - mp.mpf(1) / 2
    g = mp.mpf(gamma)
    pr, pi_, sr, si = [], [], [], []
    for k in idx:
        e = [mp.mpf(0)] * len(idx)
        e[idx.index(k)] = mp.mpf(1)
        A = dh.G_at(e, idx, L, sector, mp.mpc(d, g))
        B = dh.G_at(e, idx, L, sector, mp.mpc(-d, g))
        P, S = (A + B) / 2, (A - B) / 2
        pr.append(mp.re(P)); pi_.append(mp.im(P))
        sr.append(mp.re(S)); si.append(mp.im(S))
    return pr, pi_, sr, si


def orbit_matrix(idx, L, sector, beta, gamma, mult=1):
    pr, pi_, sr, si = orbit_vectors(idx, L, sector, beta, gamma)
    n = len(idx)
    M = mp.zeros(n, n)
    for i in range(n):
        for j in range(n):
            M[i, j] = mult * 4 * (pr[i] * pr[j] + pi_[i] * pi_[j]
                                  - sr[i] * sr[j] - si[i] * si[j])
    return M


def online_matrix(idx, L, sector, gamma, mult=1):
    """2 H(i gamma) as a rank-<=2 matrix: 2 (Re G . c)^2 + 2 (Im G . c)^2 times the
    sector sign carried by H(i gamma) = G(i gamma) G(-i gamma).  Computed by polarisation
    from dh.online_pair_contribution so no sector special-casing is needed."""
    n = len(idx)
    M = mp.zeros(n, n)
    base = []
    for i in range(n):
        e = [mp.mpf(0)] * n
        e[i] = mp.mpf(1)
        base.append(e)
    for i in range(n):
        M[i, i] = mult * dh.online_pair_contribution(base[i], idx, L, sector, gamma)
    for i in range(n):
        for j in range(i + 1, n):
            e = [mp.mpf(0)] * n
            e[i] = e[j] = mp.mpf(1)
            tot = mult * dh.online_pair_contribution(e, idx, L, sector, gamma)
            M[i, j] = M[j, i] = (tot - M[i, i] - M[j, j]) / 2
    return M


def check_orbit_matrix(idx, L, sector, beta, gamma, coeffs):
    """Polarisation check: c^T M_O c must equal dh.quartet_contribution(c)."""
    M = orbit_matrix(idx, L, sector, beta, gamma)
    n = len(idx)
    q = mp.fsum([coeffs[i] * coeffs[j] * M[i, j] for i in range(n) for j in range(n)])
    return q, dh.quartet_contribution(coeffs, idx, L, sector, beta, gamma)
