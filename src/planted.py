"""Stage 2: planted (fake) Frobenius configurations.

A functional-equation-closed, conjugation-closed multiset of "Frobenius eigenvalues"
is built from blocks. Writing  beta_j = alpha_j / sqrt(q)  (normalised Frobenius),
the functional equation alpha -> q/alpha becomes beta -> 1/beta, and RH-for-curves
becomes |beta| = 1.  The window data is

    t(n) = sum_j beta_j^n ,    t(0) = #{beta_j} = 2 g_eff .

NOTE (T1): t(n) depends ONLY on the normalised beta_j, i.e. only on (rho, phi).
q drops out entirely. So Stage 2 needs no choice of q.

Blocks:
  ONLINE pair   (rho = 1, angle phi):   {e^{+i phi}, e^{-i phi}}
        -> t contribution  2 cos(n phi).           Toeplitz block: rank 2, PSD.
  OFFLINE quartet (rho != 1, angle phi): {rho e^{+-i phi}, rho^{-1} e^{+-i phi}}
        = {beta, 1/beta, conj(beta), 1/conj(beta)} with beta = rho e^{i phi}
        -> t contribution  2 cos(n phi) (rho^n + rho^{-n}).
        Toeplitz block: rank 4, signature (2, ., 2)  -- INDEFINITE.
  OFFLINE real pair (rho != 1, phi = 0): {rho, 1/rho}   (the g=1 "|a| > 2 sqrt q" tooth)
        -> t contribution  rho^n + rho^{-n}.
        Toeplitz block: rank 2, signature (1, ., 1)  -- INDEFINITE.

The signature claims follow from  gamma^{|i-j|} + gamma^{-|i-j|} = u_i w_j + w_i u_j
with u_i = gamma^i, w_i = gamma^{-i}: an off-circle |gamma| != 1 "frequency" always
contributes a hyperbolic (1,1) plane, an on-circle one a positive rank-1 projector.
Verified numerically in stage2.py.
"""
import mpmath as mp


def t_from_blocks(blocks, nmax):
    """blocks: list of ('online', phi) / ('quartet', rho, phi) / ('realpair', rho).
    Returns t[0..nmax] as mpmath reals, and g_eff."""
    t = [mp.mpf(0)] * (nmax + 1)
    count = 0
    for b in blocks:
        if b[0] == 'online':
            phi = mp.mpf(b[1]); count += 2
            for n in range(nmax + 1):
                t[n] += 2 * mp.cos(n * phi)
        elif b[0] == 'quartet':
            rho, phi = mp.mpf(b[1]), mp.mpf(b[2]); count += 4
            for n in range(nmax + 1):
                t[n] += 2 * mp.cos(n * phi) * (rho ** n + rho ** (-n))
        elif b[0] == 'realpair':
            rho = mp.mpf(b[1]); count += 2
            for n in range(nmax + 1):
                t[n] += rho ** n + rho ** (-n)
        else:
            raise ValueError(b)
    assert t[0] == count
    return t, count // 2


def betas(blocks):
    """The normalised eigenvalues themselves (for radius/angle reporting)."""
    out = []
    for b in blocks:
        if b[0] == 'online':
            phi = mp.mpf(b[1]); out += [mp.expj(phi), mp.expj(-phi)]
        elif b[0] == 'quartet':
            rho, phi = mp.mpf(b[1]), mp.mpf(b[2])
            out += [rho * mp.expj(phi), rho * mp.expj(-phi),
                    mp.expj(phi) / rho, mp.expj(-phi) / rho]
        elif b[0] == 'realpair':
            rho = mp.mpf(b[1]); out += [mp.mpc(rho), mp.mpc(1 / rho)]
    return out


def predicted_signature(blocks, R):
    """Once R+1 exceeds the total rank, signature should be
    (2*#quartet + 1*#realpair, nullity, 2*#quartet + 1*#realpair + 2*#online)."""
    nq = sum(1 for b in blocks if b[0] == 'quartet')
    nr = sum(1 for b in blocks if b[0] == 'realpair')
    no = sum(1 for b in blocks if b[0] == 'online')
    rank = 4 * nq + 2 * nr + 2 * no
    neg = 2 * nq + nr
    pos = 2 * nq + nr + 2 * no
    if R + 1 >= rank:
        return (neg, R + 1 - rank, pos)
    return None
