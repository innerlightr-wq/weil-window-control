#!/usr/bin/env python3
r"""Background-normalized interaction between the positive and negative spectral
channels of an off-line quartet, on the Davenport-Heilbronn window form.

PROVENANCE.  Everything about the function, the window form, the quartet orbit and the
frozen witnesses is REUSED from extensions/davenport-heilbronn-2026-10-23 (commit
e2cf715): src/dh.py (the arithmetic-side matrix, the transforms G_at, kappa, b_n),
src/phaseD.py (the declared basis levels and index sets), certificates/witness*.json
(the frozen coefficients).  Nothing there is modified; this directory only reads it.

THE REDUCTION (standard finite-dimensional linear algebra, not a new theorem).
In ONE real parity sector the quartet form is a signed rank-two matrix,

    Q_O = 4 p p^T - 4 s s^T,

with p and s REAL vectors -- see `channel_vectors` for the phase removal.  Writing Q for
the complete arithmetic-side matrix on the test space, R = Q - Q_O, and u = 2s,

    B := R + 4 p p^T = Q + u u^T,        Q = B - u u^T.

B is a CANDIDATE positive background: it is the remaining-zero form plus the quartet's
POSITIVE channel only.  Nothing here assumes B > 0; it is tested.  When B > 0,

    eta := u^T B^{-1} u,      Q >= 0  <=>  eta <= 1,

and for several channels in the columns of U, with B = Q + U U^T,

    K := U^T B^{-1} U,        Q >= 0  <=>  lambda_max(K) <= 1,

with #{negative eigenvalues of Q} = #{eigenvalues of K exceeding 1}.  The last identity
also bounds what one channel can explain: a rank-one update can account for at most ONE
negative direction, so if Q has two or more, B cannot be positive definite.

Usage:
    explore_channels.py --dps 80 --windows 1.6,1.8,2.0 --levels 2,3
    explore_channels.py --dps 140 --windows 2.0 --levels 2 --enrichment --twochannel
    explore_channels.py --dps 80 --windows 2.0 --levels 2 --invariance
"""
import argparse
import csv
import json
import os
import sys

DH = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), 'davenport-heilbronn-2026-10-23', 'src')
sys.path.insert(0, DH)

import mpmath as mp                                                      # noqa: E402
import dh                                                                # noqa: E402
from phaseD import index_set, gram                                       # noqa: E402
from zeta_window import freq, norm2                                      # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, 'data')
DHEXT = os.path.dirname(DH)
SECTOR = 'odd'

# the four enumerated off-line quartets, read from the census (numerical inputs: T2,
# not independently certified -- preserved in every conclusion below)
def orbits():
    path = os.path.join(DHEXT, 'results', 'census.json')
    with open(path) as fh:
        c = json.load(fh)
    return [(o['beta'], o['gamma']) for o in c['offline']]


# ----------------------------------------------------------------- phase removal

def channel_vectors(idx, L, sector, beta_s, gamma_s, check=True):
    r"""The REAL channel vectors (p, s) of the quartet, with the phase removed.

    For real g, P = int g cosh(du) e^{i y u} du and S = int g sinh(du) e^{i y u} du, and
    A = G(d+iy), B = G(-d+iy) give P = (A+B)/2, S = (A-B)/2.

    In the ODD sector g*cosh(du) is ODD and g*sinh(du) is EVEN, so

        P = i * int g cosh(du) sin(y u) du        (purely imaginary)
        S =     int g sinh(du) cos(y u) du        (purely real)

    hence |P|^2 = (Im P)^2 and |S|^2 = (Re S)^2 and Q_O = 4(p^2 - s^2) with
    p = Im P, s = Re S both real.  In the EVEN sector the roles swap:
    P is real and S purely imaginary.  `check` verifies that the component that must
    vanish does vanish, rather than assuming it.
    """
    d = dh.mpf_hi(beta_s) - mp.mpf(1) / 2
    g = dh.mpf_hi(gamma_s)
    p, s, dead = [], [], mp.mpf(0)
    for k in idx:
        e = [mp.mpf(0)] * len(idx)
        e[list(idx).index(k)] = mp.mpf(1)
        A = dh.G_at(e, idx, L, sector, mp.mpc(d, g))
        Bv = dh.G_at(e, idx, L, sector, mp.mpc(-d, g))
        P, S = (A + Bv) / 2, (A - Bv) / 2
        if sector == 'odd':
            p.append(mp.im(P)); s.append(mp.re(S))
            dead = max(dead, abs(mp.re(P)), abs(mp.im(S)))
        else:
            p.append(mp.re(P)); s.append(mp.im(S))
            dead = max(dead, abs(mp.im(P)), abs(mp.re(S)))
    if check:
        return p, s, dead
    return p, s


def quad_check_channels(idx, L, sector, beta_s, gamma_s, nmodes=3):
    """Independent check of p, s by direct quadrature on the first few modes."""
    d = dh.mpf_hi(beta_s) - mp.mpf(1) / 2
    g = dh.mpf_hi(gamma_s)
    worst = mp.mpf(0)
    p, s, _ = channel_vectors(idx, L, sector, beta_s, gamma_s)
    for m in range(min(nmodes, len(idx))):
        k = list(idx)[m]
        w = freq(k, L, sector)
        nn = mp.sqrt(norm2(k, L, sector))
        base = (lambda u: mp.sin(w * u) / nn) if sector == 'odd' \
            else (lambda u: mp.cos(w * u) / nn)
        if sector == 'odd':
            pq = mp.quad(lambda u: base(u) * mp.cosh(d * u) * mp.sin(g * u), [-L, 0, L])
            sq = mp.quad(lambda u: base(u) * mp.sinh(d * u) * mp.cos(g * u), [-L, 0, L])
        else:
            pq = mp.quad(lambda u: base(u) * mp.cosh(d * u) * mp.cos(g * u), [-L, 0, L])
            sq = mp.quad(lambda u: base(u) * mp.sinh(d * u) * mp.sin(g * u), [-L, 0, L])
        worst = max(worst, abs(pq - p[m]), abs(sq - s[m]))
    return worst


# ----------------------------------------------------------------- linear algebra

def outer(a, b=None):
    b = a if b is None else b
    n = len(a)
    M = mp.zeros(n, n)
    for i in range(n):
        for j in range(n):
            M[i, j] = a[i] * b[j]
    return M


def madd(*Ms):
    n = Ms[0].rows
    X = mp.zeros(n, n)
    for M in Ms:
        for i in range(n):
            for j in range(n):
                X[i, j] += M[i, j]
    return X


def msub(A, Bm):
    n = A.rows
    X = mp.zeros(n, n)
    for i in range(n):
        for j in range(n):
            X[i, j] = A[i, j] - Bm[i, j]
    return X


def mnorm(A):
    return max(abs(A[i, j]) for i in range(A.rows) for j in range(A.cols))


def chol_solve(Bm, rhs):
    """Solve B x = rhs by Cholesky.  Returns (x, residual, ok).  B^{-1} is never formed.
    Raises nothing: `ok` is False when B is not positive definite."""
    try:
        Lc = mp.cholesky(Bm)
    except Exception:                                                  # noqa: BLE001
        return None, None, False
    x = mp.lu_solve(Bm, mp.matrix(rhs))
    res = mp.norm(Bm * x - mp.matrix(rhs)) / max(mp.norm(mp.matrix(rhs)), mp.mpf(1))
    return x, res, True


def inertia(A):
    ev = mp.eigsy(A, eigvals_only=True)
    ev = sorted(ev)
    neg = sum(1 for e in ev if e < 0)
    return neg, ev


# ----------------------------------------------------------------- one configuration

def run_config(Ls, level, dps, orbit_idx=0, verbose=True):
    mp.mp.dps = dps
    L = dh.mpf_hi(Ls)
    idx = index_set(L, SECTOR, level - 1)
    n = len(idx)
    beta_s, gamma_s = orbits()[orbit_idx]

    Q = dh.build_matrix(L, idx, dps, SECTOR)
    p, s, dead = channel_vectors(idx, L, SECTOR, beta_s, gamma_s)
    u = [2 * x for x in s]

    QO = msub(outer([2 * x for x in p]), outer(u))        # 4pp^T - 4ss^T
    # polarisation cross-check of the rank-two form against the scalar routine
    cc = [mp.mpf(1) / (i + 2) for i in range(n)]
    lhs = mp.fsum([cc[i] * cc[j] * QO[i, j] for i in range(n) for j in range(n)])
    rhs = dh.quartet_contribution(cc, idx, L, SECTOR,
                                  dh.mpf_hi(beta_s), dh.mpf_hi(gamma_s))
    pol = abs(lhs - rhs)

    R = msub(Q, QO)
    Bm = madd(Q, outer(u))
    recon = mnorm(msub(Q, msub(Bm, outer(u))))             # Q - (B - uu^T)
    identity = mnorm(msub(Bm, madd(R, outer([2 * x for x in p]))))   # B - (R + 4pp^T)

    # Gram matrix of the normalised basis.  This must be the identity in this sector,
    # but it is COMPUTED, from the parent extension's phaseD.gram (which pairs the basis
    # through the same cross-correlation machinery the form uses), not asserted.  An
    # earlier version of this line compared norm2(k)/norm2(k) against 1, which is a
    # tautology and could not fail; it is replaced here.
    G = gram(idx, L, SECTOR)
    gram_dev = max(abs(G[i, j] - (1 if i == j else 0))
                   for i in range(n) for j in range(n))

    negQ, evQ = inertia(Q)
    negB, evB = inertia(Bm)
    negR, evR = inertia(R)

    x, res, pd = chol_solve(Bm, u)
    eta = mp.fsum([u[i] * x[i] for i in range(n)]) if pd else None

    # The matrix-determinant lemma: det(B - uu^T) = det(B)(1 - u^T B^{-1} u), i.e.
    #     eta = 1 - det(Q)/det(B).
    # This is why eta can be O(1) even when B has an eigenvalue of order 1e-24: the
    # near-null factors cancel between the two determinants.  It also shows the
    # criterion eta <= 1 IS the sign test det(Q) >= 0 and nothing more.
    detQ, detB = mp.det(Q), mp.det(Bm)
    eta_det = 1 - detQ / detB if detB != 0 else None
    det_res = (abs(eta - eta_det) / max(abs(eta), mp.mpf(1))
               if (pd and eta_det is not None) else None)

    rec = dict(L=Ls, level=level, N=n, dps=dps, orbit='O%d' % (orbit_idx + 1),
               lmin_Q=mp.nstr(evQ[0], 16), lmin_B=mp.nstr(evB[0], 16),
               lmin_R=mp.nstr(evR[0], 16), lmax_B=mp.nstr(evB[-1], 10),
               n_neg_Q=negQ, n_neg_B=negB, n_neg_R=negR,
               B_posdef=int(pd),
               eta=(mp.nstr(eta, 20) if pd else ''),
               eta_gt_1=(int(eta > 1) if pd else ''),
               Q_psd=int(negQ == 0),
               criterion_agrees=(int((eta <= 1) == (negQ == 0)) if pd else ''),
               inertia_agrees=(int(negQ == (1 if eta > 1 else 0)) if pd else ''),
               solve_residual=(mp.nstr(res, 4) if pd else ''),
               recon_residual=mp.nstr(recon, 4),
               B_identity_residual=mp.nstr(identity, 4),
               phase_dead_component=mp.nstr(dead, 4),
               polarisation_residual=mp.nstr(pol, 4),
               gram_dev=mp.nstr(gram_dev, 4),
               cond_B=(mp.nstr(evB[-1] / evB[0], 6) if evB[0] > 0 else 'indef'),
               lam2_Q=mp.nstr(evQ[1], 12), lam2_B=mp.nstr(evB[1], 12),
               det_Q=mp.nstr(detQ, 10), det_B=mp.nstr(detB, 10),
               eta_from_det=(mp.nstr(eta_det, 20) if eta_det is not None else ''),
               det_identity_rel_residual=(mp.nstr(det_res, 4)
                                          if det_res is not None else ''))
    if verbose:
        print('  L=%-4s lvl%d N=%-3d dps=%-4d  lmin(Q)=%-22s neg(Q)=%d  '
              'lmin(B)=%-22s B>0:%s  eta=%-22s'
              % (Ls, level, n, dps, mp.nstr(evQ[0], 12), negQ,
                 mp.nstr(evB[0], 12), 'yes' if pd else 'NO',
                 mp.nstr(eta, 14) if pd else '--'))
    return rec, dict(Q=Q, B=Bm, R=R, u=u, p=p, s=s, idx=idx, L=L, evB=evB,
                     eta=eta, x=x, pd=pd)


# ----------------------------------------------------------------- eta decomposition

def eta_decomposition(ctx, topk=8, cluster_rel=mp.mpf('1e-6')):
    """eta = sum_j |v_j^T u|^2 / beta_j, with clustered eigenvalues grouped."""
    Bm, u = ctx['B'], ctx['u']
    n = len(u)
    E, V = mp.eigsy(Bm)
    order = sorted(range(n), key=lambda i: E[i])
    rows = []
    for i in order:
        ov = mp.fsum([V[k, i] * u[k] for k in range(n)])
        rows.append((E[i], ov, ov ** 2 / E[i] if E[i] != 0 else mp.inf))
    # group clusters
    groups, cur = [], [rows[0]]
    for r in rows[1:]:
        if abs(r[0] - cur[-1][0]) <= cluster_rel * max(abs(r[0]), abs(cur[-1][0])):
            cur.append(r)
        else:
            groups.append(cur); cur = [r]
    groups.append(cur)
    out = []
    for gidx, grp in enumerate(groups):
        beta_lo, beta_hi = grp[0][0], grp[-1][0]
        contrib = mp.fsum([r[2] for r in grp])
        ov2 = mp.fsum([r[1] ** 2 for r in grp])
        out.append(dict(group=gidx, size=len(grp), beta_lo=beta_lo, beta_hi=beta_hi,
                        overlap2=ov2, contribution=contrib))
    out.sort(key=lambda r: -abs(r['contribution']))
    total = mp.fsum([r['contribution'] for r in out])
    return out[:topk], total, mp.fsum([x * x for x in u])


# ----------------------------------------------------------------- enrichment

def sub(M, rows, cols):
    X = mp.zeros(len(rows), len(cols))
    for a, i in enumerate(rows):
        for b, j in enumerate(cols):
            X[a, b] = M[i, j]
    return X


def enrichment(Ls, dps, orbit_idx=0):
    r"""The nested-enrichment identity.  With V subset W at fixed L,

        B_W = [[B_V, C], [C^T, D]],   u_W = (u_V, u_new),
        S = D - C^T B_V^{-1} C,  E = u_new - C^T B_V^{-1} u_V,
        K_W - K_V = E^T S^{-1} E   (>= 0 automatically when B_V, S > 0).

    The nonnegativity is an automatic consequence of the identity.  What is reported is
    the SIZE of the increment, and the separate sizes of E and S that produce it."""
    mp.mp.dps = dps
    L = dh.mpf_hi(Ls)
    V = index_set(L, SECTOR, 1)
    W = index_set(L, SECTOR, 2)
    assert set(V) <= set(W)
    beta_s, gamma_s = orbits()[orbit_idx]
    QW = dh.build_matrix(L, W, dps, SECTOR)
    pW, sW, _ = channel_vectors(W, L, SECTOR, beta_s, gamma_s)
    uW = [2 * x for x in sW]
    BW = madd(QW, outer(uW))
    iV = [list(W).index(k) for k in V]
    iN = [i for i in range(len(W)) if i not in set(iV)]
    BV = sub(BW, iV, iV)
    C = sub(BW, iV, iN)
    D = sub(BW, iN, iN)
    uV = [uW[i] for i in iV]
    uN = [uW[i] for i in iN]

    # consistency: the V-block of B_W must equal the B built intrinsically on V
    QV = dh.build_matrix(L, V, dps, SECTOR)
    pV, sV, _ = channel_vectors(V, L, SECTOR, beta_s, gamma_s)
    BVint = madd(QV, outer([2 * x for x in sV]))
    block_res = mnorm(msub(BV, BVint))
    uV_res = max(abs(uV[i] - 2 * sV[i]) for i in range(len(V)))

    okV = True
    try:
        mp.cholesky(BV)
    except Exception:                                                  # noqa: BLE001
        okV = False
    okW = True
    try:
        mp.cholesky(BW)
    except Exception:                                                  # noqa: BLE001
        okW = False
    if not (okV and okW):
        return dict(L=Ls, dps=dps, orbit='O%d' % (orbit_idx + 1),
                    N_V=len(V), N_W=len(W), B_V_posdef=int(okV),
                    B_W_posdef=int(okW), note='background not positive definite')

    xV = mp.lu_solve(BV, mp.matrix(uV))
    KV = mp.fsum([uV[i] * xV[i] for i in range(len(V))])
    xW = mp.lu_solve(BW, mp.matrix(uW))
    KW = mp.fsum([uW[i] * xW[i] for i in range(len(W))])

    # S = D - C^T B_V^{-1} C  and  E = u_new - C^T B_V^{-1} u_V, by solves
    nV, nN = len(iV), len(iN)
    Y = mp.zeros(nV, nN)                       # Y = B_V^{-1} C
    for b in range(nN):
        col = mp.lu_solve(BV, mp.matrix([C[a, b] for a in range(nV)]))
        for a in range(nV):
            Y[a, b] = col[a]
    S = mp.zeros(nN, nN)
    for a in range(nN):
        for b in range(nN):
            S[a, b] = D[a, b] - mp.fsum([C[k, a] * Y[k, b] for k in range(nV)])
    Ev = [uN[b] - mp.fsum([C[k, b] * xV[k] for k in range(nV)]) for b in range(nN)]
    okS = True
    try:
        mp.cholesky(S)
    except Exception:                                                  # noqa: BLE001
        okS = False
    incr = None
    if okS:
        z = mp.lu_solve(S, mp.matrix(Ev))
        incr = mp.fsum([Ev[b] * z[b] for b in range(nN)])
    evS = sorted(mp.eigsy(S, eigvals_only=True))
    return dict(L=Ls, dps=dps, orbit='O%d' % (orbit_idx + 1),
                N_V=len(V), N_W=len(W), n_new=nN,
                B_V_posdef=int(okV), B_W_posdef=int(okW), S_posdef=int(okS),
                K_V=mp.nstr(KV, 20), K_W=mp.nstr(KW, 20),
                increment_direct=mp.nstr(KW - KV, 20),
                increment_identity=(mp.nstr(incr, 20) if okS else ''),
                identity_residual=(mp.nstr(abs((KW - KV) - incr), 4) if okS else ''),
                norm_E=mp.nstr(mp.sqrt(mp.fsum([x * x for x in Ev])), 10),
                lmin_S=mp.nstr(evS[0], 10), lmax_S=mp.nstr(evS[-1], 10),
                block_residual=mp.nstr(block_res, 4),
                uV_residual=mp.nstr(uV_res, 4),
                note='')


# ----------------------------------------------------------------- two channels

def two_channel(Ls, level, dps, oi=(0, 1)):
    mp.mp.dps = dps
    L = dh.mpf_hi(Ls)
    idx = index_set(L, SECTOR, level - 1)
    n = len(idx)
    Q = dh.build_matrix(L, idx, dps, SECTOR)
    cols = []
    for k in oi:
        b_s, g_s = orbits()[k]
        _, s, _ = channel_vectors(idx, L, SECTOR, b_s, g_s)
        cols.append([2 * x for x in s])
    Bm = madd(Q, *[outer(c) for c in cols])
    negQ, evQ = inertia(Q)
    negB, evB = inertia(Bm)
    pd = True
    try:
        mp.cholesky(Bm)
    except Exception:                                                  # noqa: BLE001
        pd = False
    if not pd:
        return dict(L=Ls, level=level, N=n, dps=dps,
                    orbits='+'.join('O%d' % (k + 1) for k in oi),
                    B_posdef=0, note='common background not positive definite',
                    lmin_Q=mp.nstr(evQ[0], 16), n_neg_Q=negQ)
    r = len(cols)
    K = mp.zeros(r, r)
    X = []
    for c in cols:
        X.append(mp.lu_solve(Bm, mp.matrix(c)))
    for a in range(r):
        for b in range(r):
            K[a, b] = mp.fsum([cols[a][i] * X[b][i] for i in range(n)])
    evK = sorted(mp.eigsy(K, eigvals_only=True))
    diagonly = max(K[0, 0], K[1, 1])
    return dict(L=Ls, level=level, N=n, dps=dps,
                orbits='+'.join('O%d' % (k + 1) for k in oi),
                B_posdef=1, lmin_Q=mp.nstr(evQ[0], 16), n_neg_Q=negQ,
                lmin_B=mp.nstr(evB[0], 12),
                K11=mp.nstr(K[0, 0], 18), K22=mp.nstr(K[1, 1], 18),
                K12=mp.nstr(K[0, 1], 18),
                K_symmetry=mp.nstr(abs(K[0, 1] - K[1, 0]), 4),
                lmax_K=mp.nstr(evK[-1], 18), lmin_K=mp.nstr(evK[0], 18),
                max_diag=mp.nstr(diagonly, 18),
                offdiag_effect=mp.nstr(evK[-1] - diagonly, 10),
                n_K_above_1=sum(1 for e in evK if e > 1),
                inertia_agrees=int(negQ == sum(1 for e in evK if e > 1)),
                note='')


# ----------------------------------------------------------------- invariance

def invariance(Ls, level, dps, orbit_idx=0):
    """eta must be invariant under an invertible change of coordinates on the same
    space.  Q -> T^T Q T, u -> T^T u, and the generalized eigenproblem carries
    G -> T^T G T."""
    mp.mp.dps = dps
    _, ctx = run_config(Ls, level, dps, orbit_idx, verbose=False)
    Q, u, n = ctx['Q'], ctx['u'], len(ctx['u'])
    # a fixed, reproducible, badly scaled but invertible T (no randomness)
    T = mp.zeros(n, n)
    for i in range(n):
        T[i, i] = mp.mpf(10) ** (mp.mpf(i % 7) - 3)
        if i + 1 < n:
            T[i, i + 1] = mp.mpf(1) / 3
        if i > 0:
            T[i, i - 1] = -mp.mpf(1) / 7
    Qp = T.T * Q * T
    up = T.T * mp.matrix(u)
    Bp = Qp + up * up.T
    xp = mp.lu_solve(Bp, up)
    etap = mp.fsum([up[i] * xp[i] for i in range(n)])
    Gp = T.T * T
    # generalized lmin of (Qp, Gp): with Gp = Lg Lg^T, the pencil (Qp, Gp) has the same
    # eigenvalues as the ordinary symmetric matrix Z = Lg^{-1} Qp Lg^{-T}.  mpmath's
    # lu_solve takes a VECTOR right-hand side, so the two triangular solves are done
    # column by column (an earlier version passed a matrix and raised IndexError).
    Lg = mp.cholesky(Gp)
    Y = mp.zeros(n, n)
    for j in range(n):
        col = mp.lu_solve(Lg, mp.matrix([Qp[i, j] for i in range(n)]))
        for i in range(n):
            Y[i, j] = col[i]
    Z = mp.zeros(n, n)
    for j in range(n):
        col = mp.lu_solve(Lg, mp.matrix([Y[j, i] for i in range(n)]))
        for i in range(n):
            Z[j, i] = col[i]
    ev = sorted(mp.eigsy((Z + Z.T) / 2, eigvals_only=True))
    _, evQ = inertia(Q)
    return dict(L=Ls, level=level, N=n, dps=dps, orbit='O%d' % (orbit_idx + 1),
                cond_T=mp.nstr(mp.mpf(10) ** 6, 3),
                eta_orthonormal=mp.nstr(ctx['eta'], 20),
                eta_transformed=mp.nstr(etap, 20),
                eta_rel_diff=mp.nstr(abs(etap - ctx['eta'])
                                     / abs(ctx['eta']), 4),
                lmin_Q_orthonormal=mp.nstr(evQ[0], 16),
                lmin_Q_generalized=mp.nstr(ev[0], 16),
                lmin_rel_diff=mp.nstr(abs(ev[0] - evQ[0])
                                      / max(abs(evQ[0]), mp.mpf(1)), 4))


# ----------------------------------------------------------------- witness benchmark

def witness_benchmark(dps):
    """Reproduce the two frozen witness values from their STORED coefficients, at this
    precision.  Benchmark signs only; nothing is refitted."""
    out = []
    for wf in ('witness.json', 'witness_smallest_L.json'):
        with open(os.path.join(DHEXT, 'certificates', wf)) as fh:
            c = json.load(fh)
        mp.mp.dps = dps
        L = dh.mpf_hi(c['L'])
        idx = c['idx']
        v = [dh.mpf_hi(t) for t in c['coefficients']]
        n2 = mp.fsum([x * x for x in v])
        M = dh.build_matrix(L, idx, dps, c['sector'])
        q = mp.fsum([v[i] * v[j] * M[i, j]
                     for i in range(len(idx)) for j in range(len(idx))]) / n2
        out.append(dict(witness=wf, L=c['L'], sector=c['sector'], N=len(idx), dps=dps,
                        value=mp.nstr(q, 30), stored=c['value_dps'][
                            max(c['value_dps'], key=int)],
                        negative=int(q < 0)))
    return out


# ----------------------------------------------------------------- main

def write_csv(path, rows):
    if not rows:
        return
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in keys})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dps', default='80')
    ap.add_argument('--windows', default='1.6,1.8,2.0')
    ap.add_argument('--levels', default='2,3')
    ap.add_argument('--decomp', default='')
    ap.add_argument('--enrichment', default='')
    ap.add_argument('--twochannel', default='')
    ap.add_argument('--invariance', default='')
    ap.add_argument('--witness', action='store_true')
    ap.add_argument('--out', default='channel_results.csv')
    a = ap.parse_args()
    os.makedirs(DATA, exist_ok=True)
    dpss = [int(x) for x in a.dps.split(',')]

    if a.witness:
        print('=== frozen-witness benchmark (stored coefficients, re-parsed) ===')
        for dps in dpss:
            for r in witness_benchmark(dps):
                print('  %-26s L=%-4s N=%-3d dps=%-4d  Q/||g||^2 = %-28s %s'
                      % (r['witness'], r['L'], r['N'], r['dps'], r['value'],
                         'NEGATIVE' if r['negative'] else 'POSITIVE'))

    rows = []
    if a.windows and a.levels:
        print('=== channel reduction ===')
        for dps in dpss:
            for Ls in a.windows.split(','):
                for lev in [int(x) for x in a.levels.split(',')]:
                    rec, _ = run_config(Ls, lev, dps)
                    rows.append(rec)
        write_csv(os.path.join(DATA, a.out), rows)
        print('  -> data/%s' % a.out)

    if a.decomp:
        Ls, lev = a.decomp.split(':')
        print('=== eta decomposition at L=%s level %s ===' % (Ls, lev))
        for dps in dpss:
            _, ctx = run_config(Ls, int(lev), dps, verbose=False)
            if not ctx['pd']:
                print('  dps=%d: background not positive definite, skipped' % dps)
                continue
            groups, total, u2 = eta_decomposition(ctx)
            print('  dps=%-4d  eta = %s   ||u||^2 = %s'
                  % (dps, mp.nstr(ctx['eta'], 20), mp.nstr(u2, 12)))
            print('  %-6s %-5s %-24s %-24s %-24s %-10s'
                  % ('group', 'size', 'beta (lo..hi)', '|v^T u|^2', 'contribution',
                     '% of eta'))
            for g in groups:
                print('  %-6d %-5d %-24s %-24s %-24s %-10s'
                      % (g['group'], g['size'],
                         mp.nstr(g['beta_lo'], 8) + ('' if g['size'] == 1
                                                     else '..' + mp.nstr(g['beta_hi'], 4)),
                         mp.nstr(g['overlap2'], 12), mp.nstr(g['contribution'], 12),
                         mp.nstr(100 * g['contribution'] / total, 6)))
            print('  sum of listed groups / eta = %s'
                  % mp.nstr(mp.fsum([g['contribution'] for g in groups]) / total, 12))

    if a.enrichment:
        print('=== nested enrichment V(level 2) subset W(level 3) ===')
        erows = []
        for dps in dpss:
            for Ls in a.enrichment.split(','):
                r = enrichment(Ls, dps)
                erows.append(r)
                for k, v in r.items():
                    print('    %-22s %s' % (k, v))
                print('')
        write_csv(os.path.join(DATA, 'enrichment_results.csv'), erows)
        print('  -> data/enrichment_results.csv')

    if a.twochannel:
        print('=== two channels O_1 + O_2, common background ===')
        trows = []
        for dps in dpss:
            for spec in a.twochannel.split(','):
                Ls, lev = spec.split(':')
                r = two_channel(Ls, int(lev), dps)
                trows.append(r)
                for k, v in r.items():
                    print('    %-22s %s' % (k, v))
                print('')
        write_csv(os.path.join(DATA, 'twochannel_results.csv'), trows)
        print('  -> data/twochannel_results.csv')

    if a.invariance:
        print('=== coordinate invariance ===')
        irows = []
        for dps in dpss:
            for spec in a.invariance.split(','):
                Ls, lev = spec.split(':')
                r = invariance(Ls, int(lev), dps)
                irows.append(r)
                for k, v in r.items():
                    print('    %-22s %s' % (k, v))
                print('')
        write_csv(os.path.join(DATA, 'invariance_results.csv'), irows)
    return 0


if __name__ == '__main__':
    sys.exit(main())
