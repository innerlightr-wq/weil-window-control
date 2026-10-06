#!/usr/bin/env python3
"""Task 1e: T_R as a Gram matrix in the trace form on H^1.

Let F be Frobenius on H^1 (dim 2g), char poly  Q(T) = prod_j (T - alpha_j)
= sum_i A_i T^{2g-i}  (the REVERSED L-polynomial, integer entries -> companion matrix).

The polarization gives a nondegenerate Hermitian S with  F* S F = q S, and the adjoint
is  phi^dag := S^{-1} phi* S.  Then F^dag = q F^{-1}, so U := F/sqrt(q) satisfies
U^dag = U^{-1}, and

    Gram_{ij} = Tr( U^i (U^j)^dag ) = Tr( U^{i-j} ) = sum_m beta_m^{i-j} = t(|i-j|),

i.e.  T_R = Gram matrix of {1, F/sqrt q, ..., F^R/q^{R/2}} in the trace form.  (Cor. 14
is the g = 1, R = 1 case.)

WHERE RH ENTERS.  Taking V with F = V diag(alpha) V^{-1}, the choice S = (V V*)^{-1}
satisfies F* S F = q S precisely because |alpha_j|^2 = q for EVERY j, and that S is
Hermitian POSITIVE DEFINITE.  So:

    all |alpha_j| = sqrt(q)   <=>   the q-isometry form S can be chosen positive definite
                              <=>   phi -> Tr(phi phi^dag) is a positive form
                              <=>   T_R is PSD for all R.

A planted off-line configuration still admits an S with F* S F = q S (the functional
equation guarantees that), but every such S is INDEFINITE -- which is exactly the
Stage 2 loss of positivity, seen on the H^1 side.  Checked below.
"""
import sys, os, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from weilform import toeplitz, inertia

mp.mp.dps = 40


def companion(monic_high_first):
    """Companion matrix of a monic polynomial given highest-degree-first."""
    c = [mp.mpf(x) for x in monic_high_first]
    n = len(c) - 1
    M = mp.zeros(n, n)
    for i in range(1, n):
        M[i, i - 1] = 1
    for i in range(n):
        M[i, n - 1] = -c[n - i]          # last column = -(a_0, a_1, ..., a_{n-1})
    return M


def eig_decomp(M):
    E, V = mp.eig(mp.matrix(M))
    return E, V


def pairing_in_eigenbasis(betas, q):
    """Solve  F* S F = q S  in the EIGENBASIS of F.

    With F = diag(beta_j), the constraint reads  conj(beta_j) beta_k S_jk = q S_jk,
    so S_jk = 0 unless  beta_k = q / conj(beta_j).  The map  sigma: j -> k  is an
    involution (the functional equation), so a Hermitian solution always exists:
    put 1 on each pair.

    KEY CONSEQUENCE (T1):  S_jj can be nonzero only if |beta_j|^2 = q.  If ANY
    eigenvalue is off the circle, the corresponding diagonal entry of S is forced to
    zero, and a Hermitian matrix with a zero diagonal entry is never positive definite.
    So: a POSITIVE DEFINITE q-isometry form exists  <=>  all |beta_j| = sqrt(q).
    """
    n = len(betas)
    sig = {}
    for j, bj in enumerate(betas):
        tgt = q / mp.conj(bj)
        k = min(range(n), key=lambda k: abs(betas[k] - tgt))
        assert abs(betas[k] - tgt) < mp.mpf(10) ** (-(mp.mp.dps // 2)), \
            "spectrum is not closed under beta -> q/conj(beta)"
        sig[j] = k
    S = mp.zeros(n, n)
    for j in range(n):
        S[j, sig[j]] = 1
    forced_zero_diag = [j for j in range(n) if sig[j] != j]
    return S, sig, forced_zero_diag


def build_S(V):
    """S = (V V*)^{-1}, Hermitian; positive definite iff V is invertible."""
    n = V.rows
    Vs = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            Vs[i, j] = mp.conj(V[j, i])
    return (V * Vs) ** -1


def herm_inertia(S):
    """Inertia of a Hermitian matrix via its (real) eigenvalues."""
    n = S.rows
    # Hermitian -> real 2n x 2n symmetric embedding [[Re, -Im],[Im, Re]], eigenvalues doubled
    B = mp.zeros(2 * n, 2 * n)
    for i in range(n):
        for j in range(n):
            B[i, j] = mp.re(S[i, j]); B[i, n + j] = -mp.im(S[i, j])
            B[n + i, j] = mp.im(S[i, j]); B[n + i, n + j] = mp.re(S[i, j])
    ev, _ = mp.eigsy(B)
    vals = sorted(ev[k] for k in range(2 * n))
    scale = max([abs(v) for v in vals] + [mp.mpf(1)])
    tol = scale * mp.mpf(10) ** (-(mp.mp.dps - 12))
    return (sum(1 for v in vals if v < -tol) // 2,
            sum(1 for v in vals if abs(v) <= tol) // 2,
            sum(1 for v in vals if v > tol) // 2)


def check(tag, charpoly_high_first, q, g, t_ref, Rmax, rows):
    """charpoly_high_first = [A_0=1, A_1, ..., A_2g] (char poly of F)."""
    F = companion(charpoly_high_first)
    E, V = eig_decomp(F)
    alphas = [E[k] for k in range(2 * g)]
    radii = [abs(a) / mp.sqrt(q) for a in alphas]
    Sd, sig, forced = pairing_in_eigenbasis(alphas, mp.mpf(q))
    Vi = V ** -1
    Vis = mp.matrix(2 * g, 2 * g)
    for i in range(2 * g):
        for j in range(2 * g):
            Vis[i, j] = mp.conj(Vi[j, i])
    S = Vis * Sd * Vi                 # S = V^{-*} S_diag V^{-1}
    sq = mp.sqrt(mp.mpf(q))
    U = F / sq
    # S-adjoint
    Ss = S ** -1
    def dag(M):
        n = M.rows
        Mc = mp.matrix(n, n)
        for i in range(n):
            for j in range(n):
                Mc[i, j] = mp.conj(M[j, i])
        return Ss * Mc * S
    Ud = dag(U)
    # residual of U U^dag = I
    I = mp.eye(2 * g)
    uu = U * Ud - I
    res_unitary = max(abs(uu[i, j]) for i in range(2 * g) for j in range(2 * g))
    # residual of F* S F = q S
    Fc = mp.matrix(2 * g, 2 * g)
    for i in range(2 * g):
        for j in range(2 * g):
            Fc[i, j] = mp.conj(F[j, i])
    fsf = Fc * S * F - q * S
    nS = max(abs(S[i, j]) for i in range(2 * g) for j in range(2 * g))
    res_pol = max(abs(fsf[i, j]) for i in range(2 * g) for j in range(2 * g)) / nS
    iS = herm_inertia(S)
    nforced = len(forced)

    # Gram matrix  G_ij = Tr( U^i (U^j)^dag )
    pw = [mp.eye(2 * g)]
    for k in range(1, Rmax + 1):
        pw.append(pw[-1] * U)
    pwd = [dag(P) for P in pw]
    G = mp.matrix(Rmax + 1, Rmax + 1)
    for i in range(Rmax + 1):
        for j in range(Rmax + 1):
            M = pw[i] * pwd[j]
            G[i, j] = sum(M[k, k] for k in range(2 * g))
    T = toeplitz([mp.mpf(x) for x in t_ref], Rmax)
    res_gram = max(abs(G[i, j] - T[i, j]) for i in range(Rmax + 1) for j in range(Rmax + 1))
    imag = max(abs(mp.im(G[i, j])) for i in range(Rmax + 1) for j in range(Rmax + 1))

    print(f"  {tag:18s} q={q} g={g}  max||alpha|/sqrt(q)-1| = {mp.nstr(max(abs(r-1) for r in radii),4):>10s}")
    print(f"      F* S F = q S residual        : {mp.nstr(res_pol, 4)}")
    print(f"      U U^dag = I residual         : {mp.nstr(res_unitary, 4)}")
    print(f"      inertia of S (neg,zero,pos)  : {iS}   -> "
          f"{'POSITIVE DEFINITE' if iS[0]==0 and iS[1]==0 else 'INDEFINITE'}"
          f"   (diagonal entries forced to 0 in eigenbasis: {nforced}/{2*g})")
    print(f"      max |Gram_ij - T_ij|, R<={Rmax} : {mp.nstr(res_gram, 4)}   (max |Im Gram| = {mp.nstr(imag,4)})")
    rows.append(dict(tag=tag, q=q, g=g, Rmax=Rmax,
                     max_radius_dev=mp.nstr(max(abs(r - 1) for r in radii), 6),
                     res_FSF_qS=mp.nstr(res_pol, 6),
                     res_U_Udag_I=mp.nstr(res_unitary, 6),
                     S_inertia=f"{iS}", S_posdef=(iS[0] == 0 and iS[1] == 0),
                     n_forced_zero_diag=nforced,
                     max_abs_Gram_minus_T=mp.nstr(res_gram, 6),
                     max_abs_Im_Gram=mp.nstr(imag, 6)))
    return res_gram


def main():
    from curves import count_points_hyperelliptic, newton_p_to_A, A_to_power_sums
    from stage1 import CURVES
    from exact_inertia import rational_t
    rows = []
    print("TASK 1e: T_R = Gram matrix of {1, F/sqrt q, ..., F^R/q^{R/2}} in Tr(phi psi^dag)\n")
    print("--- honest curves (RH for curves is a theorem here) ---")
    for tag, desc, f, q, g in CURVES:
        p = [None] + [1 + q ** n - count_points_hyperelliptic(f, q, n)
                      for n in range(1, g + 1)]
        A = newton_p_to_A(p, q, g)
        Rmax = 2 * g + 2
        pall = A_to_power_sums(A, Rmax)
        sq = mp.sqrt(mp.mpf(q))
        t = [mp.mpf(2 * g)] + [mp.mpf(pall[n]) / sq ** n for n in range(1, Rmax + 1)]
        check(tag, A, q, g, t, Rmax, rows)

    print("\n--- planted off-line configurations (same construction, fake spectra) ---")
    # charpoly of the planted 'Frobenius': prod (T - sqrt(q) beta_j) with q = 1 for
    # the normalised picture (beta_j are the normalised eigenvalues, so take q = 1).
    from planted import betas
    for tag, blocks in [('C0-online',   [('online', 0.7), ('online', 2.0)]),
                        ('C1-g1real',   [('realpair', '1.3')]),
                        ('C2-quartet',  [('quartet', '1.3', 0.7)]),
                        ('C4-mixed',    [('quartet', '1.3', 0.7), ('online', 2.0)])]:
        bs = betas(blocks)
        geff = len(bs) // 2
        # monic char poly with roots beta_j, highest-first
        poly = [mp.mpc(1)]
        for b in bs:
            poly = [poly[0]] + [poly[k] - b * poly[k - 1] for k in range(1, len(poly))] + \
                   [-b * poly[-1]]
        poly = [mp.re(c) for c in poly]          # real: the set is conj-closed
        from planted import t_from_blocks
        Rmax = 2 * geff + 2
        t, _ = t_from_blocks(blocks, Rmax)
        check(tag, poly, 1, geff, t, Rmax, rows)

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       'data', 'stage1e_gram_trace.csv')
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {out}")


if __name__ == '__main__':
    main()
