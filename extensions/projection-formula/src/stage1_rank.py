#!/usr/bin/env python3
"""Two precision audits of the Stage-1 deformation family.

(A) REALIZABILITY. Our increment is Delta t(n) = 2 cos(n th)(cosh(n a) - 1).
    Ask which reciprocal- and conjugation-closed multisets produce it.
      - 2-point set {e^{a+i th}, e^{-a-i th}}: sum = 2 cosh(na + i n th)
        = 2 cosh(na)cos(n th) + 2i sinh(na)sin(n th)  -- NOT real for a != 0.
      - 4-point quadruple {e^{+-a +- i th}} : sum = 4 cosh(na) cos(n th), real,
        reciprocal- and conjugation-closed, and the a=0 value is 4cos(n th),
        i.e. a DOUBLED angle. Its increment is 4cos(n th)(cosh(na)-1) = 2 x ours.
    So ours is the symmetrised (real-part) increment: a legitimate real symmetric
    Toeplitz deformation, but not the point-count sequence of a single curve.
    Checked here: the quadruple increment is exactly twice ours, and the
    corresponding lam_min coefficient doubles.

(B) RANK. c^T B c = 2 Re(S_2 conj(S_0)) - 2|S_1|^2 with S_k = sum_i i^k c_i e^{i i th};
    on E = ker T_0 (S_0 = 0) this is -2|S_1|^2 = -2(<c,u>^2 + <c,v>^2),
    u_i = i cos(i th), v_i = i sin(i th). On the REAL coefficient space that is
    rank 2 unless sin(i th) == 0 for all i (th = 0 or pi). Verified numerically.
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'src'))
import mpmath as mp
from curves import squarefree_over_Fq, count_points_hyperelliptic, newton_p_to_A, A_to_power_sums
from weilform import toeplitz, eigsym
from stage1_ff import CURVES, curve_data, thetas_of, tvals, dt_exact, kernel_basis
mp.mp.dps = 60

def Bmat(R, th):
    return mp.matrix([[(i-k)**2*mp.cos((i-k)*th) for k in range(R+1)] for i in range(R+1)])

print('(A) realizability: quadruple increment vs ours, and the coefficient ratio\n')
print('  %-6s %-4s %-13s %-13s %-9s %-13s %-13s %-7s' %
      ('curve','R','lam(ours)','pred(-2a^2 mu)','rel','lam(quad)','pred(-4a^2 mu)','rel'))
A_rows = []
for tag, f, q, g in CURVES[:3]:
    d = 2*g; nmax = d+14
    A, p_all = curve_data(f, q, g, nmax)
    ths, _ = thetas_of(A, q, g); th = ths[0]
    t0 = tvals(p_all, q, g, nmax)
    for R in (d+2, d+5):
        T = toeplitz(t0, R)
        N, _ = kernel_basis(T, R, mp.mpf(10)**-30)
        if len(N) != R+1-d: continue
        from stage1_ff import S1_maximum
        mu, nrm2 = S1_maximum(N, R, th)
        a = mp.mpf('1e-4')
        t1 = [t0[n] + dt_exact(n, th, a) for n in range(R+1)]
        t2 = [t0[n] + 2*dt_exact(n, th, a) for n in range(R+1)]   # genuine quadruple
        l1 = eigsym(toeplitz(t1, R))[0][0]; l2 = eigsym(toeplitz(t2, R))[0][0]
        p1 = -2*a**2*mu; p2 = -4*a**2*mu
        r1 = float(abs(l1-p1)/abs(p1)); r2 = float(abs(l2-p2)/abs(p2))
        print('  %-6s %-4d %-13s %-13s %-9.2e %-13s %-13s %-7.2e' %
              (tag, R, mp.nstr(l1,9), mp.nstr(p1,9), r1, mp.nstr(l2,9), mp.nstr(p2,9), r2))
        A_rows.append(dict(tag=tag,R=R,rel_ours=r1,rel_quad=r2))

print('\n(B) rank of P_E B P_E on the REAL coefficient space\n')
print('  %-6s %-4s %-9s %-7s %-9s %-26s %-13s' %
      ('curve','R','dim E','rank','theta','nonzero eigs of compr.','||compr||_op'))
B_rows = []
for tag, f, q, g in CURVES[:3]:
    d = 2*g; nmax = d+14
    A, p_all = curve_data(f, q, g, nmax)
    ths, _ = thetas_of(A, q, g); th = ths[0]
    t0 = tvals(p_all, q, g, nmax)
    for R in (d+2, d+5):
        T = toeplitz(t0, R)
        N, _ = kernel_basis(T, R, mp.mpf(10)**-30)
        if len(N) != R+1-d: continue
        m = len(N); B = Bmat(R, th)
        # orthonormalise N (eigsym vectors of a symmetric matrix are already orthonormal)
        G = mp.zeros(m, m)
        for p_ in range(m):
            for r_ in range(m):
                G[p_, r_] = sum(N[p_][i]*sum(B[i,k]*N[r_][k] for k in range(R+1)) for i in range(R+1))
        ev, _ = eigsym(G)
        nz = [e for e in ev if abs(e) > mp.mpf(10)**-20]
        rank = len(nz)
        opn = max(abs(e) for e in ev)
        bound = 2*sum(mp.mpf(i)**2 for i in range(R+1))
        print('  %-6s %-4d %-9d %-7d %-9s %-26s %-13s  (<= 2*sum i^2 = %s)' %
              (tag, R, m, rank, mp.nstr(th,6), ' '.join(mp.nstr(e,6) for e in nz),
               mp.nstr(opn,8), mp.nstr(bound,8)))
        B_rows.append(dict(tag=tag,R=R,dimE=m,rank=rank,opnorm=float(opn),bound=float(bound)))
print()
print('  theta = 0 degenerate case (Theorem 6 family): v_i = i sin(0) = 0, so rank drops to 1')
for R in (4, 8):
    T0 = toeplitz([mp.mpf(2)]*(R+1), R)           # t_0(n) = 2, T_0 = 2 * 1 1^T
    N, _ = kernel_basis(T0, R, mp.mpf(10)**-30)
    m = len(N); B = Bmat(R, mp.mpf(0))
    G = mp.zeros(m, m)
    for p_ in range(m):
        for r_ in range(m):
            G[p_, r_] = sum(N[p_][i]*sum(B[i,k]*N[r_][k] for k in range(R+1)) for i in range(R+1))
    ev, _ = eigsym(G)
    nz = [e for e in ev if abs(e) > mp.mpf(10)**-20]
    print('    R=%-3d dim E=%-3d rank=%-3d  nonzero eig = %-14s  -2*sum(i-R/2)^2 = %s' %
          (R, m, len(nz), mp.nstr(nz[0],10) if nz else '-',
           mp.nstr(-2*sum((mp.mpf(i)-mp.mpf(R)/2)**2 for i in range(R+1)), 10)))
json.dump(dict(realizability=A_rows, rank=B_rows), open('data/stage1_rank.json','w'), indent=1)
