#!/usr/bin/env python3
"""Stage 1 (GATE). Function field, exact where possible.

t(n) = p_n q^{-n/2} = sum_{j=1..g} 2 cos(n theta_j),  T_R[i,k] = t(|i-k|).
Identity (T1):  c^T T_R c = sum_j 2 |P_c(e^{i theta_j})|^2,  P_c(z) = sum_i c_i z^i.
Hence ker T_R = { c : P_c(e^{i theta_j}) = 0 for all j }, of dimension R+1-2g for R >= 2g.

Perturbation (preregistered): one conjugate pair moved off the circle with its
functional-equation partners, so
    Delta t(n) = 2 cos(n theta_j) (cosh(n a) - 1) = a^2 n^2 cos(n theta_j) + O(a^4).

Writing S_k(c) = sum_i i^k c_i e^{i i theta_j}, the a^2 form is (T1, derived in REPORT.md)
    Q_2(c) = 2 Re( S_2(c) conj(S_0(c)) ) - 2 |S_1(c)|^2,
which on ker T_R (where S_0 = 0) is  -2|S_1(c)|^2 = -2 |P_c'(e^{i theta_j})|^2.
"""
import sys, os, json, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'src'))
import mpmath as mp
from curves import squarefree_over_Fq, count_points_hyperelliptic, newton_p_to_A, A_to_power_sums
from weilform import toeplitz, eigsym

CURVES = [('E5', [0,1,0,1], 5, 1), ('H2F3', [1,0,1,0,0,1], 3, 2),
          ('G3F5', [1,1,0,0,0,0,0,1], 5, 3), ('G3F3', [1,1,0,1,0,0,0,1], 3, 3)]


def curve_data(f, q, g, nmax):
    assert squarefree_over_Fq(f, q)
    p_br = {n: 1 + q**n - count_points_hyperelliptic(f, q, n) for n in range(1, g+1)}
    A = newton_p_to_A([None] + [p_br[n] for n in range(1, g+1)], q, g)
    p_all = A_to_power_sums(A, nmax)          # exact integers
    return A, p_all


def thetas_of(A, q, g):
    rev = [mp.mpf(A[i]) for i in range(2*g+1)]
    roots = mp.polyroots(rev, maxsteps=400, extraprec=400)
    sq = mp.sqrt(q)
    out = sorted((mp.arg(r), abs(r)/sq) for r in roots)
    return [th for th, m in out if th > 1e-30], [m for th, m in out]


def tvals(p_all, q, g, nmax):
    sq = mp.sqrt(mp.mpf(q))
    return [mp.mpf(2*g)] + [mp.mpf(p_all[n])/sq**n for n in range(1, nmax+1)]


def dt_exact(n, th, a):
    return 2*mp.cos(n*th)*(mp.cosh(n*a) - 1)


def kernel_basis(T, R, tol):
    vals, vecs = eigsym(T)
    return [mp.matrix(vecs[k]) for k in range(len(vals)) if abs(vals[k]) <= tol], vals


def S1_maximum(N, R, th):
    """max_{c in span(N), ||c||=1} |S_1(c)|^2, S_1(c) = sum_i i c_i e^{i i th}.
    = largest eigenvalue of Pi_N (u u^T + v v^T) Pi_N,  u_i = i cos(i th), v_i = i sin(i th)."""
    u = mp.matrix([i*mp.cos(i*th) for i in range(R+1)])
    v = mp.matrix([i*mp.sin(i*th) for i in range(R+1)])
    m = len(N)
    G = mp.zeros(m, m)
    for p in range(m):
        for r in range(m):
            G[p, r] = (sum(N[p][i]*u[i] for i in range(R+1)) * sum(N[r][i]*u[i] for i in range(R+1))
                     + sum(N[p][i]*v[i] for i in range(R+1)) * sum(N[r][i]*v[i] for i in range(R+1)))
    ev, _ = eigsym(G)
    return ev[-1], (sum(x**2 for x in u) + sum(x**2 for x in v))


def run(dps, out):
    mp.mp.dps = dps
    rows = []
    for tag, f, q, g in CURVES:
        d = 2*g
        nmax = d + 14
        A, p_all = curve_data(f, q, g, nmax)
        ths, mods = thetas_of(A, q, g)
        th = ths[0]                                  # move the first conjugate pair
        t0 = tvals(p_all, q, g, nmax)
        # identity check: c^T T c = sum_j 2 |P_c(e^{i th_j})|^2
        R0 = d + 3
        T0 = toeplitz(t0, R0)
        import random; random.seed(0)
        c = mp.matrix([mp.mpf(random.randint(-9, 9)) for _ in range(R0+1)])
        lhs = (c.T * T0 * c)[0]
        rhs = sum(2*abs(sum(c[i]*mp.e**(1j*i*tj) for i in range(R0+1)))**2 for tj in ths)
        ident = abs(lhs - rhs)/abs(lhs)
        for R in range(d, d+11):
            T = toeplitz(t0, R)
            tol = mp.mpf(10)**(-(dps//2))
            N, vals = kernel_basis(T, R, tol)
            if len(N) != R+1-d:
                rows.append(dict(tag=tag, g=g, R=R, note='kernel dim %d != %d' % (len(N), R+1-d)))
                continue
            mu, nrm2 = S1_maximum(N, R, th)
            for a in ('1e-2', '1e-3', '1e-4'):
                aa = mp.mpf(a)
                ta = [t0[n] + dt_exact(n, th, aa) for n in range(R+1)]
                lam = eigsym(toeplitz(ta, R))[0][0]
                pred = -2*aa**2*mu
                rel = abs(lam - pred)/abs(pred)
                rows.append(dict(tag=tag, g=g, R=R, a=str(a), lam=mp.nstr(lam, 12),
                                 pred=mp.nstr(pred, 12), rel=float(rel),
                                 mu=mp.nstr(mu, 10), cos2=float(mu/nrm2),
                                 kerdim=len(N), ident=float(ident), dps=dps))
        print('  %-6s g=%d  identity rel.err %.2e  theta_1=%s' % (tag, g, float(ident), mp.nstr(th, 10)), flush=True)
    json.dump(rows, open(out, 'w'), indent=1)
    return rows


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--dps', type=int, default=60)
    ap.add_argument('--out', default='data/stage1_ff.json')
    A = ap.parse_args()
    rows = run(A.dps, A.out)
    print('\n  %-6s %-3s %-4s %-7s %-16s %-16s %-11s %-8s' %
          ('curve', 'g', 'R', 'a', 'lam_min', 'prediction', 'rel.err', 'cos^2'))
    for r in rows:
        if 'a' not in r: print('  ', r); continue
        print('  %-6s %-3d %-4d %-7s %-16s %-16s %-11.3e %-8.5f' %
              (r['tag'], r['g'], r['R'], r['a'], r['lam'], r['pred'], r['rel'], r['cos2']))
