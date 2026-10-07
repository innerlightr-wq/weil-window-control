#!/usr/bin/env python3
"""Stage 2: near-null projection (H1) and the delta-staircase (H2), zeta side.

Q_0 = Q_zeta + 2 F(g)^2  (paper eq. 3) = build_matrix + move_block(delta=0): verified.
Q_delta = build_matrix + move_block(delta).
C(delta) = (lam_min(Q_0) - lam_min(Q_delta)) / delta^2.
H1: C(delta) ~ 4 K(L,g) cos^2 theta_{N_delta},  N_delta = eigenvectors of Q_0 with
    eigenvalue <= 4K delta^2 (preregistered rule), cos^2 = ||Pi_N r||^2 / ||r||^2.
"""
import sys, os, json, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'src'))
import mpmath as mp
from zeta_window import build_matrix, F, norm2, EVEN_D
from stage3e_certify import move_block
from weilform import eigsym

# (L, N, dps) exactly as the paper's section 4.3 used them
CELLS = [('0.8', 16, 70), ('1.0', 18, 80), ('1.3', 20, 90), ('1.6', 22, 100), ('2.0', 26, 120)]
SEC = EVEN_D


def K_closed(L, g):
    L, g = mp.mpf(L), mp.mpf(g)
    return (L**3/3 - (L**2*mp.sin(2*g*L)/(2*g) + L*mp.cos(2*g*L)/(2*g**2)
                      - mp.sin(2*g*L)/(4*g**3)))


def reps(L, N, g, sec):
    """representers of c -> F(g), F'(g), F''(g) in the orthonormal basis."""
    nrm = [mp.sqrt(norm2(k, L, sec)) for k in range(N)]
    p  = [F(k, L, g, sec)/nrm[k] for k in range(N)]
    r  = [mp.diff(lambda t, k=k: F(k, L, t, sec), g)/nrm[k] for k in range(N)]
    s  = [mp.diff(lambda t, k=k: F(k, L, t, sec), g, 2)/nrm[k] for k in range(N)]
    return p, r, s


def run(out):
    g1 = None
    rows = []
    for Lx, N, dps in CELLS:
        mp.mp.dps = dps
        if g1 is None or True:
            g1 = mp.im(mp.zetazero(1))
        L = mp.mpf(Lx)
        M = build_matrix(L, N, dps, SEC)
        Q0 = M + move_block(L, range(N), g1, mp.mpf(0), SEC)
        p, r, s = reps(L, N, g1, SEC)
        K = K_closed(L, g1); fourK = 4*K
        nr2 = sum(x**2 for x in r)                 # ||r||^2 should approach K
        vals, vecs = eigsym(Q0)
        V = [mp.matrix(v) for v in vecs]
        # delta^2 coefficient matrix P:  Q_d - Q_0 = delta^2 P + O(delta^4),
        # P = -4( r r^T + sym(p s^T) )
        P = mp.zeros(N, N)
        for a in range(N):
            for b in range(N):
                P[a, b] = -4*(r[a]*r[b] + (p[a]*s[b] + p[b]*s[a])/2)
        lam0 = vals[0]
        cell = dict(L=Lx, N=N, dps=dps, K=mp.nstr(K, 12), fourK=mp.nstr(fourK, 12),
                    norm_r2=mp.nstr(nr2, 12), ratio_nr2_K=float(nr2/K),
                    lam_min_Q0=mp.nstr(lam0, 8),
                    eigs=[mp.nstr(v, 8) for v in vals[:min(N, 12)]], sweep=[])
        for e in range(1, 13):
            d = mp.mpf(10)**(-e)
            Qd = M + move_block(L, range(N), g1, d, SEC)
            lam = eigsym(Qd)[0][0]
            C = (lam0 - lam)/d**2
            thr = fourK*d**2
            idx = [k for k in range(N) if vals[k] <= thr]
            if idx:
                proj2 = sum(sum(V[k][i]*r[i] for i in range(N))**2 for k in idx)
            else:
                proj2 = mp.mpf(0)
            cos2 = proj2/nr2
            pred = fourK*cos2
            # compressed prediction: lam_min of (Lambda_N + d^2 Pi_N P Pi_N) on N
            if idx:
                m = len(idx); G = mp.zeros(m, m)
                for ii, k in enumerate(idx):
                    for jj, l in enumerate(idx):
                        G[ii, jj] = (d**2)*sum(V[k][a]*P[a, b]*V[l][b]
                                               for a in range(N) for b in range(N))
                        if ii == jj: G[ii, jj] += vals[k]
                lam_c = eigsym(G)[0][0]
                C_comp = (lam0 - lam_c)/d**2
            else:
                C_comp = mp.mpf(0)
            cell['sweep'].append(dict(
                delta=mp.nstr(d, 4), lam_min=mp.nstr(lam, 10), C=mp.nstr(C, 10),
                C_float=float(C), dimN=len(idx), cos2=float(cos2),
                pred=mp.nstr(pred, 10), pred_float=float(pred),
                C_compressed=float(C_comp),
                rel_H1=float(abs(C-pred)/abs(C)) if C != 0 else None,
                rel_comp=float(abs(C-C_comp)/abs(C)) if C != 0 else None))
        rows.append(cell)
        print('  L=%-4s N=%-3d 4K=%-12s ||r||^2/K=%.6f  lam_min(Q0)=%s'
              % (Lx, N, mp.nstr(fourK, 8), cell['ratio_nr2_K'], cell['lam_min_Q0']), flush=True)
    json.dump(rows, open(out, 'w'), indent=1)
    return rows


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--out', default='data/stage2_zeta.json')
    a = ap.parse_args()
    rows = run(a.out)
    print('\n  %-5s %-9s %-6s %-13s %-13s %-9s %-11s %-11s' %
          ('L', 'delta', 'dimN', 'C measured', 'H1 pred', 'cos^2', 'rel H1', 'rel compr'))
    for c in rows:
        for s in c['sweep']:
            print('  %-5s %-9s %-6d %-13.6g %-13.6g %-9.5f %-11.3e %-11.3e' %
                  (c['L'], s['delta'], s['dimN'], s['C_float'], s['pred_float'],
                   s['cos2'], s['rel_H1'] or 0, s['rel_comp'] or 0))
