#!/usr/bin/env python3
"""Stage 3 gate.

(A) Validate the assembled geometric-side matrix against the ZEROS side:
    under RH the explicit formula says  Q(f) = sum_rho |F(gamma_rho)|^2 = 2 sum_{gamma>0} F(gamma)^2.
    This checks the whole assembly (archimedean term included) independently of Zhu.
(B) Reproduce lambda*(0.8) and check it against Zhu's certified interval
    [8.9e-18, 2.27e-17], with basis-size convergence and variational monotonicity.
"""
import sys, os, csv, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from zeta_window import F, build_matrix
from weilform import eigsym, parity_of

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def zeros_side(L, N, gammas):
    """Z_jk = 2 sum_{gamma>0} F_j(gamma) F_k(gamma)  over the supplied ordinates."""
    Fv = [[F(k, L, g) for k in range(N)] for g in gammas]
    Z = mp.zeros(N, N)
    for row in Fv:
        for j in range(N):
            for k in range(j, N):
                Z[j, k] += 2 * row[j] * row[k]
    for j in range(N):
        for k in range(j + 1, N):
            Z[k, j] = Z[j, k]
    return Z


def tail_estimate(L, N, gmax):
    """Riemann-von Mangoldt density dN/dgamma = log(gamma/2pi)/2pi, with
    F_j F_k <= 4 a_j a_k / (gamma^2)^2 asymptotically; crude upper bound on the
    neglected tail of the zeros side."""
    from zeta_window import a_k
    out = mp.mpf(0)
    for j in range(N):
        aj = a_k(j, L)
        out = max(out, 2 * 4 * aj ** 2 / (2 * mp.pi) *
                  mp.quad(lambda g: mp.log(g / (2 * mp.pi)) / g ** 4, [gmax, mp.inf]))
    return out


def gate_A(L, N, Ks, dps):
    mp.mp.dps = dps
    L = mp.mpf(L)
    print(f"\n=== GATE A: geometric side vs zeros side, L={mp.nstr(L,6)}, N={N}, dps={dps} ===")
    M = build_matrix(L, N, dps, verbose=True)
    rows = []
    allg = [mp.im(mp.zetazero(n)) for n in range(1, max(Ks) + 1)]
    for K in Ks:
        Z = zeros_side(L, N, allg[:K]) / L        # same 1/L normalisation as build_matrix
        d = max(abs(M[j, k] - Z[j, k]) for j in range(N) for k in range(N))
        tail = tail_estimate(L, N, allg[K - 1]) / L
        print(f"  K={K:5d} zeros (gamma_K={mp.nstr(allg[K-1],8):>10s})  "
              f"max |geometric - zeros| = {mp.nstr(d,6):>12s}   "
              f"predicted tail <= {mp.nstr(tail,4)}")
        rows.append(dict(L=mp.nstr(L, 6), N=N, K=K, gamma_K=mp.nstr(allg[K - 1], 10),
                         max_abs_diff=mp.nstr(d, 8), tail_bound=mp.nstr(tail, 8)))
    return M, rows


def gate_B(L, Ns, dps_list):
    print(f"\n=== GATE B: lambda_min at L={L} vs Zhu's certified [8.9e-18, 2.27e-17] ===")
    rows, prev = [], None
    for dps in dps_list:
        print(f"  -- working precision {dps} digits --")
        prev = None
        for N in Ns:
            mp.mp.dps = dps
            M = build_matrix(mp.mpf(L), N, dps)
            vals, vecs = eigsym(M)
            lam = vals[0]
            gap = vals[1] - vals[0] if N > 1 else mp.inf
            par, pres = parity_of(vecs[0])
            mono = '' if prev is None else ('OK' if lam <= prev + mp.mpf(10) ** (-(dps - 8))
                                            else 'VIOLATED')
            inzhu = (mp.mpf('8.9e-18') <= lam <= mp.mpf('2.27e-17'))
            print(f"     N={N:3d}  lam_min={mp.nstr(lam,10):>16s}  gap={mp.nstr(gap,6):>10s}"
                  f"  sign={'+' if lam > 0 else '-'}  in Zhu interval: {str(inzhu):5s}"
                  f"  monotone:{mono}")
            rows.append(dict(L=L, N=N, dps=dps, lam_min=mp.nstr(lam, 12),
                             lam_min_sign=('pos' if lam > 0 else 'neg'),
                             gap=mp.nstr(gap, 8), basis_parity=par,
                             in_zhu_interval=inzhu, monotone=mono))
            prev = lam
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', default='0.8')
    ap.add_argument('--gateA-N', type=int, default=4)
    ap.add_argument('--gateA-K', default='50,200,800,2000')
    ap.add_argument('--N', default='2,4,6,8,10,12,14,16')
    ap.add_argument('--dps', default='40,60')
    a = ap.parse_args()
    Ks = [int(x) for x in a.gateA_K.split(',')]
    Ns = [int(x) for x in a.N.split(',')]
    dpss = [int(x) for x in a.dps.split(',')]
    _, rowsA = gate_A(mp.mpf(a.L), a.gateA_N, Ks, 30)
    rowsB = gate_B(a.L, Ns, dpss)
    os.makedirs(os.path.join(ROOT, 'results'), exist_ok=True)
    for name, data in (('stage3_gateA_explicit_formula.csv', rowsA),
                       ('stage3_gateB_lambda_min.csv', rowsB)):
        with open(os.path.join(ROOT, 'results', name), 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(data[0].keys()))
            w.writeheader(); w.writerows(data)
        print(f"wrote results/{name}")


if __name__ == '__main__':
    main()
