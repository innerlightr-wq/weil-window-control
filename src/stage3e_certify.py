#!/usr/bin/env python3
r"""Task 1: certification asymmetry, redone with a well-posed perturbation.

TWO CORRECTIONS to the first Stage 3e height sweep, both forced by this round's analysis.

(A) The on-line side must come from the GEOMETRIC assembly, not a truncated sum over the
    first K zeros.  Near the onset the negative Rayleigh quotient is ~1e-27 while
    truncation omits positive mass ~1e-6, so the truncated form cannot certify anything
    there -- and, being biased downward, it manufactures false detections.

(B) The perturbation must MOVE a zero off the line, not ADD a quartet on top of all the
    true zeros.  Writing w = gamma + i delta, for real even f,

        4 Re F(w)^2  =  4 F(gamma)^2 - 4 delta^2 (F'(gamma)^2 + F(gamma) F''(gamma)) + O(d^4)

    so simply ADDING a quartet contributes +4F(gamma)^2 > 0 unless F nearly vanishes at
    gamma.  The first sweep planted at heights 18, 50, 100, 200 -- none of which is a zeta
    zero -- so the true (untruncated) contribution there is POSITIVE, and the "detections"
    reported at those heights were artefacts of the zeros-sum truncation.

    The well-posed perturbation removes the on-line pair at +-gamma_n and inserts the
    off-line quartet in its place:

        Delta(f) = -2 F(gamma_n)^2 + 4 Re F(gamma_n + i delta)^2
                 =  2 F(gamma_n)^2 - 4 delta^2 (F'^2 + F F'') + O(delta^4),

    which is negative exactly when the trial function nearly annihilates gamma_n -- which
    is what the window minimiser does, for the leading zeros.  So we plant AT actual zeta
    zero ordinates gamma_n and sweep n.

A negative Rayleigh quotient from any trial vector certifies indefiniteness; a basis that
finds none certifies nothing.  So we grow the basis and watch whether the onset moves.
"""
import sys, os, csv, argparse, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from zeta_window import build_matrix_idx, F, norm2, EVEN_D
from ldl import ldl, first_negative_pivot, witness, rayleigh

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASES = [('B1', 8, 0), ('B3', 24, 3), ('B5', 44, 7)]


def L_pred(gs):
    """T*(L) = 2 pi e^{2L} first reaches gamma_*  <=>  L = (1/2) log(gamma_*/2pi)."""
    return mp.log(mp.mpf(gs) / (2 * mp.pi)) / 2


def index_set(L, gs, N0, m):
    kstar = int(mp.nint(mp.mpf(gs) * mp.mpf(L) / mp.pi - mp.mpf(1) / 2))
    res = [k for k in range(max(0, kstar - m), kstar + m + 1)]
    return sorted(set(list(range(N0)) + res)), kstar


def move_block(L, idx, gs, delta, sector=EVEN_D):
    """-2 F_j(g) F_k(g)  +  symmetrised quartet block at w = g + i delta."""
    L = mp.mpf(L)
    n = len(idx)
    nrm = [mp.sqrt(norm2(k, L, sector)) for k in idx]
    g = mp.mpf(gs)
    w = mp.mpc(g, mp.mpf(delta))
    Fr = [F(k, L, g, sector) / nrm[a] for a, k in enumerate(idx)]
    Fw = [F(k, L, w, sector) / nrm[a] for a, k in enumerate(idx)]
    Fm = [F(k, L, -w, sector) / nrm[a] for a, k in enumerate(idx)]
    B = mp.zeros(n, n)
    for a in range(n):
        for b in range(a, n):
            B[a, b] = B[b, a] = (-2 * Fr[a] * Fr[b]
                                 + 2 * mp.re(Fw[a] * Fm[b] + Fw[b] * Fm[a]))
    return B


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ns', default='1,5,10,30,80', help='plant at zeta zero gamma_n')
    ap.add_argument('--deltas', default='0.1,0.02')
    ap.add_argument('--factors', default='0.70,0.85,1.00,1.20,1.50,2.00')
    ap.add_argument('--out', default='stage3e_certified_onset.csv')
    ap.add_argument('--dps', type=int, default=0, help='0 = auto (30 + 40 L)')
    a = ap.parse_args()
    mp.mp.dps = 30
    ns = [int(x) for x in a.ns.split(',')]
    gammas = {n: mp.im(mp.zetazero(n)) for n in ns}
    rows = []
    print("Task 1: detection onset vs basis size, with the WELL-POSED perturbation")
    print("  on-line side = geometric (exact);  planted AT zeta zero ordinates\n")
    for n in ns:
        gs = gammas[n]
        Lp = L_pred(gs)
        print(f"=== plant at gamma_{n} = {mp.nstr(gs,10)}   "
              f"L_pred = {mp.nstr(Lp,6)} ===")
        for de in [mp.mpf(x) for x in a.deltas.split(',')]:
            print(f"  delta = {mp.nstr(de,4)}")
            print(f"    {'L':>7s} {'L/Lpred':>8s} {'basis':>6s} {'dim':>4s} {'k*':>5s} "
                  f"{'lam*(L)':>14s} {'min Rayleigh':>15s} {'verdict':>10s}")
            for fac in [mp.mpf(x) for x in a.factors.split(',')]:
                L = Lp * fac
                for tag, N0, m in BASES:
                    idx, kstar = index_set(L, gs, N0, m)
                    t0 = time.time()
                    dps = a.dps or max(60, int(30 + 40 * float(L)))
                    mp.mp.dps = dps
                    M = build_matrix_idx(L, idx, dps, EVEN_D)
                    ctrl = mp.eigsy(M)[0][0]
                    Z = M + move_block(L, idx, gs, de)
                    try:
                        Lm, D = ldl(Z)
                        j = first_negative_pivot(D)
                    except ZeroDivisionError:
                        j = None
                    rec = dict(n=n, gamma_n=mp.nstr(gs, 12), delta=mp.nstr(de, 6),
                               L=mp.nstr(L, 8), L_pred=mp.nstr(Lp, 8),
                               L_over_Lpred=mp.nstr(fac, 4), basis=tag, N0=N0, mres=m,
                               kstar=kstar, dim=len(idx), dps=dps,
                               lam_min_geom=mp.nstr(ctrl, 10))
                    floor = abs(ctrl) if ctrl < 0 else None
                    rec['noise_floor'] = mp.nstr(floor, 6) if floor else ''
                    if j is None:
                        rec.update(detected='no', min_rayleigh='', witness_support='',
                                   trustworthy='')
                    else:
                        v = witness(Lm, j)
                        r, _, _ = rayleigh(Z, v)
                        # the geometric control is PSD by construction; if it came out
                        # negative, precision is exhausted and |control| is the noise
                        # floor that a detection must clear by a wide margin.
                        trust = (floor is None) or (abs(r) > 100 * floor)
                        rec.update(detected='yes', min_rayleigh=mp.nstr(r, 10),
                                   witness_support=str(j + 1),
                                   trustworthy=('yes' if trust else
                                                'NO-below-noise-floor'))
                        if not trust:
                            j = None
                        print(f"    {mp.nstr(L,5):>7s} {mp.nstr(fac,4):>8s} {tag:>6s} "
                              f"{len(idx):>4d} {kstar:>5d} {mp.nstr(ctrl,6):>14s} "
                              f"{mp.nstr(r,6):>15s} "
                              f"{('DETECTED' if rec['trustworthy']=='yes' else 'below-noise'):>11s} "
                              f"({time.time()-t0:.0f}s)")
                    if j is None:
                        print(f"    {mp.nstr(L,5):>7s} {mp.nstr(fac,4):>8s} {tag:>6s} "
                              f"{len(idx):>4d} {kstar:>5d} {mp.nstr(ctrl,6):>14s} "
                              f"{'-':>15s} {'not found':>10s} ({time.time()-t0:.0f}s)")
                    rows.append(rec)
            print()
    with open(os.path.join(ROOT, 'data', a.out), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"wrote data/{a.out}")


if __name__ == '__main__':
    main()
