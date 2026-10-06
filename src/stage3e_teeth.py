#!/usr/bin/env python3
"""Stage 3e: teeth on the zeta side (a DIAGNOSTIC, not a certificate).

This assembles the window form from ZERO data, not from primes:

    Q(f) = sum_rho h(gamma_rho),   h(t) = F(t) F(-t),

using the first K actual zeta zeros PLUS one planted off-line quartet
{rho, 1-rho, conj(rho), 1-conj(rho)} with Re rho = 0.7.

Writing rho = 1/2 + i z, the quartet's ordinates are  {+-gamma +- i delta}, delta = 0.2
(exactly the Stage 2 off-line quartet, with e^{delta} playing the role of the radius rho
and L the role of the window R).  Its contribution to the form is

    4 Re[ F(w) F(-w) ],   w = gamma_* + i delta,

which for even f is 4 Re[F(w)^2].  Since |F(gamma + i delta)| ~ e^{delta L} |F|, the
planted term grows like e^{2 delta L} while the on-line terms stay bounded: the form must
turn indefinite once L is large enough.

BECAUSE THIS USES ZERO DATA IT PROVES NOTHING ABOUT RH.  It measures how large a window
is needed before an off-line zero at a given height and depth shows up in the sign of
lambda_min -- the zeta-side analogue of the Stage 2 detection law.
"""
import sys, os, csv, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from zeta_window import F, norm2, EVEN_D
from weilform import eigsym, poly_zeros

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def zeros_matrix(L, N, gammas, planted, sector=EVEN_D):
    """planted: list of (gamma_star, delta). Returns the (N x N) form."""
    L = mp.mpf(L)
    nrm = [mp.sqrt(norm2(k, L, sector)) for k in range(N)]
    Z = mp.zeros(N, N)
    for g in gammas:
        Fv = [F(k, L, g, sector) / nrm[k] for k in range(N)]
        for j in range(N):
            for k in range(j, N):
                Z[j, k] += 2 * Fv[j] * Fv[k]
    for gs, de in planted:
        w = mp.mpc(gs, de)
        Fw = [F(k, L, w, sector) / nrm[k] for k in range(N)]
        Fmw = [F(k, L, -w, sector) / nrm[k] for k in range(N)]
        for j in range(N):
            for k in range(j, N):
                Z[j, k] += 2 * mp.re(Fw[j] * Fmw[k] + Fw[k] * Fmw[j])
    for j in range(N):
        for k in range(j + 1, N):
            Z[k, j] = Z[j, k]
    return Z


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--K', type=int, default=400)
    ap.add_argument('--N', type=int, default=12)
    ap.add_argument('--dps', type=int, default=40)
    ap.add_argument('--gamma-star', default='14.134725141734693')
    ap.add_argument('--deltas', default='0.2,0.1,0.05,0.02')
    ap.add_argument('--Ls', default='0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.3,1.6,2.0')
    ap.add_argument('--Nscale', type=float, default=10.0,
                    help='N = max(--N, ceil(Nscale*L)); the basis must grow with L or '
                         'it under-resolves (its top frequency is N*pi/2L)')
    a = ap.parse_args()
    mp.mp.dps = a.dps
    print(f"Stage 3e: zeros-side diagnostic. K={a.K} true zeros + one planted quartet")
    print(f"  Re rho = 0.7 <-> delta = 0.2;  planted height gamma_* = {a.gamma_star}")
    print("  NOTE: built from ZERO data, so it is a diagnostic, NOT a certificate.\n")
    g = [mp.im(mp.zetazero(n)) for n in range(1, a.K + 1)]
    Ls = [mp.mpf(x) for x in a.Ls.split(',')]
    deltas = [mp.mpf(x) for x in a.deltas.split(',')]
    gs = mp.mpf(a.gamma_star)
    rows = []
    print(f"  {'L':>6s} {'control (no plant)':>20s}" +
          ''.join(f"{'delta='+str(d):>22s}" for d in deltas))
    for L in Ls:
        NL = max(a.N, int(mp.ceil(a.Nscale * L)))
        ctrl = eigsym(zeros_matrix(L, NL, g, []))[0][0]
        # the control is a sum of rank-1 PSD terms, so it is EXACTLY PSD; any negative
        # value is pure roundoff and measures the noise floor at this precision.
        floor = abs(ctrl) if ctrl < 0 else None
        line = f"  {mp.nstr(L,5):>6s} {mp.nstr(ctrl,8):>20s}"
        rec = dict(L=mp.nstr(L, 6), N=NL, K=a.K, dps=a.dps,
                   control_is_negative_roundoff=(ctrl < 0),
                   gamma_star=a.gamma_star, lam_min_control=mp.nstr(ctrl, 10))
        for d in deltas:
            lam = eigsym(zeros_matrix(L, NL, g, [(gs, d)]))[0][0]
            line += f"{mp.nstr(lam,8):>22s}"
            # only call it negative if it clears the measured noise floor
            sig = 'neg' if (lam < 0 and (floor is None or abs(lam) > 100 * floor)) else \
                  ('pos' if lam > 0 else 'below-noise')
            rec[f'lam_min_delta_{mp.nstr(d,4)}'] = mp.nstr(lam, 10)
            rec[f'sign_delta_{mp.nstr(d,4)}'] = sig
        print(line)
        rows.append(rec)
    with open(os.path.join(ROOT, 'results', 'stage3e_teeth.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print("\n  first L at which lambda_min < 0:")
    for d in deltas:
        key = f'sign_delta_{mp.nstr(d,4)}'
        first = next((r['L'] for r in rows if r[key] == 'neg'), None)
        print(f"    delta={mp.nstr(d,4):>8s} (Re rho = {mp.nstr(mp.mpf(1)/2+d,6)}) -> L = {first}")
    neg = [r['L'] for r in rows if r['control_is_negative_roundoff']]
    print(f"\n  control (no planted zero) is PSD by construction; rows where roundoff "
          f"made it negative (= noise floor reached): {neg or 'none'}")
    print("  wrote results/stage3e_teeth.csv")


if __name__ == '__main__':
    main()
