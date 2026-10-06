#!/usr/bin/env python3
r"""Final-revision checks: the exact constant K(L,gamma), and the function-field <-> zeta
dictionary."""
import sys, os, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def K_closed(L, g):
    """K(L,gamma) = int_{-L}^{L} u^2 sin^2(gamma u) du, closed form."""
    L, g = mp.mpf(L), mp.mpf(g)
    return (L ** 3 / 3
            - (L ** 2 * mp.sin(2 * g * L) / (2 * g)
               + L * mp.cos(2 * g * L) / (2 * g ** 2)
               - mp.sin(2 * g * L) / (4 * g ** 3)))


def main():
    mp.mp.dps = 40
    g1 = mp.im(mp.zetazero(1))
    print("=" * 76)
    print("(2a) K(L,gamma) closed form vs quadrature, and vs the uniform bound 2L^3/3")
    print("=" * 76)
    print(f"{'L':>5s} {'K closed':>16s} {'K quad':>16s} {'diff':>10s} "
          f"{'L^3/3':>12s} {'K/(L^3/3)':>11s} {'2L^3/3':>12s} {'K<=2L^3/3':>10s}")
    rows = []
    for Lx in ['0.6', '0.8', '1.0', '1.3', '1.6', '2.0']:
        L = mp.mpf(Lx)
        kc = K_closed(L, g1)
        kq = mp.quad(lambda u: u ** 2 * mp.sin(g1 * u) ** 2, [-L, 0, L])
        asym = L ** 3 / 3
        print(f"{Lx:>5s} {mp.nstr(kc,12):>16s} {mp.nstr(kq,12):>16s} "
              f"{mp.nstr(abs(kc-kq),4):>10s} {mp.nstr(asym,8):>12s} "
              f"{mp.nstr(kc/asym,8):>11s} {mp.nstr(2*L**3/3,8):>12s} "
              f"{str(kc <= 2*L**3/3):>10s}")
        rows.append((Lx, kc, asym))
    print("\n  correction |K/(L^3/3) - 1|, i.e. how much the sin^2 weight matters at gamma_1:")
    for Lx, kc, asym in rows:
        print(f"    L={Lx:>4s}: {mp.nstr(abs(kc/asym-1)*100, 6)} %")

    print("\n" + "=" * 76)
    print("(2b) the sup is ATTAINED by f proportional to u sin(gamma u)")
    print("=" * 76)
    for Lx in ['0.8', '2.0']:
        L = mp.mpf(Lx)
        nrm2 = mp.quad(lambda u: (u * mp.sin(g1 * u)) ** 2, [-L, 0, L])
        # F'(gamma) = -int u f(u) sin(gamma u) du  with f = u sin(gamma u)/||.||
        Fp = -mp.quad(lambda u: u * (u * mp.sin(g1 * u) / mp.sqrt(nrm2)) * mp.sin(g1 * u),
                      [-L, 0, L])
        print(f"    L={Lx}: |F'(gamma_1)|^2 = {mp.nstr(Fp**2,14)}   K(L,gamma_1) = "
              f"{mp.nstr(K_closed(L,g1),14)}   equal: "
              f"{abs(Fp**2 - K_closed(L,g1)) < mp.mpf(10)**-25}")

    print("\n" + "=" * 76)
    print("(2c) recomputed table: C(L) against 4K(L,gamma_1) and against 8L^3/3")
    print("=" * 76)
    cs = {r['L']: r for r in csv.DictReader(open(os.path.join(ROOT, 'data',
                                                             'stage3_C_vs_CS.csv')))}
    out = []
    print(f"{'L':>5s} {'measured C':>13s} {'4K(L,g1)':>13s} {'C/4K':>9s} "
          f"{'8L^3/3':>11s} {'C/(8L^3/3)':>11s} {'4L^3/3':>11s} {'C/(4L^3/3)':>11s}")
    for Lx in ['0.8', '1.0', '1.3', '1.6', '2.0']:
        L = mp.mpf(Lx)
        C = mp.mpf(cs[Lx]['C_measured'])
        fourK = 4 * K_closed(L, g1)
        print(f"{Lx:>5s} {mp.nstr(C,8):>13s} {mp.nstr(fourK,8):>13s} "
              f"{mp.nstr(C/fourK,5):>9s} {mp.nstr(8*L**3/3,7):>11s} "
              f"{mp.nstr(C/(8*L**3/3),5):>11s} {mp.nstr(4*L**3/3,7):>11s} "
              f"{mp.nstr(C/(4*L**3/3),5):>11s}")
        out.append(dict(L=Lx, C_measured=mp.nstr(C, 10),
                        K_L_gamma1=mp.nstr(K_closed(L, g1), 10),
                        four_K=mp.nstr(fourK, 10), ratio_C_over_4K=mp.nstr(C / fourK, 8),
                        proved_8L3_3=mp.nstr(8 * L ** 3 / 3, 10),
                        ratio_C_over_proved=mp.nstr(C / (8 * L ** 3 / 3), 8),
                        asymptotic_4L3_3=mp.nstr(4 * L ** 3 / 3, 10),
                        ratio_C_over_asymptotic=mp.nstr(C / (4 * L ** 3 / 3), 8),
                        sin2_correction_pct=mp.nstr(abs(K_closed(L, g1) /
                                                        (L ** 3 / 3) - 1) * 100, 6)))
    with open(os.path.join(ROOT, 'data', 'stage3_exact_constant_K.csv'), 'w',
              newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader(); w.writerows(out)
    print("  wrote data/stage3_exact_constant_K.csv")

    print("\n" + "=" * 76)
    print("(6) function-field <-> zeta dictionary, checked numerically on one curve")
    print("=" * 76)
    print("  u = n log q,  R log q = 2L,  a = delta log q,  Rayleigh quotients differ")
    print("  by a factor log q (discrete sum vs integral, du = log q dn).")
    print("  Prediction: lam_min^{ff}(R,a) = lam_min^{zeta}(L,delta) / log q.\n")
    q = 5                      # the curve y^2 = x^3 + x / F_5 lives over q = 5
    lq = mp.log(q)
    print(f"  {'R':>4s} {'L=R log q/2':>13s} {'delta=a/log q':>15s} "
          f"{'ff lam_min':>18s} {'zeta pred / log q':>19s} {'ratio':>9s}")
    for R in (8, 16, 32, 64, 128):
        a = mp.mpf('1e-6')
        rho = mp.e ** a
        u2 = (rho ** (2 * R + 2) - 1) / (rho ** 2 - 1)
        w2 = (rho ** (-2 * R - 2) - 1) / (rho ** -2 - 1)
        ff = (R + 1) - mp.sqrt(u2 * w2)
        L = R * lq / 2
        delta = a / lq
        zeta_pred = -(4 * L ** 3 / 3) * delta ** 2 / lq
        print(f"  {R:>4d} {mp.nstr(L,8):>13s} {mp.nstr(delta,8):>15s} "
              f"{mp.nstr(ff,10):>18s} {mp.nstr(zeta_pred,10):>19s} "
              f"{mp.nstr(ff/zeta_pred,8):>9s}")
    print("\n  ratio -> 1 as R grows (the gap is the O(R^2) term in binom(R+2,3) vs R^3/6).")
    print("  Factors of 2 that cancel: pair (function field) vs quartet (zeta) contributes")
    print("  a 2; the sin^2 averaging in K -> L^3/3 contributes a 1/2.")


if __name__ == '__main__':
    main()
