#!/usr/bin/env python3
"""Task 3b (C/CS across L) and Task 4b (inertia of A as N grows, both even bases)."""
import sys, os, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from delta2_law import perturbed_C
from zeta_window import build_parts, EVEN, EVEN_D
from weilform import inertia

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
g1 = None


def task3b():
    global g1
    mp.mp.dps = 30
    g1 = mp.im(mp.zetazero(1))
    rows = []
    print("TASK 3b -- measured C(L) against the Cauchy-Schwarz bound 4L^3/3\n")
    print(f"{'L':>5s} {'N':>4s} {'delta':>7s} {'C measured':>14s} "
          f"{'4(Fp^2+FFpp) pert':>19s} {'CS bound':>12s} {'C/CS':>9s} {'C/L^3':>9s}")
    for Lx, N, dps in [('0.6', 14, 60), ('0.8', 16, 70), ('1.0', 18, 80),
                       ('1.3', 20, 90), ('1.6', 22, 100), ('2.0', 26, 120)]:
        for de in ['0.02']:
            r = perturbed_C(Lx, N, dps, g1, mp.mpf(de))
            B = 4 * mp.mpf(Lx) ** 3 / 3
            print(f"{Lx:>5s} {N:>4d} {de:>7s} {mp.nstr(r['C_measured'],8):>14s} "
                  f"{mp.nstr(r['C_at_perturbed_min'],8):>19s} {mp.nstr(B,6):>12s} "
                  f"{mp.nstr(r['C_measured']/B,5):>9s} "
                  f"{mp.nstr(r['C_measured']/mp.mpf(Lx)**3,5):>9s}")
            rows.append(dict(L=Lx, N=N, dps=dps, delta=de,
                             C_measured=mp.nstr(r['C_measured'], 10),
                             C_at_perturbed_min=mp.nstr(r['C_at_perturbed_min'], 10),
                             C_at_unperturbed_gs=mp.nstr(r['C_at_unperturbed_gs'], 10),
                             cs_bound_4L3_3=mp.nstr(B, 10),
                             ratio_C_over_CS=mp.nstr(r['C_measured'] / B, 8),
                             C_over_L3=mp.nstr(r['C_measured'] / mp.mpf(Lx) ** 3, 8)))
    with open(os.path.join(ROOT, 'data', 'stage3_C_vs_CS.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print("\nwrote data/stage3_C_vs_CS.csv")


def task4b():
    rows = []
    print("\n\nTASK 4b -- inertia of A = Pole + Arch - logpi as N grows, both even bases")
    print("(expectation: A gains positive directions as higher frequencies enter, since")
    print(" the archimedean symbol Re psi(1/4+it/2) ~ log(t/2) grows)\n")
    print(f"{'L':>8s} {'basis':>8s} {'N':>4s} {'inertia A':>14s} {'inertia M':>14s} "
          f"{'inertia S':>14s} {'A neg frac':>11s}")
    for Lx in ['0.5', '0.8', '1.2824746787']:
        for sector in (EVEN_D, EVEN):
            for N in (6, 10, 14, 18, 24):
                mp.mp.dps = 40
                Pole, Arch, LogPi, Prime, _ = build_parts(mp.mpf(Lx), N, 40, sector)
                A = Pole + Arch - LogPi
                Mr = Prime / 2
                S = A - 2 * Mr
                iA, _, _ = inertia(A); iM, _, _ = inertia(Mr); iS, _, _ = inertia(S)
                print(f"{Lx:>8s} {sector:>8s} {N:>4d} {str(iA):>14s} {str(iM):>14s} "
                      f"{str(iS):>14s} {mp.nstr(mp.mpf(iA[0])/N,4):>11s}")
                rows.append(dict(L=Lx, sector=sector, N=N, dps=40,
                                 A_inertia=str(iA), M_inertia=str(iM), S_inertia=str(iS),
                                 A_neg=iA[0], A_pos=iA[2], M_neg=iM[0], M_pos=iM[2],
                                 S_psd=(iS[0] == 0),
                                 A_neg_fraction=mp.nstr(mp.mpf(iA[0]) / N, 6)))
    with open(os.path.join(ROOT, 'data', 'stage3d_inertia_vs_N.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print("\nwrote data/stage3d_inertia_vs_N.csv")


if __name__ == '__main__':
    task3b()
    task4b()
