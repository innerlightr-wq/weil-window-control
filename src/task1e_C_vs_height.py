#!/usr/bin/env python3
"""Task 1(e): how does C(L) in lam_min ~ -C delta^2 depend on the planted HEIGHT,
including heights above the window bandwidth T*(L) = 2 pi e^{2L}?"""
import sys, os, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from zeta_window import build_matrix_idx, EVEN_D
from stage3e_certify import move_block, index_set, L_pred
from weilform import eigsym

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
mp.mp.dps = 30
NS = [1, 5, 10, 20, 30, 50]
G = {n: mp.im(mp.zetazero(n)) for n in NS}
rows = []
print("Task 1(e): C(L, gamma_*) = -lam_min/delta^2 vs planted height")
print("T*(L) = 2 pi e^{2L} is the window bandwidth; heights above it should not be seen.\n")
for Lx, dps in [('1.0', 140), ('1.4', 180), ('1.8', 220)]:
    L = mp.mpf(Lx)
    Tstar = 2 * mp.pi * mp.e ** (2 * L)
    print(f"L = {Lx}   T*(L) = {mp.nstr(Tstar, 8)}   (4L^3/3 = {mp.nstr(4*L**3/3,6)})")
    print(f"  {'n':>4s} {'gamma_n':>12s} {'g/T*':>8s} {'lam*(L)':>14s} "
          f"{'lam_min pert':>15s} {'C = -lam/d^2':>14s} {'C/(4L^3/3)':>11s}")
    for n in NS:
        gs = G[n]
        de = mp.mpf('0.02')
        idx, kstar = index_set(L, gs, 24, 6)
        mp.mp.dps = dps
        M = build_matrix_idx(L, idx, dps, EVEN_D)
        ctrl = eigsym(M)[0][0]
        Z = M + move_block(L, idx, gs, de)
        lam = eigsym(Z)[0][0]
        C = -lam / de ** 2
        print(f"  {n:>4d} {mp.nstr(gs,8):>12s} {mp.nstr(gs/Tstar,4):>8s} "
              f"{mp.nstr(ctrl,6):>14s} {mp.nstr(lam,6):>15s} {mp.nstr(C,6):>14s} "
              f"{mp.nstr(C/(4*L**3/3),5):>11s}")
        rows.append(dict(L=Lx, dps=dps, n=n, gamma_n=mp.nstr(gs, 12),
                         Tstar=mp.nstr(Tstar, 10),
                         gamma_over_Tstar=mp.nstr(gs / Tstar, 8),
                         dim=len(idx), kstar=kstar, delta=str(de),
                         lam_star=mp.nstr(ctrl, 10), lam_min_perturbed=mp.nstr(lam, 10),
                         C=mp.nstr(C, 10), C_over_CS=mp.nstr(C / (4 * L ** 3 / 3), 8),
                         detected=('yes' if lam < 0 else 'no')))
    print()
with open(os.path.join(ROOT, 'data', 'stage3e_C_vs_height.csv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print("wrote data/stage3e_C_vs_height.csv")
