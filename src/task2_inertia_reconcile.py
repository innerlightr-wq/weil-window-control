#!/usr/bin/env python3
r"""Task 2: reconcile the reported inertia of A.

Round 1 (Stage 3d) tabulated columns headed "even / odd"  = (EVEN_D, ODD).
Round 2 (Task 4b) tabulated columns headed "both bases"    = (EVEN_D, EVEN).
Those are DIFFERENT second columns, which is the whole of the apparent discrepancy.
This script writes out all three bases and both bookkeeping conventions explicitly so the
question is settled by data rather than by reading two tables with different headers.

Conventions (both give the same S = A - 2M = the Weil form):
  I  (used throughout):  A = Pole + Arch - log(pi)*I,   2M = Prime
  II (alternative)    :  A = Pole + Arch,               2M = Prime + log(pi)*I
"""
import sys, os, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from zeta_window import build_parts, EVEN, EVEN_D, ODD
from weilform import inertia

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rows = []
Ls = ['0.5', '0.8', '1.0', mp.nstr(mp.log(13) / 2, 12)]
print("Task 2: inertia of A, all three bases, both conventions, N = 12, 40 digits\n")
print(f"{'L':>8s} {'basis':>8s} {'conv':>5s} {'inertia A':>13s} {'inertia M':>13s} "
      f"{'inertia S':>13s} {'A neg':>6s} {'split resid':>12s}")
for Lx in Ls:
    for sector in (EVEN_D, EVEN, ODD):
        mp.mp.dps = 40
        N = 12
        Pole, Arch, LogPi, Prime, _ = build_parts(mp.mpf(Lx), N, 40, sector)
        for conv, A, Mr in (('I', Pole + Arch - LogPi, Prime / 2),
                            ('II', Pole + Arch, (Prime + LogPi) / 2)):
            S = A - 2 * Mr
            iA, _, _ = inertia(A); iM, _, _ = inertia(Mr); iS, _, _ = inertia(S)
            resid = max(abs((A - 2 * Mr - S)[i, j]) for i in range(N) for j in range(N))
            print(f"{Lx:>8s} {sector:>8s} {conv:>5s} {str(iA):>13s} {str(iM):>13s} "
                  f"{str(iS):>13s} {iA[0]:>6d} {mp.nstr(resid,4):>12s}")
            rows.append(dict(L=Lx, sector=sector, convention=conv, N=N, dps=40,
                             A_inertia=str(iA), M_inertia=str(iM), S_inertia=str(iS),
                             A_neg=iA[0], M_neg=iM[0], S_psd=(iS[0] == 0),
                             split_residual=mp.nstr(resid, 6)))
    print()
with open(os.path.join(ROOT, 'results', 'stage3d_inertia_reconciled.csv'), 'w',
          newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print("wrote results/stage3d_inertia_reconciled.csv")
