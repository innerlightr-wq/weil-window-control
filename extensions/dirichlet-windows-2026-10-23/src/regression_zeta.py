#!/usr/bin/env python3
r"""The zeta regression suite.

The brief's working rule: "Keep known-good zeta regression values in the suite so a sector
or sign bug in the chi code is caught by the zeta control."  This script is that suite.  It
recomputes zeta-side quantities THROUGH THE SHARED CODE PATH and asserts them against the
PARENT REPOSITORY'S COMMITTED CSVs -- not against values recorded by this extension.

  (R1) Gate B lambda_min(L=0.8) at N = 16, 20, 24 in BOTH even bases
       vs  data/stage3_gateB_lambda_min.csv
  (R2) the even/odd split at L = 0.5, 0.8, 1.0, N = 12, dps = 40 -- INCLUDING the ODD
       sector, which is the value a sector/sign bug would move
       vs  data/stage3c_even_odd.csv
  (R3) the measured delta^2 coefficient C(L) and C at the perturbed minimiser
       vs  data/stage3_C_vs_CS.csv
  (R4) arch_closed_a(., ., a = 0) == the parent arch_closed, exactly
  (R5) real_block(., ., 1/2, sector) == the parent pole-term matrix, every sector

Exits non-zero if any value drifts.  Tolerances are RELATIVE and tight (1e-20), because the
parent's committed values were produced by the same mpmath at the same precision; a real
disagreement is a bug, not noise.
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp                                                      # noqa: E402
import planting as pl                                                    # noqa: E402
from zeta_window import (arch_closed, build_matrix, Csym_terms, F_at_half,  # noqa: E402
                         norm2, pole_sign, EVEN, EVEN_D, ODD)
from dirichlet_window import arch_closed_a                               # noqa: E402
from weilform import eigsym                                              # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PARENT = os.path.dirname(os.path.dirname(HERE))
TOL = mp.mpf('1e-20')
FAIL = []


def check(name, got, want, tol=TOL):
    want = mp.mpf(want)
    rel = abs(got - want) / abs(want) if want != 0 else abs(got - want)
    ok = rel < tol
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    print(f"         got  {mp.nstr(got, 14)}")
    print(f"         want {mp.nstr(want, 14)}   rel {mp.nstr(rel, 4)}")
    if not ok:
        FAIL.append(name)


def pcsv(name):
    return list(csv.DictReader(open(os.path.join(PARENT, 'data', name))))


def main():
    print("=" * 78)
    print("(R1) Gate B lambda_min at L = 0.8, two bases   [parent stage3_gateB_lambda_min]")
    print("=" * 78)
    rows = pcsv('stage3_gateB_lambda_min.csv')
    for sector in (EVEN_D, EVEN):
        for N in (16, 20, 24):
            want = [r for r in rows if r['sector'] == sector and int(r['N']) == N][0]
            dps = int(want['dps'])
            mp.mp.dps = dps
            lam = eigsym(build_matrix(mp.mpf('0.8'), N, dps, sector))[0][0]
            check(f"L=0.8 {sector} N={N} dps={dps}", lam, want['lam_min'])

    print("=" * 78)
    print("(R2) even/odd split, N = 12, dps = 40   [parent stage3c_even_odd]")
    print("=" * 78)
    rows = pcsv('stage3c_even_odd.csv')
    mp.mp.dps = 40
    for r in rows:
        L = mp.mpf(r['L'])
        N, dps = int(r['N']), int(r['dps'])
        mp.mp.dps = dps
        le = eigsym(build_matrix(L, N, dps, EVEN_D))[0][0]
        lo = eigsym(build_matrix(L, N, dps, ODD))[0][0]
        check(f"L={r['L']} N={N} EVEN_D", le, r['lam_min_even'])
        check(f"L={r['L']} N={N} ODD   (the sector-bug guard)", lo, r['lam_min_odd'])

    print("=" * 78)
    print("(R3) measured delta^2 coefficient C(L)   [parent stage3_C_vs_CS]")
    print("=" * 78)
    rows = pcsv('stage3_C_vs_CS.csv')
    for r in rows:
        if r['L'] not in ('0.8', '1.3'):
            continue
        L = mp.mpf(r['L'])
        N, dps, de = int(r['N']), int(r['dps']), mp.mpf(r['delta'])
        mp.mp.dps = dps
        M = build_matrix(L, N, dps, EVEN_D)
        g1 = mp.im(mp.zetazero(1))
        P = pl.pair_block(L, range(N), g1, EVEN_D)
        Qd = M - P + pl.quartet_block(L, range(N), g1, de, EVEN_D)
        lam = eigsym(Qd)[0][0]
        check(f"L={r['L']} N={N} dps={dps} delta={r['delta']}  C = -lam_min/delta^2",
              -lam / de ** 2, r['C_measured'], mp.mpf('1e-9'))

    print("=" * 78)
    print("(R4) arch_closed_a(a=0) == parent arch_closed, exactly")
    print("=" * 78)
    mp.mp.dps = 50
    worst = mp.mpf(0)
    for sector in (EVEN_D, EVEN, ODD):
        for Lx in ('0.5', '0.8', '1.3'):
            L = mp.mpf(Lx)
            for j in range(4):
                for k in range(j, 4):
                    tm = Csym_terms(j, k, L, sector)
                    worst = max(worst, abs(arch_closed(tm, L) - arch_closed_a(tm, L, 0)))
    print(f"  [{'PASS' if worst == 0 else 'FAIL'}] max abs diff = {mp.nstr(worst, 6)} "
          f"(must be exactly 0)")
    if worst != 0:
        FAIL.append('R4')

    print("=" * 78)
    print("(R5) real_block(delta=1/2) == parent pole-term matrix, every sector")
    print("=" * 78)
    worst = mp.mpf(0)
    for sector in (EVEN_D, EVEN, ODD):
        L = mp.mpf('0.8')
        B, _ = pl.real_block(L, range(5), mp.mpf(1) / 2, sector)
        nrm = [mp.sqrt(norm2(k, L, sector)) for k in range(5)]
        P = [F_at_half(k, L, sector) for k in range(5)]
        sg = pole_sign(sector)
        for i in range(5):
            for j in range(5):
                worst = max(worst, abs(B[i, j] - 2 * sg * P[i] * P[j]
                                       / (nrm[i] * nrm[j])))
    print(f"  [{'PASS' if worst < mp.mpf(10) ** -40 else 'FAIL'}] "
          f"max abs diff = {mp.nstr(worst, 6)}")
    if worst >= mp.mpf(10) ** -40:
        FAIL.append('R5')

    print()
    if FAIL:
        print(f"REGRESSION FAILED: {len(FAIL)}: {FAIL}")
        return 1
    print("ZETA REGRESSION SUITE: ALL PASS "
          "(every value is the parent repository's committed one)")
    return 0


if __name__ == '__main__':
    sys.exit(main())
