#!/usr/bin/env python3
r"""Gate A, third sign-flip control: the PRIME term.

phase2_gateA.py flips the conductor term and the archimedean term separately.  This script
adds the remaining one -- Prime -> -Prime, i.e. chi(n) -> -chi(n) -- and also the control
that matters most for the CHARACTER itself: replacing chi by the OTHER real primitive
character of a conductor that has two (q = 8: chi_8 vs chi_-8), and by the principal
character.  If Gate A does not blow up under each substitution, it does not discriminate
which character is being tested, and the identification would rest on the self-test alone.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp                                                      # noqa: E402
import dirichlet_window as dw                                            # noqa: E402
from zeta_window import F, norm2, EVEN_D                                 # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def zeros_side(L, sector, N, gammas, K):
    L = mp.mpf(L)
    nrm = [mp.sqrt(norm2(k, L, sector)) for k in range(N)]
    Z = mp.zeros(N, N)
    for g in gammas[:K]:
        Fv = [F(k, L, g, sector) / nrm[k] for k in range(N)]
        for j in range(N):
            for k in range(N):
                Z[j, k] += 2 * Fv[j] * Fv[k]
    return Z


def md(M, Z, N):
    return max(abs(M[j, k] - Z[j, k]) for j in range(N) for k in range(N))


def run(Lx):
    N, K, dps = 4, 160, 40
    mp.mp.dps = dps          # set BEFORE constructing L, so L is exact to `dps`
    L = mp.mpf(Lx)
    prim = [n for n, _, _ in dw.von_mangoldt_terms(2 * L)]
    print("=" * 78)
    print(f"L = {mp.nstr(L, 6)}, N = {N}, K = {K} zeros, sector = even_d, dps = {dps}")
    print(f"prime powers inside the window (log n < 2L): {prim}")
    print(f"chi(n) on them: " + ", ".join(
        f"{dw.CHARS[d]['label']}={[dw.chi(d, n) for n in prim]}"
        for d in dw.DISCRIMINANTS))
    print("=" * 78)
    print(f"  {'chi':8s} {'baseline':>12s} {'flip prime':>14s} {'ratio':>11s} "
          f"{'wrong chi (a-swap)':>20s} {'ratio':>11s}")
    # the a-swap: use the same q but the OTHER parity's archimedean shift, keeping
    # everything else -- isolates whether Gate A sees the character's PARITY
    for d in dw.DISCRIMINANTS:
        lab = dw.CHARS[d]['label']
        gs = [mp.mpf(x) for x in open(os.path.join(HERE, 'data', f'zeros_{lab}.txt'))
              if not x.startswith('#')]
        Z = zeros_side(L, EVEN_D, N, gs, min(K, len(gs)))
        base = md(dw.build_matrix_chi(L, range(N), dps, d, EVEN_D), Z, N)
        flip = md(dw.build_matrix_chi(L, range(N), dps, d, EVEN_D, flip='prime'), Z, N)
        # parity swap: temporarily flip a
        orig = dw.CHARS[d]['a']
        dw.CHARS[d]['a'] = 1 - orig
        swap = md(dw.build_matrix_chi(L, range(N), dps, d, EVEN_D), Z, N)
        dw.CHARS[d]['a'] = orig
        print(f"  {lab:8s} {mp.nstr(base, 6):>12s} {mp.nstr(flip, 6):>14s} "
              f"{mp.nstr(flip / base, 5):>11s} {mp.nstr(swap, 6):>20s} "
              f"{mp.nstr(swap / base, 5):>11s}")
    print("\n  Explicit chi-substitution control (prime weights of a DIFFERENT character,")
    print("  conductor term and archimedean term left correct):")
    print(f"  {'zeros of':8s} {'primes of':10s} {'discrepancy':>14s} {'ratio to correct':>18s}")
    for d in dw.DISCRIMINANTS:
        lab = dw.CHARS[d]['label']
        gs = [mp.mpf(x) for x in open(os.path.join(HERE, 'data', f'zeros_{lab}.txt'))
              if not x.startswith('#')]
        Z = zeros_side(L, EVEN_D, N, gs, min(K, len(gs)))
        base = md(dw.build_matrix_chi(L, range(N), dps, d, EVEN_D), Z, N)
        for d2 in dw.DISCRIMINANTS:
            if d2 == d:
                continue
            Arch, LogPi, LogQ, Prime, _ = dw.build_parts_chi(L, range(N), dps, d, EVEN_D)
            _, _, _, Prime2, _ = dw.build_parts_chi(L, range(N), dps, d2, EVEN_D)
            M = Arch - LogPi + LogQ - Prime2
            v = md(M, Z, N)
            print(f"  {lab:8s} {dw.CHARS[d2]['label']:10s} {mp.nstr(v, 6):>14s} "
                  f"{mp.nstr(v / base, 5):>18s}")
        print()


def main():
    for Lx in ('0.8', '1.0', '1.3'):
        run(Lx)
        print()


if __name__ == '__main__':
    main()
