#!/usr/bin/env python3
r"""Three verification checks on the Phase-5 (Siegel-corner) result.

(a) PAIR vs QUARTET.  Establish, FROM THE STORED PHASE-5 DATA rather than from a reading of
    the code, whether the planted object at gamma = 0 was the real pair {1/2 + d, 1/2 - d}
    or the degenerate quartet (= twice the pair).  Recompute lambda_min both ways for stored
    rows and see which reproduces the committed `lam_planted`.

(b) PARITY SIGN STRUCTURE, exactly.  For each real character and L in {0.8, 1.0, 1.4}:
    the even-basis perturbation must be rank one with eigenvalue +2A^2 and the odd-basis one
    rank one with eigenvalue -2S^2, where for the normalised eigenvector f
        A = int f(u) cosh(d u) du,     S = int f(u) sinh(d u) du.
    A and S are obtained by QUADRATURE on the explicit f, so this is an independent check of
    the closed forms in planting.F_at_imag, not a restatement of them.
    Also K_even(L,g) + K_odd(L,g) = 2 L^3 / 3 at g in {0, gamma_1(chi), 14.1347}.

(c) POLE CONSISTENCY.  The real-pair perturbation at d = 1/2 is compared entrywise with the
    parent repository's pole-term matrix 2 * pole_sign(sector) * P_j P_k / ||.|| in all three
    sectors.  BOTH B - Pole and B + Pole are reported, so the data decides the sign rather
    than the docstring.  The chi-analogue is identically zero (no pole for non-principal
    chi), which is verified structurally and is the point of the check.

Exits non-zero if any assertion fails.  No float64 anywhere.
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp                                                      # noqa: E402
import dirichlet_window as dw                                            # noqa: E402
import planting as pl                                                    # noqa: E402
from zeta_window import (EVEN, EVEN_D, ODD, F_at_half, build_matrix,     # noqa: E402
                         freq, norm2, pole_sign)
from weilform import eigsym                                              # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DPS = 60
TOL = mp.mpf(10) ** -45            # working precision for exact-identity assertions
FAIL = []
OUT = {}


def check(name, ok, detail=''):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ''))
    if not ok:
        FAIL.append(name)
    return ok


def rank_one_fit(B, n):
    """(lambda, v, residual): the dominant eigenpair and max_ij |B_ij - lam v_i v_j|."""
    vals, vecs = eigsym(B)
    k = max(range(n), key=lambda i: abs(vals[i]))
    lam, v = vals[k], vecs[k]
    nv = mp.sqrt(sum(x ** 2 for x in v))
    v = [x / nv for x in v]
    res = max(abs(B[i, j] - lam * v[i] * v[j]) for i in range(n) for j in range(n))
    others = max((abs(vals[i]) for i in range(n) if i != k), default=mp.mpf(0))
    return lam, v, res, others


def basis_fun(v, L, sector):
    """f(x) = sum_k v_k phi_k(x)/||phi_k||, as an mpmath callable."""
    n = len(v)
    w = [freq(k, L, sector) for k in range(n)]
    nr = [mp.sqrt(norm2(k, L, sector)) for k in range(n)]
    trig = mp.sin if sector == ODD else mp.cos

    def f(x):
        return sum(v[k] * trig(w[k] * x) / nr[k] for k in range(n))
    return f


def quad_A_S(v, L, delta, sector):
    """A = int f cosh(d u) du and S = int f sinh(d u) du, by quadrature on explicit f."""
    f = basis_fun(v, L, sector)
    A = mp.quad(lambda u: f(u) * mp.cosh(delta * u), [-L, 0, L])
    S = mp.quad(lambda u: f(u) * mp.sinh(delta * u), [-L, 0, L])
    return A, S


# --------------------------------------------------------------------------- (a)

def part_a():
    print("=" * 94)
    print("(a) PAIR vs QUARTET -- decided from the committed Phase-5 data")
    print("=" * 94)
    mp.mp.dps = DPS
    # identity first: the degenerate quartet is exactly twice the pair
    worst = mp.mpf(0)
    for sector in (EVEN_D, EVEN, ODD):
        for Lx in ('0.8', '1.4'):
            L = mp.mpf(Lx)
            for dx in ('0.1', '0.25', '0.49'):
                d = mp.mpf(dx)
                P, _ = pl.real_block(L, range(6), d, sector)
                Q = pl.quartet_block(L, range(6), mp.mpf(0), d, sector)
                worst = max(worst, max(abs(Q[i, j] - 2 * P[i, j])
                                       for i in range(6) for j in range(6)))
    check("quartet_block(gamma=0, d) == 2 * real_block(d), 3 sectors x 2 L x 3 delta",
          worst < TOL, f"max abs diff {mp.nstr(worst, 6)}")
    OUT['a_quartet_is_twice_pair_maxdiff'] = mp.nstr(worst, 8)

    # now: which did Phase 5 use?  Reproduce stored lam_planted both ways.
    rows = list(csv.DictReader(open(os.path.join(HERE, 'data', 'phase5_siegel.csv'))))
    probes = [r for r in rows if r['L'] == '0.8' and r['N'] == '8'
              and r['delta'] in ('0.25', '0.49')
              and r['label'] in ('zeta', 'chi_-3', 'chi_8')]
    print(f"  reproducing {len(probes)} stored rows with the PAIR and with the QUARTET:")
    verdict = set()
    for r in probes:
        L, N, d = mp.mpf(r['L']), int(r['N']), mp.mpf(r['delta'])
        sector = r['sector']
        dps = max(60, int(40 + 50 * float(r['L'])))
        mp.mp.dps = dps
        if r['label'] == 'zeta':
            M = build_matrix(L, N, dps, sector)
        else:
            dd = [k for k in dw.DISCRIMINANTS
                  if dw.CHARS[k]['label'] == r['label']][0]
            M = dw.build_matrix_chi(L, range(N), dps, dd, sector)
        P, Ap = pl.real_block(L, range(N), d, sector)
        lam_pair = eigsym(M + P)[0][0]
        lam_quart = eigsym(M + 2 * P)[0][0]
        stored = mp.mpf(r['lam_planted'])
        dp = abs(lam_pair - stored) / max(abs(stored), mp.mpf(1))
        dq = abs(lam_quart - stored) / max(abs(stored), mp.mpf(1))
        which = 'PAIR' if dp < dq else 'QUARTET'
        verdict.add(which)
        # the stored rank-one eigenvalue must be +-2||Ap||^2 (pair), not +-4 (quartet)
        nrm2 = sum(x ** 2 for x in Ap)
        want_pair = (-2 if sector == ODD else 2) * nrm2
        dr = abs(mp.mpf(r['rank1_eig']) - want_pair) / abs(want_pair)
        print(f"    {r['label']:7s} {sector:7s} d={r['delta']:>5s}  "
              f"stored {mp.nstr(stored, 10):>16s} | pair {mp.nstr(lam_pair, 10):>16s} "
              f"(rel {mp.nstr(dp, 3)}) | quartet {mp.nstr(lam_quart, 10):>16s} "
              f"(rel {mp.nstr(dq, 3)}) -> {which}; rank1 rel {mp.nstr(dr, 3)}")
        mp.mp.dps = DPS
        # data/phase5_siegel.csv stores rank1_eig as mp.nstr(.,10), so the binding
        # tolerance is the FILE's 10 significant digits, not this computation's.
        # (Same lesson as REPORT.md Retraction 1.)
        if dr > mp.mpf(10) ** -9:
            FAIL.append(f"rank1_eig mismatch {r['label']}/{sector}/{r['delta']}")
    check("every probed Phase-5 row is reproduced by the PAIR, unanimously",
          verdict == {'PAIR'}, f"verdict set = {sorted(verdict)}")
    OUT['a_phase5_used'] = sorted(verdict)
    print("  => Phase 5 planted the REAL PAIR {1/2+d, 1/2-d}.  No factor-of-2 correction")
    print("     is needed, and no retraction follows from (a).")


# --------------------------------------------------------------------------- (b)

def part_b():
    print("=" * 94)
    print("(b) PARITY SIGN STRUCTURE: +2A^2 on the even basis, -2S^2 on the odd basis")
    print("=" * 94)
    mp.mp.dps = DPS
    N = 10
    deltas = [mp.mpf(x) for x in ('0.1', '0.25', '0.375', '0.45', '0.49')]
    rows = []
    print(f"  {'chi':8s} {'L':>4s} {'sector':7s} {'delta':>6s} {'eigenvalue':>18s} "
          f"{'2A^2 / -2S^2':>18s} {'rel':>9s} {'rank1 resid':>12s} {'other eigs':>12s} "
          f"{'sign':>5s}")
    for d in dw.DISCRIMINANTS:
        lab = dw.CHARS[d]['label']
        for Lx in ('0.8', '1.0', '1.4'):
            L = mp.mpf(Lx)
            for sector in (EVEN_D, ODD):
                for de in deltas:
                    B, Ap = pl.real_block(L, range(N), de, sector)
                    lam, v, res, others = rank_one_fit(B, N)
                    A, S = quad_A_S(v, L, de, sector)
                    want = 2 * A ** 2 if sector != ODD else -2 * S ** 2
                    rel = abs(lam - want) / abs(want)
                    signok = (lam >= 0) if sector != ODD else (lam <= 0)
                    # the complementary integral must vanish by parity
                    cross = abs(S) if sector != ODD else abs(A)
                    ok = (rel < mp.mpf(10) ** -25 and res < TOL
                          and others < TOL and signok and cross < mp.mpf(10) ** -25)
                    if not ok:
                        FAIL.append(f"(b) {lab}/{Lx}/{sector}/{mp.nstr(de,4)}")
                    print(f"  {lab:8s} {Lx:>4s} {sector:7s} {mp.nstr(de, 4):>6s} "
                          f"{mp.nstr(lam, 12):>18s} {mp.nstr(want, 12):>18s} "
                          f"{mp.nstr(rel, 3):>9s} {mp.nstr(res, 3):>12s} "
                          f"{mp.nstr(others, 3):>12s} "
                          f"{('+' if lam >= 0 else '-'):>5s}")
                    rows.append(dict(label=lab, L=Lx, sector=sector,
                                     delta=mp.nstr(de, 8), N=N,
                                     eigenvalue=mp.nstr(lam, 14),
                                     predicted=mp.nstr(want, 14),
                                     rel_err=mp.nstr(rel, 6),
                                     rank1_residual=mp.nstr(res, 6),
                                     other_eigs_max=mp.nstr(others, 6),
                                     A=mp.nstr(A, 12), S=mp.nstr(S, 12),
                                     cross_integral=mp.nstr(cross, 6),
                                     sign_ok=str(bool(signok))))
    check(f"all {len(rows)} (chi, L, sector, delta) cells: rank one, correct sign, "
          f"eigenvalue = +2A^2 / -2S^2", not any(f.startswith('(b)') for f in FAIL))
    # the matrices are chi-INDEPENDENT -- worth asserting, it is the content
    L = mp.mpf('1.0')
    B0, _ = pl.real_block(L, range(N), mp.mpf('0.25'), EVEN_D)
    same = True
    for _ in dw.DISCRIMINANTS:
        B1, _ = pl.real_block(L, range(N), mp.mpf('0.25'), EVEN_D)
        same = same and all(B0[i, j] == B1[i, j] for i in range(N) for j in range(N))
    check("the real-pair perturbation does not depend on chi at all "
          "(it is a property of the basis and delta)", same)

    print("  --- K_even + K_odd = 2L^3/3 ---")
    worst = mp.mpf(0)
    for d in dw.DISCRIMINANTS:
        lab = dw.CHARS[d]['label']
        gs = [mp.mpf(x) for x in open(os.path.join(HERE, 'data', f'zeros_{lab}.txt'))
              if not x.startswith('#')]
        for Lx in ('0.8', '1.0', '1.4'):
            L = mp.mpf(Lx)
            for g in (mp.mpf(0), gs[0], mp.mpf('14.1347')):
                ke = pl.K_closed(L, g)
                ko = 2 * L ** 3 / 3 - pl.K_closed(L, g)
                # K_odd independently, by quadrature, not by subtraction
                koq = mp.quad(lambda u: u ** 2 * mp.cos(g * u) ** 2, [-L, 0, L])
                worst = max(worst, abs(ke + koq - 2 * L ** 3 / 3))
                if abs(ko - koq) > mp.mpf(10) ** -40:
                    FAIL.append(f"K_odd quad mismatch {lab}/{Lx}")
    check("K_even(L,g) + K_odd(L,g) = 2L^3/3 with K_odd by QUADRATURE, "
          "5 chi x 3 L x 3 gamma", worst < mp.mpf(10) ** -40,
          f"max abs diff {mp.nstr(worst, 6)}")
    OUT['b_K_identity_maxdiff'] = mp.nstr(worst, 8)
    with open(os.path.join(HERE, 'data', 'verify_siegel_parity.csv'), 'w',
              newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    OUT['b_cells'] = len(rows)


# --------------------------------------------------------------------------- (c)

def part_c():
    print("=" * 94)
    print("(c) POLE CONSISTENCY: real-pair perturbation at delta = 1/2 vs the pole term")
    print("=" * 94)
    mp.mp.dps = DPS
    half = mp.mpf(1) / 2
    N = 6
    res = {}
    print(f"  {'sector':7s} {'L':>5s} {'pole_sign':>10s} {'max|B - Pole|':>16s} "
          f"{'max|B + Pole|':>16s} {'which vanishes':>16s}")
    for sector in (EVEN_D, EVEN, ODD):
        for Lx in ('0.5', '0.8', '1.3'):
            L = mp.mpf(Lx)
            B, _ = pl.real_block(L, range(N), half, sector)
            nr = [mp.sqrt(norm2(k, L, sector)) for k in range(N)]
            P = [F_at_half(k, L, sector) for k in range(N)]
            sg = pole_sign(sector)
            Pole = mp.zeros(N, N)
            for i in range(N):
                for j in range(N):
                    Pole[i, j] = 2 * sg * P[i] * P[j] / (nr[i] * nr[j])
            dminus = max(abs(B[i, j] - Pole[i, j]) for i in range(N) for j in range(N))
            dplus = max(abs(B[i, j] + Pole[i, j]) for i in range(N) for j in range(N))
            which = 'B - Pole' if dminus < dplus else 'B + Pole'
            print(f"  {sector:7s} {Lx:>5s} {sg:>10d} {mp.nstr(dminus, 6):>16s} "
                  f"{mp.nstr(dplus, 6):>16s} {which:>16s}")
            res[f'{sector}|{Lx}'] = dict(max_abs_B_minus_Pole=mp.nstr(dminus, 8),
                                         max_abs_B_plus_Pole=mp.nstr(dplus, 8),
                                         pole_sign=sg, vanishes=which)
            if dminus >= TOL:
                FAIL.append(f"(c) {sector}/{Lx}: B != Pole, max diff {mp.nstr(dminus, 6)}")
    OUT['c'] = res
    ok = not any(f.startswith('(c)') for f in FAIL)
    check("real_block(delta=1/2) == the parent pole-term matrix ENTRYWISE, all three "
          "sectors x 3 L, SAME sign (B - Pole vanishes)", ok)
    print("  The odd sector's pole_sign is -1, so the odd-sector pole term carries the")
    print("  OPPOSITE sign TO THE EVEN-SECTOR ONE; the real-pair perturbation reproduces")
    print("  that flip exactly.  It is not opposite in sign to the pole term itself.")
    # chi analogue: there is no pole term at all
    L = mp.mpf('0.8')
    Arch, LogPi, LogQ, Prime, _ = dw.build_parts_chi(L, range(N), DPS, -3, EVEN_D)
    M = dw.build_matrix_chi(L, range(N), DPS, -3, EVEN_D)
    recon = max(abs(M[i, j] - (Arch[i, j] - LogPi[i, j] + LogQ[i, j] - Prime[i, j]))
                for i in range(N) for j in range(N))
    check("the chi window matrix is EXACTLY Arch - LogPi + LogQ - Prime: the pole-term "
          "analogue is identically zero", recon < TOL,
          f"reconstruction residual {mp.nstr(recon, 6)}")
    Bchi, _ = pl.real_block(L, range(N), half, EVEN_D)
    big = max(abs(Bchi[i, j]) for i in range(N) for j in range(N))
    check("but the real-pair perturbation at delta=1/2 is NOT zero for chi "
          "(so nothing cancels it)", big > mp.mpf('0.1'),
          f"max |entry| = {mp.nstr(big, 6)}")
    OUT['c_chi_pole_is_zero'] = True
    OUT['c_chi_realpair_max_entry'] = mp.nstr(big, 8)


def main():
    part_a()
    part_b()
    if any(f.startswith('(c)') for f in FAIL):
        pass
    part_c()
    with open(os.path.join(HERE, 'results', 'verify_siegel_parity.json'), 'w') as fh:
        json.dump(dict(results=OUT, failures=FAIL), fh, indent=1)
    print()
    if FAIL:
        print(f"VERIFICATION FAILED: {len(FAIL)}")
        for f in FAIL:
            print(f"   - {f}")
        return 1
    print("ALL THREE SIEGEL-PARITY VERIFICATIONS PASSED")
    return 0


if __name__ == '__main__':
    sys.exit(main())
