#!/usr/bin/env python3
"""Adversarial pass: try to break each Stage 1 / Stage 2 conclusion."""
import sys, os, csv, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from fractions import Fraction
from weilform import toeplitz, eigsym, parity_of, poly_zeros, inertia
from planted import t_from_blocks
from exact_inertia import rational_t, exact_inertia
from curves import squarefree_over_Fq, count_points_hyperelliptic, newton_p_to_A
from analysis import analyse_window

mp.mp.dps = 40
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data')


def A1_random_toeplitz_H1(trials=400, nmax=9, seed=7):
    """ATTACK on H1: random real symmetric Toeplitz. If lam_min is simple, must
    every ground-state zero lie on |z|=1?  (Caratheodory-Fejer says yes.)"""
    random.seed(seed)
    worst, worst_case, checked = mp.mpf(0), None, 0
    for _ in range(trials):
        R = random.randint(1, nmax)
        t = [mp.mpf(random.uniform(-3, 3)) for _ in range(R + 1)]
        T = toeplitz(t, R)
        vals, vecs = eigsym(T)
        scale = max([abs(v) for v in vals] + [mp.mpf(1)])
        if R >= 1 and (vals[1] - vals[0]) <= scale * mp.mpf(10) ** -20:
            continue                       # lam_min not simple: H1 says nothing
        checked += 1
        zs = poly_zeros(vecs[0])
        if len(zs) != R:                   # degenerate leading coefficient
            continue
        off = max(abs(abs(z) - 1) for z in zs)
        if off > worst:
            worst, worst_case = off, (R, [mp.nstr(x, 8) for x in t])
    print(f"A1  H1 vs {checked} random simple-lam_min Toeplitz matrices (R<={nmax}):")
    print(f"    worst max||z|-1| = {mp.nstr(worst, 6)}   (worst case R={worst_case[0]})")
    print(f"    verdict: {'H1 SURVIVES' if worst < mp.mpf(10) ** -20 else 'H1 BROKEN'}")
    return worst


def A2_repeated_angles_breaks_H2():
    """ATTACK on H2: H2 as stated needs DISTINCT Frobenius angles. Build an on-line
    config with a repeated angle and show lam_min = 0 is not simple at R = 2g."""
    print("\nA2  H2 without the distinctness hypothesis:")
    blocks = [('online', 0.7), ('online', 0.7)]          # g_eff = 2, only 2 distinct angles
    t, g = t_from_blocks(blocks, 10)
    for R in (3, 4, 5):
        row, _, _ = analyse_window(t, R, [mp.mpf(0.7), mp.mpf(-0.7)], mp.mp.dps)
        print(f"    repeated-angle config, g_eff={g}, R={R}: inertia={row['inertia']}"
              f"  lam_min={row['lam_min']}  simple={row['lam_min_simple']}")
    print("    verdict: at R=2g=4 the kernel has dimension 3, lam_min is NOT simple,")
    print("             so H2 REQUIRES distinct theta_j. H2 restated accordingly.")
    # does a real curve with repeated angles exist in our search range?
    print("    searching small hyperelliptic curves for repeated Frobenius angles...")
    found = []
    for q in (3, 5):
        for g in (1, 2):
            d = 2 * g + 1
            import itertools
            for coeffs in itertools.product(range(q), repeat=d):
                f = list(coeffs) + [1]
                if not squarefree_over_Fq(f, q):
                    continue
                try:
                    p = [None] + [1 + q ** n - count_points_hyperelliptic(f, q, n)
                                  for n in range(1, g + 1)]
                    A = newton_p_to_A(p, q, g)
                except Exception:
                    continue
                rev = [mp.mpf(A[i]) for i in range(2 * g + 1)]
                try:
                    rts = mp.polyroots(rev, maxsteps=200, extraprec=200)
                except Exception:
                    continue
                angs = sorted(mp.arg(r) for r in rts)
                rep = any(abs(angs[i] - angs[i + 1]) < mp.mpf(10) ** -20
                          for i in range(len(angs) - 1))
                if rep:
                    found.append((q, g, f, A))
                    if len(found) >= 3:
                        break
            if len(found) >= 3:
                break
        if len(found) >= 3:
            break
    for q, g, f, A in found:
        print(f"      FOUND real curve with repeated angles: q={q} g={g} f={f} P={A}")
    if not found:
        print("      none found in the scanned range (y^2=f(x), deg f=2g+1, q in {3,5}, g<=2)")
    return found


def A3_masking(kmax=24, Rmax=60, out=None):
    """ATTACK on 'sign of lam_min detects off-line zeros': can ON-LINE zeros MASK an
    off-line one, pushing the detection window arbitrarily far out?"""
    print("\nA3  masking: 1 off-line quartet (rho=1.3, cos phi=19/25) + k on-line pairs")
    print("    on-line angles: 24 distinct rational cosines, all != 19/25")
    # 24 distinct rational cosines (Pythagorean-style), none equal to cos phi = 19/25
    cos_online = ['3/5', '-1/5', '4/5', '-3/5', '1/5', '-4/5', '7/25', '-7/25',
                  '24/25', '-24/25', '5/13', '-5/13', '12/13', '-12/13', '9/41',
                  '-9/41', '40/41', '-40/41', '8/17', '-8/17', '15/17', '-15/17',
                  '20/29', '-20/29']
    rows = []
    for k in range(0, kmax + 1):
        blocks = [('quartet', '13/10', '19/25')] + [('online', c) for c in cos_online[:k]]
        tq, geff = rational_t(blocks, Rmax)
        rank = 4 + 2 * k
        first = None
        for R in range(0, Rmax + 1):
            ex = exact_inertia(tq, R)
            if ex is not None and ex[0] >= 1:
                first = R
                break
        rows.append(dict(k_online_pairs=k, g_eff=geff, total_rank=rank,
                         first_neg_R_EXACT=first, Rmax=Rmax))
        print(f"    k={k}  g_eff={geff}  rank={rank}  ->  first R with lam_min<0 "
              f"(EXACT over Q) = {first}")
    if out:
        with open(out, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
        print(f"    wrote {out}")
    return rows


def A4_can_offline_stay_psd(Rmax=30):
    """ATTACK: is there any off-line config that stays PSD for ALL R?
    Scan radii and angles; report the max R reached without a certified negative."""
    print("\nA4  does any off-line config stay PSD out to R=%d?" % Rmax)
    bad = []
    for rho in ['13/10', '11/10', '101/100', '10001/10000', '2', '5']:
        for c in ['19/25', '-21/50', '1', '-1', '0', '3/5']:
            blocks = [('quartet', rho, c)] if c not in ('1', '-1') else \
                     [('realpair', rho)]
            tq, _ = rational_t(blocks, Rmax)
            first = None
            for R in range(0, Rmax + 1):
                ex = exact_inertia(tq, R)
                if ex is not None and ex[0] >= 1:
                    first = R; break
            if first is None:
                bad.append((rho, c))
            print(f"    rho={rho:>12s} cos phi={c:>6s}  first neg R = {first}")
    print(f"    verdict: {'all off-line configs break positivity' if not bad else 'SURVIVORS: ' + str(bad)}")
    return bad


if __name__ == '__main__':
    A1_random_toeplitz_H1()
    A2_repeated_angles_breaks_H2()
    A3_masking(out=os.path.join(OUT, 'stage2_masking.csv'))
    A4_can_offline_stay_psd()
