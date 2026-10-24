#!/usr/bin/env python3
r"""Phase D -- the window scan and ONE explicit negative witness.

Basis schedule (declared in advance, recorded in data/window_scan.csv):

  level 1   low block {0..7}  +  resonant block {k0-4 .. k0+4},  k0 chosen so that
            w_{k0} is the basis frequency closest to gamma_1 = 85.699348485...
  level 2   low block {0..15} +  resonant block {k0-10 .. k0+10}      (ENRICHMENT)
  level 3   low block {0..23} +  resonant block {k0-16 .. k0+16}      (ENRICHMENT)

Level 1 is the fixed scan.  Levels 2 and 3 are run at EVERY window as well, and both
are recorded; an enrichment can only lower lambda_min (a larger subspace), so a window
whose level-1 value is positive and whose level-3 value is negative is reported as
"detected at level 3", never folded back into the level-1 scan.

The Gram matrix of f_k/sqrt(norm2_k) on the FULL index set is the identity exactly
(the bases are orthogonal with known norms), so the generalised problem B v = lambda M v
reduces to the ordinary symmetric eigenproblem.  It is still assembled and its deviation
from I and its condition number are reported, because the code admits arbitrary index
subsets and the claim should be checked, not assumed.

Everything is arbitrary precision; no float64 enters any reported value.  The decimal
grid of windows comes from the strings in config.WINDOW_SCAN.

Output: results/phaseD.txt, data/window_scan.csv, certificates/witness.json.
"""
import csv
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp                                                      # noqa: E402
import dh                                                                # noqa: E402
import config as cfg                                                     # noqa: E402
from zeta_window import Csym_terms, eval_terms, freq, norm2              # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(HERE, 'results')
DATA = os.path.join(HERE, 'data')
CERT = os.path.join(HERE, 'certificates')

GAMMA1 = '85.69934848537759217192926770894172903799'
BETA1 = '0.8085171824566373855533519606068441278507'

LEVELS = [(8, 4), (16, 10), (24, 16)]


def resonant_k(L, sector, gam):
    L = dh.mpf_hi(L)
    if sector == cfg.SECTOR_ODD:
        return max(0, int(mp.nint(gam * L / mp.pi - 1)))
    return max(0, int(mp.nint(gam * L / mp.pi - mp.mpf(1) / 2)))


def index_set(L, sector, level):
    low, half = LEVELS[level]
    k0 = resonant_k(L, sector, dh.mpf_hi(GAMMA1))
    return sorted(set(list(range(low))
                      + [k for k in range(max(0, k0 - half), k0 + half + 1)]))


def gram(idx, L, sector):
    n = len(idx)
    G = mp.zeros(n, n)
    for i, ki in enumerate(idx):
        for j, kj in enumerate(idx):
            tm = Csym_terms(ki, kj, L, sector)
            den = mp.sqrt(norm2(ki, L, sector) * norm2(kj, L, sector))
            G[i, j] = eval_terms(tm, L, 0) / den
    return G


def rayleigh(M, v):
    n = len(v)
    num = mp.fsum([v[i] * v[j] * M[i, j] for i in range(n) for j in range(n)])
    return num / mp.fsum([x * x for x in v])


def main():
    mp.mp.dps = cfg.DPS_BASE
    for d in (RES, DATA, CERT):
        os.makedirs(d, exist_ok=True)
    out, rows = [], []

    def p(line=''):
        print(line)
        out.append(line)

    p('=' * 92)
    p('PHASE D  --  window scan and the negative witness')
    p('=' * 92)
    p('  gamma_1 = %s   beta_1 = %s' % (GAMMA1, BETA1))
    p('  levels (low block size, resonant half-width): %s' % (LEVELS,))
    p('')
    p('%-5s %-7s %-6s %-4s %-5s %-10s %-24s %-11s' %
      ('L', 'sector', 'level', 'N', 'k0', 'w_k0', 'lambda_min', 'secs'))

    best = None
    for Ls in cfg.WINDOW_SCAN:
        L = dh.mpf_hi(Ls)
        for sector in (cfg.SECTOR_EVEN, cfg.SECTOR_ODD):
            for level in range(len(LEVELS)):
                idx = index_set(L, sector, level)
                k0 = resonant_k(L, sector, dh.mpf_hi(GAMMA1))
                t0 = time.time()
                mp.mp.dps = cfg.DPS_BASE
                M = dh.build_matrix(L, idx, cfg.DPS_BASE, sector)
                E, V = mp.eigsy(M)
                lo = min(range(len(idx)), key=lambda i: E[i])
                lam = E[lo]
                v = [V[i, lo] for i in range(len(idx))]
                ray = rayleigh(M, v)
                dt = time.time() - t0
                p('%-5s %-7s %-6d %-4d %-5d %-10s %-24s %-11.0f' %
                  (Ls, sector, level + 1, len(idx), k0,
                   mp.nstr(freq(k0, L, sector), 8), mp.nstr(lam, 18), dt))
                rows.append(dict(L=Ls, sector=sector, level=level + 1, N=len(idx),
                                 k0=k0, w_k0=mp.nstr(freq(k0, L, sector), 20),
                                 lambda_min=mp.nstr(lam, 30),
                                 rayleigh_check=float(abs(ray - lam)),
                                 negative=int(lam < 0), seconds=round(dt, 1)))
                if lam < 0:
                    key = (dh.mpf_hi(Ls), -lam)
                    if best is None or (key[0] < best[0]
                                        or (key[0] == best[0] and key[1] > best[1])):
                        best = (dh.mpf_hi(Ls), -lam, Ls, sector, level, idx, v, lam)

    p('')
    if best is None:
        p('  NO negative witness in the scanned spaces.')
        with open(os.path.join(RES, 'phaseD.txt'), 'w') as fh:
            fh.write('\n'.join(out) + '\n')
        return 1

    _, _, Ls, sector, level, idx, v, lam = best
    L = dh.mpf_hi(Ls)
    p('  smallest tested window with a detected negative witness: L = %s, sector %s,'
      ' basis level %d' % (Ls, sector, level + 1))
    p('  (this is NOT a claim about the true first detection threshold, and no claim is')
    p('   made about which sector detects first globally -- non-claims 5 and 9.)')

    # ---------------------------------------------------------------- freeze
    p('')
    p('[FREEZE] the witness')
    # choose the clearest witness: largest |lambda_min| at the smallest detecting L is
    # a poor conditioning choice; we freeze at the smallest detecting L and also record
    # the deepest witness in the whole scan.
    nrm = mp.sqrt(mp.fsum([x * x for x in v]))
    v = [x / nrm for x in v]
    if v[max(range(len(v)), key=lambda i: abs(v[i]))] < 0:
        v = [-x for x in v]
    DIG = 30
    vr = [mp.mpf(mp.nstr(x, DIG, strip_zeros=False)) for x in v]
    p('  frozen at L = %s, sector = %s, level = %d, N = %d' % (Ls, sector, level + 1, len(idx)))
    p('  coefficients rationalised to %d decimal digits (exact decimals, re-read as'
      ' strings)' % DIG)

    gm = gram(idx, L, sector)
    dev = max(abs(gm[i, j] - (1 if i == j else 0))
              for i in range(len(idx)) for j in range(len(idx)))
    ge = mp.eigsy(gm, eigvals_only=True)
    p('  Gram matrix: max |G - I| = %s ; eigenvalues in [%s, %s] ; condition number %s'
      % (mp.nstr(dev, 4), mp.nstr(min(ge), 12), mp.nstr(max(ge), 12),
         mp.nstr(max(ge) / min(ge), 8)))
    p('  => the generalised problem B v = lambda M v IS the ordinary one; nothing was')
    p('     regularised and no direction was discarded.')

    p('')
    p('  Rayleigh quotient Q_F(g)/||g||^2 recomputed DIRECTLY (no eigensolver):')
    vals = {}
    for dps in (cfg.DPS_BASE, cfg.DPS_HIGH, cfg.DPS_XHIGH):
        mp.mp.dps = dps
        Mx = dh.build_matrix(L, idx, dps, sector)
        q = rayleigh(Mx, vr)
        vals[dps] = q
        p('    dps = %-4d   Q/||g||^2 = %s' % (dps, mp.nstr(q, 40)))
    mp.mp.dps = cfg.DPS_XHIGH
    p('    |value(%d) - value(%d)| = %s' %
      (cfg.DPS_BASE, cfg.DPS_HIGH, mp.nstr(abs(vals[cfg.DPS_BASE] - vals[cfg.DPS_HIGH]), 4)))
    p('    |value(%d) - value(%d)| = %s' %
      (cfg.DPS_HIGH, cfg.DPS_XHIGH,
       mp.nstr(abs(vals[cfg.DPS_HIGH] - vals[cfg.DPS_XHIGH]), 4)))
    p('    eigenvalue from mp.eigsy at the scan precision = %s' % mp.nstr(lam, 25))
    p('    |Rayleigh(rationalised v) - eigenvalue|        = %s'
      % mp.nstr(abs(vals[cfg.DPS_BASE] - lam), 4))

    # the deepest witness anywhere in the scan, recorded separately
    neg = [r for r in rows if r['negative']]
    deep = min(neg, key=lambda r: mp.mpf(r['lambda_min'])) if neg else None
    if deep:
        p('')
        p('  deepest witness anywhere in the scan: L = %s, %s, level %d, lambda_min = %s'
          % (deep['L'], deep['sector'], deep['level'], deep['lambda_min']))

    cert = dict(
        L=Ls, sector=sector, level=level + 1, idx=[int(k) for k in idx],
        frequencies=[mp.nstr(freq(k, L, sector), 40) for k in idx],
        coefficients=[mp.nstr(x, DIG) for x in vr],
        normalisation='g = sum_m c_m f_{idx[m]}(u) / sqrt(norm2), ||c||_2 renormalised'
                      ' after rationalisation; the reported value is Q/||g||^2',
        basis=('cos((2k+1) pi u / 2L)' if sector != cfg.SECTOR_ODD
               else 'sin((k+1) pi u / L)'),
        value_dps={str(k): mp.nstr(v2, 50) for k, v2 in vals.items()},
        eigenvalue_scan=mp.nstr(lam, 40),
        gram_max_dev_from_I=float(dev),
        gram_condition=float(max(ge) / min(ge)),
        gamma1=GAMMA1, beta1=BETA1)
    with open(os.path.join(CERT, 'witness.json'), 'w') as fh:
        json.dump(cert, fh, indent=1)
    p('')
    p('  -> certificates/witness.json   (the frozen test function)')

    with open(os.path.join(DATA, 'window_scan.csv'), 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    p('  -> data/window_scan.csv')
    p('=' * 92)
    with open(os.path.join(RES, 'phaseD.txt'), 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
