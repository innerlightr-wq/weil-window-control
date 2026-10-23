#!/usr/bin/env python3
r"""Regenerate the zero-set summary and completeness table from data/zeros_*.txt.

Kept separate from src/zeros.py so the record is reproducible from the stored ordinates
alone, independently of how (or in how many invocations) they were produced.

COMPLETENESS CHECK.  The number of critical-line ordinates below T is
theta_chi(T)/pi + c + S(T) with S(T) = O(log qT).  A MISSED zero shows up as the printed
`found - theta/pi` column stepping down by 1 and staying there; a SPURIOUS one steps up.
A column that stays inside a narrow band across the whole ladder is the pass condition.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp                                                      # noqa: E402
import dirichlet_window as dw                                            # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GRID = [10, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275]


def main():
    mp.mp.dps = 40
    out = {}
    print(f"  {'chi':8s} {'q':>2s} {'a':>2s} {'#zeros':>7s} {'gamma_1':>26s} "
          f"{'gamma_last':>14s} {'max |found-theta/pi|':>21s} {'band':>14s}")
    for d in dw.DISCRIMINANTS:
        lab = dw.CHARS[d]['label']
        path = os.path.join(HERE, 'data', f'zeros_{lab}.txt')
        gs = [mp.mpf(x) for x in open(path) if not x.startswith('#')]
        comp = []
        for T in GRID:
            n = sum(1 for g in gs if g < T)
            th = dw.theta_chi(T, d) / mp.pi
            comp.append((T, n, float(th), float(n - th)))
        diffs = [c[3] for c in comp if c[0] <= float(gs[-1])]
        band = max(diffs) - min(diffs)
        print(f"  {lab:8s} {dw.CHARS[d]['q']:>2d} {dw.CHARS[d]['a']:>2d} "
              f"{len(gs):>7d} {mp.nstr(gs[0], 24):>26s} {mp.nstr(gs[-1], 12):>14s} "
              f"{max(abs(x) for x in diffs):>21.4f} {band:>14.4f}")
        out[lab] = dict(d=d, q=dw.CHARS[d]['q'], a=dw.CHARS[d]['a'], n_zeros=len(gs),
                        gamma_1=mp.nstr(gs[0], 40), gamma_last=mp.nstr(gs[-1], 24),
                        completeness_T_found_theta_diff=comp,
                        max_abs_diff=round(max(abs(x) for x in diffs), 5),
                        band=round(band, 5))
    with open(os.path.join(HERE, 'results', 'phase2_zeros_summary.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
    print("\n  per-character completeness ladder (found - theta(T)/pi):")
    for lab, v in out.items():
        print(f"    {lab:8s} " + " ".join(
            f"T={c[0]}:{c[3]:+.3f}" for c in v['completeness_T_found_theta_diff']
            if c[0] <= float(v['gamma_last'])))
    print("\nwrote results/phase2_zeros_summary.json")


if __name__ == '__main__':
    main()
