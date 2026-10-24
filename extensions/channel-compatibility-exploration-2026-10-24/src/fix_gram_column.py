#!/usr/bin/env python3
r"""Recompute the gram_dev column of data/channel_results.csv at each row's precision.

Why this exists.  The first version of `run_config` computed
    gram_dev = max |norm2(k)/norm2(k) - 1|
which is identically zero and could not fail.  It was replaced by a real computation of
the Gram matrix through the parent extension's phaseD.gram, but the grid and ladder rows
had already been written with the tautological value 0.0.

gram_dev is an INDEPENDENT CHECK: it is computed, stored, and never used by any other
quantity, so recomputing just this column is sound and does not require re-running the
expensive matrix assemblies.  That is the only reason this script is legitimate; a
column that fed eta would have to be regenerated from scratch.
"""
import csv
import os
import sys

DH = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), 'davenport-heilbronn-2026-10-23', 'src')
sys.path.insert(0, DH)
import mpmath as mp                                                      # noqa: E402
import dh                                                                # noqa: E402
from phaseD import index_set, gram                                       # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV = os.path.join(HERE, 'data', 'channel_results.csv')


def main():
    rows = list(csv.DictReader(open(CSV)))
    keys = list(rows[0].keys())
    for r in rows:
        dps = int(r['dps'])
        mp.mp.dps = dps
        L = dh.mpf_hi(r['L'])
        idx = index_set(L, 'odd', int(r['level']) - 1)
        G = gram(idx, L, 'odd')
        dev = max(abs(G[i, j] - (1 if i == j else 0))
                  for i in range(len(idx)) for j in range(len(idx)))
        r['gram_dev'] = mp.nstr(dev, 4)
        print('  L=%-4s lvl%s N=%-3d dps=%-4d  max|G - I| = %s'
              % (r['L'], r['level'], len(idx), dps, r['gram_dev']))
    with open(CSV, 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print('  -> %s (gram_dev column recomputed)' % os.path.relpath(CSV, HERE))
    return 0


if __name__ == '__main__':
    sys.exit(main())
