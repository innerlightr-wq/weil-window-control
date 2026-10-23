#!/usr/bin/env python3
r"""Phase 2b: cross-check the computed ordinates against LMFDB.

The L-function is NOT identified by guessing a Conrey label.  For each conductor we pull
every degree-1 primitive L-function from the LMFDB API and match its `euler_factors`
against this extension's own chi(p) at every good prime listed -- LMFDB stores the local
factor's second coefficient as chi(p) (verified on the unambiguous conductor-3 case) -- and
we additionally require the parity recorded in the label (`-r0-`/`-r1-`) to equal our `a`.
Only then are the zeros compared.

Writes results/phase2_lmfdb.json with the raw LMFDB ordinates kept, so the comparison is
re-checkable without a network.
"""
import json
import os
import sys
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp                                                      # noqa: E402
import dirichlet_window as dw                                            # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
API = ('https://www.lmfdb.org/api/lfunc_lfunctions/?degree=1&conductor={q}'
       '&_format=json&_fields=label,central_character,euler_factors,positive_zeros'
       '&_offset={off}')
PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]


def fetch(q):
    rows, off = [], 0
    while True:
        with urllib.request.urlopen(API.format(q=q, off=off), timeout=90) as fh:
            d = json.load(fh)
        rows += d['data']
        if len(d['data']) == 0 or d.get('offset', 0) + len(d['data']) <= len(rows) - 1:
            break
        if len(d['data']) < 100:
            break
        off += len(d['data'])
    return rows


def match(d, rows):
    q, a = dw.CHARS[d]['q'], dw.CHARS[d]['a']
    hits = []
    for r in rows:
        ef = r.get('euler_factors') or []
        ok = True
        n_cmp = 0
        for p, fac in zip(PRIMES, ef):
            if q % p == 0:
                continue
            if not isinstance(fac, list) or len(fac) < 2 or isinstance(fac[1], str):
                ok = False
                break
            if int(fac[1]) != dw.chi(d, p):
                ok = False
                break
            n_cmp += 1
        if ok and n_cmp >= 4 and f'-r{a}-' in r['label']:
            hits.append((r, n_cmp))
    return hits


def main():
    mp.mp.dps = 50
    out = {}
    worst_all = mp.mpf(0)
    for d in dw.DISCRIMINANTS:
        lab = dw.CHARS[d]['label']
        q = dw.CHARS[d]['q']
        path = os.path.join(HERE, 'data', f'zeros_{lab}.txt')
        mine = [mp.mpf(x) for x in open(path) if not x.startswith('#')]
        try:
            rows = fetch(q)
        except Exception as exc:                                        # noqa: BLE001
            print(f"{lab}: NETWORK FAILED ({exc}); ordinates stay T2")
            out[lab] = dict(status='network-failed', error=str(exc))
            continue
        hits = match(d, rows)
        if len(hits) != 1:
            print(f"{lab}: identification AMBIGUOUS ({len(hits)} matches) -- not used")
            out[lab] = dict(status='ambiguous', n_matches=len(hits),
                            labels=[h[0]['label'] for h in hits])
            continue
        r, n_cmp = hits[0]
        theirs = [mp.mpf(z) for z in r['positive_zeros']]
        n = min(len(mine), len(theirs))
        diffs = [abs(mine[i] - theirs[i]) for i in range(n)]
        worst = max(diffs)
        worst_all = max(worst_all, worst)
        rel = worst / theirs[n - 1]
        print(f"{lab}: LMFDB {r['label']} (chi matched at {n_cmp} primes), "
              f"{n} ordinates compared, max |diff| = {mp.nstr(worst, 5)} "
              f"(rel {mp.nstr(rel, 4)})")
        print(f"       gamma_1 mine   = {mp.nstr(mine[0], 30)}")
        print(f"       gamma_1 LMFDB  = {mp.nstr(theirs[0], 30)}")
        out[lab] = dict(status='ok', lmfdb_label=r['label'],
                        central_character=r['central_character'],
                        primes_matched=n_cmp, n_compared=n,
                        max_abs_diff=mp.nstr(worst, 8),
                        max_rel_diff=mp.nstr(rel, 6),
                        n_mine=len(mine), n_lmfdb=len(theirs),
                        lmfdb_positive_zeros=[str(z) for z in r['positive_zeros']])
    out['_worst_abs_diff_all'] = mp.nstr(worst_all, 8)
    with open(os.path.join(HERE, 'results', 'phase2_lmfdb.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
    print(f"\nworst |diff| over all characters and all compared ordinates: "
          f"{mp.nstr(worst_all, 6)}")
    print('wrote results/phase2_lmfdb.json')


if __name__ == '__main__':
    main()
