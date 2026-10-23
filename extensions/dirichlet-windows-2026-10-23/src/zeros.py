#!/usr/bin/env python3
r"""Phase 2a: critical-line zeros of L(s, chi) for the five real primitive characters.

Method.  For real primitive chi the root number is eps(chi) = 1 (selftest check 2), so
Lambda(1/2 + it, chi) is real and Z(t) = e^{i theta(t)} L(1/2 + it, chi) is a Hardy-type Z
function whose real zeros are the critical-line ordinates.  We

  1. walk an ADAPTIVE grid with step = (local mean spacing)/steps_per_gap, where the mean
     spacing is 2 pi / log(q t / 2 pi);
  2. bracket every sign change and refine by bisection then mp.findroot (secant) at the
     requested precision;
  3. COMPLETENESS CHECK: compare the number found below T against theta(T)/pi, whose
     difference is the usual O(log qT) term S(T).  A jump in that difference flags a
     missed (or spurious) zero and the run reports it rather than hiding it.

Everything this produces is T2 (REPORT.md non-claim 6): the assembly is this extension's
own code and carries no external certificate.  Cross-checking the leading ordinates
against LMFDB is a separate script (src/lmfdb_check.py).
"""
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp                                                      # noqa: E402
import dirichlet_window as dw                                            # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, 'data')


def mean_spacing(t, q):
    t = max(mp.mpf(t), mp.mpf('0.5'))
    lg = mp.log(q * t / (2 * mp.pi))
    return 2 * mp.pi / lg if lg > mp.mpf('0.2') else mp.mpf(10)


def find_zeros(d, Tmax, dps_scan=20, dps_root=60, steps_per_gap=8, tmin='0.05'):
    q = dw.CHARS[d]['q']
    out = []
    flags = []
    mp.mp.dps = dps_scan
    t = mp.mpf(tmin)
    zt = dw.Z(t, d)
    nev = 1
    while t < Tmax:
        h = mean_spacing(t, q) / steps_per_gap
        h = min(h, mp.mpf('0.4'))
        t2 = t + h
        z2 = dw.Z(t2, d)
        nev += 1
        if zt == 0:
            flags.append(('exact-zero-on-grid', float(t)))
        if zt * z2 < 0:
            lo, hi, zlo = t, t2, zt
            for _ in range(22):                       # bisect in low precision
                mid = (lo + hi) / 2
                zm = dw.Z(mid, d)
                nev += 1
                if zm == 0:
                    lo = hi = mid
                    break
                if zlo * zm < 0:
                    hi = mid
                else:
                    lo, zlo = mid, zm
            mp.mp.dps = dps_root + 10
            try:
                r = mp.findroot(lambda x: dw.Z(x, d), (lo + hi) / 2, tol=mp.mpf(10) **
                                (-2 * dps_root))
                r = mp.mpf(mp.nstr(mp.re(r), dps_root + 5))
            except Exception as exc:                  # noqa: BLE001
                flags.append((f'findroot failed: {exc}', float((lo + hi) / 2)))
                r = (lo + hi) / 2
            out.append(r)
            mp.mp.dps = dps_scan
        t, zt = t2, z2
    return out, flags, nev


def completeness(d, gammas, Tgrid):
    """found(T) - theta(T)/pi for a ladder of T; should move by O(log qT), not by 1 per step."""
    mp.mp.dps = 30
    rows = []
    for T in Tgrid:
        n = sum(1 for g in gammas if g < T)
        th = dw.theta_chi(T, d) / mp.pi
        rows.append((float(T), n, float(th), n - float(th)))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--Tmax', type=float, default=275.0)
    ap.add_argument('--dps-root', type=int, default=60)
    ap.add_argument('--dps-scan', type=int, default=20)
    ap.add_argument('--steps-per-gap', type=int, default=8)
    ap.add_argument('--discs', default=','.join(str(x) for x in dw.DISCRIMINANTS))
    a = ap.parse_args()
    os.makedirs(DATA, exist_ok=True)
    summary = {}
    for d in [int(x) for x in a.discs.split(',')]:
        lab = dw.CHARS[d]['label']
        t0 = time.time()
        gs, flags, nev = find_zeros(d, mp.mpf(a.Tmax), a.dps_scan, a.dps_root,
                                    a.steps_per_gap)
        dt = time.time() - t0
        mp.mp.dps = a.dps_root
        path = os.path.join(DATA, f'zeros_{lab}.txt')
        with open(path, 'w') as fh:
            fh.write(f"# critical-line ordinates of L(s,{lab}), q={dw.CHARS[d]['q']}, "
                     f"a={dw.CHARS[d]['a']}\n")
            fh.write(f"# T2: this extension's own root-finding, dps_root={a.dps_root}, "
                     f"Tmax={a.Tmax}\n")
            for g in gs:
                fh.write(mp.nstr(g, a.dps_root) + "\n")
        comp = completeness(d, gs, [25, 50, 100, 150, 200, 250, a.Tmax])
        print(f"=== {lab} (d={d}, q={dw.CHARS[d]['q']}, a={dw.CHARS[d]['a']}) ===")
        print(f"    {len(gs)} ordinates in (0, {a.Tmax}); {nev} Z-evaluations; {dt:.0f}s")
        print(f"    gamma_1 = {mp.nstr(gs[0], 25)}")
        print(f"    gamma_{len(gs)} = {mp.nstr(gs[-1], 15)}")
        print(f"    {'T':>8s} {'found':>6s} {'theta(T)/pi':>13s} {'S-ish':>9s}")
        for T, n, th, diff in comp:
            print(f"    {T:>8.1f} {n:>6d} {th:>13.4f} {diff:>9.4f}")
        if flags:
            print(f"    FLAGS: {flags}")
        summary[lab] = dict(d=d, q=dw.CHARS[d]['q'], a=dw.CHARS[d]['a'],
                            n_found=len(gs), Tmax=a.Tmax,
                            gamma_1=mp.nstr(gs[0], 30), gamma_last=mp.nstr(gs[-1], 20),
                            completeness=comp, flags=[str(f) for f in flags],
                            dps_root=a.dps_root, seconds=round(dt, 1))
        print()
    with open(os.path.join(HERE, 'results', 'phase2_zeros_summary.json'), 'w') as fh:
        json.dump(summary, fh, indent=1)
    print('wrote results/phase2_zeros_summary.json')


if __name__ == '__main__':
    main()
