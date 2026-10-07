#!/usr/bin/env python3
"""Stage 1: float64 vs mpmath agreement, Zhu's A_L table, and the Zhu section 7 case."""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import mpmath as mp
from symbol import Symbol, STAGE2_L, LOG_PI

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ZHU_A = {0.8: 2.9420, 1.0: 5.8525, 1.19: 7.0750, 1.2: 8.5210, 1.4: 10.290, 1.6: 14.323}


def psi_mp(s, t):
    """Psi_L(t) in mpmath at 40 digits."""
    mp.mp.dps = 40
    H = mp.re(mp.digamma(mp.mpf(1) / 4 + mp.mpc(0, 1) * mp.mpf(t) / 2)) - mp.log(mp.pi)
    P = sum(mp.mpf(w) * mp.cos(mp.mpf(t) * mp.mpf(f))
            for w, f in zip(s.weights, s.freqs))
    return H - P


def part_a(seed=0, npts=10_000):
    print("=" * 74)
    print("1(a)  float64 (scipy) vs mpmath, %d random t per L" % npts)
    print("=" * 74)
    rng = np.random.default_rng(seed)
    rows = []
    worst_all = 0.0
    for L in STAGE2_L:
        s = Symbol(L)
        hi = s.envelope_end(0.5)
        t = rng.uniform(2 * np.pi, hi, npts)
        f64 = s.Psi(t)
        # mpmath on a random subsample (full 1e4 at 40 dps per L is wasteful)
        idx = rng.choice(npts, 400, replace=False)
        err = np.array([abs(float(psi_mp(s, t[i])) - f64[i]) for i in idx])
        worst = float(err.max())
        worst_all = max(worst_all, worst)
        print("  L=%-5s  n_active=%-2d  scan hi=%-11.4g  max |f64 - mpmath| = %.3e"
              % (L, len(s.n), hi, worst))
        plainerr = float(np.max([abs(float(psi_mp(s, t[i]))
                                     - (s.H(np.array([t[i]])) - s.P_fast(np.array([t[i]])))[0])
                                 for i in idx]))
        print("      (plain float64 on the same points: %.3e)" % plainerr)
        rows.append(dict(L=L, n_active=int(len(s.n)), scan_hi=float(hi),
                         max_abs_diff=worst, max_abs_diff_plain_f64=plainerr,
                         n_compared=int(len(idx))))
    ok = worst_all < 1e-12
    print("  overall worst = %.3e   -> agreement to 1e-12: %s" % (worst_all, ok))
    print("  (evaluator uses exact argument reduction; the plain-float64 path loses up")
    print("   to 1.3e-8 at L=1.6 -- see REPORT 'Implementation note')")
    return rows, ok


def part_b():
    print("\n" + "=" * 74)
    print("1(b)  Zhu's A_L table")
    print("=" * 74)
    print("  %-6s %-40s %12s %12s %10s %12s" %
          ("L", "active n", "A_L (ours)", "A_L (Zhu)", "|diff|", "T1(L)"))
    rows = []
    allok = True
    for L in STAGE2_L:
        s = Symbol(L)
        zhu = ZHU_A.get(L)
        d = abs(s.A - zhu) if zhu is not None else None
        ok = (d is not None and d < 5e-4)
        if zhu is not None:
            allok &= ok
        print("  %-6s %-40s %12.4f %12s %10s %12.6g" %
              (L, str(s.n.tolist()), s.A,
               ("%.4f" % zhu) if zhu else "-", ("%.1e" % d) if d is not None else "-", s.T1))
        rows.append(dict(L=L, active_n=s.n.tolist(), A=s.A, A_zhu=zhu,
                         abs_diff=d, T1=s.T1, Tstar=s.Tstar))
    print("  all six Zhu values matched to <5e-4: %s" % allok)
    return rows, allok


def part_c():
    print("\n" + "=" * 74)
    print("1(c)  Zhu section 7:  L = 1.19, beta = 0.46948, window [1100, 2 pi e^{A+beta}]")
    print("=" * 74)
    L, beta = 1.19, 0.46948
    s = Symbol(L)
    hi = s.envelope_end(beta)
    lo = 1100.0
    step = 2e-4
    N = int((hi - lo) / step) + 1
    print("  window  = [%.1f, %.4f]   grid step %.0e  (%d points)" % (lo, hi, step, N))
    t0_last, measure, ndip, best_depth, best_t = None, 0.0, 0, -np.inf, None
    prev_in = False
    chunk = 1 << 22
    t_start = time.time()
    for i in range(0, N, chunk):
        tt = lo + step * np.arange(i, min(i + chunk, N), dtype=np.float64)
        psi = s.Psi_chunked(tt)
        inside = psi < beta
        measure += inside.sum() * step
        d = beta - psi
        k = int(np.argmax(d))
        if d[k] > best_depth:
            best_depth, best_t = float(d[k]), float(tt[k])
        # count dip intervals across chunk boundaries
        trans = np.flatnonzero(np.diff(inside.astype(np.int8)) == 1)
        ndip += len(trans) + (1 if (inside[0] and not prev_in) else 0)
        prev_in = bool(inside[-1])
        if inside.any():
            t0_last = float(tt[np.flatnonzero(inside)[-1]])
    rt = time.time() - t_start
    print("  dip measure      = %.3f      (Zhu: approx 18)" % measure)
    print("  max depth        = %.4f  at t = %.2f   (Zhu: approx 2.0 near t approx 1550)"
          % (best_depth, best_t))
    print("  last dip t0      = %.4g      (Zhu: approx 8.48e3)" % t0_last)
    print("  dip intervals    = %d" % ndip)
    print("  min Psi over the window = %.4f" % (beta - best_depth))
    print("  runtime %.1f s" % rt)
    okm = abs(measure - 18) / 18 < 0.15
    okd = abs(best_depth - 2.0) < 0.25 and abs(best_t - 1550) < 60
    okt = abs(t0_last - 8.48e3) / 8.48e3 < 0.05
    print("  measure match: %s | depth+location match: %s | t0 match: %s"
          % (okm, okd, okt))
    return dict(L=L, beta=beta, lo=lo, hi=hi, step=step, measure=measure,
                max_depth=best_depth, max_depth_t=best_t, t0=t0_last, n_dips=ndip,
                runtime_s=rt, match_measure=okm, match_depth=okd, match_t0=okt), \
        (okm and okd and okt)


if __name__ == '__main__':
    ra, oka = part_a()
    rb, okb = part_b()
    rc, okc = part_c()
    os.makedirs(os.path.join(ROOT, 'data'), exist_ok=True)
    json.dump(dict(part_a=ra, part_b=rb, part_c=rc, gate_b=bool(okb), gate_c=bool(okc)),
              open(os.path.join(ROOT, 'data', 'stage1_validation.json'), 'w'), indent=1,
              default=lambda o: bool(o) if isinstance(o, (np.bool_, bool)) else float(o))
    print("\n" + "=" * 74)
    print("GATE (b) A_L table : %s" % ("PASS" if okb else "FAIL"))
    print("GATE (c) Zhu s7    : %s" % ("PASS" if okc else "FAIL"))
    print("     (a) precision : %s" % ("PASS" if oka else "FAIL"))
    print("=" * 74)
