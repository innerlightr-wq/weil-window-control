#!/usr/bin/env python3
"""Stage 4 -- minima catalogue (preregistered design).

For each distinct symbol: every local minimum of Psi_L with Psi_L < 1, with
depth, curvature, phase vector, and per-prime alignment share.  Statistics:
spacing distribution, Fano clustering, depth tail, and which prime subsets
carry the deepest 1%.  All of it repeated for N2 surrogates and for the N1
torus model, which is what P3 and P4 are stated against.
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from symbol import Symbol, STAGE2_L, P_accurate, reduced_arg
from scan import scan
from stage3_nulls import surrogate

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BETA_CAT = 1.0                 # the preregistered catalogue threshold: Psi_L < 1
DT = 0.05                      # sampling step inside a dip (shortest period ~ 2pi/2L)
N_SUR = 20
SPREAD = 0.3


# ---------------------------------------------------------------- blocks ---
def blocks(sym):
    """{base prime p: (freqs, weights)} -- the prime-power block of p."""
    out = {}
    for p, f, w in zip(sym.prime, sym.freqs, sym.weights):
        out.setdefault(int(p), [[], []])
        out[int(p)][0].append(float(f)); out[int(p)][1].append(float(w))
    return {p: (np.array(fs), np.array(ws)) for p, (fs, ws) in out.items()}


def block_value(fs, ws, t):
    acc = np.zeros(np.shape(t))
    for f, w in zip(fs, ws):
        acc += w * np.cos(reduced_arg(t, f))
    return acc


def shares(sym, t, blk):
    """alignment share s_p(t) = +block_p(t)/A_p in [-1,1].

    Psi = H - P, so a DIP is where P is LARGE.  s_p = 1 therefore means the p-block
    is maximally aligned to deepen the dip.  (Sign matters: getting it backwards
    makes the real catalogue anti-correlate with the N1 prediction.)
    """
    return {p: block_value(fs, ws, t) / ws.sum() for p, (fs, ws) in blk.items()}


# --------------------------------------------------------------- minima ----
def local_minima(sym, dips, dt=DT, n_gold=60):
    """Every interior local minimum of Psi on the given dip intervals."""
    cand_lo, cand_hi = [], []
    for a, b in dips:
        n = max(5, int(np.ceil((b - a) / dt)) + 1)
        ts = np.linspace(a, b, n)
        v = sym.H(ts) - P_accurate(sym, ts)
        i = np.nonzero((v[1:-1] < v[:-2]) & (v[1:-1] <= v[2:]))[0] + 1
        if i.size == 0:                       # flat/narrow: bracket the whole cell
            cand_lo.append(a); cand_hi.append(b)
        else:
            cand_lo.extend(ts[i - 1]); cand_hi.extend(ts[i + 1])
    lo = np.array(cand_lo); hi = np.array(cand_hi)
    if lo.size == 0:
        return np.empty(0), np.empty(0)
    g = (np.sqrt(5.0) - 1.0) / 2.0            # vectorised golden section
    for _ in range(n_gold):
        c = hi - g * (hi - lo); d = lo + g * (hi - lo)
        fc = sym.H(c) - P_accurate(sym, c)
        fd = sym.H(d) - P_accurate(sym, d)
        m = fc < fd
        hi = np.where(m, d, hi); lo = np.where(m, lo, c)
    t = 0.5 * (lo + hi)
    return t, sym.H(t) - P_accurate(sym, t)


def curvature(sym, t):
    """Psi'' = H'' + sum w f^2 cos(t f).  H'' = O(1/t^2), by finite difference."""
    acc = np.zeros(np.shape(t))
    for f, w in zip(sym.freqs, sym.weights):
        acc += w * f * f * np.cos(reduced_arg(t, f))
    h = 0.05
    Hpp = (sym.H(t + h) - 2 * sym.H(t) + sym.H(t - h)) / (h * h)
    return acc + Hpp


# ----------------------------------------------------------- statistics ----
def fano(t, lo, hi, n_win=200):
    """Var/Mean of minima counts in n_win equal windows -- 1 = Poisson."""
    if t.size < 20:
        return float('nan')
    edges = np.linspace(lo, hi, n_win + 1)
    c = np.histogram(t, bins=edges)[0].astype(float)
    return float(c.var() / c.mean()) if c.mean() > 0 else float('nan')


def catalogue(sym, beta=BETA_CAT, keep_rows=False):
    sc = scan(sym, beta)
    t, d = local_minima(sym, sc['dips'])
    if t.size == 0:
        return None
    o = np.argsort(t); t, d = t[o], d[o]
    blk = blocks(sym)
    sh = shares(sym, t, blk)
    primes = sorted(blk)
    A23 = sum(blk[p][1].sum() for p in primes if p in (2, 3))
    s23 = sum(block_value(*blk[p], t) for p in primes if p in (2, 3)) / A23
    sp = np.diff(t)
    k = max(1, int(round(0.01 * t.size)))
    deep = np.argsort(d)[:k]                    # deepest 1%
    out = dict(
        L=sym.L, beta=beta, primes=primes, A=sym.A, n_minima=int(t.size),
        t_lo=float(t.min()), t_hi=float(t.max()),
        spacing_mean=float(sp.mean()) if sp.size else float('nan'),
        spacing_cv=float(sp.std() / sp.mean()) if sp.size else float('nan'),
        fano=fano(t, 2 * np.pi, sym.envelope_end(beta)),
        depth_min=float(d.min()),
        depth_q=[float(np.quantile(d, q)) for q in (0.01, 0.05, 0.25, 0.5)],
        frac_below_0=float((d < 0).mean()),
        share_all={str(p): float(sh[p].mean()) for p in primes},
        share_deep={str(p): float(sh[p][deep].mean()) for p in primes},
        s23_all=float(s23.mean()), s23_deep=float(s23[deep].mean()),
        weight={str(p): float(blk[p][1].sum()) for p in primes},
    )
    # which prime subsets carry the deepest 1%: primes aligned past 1/2
    carriers = {}
    for j in deep:
        key = ','.join(str(p) for p in primes if sh[p][j] > 0.5)
        carriers[key] = carriers.get(key, 0) + 1
    out['deep_carriers'] = sorted(carriers.items(), key=lambda kv: -kv[1])[:6]
    if keep_rows:
        out['_rows'] = (t, d, curvature(sym, t), sh, s23)
    return out


def n1_catalogue(sym, n=400_000, seed=0):
    """N1: the same statistics with the phases i.i.d. uniform on the torus."""
    rng = np.random.default_rng(seed)
    blk = blocks(sym); primes = sorted(blk)
    th = {p: rng.uniform(0, 2 * np.pi, n) for p in primes}
    bv, sh = {}, {}
    for p in primes:
        fs, ws = blk[p]
        k = np.round(fs / np.log(p)).astype(int)        # the exponent of p
        v = np.zeros(n)
        for kk, w in zip(k, ws):
            v += w * np.cos(kk * th[p])
        bv[p] = v; sh[p] = v / ws.sum()
    P = sum(bv.values())
    A23 = sum(blk[p][1].sum() for p in primes if p in (2, 3))
    s23 = sum(bv[p] for p in primes if p in (2, 3)) / A23
    # "deep" in N1 = the 1% LARGEST P (Psi = H - P), matched to the deepest-1% of minima
    kdeep = max(1, n // 100)
    deep = np.argsort(P)[-kdeep:]
    return dict(share_all={str(p): float(sh[p].mean()) for p in primes},
                share_deep={str(p): float(sh[p][deep].mean()) for p in primes},
                s23_all=float(s23.mean()), s23_deep=float(s23[deep].mean()))


# ----------------------------------------------------------------- main ----
def main():
    seen, rows = set(), []
    for L in STAGE2_L:
        sym = Symbol(L)
        key = tuple(sym.n.tolist())
        if key in seen:
            continue
        seen.add(key)
        t0 = time.time()
        real = catalogue(sym)
        n1 = n1_catalogue(sym)
        sur = []
        rng = np.random.default_rng(1234 + int(100 * L))
        for _ in range(N_SUR):
            s2 = surrogate(sym, rng, spread=SPREAD)
            s2.prime = sym.prime.copy(); s2.expo = sym.expo.copy()
            c = catalogue(s2)
            if c:
                sur.append(c)
        agg = {}
        for k in ('n_minima', 'spacing_cv', 'fano', 'depth_min', 'frac_below_0',
                  's23_all', 's23_deep'):
            v = np.array([c[k] for c in sur], dtype=float)
            agg[k] = [float(np.quantile(v, 0.025)), float(np.median(v)),
                      float(np.quantile(v, 0.975))]
        agg['share_deep'] = {p: [float(np.quantile([c['share_deep'][p] for c in sur], q))
                                 for q in (0.025, 0.5, 0.975)] for p in real['share_deep']}
        rows.append(dict(real=real, n1=n1, surrogate=agg, runtime_s=time.time() - t0))
        r = rows[-1]
        print(f"L={L:4}  primes={real['primes']}  n_min={real['n_minima']:6d} "
              f"cv={real['spacing_cv']:.3f} fano={real['fano']:7.3f} "
              f"min={real['depth_min']:.3f}  s23_deep real={real['s23_deep']:+.3f} "
              f"N1={n1['s23_deep']:+.3f} N2=[{agg['s23_deep'][0]:+.3f},{agg['s23_deep'][2]:+.3f}]"
              f"  ({r['runtime_s']:.0f}s)", flush=True)
    with open(os.path.join(ROOT, 'data', 'stage4_minima.json'), 'w') as f:
        json.dump(rows, f, indent=1, default=lambda o: o.tolist())
    print('\nwrote data/stage4_minima.json')


if __name__ == '__main__':
    main()
