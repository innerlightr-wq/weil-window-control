r"""Certified branch-and-bound scan for D_L(beta) = {t >= 2pi : Psi_L(t) < beta}.

Rigour.  A cell [a, a+h] is discarded ONLY when
    H(a) - min(A_L, P(a+h/2) + Lip_P h/2)  >=  beta,
a valid lower bound on min Psi over the cell (H increasing, notes/monotonicity.md; and
P <= A_L unconditionally).  The surviving union therefore always CONTAINS the true dip
set, so sup(surviving) is a rigorous upper bound for t0 and the last sampled point with
Psi < beta is a rigorous lower bound.

Three passes:
  A  one chunked uniform sweep at h0 over [2pi, 2pi e^{A+beta}]  -- kills the bulk;
  B  branch-and-bound on the survivors down to h_min;
  C  vectorised bisection for the actual {Psi = beta} crossings inside the survivors.
"""
import numpy as np
from kernel import cell_lower_bound, cell_bounds, psi_many


def _merge(lo, hi, tol=0.0):
    if lo.size == 0:
        return np.empty((0, 2))
    o = np.argsort(lo); lo, hi = lo[o], hi[o]
    keep_lo = [lo[0]]; keep_hi = [hi[0]]
    for i in range(1, lo.size):
        if lo[i] <= keep_hi[-1] + tol:
            if hi[i] > keep_hi[-1]:
                keep_hi[-1] = hi[i]
        else:
            keep_lo.append(lo[i]); keep_hi.append(hi[i])
    return np.column_stack([keep_lo, keep_hi])


def pass_A(freqs, weights, A, lipP, beta, lo, hi, h0, chunk_cells=8_000_000):
    """Chunked uniform sweep. Returns surviving cell left-endpoints."""
    n = int(np.ceil((hi - lo) / h0))
    out = []
    for i0 in range(0, n, chunk_cells):
        i1 = min(i0 + chunk_cells, n)
        a = lo + h0 * np.arange(i0, i1, dtype=np.float64)
        lb = cell_lower_bound(a, h0, freqs, weights, A, lipP)
        s = a[lb < beta]
        if s.size:
            out.append(s)
    return np.concatenate(out) if out else np.empty(0)


def pass_B(freqs, weights, A, lipP, beta, a, h, h_min):
    """Two-sided B&B. Returns (accepted whole-dip cells, width), (boundary cells, width)."""
    acc_lo, acc_hi = [], []
    while True:
        lo, up = cell_bounds(a, h, freqs, weights, A, lipP)
        inside = up < beta                       # certified entirely below beta
        if inside.any():
            acc_lo.append(a[inside]); acc_hi.append(a[inside] + h)
        a = a[(lo < beta) & ~inside]             # undecided only
        if a.size == 0 or h <= h_min:
            break
        a = np.concatenate((a, a + 0.5 * h)); a.sort(); h *= 0.5
    inside_iv = (_merge(np.concatenate(acc_lo), np.concatenate(acc_hi), tol=1e-12)
                 if acc_lo else np.empty((0, 2)))
    return inside_iv, a, h


def pass_C(freqs, weights, beta, iv, n_sample=512, n_bisect=52):
    """Locate the {Psi = beta} crossings. All intervals and all crossings are bisected
    in ONE batch -- per-interval numba calls dominated the runtime otherwise."""
    if iv.size == 0:
        return []
    # sample every surviving interval on a common grid resolution
    ts, owner = [], []
    for j, (a, b) in enumerate(iv):
        t = np.linspace(a, b, n_sample)
        ts.append(t); owner.append(np.full(n_sample, j))
    t = np.ascontiguousarray(np.concatenate(ts))
    own = np.concatenate(owner)
    v = psi_many(t, freqs, weights) - beta
    neg = v < 0
    # crossings only inside an interval (not across the join between two intervals)
    k = np.flatnonzero((neg[:-1] != neg[1:]) & (own[:-1] == own[1:]))
    if k.size:
        x0 = t[k].copy(); x1 = t[k + 1].copy(); f0 = v[k].copy()
        for _ in range(n_bisect):
            xm = np.ascontiguousarray(0.5 * (x0 + x1))
            fm = psi_many(xm, freqs, weights) - beta
            same = (fm < 0) == (f0 < 0)
            x0 = np.where(same, xm, x0); f0 = np.where(same, fm, f0)
            x1 = np.where(same, x1, xm)
        cross = 0.5 * (x0 + x1)
        cross_owner = own[k]
    else:
        cross = np.empty(0); cross_owner = np.empty(0, dtype=int)
    dips = []
    for j, (a, b) in enumerate(iv):
        sel = own == j
        vj = v[sel]
        pts = sorted(cross[cross_owner == j].tolist())
        edges = ([float(a)] if vj[0] < 0 else []) + pts + ([float(b)] if vj[-1] < 0 else [])
        for i in range(0, len(edges) - 1, 2):
            if edges[i + 1] > edges[i]:
                dips.append((edges[i], edges[i + 1]))
    return dips


def scan(sym, beta, h0=0.5, h_min=1e-9, verbose=False):
    freqs = np.ascontiguousarray(sym.freqs)
    weights = np.ascontiguousarray(sym.weights)
    A, lipP = sym.A, sym.lipschitz_P()
    lo, hi = 2 * np.pi, sym.envelope_end(beta)
    a = pass_A(freqs, weights, A, lipP, beta, lo, hi, h0)
    nA = a.size
    inside_iv, bnd, hb = pass_B(freqs, weights, A, lipP, beta, a, h0, h_min)
    # the surviving union = certified-inside intervals + unresolved boundary cells
    pieces = []
    if inside_iv.size:
        pieces.append(inside_iv)
    if bnd.size:
        pieces.append(np.column_stack([bnd, bnd + hb]))
    iv = _merge(np.concatenate([p[:, 0] for p in pieces]),
                np.concatenate([p[:, 1] for p in pieces]),
                tol=1e-12) if pieces else np.empty((0, 2))
    dips = pass_C(freqs, weights, beta, iv)
    # A single dip can be delivered as several touching pieces, because the
    # certified-inside cells and the refined boundary cells are separate intervals.
    # Merge pieces whose gap is below MERGE_TOL. Verified safe: the smallest genuine
    # gap between distinct dips is O(1) (0.97 at L=1.4), twelve orders above the tol.
    MERGE_TOL = 1e-9
    dips.sort()
    if dips:
        m = [list(dips[0])]
        for x, b in dips[1:]:
            if x <= m[-1][1] + MERGE_TOL:
                m[-1][1] = max(m[-1][1], b)
            else:
                m.append([x, b])
        dips = [(a_, b_) for a_, b_ in m]
    measure = float(sum(b - x for x, b in dips))
    # rigorous bracket: the dip set is contained in iv, and contains inside_iv
    meas_lo = float((inside_iv[:, 1] - inside_iv[:, 0]).sum()) if inside_iv.size else 0.0
    meas_hi = float((iv[:, 1] - iv[:, 0]).sum()) if iv.size else 0.0
    return dict(L=sym.L, beta=beta, lo=lo, hi=hi, h0=h0, h_boundary=hb,
                n_passA_cells=int(nA), n_surviving=int(iv.shape[0]),
                measure_lower=meas_lo, measure_upper=meas_hi,
                n_dips=len(dips), dip_measure=measure,
                t0_lower=max((b for _, b in dips), default=float('nan')),
                t0_upper=float(iv[:, 1].max()) if iv.size else float('nan'),
                dips=dips)
