"""Build the frozen witness cells covering [delta_new, 1/4].

Deterministic schedule: cell widths start at 1/4000 and grow by a factor 7/5 up
to a cap of 1/400; any cell whose interval enclosure fails is halved and retried.
Each cell carries ONE fixed rational witness, certified over the WHOLE closed
cell by interval arithmetic in delta. The witness is the rationalised
minimum-eigenvalue direction of M at the cell's left endpoint -- the direction
that maximises the UNIT-L2 margin, which is what the error allowance is
compared against (the H-weighted maximiser x_* = H^{-1}b maximises a different
ratio and gives a much smaller absolute margin; see VERDICT.md).
The previous cell's witness is tried first, and kept when it still certifies.
"""
import sys, os, json, resource, time
sys.path.insert(0, os.path.dirname(__file__))
resource.setrlimit(resource.RLIMIT_AS, (4 << 30, 4 << 30))
from fractions import Fraction as Q
import weil_core as C
import tail_sharp as TS
import ratlin as R

NDIM = C.NDIM
def mid(iv): return (iv.lo + iv.hi) / 2
G = [[mid(e) for e in row] for row in C.G_core()]

def M_at(delta):
    a, b = C.ab_vectors(Q(delta))
    a = [mid(t) for t in a]; b = [mid(t) for t in b]
    return [[G[i][j] + 4*a[i]*a[j] - 4*b[i]*b[j] for j in range(NDIM)] for i in range(NDIM)]

def mineig_dir(delta, iters=50):
    Mm = M_at(delta)
    lo, hi = R.min_eig_bracket(Mm, iters=iters)
    lam = (lo + hi) / 2
    shift = lam - abs(lam)/1000 - Q(1, 10**9)
    B = [[Mm[i][j] - (shift if i == j else 0) for j in range(NDIM)] for i in range(NDIM)]
    v = [Q(1)]*NDIM
    for _ in range(3):
        w = R.solve(B, v)
        if w is None: break
        s = max(abs(t) for t in w)
        v = [t/s for t in w]
    return [Q(round(t * 10**12), 10**12) for t in v], lam

def cert(x, lo, hi):
    d = TS.Q_zeta_sharp(x)
    Ad, Bd = C.quartet(x, C.Iv(Q(lo), Q(hi)))
    return d['Q'] + C.Iv(4) * Ad.sq() - C.Iv(4) * Bd.sq(), d

DELTA_NEW = Q(9, 100)
RIGHT = Q(1, 4)
WMIN, WMAX, GROW = Q(1, 4000), Q(1, 400), Q(7, 5)

cells = []
lo = DELTA_NEW
width = WMIN
prev_x = None
t0 = time.time()
print(f'target [{DELTA_NEW}, {RIGHT}]; widths {WMIN}..{WMAX}, growth {GROW}')
while lo < RIGHT:
    hi = min(lo + width, RIGHT)
    xr = prev_x
    qa = None
    if xr is not None:
        qa, d = cert(xr, lo, hi)
        if qa.hi >= 0:
            xr = None
    if xr is None:
        xr, lam = mineig_dir(lo)
        qa, d = cert(xr, lo, hi)
    else:
        lam = None
    shrinks = 0
    while qa.hi >= 0 and shrinks < 8:
        hi = lo + (hi - lo) / 2
        qa, d = cert(xr, lo, hi)
        shrinks += 1
    if qa.hi >= 0:
        print(f'  STOP: cannot certify from {float(lo):.8f}'); break
    nsq = sum(t*t for t in xr)
    cells.append(dict(lo=str(lo), hi=str(hi), x=[str(t) for t in xr],
                      Qadd_hi=str(qa.hi), Qadd_lo=str(qa.lo),
                      unitL2_margin=str(-qa.hi/nsq), normsq=str(nsq),
                      tail_halfwidth=str((d['tail_sharp'].hi - d['tail_sharp'].lo)/2),
                      lam_min_mid=(str(lam) if lam is not None else None),
                      reused_witness=(lam is None), shrinks=shrinks))
    if len(cells) % 10 == 1 or hi == RIGHT:
        print(f'  [{float(lo):.7f},{float(hi):.7f}] Qadd_hi={float(qa.hi):.4e} '
              f'margin/tail={float(-qa.hi/(d["tail_sharp"].hi-d["tail_sharp"].lo)*2):.1f} '
              f'reuse={lam is None} shrinks={shrinks}', flush=True)
    prev_x = xr
    lo = hi
    width = min(width * GROW, WMAX) if shrinks == 0 else max((hi - lo) if hi > lo else WMIN, WMIN)
    width = min(width, WMAX)
    if width < WMIN: width = WMIN

nw = len(set(tuple(c['x']) for c in cells))
blob = dict(delta_new=str(DELTA_NEW), right=str(RIGHT), cells=cells,
            n_cells=len(cells), n_distinct_witnesses=nw,
            M=C.M, L=str(C.L), gamma=str(C.GAMMA), ndim=NDIM,
            tail='sharp O(1/M^2), scripts/tail_sharp.py')
json.dump(blob, open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'results', 'certificate_new.json'), 'w'), indent=1)
print(f'\n{len(cells)} cells, {nw} distinct witnesses, built in {time.time()-t0:.1f}s')
print(f'coverage: {cells[0]["lo"]} .. {cells[-1]["hi"]}  complete = {Q(cells[-1]["hi"]) == RIGHT}')
print(f'worst unit-L2 margin = {min(float(c["unitL2_margin"]) for c in cells):.4e}')
print(f'peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024:.1f} MiB')
