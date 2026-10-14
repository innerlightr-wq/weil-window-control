#!/usr/bin/env python3
"""VERIFIER for the new certificate. Clean process, frozen inputs only.

Reads results/certificate_new.json (frozen rational witnesses and cells),
regenerates every enclosure from the mathematical input, and ASSERTS:

  1. mechanical coverage: cells are closed, consecutive, gap-free, and their
     union is exactly [delta_new, 1/4];
  2. for each cell, an enclosure of Q_add,delta(f_x) valid for EVERY delta in
     that closed cell, with a strictly negative upper bound;
  3. the join with the inherited certificate at delta = 1/4;
  4. the inherited point bound at delta = 2/5 for the published frozen vector,
     recomputed (not read);
  5. the sharper archimedean tail is CONTAINED in the published O(1/M) bound.

Exits nonzero on any failure. Relies on no precomputed verdict: the only cached
object is the x-INDEPENDENT arch sum, whose cache key carries the truncation,
the interval grid, the basis, the parameters and the source hash, and which is
recomputed if the key does not match.
"""
import sys, os, json, resource
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'scripts'))
resource.setrlimit(resource.RLIMIT_AS, (4 << 30, 4 << 30))
from fractions import Fraction as Q
import weil_core as C
import tail_sharp as TS

HERE = os.path.dirname(os.path.abspath(__file__))
fails = []
def chk(c, m):
    print(('PASS  ' if c else 'FAIL  ') + m)
    if not c: fails.append(m)

blob = json.load(open(os.path.join(HERE, 'results', 'certificate_new.json')))
cells = blob['cells']
dnew, right = Q(blob['delta_new']), Q(blob['right'])
print(f'new certificate: {len(cells)} cells covering [{dnew}, {right}]')
print(f'L = {C.L}, gamma = {C.GAMMA}, {C.NDIM} EVEN_D modes, M = {C.M}')
print()

# ---- 1. mechanical coverage ------------------------------------------------
ok = Q(cells[0]['lo']) == dnew and Q(cells[-1]['hi']) == right
for i in range(len(cells) - 1):
    if Q(cells[i]['hi']) != Q(cells[i+1]['lo']):
        ok = False
chk(ok, f'1. coverage: closed cells are consecutive with no gap, union = [{dnew}, {right}]')
chk(all(Q(c['lo']) < Q(c['hi']) for c in cells), '1. every cell is nondegenerate')

# ---- 2. per-cell certification --------------------------------------------
worst = None
nbad = 0
for i, c in enumerate(cells):
    lo, hi = Q(c['lo']), Q(c['hi'])
    x = [Q(t) for t in c['x']]
    d = TS.Q_zeta_sharp(x)
    Ad, Bd = C.quartet(x, C.Iv(lo, hi))
    qa = d['Q'] + C.Iv(4) * Ad.sq() - C.Iv(4) * Bd.sq()
    neg = qa.hi < 0
    if not neg:
        nbad += 1
        chk(False, f'2. cell {i} [{float(lo):.6f},{float(hi):.6f}]: Q_add upper bound '
                   f'{float(qa.hi):.6e} is NOT negative')
    nsq = sum(t*t for t in x)
    margin = -qa.hi / nsq
    if worst is None or margin < worst[0]:
        worst = (margin, i, lo, hi)
chk(nbad == 0,
    f'2. all {len(cells)} cells: Q_add,delta(f_x) < 0 for every delta in the closed '
    f'cell ({nbad} cells failed)')
if worst:
    print(f'      worst unit-L2 margin {float(worst[0]):.6e} on cell {worst[1]} '
          f'[{float(worst[2]):.6f}, {float(worst[3]):.6f}]')

# ---- 3. the join at 1/4 ----------------------------------------------------
xlast = [Q(t) for t in cells[-1]['x']]
d = TS.Q_zeta_sharp(xlast)
Ad, Bd = C.quartet(xlast, C.Iv(right, right))
qj = d['Q'] + C.Iv(4) * Ad.sq() - C.Iv(4) * Bd.sq()
chk(qj.hi < 0, f'3. join: the last cell\'s witness is still negative AT delta = {right} '
               f'(upper bound {float(qj.hi):.6e})')
XF = C.XFROZEN
qinh, _, _, _ = C.Q_add(XF, right)
chk(qinh.hi < 0, f'3. join: the INHERITED frozen vector is negative at delta = {right} '
                 f'(upper bound {float(qinh.hi):.6e}), so the two certified sets overlap at {right}')

# ---- 4. the inherited point bound, recomputed ------------------------------
qpt, _, _, dz = C.Q_add(XF, Q(2, 5))
chk(qpt.hi <= Q(-145107657498227, 10**16),
    f'4. inherited point bound recomputed: Q_add,2/5(frozen) <= {float(qpt.hi):.16f} '
    f'<= -0.0145107657498227')
chk(dz['Q'].lo > 0, f'4. Q_zeta(frozen) > 0, recomputed: lo = {float(dz["Q"].lo):.9e}')

# ---- 5. the sharper tail is inside the published bound ---------------------
dpub = C.Q_zeta(XF); dsh = TS.Q_zeta_sharp(XF)
chk(dpub['Q'].lo <= dsh['Q'].lo and dsh['Q'].hi <= dpub['Q'].hi,
    f'5. sharper tail: [{float(dsh["Q"].lo):.15f}, {float(dsh["Q"].hi):.15f}] is inside the '
    f'published [{float(dpub["Q"].lo):.15f}, {float(dpub["Q"].hi):.15f}]')
hw = float((dsh['tail_sharp'].hi - dsh['tail_sharp'].lo) / 2)
chk(hw < 1e-8, f'5. sharper tail allowance half-width {hw:.4e} (published 4.935034e-05)')

print()
print('NEW CERTIFICATE:', 'ALL ASSERTIONS PASS' if not fails else f'{len(fails)} FAILURES')
print(f'peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024:.1f} MiB')
sys.exit(1 if fails else 0)
