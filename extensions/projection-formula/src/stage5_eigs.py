#!/usr/bin/env python3
"""Cache the FULL spectrum of Q_0 for each cell (stage2_zeta.json stored only the 12
smallest).  Needed for a defensible K3 alignment audit: the largest eigenvalues give
delta_k = sqrt(lam_k/4K) of order 1, which is inside the sweep's top decade."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'src'))
import mpmath as mp
from zeta_window import build_matrix
from stage3e_certify import move_block
from weilform import eigsym
from stage2_zeta import CELLS, SEC, K_closed
out = []
for Lx, N, dps in CELLS:
    mp.mp.dps = dps
    g1 = mp.im(mp.zetazero(1)); L = mp.mpf(Lx)
    Q0 = build_matrix(L, N, dps, SEC) + move_block(L, range(N), g1, mp.mpf(0), SEC)
    vals, _ = eigsym(Q0)
    K = K_closed(L, g1)
    out.append(dict(L=Lx, N=N, dps=dps, fourK=mp.nstr(4*K, 16),
                    eigs_all=[mp.nstr(v, 16) for v in vals]))
    print('  L=%-4s N=%-3d full spectrum cached; lam range %s .. %s'
          % (Lx, N, mp.nstr(vals[0], 6), mp.nstr(vals[-1], 6)), flush=True)
json.dump(out, open('data/stage5_eigs_full.json', 'w'), indent=1)
print('  wrote data/stage5_eigs_full.json')
