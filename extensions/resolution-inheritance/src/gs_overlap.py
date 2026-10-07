#!/usr/bin/env python3
"""Does the inner part of the L'-ground state converge to the L-ground state?
If H_{L'} = H_L (+) H_new were reducing, the ground state would lie wholly in one
block and the overlap would be exactly 1."""
import sys
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)                                            # this note's wform.py
sys.path.insert(0, os.path.join(_HERE, '..', '..', '..', 'src'))     # the project's modules
_DATA = os.path.join(_HERE, '..', 'data')
import mpmath as mp, zeta_window as zw, wform as W, weilform as wf

print("  %-5s %-5s %-5s %10s %10s %14s" % ("L", "L'", "N", "||v_new||", "overlap", "1-overlap"))
for N in (4, 6, 8, 12, 16):
    dps = 30; mp.mp.dps = dps + 12
    L, Lp = mp.mpf('0.5'), mp.mpf('0.8')
    bi, bn = W.inner_basis(L, N, 'even'), W.annulus_basis(L, Lp, N, 'even')
    b = bi + bn; nrm = [mp.sqrt(W.l2(f, f)) for f in b]
    M = W.scale_matrix(W.Q_matrix(b, 2 * Lp, zw.von_mangoldt_terms(2 * Lp)), nrm)
    QL = W.scale_matrix(W.Q_matrix(bi, 2 * L, zw.von_mangoldt_terms(2 * L)), nrm[:N])
    mp.mp.dps = dps
    _, vecs = wf.eigsym(M); v = vecs[0]
    _, vecL = wf.eigsym(QL); u = vecL[0]
    wI = mp.sqrt(sum(v[i] ** 2 for i in range(N)))
    wN = mp.sqrt(sum(v[j] ** 2 for j in range(N, 2 * N)))
    ov = abs(sum(v[i] * u[i] for i in range(N))) / wI
    print("  %-5s %-5s %-5d %10s %10s %14s" % (L, Lp, N, mp.nstr(wN, 4), mp.nstr(ov, 6), mp.nstr(1 - ov, 4)))
