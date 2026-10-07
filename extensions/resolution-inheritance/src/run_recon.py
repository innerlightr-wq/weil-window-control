import sys, time, json
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)                                            # this note's wform.py
sys.path.insert(0, os.path.join(_HERE, '..', '..', '..', 'src'))     # the project's modules
_DATA = os.path.join(_HERE, '..', 'data')
import mpmath as mp
import zeta_window as zw
import weilform as wf
import wform as W

def build(L, Lp, Nin, Nnew, dps, kind):
    """Q_{L'} in an orthonormal basis adapted to H_{L'} = H_L (+) H_new."""
    mp.mp.dps = dps + 12
    L, Lp = mp.mpf(L), mp.mpf(Lp)
    if kind == 'even':
        bi, bn = W.inner_basis(L, Nin, 'even'), W.annulus_basis(L, Lp, Nnew, 'even')
    else:
        bi, bn = W.inner_basis(L, Nin, 'even_d'), W.sine_annulus_basis(L, Lp, Nnew)
    b = bi + bn
    nrm = [mp.sqrt(W.l2(f, f)) for f in b]
    prim = zw.von_mangoldt_terms(2 * Lp)
    M = W.scale_matrix(W.Q_matrix(b, 2 * Lp, prim), nrm)
    mp.mp.dps = dps
    return M, bi, bn, nrm

def report(L, Lp, Nin, Nnew, dps, kind):
    t = time.time()
    M, bi, bn, nrm = build(L, Lp, Nin, Nnew, dps, kind)
    n = Nin + Nnew
    old = range(Nin); new = range(Nin, n)
    # Q_L on its own, in the SAME inner basis (zero-extended): the exact comparison
    mp.mp.dps = dps + 12
    Lm = mp.mpf(L)
    primL = zw.von_mangoldt_terms(2 * Lm)
    QL = W.scale_matrix(W.Q_matrix(bi, 2 * Lm, primL), nrm[:Nin])
    mp.mp.dps = dps
    e_old = W.fro(mp.matrix([[M[i,j]-QL[i,j] for j in old] for i in old])) / W.fro(QL)
    e_cross = W.fro(M, old, new) / W.fro(M)      # single block, as defined in the brief
    # spectra
    vL, _ = wf.eigsym(QL)
    vP, _ = wf.eigsym(M)
    # for each eigenvalue of Q_L, nearest in Q_{L'}, relative
    worst = max(min(abs(a-b)/max(abs(a),mp.mpf(1)) for b in vP) for a in vL)
    return dict(L=L, Lp=Lp, Nin=Nin, Nnew=Nnew, dps=dps, kind=kind,
                e_old=e_old, e_cross=e_cross, worst_eig=worst,
                lmin_L=vL[0], lmin_Lp=vP[0], t=time.time()-t)

CASES = [
    # (L, L', Nin, Nnew, dps, kind)
    (0.5, 0.8, 4, 4, 30, 'even'),
    (0.5, 0.8, 6, 6, 30, 'even'),
    (0.5, 0.8, 8, 8, 30, 'even'),
    (0.5, 0.8, 6, 6, 50, 'even'),
    (0.5, 0.8, 6, 6, 30, 'even_d'),
    (0.8, 1.0, 6, 6, 30, 'even'),
    (0.8, 1.6, 6, 6, 30, 'even'),
    (1.0, 1.2, 6, 6, 30, 'even'),
]
rows = []
print("  %-5s %-5s %-4s %-5s %-4s %-7s %12s %12s %12s %10s" %
      ("L","L'","Nin","Nnew","dps","basis","E_old","E_cross","worst eig","t"))
for c in CASES:
    r = report(*c)
    rows.append({k: (mp.nstr(v, 6) if isinstance(v, mp.mpf) else v) for k, v in r.items()})
    print("  %-5s %-5s %-4d %-5d %-4d %-7s %12s %12s %12s %9.0fs" %
          (r['L'], r['Lp'], r['Nin'], r['Nnew'], r['dps'], r['kind'],
           mp.nstr(r['e_old'],4), mp.nstr(r['e_cross'],4), mp.nstr(r['worst_eig'],4), r['t']))
json.dump(rows, open(os.path.join(_DATA,'recon_main.json'),'w'), indent=1)
