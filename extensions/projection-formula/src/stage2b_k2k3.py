#!/usr/bin/env python3
"""Stage 2b: K2 at the paper's delta = 0.02, and K3 staircase alignment."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'src'))
import mpmath as mp
from zeta_window import build_matrix, EVEN_D
from stage3e_certify import move_block
from weilform import eigsym
from stage2_zeta import CELLS, SEC, K_closed, reps

PAPER = {'0.8': 0.3107975951, '1.0': 0.8530343312, '1.3': 2.498698441,
         '1.6': 4.921790776, '2.0': 10.08283301}
rows = []
print('K2 -- at the paper section 4.3 value delta = 0.02\n')
print('  %-5s %-13s %-13s %-13s %-9s %-10s %-10s %-10s' %
      ('L', 'C (ours)', 'C (paper 4.3)', 'H1 pred', 'cos^2', 'rel H1', 'rel compr', 'within 10%'))
for Lx, N, dps in CELLS:
    mp.mp.dps = dps
    g1 = mp.im(mp.zetazero(1)); L = mp.mpf(Lx)
    M = build_matrix(L, N, dps, SEC)
    Q0 = M + move_block(L, range(N), g1, mp.mpf(0), SEC)
    p, r, s = reps(L, N, g1, SEC)
    K = K_closed(L, g1); fourK = 4*K
    nr2 = sum(x**2 for x in r)
    vals, vecs = eigsym(Q0); V = [mp.matrix(v) for v in vecs]
    P = mp.zeros(N, N)
    for a in range(N):
        for b in range(N):
            P[a, b] = -4*(r[a]*r[b] + (p[a]*s[b] + p[b]*s[a])/2)
    d = mp.mpf('0.02')
    lam = eigsym(M + move_block(L, range(N), g1, d, SEC))[0][0]
    C = (vals[0]-lam)/d**2
    idx = [k for k in range(N) if vals[k] <= fourK*d**2]
    proj2 = sum(sum(V[k][i]*r[i] for i in range(N))**2 for k in idx) if idx else mp.mpf(0)
    cos2 = proj2/nr2; pred = fourK*cos2
    m = len(idx); G = mp.zeros(m, m)
    for ii, k in enumerate(idx):
        for jj, l in enumerate(idx):
            G[ii, jj] = (d**2)*sum(V[k][a]*P[a, b]*V[l][b] for a in range(N) for b in range(N))
            if ii == jj: G[ii, jj] += vals[k]
    Ccomp = (vals[0]-eigsym(G)[0][0])/d**2
    rel = abs(C-pred)/abs(C); relc = abs(C-Ccomp)/abs(C)
    rows.append(dict(L=Lx, C=float(C), paper=PAPER[Lx], pred=float(pred), cos2=float(cos2),
                     rel=float(rel), relc=float(relc), dimN=m,
                     rel_vs_paper=float(abs(C-PAPER[Lx])/PAPER[Lx])))
    print('  %-5s %-13.7g %-13.7g %-13.7g %-9.5f %-10.3e %-10.3e %-10s' %
          (Lx, float(C), PAPER[Lx], float(pred), float(cos2), float(rel), float(relc),
           'YES' if rel <= 0.10 else 'no'))
ok = sum(1 for r_ in rows if r_['rel'] <= 0.10)
print('\n  reconciliation with 4.3: max |C_ours - C_paper|/C_paper = %.2e' %
      max(r_['rel_vs_paper'] for r_ in rows))
print('  K2: H1 within 10%% at %d of 5 windows -> %s' % (ok, 'PASS' if ok >= 4 else 'FAIL'))
json.dump(rows, open('data/stage2b_k2.json', 'w'), indent=1)

print('\n\nK3 -- staircase: do C(delta) plateau edges align with eigenvalue scales of Q_0?\n')
cells = json.load(open('data/stage2_zeta.json'))
k3 = []
for c in cells:
    ev = [float(x) for x in c['eigs']]
    sw = c['sweep']
    print('  L=%s   eigenvalues of Q_0 (sqrt scale, delta_k := sqrt(lam_k/4K)):' % c['L'])
    fourK = float(c['fourK'])
    dk = [ (e/fourK)**0.5 for e in ev if e > 0 ]
    print('        ' + '  '.join('%.2e' % x for x in dk[:8]))
    prev = None; edges = []
    for s in sw:
        if prev is not None and s['dimN'] != prev['dimN']:
            edges.append((float(prev['delta']), float(s['delta']), prev['dimN'], s['dimN']))
        prev = s
    # a step in C is a change of more than 20% between consecutive decades
    csteps = []
    for i in range(1, len(sw)):
        a_, b_ = sw[i-1]['C_float'], sw[i]['C_float']
        if b_ > 0 and abs(a_-b_)/max(abs(a_), 1e-300) > 0.2:
            csteps.append((float(sw[i-1]['delta']), float(sw[i]['delta'])))
    # alignment: every C-step bracket should contain a delta_k within a factor 3
    aligned = 0
    for lo, hi in csteps:
        if any(min(lo, hi)/3 <= x <= max(lo, hi)*3 for x in dk): aligned += 1
    print('        C-steps: %d ; aligned with an eigenvalue scale (factor 3): %d' % (len(csteps), aligned))
    k3.append(dict(L=c['L'], steps=len(csteps), aligned=aligned))
tot, al = sum(x['steps'] for x in k3), sum(x['aligned'] for x in k3)
print('\n  K3: %d of %d C-steps align with an eigenvalue scale within a factor 3 -> %s'
      % (al, tot, 'PASS' if al >= 0.8*tot else 'FAIL'))
json.dump(k3, open('data/stage2b_k3.json', 'w'), indent=1)
