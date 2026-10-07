#!/usr/bin/env python3
"""Stage 4 (H4, confirmation only) and Stage 5 (combined formula).

Stage 4: lam*(L) = lam_min of the unperturbed Weil form, against Zhu's Landau-Widom law
         -ln lam*(L) ~ 2 pi^2 N(T*)/ln N(T*),  T* = 2 pi e^{2L}.
         Window dictionary: Connes's multiplicative window [lam^-1, lam] has additive
         log-length 2 log lam; our f is supported in [-L,L], additive length 2L, so
         log lam = L  (T1, bookkeeping).
Stage 5: lam_min(Q_delta) ~ lam_floor - 4K cos^2(theta_N) delta^2, across L and delta.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'src'))
import mpmath as mp
from zeta_window import build_matrix
from weilform import eigsym
from stage2_zeta import CELLS, SEC

def Nz(T):
    T = mp.mpf(T)
    return T/(2*mp.pi)*mp.log(T/(2*mp.pi*mp.e)) + mp.mpf(7)/8

print('STAGE 4 (confirmation only; prior art: Connes "Letter" Fig. 1, Zhu Landau-Widom)\n')
print('  window dictionary: Connes [lam^-1,lam], additive length 2 log lam; ours [-L,L],')
print('  additive length 2L  =>  log lam = L.  (No factor of 2 discrepancy.)\n')
print('  %-5s %-4s %-18s %-13s %-16s %-13s %-9s' %
      ('L', 'N', 'lam*(L)', '-ln lam*', 'T*=2 pi e^{2L}', 'Zhu law', 'ratio'))
rows4 = []
for Lx, N, dps in CELLS:
    mp.mp.dps = dps
    L = mp.mpf(Lx)
    lam = eigsym(build_matrix(L, N, dps, SEC))[0][0]
    Ts = 2*mp.pi*mp.e**(2*L)
    NT = Nz(Ts)
    zhu = 2*mp.pi**2*NT/mp.log(NT)
    mln = -mp.log(lam)
    rows4.append(dict(L=Lx, N=N, lam=mp.nstr(lam, 8), mln=float(mln),
                      Ts=float(Ts), NT=float(NT), zhu=float(zhu), ratio=float(mln/zhu)))
    print('  %-5s %-4d %-18s %-13.5f %-16.4g %-13.5f %-9.4f' %
          (Lx, N, mp.nstr(lam, 8), float(mln), float(Ts), float(zhu), float(mln/zhu)))
json.dump(rows4, open('data/stage4_floor.json', 'w'), indent=1)
print('\n  T2 confirmation: -ln lam* tracks Zhu\'s law to a slowly varying factor;')
print('  nothing new claimed. The near-null eigenspace of Q_0 is the prolate projection')
print('  Pi(lam,k) of Connes-Consani-Moscovici (arXiv:2511.22755).')

print('\n\nSTAGE 5 -- combined formula  lam_min(Q_d) ~ lam_floor - 4K cos^2(theta_N) d^2\n')
cells = json.load(open('data/stage2_zeta.json'))
print('  %-5s %-9s %-15s %-15s %-11s %-12s %-s' %
      ('L', 'delta', 'lam_min meas.', 'combined pred', 'rel err', 'C d^2 / lam*', 'regime'))
rows5 = []
for c in cells:
    lam0 = float(mp.mpf(c['lam_min_Q0']))
    fourK = float(c['fourK'])
    for s in c['sweep']:
        d = float(s['delta']); meas = float(mp.mpf(s['lam_min']))
        pred = lam0 - fourK*s['cos2']*d**2
        rel = abs(meas-pred)/abs(meas) if meas != 0 else float('nan')
        ratio = (s['C_float']*d**2)/lam0 if lam0 > 0 else float('inf')
        reg = 'C d^2 >> lam*' if ratio > 10 else ('comparable' if ratio > 0.1 else 'C d^2 << lam*')
        rows5.append(dict(L=c['L'], delta=s['delta'], meas=meas, pred=pred, rel=rel,
                          ratio=ratio, regime=reg))
        print('  %-5s %-9s %-15.6e %-15.6e %-11.3e %-12.3e %-s' %
              (c['L'], s['delta'], meas, pred, rel, ratio, reg))
json.dump(rows5, open('data/stage5_combined.json', 'w'), indent=1)
good = [r for r in rows5 if r['ratio'] > 10]
bad = [r for r in rows5 if r['ratio'] <= 0.1]
import statistics
print('\n  regime C d^2 >> lam*  : %d rows, median rel err %.3e' %
      (len(good), statistics.median([r['rel'] for r in good])))
if bad:
    print('  regime C d^2 << lam*  : %d rows, median rel err %.3e' %
          (len(bad), statistics.median([r['rel'] for r in bad])))
