#!/usr/bin/env python3
"""Audit A: the denominator in cos^2, and what the "~3%" agreement actually refers to.

PREREGISTRATION fixed  cos^2 theta_N = ||Pi_N r||^2 / ||r||^2  with ||r||^2 = K(L,gamma)
(the continuum norm).  The IMPLEMENTATION (stage2_zeta.reps + run) instead divides by
    nr2 = sum_a r[a]^2 = ||Pi_{basis} r||^2,
the squared norm of the FINITE-BASIS projection of the representer.  Per the revision
brief we RETAIN the implemented definition and report the omitted tail separately:
    omitted tail fraction  =  1 - ||Pi_{basis} r||^2 / K.
Consequence: the H1 prediction as computed is
    pred = 4K * proj2/nr2 = 4*proj2 / (nr2/K),
inflated by 1/(nr2/K) relative to the consistent finite-basis prediction 4*proj2.
Both are tabulated.  No number is altered; the two normalisations are kept distinct.
"""
import json, statistics
cells = json.load(open('data/stage2_zeta.json'))
print('  omitted-tail fraction of the representer (1 - ||Pi_basis r||^2 / K)\n')
print('  %-5s %-4s %-16s %-16s %-14s' % ('L','N','||Pi r||^2/K','omitted tail','pred inflation'))
for c in cells:
    rt = c['ratio_nr2_K']
    print('  %-5s %-4d %-16.8f %-16.4e %-14.4f' % (c['L'], c['N'], rt, 1-rt, 1/rt))
print()
print('  effect on the H1 relative error: as-implemented (4K*proj2/nr2) vs')
print('  consistent finite-basis (4*proj2 = 4K*proj2/K)\n')
print('  %-5s %-9s %-6s %-12s %-12s %-10s %-12s %-10s %-10s' %
      ('L','delta','dimN','C meas','pred impl','rel impl','pred consist','rel cons','rel compr'))
tab, agg_i, agg_c, agg_k = [], [], [], []
for c in cells:
    rt = c['ratio_nr2_K']
    for s in c['sweep']:
        if s['dimN'] == 0: continue
        C = s['C_float']; pi = s['pred_float']; pc = pi*rt
        ri = abs(C-pi)/abs(C); rc = abs(C-pc)/abs(C)
        tab.append(dict(L=c['L'], delta=s['delta'], dimN=s['dimN'], C=C,
                        pred_impl=pi, rel_impl=ri, pred_consist=pc, rel_consist=rc,
                        rel_compr=s['rel_comp']))
        agg_i.append(ri); agg_c.append(rc); agg_k.append(s['rel_comp'])
        print('  %-5s %-9s %-6d %-12.6g %-12.6g %-10.3e %-12.6g %-10.3e %-10.3e' %
              (c['L'], s['delta'], s['dimN'], C, pi, ri, pc, rc, s['rel_comp']))
print()
print('  %-34s %-12s %-12s %-12s' % ('normalisation', 'median rel', 'mean rel', 'n<=3%'))
for nm, a in (('H1 as implemented (/nr2)', agg_i), ('H1 consistent (/K)', agg_c),
              ('compressed (retains F F\'\')', agg_k)):
    print('  %-34s %-12.4f %-12.4f %d of %d' %
          (nm, statistics.median(a), statistics.fmean(a), sum(1 for x in a if x <= 0.03), len(a)))
json.dump(tab, open('data/stage5_audit_norm.json','w'), indent=1)
