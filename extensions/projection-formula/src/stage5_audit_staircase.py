#!/usr/bin/env python3
"""Audit B: what the "11 of 12 staircase steps" claim actually asserts, and how
informative the alignment test is.

EVENT DEFINITION (as implemented in stage2b_k2k3.py):
  * delta is swept on a DECADE grid 10^-1 ... 10^-12 (12 points, 11 gaps).
  * A "C-step" is a consecutive-decade pair (d_{i-1}, d_i) with
        |C(d_{i-1}) - C(d_i)| / max(|C(d_{i-1})|, 1e-300) > 0.20.
    The 12 is therefore the OBSERVED number of such pairs summed over the 5 windows
    (4+5+2+1+0), not a preregistered denominator.
  * A step is "aligned" if some delta_k = sqrt(lam_k/4K) satisfies
        min(lo,hi)/3 <= delta_k <= max(lo,hi)*3.
    With hi = lo/10 this window is [lo/30, 3*lo]: a factor 90 in delta,
    i.e. a factor 8100 in delta^2.
  * PREREGISTRATION (K3) said "within a factor 3 in delta^2", i.e. the window
    [lo_2/3, hi_2*3] in delta^2 where (lo_2,hi_2) = (lo^2, hi^2) -- a factor 300 in
    delta^2, 27x tighter than what was implemented. Both are reported here.

NULL MODEL: the alignment test is only informative if a randomly placed step would
usually FAIL it. We compute, per window, the fraction of the 11 available decade
brackets that would be scored "aligned" -- the pass rate under the null that step
locations are independent of the spectrum.
"""
import json, itertools
import mpmath as mp
mp.mp.dps = 40

sweep_cells = json.load(open('data/stage2_zeta.json'))
full = {c['L']: c for c in json.load(open('data/stage5_eigs_full.json'))}
DECADES = [mp.mpf(10)**(-e) for e in range(1, 13)]
BRACKETS = [(DECADES[i-1], DECADES[i]) for i in range(1, 12)]   # 11 available brackets

def aligned_impl(lo, hi, dk):
    a, b = min(lo, hi)/3, max(lo, hi)*3
    return any(a <= x <= b for x in dk)

def aligned_prereg(lo, hi, dk):
    """factor 3 in delta^2, as preregistered"""
    a, b = min(lo, hi)**2/3, max(lo, hi)**2*3
    return any(a <= x**2 <= b for x in dk)

print('  %-5s %-4s %-7s %-9s %-11s %-11s %-13s %-13s' %
      ('L','N','#dk','#C-steps','aligned(impl)','aligned(prereg)','null impl','null prereg'))
tot = dict(steps=0, ai=0, ap=0); nullrows = []
for c in sweep_cells:
    Lx = c['L']; f = full[Lx]
    fourK = mp.mpf(f['fourK'])
    ev = [mp.mpf(x) for x in f['eigs_all']]
    dk = [mp.sqrt(e/fourK) for e in ev if e > 0]
    sw = c['sweep']
    steps = []
    for i in range(1, len(sw)):
        a_, b_ = sw[i-1]['C_float'], sw[i]['C_float']
        if b_ > 0 and abs(a_-b_)/max(abs(a_), 1e-300) > 0.2:
            steps.append((mp.mpf(sw[i-1]['delta']), mp.mpf(sw[i]['delta'])))
    ai = sum(1 for lo, hi in steps if aligned_impl(lo, hi, dk))
    ap = sum(1 for lo, hi in steps if aligned_prereg(lo, hi, dk))
    ni = sum(1 for lo, hi in BRACKETS if aligned_impl(lo, hi, dk))
    npg = sum(1 for lo, hi in BRACKETS if aligned_prereg(lo, hi, dk))
    print('  %-5s %-4d %-7d %-9d %-11s %-11s %-13s %-13s' %
          (Lx, c['N'], len(dk), len(steps), '%d/%d' % (ai, len(steps)),
           '%d/%d' % (ap, len(steps)), '%d/11 = %.0f%%' % (ni, 100*ni/11),
           '%d/11 = %.0f%%' % (npg, 100*npg/11)))
    tot['steps'] += len(steps); tot['ai'] += ai; tot['ap'] += ap
    nullrows.append(dict(L=Lx, n_dk=len(dk), steps=len(steps), aligned_impl=ai,
                         aligned_prereg=ap, null_impl=ni, null_prereg=npg,
                         null_impl_rate=ni/11, null_prereg_rate=npg/11))
NI = sum(r['null_impl'] for r in nullrows); NP = sum(r['null_prereg'] for r in nullrows)
print()
print('  TOTAL: %d C-steps; aligned %d (implemented criterion), %d (preregistered criterion)'
      % (tot['steps'], tot['ai'], tot['ap']))
print('  NULL pass rate over all 55 decade brackets: %d/55 = %.0f%% (implemented), %d/55 = %.0f%% (prereg)'
      % (NI, 100*NI/55, NP, 100*NP/55))
print()
print('  Interpretation: the implemented criterion is passed by %.0f%% of ALL brackets, so')
print('  "%d of %d aligned" carries little information beyond the null. The preregistered' % (tot['ai'], tot['steps']))
print('  factor-3-in-delta^2 criterion has null rate %.0f%% and is the one to report.' % (100*NP/55))
json.dump(dict(per_window=nullrows, total=tot, null_impl_55=NI, null_prereg_55=NP),
          open('data/stage5_audit_staircase.json','w'), indent=1)
