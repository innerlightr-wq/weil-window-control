#!/usr/bin/env python3
"""Headline figure: Q1 (t0/T1), P3 (who carries the deep dips), P4 ({2,3} share)."""
import sys, os, json, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
t0 = list(csv.DictReader(open(os.path.join(ROOT, 'data', 'stage2_t0.csv'))))
s4 = json.load(open(os.path.join(ROOT, 'data', 'stage4_minima.json')))
p4 = json.load(open(os.path.join(ROOT, 'data', 'stage5_p4.json')))

fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.5))

# (a) Q1
for b, mk, c in [('0.0', 'o-', 'tab:red'), ('0.25', 's-', 'tab:orange'), ('0.5', '^-', 'tab:blue')]:
    r = [x for x in t0 if x['beta'] == b]
    ax[0].plot([float(x['L']) for x in r], [float(x['t0_over_T1']) for x in r],
               mk, color=c, ms=5, label=r'$\beta$ = ' + b)
ax[0].axhline(1.0, color='k', lw=1, label=r"Zhu's threshold $T_1 = 2\pi e^{A_L}$")
ax[0].axhline(0.1, color='k', ls='--', lw=1, label='K2 kill line ($T_1/10$)')
ax[0].set_ylim(0, 1.6); ax[0].set_xlabel('$L$'); ax[0].set_ylabel('$t_0(L,\\beta)\\,/\\,T_1(L)$')
ax[0].set_title('Q1: the true last dip against the worst-case\nthreshold — $T_1$ is essentially sharp',
                fontsize=10)
ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)

# (b) P3 at L=1.6
r = [x for x in s4 if x['real']['L'] == 1.6][0]
ps = sorted(r['real']['share_deep'], key=int)
x = np.arange(len(ps))
lo = [r['surrogate']['share_deep'][p][0] for p in ps]
hi = [r['surrogate']['share_deep'][p][2] for p in ps]
ax[1].fill_between(x, lo, hi, color='tab:blue', alpha=.22, label='N2 surrogate 95% band')
ax[1].plot(x, [r['real']['share_deep'][p] for p in ps], 'o-', color='tab:red', label='real primes')
ax[1].plot(x, [r['n1']['share_deep'][p] for p in ps], 'd--', color='tab:green', label='N1 torus model')
w = np.array([r['real']['weight'][p] for p in ps]); w = w / w.max()
ax[1].bar(x, w * 0.25, bottom=0.30, color='grey', alpha=.3, width=.5,
          label='block weight $A_p$ (scaled)')
ax[1].set_xticks(x); ax[1].set_xticklabels(ps)
ax[1].set_xlabel('prime $p$'); ax[1].set_ylabel('alignment share of the deepest 1%')
ax[1].set_title('P3 FAILS as stated ($L$=1.6): the deepest dips are\n'
                'not carried by the heaviest primes ($\\rho=-0.62$)', fontsize=10)
ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)

# (c) P4
Ls = [d['L'] for d in p4]
ax[2].fill_between(Ls, [d['n2_band'][0] for d in p4], [d['n2_band'][2] for d in p4],
                   color='tab:blue', alpha=.22, label='N2 surrogate 95% band')
ax[2].plot(Ls, [d['real'] for d in p4], 'o-', color='tab:red', label='real primes')
ax[2].plot(Ls, [d['n1_pred'] for d in p4], 'd--', color='tab:green', label='N1 torus prediction')
ax[2].set_xlabel('$L$'); ax[2].set_ylabel('dip measure inside $S_{23}(1/2)$')
ax[2].set_title('P4 FAILS as stated: the {2,3} share does not decay —\n'
                'but the surrogates do not decay either', fontsize=10)
ax[2].legend(fontsize=8); ax[2].grid(alpha=.3); ax[2].set_ylim(0.4, 1.0)

plt.tight_layout()
plt.savefig(os.path.join(ROOT, 'figures', 'summary.png'), dpi=140)
print('wrote figures/summary.png')
