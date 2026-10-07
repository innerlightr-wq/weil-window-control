#!/usr/bin/env python3
"""Dose-response figure with comma markers."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from symbol import Symbol
from commas import commas

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
fine = json.load(open(os.path.join(ROOT, 'data', 'stage3_dose_fine.json')))
coarse = json.load(open(os.path.join(ROOT, 'data', 'stage3_dose.json')))

NAMED = {'9/8': (9, 8), '25/24': (25, 24), '128/125': (128, 125), '81/80': (81, 80)}
s16 = Symbol(1.6)
cm = commas(sorted(set(int(p) for p in s16.prime)), 200, 0.30)
cm_by = {'%d/%d' % (a, b): (d, st) for a, b, d, s2, st in cm}

fig, ax = plt.subplots(1, 2, figsize=(12.4, 4.8))

# --- left: the dose-response ---
for cell, mk, col in [((1.6, 0.0), 'o-', 'tab:red'), ((1.6, 0.25), 's-', 'tab:orange'),
                      ((1.0, 0.25), '^-', 'tab:blue')]:
    src = fine if cell == (1.6, 0.0) else coarse
    rows = sorted([r for r in src if (r['L'], r['beta']) == cell], key=lambda r: r['spread'])
    ax[0].semilogx([r['spread'] for r in rows],
                   [100 * r['excess_low_t'] / r['sur_low_t_median'] for r in rows],
                   mk, color=col, ms=5, label='L=%s, $\\beta$=%s (low-$t$)' % cell)
ax[0].axhline(0, color='k', lw=.8)
for nm, (a, b) in NAMED.items():
    d, st = cm_by.get(nm, (abs(np.log(a / b)), None))
    ax[0].axvline(d, color='grey', ls=':', lw=1)
    ax[0].text(d, 13.2, nm, rotation=90, fontsize=7.5, ha='right', va='top', color='grey')
# effective preservation scales for the same commas
for nm, (a, b) in NAMED.items():
    st = cm_by.get(nm, (None, None))[1]
    if st:
        ax[0].axvline(st, color='tab:green', ls='--', lw=.9, alpha=.7)
ax[0].axvspan(1e-4, 1.3e-4, color='tab:purple', alpha=.18)
ax[0].text(1.1e-4, 6.5, 'trivial\nconvergence\n($\\sigma t\\lesssim1$)', fontsize=7.5,
           color='tab:purple', ha='center')
ax[0].set_xlabel('surrogate spread $\\sigma$ (log $p$ perturbed by $U(-\\sigma,\\sigma)$)')
ax[0].set_ylabel('low-$t$ dip-measure excess over surrogate median (%)')
ax[0].set_title('Dose-response: the excess fades below $\\sigma\\approx10^{-2}$\n'
                'dotted grey = comma detunings; dashed green = preservation scales $s^*$',
                fontsize=10)
ax[0].legend(fontsize=8); ax[0].grid(alpha=.3, which='both')

# --- right: the confound ---
rows = sorted(fine, key=lambda r: r['spread'])
ax[1].semilogx([r['spread'] for r in rows],
               [100 * r['excess_low_t'] / r['sur_low_t_median'] for r in rows],
               'o-', color='tab:red', label='low-$t$ excess (%)')
ax2 = ax[1].twinx()
ax2.semilogx([r['spread'] for r in rows], [r['overlap_mean'] for r in rows],
             's--', color='tab:purple', label='overlap with real dips')
ax2.set_ylabel('fraction of real dip measure covered by surrogate', color='tab:purple')
ax[1].set_xlabel('surrogate spread $\\sigma$')
ax[1].set_ylabel('low-$t$ excess (%)', color='tab:red')
ax[1].set_title('Why the test is confounded at L=1.6:\nthe excess fades exactly as the '
                'surrogate converges to the real symbol', fontsize=10)
ax[1].grid(alpha=.3, which='both')
h1, l1 = ax[1].get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax[1].legend(h1 + h2, l1 + l2, fontsize=8, loc='upper left')
plt.tight_layout()
plt.savefig(os.path.join(ROOT, 'figures', 'dose_response.png'), dpi=140)
print('wrote figures/dose_response.png')
