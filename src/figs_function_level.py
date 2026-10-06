#!/usr/bin/env python3
"""The (iv) re-scoring figure: angles converge (free), functions do not (RH content)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from planted import t_from_blocks, betas
from function_level import target_poly, normalized_value_at
from weilform import toeplitz, eigsym, poly_zeros

mp.mp.dps = 50
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
Rmax = 30
CFG = [('C1-g1real', [('realpair', '1.3')], 'tab:blue'),
       ('C2-quartet', [('quartet', '1.3', 0.7)], 'tab:red'),
       ('C4-mixed', [('quartet', '1.3', 0.7), ('online', 2.0)], 'tab:green')]

fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
for tag, blocks, col in CFG:
    t, _ = t_from_blocks(blocks, Rmax)
    Ps, uniq = target_poly(betas(blocks))
    off = [b for b in uniq if abs(abs(b) - 1) > mp.mpf(10) ** -20]
    angerr, fnerr, Rs = [], [], []
    for R in range(2, Rmax + 1):
        v = eigsym(toeplitz([mp.mpf(x) for x in t], R))[1][0]
        zs = poly_zeros(v)
        if not zs:
            continue
        za = [mp.arg(z) for z in zs]
        e = max(min(abs(mp.arg(b) - a) for a in za) for b in off)
        angerr.append(max(float(e), 1e-8))
        fnerr.append(float(min(normalized_value_at(v, b) for b in off)))
        Rs.append(R)
    ax[0].semilogy(Rs, angerr, 'o-', ms=3.5, color=col, label=tag)
    ax[1].plot(Rs, fnerr, 'o-', ms=3.5, color=col, label=tag)

ax[0].set_xlabel('window R'); ax[0].set_ylabel(r'$\max_\beta \min_z |\arg\beta - \arg z|$')
ax[0].set_title('(iv) read as ZERO ANGLES: converges\n'
                '— free, happens for off-line spectra too', fontsize=10.5)
ax[0].legend(fontsize=8); ax[0].grid(alpha=.3, which='both')

ax[1].axhline(0, color='k', lw=.8)
ax[1].set_ylim(-0.03, 0.85)
ax[1].set_xlabel('window R')
ax[1].set_ylabel(r'$|\hat c_R(\beta)| \,/\, \|c_R\|_2 (\sum_k |\beta|^{2k})^{1/2}$')
ax[1].set_title('(iv) read as FUNCTION convergence (Hurwitz): FAILS\n'
                '— the minimizer never learns to vanish at the off-circle $\\beta$',
                fontsize=10.5)
ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
ax[1].annotate('convergence would require this $\\to 0$', xy=(16, 0.03), fontsize=8.5,
               xytext=(9, 0.18), arrowprops=dict(arrowstyle='->', lw=.9))
fig.suptitle('Re-scoring step (iv): which reading you take decides whether it carries RH content',
             fontsize=11.5)
plt.tight_layout()
plt.savefig(os.path.join(ROOT, 'figures', 'stage2_iv_rescored.png'), dpi=140)
print('wrote figures/stage2_iv_rescored.png')
