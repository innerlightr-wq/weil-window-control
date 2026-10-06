#!/usr/bin/env python3
import csv, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from detection_law import lam_min_closed

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R_ = lambda *p: os.path.join(ROOT, *p)
rows = list(csv.DictReader(open(R_('data', 'stage2_spectra.csv'))))
zrows = list(csv.DictReader(open(R_('data', 'stage2_zeros.csv'))))
mask = list(csv.DictReader(open(R_('data', 'stage2_masking.csv'))))
cfgs = ['C0-online', 'C1-g1real', 'C2-quartet', 'C3-twoquartets', 'C4-mixed']

# --- 1: lam_min(R) ---------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.3))
for c in cfgs:
    rs = [r for r in rows if r['config'] == c]
    Rv = [int(r['R']) for r in rs]
    lam = [float(r['lam_min']) for r in rs]
    style = 'k-o' if c == 'C0-online' else '-o'
    ax[0].plot(Rv, lam, style, ms=3.5, label=c, lw=2 if c == 'C0-online' else 1.2)
    fn = [r for r in rs if r['lam_min_sign'] == 'neg']
    if fn:
        R0 = int(fn[0]['R'])
        ax[0].plot([R0], [float(fn[0]['lam_min'])], 'v', ms=11, mfc='none',
                   color=ax[0].lines[-1].get_color())
ax[0].set_yscale('symlog', linthresh=1e-3)
ax[0].axhline(0, color='grey', lw=.7)
ax[0].set_xlabel('window R'); ax[0].set_ylabel(r'$\lambda_{\min}(T_R)$  (symlog)')
ax[0].set_title('Sign of $\\lambda_{\\min}$ separates on-line from planted off-line\n'
                '(triangles = first certified negative)')
ax[0].legend(fontsize=7.5); ax[0].grid(alpha=.3)

ax[1].plot([int(m['k_online_pairs']) for m in mask],
           [int(m['first_neg_R_EXACT']) for m in mask], 'o-', color='tab:red',
           label='detection window $R_{\\rm det}$')
ax[1].plot([int(m['k_online_pairs']) for m in mask],
           [int(m['total_rank']) for m in mask], 's--', color='tab:grey',
           label='total rank $=4+2k$')
ax[1].set_xlabel('k = number of masking ON-LINE pairs added')
ax[1].set_ylabel('R')
ax[1].set_title('Masking is weak: one off-line quartet stays detectable\n'
                'in a window of size 13 among 26 on-line zeros')
ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
plt.tight_layout(); plt.savefig(R_('figures', 'stage2_lambda_and_masking.png'), dpi=140)
plt.close()

# --- 2: the money figure: zeros stay on the circle ------------------------
fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.6), subplot_kw={'aspect': 'equal'})
th = np.linspace(0, 2 * np.pi, 500)
for a, (cfg, planted) in zip(axs, [('C2-quartet', [(1.3, 0.7)]),
                                   ('C4-mixed', [(1.3, 0.7), (1.0, 2.0)]),
                                   ('C3-twoquartets', [(1.3, 0.7), (1.3, 2.0)])]):
    a.plot(np.cos(th), np.sin(th), 'k-', lw=.8)
    for R, col in [(6, 'tab:blue'), (12, 'tab:green'), (20, 'tab:red')]:
        zs = [r for r in zrows if r['config'] == cfg and int(r['R']) == R]
        ang = np.array([float(r['zero_angle']) for r in zs])
        rad = np.array([float(r['abs_zero']) for r in zs])
        a.plot(rad * np.cos(ang), rad * np.sin(ang), 'o', ms=5.5, mfc='none',
               color=col, label=f'ground-state zeros, R={R}')
    for rho, phi in planted:
        for s in (1, -1):
            a.plot([0, 1.62 * np.cos(s * phi)], [0, 1.62 * np.sin(s * phi)],
                   '--', lw=.9, color='goldenrod', zorder=0)
        for s in (1, -1):
            for rr in ({rho, 1 / rho} if rho != 1 else {1.0}):
                a.plot([rr * np.cos(s * phi)], [rr * np.sin(s * phi)], 'k*', ms=14,
                       mfc='gold', mec='k')
    a.set_title(f'{cfg}\n(stars = planted $\\beta_j$, off the circle)', fontsize=9.5)
    a.legend(fontsize=7, loc='lower left'); a.grid(alpha=.25)
    a.set_xlim(-1.7, 1.7); a.set_ylim(-1.7, 1.7)
fig.suptitle('The witness is free: ground-state zeros stay EXACTLY on $|z|=1$ and lock onto the '
             'planted ANGLES,\nwhile $\\lambda_{\\min}$ is hugely negative. Radial information is '
             'structurally invisible to the zeros.', fontsize=11)
plt.tight_layout(); plt.savefig(R_('figures', 'stage2_zeros_stay_on_circle.png'), dpi=140)
plt.close()

# --- 3: detection law -----------------------------------------------------
mp.mp.dps = 60
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.3))
eps = [mp.mpf(10) ** (-k) for k in range(1, 11)]
for Rv, col in [(1, 'tab:blue'), (4, 'tab:orange'), (16, 'tab:green'), (64, 'tab:red')]:
    y = [float(abs(lam_min_closed(1 + e, Rv))) for e in eps]
    ax[0].loglog([float(e) for e in eps], y, 'o-', color=col, label=f'R={Rv}')
    C = Rv * (Rv + 1) * (Rv + 2) / 6
    ax[0].loglog([float(e) for e in eps], [C * float(e) ** 2 for e in eps], ':', color=col)
for f, lab in [(1e-16, 'float64 noise $10^{-16}$'), (1.66e-17, "Zhu/Connes $\\lambda^*\\!\\approx\\!1.7{\\times}10^{-17}$")]:
    ax[0].axhline(f, ls='--', lw=1, color='k' if f > 1e-16 else 'grey')
    ax[0].text(1.2e-10, f * 1.4, lab, fontsize=7.5)
ax[0].set_xlabel(r'off-line distance $\varepsilon$   ($\rho = 1+\varepsilon$)')
ax[0].set_ylabel(r'$|\lambda_{\min}(T_R)|$')
ax[0].set_title(r'Detection law (exact):  $\lambda_{\min} = -\binom{R+2}{3}\varepsilon^2 + O(\varepsilon^3)$'
                '\ndotted = the cubic-binomial law')
ax[0].legend(fontsize=8); ax[0].grid(alpha=.3, which='both')

dl = list(csv.DictReader(open(R_('data', 'stage2_detection_law.csv'))))
for key, lab, col in [('R_for_floor_1e-16', 'noise floor $10^{-16}$ (float64)', 'tab:red'),
                      ('R_for_floor_1e-17', 'noise floor $10^{-17}$', 'tab:orange'),
                      ('R_for_floor_1e-30', 'noise floor $10^{-30}$', 'tab:green'),
                      ('R_for_floor_1e-50', 'noise floor $10^{-50}$', 'tab:blue')]:
    e = [float(d['eps']) for d in dl if d[key]]
    r = [int(d[key]) for d in dl if d[key]]
    ax[1].loglog(e, r, 'o-', color=col, label=lab)
ax[1].set_xlabel(r'off-line distance $\varepsilon$')
ax[1].set_ylabel('smallest window R that certifies a violation')
ax[1].set_title('Window size needed is governed by PRECISION, not by $\\varepsilon$:\n'
                'at infinite precision R=1 always suffices')
ax[1].legend(fontsize=8); ax[1].grid(alpha=.3, which='both')
plt.tight_layout(); plt.savefig(R_('figures', 'stage2_detection_law.png'), dpi=140)
print('stage2 figures written')
