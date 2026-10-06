#!/usr/bin/env python3
import csv, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rows = list(csv.DictReader(open(os.path.join(ROOT, 'results', 'stage1_spectra.csv'))))
zrows = list(csv.DictReader(open(os.path.join(ROOT, 'results', 'stage1_zeros.csv'))))
curves = ['E5', 'H2F3', 'G3F5', 'G3F3']
G = {'E5': 1, 'H2F3': 2, 'G3F5': 3, 'G3F3': 3}

fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for c in curves:
    rs = [r for r in rows if r['curve'] == c]
    R = [int(r['R']) for r in rs]
    lam = [max(float(r['lam_min']), 1e-55) for r in rs]
    ax[0].semilogy(R, lam, 'o-', label=f"{c} (g={G[c]}, 2g={2*G[c]})")
    ax[0].axvline(2 * G[c], ls=':', lw=.8, color=ax[0].lines[-1].get_color())
    f64 = [abs(float(r['lam_min_float64'])) for r in rs]
    ax[1].semilogy(R, f64, 's--', label=f"{c} |float64|")
ax[0].set_xlabel('window R'); ax[0].set_ylabel(r'$\lambda_{\min}(T_R)$ (clipped at 1e-55)')
ax[0].set_title(r'Weil form ground state collapses to 0 exactly at $R=2g$' '\n(mpmath, 50 dps; dotted lines at R=2g)')
ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
ax[1].set_xlabel('window R'); ax[1].set_ylabel(r'$|\lambda_{\min}|$ in float64')
ax[1].set_title('float64 evaluation: spurious $O(10^{-16})$ values\nwhere the true value is exactly 0')
ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
plt.tight_layout(); plt.savefig(os.path.join(ROOT, 'figures', 'stage1_lambda_min.png'), dpi=140)
plt.close()

# zeros on the circle, R = 2g and R = 2g-1, 2g-2
fig, axs = plt.subplots(1, 4, figsize=(16, 4.2), subplot_kw={'aspect': 'equal'})
for k, c in enumerate(curves):
    a = axs[k]; g = G[c]
    th = np.linspace(0, 2 * np.pi, 400)
    a.plot(np.cos(th), np.sin(th), 'k-', lw=.6)
    for R, mk, col in [(2 * g - 2, 'v', 'tab:orange'), (2 * g - 1, '^', 'tab:green'),
                       (2 * g, 'o', 'tab:red')]:
        zs = [r for r in zrows if r['curve'] == c and int(r['R']) == R]
        if not zs: continue
        ang = np.array([float(r['zero_angle']) for r in zs])
        rad = np.array([float(r['abs_zero']) for r in zs])
        a.plot(rad * np.cos(ang), rad * np.sin(ang), mk, ms=7, mfc='none', color=col,
               label=f'R={R}')
    a.set_title(f'{c}  (g={g})', fontsize=10); a.legend(fontsize=7, loc='upper right')
    a.set_xlim(-1.35, 1.35); a.set_ylim(-1.35, 1.35); a.grid(alpha=.25)
fig.suptitle('Ground-state polynomial zeros (H1): all on |z|=1 whenever $\\lambda_{\\min}$ is simple;'
             ' at R=2g they are exactly $e^{i\\theta_j}$ (H2)', fontsize=11)
plt.tight_layout(); plt.savefig(os.path.join(ROOT, 'figures', 'stage1_zeros_circle.png'), dpi=140)
plt.close()

# H4 convergence
fig, ax = plt.subplots(figsize=(6.2, 4.2))
for c in curves:
    rs = [r for r in rows if r['curve'] == c and r['H4_theta_to_zero'] and int(r['R']) <= 2 * G[c]]
    R = [int(r['R']) for r in rs]
    d = [max(float(r['H4_theta_to_zero']), 1e-55) for r in rs]
    ax.semilogy(R, d, 'o-', label=f"{c} (g={G[c]})")
ax.set_xlabel('window R'); ax.set_ylabel(r'$\max_j \min_z |\theta_j - \arg z|$')
ax.set_title('H4: coverage of the true Frobenius angles by\nground-state zeros, as the window grows')
ax.legend(fontsize=8); ax.grid(alpha=.3)
plt.tight_layout(); plt.savefig(os.path.join(ROOT, 'figures', 'stage1_h4_convergence.png'), dpi=140)
print('figures written')
