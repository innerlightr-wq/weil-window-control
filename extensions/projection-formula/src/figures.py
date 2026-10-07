#!/usr/bin/env python3
"""Figures for the projection-formula report."""
import json, math
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

cells = json.load(open('data/stage2_zeta.json'))
k2 = json.load(open('data/stage2b_k2.json'))
par = json.load(open('data/stage3_parity.json'))
flo = json.load(open('data/stage4_floor.json'))

# --- Fig 1: C(delta) staircase with eigenvalue markers -------------------------
fig, axes = plt.subplots(1, 5, figsize=(19, 3.6), sharey=False)
for ax, c in zip(axes, cells):
    d = [float(s['delta']) for s in c['sweep']]
    C = [s['C_float'] for s in c['sweep']]
    P = [s['pred_float'] for s in c['sweep']]
    fourK = float(c['fourK'])
    ax.loglog(d, C, 'o-', color='tab:red', ms=4, label='measured $C(\\delta)$')
    ax.loglog(d, [max(x, 1e-300) for x in P], 's--', color='tab:blue', ms=3,
              label='$4K\\cos^2\\theta_N$')
    ax.axhline(fourK, color='k', lw=.8, ls=':', label='$4K$')
    for e in c['eigs']:
        ev = float(e)
        if ev > 0:
            dk = (ev/fourK)**0.5
            if min(d) <= dk <= max(d): ax.axvline(dk, color='grey', lw=.7, alpha=.65)
    ax.set_title('$L=%s$' % c['L'], fontsize=10)
    ax.set_xlabel('$\\delta$'); ax.grid(alpha=.3, which='both')
axes[0].set_ylabel('$C(\\delta)$'); axes[0].legend(fontsize=7)
fig.suptitle('Staircase: $C(\\delta)$ against the projection prediction; grey lines are '
             'eigenvalue scales $\\delta_k=\\sqrt{\\lambda_k/4K}$ of $Q_0$', fontsize=11)
plt.tight_layout(); plt.savefig('figures/staircase.png', dpi=140); plt.close()

# --- Fig 2: cos^2 theta vs C/4K at delta = 0.02 --------------------------------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
Ls = [r['L'] for r in k2]; cos2 = [r['cos2'] for r in k2]
c4k = [r['C']/float(next(c['fourK'] for c in cells if c['L'] == r['L'])) for r in k2]
ax[0].plot([0, 1], [0, 1], 'k--', lw=.8, label='$y=x$')
ax[0].plot(cos2, c4k, 'o', ms=9, color='tab:red')
for x, y, l in zip(cos2, c4k, Ls): ax[0].annotate('$L=%s$' % l, (x, y), fontsize=8,
                                                  xytext=(5, -9), textcoords='offset points')
ax[0].set_xlabel('$\\cos^2\\theta_N$  (projection of $r_\\gamma$ on the near-null space)')
ax[0].set_ylabel('$C/4K$  (paper \\S4.3 saturation ratio)')
ax[0].set_title('H1 at $\\delta=0.02$: the saturation ratio\nis the projection fraction', fontsize=10)
ax[0].grid(alpha=.3); ax[0].legend(fontsize=8)
rel = [r['rel'] for r in k2]; relc = [r['relc'] for r in k2]
ax[1].semilogy(range(5), rel, 'o-', color='tab:red', label='scalar $4K\\cos^2\\theta_N$')
ax[1].semilogy(range(5), relc, 's-', color='tab:green', label='compressed $\\Lambda_N+\\delta^2\\Pi P\\Pi$')
ax[1].axhline(0.10, color='k', ls=':', lw=1, label='K2 threshold 10%')
ax[1].set_xticks(range(5)); ax[1].set_xticklabels(['$L=%s$' % l for l in Ls])
ax[1].set_ylabel('relative error'); ax[1].set_title('Accuracy at $\\delta=0.02$', fontsize=10)
ax[1].grid(alpha=.3, which='both'); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig('figures/cos2_vs_ratio.png', dpi=140); plt.close()

# --- Fig 3: parity -------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.4, 4.4))
L = [float(r['L']) for r in par]
ax.semilogy(L, [float(r['even']) for r in par], 'o-', label='even floor')
ax.semilogy(L, [float(r['constrained']) for r in par], '^-', label='even with $\\int f=0$')
ax.semilogy(L, [float(r['odd']) for r in par], 's-', label='odd floor')
ax.set_xlabel('$L$'); ax.set_ylabel('$\\lambda_{\\min}$')
ax.set_title('H3 FAILED: constraining $F(0)=0$ raises the even floor\nfar above the odd floor '
             '(53$\\times$ to 893$\\times$)', fontsize=10)
ax.grid(alpha=.3, which='both'); ax.legend(fontsize=9)
plt.tight_layout(); plt.savefig('figures/parity.png', dpi=140); plt.close()

# --- Fig 4: floor vs Landau-Widom ---------------------------------------------
fig, ax = plt.subplots(figsize=(6.4, 4.4))
Lf = [float(r['L']) for r in flo]
ax.plot(Lf, [r['mln'] for r in flo], 'o-', color='tab:red', label='$-\\ln\\lambda^*(L)$ measured')
ax.plot(Lf, [r['zhu'] for r in flo], 's--', color='tab:blue',
        label="Zhu $2\\pi^2 N(T^*)/\\ln N(T^*)$")
ax.set_xlabel('$L$'); ax.set_ylabel('$-\\ln\\lambda^*$')
ax.set_title('H4 (confirmation only, prior art): the floor tracks\nthe Landau–Widom law to a '
             'slowly varying factor', fontsize=10)
ax.grid(alpha=.3); ax.legend(fontsize=9)
plt.tight_layout(); plt.savefig('figures/floor.png', dpi=140); plt.close()
print('wrote 4 figures')
