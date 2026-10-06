#!/usr/bin/env python3
import csv, os, sys
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R_ = lambda *p: os.path.join(ROOT, *p)
rd = lambda n: list(csv.DictReader(open(R_('results', n))))

# ---- 1: Gate A + Gate B -------------------------------------------------
A = rd('stage3_gateA_explicit_formula.csv')
B = rd('stage3_gateB_lambda_min.csv')
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
for sec, col, lab in [('even_d', 'tab:blue', 'even (vanishes at $\\pm L$)'),
                      ('odd', 'tab:green', 'odd'),
                      ('even', 'tab:red', 'even (full cosine, jumps at $\\pm L$)')]:
    rs = [r for r in A if r['sector'] == sec]
    if not rs: continue
    g = [float(r['gamma_K']) for r in rs]
    d = [float(r['max_abs_diff']) for r in rs]
    ax[0].loglog(g, d, 'o-', color=col, label=lab)
gg = np.array([140, 1200])
ax[0].loglog(gg, 3e-2 * (gg / 140.) ** -3, 'k:', lw=1, label=r'$\gamma_K^{-3}$ (continuous f)')
ax[0].loglog(gg, 2.2e-2 * (gg / 140.) ** -1, 'k--', lw=1, label=r'$\gamma_K^{-1}$ (f jumps)')
ax[0].set_xlabel(r'$\gamma_K$ (highest zeta zero used)')
ax[0].set_ylabel(r'$\max_{jk}|$geometric $-$ zeros$|$')
ax[0].set_title('GATE A: assembled geometric side vs. $\\sum_\\rho |F(\\gamma_\\rho)|^2$\n'
                'discrepancy = the zeros-sum truncation tail, as predicted', fontsize=10)
ax[0].legend(fontsize=7.5); ax[0].grid(alpha=.3, which='both')

ax[1].axhspan(8.9e-18, 2.27e-17, color='tab:green', alpha=.18,
              label="Zhu certified $\\lambda^*(0.8)$")
for sec, col, lab in [('even_d', 'tab:blue', 'basis A: $\\cos\\frac{(2k+1)\\pi x}{2L}$'),
                      ('even', 'tab:red', 'basis B: $\\cos\\frac{k\\pi x}{L}$')]:
    rs = sorted([r for r in B if r['sector'] == sec], key=lambda r: int(r['N']))
    ax[1].semilogy([int(r['N']) for r in rs], [float(r['lam_min']) for r in rs],
                   'o-', color=col, label=lab)
ax[1].set_xlabel('basis size N'); ax[1].set_ylabel(r'$\lambda_{\min}$ at $L=0.8$')
ax[1].set_title('GATE B: two independent complete bases descend\nmonotonically onto '
                "Zhu's interval ($\\lambda_{\\min}(N)$ is an upper bound)", fontsize=10)
ax[1].legend(fontsize=8); ax[1].grid(alpha=.3, which='both')
plt.tight_layout(); plt.savefig(R_('figures', 'stage3_gates.png'), dpi=140); plt.close()

# ---- 2: Connes error profile -------------------------------------------
C = rd('stage3b_connes_zeros.csv')
fig, ax = plt.subplots(figsize=(7.2, 4.6))
n = [int(r['n']) for r in C]
e = [max(float(r['abs_err']), 1e-45) for r in C]
ax.semilogy(n, e, 'o-', ms=4, color='tab:purple')
ax.axvline(12.5, ls='--', color='grey', lw=1)
ax.text(13.2, 1e-30, 'beyond $n\\approx12$ the ground state\nno longer tracks the zeros\n'
                     '(basis exhausted at $N=24$)', fontsize=8)
ax.set_xlabel('zeta zero index $n$')
ax.set_ylabel(r'$|\gamma_n - (\text{nearest zero of } \hat f_{\rm gs})|$')
ax.set_title("Connes's experiment reproduced: support $[1,13]$, $L=\\frac{\\log 13}{2}$\n"
             "$N=24$, 110 digits, $\\lambda_{\\min}=3.65\\times10^{-43}$ (simple, even)", fontsize=10.5)
ax.grid(alpha=.3, which='both')
plt.tight_layout(); plt.savefig(R_('figures', 'stage3_connes_profile.png'), dpi=140); plt.close()

# ---- 3: zeta-side teeth -------------------------------------------------
T = rd('stage3e_teeth.csv')
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
Ls = [float(r['L']) for r in T]
ax[0].semilogy(Ls, [float(r['lam_min_control']) for r in T], 'k-o', ms=4,
               label='control: 400 true zeros only (PSD by construction)')
cols = {'0.2': 'tab:red', '0.1': 'tab:orange', '0.05': 'tab:green', '0.02': 'tab:blue'}
for d, col in cols.items():
    y = [float(r['lam_min_delta_' + d]) for r in T]
    neg = [(L, -v) for L, v in zip(Ls, y) if v < 0]
    if neg:
        ax[0].semilogy([p[0] for p in neg], [p[1] for p in neg], 'o--', color=col, ms=4,
                       label=f'$-\\lambda_{{\\min}}$, planted $\\delta={d}$ '
                             f'(Re$\\,\\rho={0.5+float(d)}$)')
ax[0].set_xlabel('window half-width L'); ax[0].set_ylabel(r'$|\lambda_{\min}|$')
ax[0].set_title('Planted off-line quartet at height $\\gamma=14.13$:\n'
                'the form turns indefinite at $L\\approx0.5$–$0.6$ for every $\\delta$', fontsize=10)
ax[0].legend(fontsize=7); ax[0].grid(alpha=.3, which='both')

ds = np.array([0.2, 0.1, 0.05, 0.02])
for Lsel, col in [(0.8, 'tab:blue'), (1.3, 'tab:green'), (2.0, 'tab:red')]:
    r = min(T, key=lambda r: abs(float(r['L']) - Lsel))
    y = [-float(r['lam_min_delta_' + f'{d:g}']) for d in ds]
    ax[1].loglog(ds, y, 'o-', color=col, label=f'L={r["L"]}')
    ax[1].loglog(ds, y[0] * (ds / ds[0]) ** 2, ':', color=col)
ax[1].set_xlabel(r'off-line depth $\delta = \mathrm{Re}\,\rho - 1/2$')
ax[1].set_ylabel(r'$-\lambda_{\min}$')
ax[1].set_title(r'$\zeta$-side detection law: $\lambda_{\min}\approx -C(L)\,\delta^2$'
                '\n(dotted = exact $\\delta^2$; fitted exponents 2.03–2.07)', fontsize=10)
ax[1].legend(fontsize=8); ax[1].grid(alpha=.3, which='both')
plt.tight_layout(); plt.savefig(R_('figures', 'stage3_teeth.png'), dpi=140); plt.close()

# ---- 4: the height sweep (the main limitation) --------------------------
H = rd('stage3e_height_sweep.csv')
hs = [('3.0', 3.0), ('14.134725141734693', 14.13), ('18.0', 18.0),
      ('50.0', 50.0), ('100.0', 100.0), ('200.0', 200.0)]
fig, ax = plt.subplots(figsize=(7.4, 4.6))
Lh = [float(r['L']) for r in H]
ax.semilogy(Lh, [abs(float(r['lam_min_control'])) for r in H], 'k-o', ms=4,
            label='control (no planted zero)')
for key, lab in hs:
    y, xs = [], []
    for r in H:
        v = float(r['lam_min_h' + key])
        if r['sign_h' + key] == 'neg':
            xs.append(float(r['L'])); y.append(-v)
    if xs:
        ax.semilogy(xs, y, 'o--', ms=4, label=f'$\\gamma_*={lab:g}$ (detected)')
    else:
        ax.plot([], [], 'o--', label=f'$\\gamma_*={lab:g}$ (never detected)')
ax.set_xlabel('window half-width L'); ax.set_ylabel(r'$-\lambda_{\min}$ where negative')
ax.set_title('The main limitation: a short window only sees LOW-LYING zeros\n'
             'planted depth $\\delta=0.1$; $\\gamma_*=100,200$ undetected out to $L=3$',
             fontsize=10.5)
ax.legend(fontsize=7.5); ax.grid(alpha=.3, which='both')
plt.tight_layout(); plt.savefig(R_('figures', 'stage3_teeth_height.png'), dpi=140)
print('stage3 figures written')
