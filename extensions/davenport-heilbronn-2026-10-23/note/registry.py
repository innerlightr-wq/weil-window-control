r"""Derive every quantity the note prints, from the extension's data files only.

Nothing here reads the note.  `verify_note_numbers.py` imports `build()` and checks the
note against it.  Each entry is (value, sig_digits_supported, source), where
sig_digits_supported is the number of significant digits the SOURCE actually carries --
the note may print that many or fewer, never more.
"""
import csv
import json
import os
import re
from decimal import Decimal, getcontext

getcontext().prec = 200

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, 'data')
RES = os.path.join(HERE, 'results')
CERT = os.path.join(HERE, 'certificates')


def _sig(s):
    """significant digits in a decimal string"""
    s = s.strip().lstrip('+-')
    s = re.sub(r'[eE].*$', '', s)
    s = s.replace('.', '').lstrip('0')
    return len(s) or 1


def trunc(s, nsig):
    """Keep nsig significant digits of a decimal string, PRESERVING the exponent.

    A plain s[:n] silently drops the 'e-21' off '-2.6618...e-21' and turns a
    1e-21 eigenvalue into a number of order 1.  That bug is the reason this exists."""
    s = s.strip()
    m = re.match(r'^([+-]?)([0-9]*\.?[0-9]+)([eE][+-]?[0-9]+)?$', s)
    if not m:
        raise ValueError('not a decimal: %r' % s)
    sign, mant, exp = m.group(1), m.group(2), m.group(3) or ''
    digits = mant.replace('.', '')
    lead = len(digits) - len(digits.lstrip('0'))
    keep = lead + nsig
    # count characters of `mant` needed to retain `keep` digit positions
    seen, out = 0, ''
    for ch in mant:
        out += ch
        if ch.isdigit():
            seen += 1
            if seen >= keep:
                break
    return sign + out + exp


class Registry(dict):
    def add(self, name, value, sig, source):
        if name in self:
            raise KeyError('duplicate registry key %s' % name)
        self[name] = dict(value=Decimal(str(value)), sig=sig, source=source)

    def add_str(self, name, s, source, sig=None):
        self.add(name, Decimal(s), sig if sig is not None else _sig(s), source)

    def add_int(self, name, n, source):
        self.add(name, Decimal(int(n)), 999, source)


def read_csv(path):
    with open(path) as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith('#')))


def read_json(path):
    with open(path) as fh:
        return json.load(fh)


_BALL = re.compile(r'\[\s*([+-]?[0-9.]+(?:[eE][+-]?[0-9]+)?)\s*\+/-\s*'
                   r'([0-9.]+(?:[eE][+-]?[0-9]+)?)\s*\]')


def ball_fields(path):
    """{label: (midpoint_string, radius_string)} from a certify.py certificate."""
    out = {}
    for line in open(path):
        if '=' not in line:
            continue
        lab, rhs = line.split('=', 1)
        m = _BALL.search(rhs)
        if m:
            out[lab.strip()] = (m.group(1), m.group(2))
        else:
            out.setdefault(lab.strip() + '#plain', rhs.strip())
    return out


def build():
    R = Registry()

    # ------------------------------------------------------------------ Phase A
    A = read_json(os.path.join(RES, 'phaseA.json'))
    R.add_str('kappa', A['kappa'], 'results/phaseA.json:kappa')
    R.add_str('kappa_closed_diff', '%.4g' % A['kappa_closed_form_diff'],
              'results/phaseA.json:kappa_closed_form_diff')
    R.add('root_number_residual', Decimal(A['root_number_identity_residual']), 999,
          'results/phaseA.json:root_number_identity_residual')
    R.add_str('fe_worst_rel_60', '%.4g' % A['fe_worst_rel_base'],
              'results/phaseA.json:fe_worst_rel_base')
    R.add_str('fe_worst_rel_100', '%.4g' % A['fe_worst_rel_high'],
              'results/phaseA.json:fe_worst_rel_high')
    R.add_str('sigma0_strip_bound', '%.12g' % A['strip']['bound'],
              'results/phaseA.json:strip.bound (zeta(2)-1)')
    off = A['offline']
    R.add_str('beta1', off['beta'], 'results/phaseA.json:offline.beta')
    R.add_str('gamma1', off['gamma'], 'results/phaseA.json:offline.gamma')
    R.add_str('delta1', off['delta'], 'results/phaseA.json:offline.delta')
    R.add_str('seed_distance', '%.2g' % off['seed_distance'],
              'results/phaseA.json:offline.seed_distance')
    R.add_str('rho1_prec_stability', '%.3g' % off['precision_stability'],
              'results/phaseA.json:offline.precision_stability')
    R.add_str('rho1_absF', '%.2g' % off['absF'], 'results/phaseA.json:offline.absF')
    R.add_str('rho1_newton', '%.2g' % off['newton_step'],
              'results/phaseA.json:offline.newton_step')
    for k in ('rel_AB', 'rel_AC'):
        worst = max(r[k] for r in A['routes'])
        R.add_str('route_%s_worst' % k[4:], '%.4g' % worst,
                  'results/phaseA.json:routes.%s (max)' % k)
    R.add_str('route_D_diff', '%.4g' % A['route_D_diff'],
              'results/phaseA.json:route_D_diff')

    # ------------------------------------------------------------------ Phase B
    B = read_json(os.path.join(RES, 'phaseB.json'))
    R.add_str('b_inverse_route_diff', '%.4g' % B['b']['inverse_route_diff'],
              'results/phaseB.json:b.inverse_route_diff')
    R.add_str('G_at_vs_quad', '%.4g' % B['G_at_vs_quad'],
              'results/phaseB.json:G_at_vs_quad')
    sr = B['scalar_routes']
    R.add_str('arch_routes_worst_diff', '%.4g' % max(r['diff'] for r in sr),
              'results/phaseB.json:scalar_routes.diff (max)')
    R.add_str('arch_routes_cut_bound', '%.4g' % max(r['cut'] for r in sr),
              'results/phaseB.json:scalar_routes.cut (max)')
    cut26 = [r for r in B['cutoff'] if r['L'] == '2.6'][0]
    R.add_int('nmax_26', cut26['nmax'], 'results/phaseB.json:cutoff[L=2.6].nmax')
    R.add_int('nterms_26', cut26['count'], 'results/phaseB.json:cutoff[L=2.6].count')
    R.add_str('cutoff_margin_26', '%.6g' % cut26['margin'],
              'results/phaseB.json:cutoff[L=2.6].margin')
    # the b_n printed in the note
    bc = read_csv(os.path.join(DATA, 'dh_coefficients.csv'))
    bmap = {int(r['n']): r['b_n'] for r in bc}
    for n in (2, 3, 6):
        R.add_str('b_%d' % n, bmap[n], 'data/dh_coefficients.csv:n=%d' % n)
    npp = [int(r['n']) for r in bc
           if int(r['is_prime_power']) == 0 and Decimal(r['b_n']) != 0
           and int(r['n']) < 50]
    R.add_int('n_nonpp_below50', len(npp),
              'data/dh_coefficients.csv (count of non-prime-powers < 50 with b_n != 0)')

    # ------------------------------------------------------------------ census
    C = read_json(os.path.join(RES, 'census.json'))
    rec = C['reconciliation']
    R.add_int('census_total', rec['contour'], 'results/census.json:reconciliation.contour')
    R.add_int('census_online', rec['online'], 'results/census.json:reconciliation.online')
    R.add_int('census_offline', rec['offline'],
              'results/census.json:reconciliation.offline')
    R.add_str('census_integrality', '%.3g' % max(r['integrality'] for r in C['ladder']),
              'results/census.json:ladder.integrality (max)')
    for i, o in enumerate(C['offline'], 1):
        R.add_str('O%d_beta' % i, o['beta'], 'results/census.json:offline[%d].beta' % (i - 1))
        R.add_str('O%d_gamma' % i, o['gamma'], 'results/census.json:offline[%d].gamma' % (i - 1))
        R.add_str('O%d_delta' % i, o['delta'], 'results/census.json:offline[%d].delta' % (i - 1))
        R.add_str('O%d_absF' % i, '%.3g' % o['absF'], 'results/census.json:offline[%d].absF' % (i - 1))
        R.add_int('O%d_winding' % i, o['multiplicity'],
                  'results/census.json:offline[%d].multiplicity' % (i - 1))
    for row in C['ladder']:
        if row['T'] == 200.0 and 'rvm_200' not in R:
            R.add_str('rvm_200', '%.6g' % row['rvm'], 'results/census.json:ladder[T=200].rvm')

    # ------------------------------------------------------------------ gate
    G = read_json(os.path.join(RES, 'phaseC.json'))
    gc = read_csv(os.path.join(DATA, 'gate_checks.csv'))
    con = [r for r in gc if r['diff_contour'] not in ('', 'nan')]
    R.add_str('gate_contour_diff', '%.3g' % max(float(r['diff_contour']) for r in con),
              'data/gate_checks.csv:diff_contour (max of finite)')
    R.add_str('gate_contour_budget', '%.4g' % max(float(r['tail_contour']) for r in con),
              'data/gate_checks.csv:tail_contour (max of finite)')
    R.add_str('gate_zeroside_diff', '%.3g' % max(float(r['diff_zeroside']) for r in gc),
              'data/gate_checks.csv:diff_zeroside (max)')
    R.add_str('gate_zeroside_budget', '%.3g' % max(float(r['tail_zeroside']) for r in gc),
              'data/gate_checks.csv:tail_zeroside (max)')
    R.add_str('ctrlA_residual', '%.3g' % G['controlA']['residual'],
              'results/phaseC.json:controlA.residual')
    R.add_str('ctrlA_QO1', G['controlA']['Q_O1'], 'results/phaseC.json:controlA.Q_O1')
    R.add_str('ctrlB_residual', '%.3g' % G['controlB']['residual'],
              'results/phaseC.json:controlB.residual')
    R.add_str('ctrlB_b6', G['controlB']['b6'], 'results/phaseC.json:controlB.b6')
    R.add_str('ctrlB_h6', G['controlB']['h6'], 'results/phaseC.json:controlB.h6')
    R.add_str('ctrlB_predicted', G['controlB']['predicted'],
              'results/phaseC.json:controlB.predicted')
    sf = read_csv(os.path.join(DATA, 'gate_signflip.csv'))
    R.add_str('signflip_worst', '%.3g' % max(float(r['residual']) for r in sf),
              'data/gate_signflip.csv:residual (max)')
    lq = [r for r in sf if r['flip'] == 'logq'][0]
    R.add_str('signflip_logq', lq['predicted'],
              'data/gate_signflip.csv:flip=logq.predicted')

    # ------------------------------------------------------------------ witnesses
    for tag, fn, wfn in (('w20', 'ball_certificate.txt', 'witness.json'),
                         ('w18', 'ball_certificate_L1.8.txt', 'witness_smallest_L.json')):
        bf = ball_fields(os.path.join(CERT, fn))
        w = read_json(os.path.join(CERT, wfn))
        R.add_str('%s_L' % tag, w['L'], 'certificates/%s:L' % wfn)
        R.add_int('%s_N' % tag, len(w['idx']), 'certificates/%s:len(idx)' % wfn)
        R.add_int('%s_level' % tag, w['level'], 'certificates/%s:level' % wfn)
        for key, lab in (('Q_F(g_*)/||g||^2', 'QF'), ('Q_{O_1}/||g||^2', 'QO1'),
                         ('Q_rest /||g||^2', 'Qrest')):
            mid, rad = bf[key]
            R.add_str('%s_%s' % (tag, lab), mid, 'certificates/%s:%s midpoint' % (fn, key))
            R.add_str('%s_%s_rad' % (tag, lab), '%.3g' % Decimal(rad),
                      'certificates/%s:%s radius' % (fn, key))
        R.add_str('%s_radius' % tag, '%.3g' % Decimal(bf['radius'][0]),
                  'certificates/%s:radius' % fn)
        txt = open(os.path.join(CERT, fn)).read()
        mrec = re.search(r'additive reconciliation:.*?\[\s*\+/-\s*'
                         r'([0-9.]+(?:[eE][+-]?[0-9]+)?)\s*\]', txt)
        R.add_str('%s_recon' % tag, '%.3g' % Decimal(mrec.group(1)),
                  'certificates/%s:additive reconciliation radius' % fn)
        R.add_int('%s_cutoff' % tag,
                  int(re.search(r'= (\d+), integer-decided',
                                open(os.path.join(CERT, fn)).read()).group(1)),
                  'certificates/%s:arithmetic cutoff' % fn)
        R.add_int('%s_nterms' % tag,
                  int(re.search(r'b_n != 0\): (\d+)',
                                open(os.path.join(CERT, fn)).read()).group(1)),
                  'certificates/%s:arithmetic terms used' % fn)
        R.add_str('%s_gram_dev' % tag, '%.3g' % w['gram_max_dev_from_I'],
                  'certificates/%s:gram_max_dev_from_I' % wfn)

    # ------------------------------------------------------------------ scan
    ws = read_csv(os.path.join(DATA, 'window_scan.csv'))
    for r in ws:
        R.add_str('scan_%s_%s_%s' % (r['L'].replace('.', ''), r['sector'], r['level']),
                  r['lambda_min'],
                  'data/window_scan.csv:L=%s,%s,level=%s' % (r['L'], r['sector'], r['level']))
    R.add_str('scan_rayleigh_worst', '%.3g' % max(float(r['rayleigh_check']) for r in ws),
              'data/window_scan.csv:rayleigh_check (max)')
    R.add_int('scan_rows', len(ws), 'data/window_scan.csv (row count)')
    R.add_int('scan_negative_rows', sum(1 for r in ws if int(r['negative'])),
              'data/window_scan.csv (rows with negative=1)')

    # ------------------------------------------------------------------ remaining form
    rf = read_csv(os.path.join(DATA, 'remaining_form.csv'))
    for r in rf:
        R.add_str('rest_%s_%s' % (r['L'].replace('.', ''), r['sector']),
                  r['lambda_min_Qrest'],
                  'data/remaining_form.csv:L=%s,%s' % (r['L'], r['sector']))
    R.add_int('rest_n_rows', len(rf), 'data/remaining_form.csv (row count)')
    R.add_int('rest_n_R4', sum(1 for r in rf if r['classification'] == 'R4'),
              'data/remaining_form.csv (rows classified R4)')

    # ------------------------------------------------------------------ masking
    mk = read_csv(os.path.join(DATA, 'masking.csv'))
    for r in mk:
        if r['kind'] == 'fixed_witness_removal':
            R.add_str('mask_%s' % r['orbit'], r['value'],
                      'data/masking.csv:fixed_witness_removal %s' % r['orbit'])
        if r['kind'] == 'reoptimised_removal' and r['orbit'] == 'O_1':
            R.add_str('reopt_without_O1', r['value'],
                      'data/masking.csv:reoptimised_removal O_1')
        if r['kind'] == 'synthetic_hybrid':
            R.add_str('hybrid_fixed', r['value'], 'data/masking.csv:synthetic_hybrid value')
            R.add_str('hybrid_reopt', r['after'], 'data/masking.csv:synthetic_hybrid after')

    # ------------------------------------------------------------------ sensitivity
    sv = read_csv(os.path.join(DATA, 'sensitivity.csv'))
    for r in sv:
        if r['kind'] == 'K' and r['L'] == '2.0':
            R.add_str('K_even_20', r['a'], 'data/sensitivity.csv:K L=2.0 K_even')
            R.add_str('K_odd_20', r['b'], 'data/sensitivity.csv:K L=2.0 K_odd')
            R.add_str('K_sum_resid_20', '%.3g' % Decimal(r['c']),
                      'data/sensitivity.csv:K L=2.0 (K_even+K_odd-2L^3/3)')
        if r['kind'] == 'overlap':
            R.add_str('overlap_witness', r['a'], 'data/sensitivity.csv:overlap a')
            R.add_str('overlap_Qrest_gs', '%.3g' % Decimal(r['c']),
                      'data/sensitivity.csv:overlap c')
    pert = [r for r in sv if r['kind'] == 'perturbation']
    for r in pert:
        R.add_str('pert_t%s' % r['t'].replace('.', ''), r['a'],
                  'data/sensitivity.csv:perturbation t=%s' % r['t'])
    smallest = min(pert, key=lambda r: float(r['t']) if float(r['t']) > 0 else 9)
    R.add_str('pert_resid_over_t4', '%.4g' % (Decimal(smallest['d'])
                                              / Decimal(smallest['t']) ** 4),
              'data/sensitivity.csv:perturbation residual/t^4 at t=%s' % smallest['t'])

    # ------------------------------------------------------------------ derived in note
    import math
    gamma1 = R['gamma1']['value']
    delta1 = R['delta1']['value']
    Lpred = Decimal(repr(0.5 * math.log(5 * float(gamma1) / (2 * math.pi))))
    R.add_str('L_pred', '%.9g' % Lpred,
              'derived: (1/2) log(5 gamma_1 / 2 pi) from gamma_1')
    Kodd = R['K_odd_20']['value']
    fourKd2 = 4 * Kodd * delta1 * delta1
    R.add_str('fourK_delta2', '%.9g' % fourKd2,
              'derived: 4 K_odd(2.0,gamma_1) delta_1^2')
    QO1 = abs(R['w20_QO1']['value'])
    R.add_str('ratio_4Kd2_over_QO1', '%.6g' % (fourKd2 / QO1),
              'derived: 4K delta^2 / |Q_{O_1}(g_*)|')
    # the collapsed quartet value Q_{O_1}(0) = Q_t at t=0 minus Q_rest
    Qt0 = R['pert_t00']['value']
    Qrest = R['w20_Qrest']['value']
    QO1_0 = Qt0 - Qrest
    R.add_str('QO1_collapsed', '%.9g' % QO1_0,
              'derived: Q_t(t=0) - Q_rest(g_*) = Q_{O_1}(0) = 4|F(gamma_1)|^2')
    drop = QO1_0 - R['w20_QO1']['value']
    R.add_str('QO1_drop', '%.9g' % drop,
              'derived: Q_{O_1}(0) - Q_{O_1}(delta_1), the actual drop')
    R.add_str('ratio_4Kd2_over_drop', '%.6g' % (fourKd2 / drop),
              'derived: 4K delta^2 / actual drop')
    # the root paper's delta^4 remainder allowance at L=2, delta=delta_1
    Lw = R['w20_L']['value']
    Ld = Lw * delta1
    rem = (Decimal('0.4') + Decimal(repr(16 / 3)) * Decimal(repr(math.exp(2 * float(Ld))))
           ) * Lw ** 5 * delta1 ** 4
    R.add_str('L_delta', '%.6g' % Ld, 'derived: L * delta_1 at L = 2.0')
    R.add_str('thm8_remainder', '%.6g' % rem,
              'derived: (2/5 + (16/3)e^{2 L delta}) L^5 delta^4, root paper Thm 8')
    R.add_str('ratio_discrepancy', '%.4g' % abs(fourKd2 - drop),
              'derived: |4K delta^2 - actual drop|')
    # masking fraction
    R.add_str('masking_fraction_pct', '%.4g' % (100 * R['w20_Qrest']['value'] / QO1),
              'derived: 100 * Q_rest(g_*) / |Q_{O_1}(g_*)|')
    QF18 = R['w18_QF']['value']
    R.add_str('masking_fraction_pct_18',
              '%.6g' % (100 * R['w18_Qrest']['value'] / abs(R['w18_QO1']['value'])),
              'derived: 100 * Q_rest / |Q_{O_1}| at L = 1.8')
    # --- strip, from src/census.py's SIGMA0 (derived in DERIVATION.md section 1)
    cen = open(os.path.join(HERE, 'src', 'census.py')).read()
    s0 = re.search(r"SIGMA0 = '([0-9.]+)'", cen).group(1)
    R.add_str('sigma0', s0, 'src/census.py:SIGMA0')
    R.add_str('strip_halfwidth', '%.6g' % (Decimal(s0) - Decimal('0.5')),
              'derived: sigma_0 - 1/2, the half-width of the zero-containing strip')
    # --- the smallest |Lambda_f| among the functional-equation test points
    fe081 = [r for r in A['fe'] if r['sigma'] == '0.81'][0]
    R.add_str('fe_scale_0p81', '%.3g' % fe081['scale'],
              "results/phaseA.json:fe[sigma=0.81,t=85.7].scale")
    R.add_str('fe_min_scale', '%.3g' % min(r['scale'] for r in A['fe']),
              'results/phaseA.json:fe.scale (min)')
    # --- magnitudes, for sentences that say "negative by X"
    for src_name, dst in (('w20_QF', 'w20_QF_abs'), ('w20_QO1', 'w20_QO1_abs'),
                          ('w18_QF', 'w18_QF_abs'), ('w18_QO1', 'w18_QO1_abs'),
                          ('scan_26_even_d_3', 'scan_26_even_d_3_abs'),
                          ('scan_20_odd_3', 'scan_20_odd_3_abs'),
                          ('scan_18_odd_3', 'scan_18_odd_3_abs')):
        R.add(dst, abs(R[src_name]['value']), R[src_name]['sig'],
              'derived: |%s|' % src_name)
    R.add_str('K_even_over_odd', '%.6g' % (R['K_even_20']['value'] / Kodd),
              'derived: K_even(2.0,gamma_1) / K_odd(2.0,gamma_1)')
    R.add_str('ratio_Lpred_over_L18', '%.6g' % (Lpred / R['w18_L']['value']),
              'derived: L_pred / 1.8 (smallest window with a CERTIFIED negative witness)')
    R.add_str('ratio_Lpred_over_L20', '%.6g' % (Lpred / Lw),
              'derived: L_pred / 2.0 (smallest window with a MACROSCOPIC witness)')
    R.add_str('remainder_over_discrepancy',
              '%.4g' % (rem / abs(fourK_delta2_minus_drop(fourKd2, drop))),
              'derived: Thm 8 delta^4 allowance / |4K delta^2 - actual drop|')
    # the five-order depth ratio between the two sectors at L = 1.8, level 3
    e18 = abs(R['scan_18_even_d_3']['value'])
    o18 = abs(R['scan_18_odd_3']['value'])
    R.add_str('sector_depth_ratio_18', '%.4g' % (o18 / e18),
              'derived: |lambda_min(odd)| / |lambda_min(even)| at L=1.8, level 3')
    R.add_str('scan_collapse_ratio', '%.4g' % (
        abs(R['scan_10_even_d_3']['value']) / abs(R['scan_16_even_d_3']['value'])),
        'derived: lambda_min(L=1.0,even,3) / lambda_min(L=1.6,even,3)')
    return R


def fourK_delta2_minus_drop(a, b):
    return a - b


if __name__ == '__main__':
    R = build()
    for k in sorted(R):
        v = R[k]
        print('%-28s %-46s sig=%-4s %s' % (k, v['value'], v['sig'], v['source']))
