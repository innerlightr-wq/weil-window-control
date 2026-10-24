r"""Derive every quantity the Siegel note prints, from this extension's data files only.

Nothing here reads the note.  `verify_note_numbers.py` imports `build()` and checks the
note against it.  Each entry is (value, sig_digits_supported, source), where
sig_digits_supported is the number of significant digits the SOURCE carries -- the note
may print that many or fewer, correctly rounded, never more.

Design copied from the companion DH note's note/registry.py (same repository), including
the exponent-preserving `trunc` and the refusal to store truncations in place of full
values.
"""
import csv
import json
import os
import re
from decimal import Decimal, getcontext

getcontext().prec = 120

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, 'data')
RES = os.path.join(HERE, 'results')
# the companion DH extension, for the one cross-note quantity
DHX = os.path.join(os.path.dirname(HERE), 'davenport-heilbronn-2026-10-23')


def _sig(s):
    s = s.strip().lstrip('+-')
    s = re.sub(r'[eE].*$', '', s)
    s = s.replace('.', '').lstrip('0')
    return len(s) or 1


class Registry(dict):
    def add(self, name, value, sig, source):
        if name in self:
            raise KeyError('duplicate registry key %s' % name)
        self[name] = dict(value=Decimal(str(value)), sig=sig, source=source)

    def add_str(self, name, s, source, sig=None):
        self.add(name, Decimal(str(s)), sig if sig is not None else _sig(str(s)), source)

    def add_int(self, name, n, source):
        self.add(name, Decimal(int(n)), 999, source)


def read_csv(p):
    with open(p) as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith('#')))


def read_json(p):
    with open(p) as fh:
        return json.load(fh)


def build():
    R = Registry()

    # ------------------------------------------------------- check (b): the 150 cells
    vs = read_csv(os.path.join(DATA, 'verify_siegel_parity.csv'))
    R.add_int('b_cells', len(vs), 'data/verify_siegel_parity.csv (row count)')
    R.add_str('b_worst_rel_err', '%.3g' % max(float(r['rel_err']) for r in vs),
              'data/verify_siegel_parity.csv:rel_err (max)')
    R.add_str('b_worst_rank1_resid',
              '%.3g' % max(float(r['rank1_residual']) for r in vs),
              'data/verify_siegel_parity.csv:rank1_residual (max)')
    R.add_str('b_worst_other_eig',
              '%.3g' % max(float(r['other_eigs_max']) for r in vs),
              'data/verify_siegel_parity.csv:other_eigs_max (max)')
    R.add_str('b_worst_cross_integral',
              '%.3g' % max(abs(float(r['cross_integral'])) for r in vs),
              'data/verify_siegel_parity.csv:cross_integral (max abs)')
    ev = [float(r['eigenvalue']) for r in vs if r['sector'] != 'odd']
    od = [float(r['eigenvalue']) for r in vs if r['sector'] == 'odd']
    R.add_int('b_even_cells', len(ev), 'data/verify_siegel_parity.csv (even rows)')
    R.add_int('b_odd_cells', len(od), 'data/verify_siegel_parity.csv (odd rows)')
    R.add_str('b_even_min', '%.6g' % min(ev), 'data/verify_siegel_parity.csv even eig min')
    R.add_str('b_even_max', '%.6g' % max(ev), 'data/verify_siegel_parity.csv even eig max')
    R.add_str('b_odd_min', '%.6g' % min(od), 'data/verify_siegel_parity.csv odd eig min')
    R.add_str('b_odd_max', '%.6g' % max(od), 'data/verify_siegel_parity.csv odd eig max')
    R.add_int('b_sign_ok', sum(1 for r in vs if r['sign_ok'] == 'True'),
              'data/verify_siegel_parity.csv:sign_ok (count True)')
    # how many distinct eigenvalues per (L, sector, delta) group -- the chi-independence
    groups = {}
    for r in vs:
        groups.setdefault((r['L'], r['sector'], r['delta']), set()).add(
            r['eigenvalue'])
    R.add_int('b_max_distinct_per_group', max(len(v) for v in groups.values()),
              'data/verify_siegel_parity.csv grouped by (L,sector,delta)')
    R.add_int('b_n_groups', len(groups),
              'data/verify_siegel_parity.csv distinct (L,sector,delta) groups')

    # ------------------------------------------------------- checks (a) and (c)
    j = read_json(os.path.join(RES, 'verify_siegel_parity.json'))['results']
    R.add_str('a_quartet_twice_pair_maxdiff',
              '%.3g' % float(j['a_quartet_is_twice_pair_maxdiff']),
              'results/verify_siegel_parity.json:a_quartet_is_twice_pair_maxdiff')
    R.add_str('K_identity_maxdiff', '%.3g' % float(j['b_K_identity_maxdiff']),
              'results/verify_siegel_parity.json:b_K_identity_maxdiff')
    cs = j['c']
    R.add_int('c_cells', len(cs), 'results/verify_siegel_parity.json:c (cell count)')
    R.add_str('c_chi_realpair_max_entry',
              '%.8g' % float(j['c_chi_realpair_max_entry']),
              'results/verify_siegel_parity.json:c_chi_realpair_max_entry')
    R.add_str('c_worst_B_minus_Pole',
              '%.3g' % max(float(v['max_abs_B_minus_Pole']) for v in cs.values()),
              'results/verify_siegel_parity.json:c.*.max_abs_B_minus_Pole (max)')
    R.add_str('c_min_B_plus_Pole',
              '%.4g' % min(float(v['max_abs_B_plus_Pole']) for v in cs.values()),
              'results/verify_siegel_parity.json:c.*.max_abs_B_plus_Pole (min)')
    R.add_str('c_max_B_plus_Pole',
              '%.4g' % max(float(v['max_abs_B_plus_Pole']) for v in cs.values()),
              'results/verify_siegel_parity.json:c.*.max_abs_B_plus_Pole (max)')
    for key, lab in ((k, k.replace('|', '_').replace('.', '')) for k in cs):
        R.add_str('c_%s_minus' % lab, '%.3g' % float(cs[key]['max_abs_B_minus_Pole']),
                  'results/verify_siegel_parity.json:c[%s].max_abs_B_minus_Pole' % key)
        R.add_str('c_%s_plus' % lab, '%.4g' % float(cs[key]['max_abs_B_plus_Pole']),
                  'results/verify_siegel_parity.json:c[%s].max_abs_B_plus_Pole' % key)
        R.add_int('c_%s_polesign' % lab, cs[key]['pole_sign'],
                  'results/verify_siegel_parity.json:c[%s].pole_sign' % key)

    # ------------------------------------------------------- the K identity table
    k4 = read_csv(os.path.join(DATA, 'phase4b_sector_extremal.csv'))
    R.add_int('K_rows', len(k4), 'data/phase4b_sector_extremal.csv (row count)')
    R.add_str('K_worst_abs_err', '%.3g' % max(float(r['abs_err']) for r in k4),
              'data/phase4b_sector_extremal.csv:abs_err (max)')
    for r in k4:
        if r['gamma'] == '0':
            tag = r['L'].replace('.', '')
            R.add_str('K_odd_0_L%s' % tag, r['K_odd'],
                      'data/phase4b_sector_extremal.csv:L=%s,gamma=0,K_odd' % r['L'])
            R.add_str('K_even_0_L%s' % tag, r['K_even'],
                      'data/phase4b_sector_extremal.csv:L=%s,gamma=0,K_even' % r['L'])
            R.add_str('two_L3_3_L%s' % tag, r['two_L3_3'],
                      'data/phase4b_sector_extremal.csv:L=%s,two_L3_3' % r['L'])

    # ------------------------------------------------------- function field
    ff = read_json(os.path.join(RES, 'ff_parity_split.json'))
    rows = ff['rows']
    R.add_int('ff_cells', len(rows), 'results/ff_parity_split.json:rows (count)')
    R.add_int('ff_failures', len(ff['failures']),
              'results/ff_parity_split.json:failures (count)')
    nondeg = [r for r in rows if r['rho'] != '1']
    R.add_int('ff_nondegenerate', len(nondeg),
              'results/ff_parity_split.json rows with rho != 1')
    for flag in ('split_exact', 'vplus_is_eigvec', 'vminus_is_eigvec',
                 'vplus_palindromic', 'vminus_antipalindromic',
                 'lam_plus_positive', 'lam_minus_nonpositive',
                 'brief_split_identity'):
        R.add_int('ff_%s' % flag, sum(1 for r in rows if r[flag]),
                  'results/ff_parity_split.json:%s (count true)' % flag)
    for flag in ('brief_upw_is_eigvec', 'brief_upw_palindromic'):
        R.add_int('ff_%s_true' % flag, sum(1 for r in rows if r[flag]),
                  'results/ff_parity_split.json:%s (count true)' % flag)
    from fractions import Fraction as _F
    R.add_int('ff_max_denom_digits',
              max(len(str(_F(r['lam_minus']).denominator)) for r in rows),
              'results/ff_parity_split.json: digits in the largest lam_minus denominator')
    R.add_int('ff_rho_values', len(set(r['rho'] for r in rows)),
              'results/ff_parity_split.json distinct rho')
    R.add_int('ff_R_values', len(set(r['R'] for r in rows)),
              'results/ff_parity_split.json distinct R')

    # ------------------------------------------------------- onset, odd sector
    on = read_csv(os.path.join(DATA, 'phase5_onset_odd.csv'))
    R.add_int('onset_rows', len(on), 'data/phase5_onset_odd.csv (row count)')
    labs, dels = [], []
    for r in on:
        if r['label'] not in labs:
            labs.append(r['label'])
        if r['delta'] not in dels:
            dels.append(r['delta'])
    R.add_int('onset_characters', len(labs), 'data/phase5_onset_odd.csv distinct label')
    R.add_int('onset_deltas', len(dels), 'data/phase5_onset_odd.csv distinct delta')
    R.add_int('onset_N', int(on[0]['N']), 'data/phase5_onset_odd.csv:N')
    Ls = sorted(set(float(r['L']) for r in on))
    R.add_int('onset_n_L', len(Ls), 'data/phase5_onset_odd.csv distinct L')
    R.add_str('onset_L_min', '%g' % min(Ls), 'data/phase5_onset_odd.csv min L')
    R.add_str('onset_L_max', '%g' % max(Ls), 'data/phase5_onset_odd.csv max L')
    for lab in labs:
        for d in dels:
            hits = [float(r['L']) for r in on
                    if r['label'] == lab and r['delta'] == d and r['detected'] == 'yes']
            key = 'onset_%s_d%s' % (lab.replace('_', '').replace('-', 'm'),
                                    d.replace('.', ''))
            if hits:
                R.add_str(key, '%g' % min(hits),
                          'data/phase5_onset_odd.csv: min L with detected=yes '
                          'for %s, delta=%s' % (lab, d))
    R.add_int('onset_detected_rows', sum(1 for r in on if r['detected'] == 'yes'),
              'data/phase5_onset_odd.csv:detected (count yes)')
    # the worst L at which every character/delta has detected
    worst = max(min(float(r['L']) for r in on
                    if r['label'] == lab and r['delta'] == d and r['detected'] == 'yes')
                for lab in labs for d in dels
                if any(r['label'] == lab and r['delta'] == d and r['detected'] == 'yes'
                       for r in on))
    R.add_str('onset_all_by_L', '%g' % worst,
              'data/phase5_onset_odd.csv: max over (label,delta) of the onset L')
    # zeta's near-singular odd baseline at L=0.8
    z08 = [r for r in on if r['label'] == 'zeta' and r['L'] == '0.8']
    if z08:
        R.add_str('zeta_odd_baseline_08', '%.6g' % float(z08[0]['lam_unperturbed']),
                  'data/phase5_onset_odd.csv: zeta, L=0.8, lam_unperturbed')
    for lab in ('chi_-4', 'chi_5'):
        rr = [r for r in on if r['label'] == lab and r['L'] == '0.8']
        if rr:
            R.add_str('%s_odd_baseline_08' % lab.replace('_', '').replace('-', 'm'),
                      '%.3g' % float(rr[0]['lam_unperturbed']),
                      'data/phase5_onset_odd.csv: %s, L=0.8, lam_unperturbed '
                      '(N = %s)' % (lab, rr[0]['N']))

    # ------------------------------------------------------- phase 5, even sector
    p5 = read_csv(os.path.join(DATA, 'phase5_siegel.csv'))
    evn = [r for r in p5 if r['sector'] != 'odd']
    odr = [r for r in p5 if r['sector'] == 'odd']
    R.add_int('p5_rows', len(p5), 'data/phase5_siegel.csv (row count)')
    R.add_int('p5_even_rows', len(evn), 'data/phase5_siegel.csv (even rows)')
    R.add_int('p5_odd_rows', len(odr), 'data/phase5_siegel.csv (odd rows)')
    R.add_int('p5_even_nonneg', sum(1 for r in evn if r['response_sign'] == '>=0'),
              "data/phase5_siegel.csv even rows with response_sign '>=0'")
    R.add_int('p5_even_detected', sum(1 for r in evn if r['detected'] == 'yes'),
              'data/phase5_siegel.csv even rows with detected=yes')
    R.add_int('p5_odd_neg', sum(1 for r in odr if r['response_sign'] == '<0'),
              "data/phase5_siegel.csv odd rows with response_sign '<0'")
    R.add_int('p5_odd_detected', sum(1 for r in odr if r['detected'] == 'yes'),
              'data/phase5_siegel.csv odd rows with detected=yes')

    # ------------------------------------------------------- companion-note quantity
    sv = read_csv(os.path.join(DHX, 'data', 'sensitivity.csv'))
    for r in sv:
        if r['kind'] == 'K' and r['L'] == '2.0':
            ke, ko = Decimal(r['a']), Decimal(r['b'])
            R.add_str('dh_K_even_over_odd', '%.6g' % (ke / ko),
                      'companion DH extension data/sensitivity.csv: '
                      'K_even/K_odd at L=2.0, gamma_1')
            break
    return R


if __name__ == '__main__':
    R = build()
    for k in sorted(R):
        print('%-34s %-28s sig=%-4s %s' % (k, R[k]['value'], R[k]['sig'],
                                           R[k]['source']))
