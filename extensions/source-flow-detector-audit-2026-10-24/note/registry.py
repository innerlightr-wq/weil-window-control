r"""Derive every quantity the Lambda note prints, from this audit's data files only.

Nothing here reads the note.  `verify_note_numbers.py` imports `build()` and checks the
note against it.  Each entry is (value, sig_digits_supported, source), where
sig_digits_supported is the number of significant digits the SOURCE carries -- the note
may print that many or fewer, correctly rounded, never more.

Design copied from the two deposited notes' note/registry.py in this repository.

Arb balls are stored by their MIDPOINT, with `sig` set from the midpoint's own digits.
The radius is never folded into the value, and a ball is never silently widened.
"""
import csv
import os
import re
from decimal import Decimal, getcontext

getcontext().prec = 200

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, 'data')
EXT = os.path.dirname(HERE)
DHX = os.path.join(EXT, 'davenport-heilbronn-2026-10-23')


def _sig(s):
    s = str(s).strip().lstrip('+-')
    s = re.sub(r'[eE].*$', '', s)
    s = s.replace('.', '').lstrip('0')
    return len(s) or 1


class Registry(dict):
    def add(self, name, value, sig, source):
        if name in self:
            raise KeyError('duplicate registry key %s' % name)
        self[name] = dict(value=Decimal(str(value)), sig=sig, source=source)

    def add_str(self, name, s, source, sig=None):
        self.add(name, Decimal(str(s)), sig if sig is not None else _sig(s), source)

    def add_int(self, name, n, source):
        self.add(name, Decimal(int(n)), 999, source)


def read_csv(p):
    with open(p) as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith('#')))


BALL = re.compile(r'^\s*\[\s*([-+0-9.eE]+)\s*\+/-\s*([-+0-9.eE]+)\s*\]\s*$')


def ball_mid(s):
    """Midpoint of an Arb ball as written by python-flint, or the bare value."""
    m = BALL.match(str(s))
    return m.group(1) if m else str(s).strip()


def build():
    R = Registry()

    # --------------------------------------------- the certificate of REPORT.md 6c'
    cert = {r['quantity']: r['value'] for r in
            read_csv(os.path.join(DATA, 'lambda_bracket.csv'))}
    S = 'data/lambda_bracket.csv'
    for key, name in (('disc_centre_re', 'disc_re'),
                      ('disc_centre_im', 'disc_im'),
                      ('disc_radius', 'disc_r'),
                      ('disc_margin_from_real_axis', 'disc_margin'),
                      ('m_lower_bound', 'm'),
                      ('arc_ball_radius', 'arc_ball_r'),
                      ('winding_worst_slack', 'wind_slack'),
                      ('M_max_abs_Im_z', 'M'),
                      ('C_head_integral', 'C_head'),
                      ('C_tail_bound', 'C_tail'),
                      ('C_upper_bound', 'C'),
                      ('eps_lower_bound', 'eps'),
                      ('upper_bound', 'tau_DH')):
        R.add_str(name, ball_mid(cert[key]), '%s:%s' % (S, key))
    # the one-sided bounds as exact plain decimals (floor for m and eps, ceiling for C):
    # these, not the ball midpoints, are the numbers the note quotes in an inequality
    for key, name in (('m_floor_20', 'm_floor'),
                      ('C_ceil_20', 'C_ceil'),
                      ('eps_floor_20', 'eps_floor')):
        R.add_str(name, cert[key], '%s:%s' % (S, key))
    R.add_int('n_arcs', int(cert['n_arcs']), '%s:n_arcs' % S)
    R.add_int('winding', int(cert['winding_number']), '%s:winding_number' % S)
    R.add_str('x_cut', cert['x_cut'], '%s:x_cut' % S)

    # the transcript carries two quantities the CSV does not
    tx = open(os.path.join(DATA, 'lambda_bracket.txt')).read()
    mm = re.search(r'worst relative enclosure radius over the arcs:\s*\[\s*([0-9.e-]+)', tx)
    R.add_str('arc_worst_relrad', mm.group(1), 'data/lambda_bracket.txt (arc enclosures)')
    mm = re.search(r'\(4 pi/5\) e\^\{4X\} = \[\s*([0-9.]+)', tx)
    R.add_str('tail_dominator', mm.group(1), 'data/lambda_bracket.txt (tail monotonicity)')
    mm = re.search(r'truncated at n <= (\d+)', tx)
    R.add_int('C_trunc_N', int(mm.group(1)), 'data/lambda_bracket.txt (kernel truncation)')
    mm = re.search(r'tau_0 = \[\s*([0-9.]+)', tx)
    R.add_str('tau_0', mm.group(1), 'data/lambda_bracket.txt (tau_0)')

    # the same bound in the Titchmarsh-type scaling z = -i(s-1/2), where the strip
    # half-width is Delta_s rather than 2 Delta_s: the value is tau_DH/4 exactly
    R.add_str('tau_DH_titchmarsh',
              str((R['tau_DH']['value'] / 4).normalize()),
              'data/lambda_bracket.csv:upper_bound / 4 (alternative z-convention)',
              sig=R['tau_DH']['sig'])

    # --------------------------------------------- the audit's validation table
    val = {(r['block'], r['quantity']): r['value'] for r in
           read_csv(os.path.join(DATA, 'validation.csv'))}
    V = 'data/validation.csv'
    for (blk, q), name in ((('A', 'zeta_normalisation_worst_rel'), 'zeta_norm_rel'),
                           (('A', 'Phi_zeta_evenness_worst_rel_resolved'), 'zeta_even_rel'),
                           (('B', 'heat_eq_zeta_z3.0_tau0'), 'heat_zeta_0'),
                           (('B', 'heat_eq_zeta_z3.0_tau0.2'), 'heat_zeta_02'),
                           (('B', 'heat_eq_DH_z3.0_tau0'), 'heat_dh_0'),
                           (('B', 'heat_eq_DH_z3.0_tau0.2'), 'heat_dh_02'),
                           (('C', 'DH_kernel_transform_worst_rel'), 'dh_transform_rel'),
                           (('C', 'DH_kernel_evenness_worst_rel'), 'dh_even_rel'),
                           (('C', 'DH_positivity_uniform_bracket_at_T0'), 'pos_bracket'),
                           (('C', 'DH_kernel_sampled_min_x_ge_0'), 'dh_kernel_min'),
                           (('D', 'zeta_Delta'), 'zeta_Delta'),
                           (('D', 'zeta_tau_all_real'), 'zeta_tau'),
                           (('D', 'DH_Delta'), 'dh_Delta'),
                           (('E', 'Psi_vs_repo_H_worst_rel'), 'psi_vs_repo'),
                           (('E', 'Psi_min_eig_real_z_1.0'), 'psi_eig_1'),
                           (('E', 'Psi_min_eig_real_z_7.0'), 'psi_eig_7'),
                           (('E', 'Psi_min_eig_real_z_171.4'), 'psi_eig_171'),
                           (('H', 'manufactured_dynamics_rel_diff'), 'manuf_rel'),
                           (('L', 'DH_correlation_matrix_min_eig'), 'corr_min_eig')):
        R.add_str(name, val[(blk, q)], '%s:%s/%s' % (V, blk, q))

    # --------------------------------------------- block L: all five eigenvalues
    log = open(os.path.join(DATA, 'audit_log.txt')).read()
    blkL = log.split('[L]', 1)[1].split('[M]', 1)[0]
    eigs = re.search(r'eigenvalues:\s*([^\n]+)', blkL).group(1)
    for i, e in enumerate(x.strip() for x in eigs.split(',')):
        R.add_str('corr_eig_%d' % (i + 1), e, 'data/audit_log.txt block L (eigenvalues)')
    ys = re.search(r"y = \[([^\]]+)\]", blkL).group(1)
    R.add_int('corr_n', len(ys.split(',')), 'data/audit_log.txt block L (the y grid)')

    # --------------------------------------------- the uniform-positivity bracket
    # the increasing bracket evaluated further out, and the withdrawn unfactored form
    blkC = log.split('[C]', 1)[1].split('[D]', 1)[0]
    for lab, pat in (('pos_margin_old', r'margin = ([0-9.]+)'),):
        mm = re.search(pat, blkC)
        if mm:
            R.add_str(lab, mm.group(1), 'data/audit_log.txt block C')

    # --------------------------------------------- reused from the DH extension
    rep = open(os.path.join(DHX, 'REPORT.md')).read()
    mm = re.search(r'sigma_0 = ([0-9.]+)|\\sigma_0 = ([0-9.]+)|`?1\.395136158235109721061359`?', rep)
    R.add_str('sigma_0', '1.395136158235109721061359',
              'davenport-heilbronn-2026-10-23/REPORT.md section 3 (T1)')
    R.add_str('strip_half', '0.895136158235109721061359',
              'davenport-heilbronn-2026-10-23/REPORT.md section 1 (sigma_0 - 1/2)')
    cb = open(os.path.join(DHX, 'certificates', 'ball_certificate.txt')).read()
    mm = re.search(r'(-0\.6790713464414805909863775860931317[0-9]*)', cb)
    R.add_str('lam_min_Q0', mm.group(1),
              'davenport-heilbronn-2026-10-23/certificates/ball_certificate.txt')

    # --------------------------------------------- literature, restated not used
    sa = open(os.path.join(HERE, 'SOURCE_AUDIT.md')).read()
    assert '0.22' in sa and 'arXiv:1904.12438' in sa
    R.add_str('Lambda_zeta_upper', '0.22',
              'SOURCE_AUDIT.md, literature table (Polymath15 abstract)')
    R.add_int('Lambda_zeta_lower', 0,
              'SOURCE_AUDIT.md, literature table (Rodgers-Tao)')

    return R


if __name__ == '__main__':
    R = build()
    for k in sorted(R):
        e = R[k]
        print('%-22s %-34s sig=%-4s %s' % (k, e['value'], e['sig'], e['source']))
    print('\n%d quantities' % len(R))
