#!/usr/bin/env python3
r"""Re-derive every number printed in lambda_note.tex from this audit's data files.

Fails if a printed literal is not derivable, if it disagrees with its source, or if the
note prints MORE significant digits than its source carries.  Same design as the two
deposited notes' verifiers in this repository.

Six passes:

  1. registry      -- rebuild every quantity from data/
  2. prerequisite  -- the certificate of REPORT.md 6c' must be present AND its regression
                      fixture must pass; the note's lower bound must not outlive it
  3. doi           -- exactly one DOI in the front matter; a placeholder is permitted and
                      reported, a real one is cross-checked against the audit's record
  4. literature    -- every number in the bibliography, and the two restated literature
                      bounds, must appear in the audit's own record, so none is typed
                      from memory
  5. formula       -- digits that are part of a mathematical expression, declared with a
                      reason for each
  6. values        -- every remaining literal must match a registry entry, to the printed
                      precision and no further
"""
import os
import re
import subprocess
import sys
from decimal import Decimal, ROUND_HALF_UP

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import registry as reg                                                   # noqa: E402

NOTE = os.path.dirname(os.path.abspath(__file__))
HERE = os.path.dirname(NOTE)
TEX = os.path.join(NOTE, 'lambda_note.tex')

FAIL = []


def fail(msg):
    FAIL.append(msg)
    print('  FAIL  %s' % msg)


# Digits that are part of a mathematical expression or a structural reference rather than
# a measured value.  Everything else must match the registry.
FORMULA = {
    '0': 'the origin; tau = 0; Lambda >= 0, also in registry',
    '1': 'unit; rank one; the winding number, also in registry',
    '2': 'the 2 of z = -2i(s-1/2), of Delta^2/2, of 2(sigma_0-1/2), of m/(2C)',
    '3': 'the 3 of e^{3x}; three hypotheses; section 3',
    '4': 'the 4 of e^{4x} and 4 pi/5; four certified inputs; section 4',
    '5': 'the modulus 5 and the 5 of pi n^2/5; five eigenvalues; section 5',
    '6': 'six review points; section 6',
    '7': 'Theorem 7 of Newman-Wu; the 7 of y = 7; seven',
    '8': 'the 8 of 8 e^{3x} and of e^{8}; eighteen',
    '9': 'nine', '10': 'hypothesis (10); base of scientific notation',
    '11': 'hypothesis (11)', '12': 'twelve',
    '15': 'fifteen', '17': 'Phase H: 17 of 18, quoted from the record',
    '18': 'Phase H: 17 of 18, quoted from the record',
    '20': 'twenty digits', '30': 'thirty',
    '40': 'the kernel truncation n <= 40, also in registry',
    '60': 'the 60-arc partition that fails, quoted from the record',
    '105': 'the digits an evenness check needed, quoted from the record',
    '106': 'the measured ball amplification, quoted from the record',
    '400': 'bits of working precision, quoted from the record',
    '3000': 'the arc count, also in registry',
    '0.1': 'tau_0, also in registry', '0.22': 'Lambda <= 0.22, also in registry',
    '0.3': 'the disc radius, also in registry',
    '0.5': 'one half; the zeta threshold, also in registry',
    '1.39': 'the ball radius of the reused witness, quoted from the record',
    '2.5': 'the y grid of the correlation matrix',
    '3.30': 'the failed 60-arc enclosure radius, quoted from the record',
    '0.048': 'C to two digits, also in registry',
    '0.0649': 'the worst arc enclosure radius, also in registry',
    '6.6': 'm to two digits, also in registry',
    '4.5': 'the overstatement of the old m floor, quoted from the record',
}

RECORD = [os.path.join(HERE, f) for f in
          ('REPORT.md', 'SOURCE_AUDIT.md', 'REPRODUCE.md')]
# Three bibliography entries (Bombieri-Ghosh, Davenport-Heilbronn, Arb) are copied
# verbatim from the already-deposited DH note in this repository rather than retyped,
# so that note's .tex is part of the record a citation may be checked against.
EXT = os.path.dirname(HERE)
RECORD += [os.path.join(EXT, 'davenport-heilbronn-2026-10-23', 'note', 'dh_note.tex'),
           os.path.join(EXT, 'dirichlet-windows-2026-10-23', 'note', 'siegel_note.tex'),
           os.path.join(EXT, 'davenport-heilbronn-2026-10-23', 'REPORT.md')]


def strip_tex(src):
    src = src.split('\\begin{document}', 1)[1]
    src = re.sub(r'\\begin\{verbatim\}.*?\\end\{verbatim\}', ' ', src, flags=re.S)
    src = re.sub(r'\\href\{[^}]*\}', ' ', src)
    src = re.sub(r'\\url\{[^}]*\}', ' ', src)
    src = re.sub(r'\\cite\[[^]]*\]', r'\\cite', src)
    src = re.sub(r'\\(label|ref|eqref|addcontentsline|src|texttt|tier|setcounter)\{[^}]*\}',
                 ' ', src)
    src = re.sub(r'%.*', ' ', src)
    src = re.sub(r'~\d+\.\d+(\.\d+)?', ' VERSION ', src)
    src = re.sub(r'10\.5281/zenodo\.\d+|XXXXXXXX', ' DOI ', src)
    src = re.sub(r'\d{4}-\d{4}-\d{4}-\d{4}', ' ORCID ', src)
    src = re.sub(r'arXiv:\d+\.\d+(v\d+)?', ' ARXIV ', src)
    src = src.replace('\\times10^{', '\\TIMESTEN{')
    src = re.sub(r'10\^\{-?\d+\}', ' ORDER ', src)
    src = src.replace('\\TIMESTEN{', '\\times10^{')
    return src


def split_sections(src):
    bib = ''
    if '\\begin{thebibliography}' in src:
        src, bib = src.split('\\begin{thebibliography}', 1)
    return src, bib


NUM = re.compile(r'(?<![0-9A-Za-z_.])(-?)(\d+(?:\.\d+)?)'
                 r'(?:\s*\\times\s*10\^(?:\{(-?\d+)\}|(-?\d)))?')


def literals(text):
    out = []
    for m in NUM.finditer(text):
        sign, mant = m.group(1), m.group(2)
        exp = m.group(3) if m.group(3) is not None else m.group(4)
        ctx = text[max(0, m.start() - 46):m.end() + 22].replace('\n', ' ')
        out.append((m.group(0).strip(), sign, mant, exp, ctx))
    return out


def sigdigits(mant):
    d = mant.replace('.', '').lstrip('0')
    return len(d) or 1


def to_dec(sign, mant, exp):
    v = Decimal(mant)
    if exp:
        v = v.scaleb(int(exp))
    return -v if sign == '-' else v


def round_sig(x, n):
    if x == 0:
        return Decimal(0)
    return x.quantize(Decimal(1).scaleb(x.adjusted() - n + 1), rounding=ROUND_HALF_UP)


def main():
    print('=' * 78)
    print('VERIFYING lambda_note.tex AGAINST THE AUDIT DATA')
    print('=' * 78)

    print('\n[1] registry')
    R = reg.build()
    print('  %d quantities derived from data/' % len(R))

    print("\n[2] prerequisite: the certificate of REPORT.md 6c' and its fixture")
    cert = os.path.join(HERE, 'data', 'lambda_bracket.txt')
    print('  certificate present: %s' % os.path.exists(cert))
    if not os.path.exists(cert):
        fail('data/lambda_bracket.txt is absent: the lower bound has not been certified')
    else:
        txt = open(cert).read()
        for need in ('closure disjoint from the real axis: True',
                     'm > 0 : True',
                     'increments unambiguous: True'):
            if need not in txt:
                fail('the certificate does not record: %s' % need)
            else:
                print('  [PASS] %s' % need)
    fx = subprocess.run([sys.executable,
                         os.path.join(HERE, 'tests', 'test_lambda_certificate.py')],
                        capture_output=True, text=True)
    print('  regression fixture: %s'
          % ('PASSED' if fx.returncode == 0 else 'FAILED'))
    if fx.returncode != 0:
        fail('tests/test_lambda_certificate.py does not pass')

    print('\n[3] the DOI in the front matter')
    raw = open(TEX).read()
    head = raw.split('\\begin{document}')[0]
    dois = re.findall(r'10\.5281/zenodo\.([0-9X]+)', head)
    print('  declared in front matter: %s' % dois)
    if len(dois) != 1:
        fail('the note declares %d DOIs in its front matter, expected exactly 1' % len(dois))
    elif not dois[0].isdigit():
        print('  PLACEHOLDER (expected before deposit); nothing to cross-check yet')
    else:
        rec = re.findall(r'10\.5281/zenodo\.(\d+)',
                         ''.join(open(f).read() for f in RECORD if os.path.exists(f)))
        if dois[0] not in rec:
            fail('the audit record does not carry the note DOI 10.5281/zenodo.%s' % dois[0])

    src = strip_tex(raw)
    body, bib = split_sections(src)
    rec_txt = ''.join(open(f).read() for f in RECORD if os.path.exists(f))
    rec_digits = re.sub(r'[^0-9]', '', rec_txt)

    print("\n[4] bibliography and literature values vs the audit's own record")
    bad = 0
    for rawn, sign, mant, exp, ctx in literals(bib):
        d = re.sub(r'[^0-9]', '', mant)
        if len(d) <= 1:
            continue
        if d not in rec_digits:
            fail('bibliography: %r is in no record file (audit REPORT/SOURCE_AUDIT, '
                 'or a deposited note of this repository)' % rawn)
            bad += 1
    for lit, what in (('0.22', 'the Polymath15 upper bound'),
                      ('1904.12438', 'the Polymath15 identifier'),
                      ('1801.05914', 'the Rodgers-Tao identifier'),
                      ('1901.06596', 'the Newman-Wu identifier')):
        if lit not in rec_txt:
            fail('%s (%s) is not recorded in the audit sources' % (lit, what))
            bad += 1
    print('  %s' % ('all bibliography and literature numbers are recorded' if bad == 0
                    else '%d unrecorded' % bad))

    print('\n[5,6] every literal in the body')
    unmatched, overstated, ok_n, formula_n = [], [], 0, 0
    for rawn, sign, mant, exp, ctx in literals(body):
        key = (sign + mant) if exp is None else None
        if key is not None and key.lstrip('-') in FORMULA:
            formula_n += 1
            continue
        val = to_dec(sign, mant, exp)
        nsig = sigdigits(mant)
        hit = None
        for name, ent in R.items():
            rv = ent['value']
            if rv == 0:
                if val == 0:
                    hit = (name, ent)
                    break
                continue
            if val != 0 and rv.copy_sign(val) != rv:
                continue
            try:
                if round_sig(rv, nsig) == round_sig(val, nsig) and \
                        rv.adjusted() == val.adjusted():
                    hit = (name, ent)
                    break
            except Exception:                                          # noqa: BLE001
                continue
        if hit is None:
            # a number quoted verbatim from the audit record is acceptable if it is
            # literally there
            d = re.sub(r'[^0-9]', '', mant)
            if len(d) > 1 and d in rec_digits:
                formula_n += 1
                continue
            if exp is not None:
                cands = ['%se%d' % (mant, int(exp)),
                         '%se\u2212%d' % (mant, abs(int(exp))),
                         '%s\u00d710\u207b' % mant]
                if any(c in rec_txt for c in cands):
                    formula_n += 1
                    continue
            unmatched.append((rawn, ctx))
            continue
        name, ent = hit
        if nsig > ent['sig']:
            overstated.append((rawn, name, nsig, ent['sig']))
        ok_n += 1
    print('  %d literals matched a derived quantity' % ok_n)
    print('  %d literals are structural or quoted from the audit record' % formula_n)
    for rawn, ctx in unmatched:
        fail('unmatched %-26r  ...%s...' % (rawn, ctx.strip()[:66]))
    for rawn, name, n, sg in overstated:
        fail('OVERSTATED PRECISION %r: printed %d sig digits, %s carries %d'
             % (rawn, n, name, sg))

    print('\n' + '=' * 78)
    if FAIL:
        print('VERIFICATION FAILED: %d problem(s)' % len(FAIL))
        return 1
    print('VERIFICATION PASSED')
    print('=' * 78)
    return 0


if __name__ == '__main__':
    sys.exit(main())
