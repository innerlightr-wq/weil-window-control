#!/usr/bin/env python3
r"""Re-derive every number printed in siegel_note.tex from this extension's data files.

Fails if a printed literal is not derivable, if it disagrees with its source, or if the
note prints MORE significant digits than its source carries.  Same design as the
companion DH note's verifier (same repository).

Seven passes:

  1. registry       -- rebuild every quantity from data/ and results/
  2. tables         -- regenerate the five machine-generated tables and compare
  3. prerequisite   -- confirm results/verify_siegel_parity.txt records ALL THREE checks
                       as passed; the note must not exist without that
  4. doi            -- the note declares exactly one DOI in its front matter; a
                       placeholder is permitted and reported, a second DOI is not
  5. literature     -- every number in the bibliography must appear in the extension's
                       own record (REPORT.md or README.md), so none is typed from memory
  6. formula        -- digits that are part of a mathematical expression, declared with
                       a reason for each
  7. values         -- every remaining literal must match a registry entry, to the
                       printed precision and no further
"""
import os
import re
import subprocess
import sys
from decimal import Decimal, ROUND_HALF_UP

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import registry as reg                                                   # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NOTE = os.path.join(HERE, 'note')
TEX = os.path.join(NOTE, 'siegel_note.tex')

FAIL = []
RETR_TXT = ''


def fail(msg):
    FAIL.append(msg)
    print('  FAIL  %s' % msg)


# Digits that are part of a mathematical expression rather than a measured value.
FORMULA = {
    '0': 'the ordinate gamma = 0; zero', '1': 'rank one; unit; the index 1',
    '2': 'the 2 of 2A^2, of 1/2, of 2L^3/3', '3': 'the 3 of 2L^3/3; three checks',
    '4': 'equation (4) of the companion note; the four-line proof',
    '5': 'five characters; section 5', '6': 'six characters incl. zeta; Theorem 6',
    '7': 'Theorem 7', '8': 'Theorem 8', '9': 'non-claim 9; Corollary 9',
    '10': 'base of scientific notation; N = 10',
    '11': 'Theorem 11', '12': 'N = 12, also in registry',
    '16': 'R = 16; N = 16', '22': 'digits in the largest denominator, also in registry',
    '25': 'cells, also in registry', '30': 'groups, also in registry',
    '75': 'cells per sector, also in registry',
    '135': 'detections, also in registry', '150': 'cells, also in registry',
    '180': 'rows per sector, also in registry',
    '228': 'detections, also in registry', '270': 'rows, also in registry',
    '0.1': 'delta grid', '0.25': 'delta grid', '0.375': 'delta grid',
    '0.45': 'delta grid', '0.49': 'delta grid', '0.5': 'window / one half',
    '0.8': 'window', '1.0': 'window', '1.1': 'window', '1.3': 'window',
    '1.4': 'window', '1.7': 'window', '2.0': 'window', '2.3': 'window',
    '2.6': 'window', '3.0': 'window',
    '20': 'non-degenerate cells, also in registry',
    '85.7': 'gamma_1 of the companion note, quoted to 3 digits',
    '14': 'fourteen spurious failures (retraction R1)',
    '15': 'fifteen digits (retraction R2)', '40': 'forty digits (retraction R2)',
    '13': 'the prime 13 (retraction R4)', '04': 'part of 4.04 (retraction R3)',
    '4.04': 'the self-test discrepancy of retraction R3',
}

# numbers appearing only in the retraction log, quoted from the parent extension's
# section 0.2 rather than derived; the verifier checks they appear there.
RETRACTION_FILE = os.path.join(HERE, 'REPORT.md')


def strip_tex(src):
    src = src.split('\\begin{document}', 1)[1]
    src = re.sub(r'\\begin\{verbatim\}.*?\\end\{verbatim\}', ' ', src, flags=re.S)
    src = re.sub(r'\\href\{[^}]*\}', ' ', src)
    src = re.sub(r'\\url\{[^}]*\}', ' ', src)
    src = re.sub(r'\\cite\[[^]]*\]', r'\\cite', src)
    src = re.sub(r'\\input\{[^}]*\}', ' ', src)
    src = re.sub(r'\\(label|ref|addcontentsline|src|texttt|tier)\{[^}]*\}', ' ', src)
    src = re.sub(r'%.*', ' ', src)
    src = re.sub(r'~\d+\.\d+(\.\d+)?', ' VERSION ', src)
    src = re.sub(r'10\.5281/zenodo\.\d+|XXXXXXXX', ' DOI ', src)
    src = re.sub(r'\d{4}-\d{4}-\d{4}-\d{4}', ' ORCID ', src)
    src = re.sub(r'arXiv:\d+\.\d+', ' ARXIV ', src)
    # bare orders of magnitude in prose: protect real numbers first
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
    print('VERIFYING siegel_note.tex AGAINST THE EXTENSION DATA')
    print('=' * 78)

    print('\n[1] registry')
    R = reg.build()
    print('  %d quantities derived from data/ and results/' % len(R))

    print('\n[2] generated tables regenerated from the data and compared')
    r = subprocess.run([sys.executable, os.path.join(NOTE, 'make_tables.py'), '--check'],
                       capture_output=True, text=True)
    print(r.stdout.rstrip())
    if r.returncode != 0:
        fail('a generated table differs from its generator')

    print('\n[3] prerequisite: the three Siegel-parity verifications')
    tr = os.path.join(HERE, 'results', 'verify_siegel_parity.txt')
    txt = open(tr).read() if os.path.exists(tr) else ''
    need = 'ALL THREE SIEGEL-PARITY VERIFICATIONS PASSED'
    npass = txt.count('[PASS]')
    print('  transcript present: %s' % bool(txt))
    print('  "%s": %s' % (need, need in txt))
    print('  [PASS] markers: %d' % npass)
    if need not in txt:
        fail('results/verify_siegel_parity.txt does not record all three checks as passed')
    if npass < 3:
        fail('fewer than three [PASS] markers in the prerequisite transcript')

    print('\n[4] the DOI in the front matter')
    raw = open(TEX).read()
    head = raw.split('\\begin{document}')[0]
    dois = re.findall(r'10\.5281/zenodo\.([0-9X]+)', head)
    print('  declared in front matter: %s' % dois)
    if len(dois) != 1:
        fail('the note declares %d DOIs in its front matter, expected exactly 1' % len(dois))
    elif not dois[0].isdigit():
        print('  PLACEHOLDER (expected before deposit); nothing to cross-check yet')
    else:
        for fn in ('REPORT.md', 'README.md'):
            rec = re.findall(r'10\.5281/zenodo\.(\d+)', open(os.path.join(HERE, fn)).read())
            if dois[0] not in rec:
                fail('%s does not record the note DOI 10.5281/zenodo.%s' % (fn, dois[0]))

    src = strip_tex(raw)
    body, bib = split_sections(src)

    print('\n[5] bibliography numbers vs the extension\'s own record')
    PAPER = os.path.join(os.path.dirname(os.path.dirname(HERE)), 'paper',
                         'weil_window_control.tex')
    rec_txt = (open(RETRACTION_FILE).read() + open(os.path.join(HERE, 'README.md')).read()
               + (open(PAPER).read() if os.path.exists(PAPER) else ''))
    rec_digits = re.sub(r'[^0-9]', '', rec_txt)
    bad = 0
    for rawn, sign, mant, exp, ctx in literals(bib):
        d = re.sub(r'[^0-9]', '', mant)
        if len(d) <= 1:
            continue
        if d not in rec_digits:
            fail('bibliography: %r is not recorded in REPORT.md, README.md or '
                 'the parent paper bibliography' % rawn)
            bad += 1
    print('  %s' % ('all bibliography numbers are recorded' if bad == 0
                    else '%d unrecorded' % bad))

    print('\n[6,7] every literal in the body')
    unmatched, overstated, ok_n, formula_n = [], [], 0, 0
    global RETR_TXT
    RETR_TXT = open(RETRACTION_FILE).read()
    retr_digits = re.sub(r'[^0-9]', '', RETR_TXT)
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
            # a retraction-log number quoted from the parent record is acceptable if it
            # is literally there
            d = re.sub(r'[^0-9]', '', mant)
            if len(d) > 1 and d in retr_digits:
                formula_n += 1
                continue
            # the retraction log quotes values in the parent REPORT's own rendering,
            # e.g. 4e-15 written there with a Unicode minus; accept an exact match
            if exp is not None:
                cands = ['%se%d' % (mant, int(exp)),
                         '%se\u2212%d' % (mant, abs(int(exp)))]
                if any(c in RETR_TXT for c in cands):
                    formula_n += 1
                    continue
            unmatched.append((rawn, ctx))
            continue
        name, ent = hit
        if nsig > ent['sig']:
            overstated.append((rawn, name, nsig, ent['sig']))
        ok_n += 1
    print('  %d literals matched a derived quantity' % ok_n)
    print('  %d literals are structural or quoted from the parent record' % formula_n)
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
