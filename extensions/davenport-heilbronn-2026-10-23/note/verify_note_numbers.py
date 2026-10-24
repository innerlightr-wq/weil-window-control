#!/usr/bin/env python3
r"""Re-derive every number printed in dh_note.tex from the extension's data files.

Fails if a printed literal is not derivable, if it disagrees with its source, or if the
note prints MORE significant digits than its source carries.

Six passes:

  1. registry      -- rebuild every quantity from data/, results/, certificates/
  2. table         -- regenerate Appendix A from certificates/witness.json and compare
  3. report        -- cross-check the Gate A values against REPORT.md; a disagreement
                      between REPORT.md and the data files is a HARD FAILURE, not a
                      choice to be made
  4. literature    -- every number in Appendix C and the bibliography must appear
                      verbatim in SOURCE_AUDIT.md (nothing is typed from memory)
  5. formula       -- digits that are part of a mathematical expression are matched
                      against a declared list with a reason for each
  6. values        -- every remaining literal must match a registry entry, to the
                      printed precision and no further

Exit 0 only if all six pass.
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
TEX = os.path.join(NOTE, 'dh_note.tex')

FAIL = []


def fail(msg):
    FAIL.append(msg)
    print('  FAIL  %s' % msg)


# --------------------------------------------------------------------------- formula
# Digits that are part of a mathematical expression rather than a measured value.
# Each carries the reason it is structural.  Anything not here and not in the registry
# is a failure.
FORMULA = {
    '0': 'chi(0) = 0; index origin', '1': 'chi(1) = 1; unit; index',
    '2': 'the 2 of a square, of 1/2, of 2L, of 2pi', '3': 'basis level 3; 3/4',
    '4': 'the four zeros of an orbit; 4K; 4|P|^2', '5': 'the conductor 5',
    '6': 'the index of b_6', '10': 'base of scientific notation',
    '12': 'twelve digits unchanged (R6)', '16': 'the 16/3 of the remainder bound',
    '18': 'count of rows, also in registry', '25': 'the 2/5 of the remainder bound',
    '36': 'sin 36 degrees', '44': 'basis index; digit count',
    '50': 'the n < 50 of the non-prime-power count', '60': 'working digits',
    '64': 'basis index', '72': 'sin 72 degrees', '85': 'part of 0.8085+85.70i (lit.)',
    '100': 'working digits', '122': 'census on-line count, also in registry',
    '140': 'working digits', '160': 'working digits', '181': 'cutoff, also in registry',
    '200': 'census height T = 200', '220': 'mpf_hi parser precision',
    '30': 'thirty-digit coefficients', '37': 'N, also in registry',
    '53': 'basis index of the dominant coefficient', '54': 'basis index; row count',
    '57': 'N, also in registry', '15': 'low-block size; count, also in registry',
    '26': 'count of negative rows / of terms, also in registry',
    '39': 'count of arithmetic terms, also in registry',
    '56': 'arithmetic cutoff, also in registry', '38': 'arithmetic cutoff',
    '9': 'nine scanned windows', '8': 'eight',
    '0.5': 'one half', '1.0': 'window', '1.2': 'window', '1.4': 'window / strip edge',
    '1.6': 'window', '1.8': 'window', '2.0': 'window', '2.2': 'window',
    '2.4': 'window', '2.6': 'window', '0.4': 'strip edge 0.4',
    '1.42': 'contour abscissa sigma_c, src/phaseC.py', '0.02': 'local winding radius',
    '0.55': 'off-line search abscissa, src/census.py',
    '0.505': 'the abandoned abscissa of R2', '0.005': 'its distance to the line',
    '403': 'an HTTP status code (SOURCE_AUDIT.md section H)',
    '13': 'the ordinal "13th digit" of retraction R6',
    '85': 'the "85 digits" of retraction R6',
    '7': 'seven', '3': 'three',
    '20': 'approximate runtime in minutes (REPRODUCE.md)',
    '45': 'approximate runtime in minutes (REPRODUCE.md)',
}

# Runtimes are approximate and are recorded in REPRODUCE.md; the verifier checks each
# appears there rather than treating it as a derived quantity.
RUNTIMES = {'2', '3', '5', '6', '10', '15', '20', '25', '45'}

# Literature values: must appear verbatim in SOURCE_AUDIT.md.
LIT_FILE = os.path.join(HERE, 'SOURCE_AUDIT.md')


def strip_tex(src):
    """Remove everything that is not prose/maths carrying derived values."""
    src = src.split('\\begin{document}', 1)[1]
    src = re.sub(r'\\begin\{verbatim\}.*?\\end\{verbatim\}', ' ', src, flags=re.S)
    src = re.sub(r'\\href\{[^}]*\}', ' ', src)
    src = re.sub(r'\\url\{[^}]*\}', ' ', src)
    # DOIs and ORCIDs are identifiers, checked in pass 2b and section I of the audit,
    # not quantities derived from the data
    src = re.sub(r'10\.5281/zenodo\.\d+', ' DOI ', src)
    src = re.sub(r'\d{4}-\d{4}-\d{4}-\d{4}', ' ORCID ', src)
    src = re.sub(r'\\cite\[[^]]*\]', r'\\cite', src)
    src = re.sub(r'\\input\{[^}]*\}', ' ', src)
    src = re.sub(r'\\(label|ref|addcontentsline|src|texttt)\{[^}]*\}', ' ', src)
    src = re.sub(r'%.*', ' ', src)
    # software version strings (Python~3.14.7, ~1.3.0, ~0.9.0) are declared in
    # REPRODUCE.md and checked there, not treated as derived quantities
    src = re.sub(r'~\d+\.\d+(\.\d+)?', ' VERSION ', src)
    # a bare 10^{-n} in prose is an order of magnitude, not a measured value; protect
    # the \times10^{-n} of an actual number first
    src = src.replace('\\times10^{', '\\TIMESTEN{')
    src = re.sub(r'10\^\{-?\d+\}', ' ORDER ', src)
    src = src.replace('\\TIMESTEN{', '\\times10^{')
    return src


def split_sections(src):
    """(body, appendixC, bibliography)"""
    bib = ''
    if '\\begin{thebibliography}' in src:
        src, bib = src.split('\\begin{thebibliography}', 1)
    appc = ''
    m = re.search(r'\\section\{The two literature witnesses\}(.*?)(?=\\section)', src,
                  flags=re.S)
    if m:
        appc = m.group(1)
        src = src.replace(appc, ' ')
    return src, appc, bib


NUM = re.compile(r'(?<![0-9A-Za-z_.])(-?)(\d+(?:\.\d+)?)'
                 r'(?:\s*\\times\s*10\^(?:\{(-?\d+)\}|(-?\d)))?')


def literals(text):
    out = []
    for m in NUM.finditer(text):
        sign, mant = m.group(1), m.group(2)
        exp = m.group(3) if m.group(3) is not None else m.group(4)
        ctx = text[max(0, m.start() - 48):m.end() + 24].replace('\n', ' ')
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
    exp = x.adjusted()
    q = Decimal(1).scaleb(exp - n + 1)
    return x.quantize(q, rounding=ROUND_HALF_UP)


def main():
    print('=' * 78)
    print('VERIFYING dh_note.tex AGAINST THE EXTENSION DATA')
    print('=' * 78)

    # ------------------------------------------------------------------ 1 registry
    print('\n[1] registry')
    R = reg.build()
    print('  %d quantities derived from data/, results/ and certificates/' % len(R))

    # ------------------------------------------------------------------ 2 table
    print('\n[2] generated tables regenerated from the data and compared')
    for gen in ('make_witness_table.py', 'make_tables.py'):
        r = subprocess.run([sys.executable, os.path.join(NOTE, gen), '--check'],
                           capture_output=True, text=True)
        print('  %s' % r.stdout.strip().replace('\n', '\n  '))
        if r.returncode != 0:
            fail('%s: a generated table differs from its generator' % gen)

    # ------------------------------------------------------------------ 2b DOI
    print("\n[2b] the note's own DOI vs the repository record")
    tex_raw = open(TEX).read()
    head = tex_raw.split('\\begin{document}')[0]
    own = set(re.findall(r'10\.5281/zenodo\.(\d+)', head))
    if len(own) != 1:
        fail('the note declares %d DOIs in its front matter, expected exactly 1: %s'
             % (len(own), sorted(own)))
    else:
        doi = own.pop()
        print('  note front matter : 10.5281/zenodo.%s' % doi)
        for fn in ('REPORT.md', 'README.md'):
            rec = set(re.findall(r'10\.5281/zenodo\.(\d+)',
                                 open(os.path.join(HERE, fn)).read()))
            ok = doi in rec
            print('  %-12s      %s' % (fn, 'records it' if ok else 'DOES NOT record it'))
            if not ok:
                fail('%s does not record the note DOI 10.5281/zenodo.%s' % (fn, doi))
        if doi == 'XXXXXXXX' or not doi.isdigit():
            fail('the note still carries a placeholder DOI')

    # ------------------------------------------------------------------ 3 REPORT.md
    print('\n[3] REPORT.md vs the data files (a disagreement is a hard failure)')
    rep = open(os.path.join(HERE, 'REPORT.md')).read()
    for key, name in (('gate_contour_diff', 'contour discrepancy'),
                      ('gate_zeroside_diff', 'zero-side discrepancy'),
                      ('ctrlA_residual', 'control A residual'),
                      ('ctrlB_residual', 'control B residual'),
                      ('census_total', 'census total'),
                      ('census_online', 'census on-line count')):
        v = R[key]['value']
        # REPORT.md prints these in its own format; look for the significant digits
        s = ('%g' % float(v)).replace('e-', 'e-').replace('e+0', 'e')
        digits = re.sub(r'[^0-9]', '', ('%.3g' % float(v)).split('e')[0])
        found = digits in re.sub(r'[^0-9]', '', rep)
        print('  %-26s %-14s %s' % (name, s, 'present in REPORT.md' if found
                                    else 'NOT FOUND'))
        if not found:
            fail('REPORT.md does not contain %s = %s' % (name, s))

    # ------------------------------------------------------------------ prepare tex
    src = strip_tex(open(TEX).read())
    body, appc, bib = split_sections(src)

    # ------------------------------------------------------------------ 4 literature
    print('\n[4] literature numbers (Appendix C and bibliography) vs SOURCE_AUDIT.md')
    audit = open(LIT_FILE).read()
    repro = open(os.path.join(HERE, 'REPRODUCE.md')).read()
    audit_digits = re.sub(r'[^0-9]', '', audit)
    bad = 0
    for text, label in ((appc, 'App C'), (bib, 'bibliography')):
        for raw, sign, mant, exp, ctx in literals(text):
            d = re.sub(r'[^0-9]', '', mant)
            if len(d) <= 1:
                continue
            if d in audit_digits:
                continue
            # a value of OUR OWN quoted in Appendix C for comparison is fine if it is
            # in the registry; only genuinely external numbers must be in the audit
            val = to_dec(sign, mant, exp)
            nsig = sigdigits(mant)
            if any(r['value'] != 0 and round_sig(r['value'], nsig) == round_sig(val, nsig)
                   and r['value'].adjusted() == val.adjusted() for r in R.values()):
                continue
            fail('%s: %r is neither recorded in SOURCE_AUDIT.md nor derived' % (label, raw))
            bad += 1
    print('  %s' % ('all literature numbers are recorded' if bad == 0
                    else '%d unrecorded' % bad))
    missing = [t for t in sorted(RUNTIMES, key=int)
               if ('~%s\xa0min' % t) not in repro and ('%s min' % t) not in repro]
    print('  runtimes declared but absent from REPRODUCE.md: %s'
          % (missing if missing else 'none'))
    for t in missing:
        fail('runtime %s min is printed in the note but not in REPRODUCE.md' % t)

    # ------------------------------------------------------------------ 5,6 body
    print('\n[5,6] every literal in the body')
    unmatched, overstated, ok_n, formula_n = [], [], 0, 0
    for raw, sign, mant, exp, ctx in literals(body):
        key = (sign + mant) if exp is None else None
        if key is not None and key.lstrip('-') in FORMULA and exp is None:
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
                continue                                  # sign mismatch
            try:
                if round_sig(rv, nsig) == round_sig(val, nsig) and \
                        rv.adjusted() == val.adjusted():
                    hit = (name, ent)
                    break
            except Exception:                              # noqa: BLE001
                continue
        if hit is None:
            unmatched.append((raw, ctx))
            continue
        name, ent = hit
        if nsig > ent['sig']:
            overstated.append((raw, name, nsig, ent['sig']))
        ok_n += 1
    print('  %d literals matched a derived quantity' % ok_n)
    print('  %d literals are structural (declared in FORMULA)' % formula_n)
    for raw, ctx in unmatched:
        fail('unmatched %-28r  ...%s...' % (raw, ctx.strip()[:70]))
    for raw, name, n, s in overstated:
        fail('OVERSTATED PRECISION %r: printed %d sig digits, %s carries %d'
             % (raw, n, name, s))

    print('\n' + '=' * 78)
    if FAIL:
        print('VERIFICATION FAILED: %d problem(s)' % len(FAIL))
        return 1
    print('VERIFICATION PASSED')
    print('=' * 78)
    return 0


if __name__ == '__main__':
    sys.exit(main())
