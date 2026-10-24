#!/usr/bin/env python3
r"""Phase H -- splitting the arithmetic side by the support of b_n.

The arithmetic side of the window form (DERIVATION.md section 3) is

    Q_f(g) = A(g) - 2 sum_{2<=n<=e^{2L}} b_n h(log n)/sqrt(n),
    A(g)   = h(0) log(5/pi) + (1/2pi) int |G(it)|^2 Re psi(3/4+it/2) dt,

and the sum splits by whether n is a prime power:

    P_pp(g)  = -2 sum_{n = p^k}       b_n h(log n)/sqrt(n)
    P_mix(g) = -2 sum_{n not a p^k}   b_n h(log n)/sqrt(n)
    Q_f = A + P_pp + P_mix.

WHAT THE SPLIT IS.  It is canonical, not arbitrary: setting
E(s) = prod_p exp(sum_k b_{p^k} p^{-ks}) and M(s) = exp(sum_{n not p^k} b_n n^{-s}),
both convergent for Re s > sigma_0, one has f = E * M with E carrying an Euler product
BY CONSTRUCTION and M carrying everything that obstructs one.  Dropping P_mix therefore
means "keep f's gamma factor, conductor and Euler-product factor, discard the rest".
Check FACT below verifies f = E*M numerically.

WHAT THE SPLIT IS NOT.  A + P_pp is NOT the Weil form of any L-function and has no
spectral side: E has an Euler product but no functional equation, so the contour
argument that produces the explicit formula does not apply to the hybrid object
(archimedean data of f) + (arithmetic data of E).  No zero of anything corresponds to it.
This is a diagnostic modification of a sum, nothing more.

TWO QUESTIONS, KEPT APART (as in Phase F).
  (A) FIXED WITNESS.  Does the certified witness stay negative with P_mix removed?
      This is additive and exact.
  (B) RE-OPTIMISED.  Is inf_g (A + P_pp)(g)/||g||^2 negative on the same space?
      (A) alone would be misleading: g_* was chosen to minimise the FULL form, so
      deleting a large negative part of the quantity it minimised is expected to leave
      something positive.  Only (B) asks whether the truncated sum can be made negative.

Output: results/phaseH.txt, data/arith_split.csv.
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp                                                      # noqa: E402
import dh                                                                # noqa: E402
import config as cfg                                                     # noqa: E402
from phaseD import index_set                                             # noqa: E402
from zeta_window import Csym_terms, eval_terms, norm2                    # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES, DATA, CERT = (os.path.join(HERE, x) for x in ('results', 'data', 'certificates'))


def is_pp(n):
    m, p = n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            return m == 1
        p += 1
    return True


def h_at(v, idx, L, sector, u):
    tot = mp.mpf(0)
    for i, ki in enumerate(idx):
        for j, kj in enumerate(idx):
            tm = Csym_terms(ki, kj, L, sector)
            den = mp.sqrt(norm2(ki, L, sector) * norm2(kj, L, sector))
            tot += v[i] * v[j] * eval_terms(tm, L, u) / den
    return tot


def split_matrices(L, idx, dps, sector):
    """(A, M_pp, M_mix) with M = A + M_pp + M_mix the full window matrix."""
    Arch, LogPi, LogQ, Prime, terms = dh.build_parts(L, idx, dps, sector)
    n = len(idx)
    A = mp.zeros(n, n)
    Mpp = mp.zeros(n, n)
    Mmix = mp.zeros(n, n)
    nrm = [mp.sqrt(norm2(k, L, sector)) for k in idx]
    for i in range(n):
        for j in range(i, n):
            A[i, j] = A[j, i] = Arch[i, j] - LogPi[i, j] + LogQ[i, j]
            tm = Csym_terms(idx[i], idx[j], L, sector)
            den = nrm[i] * nrm[j]
            pp = mp.mpf(0)
            mx = mp.mpf(0)
            for nn, bn, logn in terms:
                c = -2 * bn / mp.sqrt(mp.mpf(nn)) * eval_terms(tm, L, logn) / den
                if is_pp(nn):
                    pp += c
                else:
                    mx += c
            Mpp[i, j] = Mpp[j, i] = pp
            Mmix[i, j] = Mmix[j, i] = mx
    return A, Mpp, Mmix, terms


def quad(X, v, n2):
    n = len(v)
    return mp.fsum([v[i] * v[j] * X[i, j] for i in range(n) for j in range(n)]) / n2


def main():
    out, rows = [], []

    def p(line=''):
        print(line)
        out.append(line)

    mp.mp.dps = cfg.DPS_BASE
    p('=' * 92)
    p('PHASE H  --  prime-power versus mixed-prime split of the arithmetic side')
    p('=' * 92)

    # ------------------------------------------------------------------ FACT
    p('')
    p('[FACT] f = E * M with E = prod_p exp(sum_k (b_{p^k}/(k log p)) p^{-ks})')
    p('     From -f\'/f = sum b_n n^{-s} and log f = sum c_n n^{-s} one gets')
    p('     c_n = b_n / log n, so  log f(s) = sum_{n>=2} (b_n / log n) n^{-s}.  Splitting')
    p('     that series by the support of n gives log E and log M; log E is supported on')
    p('     prime powers, hence E = prod_p (local factor) is an Euler product BY')
    p('     CONSTRUCTION, and f = E * M.')
    p('')
    p('     %-6s %-26s %-26s %-12s' % ('s', 'log E + log M', 'log f', 'difference'))
    nmax = 20000
    b = dh.b_coeffs(nmax)
    worst = mp.mpf(0)
    for ss in ('3.0', '3.5', '4.0'):
        s = dh.mpf_hi(ss)
        lE = mp.fsum([b[n] / mp.log(n) * mp.power(mp.mpf(n), -s)
                      for n in range(2, nmax + 1) if is_pp(n) and b[n] != 0])
        lM = mp.fsum([b[n] / mp.log(n) * mp.power(mp.mpf(n), -s)
                      for n in range(2, nmax + 1) if not is_pp(n) and b[n] != 0])
        lhs = mp.log(dh.F(s))
        d = abs(lE + lM - lhs)
        worst = max(worst, d)
        p('     %-6s %-26s %-26s %-12s'
          % (ss, mp.nstr(lE + lM, 20), mp.nstr(lhs, 20), mp.nstr(d, 4)))
    p('     worst difference %s against a truncation tail of order n_max^{1-s} ~ %s'
      % (mp.nstr(worst, 4),
         mp.nstr(mp.power(mp.mpf(nmax), -2) / 2, 4)))
    p('     log E alone at s = 4: %s   log M alone: %s'
      % (mp.nstr(mp.fsum([b[n] / mp.log(n) * mp.power(mp.mpf(n), -dh.mpf_hi('4.0'))
                          for n in range(2, nmax + 1) if is_pp(n) and b[n] != 0]), 16),
         mp.nstr(mp.fsum([b[n] / mp.log(n) * mp.power(mp.mpf(n), -dh.mpf_hi('4.0'))
                          for n in range(2, nmax + 1) if not is_pp(n) and b[n] != 0]), 16)))

    # ------------------------------------------------------------------ (A)
    p('')
    p('(A) FIXED WITNESS -- additive, exact')
    for wf, tag in (('witness.json', 'g_*  (L=2.0)'),
                    ('witness_smallest_L.json', 'g_#  (L=1.8)')):
        c = json.load(open(os.path.join(CERT, wf)))
        L, sector, idx = dh.mpf_hi(c['L']), c['sector'], c['idx']
        v = [mp.mpf(s) for s in c['coefficients']]
        n2 = mp.fsum([x * x for x in v])
        A, Mpp, Mmix, terms = split_matrices(L, idx, cfg.DPS_BASE, sector)
        qa, qpp, qmx = quad(A, v, n2), quad(Mpp, v, n2), quad(Mmix, v, n2)
        full = dh.build_matrix(L, idx, cfg.DPS_BASE, sector)
        qf = quad(full, v, n2)
        npp = sum(1 for t in terms if is_pp(t[0]))
        p('')
        p('   %s   N = %d, %d arithmetic terms (%d prime-power, %d mixed)'
          % (tag, len(idx), len(terms), npp, len(terms) - npp))
        p('     A      = arch + log(5/pi) h(0)   = %s' % mp.nstr(qa, 20))
        p('     P_pp   = prime-power terms       = %s' % mp.nstr(qpp, 20))
        p('     P_mix  = mixed-prime terms       = %s' % mp.nstr(qmx, 20))
        p('     Q_f    = A + P_pp + P_mix        = %s' % mp.nstr(qa + qpp + qmx, 20))
        p('     reconciliation |Q_f - (A+P_pp+P_mix)| = %s'
          % mp.nstr(abs(qf - (qa + qpp + qmx)), 4))
        p('     A + P_pp   (mixed dropped)       = %-24s  %s'
          % (mp.nstr(qa + qpp, 20), 'NEGATIVE' if qa + qpp < 0 else 'POSITIVE'))
        p('     A + P_mix  (prime powers dropped)= %-24s  %s'
          % (mp.nstr(qa + qmx, 20), 'NEGATIVE' if qa + qmx < 0 else 'POSITIVE'))
        rows.append(dict(kind='fixed_witness', witness=tag.split()[0], L=c['L'],
                         sector=sector, N=len(idx), n_pp=npp,
                         n_mix=len(terms) - npp, A=mp.nstr(qa, 30),
                         P_pp=mp.nstr(qpp, 30), P_mix=mp.nstr(qmx, 30),
                         Q_f=mp.nstr(qa + qpp + qmx, 30),
                         A_plus_Ppp=mp.nstr(qa + qpp, 30),
                         A_plus_Pmix=mp.nstr(qa + qmx, 30)))
        # per-n contributions
        per = []
        for nn, bn, logn in terms:
            cval = -2 * bn / mp.sqrt(mp.mpf(nn)) * h_at(v, idx, L, sector, logn) / n2
            per.append((nn, is_pp(nn), cval))
        per.sort(key=lambda r: -abs(r[2]))
        p('     largest ten |contributions|:')
        for nn, pp, cv in per[:10]:
            p('        n = %-4d %-11s %s' % (nn, 'prime power' if pp else 'MIXED',
                                             mp.nstr(cv, 14)))
        for nn, pp, cv in per:
            rows.append(dict(kind='per_n', witness=tag.split()[0], L=c['L'],
                             sector=sector, N=nn, n_pp=int(pp), n_mix=0,
                             A='', P_pp='', P_mix='', Q_f=mp.nstr(cv, 30),
                             A_plus_Ppp='', A_plus_Pmix=''))

    # ------------------------------------------------------------------ (B)
    p('')
    p('(B) RE-OPTIMISED -- lambda_min of the truncated form on the SAME space')
    p('    (this, not (A), is the question "can the truncated sum be made negative")')
    p('')
    p('    %-5s %-7s %-4s %-22s %-22s %-22s' %
      ('L', 'sector', 'N', 'lmin(Q_f)', 'lmin(A + P_pp)', 'lmin(A + P_mix)'))
    for Ls in cfg.WINDOW_SCAN:
        L = dh.mpf_hi(Ls)
        for sector in (cfg.SECTOR_EVEN, cfg.SECTOR_ODD):
            idx = index_set(L, sector, 1)
            n = len(idx)
            A, Mpp, Mmix, _ = split_matrices(L, idx, cfg.DPS_BASE, sector)
            def asm(*parts):
                X = mp.zeros(n, n)
                for P in parts:
                    for i in range(n):
                        for j in range(n):
                            X[i, j] += P[i, j]
                return X
            lf = min(mp.eigsy(asm(A, Mpp, Mmix), eigvals_only=True))
            lpp = min(mp.eigsy(asm(A, Mpp), eigvals_only=True))
            lmx = min(mp.eigsy(asm(A, Mmix), eigvals_only=True))
            p('    %-5s %-7s %-4d %-22s %-22s %-22s'
              % (Ls, sector, n, mp.nstr(lf, 16), mp.nstr(lpp, 16), mp.nstr(lmx, 16)))
            rows.append(dict(kind='reoptimised', witness='', L=Ls, sector=sector,
                             N=n, n_pp=0, n_mix=0, A='', P_pp='', P_mix='',
                             Q_f=mp.nstr(lf, 30), A_plus_Ppp=mp.nstr(lpp, 30),
                             A_plus_Pmix=mp.nstr(lmx, 30)))

    p('')
    p('  NOTE.  A + P_pp is not the Weil form of any L-function: E has an Euler product')
    p('  but no functional equation, so the hybrid (archimedean data of f) + (arithmetic')
    p('  data of E) has no spectral side and no zeros correspond to it.  A sign here is')
    p('  a statement about a truncated sum, not about any zero of anything.')
    p('=' * 92)

    with open(os.path.join(DATA, 'arith_split.csv'), 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    with open(os.path.join(RES, 'phaseH.txt'), 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
