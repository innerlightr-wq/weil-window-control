#!/usr/bin/env python3
r"""Ball-arithmetic re-implementation of the Davenport-Heilbronn window form.

This module recomputes Q_F(g) for a FROZEN test function using python-flint / Arb, so
every quantity carries a rigorous enclosure.  It shares no code with src/dh.py: the
cross-correlation terms, the archimedean closed form, kappa, and the b_n are all written
again here against the Arb API.  Agreement between the two implementations is therefore
also an independent check of the mpmath assembly.

Arb supplies rigorous enclosures for exactly the special functions the closed form needs:
digamma, the Hurwitz zeta zeta(s,a), and the Lerch transcendent Phi(w,s,z).  Everything
else is elementary.

kappa is taken in the CLOSED FORM  [sqrt(10 - 2 sqrt 5) - 2]/(sqrt 5 - 1)  -- pure square
roots, hence exactly enclosable -- which Phase A [A1] proved equal to tan(arg eps_chi / 2).

What is and is not certified:
  * CERTIFIED: the value of Q_F(g) for the frozen g, as a ball, given only (i) that Arb's
    enclosures are correct and (ii) the derivation in DERIVATION.md section 3.
  * CERTIFIED: the orbit contribution Q_{O_1}(g), as a ball, given the enclosure of the
    zero rho_1 (which is supplied as an INTERVAL from the Newton/argument-principle
    localisation, and whose correctness as an enclosure of a true zero is T2, not T2-C).
  * NOT certified: that the located rho_1 is a zero of F to within the stated interval --
    that rests on the census, which is T2.  The sign of Q_F(g) does NOT depend on it.
"""
import json
import os
import sys

from flint import acb, arb, ctx                                          # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

EVEN_D, ODD = 'even_d', 'odd'


# ---------------------------------------------------------------------------- basics

def PI():
    return arb.pi()


def freq_a(k, L, sector):
    if sector == ODD:
        return (k + 1) * PI() / L
    return (2 * k + 1) * PI() / (2 * L)


def norm2_a(k, L, sector):
    return L


def kappa_a():
    """[sqrt(10 - 2 sqrt 5) - 2] / (sqrt 5 - 1)."""
    r5 = arb(5).sqrt()
    return ((10 - 2 * r5).sqrt() - 2) / (r5 - 1)


def a_coeff_a(n):
    k = kappa_a()
    r = int(n) % 5
    return (arb(1) if r == 1 else arb(-1) if r == 4 else
            k if r == 2 else -k if r == 3 else arb(0))


def b_coeffs_a(nmax):
    """b_n by the recursion of DERIVATION.md (2.1), in ball arithmetic."""
    nmax = int(nmax)
    a = [arb(0)] * (nmax + 1)
    for n in range(1, nmax + 1):
        a[n] = a_coeff_a(n)
    b = [arb(0)] * (nmax + 1)
    for n in range(2, nmax + 1):
        tot = a[n] * arb(n).log()
        d = 2
        while d * d <= n:
            if n % d == 0:
                tot -= b[d] * a[n // d]
                e = n // d
                if e != d and e < n:
                    tot -= b[e] * a[d]
            d += 1
        b[n] = tot
    return b


def b_terms_a(twoL_arb, nmax_int):
    """[(n, b_n, log n)] for 2 <= n <= nmax_int with log n < 2L.

    nmax_int is supplied as an INTEGER by the caller (computed once, in integer
    arithmetic, from floor(e^{2L})), so no float comparison decides the cutoff; the
    enclosure log n < 2L is then re-tested rigorously and a term whose enclosure
    straddles the cutoff would raise."""
    b = b_coeffs_a(nmax_int)
    out = []
    for n in range(2, nmax_int + 1):
        ln = arb(n).log()
        if (ln < twoL_arb) is True:
            if b[n] != 0:
                out.append((n, b[n], ln))
        elif (ln > twoL_arb) is not True:
            raise RuntimeError('cutoff straddles n = %d' % n)
    return out


# ------------------------------------------------------------- cross-correlation terms

def C_terms_a(j, k, L, sector):
    """Port of zeta_window.C_terms: C_{jk}(u) on [0, 2L] as a list of
       ('sin', c, alpha, beta) meaning c sin(alpha u + beta)
       ('lin', c, gamma)        meaning c (2L - u) cos(gamma u)."""
    wj, wk = freq_a(j, L, sector), freq_a(k, L, sector)
    p1, p2 = wj - wk, wj + wk
    s2 = -1 if sector == ODD else 1
    out = []
    if j == k:
        out.append(('lin', arb(1) / 2, wk))
    else:
        out.append(('sin', 1 / (2 * p1), wk, p1 * L))
        out.append(('sin', -1 / (2 * p1), wj, -p1 * L))
    out.append(('sin', s2 / (2 * p2), -wk, p2 * L))
    out.append(('sin', -s2 / (2 * p2), wj, -p2 * L))
    return out


def Csym_terms_a(j, k, L, sector):
    if j == k:
        return C_terms_a(j, k, L, sector)
    return ([(t[0], t[1] / 2) + t[2:] for t in C_terms_a(j, k, L, sector)]
            + [(t[0], t[1] / 2) + t[2:] for t in C_terms_a(k, j, L, sector)])


def eval_terms_a(terms, L, u):
    u = abs(u)
    tot = arb(0)
    for t in terms:
        if t[0] == 'lin':
            _, c, g = t
            tot += c * (2 * L - u) * (g * u).cos()
        else:
            _, c, al, be = t
            tot += c * (al * u + be).sin()
    return tot


# ------------------------------------------------------------------- archimedean term

def arch_closed_a_ball(terms, L, a):
    r"""Port of dirichlet_window.arch_closed_a at parity a, in ball arithmetic.

        Arch = -gamma_E C(0)
               + sum_{m>=0} [ C(0)/(m+1) - 2 int_0^{2L} C(u) e^{-c_m u} du ],
        c_m = 2m + 1/2 + a,

    summed in closed form via digamma / Lerch Phi / Hurwitz zeta.  The identity being
    evaluated is the parent repository's; only the arithmetic is different."""
    U = 2 * L
    w = (-2 * U).exp()
    zr = arb(1) / 4 + arb(a) / 2
    pref = (-(arb(1) / 2 + a) * U).exp()
    C0 = eval_terms_a(terms, L, arb(0))
    tot = arb(0)
    for t in terms:
        if t[0] == 'lin':
            _, c, beta = t
            z0 = acb(zr, -beta / 2)
            dg = (z0.digamma() + arb.const_euler()).real
            lp2 = acb(w).lerch_phi(acb(2), z0)
            hz = acb(2).zeta(z0)
            tot += c * U * dg
            tot -= (c / 2) * (acb(0, beta * U).exp() * pref * lp2 - hz).real
        else:
            _, c, alpha, beta = t
            z = acb(zr, -alpha / 2)
            dgg = z.digamma() + arb.const_euler()
            lp1 = acb(w).lerch_phi(acb(1), z)
            A = acb(0, beta + alpha * U).exp()
            tot += c * (acb(0, beta).exp() * dgg).imag
            tot += c * pref * (A * lp1).imag
    return -arb.const_euler() * C0 + tot


# ------------------------------------------------------------------------ the form

def Q_of(coeffs, idx, L, sector, nmax_int, drop_n=()):
    """Q_F(g) for g = sum_m coeffs[m] f_{idx[m]} / sqrt(norm2), as a ball."""
    terms = [t for t in b_terms_a(2 * L, nmax_int) if t[0] not in set(drop_n)]
    logqpi = (arb(5) / PI()).log()
    tot = arb(0)
    nn = len(idx)
    for i in range(nn):
        for j in range(nn):
            ci, cj = coeffs[i], coeffs[j]
            if ci == 0 or cj == 0:
                continue
            tm = Csym_terms_a(idx[i], idx[j], L, sector)
            den = (norm2_a(idx[i], L, sector) * norm2_a(idx[j], L, sector)).sqrt()
            v = arch_closed_a_ball(tm, L, 1) / den
            if i == j:
                v += logqpi * eval_terms_a(tm, L, arb(0)) / den
            pr = arb(0)
            for n, bn, logn in terms:
                pr += 2 * bn / arb(n).sqrt() * eval_terms_a(tm, L, logn)
            v -= pr / den
            tot += ci * cj * v
    return tot


def G_ball(coeffs, idx, L, sector, z):
    """G(z) = int g(u) e^{zu} du as an acb ball."""
    tot = acb(0)
    for m, k in enumerate(idx):
        c = coeffs[m]
        if c == 0:
            continue
        w = freq_a(k, L, sector)
        den = z * z + w * w
        if sector == ODD:
            num = 2 * (z * (w * L).sin() * (z * L).cosh()
                       - w * (w * L).cos() * (z * L).sinh())
        else:
            num = 2 * (z * (w * L).cos() * (z * L).sinh()
                       + w * (w * L).sin() * (z * L).cosh())
        tot += c * num / den / (norm2_a(k, L, sector)).sqrt()
    return tot


def quartet_ball(coeffs, idx, L, sector, beta, gamma, mult=1):
    """4 Re(A conj B); beta and gamma are acb/arb BALLS enclosing the zero."""
    d = beta - arb(1) / 2
    A = G_ball(coeffs, idx, L, sector, acb(d, gamma))
    B = G_ball(coeffs, idx, L, sector, acb(-d, gamma))
    return mult * 4 * (A * B.conjugate()).real


def online_ball(coeffs, idx, L, sector, gamma, mult=1):
    z = acb(0, gamma)
    return mult * 2 * (G_ball(coeffs, idx, L, sector, z)
                       * G_ball(coeffs, idx, L, sector, -z)).real


def load_witness(path=None):
    path = path or os.path.join(HERE, 'certificates', 'witness.json')
    with open(path) as fh:
        c = json.load(fh)
    L = arb(c['L'])
    coeffs = [arb(s) for s in c['coefficients']]
    return c, L, coeffs


def nmax_for(Lstr):
    """floor(e^{2L}) + 2, computed in integer arithmetic from a high-precision e^{2L}."""
    old = ctx.prec
    ctx.prec = 400
    try:
        v = (2 * arb(Lstr)).exp()
        return int(float(v.mid())) + 2
    finally:
        ctx.prec = old


def fmt(x):
    return str(x)


def main():
    # usage: certify.py [prec_bits] [witness.json] [out_name]
    ctx.prec = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    wpath = sys.argv[2] if len(sys.argv) > 2 else None
    oname = sys.argv[3] if len(sys.argv) > 3 else 'ball_certificate.txt'
    out = []

    def p(line=''):
        print(line)
        out.append(line)

    c, L, coeffs = load_witness(wpath)
    idx = c['idx']
    sector = c['sector']
    nmax = nmax_for(c['L'])
    p('=' * 78)
    p('BALL-ARITHMETIC CERTIFICATE FOR THE FROZEN WITNESS')
    p('=' * 78)
    p('  implementation : python-flint / Arb, working precision %d bits' % ctx.prec)
    p('  every number below is a BALL [midpoint +/- radius] that rigorously encloses the')
    p('  quantity named, given only Arb\'s enclosures and the derivation in DERIVATION.md.')
    p('')
    p('  THE TEST FUNCTION')
    p('    L          = %s   (window half-length of g; h = g * g lives on [-2L, 2L])' % c['L'])
    p('    sector     = %s' % sector)
    p('    basis      = %s, k in idx, normalised by 1/sqrt(L)' % c['basis'])
    p('    N          = %d   basis level %d' % (len(idx), c['level']))
    p('    idx        = %s' % c['idx'])
    p('    arithmetic cutoff: n <= floor(e^{2L}) + 2 = %d, integer-decided' % nmax)
    terms = b_terms_a(2 * L, nmax)
    p('    arithmetic terms actually used (log n < 2L, b_n != 0): %d' % len(terms))
    p('    coefficients (exact finite decimals, 30 digits, read back as strings):')
    for m, k in enumerate(idx):
        p('      c[%-2d] (k=%-3d, w = %-24s) = %s'
          % (m, k, freq_a(k, L, sector).str(18), c['coefficients'][m]))
    p('')
    n2 = sum((x * x for x in coeffs), arb(0))
    p('  ||g||^2 = sum_m c_m^2 = %s' % n2.str(40))
    p('')
    p('  THE VALUE')
    Q = Q_of(coeffs, idx, L, sector, nmax)
    R = Q / n2
    p('    Q_F(g_*)          = %s' % Q)
    p('    Q_F(g_*)/||g||^2  = %s' % R)
    p('    radius            = %s' % R.rad())
    p('    STRICTLY NEGATIVE : %s' % ((R < 0) is True))
    p('')
    p('  THE DECOMPOSITION  (O_1 is the selected off-line quartet)')
    beta = arb(c['beta1'])
    gamma = arb(c['gamma1'])
    p('    beta_1  = %s' % c['beta1'])
    p('    gamma_1 = %s' % c['gamma1'])
    QO = quartet_ball(coeffs, idx, L, sector, beta, gamma)
    QR = Q - QO
    p('    Q_{O_1}(g_*)      = %s' % QO)
    p('    Q_rest(g_*)       = %s' % QR)
    p('    Q_{O_1}/||g||^2   = %s' % (QO / n2))
    p('    Q_rest /||g||^2   = %s' % (QR / n2))
    p('    additive reconciliation:  Q_F - (Q_{O_1} + Q_rest) = %s' % (Q - (QO + QR)))
    p('    sign of Q_rest    : %s'
      % ('POSITIVE' if (QR > 0) is True else
         'NEGATIVE' if (QR < 0) is True else 'UNRESOLVED'))
    p('')
    p('  CAVEAT ON THE ORBIT TERM.  Q_{O_1} is certified as a FUNCTION of the supplied')
    p('    (beta_1, gamma_1).  That those numbers enclose an actual zero of F rests on the')
    p('    census, which is T2, not T2-C.  The sign of Q_F(g_*) does not depend on it:')
    p('    Q_F is computed from the arithmetic side alone and knows nothing about where')
    p('    any zero is.')
    p('')
    p('  ERROR BUDGET for Q_F(g_*)/||g||^2')
    p('    Everything below is inside the single Arb radius %s.' % R.rad())
    p('      rounding of every arithmetic operation      : Arb ball arithmetic')
    p('      kappa                                       : sqrt only, [sqrt(10-2 sqrt5)-2]/(sqrt5-1)')
    p('      b_n (n <= %-3d)                               : exact recursion over kappa and log n' % nmax)
    p('      log(5/pi), log n, sqrt n                     : Arb elementary functions')
    p('      archimedean closed form                     : Arb digamma, Hurwitz zeta, Lerch Phi')
    p('      the cross-correlations C_jk(u)              : Arb sin/cos of exact arguments')
    p('      truncation of the arithmetic sum            : NONE -- h vanishes beyond 2L,')
    p('                                                    the sum is finite and complete')
    p('      spectral tail                               : NONE -- no zero is used')
    p('      quadrature                                  : NONE -- no quadrature is used')
    p('      root uncertainty                            : NONE for Q_F; enters Q_{O_1} only')
    p('    The coefficients are exact finite decimals, so the test function carries no')
    p('    error of its own: the ball is an enclosure for THIS g, not for an approximation')
    p('    to some other g.')
    p('')
    p('  TIER: T2-C.  This is an Arb ball enclosure and a reproducible certificate.')
    p('        It is NOT a Lean verification and is nowhere described as one.')
    p('=' * 78)
    with open(os.path.join(HERE, 'certificates', oname), 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
