#!/usr/bin/env python3
r"""Phase C unit tests -- orbit conventions and multiplicity bookkeeping.

Each test compares the orbit routines in dh.py against a BRUTE-FORCE sum of
H(rho - 1/2) over an explicitly written-out list of zeros, with H computed by direct
quadrature of g rather than by the closed-form G_at.  Nothing here touches F(s).
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), 'src'))
import mpmath as mp                                                      # noqa: E402
import dh                                                                # noqa: E402
from zeta_window import freq, norm2                                      # noqa: E402

FAIL = []


def check(name, got, want, tol):
    ok = abs(got - want) <= tol
    print('  %-58s %-14s %s' % (name, 'OK' if ok else 'FAIL',
                                mp.nstr(abs(got - want), 4)))
    if not ok:
        FAIL.append((name, mp.nstr(got, 25), mp.nstr(want, 25)))
    return ok


def g_fn(coeffs, idx, L, sector):
    def g(u):
        tot = mp.mpf(0)
        for c, k in zip(coeffs, idx):
            w = freq(k, L, sector)
            tot += c * ((mp.sin(w * u) if sector == 'odd' else mp.cos(w * u))
                        / mp.sqrt(norm2(k, L, sector)))
        return tot
    return g


def H_brute(coeffs, idx, L, sector, z):
    g = g_fn(coeffs, idx, L, sector)
    Gp = mp.quad(lambda u: g(u) * mp.e ** (z * u), [-L, 0, L])
    Gm = mp.quad(lambda u: g(u) * mp.e ** (-z * u), [-L, 0, L])
    return Gp * Gm


def main():
    mp.mp.dps = 40
    L = mp.mpf('1.7')
    print('Phase C unit tests: orbit conventions')
    for sector in ('even_d', 'odd'):
        print(' sector = %s' % sector)
        idx = [0, 2, 5]
        coeffs = [mp.mpf('0.6'), mp.mpf('-0.3'), mp.mpf('0.75')]
        gam = mp.mpf('3.4')
        beta = mp.mpf('0.77')
        dlt = beta - mp.mpf(1) / 2

        # 1. on-line conjugate pair
        want = mp.fsum([H_brute(coeffs, idx, L, sector, mp.mpc(0, s * gam))
                        for s in (1, -1)])
        got = dh.online_pair_contribution(coeffs, idx, L, sector, gam)
        check('on-line pair {1/2 +- i gamma} vs brute force', got, mp.re(want),
              mp.mpf('1e-28'))
        check('  ... and it is >= 0', mp.sign(got), 1, 0)

        # 2. full off-line quartet, brute force over the four shifts
        shifts = [mp.mpc(dlt, gam), mp.mpc(-dlt, gam),
                  mp.mpc(dlt, -gam), mp.mpc(-dlt, -gam)]
        want = mp.fsum([H_brute(coeffs, idx, L, sector, z) for z in shifts])
        got = dh.quartet_contribution(coeffs, idx, L, sector, beta, gam)
        check('off-line quartet vs brute force over 4 shifts', got, mp.re(want),
              mp.mpf('1e-27'))
        check('  ... imaginary part of the brute-force sum vanishes',
              mp.im(want), 0, mp.mpf('1e-28'))

        # 3. the |P|^2 - |S|^2 form
        g = g_fn(coeffs, idx, L, sector)
        P = mp.quad(lambda u: g(u) * mp.cosh(dlt * u) * mp.e ** (1j * gam * u),
                    [-L, 0, L])
        S = mp.quad(lambda u: g(u) * mp.sinh(dlt * u) * mp.e ** (1j * gam * u),
                    [-L, 0, L])
        check('quartet = 4(|P|^2 - |S|^2)', got, 4 * (abs(P) ** 2 - abs(S) ** 2),
              mp.mpf('1e-27'))

        # 4. NOT |G(delta + i gamma)|^2
        bad = abs(dh.G_at(coeffs, idx, L, sector, mp.mpc(dlt, gam))) ** 2
        print('  %-58s %s' % ('(control) 4|G(delta+i gamma)|^2 differs from the orbit',
                              'YES' if abs(got - 4 * bad) > mp.mpf('1e-6') else 'NO'))
        if abs(got - 4 * bad) <= mp.mpf('1e-6'):
            FAIL.append(('quartet-vs-|G|^2 control', 0, 0))

        # 5. collapse delta -> 0 gives TWO DOUBLE on-line zeros, not a simple pair
        lim = dh.online_pair_contribution(coeffs, idx, L, sector, gam)
        for e in ('1e-4', '1e-6', '1e-8'):
            q = dh.quartet_contribution(coeffs, idx, L, sector,
                                        mp.mpf(1) / 2 + mp.mpf(e), gam)
            check('collapse delta=%s -> 2 x (double on-line pair)' % e, q, 2 * lim,
                  mp.mpf(e) ** 2 * 100 + mp.mpf('1e-25'))

        # 6. simple real off-line pair at gamma = 0 (TWO zeros, not four)
        want = mp.fsum([H_brute(coeffs, idx, L, sector, s * dlt) for s in (1, -1)])
        got = dh.real_pair_contribution(coeffs, idx, L, sector, beta)
        check('real pair {beta, 1-beta} vs brute force (2 zeros)', got, mp.re(want),
              mp.mpf('1e-27'))
        # sign: + in the even sector, - in the odd sector (dirichlet-windows result)
        sgn = 1 if sector != 'odd' else -1
        check('  ... sign is %+d in the %s sector' % (sgn, sector),
              mp.sign(got), sgn, 0)

    # 7. mixed even-odd terms vanish for the full real self-reciprocal form:
    #    an even g and an odd g' have <g, g'> = 0 and their window matrices never mix,
    #    which we check by evaluating the orbit of a quartet on g + g' and comparing
    #    with the sum of the two separate evaluations.
    print(' parity decoupling')
    L = mp.mpf('1.3')
    ce, io = [mp.mpf('0.8'), mp.mpf('0.2')], [0, 3]
    co, jo = [mp.mpf('-0.4'), mp.mpf('0.9')], [1, 2]
    beta, gam = mp.mpf('0.65'), mp.mpf('5.1')

    def g_mix(u):
        return (g_fn(ce, io, L, 'even_d')(u) + g_fn(co, jo, L, 'odd')(u))

    def Hmix(z):
        Gp = mp.quad(lambda u: g_mix(u) * mp.e ** (z * u), [-L, 0, L])
        Gm = mp.quad(lambda u: g_mix(u) * mp.e ** (-z * u), [-L, 0, L])
        return Gp * Gm
    d = beta - mp.mpf(1) / 2
    tot = mp.re(mp.fsum([Hmix(z) for z in (mp.mpc(d, gam), mp.mpc(-d, gam),
                                           mp.mpc(d, -gam), mp.mpc(-d, -gam))]))
    sep = (dh.quartet_contribution(ce, io, L, 'even_d', beta, gam)
           + dh.quartet_contribution(co, jo, L, 'odd', beta, gam))
    check('quartet on (even + odd) = quartet(even) + quartet(odd)', tot, sep,
          mp.mpf('1e-25'))

    print('')
    if FAIL:
        print('FAILURES: %d' % len(FAIL))
        for f in FAIL:
            print('  ', f)
        return 1
    print('ALL ORBIT UNIT TESTS PASS')
    return 0


if __name__ == '__main__':
    sys.exit(main())
