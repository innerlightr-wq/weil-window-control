#!/usr/bin/env python3
r"""A CERTIFIED lower bound for Lambda_DH by local persistence (Rouche), in Arb.

Follows the review's section 3 argument exactly.  Nothing here touches the global
differentiability material of block [F] or runs any flow numerics: the only tau-dependence
used is the explicit perturbation bound below, evaluated once.

THE ARGUMENT.
  H_tau(z) = int_0^inf e^{tau x^2} Phi_F(x) cos(zx) dx,   H_0(z) = (1/2) Lambda_F(1/2+iz/2).
  Let D be a closed disc about the first off-line zero z_1 of H_0, with closure disjoint
  from the real axis.  Put

      m = min_{dD} |H_0|,     M = max_{dD} |Im z|,
      C = int_0^inf x^2 |Phi_F(x)| e^{tau_0 x^2 + M x} dx.

  For |tau| <= tau_0 and z on dD,

      |H_tau(z) - H_0(z)| = | int_0^inf (e^{tau x^2} - 1) Phi_F(x) cos(zx) dx |
                         <= |tau| int_0^inf x^2 |Phi_F(x)| e^{tau_0 x^2} e^{M x} dx
                          = |tau| C,

  using  |e^{tau x^2} - 1| <= |tau| x^2 e^{tau_0 x^2}  for |tau| <= tau_0,  which holds
  because |e^u - 1| = |int_0^u e^t dt| <= |u| e^{|u|} and |u| = |tau| x^2 <= tau_0 x^2;
  and |cos(zx)| <= e^{|Im z| x} <= e^{M x}.

  So if |tau| C < m, Rouche gives H_tau the same zero count in D as H_0.  With a zero of
  H_0 in D (certified below by the winding number), H_tau then has a zero in D, which is
  NONREAL because D misses the real axis, hence tau < Lambda_DH.  Therefore

      Lambda_DH >= eps := min(tau_0, m/(2C)).

WHAT IS CERTIFIED.  m, C and the winding number are Arb enclosures (T2-C).  The identity
H_0 = (1/2) Lambda_F(1/2 + iz/2) is the Mellin derivation of REPORT.md section 2 (T1).
The disc centre is a numerical input from the parent census, but NOTHING depends on it
being a zero: the winding number certifies a zero inside D by itself.
"""
import argparse
import csv
import os
import sys
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, getcontext

from flint import acb, arb, ctx                                          # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROWS = []


def rec(name, value, note=''):
    ROWS.append(dict(quantity=name, value=str(value), note=note))


# --------------------------------------------------- printing a ONE-SIDED bound honestly
#
# `arb.str(n)` rounds the midpoint to nearest, so "m >= <that>" can assert slightly MORE
# than the enclosure supports -- for this run by 4.5e-50, harmless numerically and wrong
# in principle.  A certified floor must round DOWN and a certified ceiling UP, from the
# ball's own endpoint.  These return plain decimals, with no +/- part, that are true on
# every point of the ball.

def _endpoint(v, upper):
    getcontext().prec = 220
    mid = Decimal(v.mid().str(90, radius=False))
    rad = Decimal(v.rad().str(30, radius=False))
    return mid + rad if upper else mid - rad


def floor_str(v, n):
    """Plain decimal with n significant digits, <= every point of the ball `v`."""
    lo = _endpoint(v, False)
    if lo == 0:
        return '0'
    return str(lo.quantize(Decimal(1).scaleb(lo.adjusted() - n + 1),
                           rounding=ROUND_FLOOR)).replace('E', 'e')


def ceil_str(v, n):
    """Plain decimal with n significant digits, >= every point of the ball `v`."""
    hi = _endpoint(v, True)
    if hi == 0:
        return '0'
    return str(hi.quantize(Decimal(1).scaleb(hi.adjusted() - n + 1),
                           rounding=ROUND_CEILING)).replace('E', 'e')


# ------------------------------------------------------------------ the function, in Arb

def kappa():
    """[sqrt(10 - 2 sqrt 5) - 2]/(sqrt 5 - 1); pure square roots, exactly enclosable."""
    r5 = arb(5).sqrt()
    return ((10 - 2 * r5).sqrt() - 2) / (r5 - 1)


def a_vec():
    k = kappa()
    return [arb(1), k, -k, arb(-1), arb(0)]


def Lambda_F(s):
    """(5/pi)^{(s+1)/2} Gamma((s+1)/2) 5^{-s} sum_{r=1}^{4} a_r zeta(s, r/5)."""
    av = a_vec()
    tot = acb(0)
    for r in range(1, 5):
        tot += av[r - 1] * acb(s).zeta(acb(arb(r) / 5))
    pre = (acb(5) / acb.pi()) ** ((s + 1) / 2)
    return pre * ((s + 1) / 2).gamma() * acb(5) ** (-s) * tot


def H0(z):
    """H_0(z) = (1/2) Lambda_F(1/2 + i z/2)."""
    s = acb(arb(1) / 2) + acb(0, 1) * z / 2
    return Lambda_F(s) / 2


# ------------------------------------------------------------------ the kernel bound

NKER = 40


def Phi_F_bound_terms(x):
    """Phi_F(x) as an enclosure: 4 e^{3x} [ sum_{n<=N} n a_n e^{-pi n^2 e^{4x}/5} +- tail ].

    For x >= 0, e^{4x} >= 1, so the omitted terms are bounded by
      sum_{n>N} n e^{-pi n^2/5} <= 2 (N+1) e^{-pi (N+1)^2/5}
    (consecutive terms fall by a factor at most e^{-pi} < 1/20, so the geometric
    majorant costs a factor 20/19 < 2, and the linear growth of n is absorbed in it)."""
    av = a_vec()
    E = (4 * x).exp()
    tot = acb(0)
    for n in range(1, NKER + 1):
        tot += n * av[(n - 1) % 5] * (-acb.pi() * n * n * E / 5).exp()
    tail = 2 * (NKER + 1) * (-arb.pi() * (NKER + 1) ** 2 / 5).exp()
    tot += acb(arb(0, tail.mid() + tail.rad()))
    return 4 * (3 * x).exp() * tot


def Phi_F_upper(x):
    """An upper bound for |Phi_F(x)|, used only in the tail argument:
       |Phi_F(x)| <= 8 e^{3x} e^{-pi e^{4x}/5}  for x >= 1,
    since sum_{n>=1} n e^{-n^2 T} <= 2 e^{-T} once T >= 34 (and T = pi e^{4x}/5 >= 34
    for x >= 1)."""
    return 8 * (3 * x).exp() * (-arb.pi() * (4 * x).exp() / 5).exp()


# ------------------------------------------------------------------ main certification

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--prec', type=int, default=400)
    ap.add_argument('--arcs', type=int, default=400)
    ap.add_argument('--radius', default='0.3')
    ap.add_argument('--tau0', default='0.1')
    ap.add_argument('--xcut', default='2')
    a = ap.parse_args()
    ctx.prec = a.prec
    out = []

    def p(line=''):
        print(line)
        out.append(line)

    gam1 = arb('85.69934848537759217192926770894172903799')
    bet1 = arb('0.8085171824566373855533519606068441278507')
    del1 = bet1 - arb(1) / 2
    cre, cim = 2 * gam1, -2 * del1
    r = arb(a.radius)
    tau0 = arb(a.tau0)
    X = arb(a.xcut)

    p('=' * 86)
    p('CERTIFIED LOWER BOUND FOR Lambda_DH BY LOCAL PERSISTENCE   (Arb, %d bits)' % a.prec)
    p('=' * 86)
    p('')
    p('[1] THE DISC')
    p('    z_1 = 2 gamma_1 - 2 i delta_1  (critical line = real axis in this variable)')
    p('      Re z_1 = %s' % cre.str(30))
    p('      Im z_1 = %s' % cim.str(30))
    p('    radius r              = %s' % r.str(10))
    p('    margin from real axis = |Im z_1| - r = %s' % (-cim - r).str(20))
    p('    closure disjoint from the real axis: %s' % bool((-cim - r) > 0))
    rec('disc_centre_re', cre.str(30))
    rec('disc_centre_im', cim.str(30))
    rec('disc_radius', r.str(10))
    rec('disc_margin_from_real_axis', (-cim - r).str(20),
        'must be > 0 for the disc to miss the real axis')
    if not bool((-cim - r) > 0):
        p('    STOP: the disc is not disjoint from the real axis.')
        return 1

    # ---------------------------------------------------------------- m and winding
    N = a.arcs
    dth = 2 * arb.pi() / N
    arcrad = 2 * r * (dth / 4).sin()
    p('')
    p('[2] m = min_{dD} |H_0|, AND THE WINDING NUMBER, BY ARC ENCLOSURES')
    p('    partition: %d arcs of angular width 2 pi/%d' % (N, N))
    p('    each arc {c + r e^{i th} : |th - th_m| <= dth/2} is contained in the ball of')
    p('    centre c + r e^{i th_m} and radius 2 r sin(dth/4) = %s,' % arcrad.str(12))
    p('    since |r e^{i th} - r e^{i th_m}| = 2 r |sin((th - th_m)/2)|.')
    mlow = None
    args, uncs, rels = [], [], []
    for k in range(N):
        thm = (2 * k + 1) * arb.pi() / N
        mid = acb(cre + r * thm.cos(), cim + r * thm.sin())
        zb = acb(arb(mid.real.mid(), arcrad.mid() + arcrad.rad()),
                 arb(mid.imag.mid(), arcrad.mid() + arcrad.rad()))
        v = H0(zb)
        lo = v.abs_lower()
        if mlow is None or lo < mlow:
            mlow = lo
        # argument interval: for a ball of midpoint mu and radius rho < |mu|,
        # arg lies within arg(mu) +- asin(rho/|mu|)
        mu = acb(v.real.mid(), v.imag.mid())
        rho = arb(max(v.real.rad(), v.imag.rad())) * arb(2).sqrt()
        amu = mu.abs_lower()
        args.append(mu.arg())
        rels.append(rho / amu)
        uncs.append((rho / amu).asin() if bool(rho < amu) else arb.pi())
    p('    worst relative enclosure radius over the arcs: %s'
      % max(rels, key=lambda t: t.mid()).str(8))
    p('')
    p('    m enclosure = %s' % mlow.str(20))
    p('    m >= %s      (plain decimal floor: the ball\'s lower endpoint,' % floor_str(mlow, 20))
    p('                                       rounded DOWN, so the inequality is exact)')
    p('    m > 0 : %s   (so dD passes through no zero)' % bool(mlow > 0))
    rec('m_lower_bound', mlow.str(25), 'min over dD of |H_0|, Arb arc enclosures')
    rec('m_floor_20', floor_str(mlow, 20), 'plain decimal floor for m, rounded DOWN')
    rec('n_arcs', N)
    rec('arc_ball_radius', arcrad.str(12))
    if not bool(mlow > 0):
        p('    STOP: m could not be certified positive.')
        return 1

    # winding number from the certified argument intervals
    tot = arb(0)
    worst = arb(0)
    ok = True
    for k in range(N):
        d = args[(k + 1) % N] - args[k]
        while bool(d > arb.pi()):
            d = d - 2 * arb.pi()
        while bool(d < -arb.pi()):
            d = d + 2 * arb.pi()
        slack = d.abs_upper() + uncs[k] + uncs[(k + 1) % N]
        if bool(slack >= arb.pi()):
            ok = False
        if bool(slack > worst):
            worst = slack
        tot = tot + d
    wind = tot / (2 * arb.pi())
    p('[3] WINDING NUMBER OF H_0 ON dD')
    p('    sum of principal argument increments / 2 pi = %s' % wind.str(15))
    p('    worst (|increment| + both arc uncertainties) = %s  (must be < pi = %s)'
      % (worst.str(10), arb.pi().str(10)))
    p('    increments unambiguous: %s' % ok)
    nz = round(float(wind.mid()))
    p('    => H_0 has exactly %d zero(s) in D, counted with multiplicity' % nz)
    rec('winding_number', nz, 'certified: all increments unambiguous' if ok else 'AMBIGUOUS')
    rec('winding_worst_slack', worst.str(10), 'must be < pi')
    if not ok or nz < 1:
        p('    STOP: the winding number is not certified to be >= 1.')
        return 1

    # ---------------------------------------------------------------- M and C
    M = -cim + r
    p('')
    p('[4] M = max_{dD} |Im z| = |Im z_1| + r = %s' % M.str(20))
    rec('M_max_abs_Im_z', M.str(20))
    p('')
    p('[5] C = int_0^inf x^2 |Phi_F(x)| e^{tau_0 x^2 + M x} dx,  tau_0 = %s'
      % tau0.str(6))
    p('    split at x = %s.  On [0, X] by rigorous Arb integration of the explicit' % X.str(4))
    p('    integrand, with the kernel truncated at n <= %d and its omitted terms' % NKER)
    p('    carried as an error ball.')

    def integrand(x, _):
        return x * x * Phi_F_bound_terms(x) * (tau0 * x * x + M * x).exp()

    head = acb.integral(integrand, 0, float(X.mid()))
    p('      int_0^X = %s' % head.real.str(20))

    p('    Tail x > X.  For x >= 1, |Phi_F(x)| <= 8 e^{3x} e^{-pi e^{4x}/5}, so with')
    p('      f(x) = x^2 8 e^{3x} e^{-pi e^{4x}/5} e^{tau_0 x^2 + M x},')
    p('      (log f)\'(x) = 2/x + 3 - (4 pi/5) e^{4x} + 2 tau_0 x + M <= -1  for x >= X,')
    p('    because (4 pi/5) e^{4X} = %s dominates 2/X + 3 + 2 tau_0 x + M + 1 on x >= X'
      % (4 * arb.pi() / 5 * (4 * X).exp()).str(8))
    p('    and grows far faster; hence f(x) <= f(X) e^{-(x-X)} and int_X^inf f <= f(X).')
    fX = X * X * Phi_F_upper(X) * (tau0 * X * X + M * X).exp()
    p('      f(X) = %s   (the tail bound)' % fX.str(10))
    Cup = arb(head.real.mid(), head.real.rad()) + fX
    Cup = arb(Cup.mid() + Cup.rad())
    p('    C enclosure = %s' % Cup.str(20))
    p('    C <= %s      (plain decimal ceiling, rounded UP)' % ceil_str(Cup, 20))
    rec('C_head_integral', head.real.str(20), 'rigorous Arb integration on [0,X]')
    rec('C_tail_bound', fX.str(10), 'f(X), via (log f)\' <= -1 on x >= X')
    rec('C_upper_bound', Cup.str(20))
    rec('C_ceil_20', ceil_str(Cup, 20), 'plain decimal ceiling for C, rounded UP')
    rec('x_cut', X.str(4))

    # ---------------------------------------------------------------- eps
    p('')
    p('[6] THE BOUND')
    eps_r = arb(mlow.mid() - mlow.rad()) / (2 * Cup)
    eps_r = arb(eps_r.mid() - eps_r.rad())
    eps = eps_r if bool(eps_r < tau0) else tau0
    p('    m/(2C) >= %s' % floor_str(eps_r, 20))
    p('    tau_0   = %s' % tau0.str(6))
    p('    eps = min(tau_0, m/(2C)) >= %s' % floor_str(eps, 20))
    p('')
    p('    CERTIFIED INTERVAL:')
    p('      %s <= Lambda_DH <= 1.602537483559822777787819' % floor_str(eps, 20))
    p('')
    p('    Tiers: eps is T2-C (m, C and the winding number are Arb enclosures; the')
    p('    identity H_0 = Lambda_F(1/2+iz/2)/2 is the T1 Mellin derivation of')
    p('    REPORT.md section 2).  The upper bound is T1 GIVEN the strip bound')
    p('    sigma_0 = 1.395136158235109721061359, itself T1 in the parent extension.')
    rec('eps_lower_bound', eps.str(25), 'T2-C: Lambda_DH >= eps')
    rec('eps_floor_20', floor_str(eps, 20),
        'plain decimal floor for eps, rounded DOWN; this is the number to quote')
    rec('upper_bound', '1.602537483559822777787819', 'T1 given the strip bound')
    p('=' * 86)

    os.makedirs(os.path.join(HERE, 'data'), exist_ok=True)
    with open(os.path.join(HERE, 'data', 'lambda_bracket.csv'), 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=['quantity', 'value', 'note'])
        w.writeheader()
        for row in ROWS:
            w.writerow(row)
    with open(os.path.join(HERE, 'data', 'lambda_bracket.txt'), 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    p('  -> data/lambda_bracket.csv, data/lambda_bracket.txt')
    return 0


if __name__ == '__main__':
    sys.exit(main())
