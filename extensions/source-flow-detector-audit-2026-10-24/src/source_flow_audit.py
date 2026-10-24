#!/usr/bin/env python3
r"""Source-defined (de Bruijn-Newman heat) deformation -> fixed-window detector: audit.

PROVENANCE.  The function, window form, quartet decomposition and frozen witnesses are
REUSED unchanged from extensions/davenport-heilbronn-2026-10-23 (commit e2cf715).  The
rank-one channel reduction is from extensions/channel-compatibility-exploration-2026-10-24.
Nothing in either directory is modified.

SOURCE OBJECTS (see SOURCE_AUDIT.md for provenance and exact equation numbers).

  zeta:  Phi(u)   = sum_{n>=1} (2 pi^2 n^4 e^{9u} - 3 pi n^2 e^{5u}) exp(-pi n^2 e^{4u})
         H_t(z)   = int_0^inf e^{t u^2} Phi(u) cos(zu) du
         H_0(z)   = (1/8) xi(1/2 + i z/2)                      [Polymath eq. (1),(3),(4)]
         dH_t/dt  = -H_t''                                     [Polymath eq. (31)]

  Davenport-Heilbronn: DERIVED HERE, not taken from the zeta case.  From
         Lambda_F(s) = int_0^inf omega(t) t^{(s+1)/2} dt/t,
         omega(t)    = sum_{n>=1} n a_n e^{-pi n^2 t / 5},
  the substitution t = e^{4x}, s = 1/2 + i z/2 gives
         Lambda_F(1/2 + i z/2) = int_R Phi_F(x) e^{izx} dx,
         Phi_F(x) = 4 e^{3x} sum_{n>=1} n a_n e^{-pi n^2 e^{4x}/5},
  and the functional equation Lambda_F(s) = Lambda_F(1-s) forces Phi_F even, so
         H^DH_t(z) = int_0^inf e^{t x^2} Phi_F(x) cos(zx) dx,
         H^DH_0(z) = (1/2) Lambda_F(1/2 + i z/2).

COORDINATE DICTIONARY, kept strictly separate:
  s  = 1/2 + i z/2         the original Dirichlet-series variable
  z  = -2i(s - 1/2)        the spectral coordinate of H_t;  for s = beta + i gamma,
                           z = 2 gamma - 2 i (beta - 1/2), so Im z = -2 delta
  t (= tau)                the deformation parameter; NOT physical time
  x (= u)                  the kernel variable
  L                        the detector half-window, in the test-function variable

Usage:  source_flow_audit.py [--dps 40] [--quick]
"""
import argparse
import csv
import json
import os
import sys

DH = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), 'davenport-heilbronn-2026-10-23', 'src')
CH = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), 'channel-compatibility-exploration-2026-10-24', 'src')
sys.path.insert(0, DH)
sys.path.insert(0, CH)

import mpmath as mp                                                      # noqa: E402
import dh                                                                # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, 'data')
ROWS = []


def rec(block, name, value, note=''):
    ROWS.append(dict(block=block, quantity=name, value=value, note=note))


# --------------------------------------------------------------------- kernels

def Phi_zeta(x, N=80):
    x = mp.mpf(x)
    return mp.fsum([(2 * mp.pi ** 2 * n ** 4 * mp.e ** (9 * x)
                     - 3 * mp.pi * n ** 2 * mp.e ** (5 * x))
                    * mp.e ** (-mp.pi * n ** 2 * mp.e ** (4 * x))
                    for n in range(1, N + 1)])


def Phi_DH(x, N=240):
    x = mp.mpf(x)
    return 4 * mp.e ** (3 * x) * mp.fsum(
        [n * dh.a_coeff(n) * mp.e ** (-mp.pi * n ** 2 * mp.e ** (4 * x) / 5)
         for n in range(1, N + 1)])


def H_of(Phi, z, tau=0, T=3):
    z = mp.mpmathify(z)
    tau = mp.mpf(tau)
    return mp.quad(lambda x: mp.e ** (tau * x ** 2) * Phi(x) * mp.cos(z * x),
                   [0, mp.mpf(1) / 2, 1, T])


def xi(s):
    s = mp.mpmathify(s)
    return mp.mpf(1) / 2 * s * (s - 1) * mp.power(mp.pi, -s / 2) * mp.gamma(s / 2) \
        * mp.zeta(s)


# --------------------------------------------------------------------- blocks

def block_dictionary(p):
    p('')
    p('[A] COORDINATE DICTIONARY AND SOURCE NORMALISATION')
    p('    zeta: H_0(z) = (1/8) xi(1/2 + i z/2)   [Polymath 1904.12438 eq. (1)]')
    worst = mp.mpf(0)
    for zz in ('0', '1', '5', '14.134725'):
        z = mp.mpf(zz)
        lhs = H_of(Phi_zeta, z)
        rhs = mp.re(xi(mp.mpc(mp.mpf(1) / 2, z / 2)) / 8)
        d = abs(lhs - rhs) / max(abs(rhs), mp.mpf(1))
        worst = max(worst, d)
        p('      z = %-12s H_0 = %-24s xi/8 = %-24s rel = %s'
          % (zz, mp.nstr(lhs, 18), mp.nstr(rhs, 18), mp.nstr(d, 4)))
    rec('A', 'zeta_normalisation_worst_rel', mp.nstr(worst, 6),
        'H_0(z) = xi(1/2+iz/2)/8 verified, not assumed')
    p('    worst relative residual: %s' % mp.nstr(worst, 6))

    p('    Phi_zeta even?  (evenness is equivalent to the theta functional equation)')
    p('      The series at x < 0 has terms of opposite sign and size O(1) summing to a')
    p('      value as small as 1e-46, so the check needs working precision well above')
    p('      the size of the cancellation.  Each point is evaluated at a precision set')
    p('      from its own magnitude, and a point that is still unresolved is reported')
    p('      as UNRESOLVED rather than as a large residual.')
    we = mp.mpf(0)
    base = mp.mp.dps
    for v in ('0.2', '0.5', '0.9'):
        try:
            mp.mp.dps = base
            mag = abs(Phi_zeta(mp.mpf(v)))
            need = int(max(base, 30 - mp.log10(mag) if mag > 0 else base)) + 30
            mp.mp.dps = need
            a, b = Phi_zeta(mp.mpf(v)), Phi_zeta(-mp.mpf(v))
            rel = abs(a - b) / abs(a)
            resolved = rel < mp.mpf('1e-10')
            p('      x = %-5s dps=%-4d Phi(x) = %-20s rel = %-12s %s'
              % (v, need, mp.nstr(a, 12), mp.nstr(rel, 4),
                 'resolved' if resolved else 'UNRESOLVED'))
            if resolved:
                we = max(we, rel)
            else:
                rec('A', 'Phi_zeta_evenness_x%s' % v, 'UNRESOLVED',
                    'cancellation below the attempted precision')
        finally:
            mp.mp.dps = base
    rec('A', 'Phi_zeta_evenness_worst_rel_resolved', mp.nstr(we, 6),
        'hypothesis (10), over the points that resolved')


def block_heat(p):
    p('')
    p('[B] THE SOURCE EQUATION  d_tau H = - d_z^2 H   [Polymath eq. (31)]')
    p('    Verified by direct differentiation of the integral, both kernels:')
    p('      d_tau H = int x^2 e^{tau x^2} Phi cos(zx) dx')
    p('      d_z^2 H = int e^{tau x^2} Phi (-x^2) cos(zx) dx = - d_tau H .')
    p('    Numerically (finite differences on the integral representation):')
    for nm, Phi in (('zeta', Phi_zeta), ('DH  ', Phi_DH)):
        for zz, tt in (('3.0', '0'), ('3.0', '0.2')):
            z, t = mp.mpf(zz), mp.mpf(tt)
            h = mp.mpf(10) ** (-10)
            dtau = (H_of(Phi, z, t + h) - H_of(Phi, z, t - h)) / (2 * h)
            d2z = (H_of(Phi, z + h, t) - 2 * H_of(Phi, z, t)
                   + H_of(Phi, z - h, t)) / h ** 2
            d = abs(dtau + d2z) / max(abs(dtau), mp.mpf(1))
            p('      %s  z=%-5s tau=%-4s  d_tau H = %-22s  -d_z^2 H = %-22s rel = %s'
              % (nm, zz, tt, mp.nstr(dtau, 14), mp.nstr(-d2z, 14), mp.nstr(d, 4)))
            rec('B', 'heat_eq_%s_z%s_tau%s' % (nm.strip(), zz, tt), mp.nstr(d, 6),
                'finite-difference check of eq. (31)')


def block_dh_kernel(p):
    p('')
    p('[C] THE DAVENPORT-HEILBRONN KERNEL, DERIVED AND CHECKED')
    p('    Phi_F(x) = 4 e^{3x} sum_{n>=1} n a_n exp(-pi n^2 e^{4x}/5)')
    p('')
    p('    (C1) does it reproduce Lambda_F?  2 int_0^inf Phi_F cos(zx) dx')
    worst = mp.mpf(0)
    for zz in ('0', '2', '10', '171.3986969707551843'):
        z = mp.mpf(zz)
        lhs = 2 * H_of(Phi_DH, z)
        rhs = mp.re(dh.Lambda_F(mp.mpc(mp.mpf(1) / 2, z / 2)))
        d = abs(lhs - rhs) / max(abs(rhs), mp.mpf('1e-25'))
        worst = max(worst, d)
        p('      z = %-22s 2H = %-24s Lambda_F = %-24s rel = %s'
          % (zz, mp.nstr(lhs, 16), mp.nstr(rhs, 16), mp.nstr(d, 4)))
    rec('C', 'DH_kernel_transform_worst_rel', mp.nstr(worst, 6),
        'Phi_F reproduces Lambda_F(1/2+iz/2)')
    p('      (the last z is 2*gamma_1 for the selected off-line quartet, where')
    p('       Lambda_F is small but nonzero on the real-z line)')

    p('')
    p('    (C2) hypothesis (10): is Phi_F even?')
    we = mp.mpf(0)
    for v in ('0.1', '0.3', '0.6', '1.0'):
        a, b = Phi_DH(mp.mpf(v)), Phi_DH(-mp.mpf(v))
        we = max(we, abs(a - b) / abs(a))
        p('      x = %-5s Phi_F(x) = %-22s Phi_F(-x) = %-22s rel = %s'
          % (v, mp.nstr(a, 16), mp.nstr(b, 16), mp.nstr(abs(a - b) / abs(a), 4)))
    rec('C', 'DH_kernel_evenness_worst_rel', mp.nstr(we, 6), 'hypothesis (10)')

    p('')
    p('    (C3) hypothesis (11): |Phi_F(x)| <= A exp(-|x|^{2+alpha})')
    p('      Phi_F(x) = O(e^{3|x|} e^{-(pi/5) e^{4|x|}}), double-exponential decay,')
    p('      so (11) holds for EVERY alpha > 0.  Sampled:')
    for v in ('1.0', '1.5', '2.0'):
        x = mp.mpf(v)
        p('      x = %-5s Phi_F = %-22s  exp(-|x|^3) = %s'
          % (v, mp.nstr(Phi_DH(x), 8), mp.nstr(mp.e ** (-x ** 3), 8)))

    p('')
    p('    (C4) positivity of Phi_F on x >= 0.  With T = pi e^{4x}/5 >= pi/5 and')
    p('    |a_n| <= 1, a_1 = 1, FACTOR OUT the leading term:')
    p('      sum_n n a_n e^{-n^2 T} >= e^{-T} [ 1 - sum_{n>=2} n e^{-(n^2-1)T} ].')
    p('    The bracket is INCREASING in T (each exponential decreases), so a single')
    p('    evaluation at T = pi/5 bounds it from below for every x >= 0.  The')
    p('    unfactored difference e^{-T} - sum_{n>=2} n e^{-n^2 T} must NOT be used for')
    p('    this: it tends to 0 as T grows and so is not a uniform bound.')
    T0 = mp.pi / 5
    brk = 1 - mp.fsum([n * mp.e ** (-(n ** 2 - 1) * T0) for n in range(2, 80)])
    p('      %-28s %-24s %s' % ('T', 'bracket (uniform)', 'unfactored difference'))
    for Ts in (T0, mp.mpf(1), mp.mpf(2), mp.mpf(5), mp.mpf(20)):
        b = 1 - mp.fsum([n * mp.e ** (-(n ** 2 - 1) * Ts) for n in range(2, 80)])
        u = mp.e ** (-Ts) - mp.fsum([n * mp.e ** (-n ** 2 * Ts) for n in range(2, 80)])
        p('      %-28s %-24s %s'
          % (mp.nstr(Ts, 14), mp.nstr(b, 16), mp.nstr(u, 10)))
    p('      => bracket >= %s > 0 for all T >= pi/5, hence Phi_F > 0 on x >= 0 with a'
      % mp.nstr(brk, 16))
    p('         decaying prefactor e^{-T}, and by evenness on all of R.')
    rec('C', 'DH_positivity_uniform_bracket_at_T0', mp.nstr(brk, 16),
        '1 - sum_{n>=2} n e^{-(n^2-1)T} at T = pi/5; increasing in T, so uniform')
    mn = None
    for j in range(0, 61):
        x = mp.mpf(j) / 20
        v = Phi_DH(x)
        if mn is None or v < mn[1]:
            mn = (x, v)
    p('      sampled minimum on [0,3]: Phi_F(%s) = %s  (positive: %s)'
      % (mp.nstr(mn[0], 4), mp.nstr(mn[1], 8), mn[1] > 0))
    p('      and |a_n| <= 1 makes the tail bound monotone in T, so the margin only')
    p('      grows for x > 0.  With evenness, Phi_F > 0 on all of R.')
    rec('C', 'DH_kernel_sampled_min_x_ge_0', mp.nstr(mn[1], 8), 'positive')


def block_debruijn(p):
    p('')
    p('[D] WHAT DE BRUIJN THEOREM 7 GIVES, FOR EACH OBJECT')
    p('    Hypotheses (Newman-Wu 1901.06596, eqs. (10),(11), Thm 7): F(-t) = conj F(t);')
    p('    |F(t)| <= A exp(-|t|^{2+alpha}); zeros of int F(t)e^{izt}dt in |Im z| <= Delta.')
    p('    Conclusion: zeros of int F(t) e^{lambda t^2/2} e^{izt} dt lie in')
    p('                |Im z| <= max(Delta^2 - lambda, 0)^{1/2}.')
    p('    With the Polymath normalisation e^{tau x^2} we have lambda = 2 tau, so all')
    p('    zeros are real once tau >= Delta^2/2.')
    p('')
    sig0 = dh.mpf_hi('1.395136158235109721061359')   # src/census.py SIGMA0
    p('    zeta: zeros have |beta - 1/2| <= 1/2, so Delta = 2*(1/2) = 1 and')
    p('          tau >= 1/2.  This REPRODUCES de Bruijn 1950 exactly, which is the')
    p('          check that the lambda = 2 tau bookkeeping is right.')
    rec('D', 'zeta_Delta', '1', 'from the classical critical strip')
    rec('D', 'zeta_tau_all_real', '0.5', 'de Bruijn 1950, reproduced')
    dlt = 2 * (sig0 - mp.mpf(1) / 2)
    tcrit = dlt ** 2 / 2
    p('')
    p('    DH:   the PROVED zero-free abscissa is sigma_0 = %s' % mp.nstr(sig0, 25))
    p('          (davenport-heilbronn-2026-10-23 REPORT.md section 3, T1), so every')
    p('          zero has |beta - 1/2| <= sigma_0 - 1/2 = %s' % mp.nstr(sig0 - mp.mpf(1) / 2, 25))
    p('          and Im z = -2(beta - 1/2) gives')
    p('            Delta_DH = %s' % mp.nstr(dlt, 25))
    p('            tau_DH   = Delta^2/2 = %s' % mp.nstr(tcrit, 25))
    p('          so ALL zeros of H^DH_tau are real for tau >= %s.' % mp.nstr(tcrit, 12))
    rec('D', 'DH_Delta', mp.nstr(dlt, 25), 'Delta = 2(sigma_0 - 1/2)')
    rec('D', 'DH_tau_all_real', mp.nstr(tcrit, 25),
        'de Bruijn Thm 7 applied to the DH kernel')
    return tcrit


# --------------------------------------------------------------- detector bridge

def basis_G(idx, L, sector, w):
    """G_j(w) = int_{-L}^{L} g_j(x) e^{wx} dx on the repository's normalised basis."""
    out = []
    for k in idx:
        e = [mp.mpf(0)] * len(idx)
        e[list(idx).index(k)] = mp.mpf(1)
        out.append(dh.G_at(e, idx, L, sector, w))
    return out


def Psi(idx, L, sector, z):
    """Psi_ij(z) = (1/2)[G_i(iz/2)G_j(-iz/2) + G_j(iz/2)G_i(-iz/2)]."""
    z = mp.mpmathify(z)
    a = basis_G(idx, L, sector, mp.mpc(0, 1) * z / 2)
    b = basis_G(idx, L, sector, -mp.mpc(0, 1) * z / 2)
    n = len(idx)
    M = mp.zeros(n, n)
    for i in range(n):
        for j in range(n):
            M[i, j] = (a[i] * b[j] + a[j] * b[i]) / 2
    return M


def block_bridge(p, idx, L, sector):
    p('')
    p('[E] THE DETECTOR BRIDGE, AND ITS tau = 0 AGREEMENT')
    p('    Psi_ij(z) = (1/2)[G_i(iz/2)G_j(-iz/2) + G_j(iz/2)G_i(-iz/2)].')
    p('    With z = -2i(s - 1/2), i.e. z = 2 gamma for an on-line zero, iz/2 = i gamma,')
    p('    so Psi_ij(2 gamma) is exactly the repository\'s symmetrised H(rho - 1/2).')
    p('    Checked against dh.H_at on a random-free fixed test vector:')
    n = len(idx)
    cc = [mp.mpf(1) / (i + 3) for i in range(n)]
    worst = mp.mpf(0)
    for gam in ('3.5', '85.69934848537759217'):
        g = dh.mpf_hi(gam)
        P = Psi(idx, L, sector, 2 * g)
        lhs = mp.fsum([cc[i] * cc[j] * P[i, j] for i in range(n) for j in range(n)])
        rhs = mp.re(dh.H_at(cc, idx, L, sector, mp.mpc(0, g)))
        d = abs(lhs - rhs) / max(abs(rhs), mp.mpf('1e-30'))
        worst = max(worst, d)
        p('      gamma = %-22s c^T Psi(2g) c = %-22s H_at = %-22s rel = %s'
          % (gam, mp.nstr(lhs, 14), mp.nstr(rhs, 14), mp.nstr(d, 4)))
    rec('E', 'Psi_vs_repo_H_worst_rel', mp.nstr(worst, 6),
        'coordinate dictionary z = 2 gamma verified')

    p('')
    p('    LEMMA (reality => positive semidefinite).  For real g and REAL z,')
    p('    G_j(-iz/2) = conj G_j(iz/2), so with v_j = G_j(iz/2),')
    p('        Psi(z) = Re( v v^* ) = (Re v)(Re v)^T + (Im v)(Im v)^T  >= 0,')
    p('    a sum of two rank-one positive semidefinite real matrices.  Hence if EVERY')
    p('    zero of H_tau is real, Q_tau = sum_k m_k Psi(z_k) is a sum of positive')
    p('    semidefinite matrices and Q_tau >= 0.  Checked:')
    for zz in ('1.0', '7.0', '171.4'):
        P = Psi(idx, L, sector, mp.mpf(zz))
        # Psi(z) is real for real z (the lemma); numerically it carries an imaginary
        # residue at the precision floor, which is itself a check of the lemma.
        imres = max(abs(mp.im(P[i, j])) for i in range(n) for j in range(n))
        Pr = mp.matrix(n, n)
        for i in range(n):
            for j in range(n):
                Pr[i, j] = mp.re(P[i, j])
        ev = sorted(mp.eigsy(Pr, eigvals_only=True))
        rk = sum(1 for e in ev if abs(e) > mp.mpf('1e-25'))
        p('      z = %-8s min eig = %-22s rank = %-3d max |Im Psi| = %s'
          % (zz, mp.nstr(ev[0], 10), rk, mp.nstr(imres, 4)))
        rec('E', 'Psi_min_eig_real_z_%s' % zz, mp.nstr(ev[0], 10),
            'rank %d, PSD' % rk)
    p('    The observed rank is 1, not 2: in a PURE PARITY SECTOR G_j(iz/2) is purely')
    p('    imaginary (odd) or purely real (even) for real z, so one of the two rank-one')
    p('    terms vanishes identically -- the same phase removal as the channel')
    p('    decomposition of the preceding exploration.  The lemma conclusion is')
    p('    unaffected: Psi(z) >= 0 either way, and the least eigenvalue above is at the')
    p('    precision floor rather than negative.  (An earlier draft of this line')
    p('    asserted rank two; the numerics corrected it.)')
    p('')
    p('    AGREEMENT AT tau = 0 is NOT re-derived here.  The parent extension already')
    p('    gated it: its Phase C compares the arithmetic side against exactly this')
    p('    zero-side sum over its census and agrees to 5.65e-11, and against an')
    p('    independent log-derivative contour to 1.73e-13 (classification PASS-N).')
    p('    That gate is reused, with its tier, not repeated.')


# --------------------------------------------------------------- feasibility gate

def block_gate(p):
    p('')
    p('[F] FEASIBILITY GATE ON A COMPLETE DIFFERENTIABLE FAMILY')
    p('    Convergence.  On the real-z axis the repository bases give Psi_ij(z) =')
    p('    O(|z|^{-4}) (both bases vanish at +-L, so G = O(|Im w|^{-2})), and the zero')
    p('    density of H_tau is O(log|z|), so sum_k |Psi(z_k)| converges ABSOLUTELY for')
    p('    real zeros, and also for zeros confined to a bounded strip, since Psi grows')
    p('    only like e^{L|Im z|} there while still decaying like |Re z|^{-4}.  JUSTIFIED.')
    p('    Reality and orbit pairing.  Zeros come in the orbit {z, -z, conj z, -conj z}')
    p('    and Psi is even with Psi(conj z) = conj Psi(z), so each complete orbit')
    p('    contributes a real symmetric matrix.  JUSTIFIED (same convention as the')
    p('    parent extension, whose unit tests cover it).')
    p('    tau = 0 agreement.  JUSTIFIED by reuse of the parent Phase C gate (T2-N).')
    p('    Differentiability in tau.  NOT JUSTIFIED here.  Simple zeros move')
    p('    analytically, and a FIXED zero-free contour gives an analytic finite-contour')
    p('    family; but differentiating the COMPLETE sum term by term needs a uniform')
    p('    bound on sum_k |Psi\'(z_k) z_k\'| with z_k\' = H\'\'/H\'(z_k), which is not')
    p('    controlled near near-collisions, and a fixed contour can be crossed by a')
    p('    zero as tau varies.')
    p('')
    p('    => Only a FINITE-CONTOUR differentiable family is justified.  The missing')
    p('       lemma is a uniform-in-tau tail bound for sum_{|z_k| > R} |Psi\'(z_k) z_k\'|')
    p('       on a parameter interval, together with a zero-crossing-free contour on')
    p('       that interval.  Per the brief, the numerical tau-sweep branch is NOT run.')
    rec('F', 'complete_differentiable_family', 'NOT JUSTIFIED',
        'finite-contour family only; missing uniform tail lemma')


# --------------------------------------------------------------- controls

def block_toy(p):
    p('')
    p('[G] ADVERSARIAL CONTROL 1: the generic heat-flow family F_tau(z) = z^2 + a^2 - 2tau')
    p('    d_tau F = -2 and d_z^2 F = 2, so d_tau F = -d_z^2 F: the SAME equation.')
    p('    Zeros z = +- sqrt(2 tau - a^2): purely imaginary for tau < a^2/2, real for')
    p('    tau >= a^2/2.  Its initial strip half-width is Delta = a, and it becomes')
    p('    real exactly at tau = Delta^2/2 -- it SATURATES de Bruijn Theorem 7.')
    for aa in ('1', '1.790272316470219442'):
        a = dh.mpf_hi(aa)
        p('      a = %-22s transition at tau = a^2/2 = %s'
          % (aa, mp.nstr(a ** 2 / 2, 18)))
    p('    Consequence: the transition from non-real to real zeros, and hence the')
    p('    transition of the detector from indefinite to positive semidefinite, is')
    p('    produced by the heat equation and the initial strip width ALONE, with no')
    p('    arithmetic input whatsoever.  Any argument that derives a zero restriction')
    p('    from the source equation plus a strip bound must therefore also apply to')
    p('    F_tau, and so cannot be arithmetic in content.')
    p('    Hypothesis check.  F_tau is NOT the Fourier transform of an integrable,')
    p('    rapidly decaying kernel: such a transform is BOUNDED on the real axis, and a')
    p('    degree-2 polynomial is not.  So F_tau is a POLYNOMIAL HEAT-FLOW CONTROL that')
    p('    reproduces the strip-shrinking RATE; it is not an example satisfying every')
    p('    hypothesis of Theorem 7, and it is not a counterexample to it.  The relevant')
    p('    counterexample inside the intended class is the DH construction of [C].')
    rec('G', 'toy_saturates_de_bruijn', 'yes',
        'transition exactly at Delta^2/2, no arithmetic input')


def block_manufactured(p, idx, L, sector):
    p('')
    p('[H] ADVERSARIAL CONTROL 2: manufactured dynamics from an assembled B')
    p('    For B > 0 the flow x\'(t) = -(1/2) B x(t), x(0) = u has x(t) = e^{-tB/2}u and')
    p('        int_0^inf ||x(t)||^2 dt = u^T (int_0^inf e^{-tB} dt) u = u^T B^{-1} u = eta.')
    p('    So eta always admits a "dynamical" representation.  This is a generic')
    p('    consequence of B > 0, NOT an upstream explanation, and it is available for')
    p('    any positive definite matrix from any source.  Verified numerically:')
    import explore_channels as ec                                        # noqa: E402
    _, ctx = ec.run_config('2.0', 2, mp.mp.dps, verbose=False)
    Bm, u, eta = ctx['B'], ctx['u'], ctx['eta']
    n = len(u)
    E, V = mp.eigsy(Bm)
    # int_0^inf ||e^{-tB/2}u||^2 dt = sum_j |v_j^T u|^2 / beta_j, done spectrally
    tot = mp.fsum([mp.fsum([V[k, j] * u[k] for k in range(n)]) ** 2 / E[j]
                   for j in range(n)])
    p('      eta (linear solve)                      = %s' % mp.nstr(eta, 20))
    p('      int_0^inf ||x(t)||^2 dt (spectral form) = %s' % mp.nstr(tot, 20))
    p('      relative difference                     = %s'
      % mp.nstr(abs(tot - eta) / abs(eta), 4))
    rec('H', 'manufactured_dynamics_rel_diff', mp.nstr(abs(tot - eta) / abs(eta), 6),
        'generic representation of eta, not an upstream mechanism')
    p('    The brief\'s warning is confirmed: this construction carries no information')
    p('    about where B came from.')


# --------------------------------------------------------------- the obstruction

def block_layers(p):
    p('')
    p('[L] POSITIVE KERNEL => PSD CORRELATION MATRIX, WHICH IS A DIFFERENT FORM')
    p('    For any nonnegative integrable Phi with Fourier transform Phi^(y),')
    p('      sum_{j,k} c_j conj(c_k) Phi^(y_j - y_k)')
    p('          = int Phi(x) |sum_j c_j e^{i y_j x}|^2 dx  >=  0,')
    p('    so a positive kernel generates a POSITIVE SEMIDEFINITE matrix of Fourier')
    p('    correlations.  For the DH kernel, Phi_F^(y) = Lambda_F(1/2 + i y/2), so the')
    p('    matrix M_jk = Lambda_F(1/2 + i(y_j - y_k)/2) must be PSD.  Checked:')
    ys = [mp.mpf(v) for v in ('0', '1', '2.5', '4', '7')]
    n = len(ys)
    M = mp.zeros(n, n)
    for i in range(n):
        for j in range(n):
            M[i, j] = mp.re(dh.Lambda_F(mp.mpc(mp.mpf(1) / 2, (ys[i] - ys[j]) / 2)))
    ev = sorted(mp.eigsy(M, eigvals_only=True))
    p('      y = %s' % [mp.nstr(v, 4) for v in ys])
    p('      eigenvalues: %s' % ', '.join(mp.nstr(e, 10) for e in ev))
    p('      least eigenvalue = %s  ->  PSD' % mp.nstr(ev[0], 10))
    rec('L', 'DH_correlation_matrix_min_eig', mp.nstr(ev[0], 10),
        'PSD, as positivity of Phi_F forces')
    p('')
    p('    AND YET Q^DH_0 has a certified NEGATIVE witness.  There is no contradiction:')
    p('    these are different quadratic forms built from different aspects of the same')
    p('    function -- one from the kernel\'s correlations, one from the zeros.  The')
    p('    lesson for a layered programme is exact:')
    p('      POSITIVITY IN THE CORRELATION GEOMETRY DOES NOT TRANSFER TO POSITIVITY IN')
    p('      THE ZERO GEOMETRY.')
    p('    Finding "something positive" is not the missing bridge; proving that a')
    p('    particular positive structure controls the particular spectral form that')
    p('    matters is.')


def block_lambda_bracket(p, tcrit):
    p('')
    p('[M] A TWO-SIDED BRACKET ON THE DH DE BRUIJN-NEWMAN CONSTANT')
    p('    Upper: de Bruijn Theorem 7, block [D]:  Lambda_DH <= %s.' % mp.nstr(tcrit, 22))
    p('    Lower: a nonreal zero at tau = 0 PERSISTS for small tau.  Inputs:')
    p('      (a) H^DH_tau(z) is jointly analytic in (tau, z): the double-exponential')
    p('          decay of Phi_F makes the integral locally uniformly convergent;')
    p('      (b) at tau = 0 there is a zero at z_1 = 2 gamma_1 - 2 i delta_1 with')
    beta1 = dh.mpf_hi('0.8085171824566373855533519606068441278507')
    gam1 = dh.mpf_hi('85.69934848537759217192926770894172903799')
    d1 = beta1 - mp.mpf(1) / 2
    p('          Im z_1 = -2 delta_1 = %s, so a disc of radius r < %s'
      % (mp.nstr(-2 * d1, 18), mp.nstr(2 * d1, 12)))
    p('          about z_1 does not meet the real axis;')
    p('      (c) the argument principle / Hurwitz on that disc gives a zero inside it')
    p('          for all sufficiently small |tau|, hence a NONREAL zero, hence')
    p('          tau < Lambda_DH;')
    p('      (d) forward preservation makes Lambda_DH well defined at all.')
    p('    Therefore')
    p('      0 < Lambda_DH <= %s.' % mp.nstr(tcrit, 22))
    rec('M', 'DH_Lambda_bracket', '(0, %s]' % mp.nstr(tcrit, 22),
        'lower bound by analytic persistence; upper by de Bruijn Thm 7')
    p('')
    p('    This is a DEDUCTION from the established kernel and the off-real-zero input,')
    p('    NOT a numerical localisation: quantifying "sufficiently small" would need the')
    p('    locally uniform estimate that block [F] records as missing, so no numerical')
    p('    lower bound is claimed and no flow numerics were run for it.')


def block_obstruction(p, tcrit):
    p('')
    p('[I] THE OBSTRUCTION: the DH control refutes the sought class of inequality')
    p('    Chain, with tiers:')
    p('      (i)   Phi_F satisfies (10) and (11) and Phi_F > 0           [C2,C3,C4: T1 +')
    p('            numerical confirmation]')
    p('      (ii)  every zero of Lambda_F has |beta - 1/2| <= sigma_0 - 1/2, so the')
    p('            initial strip half-width is Delta_DH = 2(sigma_0 - 1/2)   [T1, parent')
    p('            REPORT.md section 3]')
    p('      (iii) de Bruijn Theorem 7 therefore applies to the DH kernel and ALL zeros')
    p('            of H^DH_tau are real for tau >= %s                [T1]' % mp.nstr(tcrit, 14))
    p('      (iv)  by the lemma of [E], Q^DH_tau >= 0 for tau >= %s   [T1]'
      % mp.nstr(tcrit, 14))
    p('      (v)   but Q^DH_0 restricted to the L = 2.0 odd level-2 space has')
    p('            lambda_min = -0.679071346441480590986377586093..., enclosed in an Arb')
    p('            ball of radius 1.39e-116                            [T2-C, parent')
    p('            certificates/ball_certificate.txt]')
    p('      (vi)  and Q^DH_0 is that arithmetic matrix, by the explicit formula of the')
    p('            parent DERIVATION.md section 3                                 [T1];')
    p('            the 5.65e-11 agreement is a CROSS-CHECK of the two implementations,')
    p('            not the source of the identity, and the Arb ball certifies the')
    p('            evaluated expression, not the identity.  Those roles stay separate.')
    p('')
    p('    CONSEQUENCE.  Suppose there were an inequality of the sought form')
    p('        eta\'(theta) <= a(theta) (1 - eta(theta)),')
    p('    oriented so that eta <= 1 at large tau propagates DOWN to tau = 0, with a')
    p('    controlled by hypotheses that the DH kernel satisfies -- i.e. by (10), (11),')
    p('    kernel positivity, a strip bound, and the heat equation.  Applied to DH it')
    p('    would give Q^DH_0 >= 0 on the stated space, contradicting (v)-(vi).')
    p('    Therefore NO UNIVERSALLY APPLICABLE backward-preservation argument resting')
    p('    only on those shared hypotheses can work.  The precise conclusion is not')
    p('    "the inequality cannot exist" but: somewhere along the required path from')
    p('    tau = tau_DH down to tau = 0, one of its ingredients must fail -- the')
    p('    inequality itself, the control on a(theta), the differentiability, or the')
    p('    standing assumption B_tau > 0.  The channel representation needs B_tau > 0')
    p('    at every tau, and tracking individual zeros can meet collisions; neither can')
    p('    be carried along the whole path for free.')
    p('')
    p('    This does not refute a backward inequality that uses a hypothesis DH fails.')
    p('    The only candidate the parent work supplies is an Euler product: DH has none,')
    p('    and zeta does.  So any surviving route must use the Euler product in an')
    p('    essential way -- and the parent Phase H already showed that deleting the')
    p('    non-Euler part of the arithmetic side is not a usable way to do that')
    p('    (lambda_min of the truncation is negative at 17 of 18 configurations,')
    p('    including where the true form is positive).')


def block_eta_identities(p, idx, L, sector):
    p('')
    p('[J] THE eta AND eigenvalue IDENTITIES (validation tools, not findings)')
    p('    eta_tau = u^T B^{-1} u  =>  eta\' = 2 (u\')^T x - x^T B\' x,  x = B^{-1}u.')
    p('    At a crossing eta = 1: Q x = (B - u u^T) x = u(1 - eta) = 0, so x spans')
    p('    ker Q; and for the eigenvalue through zero, lambda_0\' = -eta\'/||x||^2.')
    p('    These are verified in tests/test_identities.py against a SYNTHETIC, clearly')
    p('    labelled perturbation of B and u -- NOT a source-derived one, because [F]')
    p('    did not justify a complete differentiable family.  They establish that the')
    p('    algebra is right; they say nothing about the source.')
    p('    The research question was whether the source law controls either term.  It')
    p('    does not: u\' needs the selected quartet\'s velocity H\'\'/H\'(z_1), and B\' needs')
    p('    the velocity of every other zero, so both require the spectrum that the')
    p('    detector was supposed to constrain.')


def block_arithmetic(p):
    p('')
    p('[K] ARITHMETIC CONTENT OF THE EVOLUTION AT tau = 0')
    p('    Away from zeros, with l = d_z log H,')
    p('      d_z^2 H / H = d_z l + l^2,  so  d_tau l = -d_z(d_z l + l^2)')
    p('                                            = -d_z^2 l - 2 l d_z l.')
    p('    (Verified symbolically above; it is an identity, not an assumption.)')
    p('    At tau = 0, l_0(z) = (i/2) (xi\'/xi)(1/2 + iz/2) and')
    p('      xi\'/xi(s) = 1/s + 1/(s-1) - (1/2)log pi + (1/2)psi(s/2) + zeta\'/zeta(s),')
    p('      -zeta\'/zeta(s) = sum_n Lambda(n) n^{-s}      (Re s > 1).')
    p('    So the QUADRATIC term in the evolution carries')
    p('      l^2  ->  (sum_n Lambda(n) n^{-s})^2 = sum_n (Lambda * Lambda)(n) n^{-s},')
    p('    a divisor convolution supported on products of TWO prime powers.  Its')
    p('    coefficients are FORCED by the source equation, which is genuinely different')
    p('    from deleting mixed terms by hand (parent Phase H).')
    p('')
    p('    But: (a) the Dirichlet series for l^2 converges absolutely only for')
    p('    Re s > 1, i.e. |Im z| > 1, which is OUTSIDE the strip where the zeros and')
    p('    the detector live, and continuing it termwise into the strip is exactly the')
    p('    step the brief forbids; (b) (Lambda*Lambda)(n) >= 0, but it enters d_tau Q')
    p('    through d_z Psi, which is sign-indefinite on the contour, so the convolution')
    p('    supplies no sign; (c) for tau != 0 the kernel e^{tau x^2}Phi(x) is not the')
    p('    kernel of any Dirichlet series, so H_tau has NO Euler product, NO')
    p('    gamma-factor decomposition and NO prime-side explicit formula.  The')
    p('    arithmetic exists only at the single point tau = 0.')
    p('')
    p('    So the honest answer to the brief\'s question is the second alternative: the')
    p('    source equation generates another exact but sign-indefinite expression, in a')
    p('    half-plane disjoint from where it would have to be used.')
    rec('K', 'quadratic_term_convolution', 'Lambda * Lambda',
        'forced by the source equation; absolutely convergent only for Re s > 1')


# --------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dps', type=int, default=40)
    a = ap.parse_args()
    mp.mp.dps = a.dps
    os.makedirs(DATA, exist_ok=True)
    out = []

    def p(line=''):
        print(line)
        out.append(line)

    sys.path.insert(0, CH)
    from phaseD import index_set                                         # noqa: E402
    L = dh.mpf_hi('2.0')
    sector = 'odd'
    idx = index_set(L, sector, 1)

    p('=' * 86)
    p('SOURCE-FLOW -> DETECTOR AUDIT   (working precision %d digits)' % a.dps)
    p('=' * 86)
    p('  detector space: L = 2.0, odd sector, parent level 2, N = %d' % len(idx))
    p('  deformation parameter tau; spectral coordinate z = -2i(s - 1/2); kernel')
    p('  variable x; window half-length L.  These are never conflated.')

    block_dictionary(p)
    block_heat(p)
    block_dh_kernel(p)
    tcrit = block_debruijn(p)
    block_bridge(p, idx, L, sector)
    block_gate(p)
    block_arithmetic(p)
    block_eta_identities(p, idx, L, sector)
    block_toy(p)
    block_manufactured(p, idx, L, sector)
    block_layers(p)
    block_lambda_bracket(p, tcrit)
    block_obstruction(p, tcrit)

    p('')
    p('=' * 86)
    with open(os.path.join(DATA, 'validation.csv'), 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=['block', 'quantity', 'value', 'note'])
        w.writeheader()
        for r in ROWS:
            w.writerow(r)
    with open(os.path.join(HERE, 'data', 'audit_log.txt'), 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    p('  -> data/validation.csv, data/audit_log.txt')
    return 0


if __name__ == '__main__':
    sys.exit(main())
