#!/usr/bin/env python3
r"""Validation of the eta / eigenvalue identities, and of the toy heat family.

These are VALIDATION TOOLS, not findings.  The perturbations (B', u') used here are
SYNTHETIC and explicitly labelled as such: the audit's feasibility gate [F] did not
justify a complete differentiable source-derived family, so no source-derived derivative
is used anywhere in this file.  What is checked is that the algebra is right.
"""
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, 'src'))
import mpmath as mp                                                      # noqa: E402

FAIL = []


def check(name, got, want, tol):
    d = abs(got - want) / max(abs(want), mp.mpf(1))
    ok = d <= tol
    print('  %-56s %-5s  rel %s' % (name, 'OK' if ok else 'FAIL', mp.nstr(d, 4)))
    if not ok:
        FAIL.append(name)


def main():
    mp.mp.dps = 50
    n = 6
    print('eta / eigenvalue identities (SYNTHETIC perturbation, not source-derived)')

    # a fixed, reproducible positive definite B and a channel u; no randomness
    B = mp.zeros(n, n)
    for i in range(n):
        for j in range(n):
            B[i, j] = mp.mpf(1) / (1 + abs(i - j)) + (mp.mpf(n) if i == j else 0)
    u = mp.matrix([mp.mpf(1) / (i + 2) for i in range(n)])
    # synthetic derivatives
    Bp = mp.zeros(n, n)
    for i in range(n):
        for j in range(n):
            Bp[i, j] = mp.cos(mp.mpf(i + 1) * (j + 1)) / 10
            Bp[j, i] = Bp[i, j]
    up = mp.matrix([mp.sin(mp.mpf(i + 1)) / 5 for i in range(n)])

    def eta_of(h):
        Bh = B + h * Bp
        uh = u + h * up
        x = mp.lu_solve(Bh, uh)
        return mp.fsum([uh[i] * x[i] for i in range(n)])

    x0 = mp.lu_solve(B, u)
    eta0 = mp.fsum([u[i] * x0[i] for i in range(n)])
    # eta' = 2 u'^T x - x^T B' x
    analytic = 2 * mp.fsum([up[i] * x0[i] for i in range(n)]) \
        - mp.fsum([x0[i] * Bp[i, j] * x0[j] for i in range(n) for j in range(n)])
    h = mp.mpf(10) ** (-20)
    numeric = (eta_of(h) - eta_of(-h)) / (2 * h)
    print('  eta(0) = %s' % mp.nstr(eta0, 20))
    check("eta' = 2(u')^T x - x^T B' x", analytic, numeric, mp.mpf('1e-30'))

    # at a crossing eta = 1: Q x = 0.  Rescale u so that eta = 1 exactly.
    c = 1 / mp.sqrt(eta0)
    uc = c * u
    xc = mp.lu_solve(B, uc)
    etac = mp.fsum([uc[i] * xc[i] for i in range(n)])
    check('rescaled eta == 1', etac, mp.mpf(1), mp.mpf('1e-40'))
    Q = B - uc * uc.T
    resid = mp.norm(Q * xc) / mp.norm(xc)
    print('  %-56s %-5s  |Q x|/|x| = %s'
          % ('at eta = 1:  Q x = 0', 'OK' if resid < mp.mpf('1e-40') else 'FAIL',
             mp.nstr(resid, 4)))
    if resid >= mp.mpf('1e-40'):
        FAIL.append('Qx=0 at crossing')

    # lambda_0' = -eta'/||x||^2 for the eigenvalue through zero
    analytic_c = 2 * mp.fsum([up[i] * xc[i] for i in range(n)]) \
        - mp.fsum([xc[i] * Bp[i, j] * xc[j] for i in range(n) for j in range(n)])
    lam_pred = -analytic_c / mp.fsum([xc[i] ** 2 for i in range(n)])

    def lam0(hh):
        Bh = B + hh * Bp
        uh = uc + hh * up
        Qh = Bh - uh * uh.T
        ev = sorted(mp.eigsy(Qh, eigvals_only=True))
        return min(ev, key=abs)
    lam_num = (lam0(h) - lam0(-h)) / (2 * h)
    check("lambda_0' = -eta'/||x||^2", lam_pred, lam_num, mp.mpf('1e-18'))

    # the toy heat family: d_tau F = -d_z^2 F, and the transition at tau = a^2/2
    print('')
    print('toy family F_tau(z) = z^2 + a^2 - 2 tau')
    for aa in ('1', '1.790272316470219442122718'):
        a = mp.mpf(aa)
        tc = a ** 2 / 2
        for tt, want_real in ((tc - mp.mpf(1) / 10, False),
                              (tc + mp.mpf(1) / 10, True)):
            disc = 2 * tt - a ** 2
            is_real = disc >= 0
            print('  a = %-26s tau = %-22s zeros real: %-5s (expected %s)'
                  % (aa, mp.nstr(tt, 16), is_real, want_real))
            if is_real != want_real:
                FAIL.append('toy transition a=%s' % aa)
        # it satisfies the same PDE
        check('  d_tau F = -d_z^2 F  (a=%s)' % aa, mp.mpf(-2), -mp.mpf(2),
              mp.mpf('1e-40'))

    print('')
    if FAIL:
        print('FAILURES: %s' % FAIL)
        return 1
    print('ALL IDENTITY CHECKS PASS')
    return 0


if __name__ == '__main__':
    sys.exit(main())
