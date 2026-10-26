"""H_tau from the s-derivatives of Lambda_F -- the representation the fixed-parameter
certificate needs.

Since  d_tau H = -d_z^2 H  and  s = 1/2 + i z/2  (so d/dz = (i/2) d/ds, d_z^2 = -(1/4) d_s^2):

    H_tau(z) = sum_{m>=0} (tau^m/m!) (-1)^m d_z^{2m} H_0(z)
             = (1/2) sum_{m>=0} (tau/4)^m Lambda_F^{(2m)}(s) / m!            (*)

so H_tau at a FIXED tau is a power series in tau whose coefficients are even s-derivatives
of Lambda_F -- each of which Arb encloses rigorously via acb_series.  No oscillatory
quadrature and no comparison with H_0 is needed, which is exactly what made the parent
extension's perturbation bound weak (it replaced cos(zx) by e^{Mx} and discarded every bit
of cancellation).

THE TAIL IS MAJORISED IN CLOSED FORM.  From  H_0^{(2m)}(z) = (-1)^m int_0^inf Phi_F(u)
u^{2m} cos(zu) du,
    |H_0^{(2m)}(z)|  <=  int_0^inf |Phi_F(u)| u^{2m} cosh(u |Im z|) du  =:  M_{2m} ,
and  sum_m (tau^m/m!) M_{2m} = int_0^inf |Phi_F(u)| cosh(u|Im z|) e^{tau u^2} du  <  inf.
So the tail beyond m = M is bounded by

    T_M(tau, Y)  =  (tau^{M+1}/(M+1)!) int_0^inf |Phi_F(u)| u^{2M+2} cosh(uY) e^{tau u^2} du ,

a NON-OSCILLATORY integral of a positive function, with Y >= |Im z| on the contour.  (The
bound uses  sum_{m>M} a_m <= a_{M+1} sum_{k>=0} (tau u^2)^k/k! absorbed into the e^{tau u^2}
factor, i.e. m!/(M+1)! >= (m-M-1)! termwise.)
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from mpmath import mp, mpf, mpc, exp, cos, cosh, pi, quad, factorial
import kernels as K


def Lambda_derivs(s, M, kap=None):
    """Lambda_F^{(0)}, ..., Lambda_F^{(2M)} at s, by mpmath numerical differentiation of
    the analytic assembly.  Used for NAVIGATION only; the certificate uses Arb."""
    from dh_flow import Lambda_log_derivs
    if kap is None:
        kap = K.kappa_closed()
    f = lambda w: K.Lambda_DH(w, kap)
    return [mp.diff(f, s, n) for n in range(0, 2*M + 1, 2)]


def H_taylor(z, tau, M=8, kap=None):
    """H_tau(z) by (*), truncated at m = M.  Navigation only."""
    s = mpf(1)/2 + mpc(0, 1)*z/2
    d = Lambda_derivs(s, M, kap)
    acc = mpc(0)
    for m in range(M + 1):
        acc += (mpf(tau)/4)**m * d[m] / factorial(m)
    return acc/2


def tail_majorant(tau, Y, M, kap=None, xmax=3):
    """T_M(tau, Y) as above, by quadrature (navigation); positive, non-oscillatory."""
    if kap is None:
        kap = K.kappa_closed()
    f = lambda u: abs(K.Phi_DH(u, kap)) * u**(2*M + 2) * cosh(u*mpf(Y)) * exp(mpf(tau)*u*u)
    I = quad(f, [0, mpf(1)/2, 1, mpf(3)/2, mpf(xmax)])
    return mpf(tau)**(M + 1)/factorial(M + 1) * I
