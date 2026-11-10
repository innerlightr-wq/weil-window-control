"""R_T(f)/||f||^2 for a FROZEN vector, by direct scalar quadrature with a
certified Bernstein remainder.  No matrix of any size is assembled.

    R_T(f) = beta ||f||^2 + 2 s F(i/2)^2
             + (1/pi) int_0^{T} (Psi_L(t) - beta) F(t)^2 dt ,
    F(t)   = sum_n y_n v_n(t) ,   v_n(t) = sqrt(2L(2n+1)) j_n(Lt) ,

with y the y-coordinates (the same convention as the matrices).  Even sector.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "src"))
sys.path.insert(0, HERE)
from flint import arb, acb
import block as BK
from prec import owned, require_tight

def R_scalar(L_num, L_den, y_exact, degs, T, beta, panels, order, prec,
             b_e="0.30", progress=None):
    """Certified enclosure of R_T(f)/||f||^2.

    INPUTS ARE REBUILT INSIDE THE SCOPE.  Passing a pre-built `L` or `Symbol`
    is unsafe: an arb constructed at the ambient precision carries that
    precision's radius, and entering a higher-precision scope afterwards does
    NOT recover the lost information.  That is how a 53-bit `arb(4)/5` turned
    the integral radius into 1e105 during this run's first attempt -- the same
    defect class as the ctx.prec leak, in input construction rather than scope.
    `y_exact` must be EXACT data (Fraction or a decimal/rational string), not
    pre-built arbs: a coefficient vector built at the ambient precision carries
    that precision's radius into a quantity that needs ~17 digits of
    cancellation.  That cost two attempts in this run -- first a 53-bit
    `arb(4)/5` (radius 1e105), then 53-bit coefficients (radius 1.8e-15).
    `beta` is checked for tightness rather than rebuilt, since it is a certified
    constant supplied by `constants()`.
    """
    from fractions import Fraction
    with owned(prec):
        require_tight(arb(beta), "beta (input)", 1e-25)
        L = arb(L_num) / arb(L_den)
        S = BK.Symbol(L_num, L_den, T)
        y = []
        for v in y_exact:
            f = v if isinstance(v, Fraction) else Fraction(str(v))
            y.append(arb(f.numerator) / arb(f.denominator))
        nodes, wts = BK.gl_rule(order)
        hw = arb(T) / (2 * panels)
        rho = BK.panel_rho(hw, arb(b_e))
        efac = BK.panel_error_factor(hw, rho, order)
        rad_re = hw * (rho + 1 / rho) / 2
        VMAX = BK.v_hol_absmax_uniform(L, arb(b_e))
        nrm2 = sum((v * v for v in y), arb(0))
        ynorm = nrm2.sqrt()
        # |F(z)| <= ||y|| * sqrt(sum v_n^2) <= ||y|| * sqrt(#modes) * VMAX on the strip
        Fbound = ynorm * arb(len(degs)).sqrt() * VMAX
        maxdeg = max(degs)
        tot = arb(0)
        errsum = arb(0)
        for ip in range(panels):
            c = hw * (2 * ip + 1)
            box = acb(arb(c.mid(), rad_re.upper()), arb(0, arb(b_e).upper()))
            Wm = arb(S.W_hol_at(box, beta).abs_upper())
            errsum += efac * Wm * Fbound * Fbound
            s = arb(0)
            for q in range(order):
                t = c + hw * nodes[q]
                vv = BK.v_real(L, t, maxdeg + 1)
                F = sum((y[i] * vv[degs[i]] for i in range(len(degs))), arb(0))
                s += wts[q] * (S.W_real_at(t, beta)) * F * F
            tot += hw * s
            if progress and ip % progress == 0:
                print("      panel %d/%d" % (ip, panels), flush=True)
        integral = arb((tot / arb.pi()).mid(),
                       (arb((tot / arb.pi()).rad()) + errsum / arb.pi()).upper())
        p = BK.pole_vector(L, maxdeg + 1)          # y-coordinates
        Fh = sum((y[i] * p[degs[i]] for i in range(len(degs))), arb(0))
        R = beta * nrm2 + 2 * Fh * Fh + integral
        require_tight(R / nrm2, "R_T/||f||^2", 1e-22)
        return R / nrm2, dict(quad_err=(errsum / arb.pi()).str(10), rho=rho.str(8),
                              panels=panels, order=order)
