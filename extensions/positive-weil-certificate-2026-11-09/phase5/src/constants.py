"""Certificate constants with OWNED precision, replacing the leaky sector_dk."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "src"))
sys.path.insert(0, HERE)
from flint import arb
import baseline as B
from prec import owned, PrecCache, require_tight

_CACHE = PrecCache()

def constants(L_num, L_den, T, N, prec, tail_extra=200):
    """All precision-dependent certificate constants, computed inside owned(prec).

    Returns a dict of Arb balls.  Cached by (L, T, N, tail_extra, prec).
    """
    key = (L_num, L_den, T, N, tail_extra)

    def build():
        L = arb(L_num) / arb(L_den)
        A_L = B.A_L(L_num, L_den)
        beta = B.H(T) - A_L
        M = B.M_sup(L_num, L_den, beta, T)
        x = arb(T) * L
        S = B.S_tail(N, x, top=N + tail_extra)
        require_tight(S, "S_N", 1e-50)
        g = (2 * L * arb(T) / arb.pi()) * S
        pt2 = B.pole_vec_tail(N, L)
        pn = (2 * L.sinh()).sqrt()
        kappa = 2 * pn * pt2.sqrt() + M * g.sqrt()
        d_even = beta - M * g
        d_odd = beta - M * g - 2 * pt2
        for nm, v, r in (("beta", beta, 1e-20), ("M", M, 1e-20),
                         ("kappa", kappa, 1e-30),
                         ("d_even", d_even, 1e-20), ("d_odd", d_odd, 1e-20)):
            require_tight(v, nm, r)
        return dict(L=L, A_L=A_L, beta=beta, M=M, S=S, g=g, p_perp2=pt2,
                    p_norm=pn, kappa=kappa, d_even=d_even, d_odd=d_odd)

    return _CACHE.get(key, prec, build)
