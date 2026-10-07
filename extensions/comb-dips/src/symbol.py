r"""Zhu's prime-comb symbol (arXiv:2608.24827, eq. (3)).

    Psi_L(t) = Re psi(1/4 + i t/2) - log pi - P_L(t)
    P_L(t)   = sum_{n : log n < 2L} (2 Lambda(n)/sqrt n) cos(t log n)
    A_L      = sum_{log n < 2L} 2 Lambda(n)/sqrt n
    T1(L)    = 2 pi e^{A_L}

Psi_L depends on L only through the active set {n : log n < 2L}, so it is piecewise
constant in L and changes only at L = (1/2) log n.
"""
import numpy as np
from scipy.special import digamma

LOG_PI = np.log(np.pi)


def von_mangoldt_terms(twoL):
    """[(n, Lambda(n), log n)] for prime powers n >= 2 with log n < 2L."""
    nmax = int(np.floor(np.exp(twoL))) + 2
    out = []
    for n in range(2, nmax + 1):
        m, p, lam = n, 2, None
        while p * p <= m:
            if m % p == 0:
                q = m
                while q % p == 0:
                    q //= p
                lam = np.log(p) if q == 1 else None
                break
            p += 1
        else:
            lam = np.log(m)
        if lam is not None and np.log(n) < twoL:
            out.append((n, lam, np.log(n)))
    return out


class Symbol:
    """The comb at a fixed L (i.e. a fixed active set)."""

    def __init__(self, L, freqs=None, weights=None, label=None):
        self.L = float(L)
        self.label = label or f"L={L}"
        if freqs is None:
            terms = von_mangoldt_terms(2 * self.L)
            self.n = np.array([t[0] for t in terms], dtype=np.int64)
            self.weights = np.array([2 * t[1] / np.sqrt(t[0]) for t in terms])
            self.freqs = np.array([t[2] for t in terms])          # log n
            # which prime each prime power belongs to, and its exponent
            self.prime = np.array([int(round(np.exp(t[1]))) for t in terms], dtype=np.int64)
            self.expo = np.array([int(round(t[2] / t[1])) for t in terms], dtype=np.int64)
        else:
            self.freqs = np.asarray(freqs, dtype=float)
            self.weights = np.asarray(weights, dtype=float)
            self.n = np.arange(len(self.freqs))
            self.prime = self.n
            self.expo = np.ones(len(self.freqs), dtype=np.int64)

    # --- scalar quantities -------------------------------------------------
    @property
    def A(self):
        return float(self.weights.sum())

    @property
    def T1(self):
        return 2 * np.pi * np.exp(self.A)

    @property
    def Tstar(self):
        return 2 * np.pi * np.exp(2 * self.L)

    def envelope_end(self, beta):
        """No dip beyond this: H(t) - A_L >= beta forces Psi_L >= beta."""
        return 2 * np.pi * np.exp(self.A + beta)

    def lipschitz_P(self):
        """sum |w_k| * freq_k  bounds |P_L'|."""
        return float(np.sum(np.abs(self.weights) * np.abs(self.freqs)))

    # --- the symbol --------------------------------------------------------
    def H(self, t):
        """Re psi(1/4 + i t/2) - log pi  (the smooth trend)."""
        t = np.asarray(t, dtype=float)
        return np.real(digamma(0.25 + 0.5j * t)) - LOG_PI

    def P_fast(self, t):
        """Plain float64. Loses ~|t log n| 2^-53 (1.3e-8 at L=1.6). Benchmarks only."""
        t = np.asarray(t, dtype=float)
        return np.cos(np.multiply.outer(t, self.freqs)) @ self.weights

    def P(self, t):
        """P_L(t) with exact argument reduction -- the default (see reduced_arg)."""
        return P_accurate(self, t)

    def Psi(self, t):
        return self.H(t) - self.P(t)

    def Psi_chunked(self, t, chunk=1 << 20):
        t = np.asarray(t, dtype=float)
        out = np.empty(t.shape, dtype=float)
        for i in range(0, t.size, chunk):
            s = slice(i, i + chunk)
            out[s] = self.H(t[s]) - P_accurate(self, t[s])
        return out

    def describe(self):
        return dict(label=self.label, L=self.L,
                    active_n=self.n.tolist(), A=self.A, T1=self.T1, Tstar=self.Tstar,
                    lipschitz_P=self.lipschitz_P())


# L values that sit strictly inside distinct prime-power sets
STAGE2_L = [0.55, 0.6, 0.7, 0.8, 0.9, 1.0, 1.19, 1.2, 1.4, 1.6]


# ---------------------------------------------------------------------------
# Accurate argument reduction.
#
# Plain float64 loses ~|t * log n| * 2^-53 in the cosine argument, i.e. ~1e-8 at
# t ~ 1.7e7 (L = 1.6).  That is harmless for dip detection (depths are O(1)) but it
# breaks the 1e-12 validation target and would have to be carried as an error term
# through the branch-and-bound.  Dekker two-product + two-term reduction mod 2*pi
# restores ~1e-13 at a small constant cost.
# ---------------------------------------------------------------------------
_SPLIT = 134217729.0          # 2^27 + 1
TWOPI_HI = 6.283185307179586
TWOPI_LO = 2.4492935982947064e-16      # 2*pi - TWOPI_HI, to double precision


def _two_prod(a, b):
    """Exact product: returns (p, e) with p + e == a*b exactly."""
    p = a * b
    ca = _SPLIT * a; ahi = ca - (ca - a); alo = a - ahi
    cb = _SPLIT * b; bhi = cb - (cb - b); blo = b - bhi
    e = ((ahi * bhi - p) + ahi * blo + alo * bhi) + alo * blo
    return p, e


def reduced_arg(t, f):
    """(t*f) mod 2pi to full double accuracy.

    Both products must be done exactly.  Taking q*TWOPI_HI in plain float64 is the
    trap: for t*f ~ 5e7 we have q ~ 8e6, and q*TWOPI_HI then needs ~77 bits, leaking
    ~|t f| * 2^-53 ~ 5e-9 -- the whole error being chased.  So two_prod it as well.
    """
    ph, pe = _two_prod(t, f)                    # t*f = ph + pe   exactly
    q = np.rint(ph / TWOPI_HI)
    qh, qe = _two_prod(q, TWOPI_HI)             # q*TWOPI_HI = qh + qe  exactly
    return (((ph - qh) - qe) - q * TWOPI_LO) + pe


def P_accurate(sym, t):
    """P_L(t) with accurate argument reduction (~1e-13 even at t ~ 1e7)."""
    t = np.asarray(t, dtype=float)
    acc = np.zeros(t.shape, dtype=float)
    for w, f in zip(sym.weights, sym.freqs):
        acc += w * np.cos(reduced_arg(t, f))
    return acc


def Psi_accurate(sym, t):
    return sym.H(t) - P_accurate(sym, t)


def float64_arg_error_bound(sym, t):
    """Rigorous bound on the plain-float64 evaluation error of P_L at |t| <= t."""
    return float(np.sum(np.abs(sym.weights) * np.abs(sym.freqs)) * abs(t) * 2.0 ** -53)
