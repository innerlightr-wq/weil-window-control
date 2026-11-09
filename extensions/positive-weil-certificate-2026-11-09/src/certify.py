"""Certified (Arb ball) evaluation of the Weil form on the EVEN window basis.

This is a SECOND IMPLEMENTATION of the repository's closed form, in a different
arithmetic backend (python-flint / Arb ball arithmetic instead of mpmath
floating point).  DISCLOSURE: it uses the SAME closed-form derivation as
`src/zeta_window.py` (notes/stage3_assembly.md); it is an independent
*implementation and arithmetic*, NOT an independent derivation.  It is used to
turn NUMERICAL values from the mpmath path into CERTIFIED enclosures.

Basis (EVEN sector of src/zeta_window.py):  f_k(x) = cos(k pi x / L), k >= 0,
on [-L, L], extended by zero.  int f_j f_k = L delta_jk for k >= 1 and 2L for k = 0.

Weil form:  M = Pole + Arch - LogPi - Prime, in the ORTHONORMAL basis.
"""
from fractions import Fraction
from flint import arb, acb, ctx
from wcore import von_mangoldt, prime_powers_below

def _A(x):
    if isinstance(x, Fraction):
        return arb(x.numerator) / arb(x.denominator)
    return arb(x)

class Win:
    def __init__(self, Lnum, Lden):
        self.Lf = Fraction(Lnum, Lden)
        self.L = _A(self.Lf)
        self.U = 2 * self.L
        self._cache = {}

    def freq(self, k):
        return arb(k) * arb.pi() / self.L

    def norm2(self, k):
        return 2 * self.L if k == 0 else self.L

    def F_at_half(self, k):
        """P_k = int f_k(x) e^{-x/2} dx.  EVEN: cos(wL)=(-1)^k, sin(wL)=0."""
        w = self.freq(k)
        q = w ** 2 + arb(1) / 4
        return arb(-1) ** k * (self.L / 2).sinh() / q

    def C_terms(self, j, k):
        """(c, alpha, beta) meaning c*sin(alpha u + beta); ('lin', gamma) meaning
        c*(2L-u)*cos(gamma u).  Mirrors zeta_window.C_terms for sector EVEN."""
        wj, wk = self.freq(j), self.freq(k)
        p1d, p2d = j - k, j + k           # exact integer tests, no float comparison
        out = []
        if p1d == 0:
            out.append((arb(1) / 2, 'lin', wk))
        else:
            p1 = wj - wk
            out.append((1 / (2 * p1), wk, p1 * self.L))
            out.append((-1 / (2 * p1), wj, -p1 * self.L))
        if p2d == 0:
            out.append((arb(1) / 2, 'lin', wk))
        else:
            p2 = wj + wk
            out.append((1 / (2 * p2), -wk, p2 * self.L))
            out.append((-1 / (2 * p2), wj, -p2 * self.L))
        return out

    def Csym(self, j, k):
        if j == k:
            return self.C_terms(j, k)
        return ([(c / 2, a, b) for c, a, b in self.C_terms(j, k)] +
                [(c / 2, a, b) for c, a, b in self.C_terms(k, j)])

    def eval_terms(self, terms, u):
        u = _A(u)
        tot = arb(0)
        for c, a, b in terms:
            tot += c * ((2 * self.L - u) * (b * u).cos() if a == 'lin'
                        else (a * u + b).sin())
        return tot

    def _sf(self, alpha, is_lin):
        """Cached special-function values keyed on the frequency only.

        The matrix has O(N) distinct frequencies but O(N^2) entries, so caching
        here turns O(N^2) Lerch/digamma evaluations into O(N).  Keyed on
        (precision, kind, frequency) so a precision change never reuses a value.
        """
        key = (ctx.prec, is_lin, alpha.str(40))
        hit = self._cache.get(key)
        if hit is not None:
            return hit
        U = self.U
        w = (-2 * U).exp()
        eul = arb.const_euler()
        if is_lin:
            z0 = acb(arb(1) / 4, -alpha / 2)
            val = ((z0.digamma() + eul).real,
                   acb.lerch_phi(acb(w), acb(2), z0), acb(2).zeta(z0))
        else:
            z = acb(arb(1) / 4, -alpha / 2)
            val = (z.digamma() + eul, acb.lerch_phi(acb(w), acb(1), z))
        self._cache[key] = val
        return val

    def arch(self, terms):
        """Archimedean term in closed form (digamma / Lerch Phi / Hurwitz zeta)."""
        U, eul = self.U, arb.const_euler()
        w = (-2 * U).exp()
        C0 = self.eval_terms(terms, 0)
        tot = arb(0)
        for c, alpha, beta in terms:
            if alpha == 'lin':
                dg, lp2, hz = self._sf(beta, True)
                tot += c * U * dg
                ph = acb(0, beta * U).exp() * (-U / 2).exp()
                tot -= (c / 2) * (ph * lp2 - hz).real
            else:
                dgg, lp1 = self._sf(alpha, False)
                tot += c * (acb(0, beta).exp() * dgg).imag
                tot += c * (-U / 2).exp() * (acb(0, beta + alpha * U).exp() * lp1).imag
        return -eul * C0 + tot

    def matrix(self, N):
        """Pole + Arch - LogPi - Prime on modes 0..N-1, orthonormalised. Arb balls."""
        pp = prime_powers_below(2 * self.Lf.numerator, self.Lf.denominator)
        P = [self.F_at_half(k) for k in range(N)]
        nrm = [self.norm2(k).sqrt() for k in range(N)]
        logpi = arb.pi().log()
        M = [[arb(0)] * N for _ in range(N)]
        for a in range(N):
            for b in range(a, N):
                tm = self.Csym(a, b)
                d = nrm[a] * nrm[b]
                v = 2 * P[a] * P[b] / d + self.arch(tm) / d
                if a == b:
                    v -= logpi * self.eval_terms(tm, 0) / d
                pr = arb(0)
                for n, p, k in pp:
                    pr += 2 * arb(p).log() / arb(n).sqrt() * self.eval_terms(tm, arb(n).log())
                v -= pr / d
                M[a][b] = M[b][a] = v
        return M

    def quad_form(self, M, x):
        """x^T M x / (x^T x) with x a list of Fractions -- a CERTIFIED enclosure."""
        N = len(x)
        xs = [_A(v) for v in x]
        num = arb(0)
        for i in range(N):
            for j in range(N):
                num += xs[i] * M[i][j] * xs[j]
            
        den = arb(0)
        for i in range(N):
            den += xs[i] * xs[i]
        return num / den
