"""weil_core.py -- parameterised reuse of the PUBLISHED certified functional.

The interval class, every proved transcendental, the canonical autocorrelation,
the pole / archimedean / log-pi / prime terms, the archimedean tail bound, the
gamma_E elimination and the exact quartet closed forms are taken UNCHANGED, by
import, from the published verifier

    extensions/planted-witness-closure-2026-10-08/verify_certificate.py
    sha256 c49858a1d50afbfefaaa63a4914c0c8bd2eb88269c21b4dd15e24b03deb510f0

which lives in this repository and is imported from there -- not copied, not
edited.  Its sha256 is checked at import time, so a fresh checkout either
imports the identical published file or fails loudly.  Nothing in the certified
path is reimplemented here.  What this module adds is only:

  (1) the vector X is a parameter instead of a module constant;
  (2) the x-INDEPENDENT part of the archimedean sum is precomputed once,
          SA_k = sum_{m<M} w_k (1+E_m)/d_{k,m},
          SB_k = sum_{m<M} [ s_m U/d_{k,m} - (1+E_m)(s_m^2-w_k^2)/d_{k,m}^2 ],
      d_{k,m} = s_m^2 + w_k^2, s_m = 2m+1/2,
      so that  sum_{m<M} laplace(A,B,s_m,E_m) = sum_k [A_k SA_k + B_k SB_k]
      exactly (laplace is linear in (A,B)).  This is an algebraic
      rearrangement of the published loop, not a new bound;
  (3) the resulting split
          Q_zeta(x)  in  Qcore(x) + tail(x),
      Qcore quadratic in x with interval coefficients, tail(x) the published
      bound  +-[ (3/4)C0(x) + (1/2)sup|C'|(x) ] / (M - 3/4).
      The ONLY x-dependent error is tail(x), given by that explicit formula.

Validated in tests/test_core_vs_published.py against the published run.
"""
import sys, os, json, hashlib
from fractions import Fraction as Q

HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.dirname(HERE)                       # extensions/<this extension>/
REPO = os.path.dirname(os.path.dirname(EXT))      # repository root

PUBLISHED_SHA = 'c49858a1d50afbfefaaa63a4914c0c8bd2eb88269c21b4dd15e24b03deb510f0'


def _sha(path):
    with open(path, 'rb') as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def _locate_published():
    """Find the published verifier relative to this file -- no absolute paths,
    so a fresh clone works unchanged -- and refuse to run on anything whose
    sha256 differs from the published one."""
    cands = [os.path.join(REPO, 'extensions', 'planted-witness-closure-2026-10-08'),
             os.path.join(EXT, 'sources', 'from_repo')]
    seen = []
    for d in cands:
        f = os.path.join(d, 'verify_certificate.py')
        if os.path.exists(f):
            h = _sha(f)
            if h == PUBLISHED_SHA:
                return d
            seen.append((f, h))
    raise RuntimeError(
        'the published verifier verify_certificate.py (sha256 %s) was not found '
        'in any of %r; candidates present but with a different hash: %r'
        % (PUBLISHED_SHA, cands, seen))


SRC = _locate_published()
sys.path.insert(0, SRC)
import verify_certificate as V            # the published module, unchanged

Iv, mk = V.Iv, V.mk
L, U, GAMMA, NDIM = V.L, V.U, V.GAMMA, V.NDIM
M, MEXP = V.M, V.MEXP
LOG2, LOGPI, PI = V.LOG2, V.LOGPI, V.PI
wk = V.wk
CertificationFailure = V.CertificationFailure
XFROZEN = list(V.X)


def _set_X(x):
    V.X = [Q(t) for t in x]


# ---------------------------------------------------------------- arch sums (cache)
_CACHE = {}
CACHE_PATH = os.path.join(EXT, 'results', 'arch_sums_cache.json')


def source_sha():
    return _sha(os.path.join(SRC, 'verify_certificate.py'))


def cache_key():
    """Cache key: truncation, interval grid, basis/parameters, source hash."""
    return dict(M=M, GRID_bits=V.GRID.bit_length() - 1, L=str(L), GAMMA=str(GAMMA),
                NDIM=NDIM, basis='EVEN_D cos((2k+1)pi u/(2L))/sqrt(L)',
                verifier_sha256=source_sha(),
                quantity='SA_k, SB_k : x-independent arch sums')


def arch_sums(progress=False):
    """SA_k, SB_k: the x-independent part of sum_{m<M} laplace(.,.,s_m,E_m)."""
    if 'SA' in _CACHE:
        return _CACHE['SA'], _CACHE['SB']
    if os.path.exists(CACHE_PATH):
        with open(CACHE_PATH) as fh:
            blob = json.load(fh)
        if blob.get('key') == cache_key():
            _CACHE['SA'] = [Iv(Q(t[0]), Q(t[1])) for t in blob['SA']]
            _CACHE['SB'] = [Iv(Q(t[0]), Q(t[1])) for t in blob['SB']]
            return _CACHE['SA'], _CACHE['SB']
        print('  arch_sums cache key MISMATCH: recomputing', file=sys.stderr)
    small = Iv(0, Q(1, 1 << 90))
    SA = [Iv(0) for _ in range(NDIM)]
    SB = [Iv(0) for _ in range(NDIM)]
    for m in range(M):
        s = Q(4*m + 1, 2)
        E = V.iv_exp(Iv(-s * U)) if m < 20 else small
        one_pE = Iv(1) + E
        sI = Iv(s)
        for k in range(NDIM):
            w = wk(k)
            d = sI.sq() + w.sq()
            SA[k] = SA[k] + w * one_pE / d
            SB[k] = SB[k] + (sI * Iv(U) / d - one_pE * (sI.sq() - w.sq()) / (d * d))
        if progress and m % 20000 == 0:
            print(f'    arch_sums m={m}/{M}', file=sys.stderr)
    _CACHE['SA'], _CACHE['SB'] = SA, SB
    with open(CACHE_PATH, 'w') as fh:
        json.dump(dict(key=cache_key(),
                       SA=[[str(v.lo), str(v.hi)] for v in SA],
                       SB=[[str(v.lo), str(v.hi)] for v in SB]), fh, indent=1)
    return SA, SB


# ---------------------------------------------------------------- the functional
def Qcore(x):
    """The tail-free part of Q_zeta(f_x): an exact quadratic form in x, enclosed."""
    _set_X(x)
    A, B = V.C_canonical()
    C0 = sum(Q(t) * Q(t) for t in x)
    SA, SB = arch_sums()
    T = Iv(Q(1, 2*(M + 1)) - Q(1, 6*M*M), Q(1, 2*M))
    acc = Iv(0)
    for k in range(NDIM):
        acc = acc + A[k] * SA[k] + B[k] * SB[k]
    arch_core = Iv(C0) * (Iv(17) * LOG2 - T) - Iv(2) * acc
    pole = V.pole_term()
    prime = V.prime_term(A, B)
    logpi = LOGPI * Iv(C0)
    return dict(core=pole + arch_core - logpi - prime, pole=pole, arch_core=arch_core,
                logpi=logpi, prime=prime, C0=C0, A=A, B=B,
                dCsup=V.dC_sup_bound(A, B))


def tail_of(x):
    """The published archimedean tail allowance, recomputed for THIS x."""
    d = Qcore(x)
    tailc = Iv(Q(3, 4)) * Iv(d['C0']) + Iv(Q(1, 2)) * d['dCsup']
    mag = tailc / Iv(Q(M) - Q(3, 4))
    return Iv(-mag.hi, mag.hi)


def Q_zeta(x):
    """Enclosure of Q_zeta(f_x) = Qcore(x) + tail(x)."""
    d = Qcore(x)
    tailc = Iv(Q(3, 4)) * Iv(d['C0']) + Iv(Q(1, 2)) * d['dCsup']
    mag = tailc / Iv(Q(M) - Q(3, 4))
    tail = Iv(-mag.hi, mag.hi)
    d['tail'] = tail
    d['Q'] = d['core'] + tail
    return d


def quartet(x, delta):
    _set_X(x)
    return V.quartet(delta)


def Q_add(x, delta):
    """Enclosure of Q_add,delta(f_x) = Q_zeta(f_x) + 4[A_delta^2 - B_delta^2]."""
    d = Q_zeta(x)
    Ad, Bd = quartet(x, delta)
    return d['Q'] + Iv(4) * Ad.sq() - Iv(4) * Bd.sq(), Ad, Bd, d


def F_at_gamma(x):
    _set_X(x); return V.F_at_gamma()


# ---------------------------------------------------------------- G by polarization
def G_core():
    """The 4x4 interval matrix of the TAIL-FREE quadratic form Qcore.

    Qcore is exactly a quadratic form in x (every term is), so polarization is
    exact up to interval widths:
        G[j][j] = Qcore(e_j),
        G[j][l] = (Qcore(e_j+e_l) - Qcore(e_j) - Qcore(e_l))/2 .
    G_core is NOT Q_zeta's matrix: Q_zeta(x) = x^T G_core x + tail(x), and
    tail(x) is NOT a quadratic form.  Any statement about Q_zeta from G_core
    must carry tail(x) explicitly.
    """
    if 'G' in _CACHE:
        return _CACHE['G']
    e = lambda j: [Q(1) if i == j else Q(0) for i in range(NDIM)]
    diag = [Qcore(e(j))['core'] for j in range(NDIM)]
    G = [[None]*NDIM for _ in range(NDIM)]
    for j in range(NDIM):
        G[j][j] = diag[j]
    for j in range(NDIM):
        for l in range(j+1, NDIM):
            v = [Q(1) if i in (j, l) else Q(0) for i in range(NDIM)]
            off = (Qcore(v)['core'] - diag[j] - diag[l]) * Iv(Q(1, 2))
            G[j][l] = G[l][j] = off
    _CACHE['G'] = G
    return G


def ab_vectors(delta):
    """a_k(delta), b_k(delta) with  A_delta(f_x) = a.x,  B_delta(f_x) = b.x."""
    e = lambda j: [Q(1) if i == j else Q(0) for i in range(NDIM)]
    a = []; b = []
    for j in range(NDIM):
        A_, B_ = quartet(e(j), delta)
        a.append(A_); b.append(B_)
    return a, b


def quad(Gm, x):
    tot = Iv(0)
    for j in range(NDIM):
        for l in range(NDIM):
            tot = tot + Gm[j][l] * Iv(Q(x[j]) * Q(x[l]))
    return tot


def dot(v, x):
    tot = Iv(0)
    for j in range(NDIM):
        tot = tot + v[j] * Iv(Q(x[j]))
    return tot
