"""INDEPENDENT cross-check of a NEW witness, using the repository's
independently coded analytic route.

What is independent: the archimedean term is evaluated by a COMPLETELY
DIFFERENT closed form -- digamma/trigamma at the four complex points
z_k = 1/4 - i w_k/2, with the Euler constant cancelling identically and a
g-series/midpoint error budget -- instead of 131071 Laplace integrals with a
1/M (or the new 1/M^2) tail.  The module is not imported by, and does not read
any output of, the primary verifier.

What is SHARED, stated plainly: the frozen mathematical input (x, L, gamma);
the classical explicit formula; the canonical (A_k,B_k) reduction of the
autocorrelation (recoded there, not imported); and the idea of exact rational
interval arithmetic (re-implemented there).  This is a second implementation,
not external peer review.

Its error budget (~2.5e-5, set by JG and the trigamma midpoint bound) is four
orders LOOSER than the sharp tail of tail_sharp.py, so it can confirm the new
witnesses only where the margin exceeds that -- i.e. at the right part of the
new interval, NOT at the thin left endpoint.  That limitation is reported, not
hidden.
"""
import sys, os, json, resource
HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.dirname(HERE)
REPO = os.path.dirname(os.path.dirname(EXT))
IP_SHA = 'bb72086ecb6e41c7450ac2e685149942fdff9785652de4744e81f5bc164ef4a8'


def _locate():
    """Resolve the repository's independent implementation relative to this
    file -- no absolute paths -- and refuse anything whose sha256 differs."""
    import hashlib
    cands = [os.path.join(REPO, 'extensions', 'planted-witness-consolidation-2026-10-08'),
             os.path.join(EXT, 'sources', 'from_repo')]
    for d in cands:
        f = os.path.join(d, 'independent_point_check.py')
        if os.path.exists(f):
            with open(f, 'rb') as fh:
                if hashlib.sha256(fh.read()).hexdigest() == IP_SHA:
                    return d
    raise RuntimeError('independent_point_check.py (sha256 %s) not found in %r'
                       % (IP_SHA, cands))


sys.path.insert(0, _locate())
resource.setrlimit(resource.RLIMIT_AS, (4 << 30, 4 << 30))
from fractions import Fraction as Q
import independent_point_check as IP

blob = json.load(open(os.path.join(EXT, 'results', 'certificate_new.json')))
cells = blob['cells']

def check(x, delta, label):
    IP.X = [Q(t) for t in x]
    A, B = IP.canonical()
    C0 = sum((Q(t)*Q(t) for t in IP.X), Q(0))
    c0 = IP.R(0)
    for k in range(4):
        c0 = c0 + B[k]*IP.R(IP.U)
    assert c0.a <= C0 <= c0.b, 'C(0) != sum x_k^2 (recoded reduction disagrees)'
    pole, prime, logpi, _ = IP.pole_prime_logpi(A, B)
    ar = IP.arch_via_digamma(A, B)
    Qz = pole + ar - logpi - prime
    Ad, Bd = IP.quartet_complex_sinc(Q(delta))
    Qa = Qz + IP.R(4)*Ad.sq() - IP.R(4)*Bd.sq()
    print(f'{label}: delta = {delta}')
    print(f'   Q_zeta  = [{float(Qz.a):.12f}, {float(Qz.b):.12f}]  width {float(Qz.b-Qz.a):.3e}')
    print(f'   Q_add   = [{float(Qa.a):.12f}, {float(Qa.b):.12f}]  upper {float(Qa.b):.9e}')
    print(f'   NEGATIVE UPPER BOUND: {Qa.b < 0}')
    return dict(label=label, delta=str(delta), x=[str(t) for t in x],
                Q_zeta=[str(Qz.a), str(Qz.b)], Q_add=[str(Qa.a), str(Qa.b)],
                negative=bool(Qa.b < 0),
                truncations=dict(JREC=IP.JREC, JG=IP.JG, NG=IP.NG, MEXP=IP.MEXP, gridbits=IP.GB))

out = []
# the witness of the last cell (the one adjacent to delta = 1/4), at two points
xlast = cells[-1]['x']
for dd in (Q(1, 5), Q(1, 4)):
    if Q(cells[-1]['lo']) <= dd <= Q(cells[-1]['hi']) or dd == Q(1, 5):
        pass
    out.append(check(xlast, dd, 'last-cell witness'))
# the witness of the FIRST cell, at the right end where its margin is large
out.append(check(cells[0]['x'], Q(1, 4), 'first-cell witness'))
# ---------------------------------------------------------------------------
# The comparison that makes this a cross-check rather than two monologues:
# the two implementations must enclose the SAME numbers.  The primary is
# imported only here, after the independent values are already fixed above.
sys.path.insert(0, HERE)
import weil_core as PC, tail3 as PT3            # primary route

print()
print('CROSS-CHECK: the two enclosures of the same real numbers must overlap,')
print('and the (much sharper) primary one must sit inside the independent one.')
agree = True
for rec in out:
    x = [Q(t) for t in rec['x']]; d = Q(rec['delta'])
    ilo, ihi = Q(rec['Q_zeta'][0]), Q(rec['Q_zeta'][1])
    pz = PT3.Q_zeta3(x)['Q']
    inside = ilo <= pz.lo and pz.hi <= ihi
    alo, ahi = Q(rec['Q_add'][0]), Q(rec['Q_add'][1])
    pa, _, _, _ = PT3.Q_add3(x, d)
    ov = alo <= pa.hi and pa.lo <= ahi
    both_neg = bool(pa.hi < 0 and ahi < 0)
    rec['primary_Q_zeta'] = [str(pz.lo), str(pz.hi)]
    rec['primary_Q_add'] = [str(pa.lo), str(pa.hi)]
    rec['primary_Qzeta_inside_independent'] = bool(inside)
    rec['Q_add_overlap'] = bool(ov)
    rec['both_negative'] = both_neg
    agree = agree and inside and ov and both_neg
    print(f"  {rec['label']}, delta = {rec['delta']}:")
    print(f"     Q_zeta  primary [{float(pz.lo):.12e}, {float(pz.hi):.12e}]")
    print(f"             indep   [{float(ilo):.12e}, {float(ihi):.12e}]   primary inside: {inside}")
    print(f"     Q_add   primary [{float(pa.lo):.9e}, {float(pa.hi):.9e}]")
    print(f"             indep   [{float(alo):.9e}, {float(ahi):.9e}]   overlap: {ov}   both < 0: {both_neg}")
print()
print('CROSS-CHECK AGREES:' if agree else 'CROSS-CHECK DISAGREES:', agree)
print('Limitation, stated: the independent budget (~1.1e-04 on Q_zeta here) is')
print('about 90x the worst margin of the new certificate (1.28e-06 at the left')
print('endpoint), so it corroborates the RIGHT part of [9/100, 1/4] only.  It')
print('does NOT reach delta = 9/100, and is not claimed to.')
json.dump(out, open(os.path.join(EXT, 'results', 'indep_crosscheck.json'), 'w'), indent=1)
print()
print(f'peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024:.1f} MiB')
sys.exit(0 if agree else 1)
