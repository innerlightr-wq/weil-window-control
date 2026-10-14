"""ROUTE O, step 2: certify  M(delta) >= 0  on cells below the onset,
i.e. CERTIFIED NON-DETECTION on the four-mode subspace.

For a cell [d1, d2], w = d2-d1, write A_d = a(delta).x, A_1 = a(d1).x (same for b).
Since a(delta) - a(d1) = int_{d1}^{delta} a'(s) ds,
   |A_d - A_1| <= w * Ahat * ||x||,     Ahat := sup_{cell} ||a'(s)|| ,
so  |A_d^2 - A_1^2| <= 2 w Ahat ||x|| |A_1| + w^2 Ahat^2 ||x||^2
                    <= tau A_1^2 + w^2 Ahat^2 (1 + 1/tau) ||x||^2      (AM-GM)
for any tau > 0.  The AM-GM step is what makes this work: the loss is charged
against the term 4 A_1^2 that is ALREADY in M, instead of against ||a|| ||a'||,
and ||a|| = 0.88 while the component of a along the near-null direction is
~4e-5.  Hence, for every delta in the cell and every x,

   x^T M(delta) x  >=  x^T N x  -  Loss ||x||^2 ,
   N    := Gt3_mid + 4(1-tau) a1 a1^T - 4(1+tau') b1 b1^T          (rational)
   Loss := eps + 4 w^2 [ Ahat^2 (1+1/tau) + Bhat^2 (1+1/tau') ] + tiny ,
so  lambda_min(N) >= Loss   (certified by exact rational LDL)  implies
M(delta) >= 0 throughout the cell.

Rigorous derivative bounds, elementary:
   |a'_k(d)| = |int phi_k(u) cos(gamma u) u sinh(d u) du| <= (2/3) L^{3/2} sinh(dL),
       using |sinh(d u)| <= |u| sinh(dL)/L for |u| <= L (sinh(t)/t increasing)
       and int |phi_k| u^2 du <= (2/3) L^{5/2};
   |b'_k(d)| = |int phi_k(u) sin(gamma u) u cosh(d u) du| <= L^{3/2} cosh(dL),
       using int |phi_k| |u| du <= L^{3/2}.
With four components, Ahat <= 2 (2/3) L^{3/2} sinh(d2 L), Bhat <= 2 L^{3/2} cosh(d2 L).
"""
import sys, os, json, resource, time
sys.path.insert(0, os.path.dirname(__file__))
resource.setrlimit(resource.RLIMIT_AS, (4 << 30, 4 << 30))
from fractions import Fraction as Q
import weil_core as C
import ratlin as R
V = C.V; Iv = C.Iv; NDIM = C.NDIM; L = C.L

blob = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'results',
                                   'route_o_step1.json')))
Gt = [[Q(v) for v in row] for row in blob['Gt3_mid']]
EPS = Q(blob['eps'][1])
M0 = Q(blob['m0_lo'])
mid = lambda iv: (iv.lo + iv.hi) / 2
rad = lambda iv: (iv.hi - iv.lo) / 2

def up(iv): return iv.hi
def sq_up(x): return Q(x) * Q(x)

L32 = None
def Lpow32_up():
    """an upper rational bound for L^{3/2}"""
    global L32
    if L32 is None:
        s = V.iv_sqrt(Iv(L))
        L32 = (Iv(L) * s).hi
    return L32

def derivative_bounds(d2):
    """Ahat, Bhat : rational upper bounds on sup_cell ||a'||, ||b'||."""
    sh = V.iv_sinh(Iv(Q(d2) * L)).hi
    ch = V.iv_cosh(Iv(Q(d2) * L)).hi
    l32 = Lpow32_up()
    Ahat = Q(2) * Q(2, 3) * l32 * sh
    Bhat = Q(2) * l32 * ch
    return Ahat, Bhat

def ab_rational(d1):
    a, b = C.ab_vectors(Q(d1))
    return ([mid(t) for t in a], [mid(t) for t in b],
            max(rad(t) for t in a), max(rad(t) for t in b))

def try_cell(d1, w, taus):
    """Returns (ok, best_info). Certifies M(delta)>=0 on [d1, d1+w]."""
    d2 = d1 + w
    a1, b1, ra, rb = ab_rational(d1)
    Ahat, Bhat = derivative_bounds(d2)
    na = sum(t*t for t in a1); nb = sum(t*t for t in b1)
    # the (negligible) allowance for the 1e-70 widths of a1, b1 at rational d1
    tiny = 8 * (ra * Q(2) + rb * Q(2)) + Q(1, 10**60)
    best = None
    for tau in taus:
        for taup in taus:
            N = [[Gt[j][l] + 4*(1-tau)*a1[j]*a1[l] - 4*(1+taup)*b1[j]*b1[l]
                  for l in range(NDIM)] for j in range(NDIM)]
            Loss = EPS + 4*w*w*(Ahat*Ahat*(1 + 1/tau) + Bhat*Bhat*(1 + 1/taup)) + tiny
            B = [[N[j][l] - (Loss if j == l else 0) for l in range(NDIM)]
                 for j in range(NDIM)]
            piv = R.cholesky_pd(B)
            if piv is not None and all(p > 0 for p in piv):
                if best is None:
                    best = dict(tau=str(tau), taup=str(taup), Loss=str(Loss),
                                pivots=[str(p) for p in piv],
                                Ahat=str(Ahat), Bhat=str(Bhat))
                return True, best
    return False, None

TAUS = [Q(1, 10**k) for k in range(1, 9)] + [Q(1, 2)]
RIGHT = Q(13, 200)
BUDGET_S = float(os.environ.get('ROUTEO_BUDGET', '1500'))

cells = []
d1 = Q(0)
w = Q(1, 1000)
t0 = time.time()
while d1 < RIGHT and time.time() - t0 < BUDGET_S:
    w = min(w, RIGHT - d1)
    ok, info = try_cell(d1, w, TAUS)
    shr = 0
    while not ok and shr < 30:
        w = w / 2
        ok, info = try_cell(d1, w, TAUS)
        shr += 1
    if not ok:
        print(f'  STOP at {float(d1):.8f}: cannot certify even with w = {float(w):.3e}')
        break
    # store ONLY what verify_route_o.py consumes.  Loss and N are recomputed
    # there from (lo, hi, w, tau, tau'), so storing the exact-rational Loss adds
    # 50+ MB of digits and no verification value.
    cells.append(dict(lo=str(d1), hi=str(d1 + w), w=str(w), shrinks=shr,
                      tau=info['tau'], taup=info['taup']))
    if len(cells) % 25 == 1 or d1 + w >= RIGHT:
        print(f'  [{float(d1):.8f}, {float(d1+w):.8f}] w={float(w):.3e} '
              f'tau={info["tau"]} tau\'={info["taup"]} Loss={float(Q(info["Loss"])):.3e} '
              f'shrinks={shr}  n={len(cells)}', flush=True)
    d1 = d1 + w
    if shr == 0:
        w = w * 2
json.dump(dict(right=str(RIGHT), reached=str(d1), n_cells=len(cells),
               note=("Loss and N are NOT stored: verify_route_o.py recomputes both "
                     "from (lo, hi, w, tau, taup) and from route_o_step1.json, so "
                     "nothing here is trusted."),
               cells=cells),
          open(os.path.join(os.path.dirname(__file__), '..', 'results',
                            'route_o_step2.json'), 'w'), indent=1)
print(f'\n{len(cells)} cells; reached delta = {float(d1):.8f} of {float(RIGHT):.6f} '
      f'in {time.time()-t0:.0f}s')
print(f'complete = {d1 >= RIGHT}')
print(f'peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024:.1f} MiB')
