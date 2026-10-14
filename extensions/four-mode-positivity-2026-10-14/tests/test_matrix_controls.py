"""D6, D7: software controls for the generic matrix statements.
These are NOT RH results. Exact rationals throughout."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))
from fractions import Fraction as Q
import ratlin as R
fails = []
def chk(c, m):
    print(('PASS  ' if c else 'FAIL  ') + m)
    if not c: fails.append(m)

def theta_and_verdict(H, b):
    """Theta = 4 b^T H^-1 b by exact solve; and the direct inertia of M = H - 4bb^T."""
    y = R.solve(H, b)
    z = R.dot(b, y) if y is not None else None
    n = len(b)
    M = [[Q(H[i][j]) - 4*Q(b[i])*Q(b[j]) for j in range(n)] for i in range(n)]
    d = R.cholesky_pd(M)
    cs = R.charpoly(M)
    # M has a negative eigenvalue iff p(t)=det(tI-M) has a root < 0 iff
    # sign pattern of c[0]=det(-M)... use LDL pivots when available
    neg = None
    if d is not None:
        neg = any(t < 0 for t in d)
    return z, neg, M, cs

# ---- Theta < 1, = 1, > 1 on exact rationals -------------------------------
H = [[Q(2),Q(0)],[Q(0),Q(1)]]
for lab, b, want in (('Theta<1', [Q(1,4),Q(1,4)], False),
                     ('Theta>1', [Q(1),Q(1)],    True)):
    z, neg, M, cs = theta_and_verdict(H, b)
    chk(z is not None and ((4*z > 1) == want) and (neg == want),
        f'D6 {lab}: Theta = {float(4*z):.6f}, M has a negative direction = {neg} (expected {want})')
# exact Theta = 1 : choose b = (t,0) with 4 t^2/2 = 1 -> t^2 = 1/2, irrational.
# use H = [[8,0],[0,1]], b = (1,0): Theta = 4*1/8 = 1/2. Instead scale: H=[[4,0],[0,1]], b=(1,0) -> 4/4=1
H1 = [[Q(4),Q(0)],[Q(0),Q(1)]]; b1 = [Q(1),Q(0)]
z, neg, M, cs = theta_and_verdict(H1, b1)
chk(4*z == 1, f'D6 Theta = 1 exactly: 4z = {4*z}')
chk(R.cholesky_pd(M) is None, 'D6 Theta = 1: M is singular (LDL hits a zero pivot)')
xs = R.solve(H1, b1)
chk(all(v == 0 for v in R.matvec(M, xs)), 'D6 Theta = 1: x_* = H^-1 b spans a kernel direction (Cor 3)')
chk(R.quadform(M, xs) == 0, 'D6 Theta = 1: x_*^T M x_* = z - 4z^2 = 0')
# Corollary 2 identity on a random-ish exact example
H2 = [[Q(3),Q(1),Q(0)],[Q(1),Q(5),Q(2)],[Q(0),Q(2),Q(7)]]
b2 = [Q(1,3),Q(-2,5),Q(4,7)]
y2 = R.solve(H2, b2); z2 = R.dot(b2, y2)
M2 = [[H2[i][j] - 4*b2[i]*b2[j] for j in range(3)] for i in range(3)]
chk(R.quadform(M2, y2) == z2 - 4*z2*z2, f'D7 Cor 2: x_*^T M x_* = z - 4z^2 exactly (z = {z2})')
chk(R.dot(b2, y2) == z2 and R.quadform(H2, y2) == z2, 'D7 Cor 2: b.x_* = x_*^T H x_* = z')
# b = 0
z3, neg3, M3, _ = theta_and_verdict(H2, [Q(0)]*3)
chk(z3 == 0 and neg3 is False, 'D6 b = 0: Theta = 0 and M = H has no negative direction')
# nearly singular positive H
eps = Q(1, 10**9)
H4 = [[Q(1),Q(0)],[Q(0),eps]]
b4 = [Q(0), eps]                  # Theta = 4 eps
z4, neg4, _, _ = theta_and_verdict(H4, b4)
chk(4*z4 == 4*eps and neg4 is False,
    f'D6 nearly singular positive H: Theta = {float(4*z4):.3e} < 1, no negative direction')
b5 = [Q(0), Q(1)]                 # Theta = 4/eps >> 1
z5, neg5, _, _ = theta_and_verdict(H4, b5)
chk(4*z5 > 1 and neg5 is True, f'D6 nearly singular positive H, large Theta = {float(4*z5):.3e}: detection')
# outside the hypothesis: H indefinite -> Proposition 1 unavailable
H6 = [[Q(1),Q(0)],[Q(0),Q(-1)]]
b6 = [Q(0), Q(1)]
y6 = R.solve(H6, b6); z6 = R.dot(b6, y6)
M6 = [[H6[i][j] - 4*b6[i]*b6[j] for j in range(2)] for i in range(2)]
chk(4*z6 < 1 and R.cholesky_pd(M6) is not None and any(t < 0 for t in R.cholesky_pd(M6)),
    f'D6 H indefinite: 4z = {float(4*z6)} <= 1 yet M HAS a negative direction '
    f'-> Proposition 1 is genuinely unavailable without [H1]')
# congruence invariance
S = [[Q(1),Q(2),Q(0)],[Q(0),Q(1),Q(-1)],[Q(3),Q(0),Q(1)]]
def cong(A):
    n = 3
    return [[sum(S[k][i]*Q(A[k][l])*S[l][j] for k in range(n) for l in range(n))
             for j in range(n)] for i in range(n)]
Hc = cong(H2); bc = [sum(S[k][i]*b2[k] for k in range(3)) for i in range(3)]
yc = R.solve(Hc, bc); zc = R.dot(bc, yc)
chk(zc == z2, f'D6 congruence invariance: z is unchanged under x = S y (z = {z2})')
# D7 residual enclosure
m = min(R.cholesky_pd(H2))        # a valid lower bound on the spectrum? no: pivots != eigenvalues
# use a genuine m: H2 - m I must be PSD; find the largest such m by bisection on LDL pivots
lo, hi = Q(0), Q(10)
for _ in range(60):
    mid = (lo+hi)/2
    Hm = [[H2[i][j] - (mid if i==j else 0) for j in range(3)] for i in range(3)]
    d = R.cholesky_pd(Hm)
    if d is not None and all(t > 0 for t in d): lo = mid
    else: hi = mid
mval = lo
for y in ([Q(0)]*3, [Q(1,10),Q(0),Q(1,10)], y2):
    r = [b2[i] - R.matvec(H2, y)[i] for i in range(3)]
    s = 2*R.dot(b2, y) - R.quadform(H2, y)
    rn2 = R.dot(r, r)
    chk(s <= z2 <= s + rn2/mval,
        f'D7 residual enclosure s <= z <= s + ||r||^2/m : [{float(s):.9f}, '
        f'{float(s + rn2/mval):.9f}] contains z = {float(z2):.9f}  (m = {float(mval):.6f})')

print()
print('MATRIX CONTROLS:', 'ALL PASS' if not fails else f'{len(fails)} FAILURES')
sys.exit(1 if fails else 0)
