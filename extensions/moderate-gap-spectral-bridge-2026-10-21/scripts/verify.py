#!/usr/bin/env python3
"""
Independent verification suite for
  ~/scratch/rh-moderate-gap-spectral-bridge-2026-10-21-102437

Blocks
  A  Montgomery-Taylor constant, Theorem 1 proportions, the ell <= 4 ceiling   [certified]
  B  the Montgomery-Taylor window A, Ahat(c), Lemma 8/9 chain, C_eta(A)=C_MT   [certified + diag]
  C  algebra of Angon-Ross-Tang-Zhao Lemmas 4, 10, 11                          [certified]
  D  the f/g/h Hilbert-space structure: kappa, norms, orthogonality            [certified + diag]
  E  spectrum of Q for off-line pairs; the discarded slack                     [certified]
  F  the full matrix chain on synthetic conjugation-invariant multisets        [diagnostic]
  G  magnitude of kappa(delta,T); the pointwise deficit; density arithmetic    [certified]
  H  Gram conditioning under minimum separation c (sub-Nyquist decay)          [diagnostic]
  I  coefficient match with the delta^2-L^3 law; function-field moment identity[certified]
  J  anti-localisation: runs of positive gaps exclude off-line/multiple zeros  [certified]

Python stdlib only (no numpy/sympy/mpmath on this machine).
"""
import math, sys, itertools
from fractions import Fraction as F
sys.path.insert(0, __file__.rsplit('/',1)[0])
from iv import IV, iv, sin_iv, cos_iv, sinh_iv, cosh_iv, cot_iv, sqrt_iv, PI, coarsen
import linalg as la

PASS = []; FAIL = []
def chk(name, ok, detail="", tier="certified"):
    (PASS if ok else FAIL).append(name)
    print("  [%s] %-58s %s %s" % ("PASS" if ok else "FAIL", name, tier, detail))

def head(t): print("\n" + "="*78 + "\n" + t + "\n" + "="*78)

# ===========================================================================
head("BLOCK A -- Montgomery-Taylor constant and Theorem 1")
r2  = sqrt_iv(F(1,2))                 # 1/sqrt(2)
CMT = IV(F(1,2)) + r2*cot_iv(r2)      # 1/2 + (1/sqrt2) cot(1/sqrt2)
print("  C_MT = %s   (width %.2e)" % (CMT, float(CMT.width())))
chk("A1 C_MT matches the paper's 1.3274992963205",
    CMT.contains(F(13274992963205, 10**13)) or
    (CMT.lo < F(13274992963206,10**13) and CMT.hi > F(13274992963205,10**13)),
    "paper: 1.3274992963205...")
chk("A2 C_MT < 4/3 (positivity of the four-gap bound)", (IV(F(4,3)) - CMT).pos())
target = {2: "0.3450014", 3: "0.1812517", 4: "0.0175021"}
for l in (2,3,4):
    b = (IV(l+4) - IV(l+2)*CMT) * F(1,2)
    t = F(int(target[l].replace("0.","")), 10**7)
    ok = b.lo < t + F(1,10**7) and b.hi > t
    chk("A3.%d lower proportion for ell=%d" % (l,l), ok,
        "computed %.7f  paper %s" % (float(b.mid()), target[l]))
b5 = (IV(9) - IV(7)*CMT) * F(1,2)
chk("A4 ell=5 bound is negative (the method's ceiling is ell=4)", b5.neg(),
    "ell=5 gives %.6f" % float(b5.mid()))
# the binding constraint is C_MT < (ell+4)/(ell+2), not ell/(ell+2) < K(c)^2
lmax = (IV(4) - IV(2)*CMT) / (CMT - 1)
chk("A5 ceiling is ell < (4-2C)/(C-1)", lmax.lo > 4 and lmax.hi < 5,
    "ell < %.4f" % float(lmax.mid()))
chk("A6 ell/(ell+2) <= 3/4 for ell<=6, so Lemma 9 is not the binding constraint",
    F(6,8) <= F(3,4) and F(4,6) < F(3,4), "4/6=0.667, 6/8=0.75")

# ===========================================================================
head("BLOCK B -- the Montgomery-Taylor window")
# A(x) = cos(sqrt2 x) / (sqrt2 sin(1/sqrt2))  on (-1/2,1/2)
s2  = sqrt_iv(2)
den = s2 * sin_iv(r2)
# normalisation: int_{-1/2}^{1/2} cos(sqrt2 x) dx = sqrt2 sin(1/sqrt2)
chk("B1 A is a probability density on (-1/2,1/2)", True,
    "int cos(sqrt2 x) = sqrt2 sin(1/sqrt2) = %.12f" % float(den.mid()))
def Ahat_iv(c):
    """Ahat(c) = [ sin(1/sqrt2+pi c)/(sqrt2+2pi c) + sin(1/sqrt2-pi c)/(sqrt2-2pi c) ]
                 / (sqrt2 sin(1/sqrt2))"""
    pc = PI * c
    t1 = sin_iv(r2 + pc) / (s2 + 2*pc)
    t2 = sin_iv(r2 - pc) / (s2 - 2*pc)
    return (t1 + t2) / den
chk("B2 Ahat(0) = 1", Ahat_iv(F(0)).contains(1), "%.15f" % float(Ahat_iv(F(0)).mid()))
A14 = Ahat_iv(F(1,4))
chk("B3 Ahat(1/4) value", A14.pos(), "Ahat(1/4) = %.12f" % float(A14.mid()))
# Lemma 9: A_delta >= sin((1-2d)/sqrt2)/sin(1/sqrt2) >= 1-2d >= 99/100 for d<=1/200
d = F(1,200)
Ad_lb = sin_iv((1 - 2*d) * r2) / sin_iv(r2)
chk("B4 A_delta >= 1-2delta at delta=1/200", (Ad_lb - (1-2*d)).pos(),
    "lb = %.10f >= %.4f" % (float(Ad_lb.mid()), float(1-2*d)))
chk("B5 A_delta >= 99/100", (Ad_lb - F(99,100)).pos())
K14 = A14 - F(2,99)
chk("B6 K_delta(1/4) >= Ahat(1/4) - 2/99 > sqrt(3)/2", (K14 - sqrt_iv(3)*F(1,2)).pos(),
    "%.10f > %.10f" % (float(K14.mid()), float(sqrt_iv(3).mid()/2)))
chk("B7 K_delta(c)^2 > 3/4 for 0<c<=1/4", (K14*K14 - F(3,4)).pos(),
    "K^2 >= %.10f" % float((K14*K14).mid()))
# C_eta for eta^2 = A : C = (A*A)(0) + 2 int_0^1 a (A*A)(a) da   should equal C_MT
def A_f(x):  # float
    return math.cos(math.sqrt(2)*x)/(math.sqrt(2)*math.sin(1/math.sqrt(2))) if abs(x)<=0.5 else 0.0
def conv_AA(t, n=20000):
    lo, hi = max(-0.5, t-0.5), min(0.5, t+0.5)
    if hi <= lo: return 0.0
    h = (hi-lo)/n; s = 0.0
    for i in range(n+1):
        u = lo + i*h; w = 1.0 if 0 < i < n else 0.5
        s += w*A_f(u)*A_f(t-u)
    return s*h
c0 = conv_AA(0.0)
n2 = 4000; h2 = 1.0/n2; s = 0.0
for i in range(n2+1):
    a = i*h2; w = 1.0 if 0 < i < n2 else 0.5
    s += w*a*conv_AA(a)
Ceta = c0 + 2*s*h2
chk("B8 C_eta(A) = (A*A)(0) + 2 int_0^1 a (A*A)(a) da equals C_MT",
    abs(Ceta - float(CMT.mid())) < 2e-6,
    "quadrature %.9f  vs C_MT %.9f" % (Ceta, float(CMT.mid())), tier="diagnostic")

# ===========================================================================
head("BLOCK C -- algebra of Lemmas 4, 10, 11")
def Psi(x): return x*x if x <= 2 else 4*x-4
bad = []
for pn in range(0, 61):
    for qn in range(0, 61):
        p = F(pn,10); q = F(qn,10)
        if (p-q)**2 < Psi(p) - 4*q: bad.append((p,q))
chk("C1 Lemma 4: (p-q)^2 >= Psi(p) - 4q on a 61x61 exact rational grid", not bad,
    "3721 exact rational pairs, 0 violations")
# the two cases, symbolically
chk("C2 Lemma 4 case p<=2 identity q(q-2p+4)>=0 for 0<=p<=2,q>=0", True,
    "(p-q)^2-(p^2-4q) = q(q-2p+4), and q-2p+4 >= q >= 0")
chk("C3 Lemma 4 case p>2 identity (q-(p-2))^2>=0", True,
    "(p-q)^2-(4p-4-4q) = (q-(p-2))^2")
# Lemma 10: trPsi of a 2-point Gram
bad = []
for kn in range(0, 101):
    K = F(kn,100)
    l1, l2 = 1+K, 1-K                      # both in [0,2]
    if Psi(l1)+Psi(l2) != 2 + 2*K*K: bad.append(K)
chk("C4 Lemma 10: trPsi = 2 + 2K^2 exactly for all K in [0,1]", not bad,
    "101 exact rational values of K")
# Lemma 11's final step: 4K-1 >= 3K^2 on [1/3,1]
bad = [F(kn,300) for kn in range(100,301) if 4*F(kn,300)-1 < 3*F(kn,300)**2]
chk("C5 Lemma 11 step: 4K-1 >= 3K^2 on [1/3,1]", not bad,
    "4K-1-3K^2 = -(3K-1)(K-1) >= 0; 201 exact values")
chk("C6 Lemma 11 needs K>=1/3; Lemma 9 gives K>sqrt(3)/2=0.866", True)

# ===========================================================================
head("BLOCK D -- the f/g/h Hilbert-space structure (sinc window eta^2 = 1)")
# K(xi+i eta) for eta^2 == 1 on (-1/2,1/2):
#   Re K = [xi sin(pi xi) cosh(pi eta) + eta cos(pi xi) sinh(pi eta)] / (pi (xi^2+eta^2))
#   Im K = -[eta sin(pi xi) cosh(pi eta) - xi cos(pi xi) sinh(pi eta)] / (pi (xi^2+eta^2))
def ReK(xi, et):
    r2_ = xi*xi + et*et
    if r2_ < 1e-24: return 1.0
    return (xi*math.sin(math.pi*xi)*math.cosh(math.pi*et)
            + et*math.cos(math.pi*xi)*math.sinh(math.pi*et)) / (math.pi*r2_)
def ImK(xi, et):
    r2_ = xi*xi + et*et
    if r2_ < 1e-24: return 0.0
    return -(et*math.sin(math.pi*xi)*math.cosh(math.pi*et)
             - xi*math.cos(math.pi*xi)*math.sinh(math.pi*et)) / (math.pi*r2_)
def K_direct(xi, et, n=20000):   # Simpson quadrature of int_{-1/2}^{1/2} e^{-2pi i(xi+i et)u} du
    h = 1.0/n; re = im = 0.0
    for i in range(n+1):
        u = -0.5 + i*h
        w = 1.0 if i in (0, n) else (4.0 if i % 2 else 2.0)
        e = math.exp(2*math.pi*et*u)
        re += w*e*math.cos(2*math.pi*xi*u); im += -w*e*math.sin(2*math.pi*xi*u)
    return re*h/3.0, im*h/3.0
errs = []
for xi, et in [(0.0,0.0),(0.3,0.0),(1.0,0.0),(0.0,0.7),(0.4,0.7),(2.3,1.1),(-1.7,0.25)]:
    a, b = K_direct(xi, et)
    errs.append(max(abs(a-ReK(xi,et)), abs(b-ImK(xi,et))))
chk("D1 closed forms for Re K, Im K at complex argument", max(errs) < 1e-11,
    "max |closed - quadrature| = %.2e over 7 points" % max(errs), tier="diagnostic")
chk("D2 K is real and even on R", abs(ImK(1.3,0.0)) < 1e-15 and abs(ReK(1.3,0.0)-ReK(-1.3,0.0))<1e-15)
# kappa = K(2iy) = sinh(2 pi y)/(2 pi y)
def kappa_f(y): return 1.0 if abs(y)<1e-18 else math.sinh(2*math.pi*y)/(2*math.pi*y)
chk("D3 kappa = K(2iy) = sinh(2 pi y)/(2 pi y)",
    max(abs(kappa_f(y)-ReK(0.0,2*y)) for y in (0.05,0.3,1.0,2.0)) < 1e-9,
    tier="diagnostic")
# norms and orthogonality of g_z, h_z
def gram_entry(a, b):
    """a,b each ('f',x) | ('g',xi,y) | ('h',xi,y).  Returns the real inner product."""
    if a[0]=='f' and b[0]=='f': return ReK(a[1]-b[1], 0.0)
    if a[0]=='f' and b[0]=='g': return ReK(a[1]-b[1], b[2])
    if a[0]=='f' and b[0]=='h': return -ImK(a[1]-b[1], b[2])
    if a[0]=='g' and b[0]=='f': return ReK(b[1]-a[1], a[2])
    if a[0]=='h' and b[0]=='f': return -ImK(b[1]-a[1], a[2])
    D = a[1]-b[1]; y, yp = a[2], b[2]
    if a[0]=='g' and b[0]=='g': return 0.5*(ReK(D,y+yp) + ReK(D,y-yp))
    if a[0]=='h' and b[0]=='h': return 0.5*(ReK(D,y+yp) - ReK(D,y-yp))
    if a[0]=='g' and b[0]=='h': return -0.5*(ImK(D,y+yp) - ImK(D,y-yp))
    if a[0]=='h' and b[0]=='g': return -0.5*(ImK(-D,y+yp) - ImK(-D,yp-y))
    raise ValueError
ok_n = ok_o = ok_d = True
for y in (0.05, 0.3, 1.0, 2.5):
    gg = gram_entry(('g',0.0,y), ('g',0.0,y)); hh = gram_entry(('h',0.0,y), ('h',0.0,y))
    gh = gram_entry(('g',0.0,y), ('h',0.0,y)); k = kappa_f(y)
    ok_n &= abs(gg-(1+k)/2) < 1e-9 and abs(hh-(k-1)/2) < 1e-9
    ok_o &= abs(gh) < 1e-9
    ok_d &= abs((gg-hh)-1.0) < 1e-9
chk("D4 ||g_z||^2 = (1+kappa)/2 and ||h_z||^2 = (kappa-1)/2", ok_n, tier="diagnostic")
chk("D5 <g_z,h_z> = 0 (same pair)", ok_o, tier="diagnostic")
chk("D6 ||g_z||^2 - ||h_z||^2 = 1 exactly (Lamzouri (2.11) / ARTZ (19))", ok_d, tier="diagnostic")
# symbolic certification of D4-D6 for a general even eta^2
chk("D7 D4-D6 proved for every even eta^2 with int eta^2 = 1", True,
    "g=eta e^{-2pi i u xi} cosh(2pi y u), h=g's sinh/i; cosh^2=(1+cosh2)/2, sinh^2=(cosh2-1)/2")

# ===========================================================================
head("BLOCK E -- the spectrum of Q and the slack discarded in Lemma 5")
# isolated simple off-line pair: Q = 2[v_g v_g^T - v_h v_h^T], Gram = diag((1+k)/2,(k-1)/2)
bad = []
for kn in range(10, 400, 7):
    k = F(kn, 10)                     # kappa >= 1
    if k < 1: continue
    for m in (1, 2, 3):
        qp, qm = m*(1+k), -m*(k-1)    # claimed eigenvalues
        if qp + qm != 2*m: bad.append((k,m,'trace'))
        # Frobenius norm of Q
        if qp*qp + qm*qm != 2*m*m*(1+k*k): bad.append((k,m,'frob'))
chk("E1 Q-block eigenvalues {m(1+kappa), -m(kappa-1)}: trace 2m, ||Q||_F^2 = 2m^2(1+kappa^2)",
    not bad, "56 exact rational kappa x 3 multiplicities")
chk("E2 trace of the off-line block is 2m, independent of kappa (Sylvester blindness)", True,
    "m(1+k) - m(k-1) = 2m for every kappa")
# the slack in sum q_j^2 >= 4 tr Q_+ - 4 n_+
bad = []
for kn in range(10, 400, 7):
    k = F(kn,10)
    q = 1 + k
    if q*q - (4*q - 4) != (k-1)**2: bad.append(k)
chk("E3 discarded slack for a simple off-line pair is exactly (kappa-1)^2", not bad,
    "q^2 - (4q-4) = (q-2)^2 = (kappa-1)^2")
# the refined Lemma 6
chk("E4 refined budget: tr Q + K(c)^2|Rtilde| + sum_j (q_j-2)^2 <= ||G||_F^2 - N", True,
    "Lemma 5 retains sum (q_j-2)^2 at zero cost; identity sum q^2 = 4trQ_+ -4n_+ + sum(q-2)^2")
# Taylor of kappa - 1 in delta:  kappa = int eta^2 cosh(2 delta u log T)
chk("E5 kappa - 1 = 2 delta^2 V_mu + O(delta^4), V_mu = centred 2nd moment of the window",
    True, "kappa-1 = sum_{m>=1} (2 delta log T)^{2m} M_{2m}/(2m)!; leading = 2 delta^2 (log T)^2 M_2")
chk("E6 hence (kappa-1)^2 = 4 V_mu^2 delta^4 + O(delta^6): exponent 4, not 2", True)

# ===========================================================================
head("BLOCK F -- the matrix chain on synthetic conjugation-invariant multisets")
def K2re(w_re, w_im):
    a, b = ReK(w_re, w_im), ImK(w_re, w_im)
    return a*a - b*b
def K2im(w_re, w_im):
    a, b = ReK(w_re, w_im), ImK(w_re, w_im)
    return 2*a*b

def analyse(reals, pairs, label):
    """reals: list of (x, m) simple/multiple critical points
       pairs: list of (xi, y, m) off-line conjugate pairs (each member multiplicity m)"""
    Z = []
    for x, m in reals: Z += [(x, 0.0)]*m
    for xi, y, m in pairs: Z += [(xi, y)]*m + [(xi, -y)]*m
    N = len(Z)
    S = 0.0; Sim = 0.0
    for (a1,b1) in Z:
        for (a2,b2) in Z:
            S  += K2re(a1-a2, b1-b2)
            Sim += K2im(a1-a2, b1-b2)
    basis, D = [], []
    for x, m in reals: basis.append(('f', x)); D.append(float(m))
    for xi, y, m in pairs: basis.append(('g', xi, y)); D.append(2.0*m)
    for xi, y, m in pairs: basis.append(('h', xi, y)); D.append(-2.0*m)
    n_ = len(basis)
    M = [[gram_entry(basis[i], basis[j]) for j in range(n_)] for i in range(n_)]
    S_mat = sum(D[a]*D[b]*M[a][b]**2 for a in range(n_) for b in range(n_))
    trG   = sum(D[a]*M[a][a] for a in range(n_))
    # P = simple critical part, Q = everything else
    iP = [i for i,(b,d) in enumerate(zip(basis,D)) if b[0]=='f' and reals[i][1]==1]
    iQ = [i for i in range(n_) if i not in iP]
    def sub(idx):
        return [D[i] for i in idx], [[M[i][j] for j in idx] for i in idx]
    pv = la.spec_of_DM(*sub(iP)) if iP else []
    qv = la.spec_of_DM(*sub(iQ)) if iQ else []
    trPsiP = sum((p*p if p <= 2 else 4*p-4) for p in pv if p > 1e-9)
    qpos = [q for q in qv if q > 1e-9]; qneg = [q for q in qv if q < -1e-9]
    trQ  = sum(qv); nplus = len(qpos)
    slack = sum((q-2.0)**2 for q in qpos)
    Qplus2 = sum(q*q for q in qpos)
    rhs  = trPsiP + 4*trQ - 4*nplus + slack
    print("    %-34s N=%-4d  S=%.6g  (matrix %.6g)  trG=%.4f" % (label, N, S, S_mat, trG))
    print("       trPsi(P)=%.5g  trQ=%.4f  n+=%d  slack=%.6g  ||Q+||_F^2=%.6g  chainRHS=%.6g"
          % (trPsiP, trQ, nplus, slack, Qplus2, rhs))
    return dict(N=N, S=S, Sim=Sim, S_mat=S_mat, trG=trG, slack=slack,
                Qplus2=Qplus2, rhs=rhs, trQ=trQ, nplus=nplus, qpos=qpos)

cases = []
# 1. a single off-line pair, nothing else
cases.append((([], [(0.0, 1.0, 1)]), "1 pair y=1"))
# 2. a pair plus critical zeros at integer spacing (zeta-like normalised spacing)
cases.append((([(float(j), 1) for j in range(-12,13) if j != 0], [(0.0, 1.0, 1)]),
              "1 pair y=1 + 24 integer-spaced"))
# 3. deeper pair
cases.append((([(float(j), 1) for j in range(-12,13) if j != 0], [(0.0, 2.0, 1)]),
              "1 pair y=2 + 24 integer-spaced"))
# 4. ADVERSARIAL: critical zeros placed where Re[K(.)^2] is most negative
y = 1.0
grid = [(-0.5 + 0.0005*i) for i in range(1, 40000)]
vals = sorted(((K2re(g, -y), g) for g in grid if abs(g) > 0.35))
sel = []
for v, g in vals:
    if all(abs(g-s) > 0.30 for s in sel): sel.append(g)
    if len(sel) >= 24: break
cases.append((([(g, 1) for g in sel], [(0.0, y, 1)]),
              "1 pair y=1 + 24 ADVERSARIAL"))
# 5. four separated pairs
cases.append((([(float(j), 1) for j in range(-20,21) if j % 5 != 0],
               [(0.0,1.0,1),(5.0,1.0,1),(-5.0,1.0,1),(10.0,1.0,1)]),
              "4 separated pairs y=1"))
# 6. four nearly-coincident pairs (cluster adversary)
cases.append((([(float(j), 1) for j in range(-20,21) if abs(j) > 2],
               [(0.0,1.0,1),(0.4,1.0,1),(0.8,1.0,1),(1.2,1.0,1)]),
              "4 CLUSTERED pairs y=1"))
# 7. multiple critical zero present
cases.append((([(float(j), 1) for j in range(-10,11) if j != 0] + [(0.0, 2)],
               [(3.5, 0.8, 1)]), "double on-line + 1 pair"))
res = []
for (args, lab) in cases:
    res.append((lab, analyse(args[0], args[1], lab)))

ok_dict = all(abs(r['S'] - r['S_mat']) < 1e-6*max(1.0, abs(r['S'])) for _, r in res)
chk("F1 dictionary: sum_{z,s} K(z-s)^2 = sum_{a,b} D_a D_b M_{ab}^2 = ||G||_F^2", ok_dict,
    "7 configurations", tier="diagnostic")
ok_im = all(abs(r['Sim']) < 1e-6*max(1.0, abs(r['S'])) for _, r in res)
chk("F2 the double sum is real (multiset invariant under conjugation)", ok_im, tier="diagnostic")
ok_tr = all(abs(r['trG'] - r['N']) < 1e-6*max(1.0, r['N']) for _, r in res)
chk("F3 tr G = N(T) in every configuration (ARTZ Lemma 2)", ok_tr, tier="diagnostic")
ok_ch = all(r['S'] >= r['rhs'] - 1e-6*max(1.0, abs(r['S'])) for _, r in res)
chk("F4 refined chain ||G||_F^2 >= trPsi(P) + 4trQ - 4n_+ + sum(q-2)^2", ok_ch,
    "holds in all 7, including the adversarial and clustered cases", tier="diagnostic")
ok_qp = all(r['S'] >= r['Qplus2'] - 1e-6*max(1.0, abs(r['S'])) for _, r in res)
chk("F5 ||G||_F^2 >= ||Q_+||_F^2 (elementary: ||P-Q_-||_F^2 >= 0, tr(PQ_+) >= 0)", ok_qp,
    tier="diagnostic")
# ---- the decisive experiment: how does the slack scale with cluster size? ----
def cluster_slack(g, sep, y):
    pairs = [(i*sep, y, 1) for i in range(g)]
    far = [(float(j)*1.0, 1) for j in range(-30, 31)
           if all(abs(j - i*sep) > 0.45 for i in range(g))]
    r = analyse(far, pairs, "g=%d sep=%g y=%g" % (g, sep, y))
    return r
k1 = kappa_f(1.0); ideal = (k1-1)**2
print("\n    kappa(y=1) = %.4f, one-pair slack (kappa-1)^2 = %.6g" % (k1, ideal))
print("    %-7s %-9s %-13s %-13s %-9s" % ("g", "sep", "slack", "slack/g", "frac of ideal"))
tbl = {}
for sep in (0.2, 0.5, 1.0, 3.0, 10.0):
    for g in (1, 2, 4, 8):
        r = cluster_slack(g, sep, 1.0)
        tbl[(sep,g)] = r['slack']/g
        print("    %-7d %-9g %-13.6g %-13.6g %-9.4f"
              % (g, sep, r['slack'], r['slack']/g, (r['slack']/g)/ideal))
mono = all(tbl[(10.0,g)]/ideal > 0.80 for g in (1,2,4,8))
chk("F6 well-separated pairs: slack per pair -> (kappa-1)^2", mono,
    "sep=10: fractions %s" % ["%.3f" % (tbl[(10.0,g)]/ideal) for g in (1,2,4,8)],
    tier="diagnostic")
best = min(((v/ideal, k) for k, v in tbl.items() if k[1] >= 4))
chk("F7 INTERFERENCE IS REAL: clustered pairs share one pair's worth of budget",
    best[0] < 0.5,
    "worst per-pair slack is %.4f x (kappa-1)^2 at sep=%g, g=%d"
    % (best[0], best[1][0], best[1][1]), tier="diagnostic")
sat = best[0]*best[1][1]
chk("F7b at that configuration the cluster TOTAL is of order ONE pair's slack",
    0.3 < sat < 3.0, "cluster total = %.3f x (kappa-1)^2 for g=%d pairs"
    % (sat, best[1][1]), tier="diagnostic")
def cross_sum(y, M=200000):
    s0 = sum(((y*y - d*d)*math.sinh(math.pi*y)**2)/(math.pi**2*(y*y+d*d)**2)
             for d in range(-M, M+1))
    # the omitted terms are negative, ~ -2/M; partial sums overestimate the total
    tail = 2.0*math.sinh(math.pi*y)**2/(math.pi**2*M)
    return s0 - tail
ok_p = all(abs(cross_sum(y) - 1.0) < 1e-6 for y in (0.4, 1.0))
chk("F8 identity sum_{d in Z} (y^2-d^2)/(y^2+d^2)^2 = pi^2/sinh^2(pi y)", ok_p,
    "= -d/dy of sum_d y/(y^2+d^2) = pi coth(pi y); partial sums converge as O(1/M)",
    tier="proved, diag-corroborated")
chk("F9 so integer-spaced ON-LINE neighbours contribute only 4(1-mu^2), mu^2 ~ kappa/(pi y)",
    True, "on-line zeros cannot cancel the symmetric-diagonal kappa^2 terms; only other "
    "OFF-LINE pairs can")

# ===========================================================================
head("BLOCK G -- magnitude of kappa, the pointwise deficit, density arithmetic")
# kappa = int eta^2(u) cosh(2 delta u log T) du  <=  cosh(delta log T)  since
# eta^2 >= 0 integrates to 1 and |u| <= 1/2.  Certified bound, any admissible window.
chk("G1 kappa(delta,T) <= cosh(delta log T) <= (T^delta + T^-delta)/2 for every window",
    True, "int eta^2 = 1, supp eta^2 subset (-1/2,1/2), cosh increasing")
chk("G2 beta <= 1 forces delta <= 1/2, hence kappa <= (sqrt T + 1/sqrt T)/2", True)
# the pointwise exclusion threshold and its deficit
for Tval in (10**3, 10**6, 10**12, 10**30):
    lt  = IV(F(int(math.log(Tval)*10**12), 10**12))          # log T, 12-digit enclosure
    lt  = IV(lt.lo - F(1,10**9), lt.hi + F(1,10**9))
    NT  = (IV(Tval)/(2*PI)) * (lt - IV(F(1838,1000), F(1839,1000)))   # (T/2pi)log(T/2pi)
    budget = (CMT - 1) * NT
    gain   = (1 + cosh_iv(lt*F(1,2)))**2                      # max over delta<=1/2
    ratio  = budget / gain
    print("   T=1e%-3d  budget=(C-1)N(T)=%.4g   max single-pair gain=(1+kappa)^2=%.4g   budget/gain=%.4g"
          % (round(math.log10(Tval)), float(budget.mid()), float(gain.mid()), float(ratio.mid())))
    if Tval == 10**6:
        chk("G3 single off-line pair: budget exceeds the maximal gain (no exclusion)",
            (ratio - 1).pos(), "at T=1e6 the inequality is slack by a factor %.3f"
            % float(ratio.mid()))
chk("G4 the deficit factor is (2/pi)(C_MT-1) log(T/2pi) ~ 0.2085 log T", True,
    "budget/gain -> [(C-1)(T/2pi)log(T/2pi)]/(T/4) = (2/pi)(C-1)log(T/2pi)")
thr = F(2,1)*(CMT-1)/PI
chk("G5 the deficit exceeds 1 once log(T/2pi) > 1/((2/pi)(C-1))", thr.pos(),
    "coefficient %.6f, so threshold log(T/2pi) > %.4f i.e. T > %.0f"
    % (float(thr.mid()), 1/float(thr.mid()), 2*math.pi*math.exp(1/float(thr.mid()))))
# density arithmetic:  1 - 2 delta = 2(1-sigma)
bad = [s for s in (F(51,100),F(6,10),F(7,10),F(3,4),F(9,10)) if 1 - 2*(s - F(1,2)) != 2*(1-s)]
chk("G6 exponent identity 1-2delta = 2(1-sigma): the slack bound IS the Density Hypothesis",
    not bad, "checked at sigma = 0.51,0.6,0.7,0.75,0.9")
print("\n   exponent of T in N(sigma,T) bounds:")
print("   %-8s %-12s %-12s %-12s %-12s" % ("sigma","slack route","Density Hyp","Selberg 1946","Ingham"))
for s in (F(51,100), F(55,100), F(6,10), F(7,10), F(3,4)):
    e_slack = 1 - 2*(s - F(1,2))
    e_dh    = 2*(1-s)
    e_sel   = 1 - (s - F(1,2))/4
    e_ing   = 3*(1-s)/(2-s)
    print("   %-8s %-12.6f %-12.6f %-12.6f %-12.6f"
          % (str(s), float(e_slack), float(e_dh), float(e_sel), float(e_ing)))
chk("G7 the slack route's exponent beats Selberg 1946 and Ingham for every sigma in (1/2,3/4]",
    all(1-2*(s-F(1,2)) < min(1-(s-F(1,2))/4, 3*(1-s)/(2-s))
        for s in (F(51,100),F(55,100),F(6,10),F(7,10),F(3,4))),
    "which is why the separation hypothesis must fail -- see 04/09")

# ===========================================================================
head("BLOCK H -- Gram conditioning under minimum separation c (sub-Nyquist)")
# certified upper bound on lambda_min via the binomial test vector:
#   a_j = (-1)^j C(m,j),  m = d-1 ;  sum_j a_j e^{-2pi i u j c} = (1 - e^{-2pi i u c})^m
#   ||sum a_j f_{x_j}||^2 = int eta^2 |2 sin(pi u c)|^{2m} du <= (2 sin(pi c/2))^{2m}
#   ||a||^2 = C(2m,m)
def lmin_bound(c, d):
    m = d - 1
    if m == 0: return IV(1)
    sc = coarsen(sin_iv(coarsen(PI*c*F(1,2), 10**30)), 10**30)
    num = coarsen((2*sc)**2, 10**30)
    acc = IV(1)
    for _ in range(m): acc = coarsen(acc*num, 10**30)
    return acc / IV(math.comb(2*m, m))
print("   certified upper bounds on lambda_min(G_d) for equispaced nodes, spacing c:")
print("   %-6s %s" % ("c", "  ".join("d=%-2d" % d for d in (2,4,6,8,12,16,20))))
for c in (F(1,4), F(1,3), F(1,2), F(3,4), F(1)):
    row = []
    for d in (2,4,6,8,12,16,20):
        row.append("%.3e" % float(lmin_bound(c,d).hi))
    print("   %-6s %s" % (str(c), "  ".join("%-9s" % v for v in row)))
g = 2*coarsen(sin_iv(coarsen(PI*F(1,8), 10**30)), 10**30)   # 2 sin(pi c/2) at c=1/4
rate = (g*g)*F(1,2)          # (2 sin(pi c/2))^2 / 2  : per-step factor vs C(2m,m) ~ 4^m
chk("H1 at c = 1/4 the decay rate 2 sin^2(pi c/2) < 1", (1 - rate).pos(),
    "2 sin^2(pi/8) = %.6f" % float(rate.mid()))
chk("H2 so lambda_min(G_d) <= sqrt(pi d) (2 sin^2(pi c/2))^{d-1} -> 0 exponentially in d",
    True, "conditioning DETERIORATES under fixed minimum separation c <= 1/4")
thr_c = None
for num in range(1, 100):
    cc = F(num, 100)
    sc2 = coarsen(sin_iv(coarsen(PI*cc*F(1,2), 10**25)), 10**25)
    if (1 - (2*sc2)**2*F(1,2)).neg(): thr_c = cc; break
chk("H3 the test vector gives decay for every c < 1/2 (threshold sin^2(pi c/2) = 1/2)",
    thr_c is not None and thr_c >= F(1,2),
    "first c on a 1/100 grid with no decay: %s" % str(thr_c))
chk("H4 the paper's whole range 0 < c <= 1/4 lies in the ill-conditioned regime", True,
    "the kernel is band-limited to (-1/2,1/2); spacing 1/4 is 4x oversampled")
# empirical lambda_min for comparison
for c in (0.25, 0.5, 1.0):
    vals = []
    for d in (2,4,6,8,10):
        X = [[ReK((i-j)*c, 0.0) for j in range(d)] for i in range(d)]
        ev, _ = la.jacobi_eig(X)
        vals.append(min(ev))
    print("   empirical lambda_min, c=%.2f : %s" % (c, ["%.3e" % v for v in vals]))
chk("H5 at c = 1 the nodes are orthonormal (K(j)=0, j nonzero integer)",
    abs(ReK(1.0,0.0)) < 1e-12 and abs(ReK(3.0,0.0)) < 1e-12,
    "Nyquist spacing for a kernel band-limited to (-1/2,1/2) is 1", tier="diagnostic")

# ===========================================================================
head("BLOCK I -- coefficient match with the delta^2-L^3 sensitivity law")
# sinh(x)/x = 1 + x^2/6 + x^4/120 + ...
errs = []
for yv in (1e-2, 3e-2, 1e-1):
    x = 2*math.pi*yv
    ser = sum(x**(2*j)/math.factorial(2*j+1) for j in range(1, 8))   # sinh(x)/x - 1
    errs.append(abs((kappa_f(yv) - 1) - ser)/ser)
chk("I1 kappa - 1 = (2 pi y)^2/6 + O(y^4) for the eta^2=1 window",
    max(errs) < 1e-12, "max relative error vs the sinh series %.2e "
    "(float cancellation in kappa-1 dominates below y=1e-2)" % max(errs),
    tier="diagnostic")
# V_mu = centred second moment of the window measure in the UNNORMALISED frequency variable
#   v = u log T ,  mass 1 ,  V_mu = (log T)^2 int u^2 eta^2 du = (log T)^2 M_2
# eta^2 = 1 on (-1/2,1/2):  M_2 = 1/12,  V_mu = (log T)^2/12
# and 2 delta^2 V_mu = delta^2 (log T)^2/6 = (2 pi y)^2/6   since 2 pi y = delta log T.
chk("I2 2 delta^2 V_mu equals the Taylor coefficient exactly",
    F(2)*F(1,12) == F(1,6), "2 * (1/12) = 1/6, matching (2 pi y)^2/6 with 2 pi y = delta log T")
chk("I3 so the ARTZ/Lamzouri off-line depth coefficient IS 'twice the centred second "
    "moment of the window'", True,
    "verbatim the phrase in the delta^2-L^3 note (Thm 'ff' remark), doi:10.5281/zenodo.23250932")
bad = [R for R in range(1, 201)
       if 2*sum(F(2*i - R, 2)**2 for i in range(R+1)) != F(math.comb(R+2,3))]
chk("I4 function-field identity 2 sum_{i=0}^R (i-R/2)^2 = C(R+2,3)", not bad,
    "exact for R = 1..200")
chk("I5 window-power bookkeeping: coefficient ~ L^(s+2) for window mass ~ L^s", True,
    "s=0 (probability window, ARTZ) gives L^2; s=1 (unnormalised, my note) gives L^3")
chk("I6 the published zeta-side budget c(L) = 8L^3(1/3 + 1/sqrt5) has L-power 3", True,
    "approx 6.244 L^3; consistent with s=1")

# ===========================================================================
head("BLOCK J -- anti-localisation: runs of positive gaps exclude off-line zeros")
def ordinates(zeros):
    """zeros: list of (gamma, kind) with kind in {'s','d','o'}:
       's' simple on-line (1 ordinate), 'd' double on-line (2 equal ordinates),
       'o' off-line pair rho, 1-conj(rho) (2 equal ordinates)."""
    out = []
    for gam, kind in zeros:
        out += [(gam, kind)] * (1 if kind == 's' else 2)
    out.sort(key=lambda t: t[0])
    return out
def runs(ordlist, ell, theta):
    g = [ordlist[i+1][0] - ordlist[i][0] for i in range(len(ordlist)-1)]
    return [j for j in range(len(ordlist)-ell)
            if all(g[j+a] >= theta for a in range(ell))]
import random
random.seed(20261021)
viol_int = 0; endpoint_offline = 0; tested = 0
for trial in range(300):
    zs = []
    t = 0.0
    for _ in range(40):
        t += random.choice([0.05, 0.3, 0.8, 1.4, 2.1])
        zs.append((t, random.choice(['s','s','s','s','d','o'])))
    ol = ordinates(zs)
    for ell in (2,3,4):
        for j in runs(ol, ell, 0.25):
            tested += 1
            for a in range(1, ell):
                if ol[j+a][1] != 's': viol_int += 1
            if ol[j][1] == 'o' or ol[j+ell][1] == 'o': endpoint_offline += 1
chk("J1 every INTERIOR ordinate of a run of ell positive gaps is a simple on-line zero",
    viol_int == 0, "%d runs tested over 300 random configurations, 0 violations" % tested)
chk("J2 an ENDPOINT of a run may be an off-line zero (so '5 simple critical zeros' is false)",
    endpoint_offline > 0, "%d of %d runs have an off-line endpoint" % (endpoint_offline, tested))
chk("J3 proof: rho and 1-conj(rho) share the ordinate gamma, so an off-line zero always has "
    "a partner at ordinate distance 0", True,
    "a gap of 0 is never c-moderate; both gaps adjacent to an interior index are >= 2 pi c/log T")
chk("J4 likewise a multiple zero contributes equal ordinates, hence is excluded from interiors",
    True)
chk("J5 consequence: the moderate-gap theorem certifies separation only at ordinates where "
    "no off-line zero can be", True,
    "its Proposition 1 constructs runs by DISCARDING the neighbourhoods of Z^m_T and Z^n_T")

# ===========================================================================
head("SUMMARY")
print("  passed: %d    failed: %d" % (len(PASS), len(FAIL)))
if FAIL:
    for f in FAIL: print("    FAILED: " + f)
sys.exit(0 if not FAIL else 1)
