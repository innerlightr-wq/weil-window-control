import sys, json; sys.path.insert(0,'src')
import mpmath as mp
mp.mp.dps = 50
g1 = mp.im(mp.zetazero(1))
print('  gamma_1 =', mp.nstr(g1, 10), ' (paper states 14.1347)')

def K(L, g):
    L, g = mp.mpf(L), mp.mpf(g)
    return L**3/3 - (L**2*mp.sin(2*g*L)/(2*g) + L*mp.cos(2*g*L)/(2*g**2) - mp.sin(2*g*L)/(4*g**3))
def Kquad(L, g):   # independent check: direct quadrature of int u^2 sin^2(gu)
    L, g = mp.mpf(L), mp.mpf(g)
    return mp.quad(lambda u: u**2*mp.sin(g*u)**2, [-L, L])

print('\n=== (a) K closed form vs direct quadrature (eq. 5) ===')
for Lx in ('0.6','0.8','1.0','1.3','1.6','2.0'):
    a, b = K(Lx, g1), Kquad(Lx, g1)
    print('  L=%-4s K=%-20s quad=%-20s rel.diff=%.2e' % (Lx, mp.nstr(a,12), mp.nstr(b,12), float(abs(a-b)/a)))

# measured C(L): the paper's own values, cross-checked against the projection-formula run
PAPER_C = {'0.8':'0.3108','1.0':'0.8530','1.3':'2.4987','1.6':'4.9218','2.0':'10.083'}
K2 = {r['L']: r for r in json.load(open('extensions/projection-formula/data/stage2b_k2.json'))}
print('\n=== (b) section 4.3 table ===')
print('  %-5s %-10s %-10s %-12s %-12s %-10s %-10s %-10s %-10s %-10s' %
      ('L','C paper','C rerun','4K paper','4K computed','C/4K pap','C/4K comp','C/(8L3/3)','C/(4L3/3)','ok'))
PAP_4K = {'0.8':'0.7419','1.0':'1.3427','1.3':'3.1158','1.6':'5.1130','2.0':'10.652'}
PAP_R1 = {'0.8':'0.419','1.0':'0.635','1.3':'0.802','1.6':'0.963','2.0':'0.947'}
PAP_R2 = {'0.8':'0.228','1.0':'0.320','1.3':'0.427','1.6':'0.451','2.0':'0.473'}
PAP_R3 = {'0.8':'0.455','1.0':'0.640','1.3':'0.853','1.6':'0.901','2.0':'0.945'}
bad = []
for Lx in ('0.8','1.0','1.3','1.6','2.0'):
    L = mp.mpf(Lx); fk = 4*K(Lx, g1)
    Cp = mp.mpf(PAPER_C[Lx]); Cr = mp.mpf(K2[Lx]['C'])
    r1 = Cr/fk; r2 = Cr/(8*L**3/3); r3 = Cr/(4*L**3/3)
    chk = (abs(fk-mp.mpf(PAP_4K[Lx]))/fk < 1e-3 and abs(r1-mp.mpf(PAP_R1[Lx])) < 1e-3
           and abs(r2-mp.mpf(PAP_R2[Lx])) < 1e-3 and abs(r3-mp.mpf(PAP_R3[Lx])) < 1e-3
           and abs(Cr-Cp)/Cp < 1e-3)
    if not chk: bad.append(Lx)
    print('  %-5s %-10s %-10s %-12s %-12s %-12s %-10s %-10s %-10s %-10s' %
          (Lx, PAPER_C[Lx], mp.nstr(Cr,6), PAP_4K[Lx], mp.nstr(fk,6),
           PAP_R1[Lx], mp.nstr(r1,4), mp.nstr(r2,4), mp.nstr(r3,4), 'OK' if chk else 'MISMATCH'))
print('  4.3 table:', 'all rows verified' if not bad else 'MISMATCH at %s' % bad)

print('\n=== (c) the sin^2 deviation list |K/(L^3/3) - 1| ===')
PAP_DEV = {'0.6':'17.3','0.8':'8.68','1.0':'0.699','1.3':'6.37','1.6':'6.38','2.0':'0.136'}
bad=[]
for Lx in ('0.6','0.8','1.0','1.3','1.6','2.0'):
    L=mp.mpf(Lx); dev = abs(K(Lx,g1)/(L**3/3) - 1)*100
    ok = abs(dev-mp.mpf(PAP_DEV[Lx])) <= mp.mpf('0.01')*max(mp.mpf(1),mp.mpf(PAP_DEV[Lx]))
    if not ok: bad.append((Lx, mp.nstr(dev,6), PAP_DEV[Lx]))
    print('  L=%-4s computed %-10s%%  paper %-8s%%  %s' % (Lx, mp.nstr(dev,4), PAP_DEV[Lx], 'OK' if ok else 'MISMATCH'))
print('  sin^2 list:', 'all verified' if not bad else 'MISMATCH %s' % bad)

print('\n=== (d) dictionary ratios, y^2=x^3+x over F_5, a=1e-6 ===')
from curves import squarefree_over_Fq, count_points_hyperelliptic, newton_p_to_A, A_to_power_sums
from weilform import toeplitz, eigsym
PAP = {8:'1.406',16:'1.195',32:'1.096',64:'1.047',128:'1.024'}
q=5; a=mp.mpf('1e-6')
bad=[]
for R in (8,16,32,64,128):
    # normalised spectrum = single off-circle real pair {rho, rho^-1}, rho = e^a (Thm 6 family)
    t=[2*mp.cosh(n*a) for n in range(R+1)]
    lam=eigsym(toeplitz(t,R))[0][0]
    L=mp.mpf(R)*mp.log(q)/2; delta=a/mp.log(q)
    pred=-mp.mpf(4)/3*L**3*delta**2/mp.log(q)
    ratio=lam/pred
    closed=(1+mp.mpf(1)/R)*(1+mp.mpf(2)/R)
    ok=abs(ratio-mp.mpf(PAP[R]))<mp.mpf('0.001')
    if not ok: bad.append((R,mp.nstr(ratio,6),PAP[R]))
    print('  R=%-4d ratio %-10s paper %-8s (1+1/R)(1+2/R)=%-10s %s'
          % (R, mp.nstr(ratio,6), PAP[R], mp.nstr(closed,6), 'OK' if ok else 'MISMATCH'))
print('  dictionary:', 'all verified' if not bad else 'MISMATCH %s' % bad)

print('\n=== (e) Theorem 6: binom(R+2,3) and the R=1 exact value ===')
for R in (1,2,4,8,16,32):
    t=[2*mp.cosh(n*a) for n in range(R+1)]
    lam=eigsym(toeplitz(t,R))[0][0]
    b=mp.binomial(R+2,3); mom=2*sum((mp.mpf(i)-mp.mpf(R)/2)**2 for i in range(R+1))
    print('  R=%-4d -lam/a^2 = %-18s binom(R+2,3) = %-10s 2*centred moment = %-10s'
          % (R, mp.nstr(-lam/a**2,12), mp.nstr(b,8), mp.nstr(mom,8)))
r1=-4*mp.sinh(a/2)**2; t=[2*mp.cosh(n*a) for n in range(2)]
print('  R=1 exact -4sinh^2(a/2) = %s  vs eig = %s' % (mp.nstr(r1,14), mp.nstr(eigsym(toeplitz(t,1))[0][0],14)))

print('\n=== (f) c(L) = 8L^3(1/3 + 1/sqrt5) ~ 6.244 L^3 ===')
c=8*(mp.mpf(1)/3+1/mp.sqrt(5)); print('  coefficient 8(1/3+1/sqrt5) =', mp.nstr(c,8), ' paper states 6.244')
print('\n=== (g) K(L,gamma) -> L^3/3 as gamma -> inf ; <= 2L^3/3 always ===')
for Lx in ('0.8','2.0'):
    L=mp.mpf(Lx)
    print('  L=%-4s K(L,1e6)=%-14s L^3/3=%-14s  2L^3/3=%s' %
          (Lx, mp.nstr(K(Lx,mp.mpf('1e6')),8), mp.nstr(L**3/3,8), mp.nstr(2*L**3/3,8)))
