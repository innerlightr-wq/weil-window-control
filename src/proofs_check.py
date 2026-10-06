#!/usr/bin/env python3
r"""Task 3: numerical verification of every step of the two proofs."""
import sys, os, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from zeta_window import build_matrix, F, norm2, EVEN_D
from weilform import toeplitz, eigsym
from delta2_law import Fv_and_derivs

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
mp.mp.dps = 60
print("=" * 78)
print("TASK 3a -- function field")
print("=" * 78)
print("\n(i) T_R = u w^T + w u^T  with u_i = rho^i, w_i = rho^-i")
worst = mp.mpf(0)
for rho in ['3/2', '13/10', '11/10', '101/100']:
    r = mp.mpf(rho)
    for R in range(1, 10):
        t = [r ** n + r ** (-n) for n in range(R + 1)]
        T = toeplitz(t, R)
        u = [r ** i for i in range(R + 1)]
        w = [r ** (-i) for i in range(R + 1)]
        d = max(abs(T[i, j] - (u[i] * w[j] + w[i] * u[j]))
                for i in range(R + 1) for j in range(R + 1))
        worst = max(worst, d)
print(f"    max |T_ij - (u_i w_j + w_i u_j)| over rho in 4 values, R=1..9 : {mp.nstr(worst,4)}")

print("\n(ii) nonzero eigenvalues are (R+1) +- |u||w|;  all others are 0")
worst2 = worst3 = mp.mpf(0)
for rho in ['3/2', '13/10', '101/100']:
    r = mp.mpf(rho)
    for R in range(1, 9):
        t = [r ** n + r ** (-n) for n in range(R + 1)]
        vals, _ = eigsym(toeplitz(t, R))
        nu = mp.sqrt(sum((r ** i) ** 2 for i in range(R + 1)))
        nw = mp.sqrt(sum((r ** (-i)) ** 2 for i in range(R + 1)))
        pred = sorted([(R + 1) - nu * nw, (R + 1) + nu * nw])
        worst2 = max(worst2, abs(vals[0] - pred[0]), abs(vals[-1] - pred[1]))
        worst3 = max(worst3, max(abs(vals[k]) for k in range(1, R)) if R >= 2 else mp.mpf(0))
print(f"    max |computed - ((R+1) -+ |u||w|)| : {mp.nstr(worst2,4)}")
print(f"    max |middle eigenvalues| (should be 0) : {mp.nstr(worst3,4)}")

print("\n(iii) Cauchy-Schwarz: |u||w| >= u.w = R+1, equality iff rho = 1")
print(f"    u.w computed vs R+1 : ", end='')
ok = all(abs(sum((mp.mpf('13/10') ** i) * (mp.mpf('13/10') ** (-i))
                 for i in range(R + 1)) - (R + 1)) < mp.mpf(10) ** -50
         for R in range(1, 10))
print(ok)

print("\n(iv) centred-moment identity  S2 - S1^2/(R+1) = R(R+1)(R+2)/12  (exact integers)")
bad = []
for R in range(1, 200):
    S1 = sum(range(R + 1)); S2 = sum(i * i for i in range(R + 1))
    lhs = mp.mpf(S2) - mp.mpf(S1) ** 2 / (R + 1)
    rhs = mp.mpf(R * (R + 1) * (R + 2)) / 12
    if abs(lhs - rhs) > mp.mpf(10) ** -40:
        bad.append(R)
print(f"    failures for R = 1..199 : {bad or 'none'}")

print("\n(v) expansion in a = log(rho) vs in eps = rho-1")
print(f"    {'R':>4s} {'-lam/a^2':>18s} {'-lam/eps^2':>18s} {'binom(R+2,3)':>14s}")
for R in (1, 2, 4, 8, 16, 32):
    eps = mp.mpf(10) ** -8
    r = 1 + eps
    a = mp.log(r)
    u2 = (r ** (2 * R + 2) - 1) / (r ** 2 - 1)
    w2 = (r ** (-2 * R - 2) - 1) / (r ** -2 - 1)
    lam = (R + 1) - mp.sqrt(u2 * w2)
    b = mp.mpf(R * (R + 1) * (R + 2)) / 6
    print(f"    {R:>4d} {mp.nstr(-lam/a**2,12):>18s} {mp.nstr(-lam/eps**2,12):>18s} "
          f"{mp.nstr(b,10):>14s}")

print("\n" + "=" * 78)
print("TASK 3b -- zeta side")
print("=" * 78)
g1 = mp.im(mp.zetazero(1))
print(f"\n(i) normalisation: which convention gives  dQ = -4 d^2 (F'^2 + F F'') ?")
print("    (count-preserving: a DOUBLE on-line pair at gamma -> one off-line quartet)")
L = mp.mpf('0.8'); N = 16
M = build_matrix(L, N, 60, EVEN_D)
vals, vecs = eigsym(M)
v = vecs[0]; nv = mp.sqrt(sum(x ** 2 for x in v)); v = [x / nv for x in v]
F0, F1, F2 = Fv_and_derivs(v, L, g1)
nrm = [mp.sqrt(norm2(k, L, EVEN_D)) for k in range(N)]
def quartet(d):
    w = mp.mpc(g1, d)
    Sw = sum(v[k] * F(k, L, w, EVEN_D) / nrm[k] for k in range(N))
    Sm = sum(v[k] * F(k, L, -w, EVEN_D) / nrm[k] for k in range(N))
    return 4 * mp.re(Sw * Sm)
print(f"    {'delta':>10s} {'dQ = -4F^2 + quartet':>24s} {'-4 d^2 (Fp^2+FFpp)':>22s} "
      f"{'ratio':>12s}")
for d in ['1e-2', '1e-3', '1e-4']:
    d = mp.mpf(d)
    dQ = -4 * F0 ** 2 + quartet(d)
    pred = -4 * d ** 2 * (F1 ** 2 + F0 * F2)
    print(f"    {mp.nstr(d,4):>10s} {mp.nstr(dQ,14):>24s} {mp.nstr(pred,14):>22s} "
          f"{mp.nstr(dQ/pred,12):>12s}")

print("\n(ii) Cauchy-Schwarz bounds (crude, valid for ALL gamma)")
for Lx in ['0.5', '0.8', '1.0', '1.3']:
    Lv = mp.mpf(Lx)
    b1 = 2 * Lv ** 3 / 3
    b2 = 2 * Lv ** 5 / 5
    b0 = 2 * Lv
    num1 = mp.quad(lambda u: u ** 2 * mp.sin(g1 * u) ** 2, [-Lv, 0, Lv])
    num2 = mp.quad(lambda u: u ** 4 * mp.cos(g1 * u) ** 2, [-Lv, 0, Lv])
    num0 = mp.quad(lambda u: mp.cos(g1 * u) ** 2, [-Lv, 0, Lv])
    print(f"    L={Lx}: int u^2 sin^2 = {mp.nstr(num1,6)} <= 2L^3/3 = {mp.nstr(b1,6)} : "
          f"{num1 <= b1};  int u^4 cos^2 = {mp.nstr(num2,6)} <= 2L^5/5 = {mp.nstr(b2,6)} : "
          f"{num2 <= b2};  int cos^2 = {mp.nstr(num0,6)} <= 2L = {mp.nstr(b0,6)} : {num0 <= b0}")

print("\n(iii) c(L) = 8 L^3 (1/3 + 1/sqrt 5)  vs the MEASURED C(L)")
meas = {'0.8': '0.3107976', '1.0': '0.85303433', '1.3': '2.4986984',
        '1.6': '4.9217908', '2.0': '10.082833'}
print(f"    {'L':>5s} {'measured C':>13s} {'c(L) bound':>13s} {'C <= c(L)':>10s} "
      f"{'C/c(L)':>9s}")
for Lx, Cm in meas.items():
    Lv = mp.mpf(Lx); Cm = mp.mpf(Cm)
    c = 8 * Lv ** 3 * (mp.mpf(1) / 3 + 1 / mp.sqrt(5))
    print(f"    {Lx:>5s} {mp.nstr(Cm,8):>13s} {mp.nstr(c,8):>13s} {str(Cm <= c):>10s} "
          f"{mp.nstr(Cm/c,6):>9s}")

print("\n(iv) remainder bound  |R_4| <= (16/3) L^5 delta^4 e^{2 L delta}")
for Lx in ['0.8']:
    Lv = mp.mpf(Lx)
    for d in ['1e-1', '1e-2', '1e-3']:
        d = mp.mpf(d)
        actual = abs((-4 * F0 ** 2 + quartet(d)) - (-4 * d ** 2 * (F1 ** 2 + F0 * F2)))
        bound = mp.mpf(16) / 3 * Lv ** 5 * d ** 4 * mp.e ** (2 * Lv * d)
        print(f"    L={Lx} delta={mp.nstr(d,4)}: |actual remainder| = {mp.nstr(actual,6)}"
              f"  <= bound = {mp.nstr(bound,6)} : {actual <= bound}")
