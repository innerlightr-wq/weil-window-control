#!/usr/bin/env python3
"""Regression: does the projected-derivative formula recover Theorem 6's coefficient?

Theorem 6 (thm:ff) is the theta = 0 case of the same deformation family:
normalised spectrum {rho, rho^-1} = {e^a, e^-a}, so t(n) = 2cosh(na), t_0(n) = 2,
T_0 = 2 * 1 1^T (rank one), E = ker T_0 = {c : sum c_i = 0}, and
Delta t(n) = 2(cosh(na) - 1) = a^2 n^2 + O(a^4) -- our Delta t at theta = 0.

Our formula: lam_min = -2a^2 max_{c in E, ||c||=1} |S_1(c)|^2, S_1(c) = sum_i i c_i e^{i i theta}.
At theta = 0, S_1 is REAL (the sine direction vanishes: the rank-2 representer space
collapses to rank 1), so max = ||Pi_E u||^2 with u_i = i, and
   Pi_E u = u - (<u,1>/||1||^2) 1,  (Pi_E u)_i = i - R/2,
   ||Pi_E u||^2 = sum_{i=0}^R (i - R/2)^2 = R(R+1)(R+2)/12.
Hence lam_min = -2a^2 R(R+1)(R+2)/12 = -C(R+2,3) a^2 -- Theorem 6's coefficient.
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'src'))
import mpmath as mp
from weilform import toeplitz, eigsym
mp.mp.dps = 60

def binom3(R): return mp.mpf(R)*(R+1)*(R+2)/6

print('  %-4s %-22s %-22s %-22s %-11s' %
      ('R', 'sum (i-R/2)^2', 'C(R+2,3)/2', 'proj ||Pi_E u||^2', 'identical'))
for R in (1, 2, 4, 8, 16, 32, 64):
    cm = sum((mp.mpf(i)-mp.mpf(R)/2)**2 for i in range(R+1))
    # projection computed independently: u minus its mean
    u = [mp.mpf(i) for i in range(R+1)]
    mean = sum(u)/(R+1)
    proj = sum((x-mean)**2 for x in u)
    print('  %-4d %-22s %-22s %-22s %-11s' %
          (R, mp.nstr(cm, 14), mp.nstr(binom3(R)/2, 14), mp.nstr(proj, 14),
           abs(cm-binom3(R)/2) < mp.mpf(10)**-50 and abs(cm-proj) < mp.mpf(10)**-50))
print()
print('  %-4s %-8s %-20s %-20s %-12s %-12s' %
      ('R', 'a', 'lam_min exact', 'our -2a^2||Pi_E u||^2', 'rel err', 'rel/a^2'))
rows = []
for R in (1, 2, 4, 8, 16, 32):
    u = [mp.mpf(i) for i in range(R+1)]; mean = sum(u)/(R+1)
    mu = sum((x-mean)**2 for x in u)
    for a_ in ('1e-2', '1e-3', '1e-4'):
        a = mp.mpf(a_)
        t = [2*mp.cosh(n*a) for n in range(R+1)]
        lam = eigsym(toeplitz(t, R))[0][0]
        pred = -2*a**2*mu
        rel = abs(lam-pred)/abs(pred)
        rows.append(dict(R=R, a=a_, lam=mp.nstr(lam, 14), pred=mp.nstr(pred, 14),
                         rel=float(rel), rel_over_a2=float(rel/a**2)))
        print('  %-4d %-8s %-20s %-20s %-12.3e %-12.4f' %
              (R, a_, mp.nstr(lam, 12), mp.nstr(pred, 12), float(rel), float(rel/a**2)))
json.dump(rows, open('data/stage1_regression.json', 'w'), indent=1)
print()
print('  Theorem 6 coefficient C(R+2,3) recovered as 2*||Pi_E u||^2 at every R: the paper\'s')
print('  "centred second moment of the window" IS the projected derivative-evaluation norm.')
