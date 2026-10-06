#!/usr/bin/env bash
# Fast self-test: exact arithmetic and the function-field proof, no high-precision runs.
# Target: under 5 minutes.
set -euo pipefail
cd "$(dirname "$0")/.."
echo "== finite fields =="
python3 -c "
import sys; sys.path.insert(0,'src')
from ffield import GF
for q,n in [(5,1),(5,2),(3,3)]:
    F=GF(q,n)
    assert all(F.pow(x,F.Q-1)==F.one for x in F.elements() if x!=F.zero)
print('  Fermat check passed for F_5, F_25, F_27')
"
echo "== LDL inertia + witness self-test =="
python3 -c "
import sys; sys.path.insert(0,'src')
import mpmath as mp; mp.mp.dps=40
from ldl import ldl, first_negative_pivot, witness, rayleigh
import random; random.seed(3)
bad=0
for _ in range(60):
    n=random.randint(2,7)
    A=mp.zeros(n,n)
    for i in range(n):
        for j in range(i,n): A[i,j]=A[j,i]=mp.mpf(random.uniform(-3,3))
    try: Lm,D=ldl(A)
    except ZeroDivisionError: continue
    for k in range(1,n+1):
        sub=mp.matrix([[A[i,j] for j in range(k)] for i in range(k)])
        ev,_=mp.eigsy(sub)
        if sum(1 for t in range(k) if ev[t]<0)!=sum(1 for t in range(k) if D[t]<0): bad+=1
assert bad==0, bad
print('  inertia of every leading block matches eigsy on 60 random matrices')
"
echo "== Stage 1 at reduced precision (regression values) =="
python3 src/stage1.py --dps 30 | grep -E "H3 \(|variational|max \| \|alpha"
echo "== Stage 2 exact-over-Q first-negative windows =="
python3 -c "
import sys; sys.path.insert(0,'src')
from exact_inertia import rational_t, exact_inertia
exp={'realpair':1,'quartet':2,'mixed':4}
cfg={'realpair':[('realpair','13/10')],
     'quartet':[('quartet','13/10','19/25')],
     'mixed':[('quartet','13/10','19/25'),('online','-21/50')]}
for k,blocks in cfg.items():
    tq,_=rational_t(blocks,14)
    first=next((R for R in range(15) if (lambda e: e and e[0]>=1)(exact_inertia(tq,R))),None)
    assert first==exp[k], (k,first,exp[k])
    print('  %-9s first negative window R=%d  (exact over Q)'%(k,first))
"
echo "== function-field delta^2 law =="
python3 -c "
import sys; sys.path.insert(0,'src')
import mpmath as mp; mp.mp.dps=60
def lam(a,R):
    r=mp.e**mp.mpf(a)
    u2=(r**(2*R+2)-1)/(r**2-1); w2=(r**(-2*R-2)-1)/(r**-2-1)
    return (R+1)-mp.sqrt(u2*w2)
for R in (1,2,4,8,16,32):
    b=mp.mpf(R*(R+1)*(R+2))/6
    got=-lam('1e-8',R)/mp.mpf('1e-8')**2
    assert abs(got-b)/b < 1e-10, (R,got,b)
print('  lam_min = -binom(R+2,3) a^2 confirmed for R = 1,2,4,8,16,32')
"
echo
echo "QUICK CHECK PASSED"
