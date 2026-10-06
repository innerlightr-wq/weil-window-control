#!/usr/bin/env bash
# Connes's [1,13] window at high precision, plus even/odd gap and recorder split.
# ~60 min at N=36/140 digits.
set -euo pipefail
cd "$(dirname "$0")/.."
python3 -c "
import sys; sys.path.insert(0,'src')
import mpmath as mp, csv, os
from stage3_connes import part3b
mp.mp.dps = 140
rows = []; part3b(140, 36, 50, rows)
with open(os.path.join('data','stage3b_connes_zeros_N36.csv'),'w',newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print('wrote data/stage3b_connes_zeros_N36.csv')
"
python3 src/stage3_connes.py --N 12 --dps 40 --parts c,d
python3 src/task2_inertia_reconcile.py
