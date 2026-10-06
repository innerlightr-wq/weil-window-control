#!/usr/bin/env bash
# Detection onset vs L_pred and basis size.  HOURS at these precisions.
set -euo pipefail
cd "$(dirname "$0")/.."
python3 src/stage3e_certify.py --ns 1,5,10,30,80 --deltas 0.1 \
        --factors 0.70,1.00,1.30,1.70,2.20
python3 src/stage3e_certify.py --ns 30,80 --factors 0.85,1.00,1.15 --dps 220 \
        --out stage3e_onset_highprec.csv
python3 src/task1e_C_vs_height.py
