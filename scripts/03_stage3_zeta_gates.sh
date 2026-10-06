#!/usr/bin/env bash
# Stage 3 -- zeta window form: Gate A (assembly vs zeros side) and Gate B (lambda*).
# ~45 min; the zetazero calls dominate the first few minutes.
set -euo pipefail
cd "$(dirname "$0")/.."
python3 src/stage3_converge.py --N 4,6,8,10,12,14,16,18,20,24 --dps 50 --zeros 800
python3 src/figs_stage3.py
