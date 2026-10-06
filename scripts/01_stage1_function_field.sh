#!/usr/bin/env bash
# Stage 1 -- function-field control on genuine spectra.  ~2 min.
set -euo pipefail
cd "$(dirname "$0")/.."
python3 src/stage1.py --dps 50
python3 src/figs_stage1.py
