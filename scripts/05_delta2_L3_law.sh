#!/usr/bin/env bash
# The delta^2--L^3 law: proof checks, C(L) vs Cauchy-Schwarz, inertia vs N.  ~25 min.
set -euo pipefail
cd "$(dirname "$0")/.."
python3 src/proofs_check.py
python3 -c "import sys; sys.path.insert(0,'src'); import proofs_check as p; p.task3b_sharp()"
python3 src/delta2_law.py
python3 src/task34_extra.py
