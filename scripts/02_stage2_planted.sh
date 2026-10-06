#!/usr/bin/env bash
# Stage 2 -- planted off-circle spectra, exact over Q, + detection law + adversarial pass.
# ~6 min.
set -euo pipefail
cd "$(dirname "$0")/.."
python3 src/stage2.py --dps 50 --Rmax 20
python3 src/detection_law.py
python3 src/adversarial.py
python3 src/function_level.py --dps 50 --Rmax 14
python3 src/gram_trace.py
python3 src/figs_stage2.py
python3 src/figs_function_level.py
