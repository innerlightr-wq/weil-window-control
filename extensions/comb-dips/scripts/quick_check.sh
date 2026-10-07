#!/usr/bin/env bash
# comb-dips quick check -- Stage 1 gate + two recomputed Stage 2 rows.  Under 5 minutes.
set -euo pipefail
cd "$(dirname "$0")/.."
exec python3 -u src/quick_check.py
