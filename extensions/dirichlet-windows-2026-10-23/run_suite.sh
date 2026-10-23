#!/usr/bin/env bash
# Full test suite for extensions/dirichlet-windows-2026-10-23.
# Every script exits non-zero on failure; this runner propagates the first failure.
# Run from the repository root:   ./extensions/dirichlet-windows-2026-10-23/run_suite.sh
set -u
cd "$(dirname "$0")/../.." || exit 2
PY=${PY:-.venv/bin/python}
E=extensions/dirichlet-windows-2026-10-23
rc=0
declare -a NAMES=(
  "Phase 1 self-test (characters, Gauss sums, class-number formula, functional equation)"
  "planting-block regressions (F_at_half, pole-term anchor, degeneration identities)"
  "Step 1: Siegel parity verification (a) pair-vs-quartet (b) sign structure (c) pole"
  "Step 2: function-field parity split, exact over Q"
  "Phase 4d crossover criterion R = 4d^2 K / lambda_2"
  "zeta regression suite (18 values, all from the parent's committed CSVs)"
)
declare -a SCRIPTS=(
  "$E/src/selftest.py"
  "$E/src/selftest_planting.py"
  "$E/src/verify_siegel_parity.py"
  "$E/src/ff_parity_split.py"
  "$E/src/phase4d_crossover.py"
  "$E/src/regression_zeta.py"
)
for i in "${!SCRIPTS[@]}"; do
  echo "=============================================================================="
  echo ">>> ${NAMES[$i]}"
  echo "    $PY ${SCRIPTS[$i]}"
  echo "=============================================================================="
  if $PY -u "${SCRIPTS[$i]}" > "/tmp/suite_$i.log" 2>&1; then
    echo "    PASS   (tail:)"; tail -2 "/tmp/suite_$i.log" | sed 's/^/      /'
  else
    echo "    *** FAIL (exit $?) ***"; tail -25 "/tmp/suite_$i.log" | sed 's/^/      /'
    rc=1
  fi
done
echo "=============================================================================="
if [ "$rc" = "0" ]; then echo "SUITE: ALL PASS"; else echo "SUITE: FAILURES PRESENT"; fi
exit $rc
