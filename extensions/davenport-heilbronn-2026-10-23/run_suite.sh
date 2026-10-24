#!/usr/bin/env bash
# Davenport-Heilbronn extension: the whole suite, in dependency order.
# Stops at the first failing gate rather than continuing past it.
set -euo pipefail
cd "$(dirname "$0")/../.."
PY=.venv/bin/python
E=extensions/davenport-heilbronn-2026-10-23
mkdir -p "$E/logs"

run () {
  echo "=== $* ==="
  "$PY" "$@" 2>&1 | tee -a "$E/logs/suite.log"
  return "${PIPESTATUS[0]}"
}

run "$E/tests/regression.py"
run "$E/src/phaseA.py"
run "$E/src/phaseB.py"
run "$E/tests/test_orbits.py"
run "$E/src/census.py" --Tmax 200
run "$E/src/phaseC.py"
run "$E/src/phaseD.py"
run "$E/src/phaseE.py"
run "$E/src/phaseF.py"
run "$E/src/phaseG.py"
run "$E/src/phaseH.py"
run "$E/src/gate_signflip.py"
run "$E/src/certify.py" 400
run "$E/tests/regression.py"
echo "SUITE COMPLETE"
