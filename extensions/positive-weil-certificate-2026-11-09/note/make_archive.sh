#!/bin/sh
# Build the reproducibility archive for the technical note.
#
# Redistributable files only: everything below is authored in this repository
# and covered by its MIT (code) / CC BY 4.0 (documents) split.  No third-party
# material is included -- in particular no copy of arXiv:2608.24827 or
# arXiv:2611.13129, whose texts are cited, not redistributed.
#
# LAYOUT.  Entries are written with their REPOSITORY-RELATIVE paths, so
# extracting the archive into an empty directory reproduces the layout the
# release gate expects:
#
#     <dir>/extensions/positive-weil-certificate-2026-11-09/...
#     <dir>/src/zeta_window.py          (the only repo-level module imported)
#     <dir>/LICENSE  <dir>/LICENSE-CC-BY-4.0
#
# and then
#
#     cd <dir>/extensions/positive-weil-certificate-2026-11-09
#     python phase7/acceptance.py
#
# runs the full gate with no $HOME path and no hidden cache.  The only external
# dependencies are Python 3, python-flint (FLINT/Arb) and mpmath; versions are
# recorded in note/README.md.
#
# The archive is built DETERMINISTICALLY: sorted entries, fixed mtime, numeric
# owner 0:0, and gzip without a timestamp.  Two builds of an unchanged tree give
# byte-identical output, so its checksum is reproducible.
#
# LaTeX build products under note/, caches and the archive itself are excluded.
# The research logs under */logs/ ARE included: CLAIM_INDEX.md cites them for
# the Rump residuals and the measured timings.
set -e
cd "$(dirname "$0")/.."
EXT=$(basename "$PWD")
REL="extensions/$EXT"
VER=2.0.0
OUT="note/weil-two-window-certificate-$VER.tar.gz"
cd ../..

tar --create --file - \
    --sort=name --mtime="UTC 2026-11-12" \
    --owner=0 --group=0 --numeric-owner --format=gnu \
    --exclude='note/*.aux' --exclude='note/*.log' --exclude='note/*.out' \
    --exclude='note/*.fls' --exclude='note/*.fdb_latexmk' --exclude='note/*.synctex.gz' \
    --exclude='__pycache__' --exclude='*.pyc' \
    --exclude='weil-L45-certificate-*.tar.gz' \
    --exclude='weil-two-window-certificate-*.tar.gz' \
    --exclude='.gitignore' --exclude='SHA256SUMS.txt' \
    "LICENSE" "LICENSE-CC-BY-4.0" "src/zeta_window.py" \
    "$REL/README.md" "$REL/HYPOTHESES.md" "$REL/DERIVATION.md" \
    "$REL/SOURCE_AUDIT.md" "$REL/REPORT.md" "$REL/RETRACTIONS.md" \
    "$REL/PROVENANCE.txt" "$REL/verify.py" \
    "$REL/src" "$REL/certs" "$REL/data" "$REL/logs" \
    "$REL/phase2/src" "$REL/phase2/logs" "$REL/phase2/certs" \
    "$REL/phase2/verify_phase2.py" "$REL/phase2/verify_phase3.py" \
    "$REL/phase3/src" "$REL/phase3/logs" "$REL/phase3/certs" "$REL/phase3/verify_audit.py" \
    "$REL/phase4/src" "$REL/phase4/logs" "$REL/phase4/certs" "$REL/phase4/verify_t400.py" \
    "$REL/phase5/src" "$REL/phase5/logs" "$REL/phase5/certs" "$REL/phase5/test_prec.py" \
    "$REL/phase7" \
    "$REL/note" \
  | gzip -9 -n > "$REL/$OUT"

echo "built $REL/$OUT"
cd "$REL" && ls -l "$OUT"
