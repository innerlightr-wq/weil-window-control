#!/bin/sh
# Build the reproducibility archive for the technical note.
#
# Redistributable files only: everything below is authored in this repository
# and covered by its MIT (code) / CC BY 4.0 (documents) split.  No third-party
# material is included -- in particular no copy of arXiv:2608.24827, whose
# reproducibility materials were never obtained and whose matrix step is not
# reproduced here.
#
# The archive is built DETERMINISTICALLY: sorted entries, fixed mtime, numeric
# owner 0:0, and gzip without a timestamp.  Two builds of an unchanged tree give
# byte-identical output, so the checksum in SHA256SUMS.txt is reproducible.
#
# LaTeX build products under note/, caches and the archive itself are excluded.
# The research logs under */logs/ ARE included: CLAIM_INDEX.md cites them for
# the Rump residuals and the measured timings.
set -e
cd "$(dirname "$0")/.."
EXT=$(basename "$PWD")
VER=1.2.0
OUT="note/weil-L45-certificate-$VER.tar.gz"

cd ..
tar --create --file - \
    --sort=name --mtime="UTC 2026-11-10" \
    --owner=0 --group=0 --numeric-owner --format=gnu \
    --exclude='note/*.aux' --exclude='note/*.log' --exclude='note/*.out' \
    --exclude='note/*.fls' --exclude='note/*.fdb_latexmk' --exclude='note/*.synctex.gz' \
    --exclude='__pycache__' --exclude='*.pyc' \
    --exclude='weil-L45-certificate-*.tar.gz' \
    --exclude='.gitignore' --exclude='SHA256SUMS.txt' \
    "$EXT/README.md" "$EXT/HYPOTHESES.md" "$EXT/DERIVATION.md" \
    "$EXT/SOURCE_AUDIT.md" "$EXT/REPORT.md" "$EXT/RETRACTIONS.md" \
    "$EXT/PROVENANCE.txt" "$EXT/verify.py" \
    "$EXT/src" "$EXT/certs" "$EXT/data" "$EXT/logs" \
    "$EXT/phase2/src" "$EXT/phase2/logs" "$EXT/phase2/certs" "$EXT/phase2/verify_phase2.py" "$EXT/phase2/verify_phase3.py" \
    "$EXT/phase3/src" "$EXT/phase3/logs" "$EXT/phase3/certs" "$EXT/phase3/verify_audit.py" \
    "$EXT/phase4/src" "$EXT/phase4/logs" "$EXT/phase4/certs" "$EXT/phase4/verify_t400.py" \
    "$EXT/phase5/src" "$EXT/phase5/logs" "$EXT/phase5/certs" "$EXT/phase5/test_prec.py" \
    "$EXT/note" \
  | gzip -9 -n > "$EXT/$OUT"

echo "built $EXT/$OUT"
cd "$EXT" && ls -l "$OUT"
