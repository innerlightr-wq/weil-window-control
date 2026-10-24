#!/usr/bin/env bash
# Runs the exploration's configurations in order.  Each invocation is independent and
# restartable; the loop lives here rather than in an inline shell command.
set -u
cd "$(dirname "$0")/../../.."
PY=.venv/bin/python
E=extensions/channel-compatibility-exploration-2026-10-24
D=$E/data

run () { echo "### $*"; timeout 7200 "$PY" "$E/src/explore_channels.py" "$@" 2>&1; echo; }

# precision ladder on the two most delicate configurations (the L=1.8 crossing and the
# worst-conditioned background), in the SAME run so the comparison is apples to apples
run --dps 80,140,220 --windows 1.8 --levels 3 --out channel_ladder_18.csv   > "$D/ladder_18.log"
run --dps 80,140,220 --windows 2.0 --levels 2 --out channel_ladder_20.csv   > "$D/ladder_20.log"

# eta decomposition: where does the normalized interaction live?
run --dps 80,140 --windows '' --levels '' --decomp 2.0:2                    > "$D/decomp_20_lvl2.log"
run --dps 80     --windows '' --levels '' --decomp 1.8:3                    > "$D/decomp_18_lvl3.log"
run --dps 80     --windows '' --levels '' --decomp 1.6:2                    > "$D/decomp_16_lvl2.log"

# the nested enrichment identity, at the window where enrichment crosses the threshold
run --dps 80,140 --windows '' --levels '' --enrichment 1.8                  > "$D/enrichment_18.log"

# coordinate invariance, and the two-channel case
run --dps 80 --windows '' --levels '' --invariance 2.0:2                    > "$D/invariance.log"
run --dps 80 --windows '' --levels '' --twochannel 2.0:2                    > "$D/twochannel_20.log"
echo "ALL CONFIGURATIONS DONE"
