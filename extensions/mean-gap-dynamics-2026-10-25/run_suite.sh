#!/bin/sh
# Full suite, in dependency order.  Exits non-zero on the first failure.
# Runtimes in REPRODUCE.md; the whole thing is about an hour without gmpy2.
set -e
cd "$(dirname "$0")"
P=${PYTHON:-python3}

run() { echo; echo "### $1"; $P -u "$1" > "results/$2" 2>&1 || { tail -40 "results/$2"; exit 1; }
        tail -2 "results/$2"; }

run src/check_dictionary.py   stage0_dictionary.txt
run src/means.py              part2_means.txt
run src/toy.py                part4_toy.txt
run src/separation_bound.py   part3_separation.txt
run src/dh_flow.py            part4_dh.txt
run src/certify.py            certificates.txt
run src/alt_evolution.py      part5_alt.txt
run src/zeta_lehmer.py        part4_zeta_lehmer.txt
run src/dh_interaction.py     part3_interaction.txt
run src/zeta_flow.py          part4_zeta.txt

# --- the follow-up constraint audit (AUDIT.md); audit4 needs dh_interaction's zero list
run src/audit1_claims.py      audit1.txt
run src/audit2_constraints.py audit2.txt
run src/audit4_lemma.py       audit4.txt
run src/dh_fixed_tau.py       dh_fixed_tau_nav.txt
run src/check_taylor_vs_integral.py taylor_vs_integral.txt
run src/certify_fixed_tau.py  certificate_fixed_tau.txt

echo; echo "### tests/regression.py"
$P -u tests/regression.py | tail -2
echo; echo "### tests/verify_numbers.py"
$P -u tests/verify_numbers.py | tail -2
echo; echo "### SUITE: ALL PASS"
