#!/bin/bash
# Voxel-size and time-step test: one e_coli_core cell in a closed 16 um box, dx in {0.8,1.6,3.2} x dt in {0.005,0.01,0.02},
# an "inactive" series (cap not limiting; outputs must agree across dx and dt) and a "capped" series (glucose D=0, cap must act
# and empty the cell's voxel without overdraw), plus ledger residual and carbon balance, and the elemental balance of the model.
# Needs a staged dfba-unit-test binary. Single thread per run, runs in parallel (JOBS, default 4).
# usage (from the staged root): scripts/run_dxdt_test.sh [./dfba-unit-test] [max_time_min]
# Output goes to output_dxdt_<kind>_dx<dx>_dt<dt>/; existing results are moved aside, never deleted.
BIN=${1:-./dfba-unit-test}
T=${2:-120}
PY=${PYTHON:-python3}
JOBS=${JOBS:-4}
names=$($PY scripts/gen_dxdt_configs.py config $T) || { echo "could not generate configs"; exit 2; }
run_one() {
    name=$1; out=output_$name
    if [ -d "$out" ] && [ -n "$(ls -A $out 2>/dev/null)" ]; then mv "$out" "${out}_$(date +%Y%m%d_%H%M%S)"; fi
    mkdir -p "$out"
    "$BIN" config/Unit_test_$name.xml > "$out.log" 2>&1 || echo "simulation failed for $name, see $out.log"
}
export -f run_one; export BIN
echo "$names" | xargs -P "$JOBS" -I{} bash -c 'run_one {}'
rc=0
$PY scripts/check_dxdt.py || rc=$?
$PY scripts/elemental_check.py config/Ecoli_core.xml | tee "${TMPDIR:-/tmp}/dxdt_elemental.txt"
grep -q "elementally unbalanced 0" "${TMPDIR:-/tmp}/dxdt_elemental.txt" && echo "PASS E model elemental balance" || { echo "FAIL E model elemental balance"; [ $rc -ne 0 ] || rc=1; }
exit $rc
