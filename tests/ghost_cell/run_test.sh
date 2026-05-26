#!/usr/bin/env bash
# Ghost-cell MPI validation test.
#
# Runs the WTM on a 102x3 heterogeneous-ksat domain with 1 and then 2 MPI
# processes, then compares the output TIFs.  With the ghost-cell fix the two
# runs agree; without it they diverge at the MPI processor boundary.
#
# Usage:
#   cd tests/ghost_cell
#   ./run_test.sh [path/to/wtm.x]
#
# Default binary: ../../build/wtm.x

set -euo pipefail
cd "$(dirname "$0")"

WTM=${1:-../../build/wtm.x}

if [[ ! -x "$WTM" ]]; then
    echo "ERROR: WTM binary not found at $WTM" >&2
    echo "Build the project first (cmake --build build) or supply the path as \$1." >&2
    exit 1
fi

echo "=== Ghost-cell MPI validation test ==="
echo "WTM binary: $WTM"
echo

# 1. Generate synthetic inputs
echo "--- Generating inputs ---"
python3 make_inputs.py

# 2. One-process reference run
echo
echo "--- 1-process reference run ---"
rm -rf out_1p run_1p.txt
mkdir -p out_1p
# outfile_prefix and textfilename are relative to CWD when wtm.x is called.
# Override them with sed-generated temp configs so both runs share the same
# base config without modifying it.
CFG_1P=$(mktemp /tmp/ghost_cell_1p_XXXXXX.cfg)
sed 's|^outfile_prefix.*|outfile_prefix     out_1p/out_|;
     s|^textfilename.*|textfilename       run_1p.txt|' ghost_cell.cfg > "$CFG_1P"
mpirun -n 1 "$WTM" "$CFG_1P" \
    -snes_mf \
    -snes_stol 1e-6 \
    2>&1 | grep -E 'SNES|converged|norm|Error|error' || true
echo "1-process run complete."
rm -f "$CFG_1P"

# 3. Two-process run (processor boundary at the ksat discontinuity)
echo
echo "--- 2-process run (split at ksat boundary) ---"
rm -rf out_2p run_2p.txt
mkdir -p out_2p
CFG_2P=$(mktemp /tmp/ghost_cell_2p_XXXXXX.cfg)
sed 's|^outfile_prefix.*|outfile_prefix     out_2p/out_|;
     s|^textfilename.*|textfilename       run_2p.txt|' ghost_cell.cfg > "$CFG_2P"
mpirun -n 2 "$WTM" "$CFG_2P" \
    -snes_mf \
    -snes_stol 1e-6 \
    -da_processors_x 2 -da_processors_y 1 \
    2>&1 | grep -E 'SNES|converged|norm|Error|error' || true
echo "2-process run complete."
rm -f "$CFG_2P"

# 4. Compare outputs
echo
echo "--- Comparing outputs ---"
python3 check_results.py
