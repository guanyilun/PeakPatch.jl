#!/bin/bash
# Resubmit some or all of the v5fs analysis steps (run_analysis.slurm) without the rest of the chain.
# Works from any directory:  bash validation/websky_6144/production_v5/rerun_analysis.sh [step ...]
# steps: assemble seams autos cross theory tails (default: all)
set -euo pipefail
export REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
source "${REPO}/validation/websky_6144/production_v5/env.sh"
export STEPS="${*:-assemble seams autos cross theory tails}"
sbatch --export=ALL --chdir="${REPO}" --output="${WS_ROOT}/logs/v5fs_analysis_%j.out" \
    "${REPO}/validation/websky_6144/production_v5/run_analysis.slurm"
echo "steps: ${STEPS}"
