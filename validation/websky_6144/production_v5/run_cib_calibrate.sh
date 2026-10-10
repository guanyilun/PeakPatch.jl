#!/bin/bash
# Submit the CIB calibration (cib_calibrate.jl). CPU only; ~12 full-sky Nside 4096 transforms to lmax 3085,
# a few GB of memory (not measured). Works from any directory:
#   bash validation/websky_6144/production_v5/run_cib_calibrate.sh
set -euo pipefail
export REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
source "${REPO}/validation/websky_6144/production_v5/env.sh"
SCRIPT=validation/websky_6144/production_v5/cib_calibrate.jl sbatch --export=ALL --chdir="${REPO}" \
    --job-name=v5-cib-calib --cpus-per-task=32 --mem=48G --time=1:30:00 \
    --output="${WS_ROOT}/logs/v5_cib_calibrate_%j.out" "${REPO}/validation/websky_6144/production_v5/run_jl.slurm"
