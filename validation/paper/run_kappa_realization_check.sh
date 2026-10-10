#!/bin/bash
# Submit kappa_realization_check.jl (CPU; ~2.3e11 noise draws + a 512^3 FFT; memory a few GB, not measured).
# Works from any directory:  bash validation/paper/run_kappa_realization_check.sh
set -euo pipefail
export REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "${REPO}/validation/websky_6144/production_v5/env.sh"
SCRIPT=validation/paper/kappa_realization_check.jl sbatch --export=ALL --chdir="${REPO}" \
    --job-name=kappa-realization --cpus-per-task=32 --mem=32G --time=1:00:00 \
    --output="${WS_ROOT}/logs/kappa_realization_%j.out" "${REPO}/validation/websky_6144/production_v5/run_jl.slurm"
