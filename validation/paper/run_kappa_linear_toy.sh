#!/bin/bash
# Submit kappa_linear_toy.jl (CPU; regenerates the seed-12345 coarse noise, ~5 min, then caches it under
# $WS_ROOT/scratch; 9 linear κ maps at Nside 256 plus one Nside-4096 map read; memory a few GB, not measured).
# Works from any directory:  bash validation/paper/run_kappa_linear_toy.sh
set -euo pipefail
export REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "${REPO}/validation/websky_6144/production_v5/env.sh"
SCRIPT=validation/paper/kappa_linear_toy.jl sbatch --export=ALL --chdir="${REPO}" \
    --job-name=kappa-toy --cpus-per-task=32 --mem=16G --time=1:00:00 \
    --output="${WS_ROOT}/logs/kappa_toy_%j.out" "${REPO}/validation/websky_6144/production_v5/run_jl.slurm"
