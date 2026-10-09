#!/bin/bash
#SBATCH --job-name=v5-jl
#SBATCH --partition=genx
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1:00:00
# Generic CPU runner for the post-catalog steps (AM, halo paint, CIB, Tier-A). None of them load
# CUDA, so they run on genx (shared nodes, no GPU). submit_all.sh overrides cpus/mem/time per step.
#   sbatch --export=ALL,REPO=...,SCRIPT=validation/x.jl,SARGS="...",[JPROJECT=...] run_jl.slurm
set -euo pipefail
: "${REPO:?}" "${SCRIPT:?}"
source "${REPO}/validation/websky_6144/production_v5/env.sh"
cd "${REPO}"
echo "Node: $(hostname) commit: $(git rev-parse --short HEAD) Start: $(date)"
julia --project="${JPROJECT:-validation}" -t ${JULIA_NUM_THREADS} ${SCRIPT} ${SARGS:-}
echo "End: $(date)"
