#!/bin/bash
#SBATCH --job-name=v5-cat
#SBATCH --partition=gpu
#SBATCH --constraint=a100-80gb
#SBATCH --nodes=1
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=32
#SBATCH --mem=160G
#SBATCH --time=4:00:00
# Production-v5: one octant catalog (halo finding, merge, finalize, write); AM follows with one
# full-sky table (submit_all.sh). Pass OCT and REPO via --export.
# Sizing from Killarney (4×L40S): v4 job 1.62-1.95 h, v5test +9%; host MaxRSS 89.7 GiB (v4),
# 96 GB (v5test), 134.5 GiB (v3, before the per-tile threshold). GPU memory per tile ~22 GiB.
# The last ~45-60 min (merge + finalize) is CPU-only with the GPUs idle.
# Do NOT unset CUDA_VISIBLE_DEVICES (the Killarney script did): gpu-partition nodes are shared,
# and Slurm's value is what limits run_gpu_octant.jl to this job's 4 GPUs.
set -euo pipefail
: "${REPO:?}" "${OCT:?}"
source "${REPO}/validation/websky_6144/production_v5/env.sh"
cd "${REPO}"
CONFIG="validation/websky_6144/production_v5/config_v5_oct${OCT}.toml"
echo "Node: $(hostname)  Octant: ${OCT}  commit: $(git rev-parse --short HEAD)  GPUs: ${CUDA_VISIBLE_DEVICES:-?}  Start: $(date)"
julia --project=validation --threads=${JULIA_NUM_THREADS} validation/websky_6144/run_gpu_octant.jl "${CONFIG}"
echo "End: $(date)"
