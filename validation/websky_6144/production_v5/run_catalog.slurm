#!/bin/bash
#SBATCH --job-name=v5-cat
#SBATCH --partition=gpu
#SBATCH --constraint=a100-80gb&rocky9
#SBATCH --nodes=1
#SBATCH --gpus-per-node=2
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=3:30:00
# Production-v5: one octant catalog (halo finding, merge, finalize, write); AM follows with one
# full-sky table (submit_all.sh). Pass OCT and REPO via --export.
# Sizing, measured on Rusty (v5, 2×A100, jobs 7206249+2k): 1.75-1.82 h, MaxRSS 90.1-93.5 GiB.
# (Killarney 4×L40S: v4 1.62-1.95 h, MaxRSS 89.7 GiB; v3 134.5 GiB before the per-tile threshold.)
# The last ~45-60 min (merge + finalize) is CPU-only with the GPUs idle.
# 2 GPUs, not 4 (2026-10-09): 4-GPU jobs waited many hours for a fully free node on a busy cluster,
# while 2 free GPUs on a node are common. 8 jobs x 2 GPUs = the 16-GPU gpu QoS cap, so all octants
# can run at once. The multi-GPU
# dispatcher matches 1-GPU results (test_multigpu_dispatch.jl, passed on Rusty). Host memory is set by the merge.
# Do NOT unset CUDA_VISIBLE_DEVICES (the Killarney script did): gpu-partition nodes are shared,
# and Slurm's value is what limits run_gpu_octant.jl to this job's GPUs.
# rocky9: env.sh loads modules from the Rocky 9 tree; the 3 Rocky 8 A100 nodes (workergpu038-040) are untested.
set -euo pipefail
: "${REPO:?}" "${OCT:?}"
source "${REPO}/validation/websky_6144/production_v5/env.sh"
cd "${REPO}"
CONFIG="validation/websky_6144/production_v5/config_v5_oct${OCT}.toml"
echo "Node: $(hostname)  Octant: ${OCT}  commit: $(git rev-parse --short HEAD)  GPUs: ${CUDA_VISIBLE_DEVICES:-?}  Start: $(date)"
julia --project=validation --threads=${JULIA_NUM_THREADS} validation/websky_6144/run_gpu_octant.jl "${CONFIG}"
echo "End: $(date)"
