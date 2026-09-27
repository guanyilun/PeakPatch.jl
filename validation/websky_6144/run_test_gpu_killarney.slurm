#!/bin/bash
#SBATCH --job-name=pp-gpu-test
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_interac
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --output=/home/yguan/scratch/websky_6144/logs/gpu_test_%j.out
#SBATCH --error=/home/yguan/scratch/websky_6144/logs/gpu_test_%j.err

# Small-scale GPU pipeline test on Killarney L40S
set -euo pipefail

module load julia/1.12.5 cuda/12.6

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"
unset CUDA_VISIBLE_DEVICES

cd ${SLURM_SUBMIT_DIR}

echo "=== PeakPatch.jl GPU Test ==="
echo "Node:    $(hostname)"
echo "GPU:     $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'unknown')"
echo "Threads: ${JULIA_NUM_THREADS}"
echo "Start:   $(date)"

julia --project=validation --threads=${JULIA_NUM_THREADS} \
    validation/websky_6144/run_test_gpu.jl

echo "End: $(date)"
