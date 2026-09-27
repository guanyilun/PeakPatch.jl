#!/bin/bash
#SBATCH --job-name=pp6144-oct${OCTANT:-0}
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b3
#SBATCH --nodes=1
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=32
#SBATCH --mem=400G
#SBATCH --time=06:00:00
#SBATCH --output=/home/yguan/scratch/websky_6144/logs/oct${OCTANT:-0}_%j.out
#SBATCH --error=/home/yguan/scratch/websky_6144/logs/oct${OCTANT:-0}_%j.err

# Websky 6144^3 GPU multi-resolution run on Killarney.
# Runs one octant of the full-sky lightcone.
#
# Usage:
#   sbatch --export=OCTANT=0 run_websky_6144_killarney.slurm
#
# All 8 octants in parallel:
#   for oct in 000 001 010 011 100 101 110 111; do
#       sbatch --export=OCTANT=$oct run_websky_6144_killarney.slurm
#   done

set -euo pipefail

module load julia/1.12.5 cuda/12.6

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
# Keep all 4 L40S visible so devices=[0,1,2,3] works
unset CUDA_VISIBLE_DEVICES

cd /home/yguan/work/PeakPatch.jl

# OCTANT is a 3-digit binary name (000-111).  Force decimal base to avoid
# bash treating leading zeros as octal (e.g. 010 → 8 instead of ten).
OCT=${OCTANT:-000}

# Config file for this octant
CONFIG="validation/websky_6144/config_websky_6144_oct${OCT}.toml"

if [ ! -f "$CONFIG" ]; then
    echo "ERROR: Config file not found: $CONFIG"
    echo "Available configs:"
    ls validation/websky_6144/config_websky_6144_oct*.toml 2>/dev/null || echo "  (none)"
    exit 1
fi

OUTDIR="/home/yguan/scratch/websky_6144"
mkdir -p ${OUTDIR}/logs

echo "=== PeakPatch.jl Websky 6144^3 — GPU Multi-Resolution ==="
echo "Octant:  ${OCT}"
echo "Config:  ${CONFIG}"
echo "Node:    $(hostname)"
echo "GPUs:    ${SLURM_GPUS_PER_NODE:-?}"
echo "Threads: ${JULIA_NUM_THREADS}"
echo "Start:   $(date)"

julia --project=validation --threads=${JULIA_NUM_THREADS} \
    validation/websky_6144/run_gpu_octant.jl "${CONFIG}"

echo "End: $(date)"
