#!/bin/bash
#SBATCH --job-name=pp-composite
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b1
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1:30:00
#SBATCH --output=/home/yguan/scratch/websky_6144/logs/composite_%j.out
#SBATCH --error=/home/yguan/scratch/websky_6144/logs/composite_%j.err

# CPU-only composite kappa analysis (the GPU is just the smallest schedulable
# allocation unit on this cluster). Minimal request -> fast backfill.
set -euo pipefail
module load julia/1.12.5
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"

cd /home/yguan/work/PeakPatch.jl
FMAPS=/home/yguan/scratch/websky_6144/fieldmaps
FIELD="${1:-${FMAPS}/kappa_websky_6144_oct000_finecell_field_nside2048.fits}"
julia --project=. --threads=${JULIA_NUM_THREADS} \
    validation/websky_6144/compare_composite_kappa.jl \
    "${FIELD}" \
    "${FMAPS}/kappa_websky_6144_oct000_finecell_nside2048.fits"
