#!/bin/bash
#SBATCH --job-name=pp-ksz-batt
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b1
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1:30:00
#SBATCH --output=/home/yguan/scratch/websky_6144/logs/composite_kszbatt_%j.out
#SBATCH --error=/home/yguan/scratch/websky_6144/logs/composite_kszbatt_%j.err

# CPU-only composite kSZ v2 (Battaglia halo model, Fortran-verbatim); minimal
# 1-GPU slice for fast backfill.
set -euo pipefail
module load julia/1.12.5
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"

cd /home/yguan/work/PeakPatch.jl
FMAPS=/home/yguan/scratch/websky_6144/fieldmaps
FIELD="${1:-${FMAPS}/ksz_websky_6144_oct000_finecell_cf32_nside4096.fits}"
REF="${2:-/home/yguan/scratch/websky_6144/websky_ref/ksz.fits}"
CATA="${3:-}"
CONV="${4:-}"
SHUF="${5:-}"
julia --project=. --threads=${JULIA_NUM_THREADS} \
    validation/websky_6144/compare_composite_ksz_battaglia.jl \
    "${FIELD}" "${REF}" ${CATA:+"${CATA}"} ${CONV:+"${CONV}"} ${SHUF:+"${SHUF}"}
