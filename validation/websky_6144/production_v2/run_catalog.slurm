#!/bin/bash
#SBATCH --job-name=v2-cat
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b2
#SBATCH --nodes=1
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=32
#SBATCH --mem=480G
#SBATCH --time=8:00:00
# Production-v2 campaign: one octant catalog (periodic cores, cf32, z_max=4.5), then
# abundance matching (fixed top-halo mapping). Pass OCT, OBSX/OBSY/OBSZ via --export
# (submit_all.sh). Expected ~2.6 h with the data-sized merge hash (was ~7 h).
set -euo pipefail
module load julia/1.12.5 cuda/12.6
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"
unset CUDA_VISIBLE_DEVICES
: "${OCT:?}"; : "${OBSX:?}"; : "${OBSY:?}"; : "${OBSZ:?}"
cd /home/yguan/work/PeakPatch.jl
CONFIG="validation/websky_6144/production_v2/config_v2_oct${OCT}.toml"
D=/home/yguan/projects/aip-aspuru-ab/yguan/websky
RAW=${D}/catalog_websky_6144_v2_oct${OCT}.pksc
AM=${D}/catalog_websky_6144_v2_oct${OCT}_AM.pksc
echo "Node: $(hostname)  Octant: ${OCT}  obs=(${OBSX},${OBSY},${OBSZ})  commit: $(git rev-parse --short HEAD)  Start: $(date)"
echo "=== [1/2] halo finding $(date) ==="
julia --project=validation --threads=${JULIA_NUM_THREADS} validation/websky_6144/run_gpu_octant.jl "${CONFIG}"
echo "=== [2/2] abundance matching $(date) ==="
julia --project=validation -t ${SLURM_CPUS_PER_TASK} validation/websky_6144/apply_abundance_match.jl \
    "${RAW}" "${AM}" "${OBSX},${OBSY},${OBSZ}" 4.5
echo "End: $(date)"
