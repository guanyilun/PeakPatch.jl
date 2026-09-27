#!/bin/bash
#SBATCH --job-name=pp-prod-cat
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b2
#SBATCH --nodes=1
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=32
#SBATCH --mem=480G
#SBATCH --time=10:00:00

# FROZEN production campaign (paper): one octant catalog at coarse_factor=32,
# z_max=4.5, then abundance matching. Generic over octants: pass OCT (e.g. 010)
# and OBSX/OBSY/OBSZ via sbatch --export (see submit_all.sh, which also sets
# per-octant log paths — #SBATCH directives cannot expand variables).
# Timing anchor: oct000 cf32 z<4.6 took 7h13m on 4x L40S.
set -euo pipefail
module load julia/1.12.5 cuda/12.6
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"
unset CUDA_VISIBLE_DEVICES

: "${OCT:?need OCT=zyx octant bits}"; : "${OBSX:?}"; : "${OBSY:?}"; : "${OBSZ:?}"

cd /home/yguan/work/PeakPatch.jl
CONFIG="validation/websky_6144/production/config_prod_oct${OCT}.toml"
D=/home/yguan/projects/aip-aspuru-ab/yguan/websky
RAW=${D}/catalog_websky_6144_prod_oct${OCT}.pksc
AM=${D}/catalog_websky_6144_prod_oct${OCT}_AM.pksc

echo "Node: $(hostname)  Octant: ${OCT}  obs=(${OBSX},${OBSY},${OBSZ})  Start: $(date)"

echo "=== [1/2] halo finding, coarse_factor=32, z_max=4.5 ==="
julia --project=validation --threads=${JULIA_NUM_THREADS} \
    validation/websky_6144/run_gpu_octant.jl "${CONFIG}"

echo "=== [2/2] abundance matching (per-axis observer, z_max=4.5) ==="
julia --project=validation -t ${SLURM_CPUS_PER_TASK} \
    validation/websky_6144/apply_abundance_match.jl "${RAW}" "${AM}" \
    "${OBSX},${OBSY},${OBSZ}" 4.5

echo "End: $(date)"
