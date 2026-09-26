#!/bin/bash
#SBATCH --job-name=tiling-cat
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b1
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=2:50:00
#SBATCH --output=/home/yguan/scratch/websky_6144/logs/tiling_cat_%j.out
#SBATCH --error=/home/yguan/scratch/websky_6144/logs/tiling_cat_%j.err
# Catalog-level tiling-independence test (A5a) on one L40S. Env: CAT_N (480), CAT_NTILES (1,2,4,8).
set -euo pipefail
module load julia/1.12.5 cuda/12.6
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"
cd /home/yguan/work/PeakPatch.jl
OUT=validation/tiling/results/${OUTSUB:-catalog_N${CAT_N:-480}}
echo "Node: $(hostname) commit: $(git rev-parse --short HEAD) N=${CAT_N:-480} ntiles=${CAT_NTILES:-1,2,4,8} Start: $(date)"
julia --project=validation -t ${SLURM_CPUS_PER_TASK} validation/tiling/catalog_tiling_test.jl $OUT
echo "End: $(date)"
