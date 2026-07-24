#!/bin/bash
#SBATCH --job-name=pp-tsz
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b1
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1:30:00
#SBATCH --output=/home/yguan/scratch/websky_6144/logs/tsz_%j.out
#SBATCH --error=/home/yguan/scratch/websky_6144/logs/tsz_%j.err

# CPU-only tSZ halo painting vs tsz_2048.fits; minimal 1-GPU slice for backfill.
set -euo pipefail
module load julia/1.12.5
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"

cd /home/yguan/work/PeakPatch.jl
julia --project=. --threads=${JULIA_NUM_THREADS} validation/websky_6144/compare_tsz.jl
