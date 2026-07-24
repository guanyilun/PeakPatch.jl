#!/bin/bash
#SBATCH --job-name=pp-cib
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b1
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=2:30:00
#SBATCH --output=/home/yguan/scratch/websky_6144/logs/cib_%j.out
#SBATCH --error=/home/yguan/scratch/websky_6144/logs/cib_%j.err

# CIB painting via XGPaint (CIB_Planck2013 = Websky CIB model), full AM catalog.
set -euo pipefail
module load julia/1.12.5
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"

cd /home/yguan/work/PeakPatch.jl
FREQ="${1:-545}"
julia --project=/home/yguan/work/XGPaint.jl --threads=${JULIA_NUM_THREADS} \
    validation/websky_6144/run_cib.jl "${FREQ}" 4096 0
