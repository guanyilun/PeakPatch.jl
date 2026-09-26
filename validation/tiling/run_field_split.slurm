#!/bin/bash
#SBATCH --job-name=tiling-field
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b1
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=2:50:00
#SBATCH --output=/home/yguan/scratch/websky_6144/logs/tiling_field_%j.out
#SBATCH --error=/home/yguan/scratch/websky_6144/logs/tiling_field_%j.err
# Field-level split test (A5b). CPU work only (GPU partition is the only one available).
# Env: FIELD_N (480), FIELD_BLOCK (12), FIELD_NTILES (1,2,4,8), OUTSUB (results subdir).
set -euo pipefail
module load julia/1.12.5 cuda/12.6
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"
export FIELD_FFTW_THREADS=${FIELD_FFTW_THREADS:-1}
cd /home/yguan/work/PeakPatch.jl
OUT=validation/tiling/results/${OUTSUB:-field_N${FIELD_N:-480}_b${FIELD_BLOCK:-12}}
echo "Node: $(hostname) commit: $(git rev-parse --short HEAD) N=${FIELD_N:-480} block=${FIELD_BLOCK:-12} ntiles=${FIELD_NTILES:-1,2,4,8} Start: $(date)"
julia --project=validation -t ${SLURM_CPUS_PER_TASK} validation/tiling/field_split_test.jl $OUT
echo "End: $(date)"
