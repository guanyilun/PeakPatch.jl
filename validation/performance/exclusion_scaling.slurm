#!/bin/bash
#SBATCH --job-name=pp-excl
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b1
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --time=2:30:00
#SBATCH --output=/home/yguan/scratch/websky_6144/logs/bench_exclusion_%j.out
#SBATCH --error=/home/yguan/scratch/websky_6144/logs/bench_exclusion_%j.err
# CPU-only microbenchmark (the GPU is unused; the b1 GPU slice is the fastest-backfilling
# allocation on Killarney). Single-threaded, like merge_catalog.
set -euo pipefail
module load julia/1.12.5
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"
cd /home/yguan/work/PeakPatch.jl
echo "Node: $(hostname)  Start: $(date)"
julia --project=validation -t 1 validation/performance/exclusion_scaling.jl
echo "End: $(date)"
