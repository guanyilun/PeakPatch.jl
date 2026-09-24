#!/bin/bash
#SBATCH --job-name=rng-null
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b1
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=2:30:00
#SBATCH --output=/home/yguan/scratch/websky_6144/rng_test/rng_%j.out
#SBATCH --error=/home/yguan/scratch/websky_6144/rng_test/rng_%j.err
# Large-scale RNG null test (see rng_large_scale_test.jl). CPU-only work; the GPU
# is requested only because every Killarney partition is a GPU partition.
set -euo pipefail
module load julia/1.12.5 cuda/12.6
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"
cd /home/yguan/work/PeakPatch.jl
echo "Node: $(hostname)  Start: $(date)"
RNG_PART=${RNG_PART:-AB} RNG_FFTW_THREADS=${RNG_FFTW_THREADS:-32} julia --project=validation -t ${SLURM_CPUS_PER_TASK} validation/rng/rng_large_scale_test.jl \
    /home/yguan/scratch/websky_6144/rng_test ${K:-40}
echo "End: $(date)"
