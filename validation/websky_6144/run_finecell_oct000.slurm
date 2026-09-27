#!/bin/bash
#SBATCH --job-name=pp6144-finecell
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b3
#SBATCH --nodes=1
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=32
#SBATCH --mem=480G
#SBATCH --time=23:00:00
#SBATCH --output=/home/yguan/scratch/websky_6144/logs/finecell_oct000_%j.out
#SBATCH --error=/home/yguan/scratch/websky_6144/logs/finecell_oct000_%j.err

# Factor-h-CORRECTED resolution re-run of Websky octant 000.
# box = 5236 Mpc/h (= 7700 Mpc), cellsize 0.852 Mpc/h, finecell filter bank.
# Confirms the cellsize root cause: raw low-mass floor should drop ~3e12 -> ~1e12.
set -euo pipefail
module load julia/1.12.5 cuda/12.6
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"
unset CUDA_VISIBLE_DEVICES

cd /home/yguan/work/PeakPatch.jl
CONFIG="validation/websky_6144/config_websky_6144_oct000_finecell.toml"

echo "=== Websky 6144^3 finecell (factor-h corrected) — octant 000 ==="
echo "Config: ${CONFIG}"; echo "Node: $(hostname)"; echo "Start: $(date)"

julia --project=validation --threads=${JULIA_NUM_THREADS} \
    validation/websky_6144/run_gpu_octant.jl "${CONFIG}"

echo "=== raw low-mass MF check (floor should reach ~1e12) ==="
julia --project=validation validation/websky_6144/compare_lowmass_mf.jl \
    /home/yguan/projects/aip-aspuru-ab/yguan/websky/catalog_websky_6144_oct000_finecell.pksc 0.85221

echo "End: $(date)"
