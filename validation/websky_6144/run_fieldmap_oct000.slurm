#!/bin/bash
#SBATCH --job-name=pp-fieldmap
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b2
#SBATCH --nodes=1
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=32
#SBATCH --mem=480G
#SBATCH --time=11:00:00
#SBATCH --output=/home/yguan/scratch/websky_6144/logs/fieldmap_oct000_%j.out
#SBATCH --error=/home/yguan/scratch/websky_6144/logs/fieldmap_oct000_%j.err

# Field-matter lightcone kappa/mass maps for octant 000 (Phase A,
# docs/field_lightcone_plan.md). Full matter (halo+field in one map), Nside 2048.
set -euo pipefail
module load julia/1.12.5 cuda/12.6
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"
unset CUDA_VISIBLE_DEVICES

cd /home/yguan/work/PeakPatch.jl
CONFIG="validation/websky_6144/config_websky_6144_oct000_finecell.toml"

echo "=== fieldmap GPU smoke test ==="
julia --project=. --threads=8 validation/websky_6144/test_fieldmap_smoke.jl --gpu

echo "=== fieldmap octant 000 (Nside 2048) ==="
echo "Node: $(hostname)"; echo "Start: $(date)"
julia --project=. --threads=${JULIA_NUM_THREADS} \
    validation/websky_6144/run_fieldmap_octant.jl "${CONFIG}" 2048

echo "End: $(date)"
