#!/bin/bash
#SBATCH --job-name=pp-cat-cf32
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b2
#SBATCH --nodes=1
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=32
#SBATCH --mem=480G
#SBATCH --time=10:00:00
#SBATCH --output=/home/yguan/scratch/websky_6144/logs/finecell_oct000_cf32_%j.out
#SBATCH --error=/home/yguan/scratch/websky_6144/logs/finecell_oct000_cf32_%j.err

# Octant catalog regeneration at coarse_factor=32: fixes the large-scale velocity
# decoherence in halo velocities (KSZ_COMPOSITE_2026-07-19.md — the cf4 catalog's
# halo kSZ term is ~7x Websky's at ell~500-1200 from exactly this). New-convention
# output (Eulerian + km/s via finalize_eulerian), then abundance matching.
set -euo pipefail
module load julia/1.12.5 cuda/12.6
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"
unset CUDA_VISIBLE_DEVICES

cd /home/yguan/work/PeakPatch.jl
CONFIG="validation/websky_6144/config_websky_6144_oct000_finecell_cf32.toml"
D=/home/yguan/projects/aip-aspuru-ab/yguan/websky
RAW=${D}/catalog_websky_6144_oct000_finecell_cf32.pksc
AM=${D}/catalog_websky_6144_oct000_finecell_cf32_AM.pksc

echo "=== [1/2] halo finding, coarse_factor=32 ==="
julia --project=. --threads=${JULIA_NUM_THREADS} \
    validation/websky_6144/run_gpu_octant.jl "${CONFIG}"

echo "=== [2/2] abundance matching ==="
julia --project=. -t ${SLURM_CPUS_PER_TASK} \
    validation/websky_6144/apply_abundance_match.jl "${RAW}" "${AM}" -2618

echo "End: $(date)"
