#!/bin/bash
#SBATCH --job-name=pp-field4096
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b2
#SBATCH --nodes=1
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=32
#SBATCH --mem=480G
#SBATCH --time=11:00:00
#SBATCH --output=/home/yguan/scratch/websky_6144/logs/fieldmap_oct000_field4096_%j.out
#SBATCH --error=/home/yguan/scratch/websky_6144/logs/fieldmap_oct000_field4096_%j.err

# Field-ONLY octant maps at Nside 4096 with on-device painting (gpu_paint -> subdiv 5):
# removes the Nside-2048 pixel window + subdiv-3 granularity from the composite's
# high-ell deficit diagnosis (COMPOSITE_KAPPA_2026-07-18.md follow-up plan, step 1).
set -euo pipefail
module load julia/1.12.5 cuda/12.6
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"
unset CUDA_VISIBLE_DEVICES

cd /home/yguan/work/PeakPatch.jl
CONFIG="validation/websky_6144/config_websky_6144_oct000_finecell.toml"
EXCL="/home/yguan/projects/aip-aspuru-ab/yguan/websky/catalog_websky_6144_oct000_finecell_AM.pksc"

echo "=== fieldmap GPU smoke test (incl. gpu_paint x exclusion) ==="
julia --project=validation --threads=8 validation/websky_6144/test_fieldmap_smoke.jl --gpu

echo "=== fieldmap octant 000, field-only, Nside 4096, gpu_paint/subdiv 5 ==="
echo "Node: $(hostname)"; echo "Start: $(date)"
julia --project=validation --threads=${JULIA_NUM_THREADS} \
    validation/websky_6144/run_fieldmap_octant.jl "${CONFIG}" 4096 \
    /home/yguan/scratch/websky_6144/fieldmaps gpu_paint "exclude=${EXCL}"

echo "End: $(date)"
