#!/bin/bash
#SBATCH --job-name=v2-fld
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b1
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=2:55:00
# Production-v2 campaign: full-matter field maps (kappa/mass/tau/ksz/isw, Nside 4096,
# GPU pixelization, cf32, z_max=4.5, no exclusion) for one octant. Pass OCT.
set -euo pipefail
module load julia/1.12.5 cuda/12.6
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"
unset CUDA_VISIBLE_DEVICES
: "${OCT:?}"
cd /home/yguan/work/PeakPatch.jl
CONFIG="validation/websky_6144/production_v2/config_v2_oct${OCT}.toml"
OUTDIR=/home/yguan/projects/aip-aspuru-ab/yguan/websky/fieldmaps_v2
mkdir -p "${OUTDIR}"
echo "Node: $(hostname)  Octant: ${OCT}  commit: $(git rev-parse --short HEAD)  Start: $(date)"
julia --project=validation --threads=${JULIA_NUM_THREADS} \
    validation/websky_6144/run_fieldmap_octant.jl "${CONFIG}" 4096 "${OUTDIR}" gpu_paint
echo "End: $(date)"
