#!/bin/bash
#SBATCH --job-name=pp-prod-field
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b1
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=2:55:00

# FROZEN production campaign (paper): full-matter field maps for one octant —
# kernels kappa/mass/tau/ksz/isw, Nside 4096, GPU pixelization, cf32, z_max=4.5.
# NO halo exclusion (the production kappa/kSZ composites use compensated halo
# painting on top of the full-matter field — TIER_B_SUMMARY recipes).
# Independent of the catalog job. Pass OCT via sbatch --export (submit_all.sh).
# Timing anchor: cf32 4-kernel 4096 run took ~1.3h on one L40S.
set -euo pipefail
module load julia/1.12.5 cuda/12.6
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"
unset CUDA_VISIBLE_DEVICES

: "${OCT:?need OCT=zyx octant bits}"

cd /home/yguan/work/PeakPatch.jl
CONFIG="validation/websky_6144/production/config_prod_oct${OCT}.toml"
OUTDIR=/home/yguan/scratch/websky_6144/fieldmaps_prod

echo "Node: $(hostname)  Octant: ${OCT}  Start: $(date)"
julia --project=. --threads=${JULIA_NUM_THREADS} \
    validation/websky_6144/run_fieldmap_octant.jl "${CONFIG}" 4096 "${OUTDIR}" gpu_paint

echo "End: $(date)"
