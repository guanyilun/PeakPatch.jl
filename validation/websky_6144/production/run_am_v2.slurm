#!/bin/bash
#SBATCH --job-name=prod-amv2
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b1
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=480G
#SBATCH --time=1:30:00
# Re-run ONLY abundance matching on a frozen production octant with the fixed
# build_abundance_table (top-halo fix, 9f92778; see ../AM_TOPHALO_FIX_2026-09-25.md).
# Same settings as run_catalog.slurm stage 2. Writes *_AMv2.pksc; originals untouched.
# Pass OCT, OBSX, OBSY, OBSZ via --export (see submit_am_v2.sh).
set -euo pipefail
module load julia/1.12.5 cuda/12.6
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"
: "${OCT:?}"; : "${OBSX:?}"; : "${OBSY:?}"; : "${OBSZ:?}"
cd /home/yguan/work/PeakPatch.jl
D=/home/yguan/projects/aip-aspuru-ab/yguan/websky
echo "Node: $(hostname)  Octant: ${OCT}  commit: $(git rev-parse --short HEAD)  Start: $(date)"
julia --project=validation -t ${SLURM_CPUS_PER_TASK} validation/websky_6144/apply_abundance_match.jl \
    ${D}/catalog_websky_6144_prod_oct${OCT}.pksc ${D}/catalog_websky_6144_prod_oct${OCT}_AMv2.pksc \
    "${OBSX},${OBSY},${OBSZ}" 4.5
echo "End: $(date)"
