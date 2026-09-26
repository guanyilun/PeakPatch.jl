#!/bin/bash
#SBATCH --job-name=prod-cib
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b1
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=180G
#SBATCH --time=2:55:00
# CIB (XGPaint CIB_Planck2013, 6 Planck frequencies, Websky completeness cut) of one
# frozen-campaign octant from its AMv2 catalog (paint_cib.jl).
set -euo pipefail
module load julia/1.12.5
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"
: "${OCT:?}"
cd /home/yguan/work/PeakPatch.jl
D=/home/yguan/projects/aip-aspuru-ab/yguan/websky
echo "Node: $(hostname)  Octant: ${OCT}  commit: $(git rev-parse --short HEAD)  Start: $(date)"
CAT="${CAT:-${D}/catalog_websky_6144_prod_oct${OCT}_AMv2.pksc}"
OUTD="${OUTD:-${D}/cibmaps_prod}"; TAG="${TAG:-prod_oct${OCT}_AMv2}"
julia --project=/home/yguan/work/XGPaint.jl -t ${SLURM_CPUS_PER_TASK} validation/websky_6144/production/paint_cib.jl \
    ${OCT} ${CAT} ${OUTD} 4096 100,143,217,353,545,857 wcut ${TAG}
echo "End: $(date)"
