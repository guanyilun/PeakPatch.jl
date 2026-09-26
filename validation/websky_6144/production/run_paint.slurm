#!/bin/bash
#SBATCH --job-name=prod-paint
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b1
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=2:55:00
# CPU halo painting (κ/tSZ/kSZ, Nside 4096) of one frozen-campaign octant from its AMv2
# catalog (paint_octant.jl). The GPU is only the smallest schedulable unit here.
# Pass OCT via --export (see submit_paint.sh); optional PRODS (default kappa,tsz,ksz).
set -euo pipefail
module load julia/1.12.5
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"
: "${OCT:?}"
PRODS="${PRODS:-kappa,tsz,ksz}"
cd /home/yguan/work/PeakPatch.jl
D=/home/yguan/projects/aip-aspuru-ab/yguan/websky
echo "Node: $(hostname)  Octant: ${OCT}  commit: $(git rev-parse --short HEAD)  Start: $(date)"
CAT="${CAT:-${D}/catalog_websky_6144_prod_oct${OCT}_AMv2.pksc}"
OUTD="${OUTD:-${D}/halomaps_prod}"; TAG="${TAG:-prod_oct${OCT}_AMv2}"
julia --project=validation -t ${SLURM_CPUS_PER_TASK} validation/websky_6144/production/paint_octant.jl \
    ${OCT} ${CAT} ${OUTD} 4096 ${PRODS} ${TAG}
echo "End: $(date)"
