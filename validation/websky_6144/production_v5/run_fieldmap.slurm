#!/bin/bash
#SBATCH --job-name=v5-fld
#SBATCH --partition=gpu
#SBATCH --constraint=a100-80gb&rocky9
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2:30:00
# Production-v5 full-matter field maps (kappa/mass/tau/ksz/isw, Nside 4096, GPU pixelization,
# cf32, z_max 4.5) for one octant. Not reusable from v3/v4: gaussian_split also changes the field
# (FieldMap.jl, c880c24). Halo-independent, so it runs in parallel with the catalogs.
# Sizing (v5 on Rusty A100, jobs 7206250+2k): 1.26-1.31 h, MaxRSS 35-40 GiB (v3 on L40S: ~1.35 h).
# Pass OCT and REPO via --export.
# rocky9: env.sh loads modules from the Rocky 9 tree; the 3 Rocky 8 A100 nodes (workergpu038-040) are untested.
set -euo pipefail
: "${REPO:?}" "${OCT:?}"
source "${REPO}/validation/websky_6144/production_v5/env.sh"
cd "${REPO}"
CONFIG="validation/websky_6144/production_v5/config_v5_oct${OCT}.toml"
OUTDIR="${WS_ROOT}/fieldmaps_v5fs"
mkdir -p "${OUTDIR}"
echo "Node: $(hostname)  Octant: ${OCT}  commit: $(git rev-parse --short HEAD)  Start: $(date)"
julia --project=validation --threads=${JULIA_NUM_THREADS} \
    validation/websky_6144/run_fieldmap_octant.jl "${CONFIG}" 4096 "${OUTDIR}" gpu_paint
# Files are tagged by config name (v5_octZYX); the analysis scripts expect the campaign tag v5fs.
cd "${OUTDIR}"
for f in *_v5_oct${OCT}_nside4096.fits; do ln -sf "$f" "${f/_v5_oct/_v5fs_oct}"; done
echo "End: $(date)"
