#!/bin/bash
#SBATCH --job-name=pp-field-isw
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b1
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=2:55:00
#SBATCH --output=/home/yguan/scratch/websky_6144/logs/fieldmap_oct000_isw_%j.out
#SBATCH --error=/home/yguan/scratch/websky_6144/logs/fieldmap_oct000_isw_%j.err

# ISW octant map (new :isw kernel) at coarse_factor=32 (potential is dominated by
# the largest scales), Nside 1024 (ISW dies by ell~100-300). GPU smoke test first:
# the pot-convolve (kernel id 5) + pwmask painting GPU paths are new.
set -euo pipefail
module load julia/1.12.5 cuda/12.6
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"
unset CUDA_VISIBLE_DEVICES

cd /home/yguan/work/PeakPatch.jl
echo "=== GPU smoke test (incl. new :isw kernel) ==="
julia --project=. --threads=8 validation/websky_6144/test_fieldmap_smoke.jl --gpu

echo "=== ISW octant, cf32, Nside 1024 ==="
CONFIG="validation/websky_6144/config_websky_6144_oct000_finecell_cf32.toml"
julia --project=. --threads=${JULIA_NUM_THREADS} \
    validation/websky_6144/run_fieldmap_octant.jl "${CONFIG}" 1024 \
    /home/yguan/scratch/websky_6144/fieldmaps gpu_paint
echo "End: $(date)"
