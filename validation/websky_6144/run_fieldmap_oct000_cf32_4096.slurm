#!/bin/bash
#SBATCH --job-name=pp-field-cf32
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b1
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=2:55:00
#SBATCH --output=/home/yguan/scratch/websky_6144/logs/fieldmap_oct000_cf32_4096_%j.out
#SBATCH --error=/home/yguan/scratch/websky_6144/logs/fieldmap_oct000_cf32_4096_%j.err

# kSZ low-ell Doppler diagnosis: rerun the FULL-matter octant field maps with
# coarse_factor=16 (coarse grid 256^3, Nyquist 0.154 h/Mpc, vs production 64^3 /
# 0.038 h/Mpc). Hypothesis: the k~0.03-0.1 coarse<->fine transition band carries
# excess/decoherent VELOCITY power -> 2-3x kSZ Doppler excess at ell<350 (measured
# octant-wide vs ksz.fits). Density/kappa should be nearly unchanged (control).
set -euo pipefail
module load julia/1.12.5 cuda/12.6
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"
unset CUDA_VISIBLE_DEVICES

cd /home/yguan/work/PeakPatch.jl
CONFIG="validation/websky_6144/config_websky_6144_oct000_finecell_cf32.toml"

echo "Node: $(hostname)"; echo "Start: $(date)"
julia --project=. --threads=${JULIA_NUM_THREADS} \
    validation/websky_6144/run_fieldmap_octant.jl "${CONFIG}" 4096 \
    /home/yguan/scratch/websky_6144/fieldmaps gpu_paint

echo "End: $(date)"
