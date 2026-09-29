#!/bin/bash
#SBATCH --job-name=v4-cat
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b2
#SBATCH --nodes=1
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=32
#SBATCH --mem=480G
#SBATCH --time=8:00:00
# Production-v4: one octant catalog (v3 + volume reduction + per-tile threshold). Halo finding,
# merge and finalize only; AM is done afterwards with one full-sky table (submit_all.sh).
# Pass OCT via --export. Expected ~2.4 h (v3: 142-153 min for this step).
set -euo pipefail
module load julia/1.12.5 cuda/12.6
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"
unset CUDA_VISIBLE_DEVICES
: "${OCT:?}"
cd /home/yguan/work/PeakPatch.jl
CONFIG="validation/websky_6144/production_v4/config_v4_oct${OCT}.toml"
echo "Node: $(hostname)  Octant: ${OCT}  commit: $(git rev-parse --short HEAD)  Start: $(date)"
julia --project=validation --threads=${JULIA_NUM_THREADS} validation/websky_6144/run_gpu_octant.jl "${CONFIG}"
echo "End: $(date)"
