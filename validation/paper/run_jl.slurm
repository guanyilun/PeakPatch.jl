#!/bin/bash
#SBATCH --job-name=paper-jl
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=2:00:00
# Generic runner: sbatch --export=ALL,SCRIPT=validation/paper/x.jl,SARGS="..." run_jl.slurm
module load julia/1.12.5
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"
cd /home/yguan/work/PeakPatch.jl
echo "Node: $(hostname) commit: $(git rev-parse --short HEAD) Start: $(date)"
julia --project=validation -t ${SLURM_CPUS_PER_TASK} ${SCRIPT} ${SARGS:-}
echo "End: $(date)"
