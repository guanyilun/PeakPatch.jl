#!/bin/bash
#SBATCH --job-name=v2-analysis
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=160G
#SBATCH --time=2:55:00
# v2 full-sky pipeline: assemble -> seam check -> all auto spectra -> cross spectra
module load julia/1.12.5
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"
export CAMPAIGN=v2
cd /home/yguan/work/PeakPatch.jl
J="julia --project=validation -t ${SLURM_CPUS_PER_TASK}"
echo "=== assemble $(date)"; $J validation/paper/assemble_fullsky.jl
echo "=== seam check $(date)"; $J validation/paper/diag_holes_tsz.jl
echo "=== autos $(date)"; bash validation/paper/run_fullsky_compare.sh tsz ksz isw kappa cib
echo "=== cross $(date)"; $J validation/paper/compare_cross.jl
echo "=== end $(date)"
