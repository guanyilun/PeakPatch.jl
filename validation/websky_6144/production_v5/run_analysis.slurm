#!/bin/bash
#SBATCH --job-name=v5fs-analysis
#SBATCH --partition=genx
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=1:30:00
# v5fs full-sky pipeline: assemble -> seam check -> auto spectra -> cross spectra -> theory, and
# the high-mass AM tail per octant. Same steps as validation/paper/v3_analysis.slurm (v4fs: ~26 min,
# plus ~9 min of tail checks); memory was never measured. Pass REPO via --export.
set -euo pipefail
: "${REPO:?}"
source "${REPO}/validation/websky_6144/production_v5/env.sh"
export CAMPAIGN=v5fs
cd "${REPO}"
J="julia --project=validation -t ${JULIA_NUM_THREADS}"
echo "=== assemble $(date)"; $J validation/paper/assemble_fullsky.jl
echo "=== seam check $(date)"; $J validation/paper/diag_holes_tsz.jl
echo "=== autos $(date)"; bash validation/paper/run_fullsky_compare.sh tsz ksz isw kappa cib
echo "=== cross $(date)"; $J validation/paper/compare_cross.jl
echo "=== theory $(date)"; $J validation/paper/compare_theory_fullsky.jl
echo "=== AM tails $(date)"
for o in 000 001 010 011 100 101 110 111; do TAG=v5fs OCT=$o $J validation/paper/check_am_tail.jl; done
echo "=== end $(date)"
