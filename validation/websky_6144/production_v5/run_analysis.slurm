#!/bin/bash
#SBATCH --job-name=v5fs-analysis
#SBATCH --partition=genx
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=1:30:00
# v5fs full-sky pipeline: assemble -> seam check -> auto spectra -> cross spectra -> theory, and
# the high-mass AM tail per octant. Same steps as validation/paper/v3_analysis.slurm (v4fs: ~26 min,
# plus ~9 min of tail checks). v5 (job 7206290): MaxRSS 35.7 GiB through the cross spectra. Pass REPO via --export.
set -euo pipefail
: "${REPO:?}"
source "${REPO}/validation/websky_6144/production_v5/env.sh"
export CAMPAIGN=v5fs
cd "${REPO}"
J="julia --project=validation -t ${JULIA_NUM_THREADS}"
# STEPS (optional, space-separated) runs a subset, e.g. STEPS="theory tails"; see rerun_analysis.sh
STEPS="${STEPS:-assemble seams autos cross theory tails}"
for step in ${STEPS}; do
    echo "=== ${step} $(date)"
    case ${step} in
        assemble) $J validation/paper/assemble_fullsky.jl ;;
        seams)    $J validation/paper/diag_holes_tsz.jl ;;
        autos)    bash validation/paper/run_fullsky_compare.sh tsz ksz isw kappa cib ;;
        cross)    $J validation/paper/compare_cross.jl ;;
        theory)   $J validation/paper/compare_theory_fullsky.jl ;;   # "before" column = v4fs (fullsky_v4fs/ksz_field_uK_*)
        tails)    for o in 000 001 010 011 100 101 110 111; do TAG=v5fs OCT=$o $J validation/paper/check_am_tail.jl; done ;;
        *)        echo "unknown step ${step}" >&2; exit 1 ;;
    esac
done
echo "=== end $(date)"
