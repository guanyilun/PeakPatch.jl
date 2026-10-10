#!/bin/bash
#SBATCH --job-name=repro-websky
#SBATCH --partition=genx
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=48G
#SBATCH --time=3:00:00
# Reproduce Websky science results on Websky and Huntian v0.1 (validation/paper/repro/):
#   rodriguez2026.jl  Table III, Fig. 2, Fig. 3 of arXiv:2509.03458   (streams both catalogs, ~155 GB)
#   lee2024.jl        Fig. 5 band statistics of arXiv:2304.07283      (6 CIB maps)
# A 200k-halo smoke run of rodriguez2026.jl and the short lee2024.jl go first, so a bug fails early.
# Memory is not measured: maps are a few GB; 48G is a guess. Submit from anywhere:
#   sbatch --chdir=/mnt/home/yguan/work/uoft/peakpatch/PeakPatch.jl validation/paper/repro/run_repro.slurm
set -euo pipefail
source /mnt/home/yguan/work/uoft/peakpatch/PeakPatch.jl/validation/websky_6144/production_v5/env.sh
cd "${REPO}"
J="julia --project=validation -t ${JULIA_NUM_THREADS}"
echo "=== smoke $(date)"; MAXH=200000 $J validation/paper/repro/rodriguez2026.jl
echo "=== lee2024 $(date)"; $J validation/paper/repro/lee2024.jl
echo "=== rodriguez2026 $(date)"; $J validation/paper/repro/rodriguez2026.jl
echo "=== end $(date)"
