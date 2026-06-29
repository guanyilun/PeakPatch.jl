#!/bin/bash
#SBATCH --job-name=pp-am-finecell
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b3
#SBATCH --cpus-per-task=16
#SBATCH --mem=480G
#SBATCH --time=06:00:00
#SBATCH --output=/home/yguan/scratch/websky_6144/logs/am_finecell_%j.out
#SBATCH --error=/home/yguan/scratch/websky_6144/logs/am_finecell_%j.err

# Apply abundance matching to the factor-h-corrected FINECELL octant (197.5M halos),
# then compare to Websky across the full mass range. This is the end-to-end confirmation:
# finecell (fixes M<3e12 resolution floor) + AM (fixes M>3e12 mass calibration) -> full match.
# CPU-only; needs big memory (loads 197.5M-halo / 26GB catalog + AM copy). obs = -2618 (box 5236 Mpc/h).
set -euo pipefail
module load julia/1.12.5
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"

REPO=/home/yguan/work/PeakPatch.jl
D=/home/yguan/projects/aip-aspuru-ab/yguan/websky
cd ${REPO}

FINE=${D}/catalog_websky_6144_oct000_finecell.pksc
FINE_AM=${D}/catalog_websky_6144_oct000_finecell_AM.pksc

echo "=== [1/3] apply AM to finecell octant (obs -2618) ==="
julia --project=. -t ${SLURM_CPUS_PER_TASK} \
    validation/websky_6144/apply_abundance_match.jl "${FINE}" "${FINE_AM}" -2618

echo ""; echo "=== [2/3] cumulative N(>M) per deg²: finecell RAW + AM vs Websky ==="
julia --project=. validation/websky_6144/compare_am_vs_websky.jl "${FINE}" "${FINE_AM}"

echo ""; echo "=== [3/3] low-mass MF shape of the AM'd finecell catalog ==="
julia --project=. validation/websky_6144/compare_lowmass_mf.jl "${FINE_AM}" 0.85221

echo ""; echo "=== DONE ==="
