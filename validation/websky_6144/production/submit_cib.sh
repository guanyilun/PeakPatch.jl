#!/bin/bash
# Submit production CIB painting for the given octants (default: all 8), e.g.
#   bash submit_paint.sh 000        |   bash submit_paint.sh
set -euo pipefail
REPO=/home/yguan/work/PeakPatch.jl/validation/websky_6144/production
SUB=/home/yguan/scratch/websky_6144/slurm_prod; LOGS=/home/yguan/scratch/websky_6144/logs
cp "${REPO}/run_cib_prod.slurm" "${SUB}/"; cd "${SUB}"
OCTS=("$@"); [ ${#OCTS[@]} -eq 0 ] && OCTS=(000 001 010 011 100 101 110 111)
for OCT in "${OCTS[@]}"; do
    jid=$(sbatch --parsable --job-name="cib-${OCT}" --export=ALL,OCT=${OCT} \
        --output="${LOGS}/prod_cib_oct${OCT}_%j.out" --error="${LOGS}/prod_cib_oct${OCT}_%j.err" run_cib_prod.slurm)
    echo "oct${OCT}: job ${jid}"
done
