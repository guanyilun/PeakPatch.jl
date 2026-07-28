#!/bin/bash
# Submit the frozen production campaign: 8 octants x (catalog+AM job, fieldmap job).
# sbatch cannot run from /home on Killarney — scripts are copied to scratch first.
# Usage: bash submit_all.sh [catalog|fieldmap|all]   (default: all)
set -euo pipefail

STAGE="${1:-all}"
REPO=/home/yguan/work/PeakPatch.jl/validation/websky_6144/production
SUB=/home/yguan/scratch/websky_6144/slurm_prod
LOGS=/home/yguan/scratch/websky_6144/logs
mkdir -p "${SUB}" "${LOGS}" /home/yguan/scratch/websky_6144/fieldmaps_prod
cp "${REPO}"/run_catalog.slurm "${REPO}"/run_fieldmap.slurm "${SUB}/"
cd "${SUB}"

for OCT in 000 001 010 011 100 101 110 111; do
    Z=${OCT:0:1}; Y=${OCT:1:1}; X=${OCT:2:1}
    [ "$X" = "1" ] && OBSX=2618.0  || OBSX=-2618.0
    [ "$Y" = "1" ] && OBSY=2618.0  || OBSY=-2618.0
    [ "$Z" = "1" ] && OBSZ=2618.0  || OBSZ=-2618.0

    if [ "${STAGE}" = "all" ] || [ "${STAGE}" = "catalog" ]; then
        jid=$(sbatch --parsable \
            --job-name="prod-cat-${OCT}" \
            --export=ALL,OCT=${OCT},OBSX=${OBSX},OBSY=${OBSY},OBSZ=${OBSZ} \
            --output="${LOGS}/prod_cat_oct${OCT}_%j.out" \
            --error="${LOGS}/prod_cat_oct${OCT}_%j.err" \
            run_catalog.slurm)
        echo "oct${OCT} catalog:  job ${jid}  obs=(${OBSX},${OBSY},${OBSZ})"
    fi
    if [ "${STAGE}" = "all" ] || [ "${STAGE}" = "fieldmap" ]; then
        jid=$(sbatch --parsable \
            --job-name="prod-fld-${OCT}" \
            --export=ALL,OCT=${OCT} \
            --output="${LOGS}/prod_fld_oct${OCT}_%j.out" \
            --error="${LOGS}/prod_fld_oct${OCT}_%j.err" \
            run_fieldmap.slurm)
        echo "oct${OCT} fieldmap: job ${jid}"
    fi
done
