#!/bin/bash
# Submit the production-v3 campaign: 8 octants x (catalog+AM job, fieldmap job).
# sbatch cannot run from /home on Killarney — scripts are copied to scratch first.
# Usage: bash submit_all.sh [catalog|fieldmap|all] [OCT ...]   (default: all, all octants)
# "catalog"/"all" also chain, per octant, halo painting (paint_octant.jl) and CIB
# (paint_cib.jl) on the v3 AM catalog with --dependency=afterok on the catalog job.
set -euo pipefail
STAGE="${1:-all}"; shift || true
OCTS=("$@"); [ ${#OCTS[@]} -eq 0 ] && OCTS=(000 001 010 011 100 101 110 111)
REPO=/home/yguan/work/PeakPatch.jl/validation/websky_6144/production_v3
SUB=/home/yguan/scratch/websky_6144/slurm_v3; LOGS=/home/yguan/scratch/websky_6144/logs
mkdir -p "${SUB}" "${LOGS}" /home/yguan/scratch/websky_6144/catalogs_v3
cp "${REPO}"/run_catalog.slurm "${REPO}"/run_fieldmap.slurm "${REPO}"/../production/run_paint.slurm \
   "${REPO}"/../production/run_cib_prod.slurm "${SUB}/"; cd "${SUB}"
D=/home/yguan/projects/aip-aspuru-ab/yguan/websky
for OCT in "${OCTS[@]}"; do
    Z=${OCT:0:1}; Y=${OCT:1:1}; X=${OCT:2:1}
    [ "$X" = "1" ] && OBSX=2618.0 || OBSX=-2618.0
    [ "$Y" = "1" ] && OBSY=2618.0 || OBSY=-2618.0
    [ "$Z" = "1" ] && OBSZ=2618.0 || OBSZ=-2618.0
    if [ "${STAGE}" = "all" ] || [ "${STAGE}" = "catalog" ]; then
        jid=$(sbatch --parsable --job-name="v3-cat-${OCT}" \
            --export=ALL,OCT=${OCT},OBSX=${OBSX},OBSY=${OBSY},OBSZ=${OBSZ} \
            --output="${LOGS}/v3_cat_oct${OCT}_%j.out" --error="${LOGS}/v3_cat_oct${OCT}_%j.err" run_catalog.slurm)
        echo "oct${OCT} catalog:  job ${jid}  obs=(${OBSX},${OBSY},${OBSZ})"
        C=${D}/catalog_websky_6144_v3_oct${OCT}_AM.pksc; T=v3_oct${OCT}_AM
        pj=$(sbatch --parsable --dependency=afterok:${jid} --job-name="v3-paint-${OCT}" \
            --export=ALL,OCT=${OCT},CAT=${C},OUTD=${D}/halomaps_v3,TAG=${T} \
            --output="${LOGS}/v3_paint_oct${OCT}_%j.out" --error="${LOGS}/v3_paint_oct${OCT}_%j.err" run_paint.slurm)
        cj=$(sbatch --parsable --dependency=afterok:${jid} --job-name="v3-cib-${OCT}" \
            --export=ALL,OCT=${OCT},CAT=${C},OUTD=${D}/cibmaps_v3,TAG=${T} \
            --output="${LOGS}/v3_cib_oct${OCT}_%j.out" --error="${LOGS}/v3_cib_oct${OCT}_%j.err" run_cib_prod.slurm)
        echo "oct${OCT} paint:    job ${pj} (after ${jid});  cib: job ${cj}"
    fi
    if [ "${STAGE}" = "all" ] || [ "${STAGE}" = "fieldmap" ]; then
        jid=$(sbatch --parsable --job-name="v3-fld-${OCT}" --export=ALL,OCT=${OCT} \
            --output="${LOGS}/v3_fld_oct${OCT}_%j.out" --error="${LOGS}/v3_fld_oct${OCT}_%j.err" run_fieldmap.slurm)
        echo "oct${OCT} fieldmap: job ${jid}"
    fi
done
