#!/bin/bash
# Submit the AM-v2 rerun (fixed top-halo mapping) for all 8 frozen octants.
set -euo pipefail
REPO=/home/yguan/work/PeakPatch.jl/validation/websky_6144/production
SUB=/home/yguan/scratch/websky_6144/slurm_prod; LOGS=/home/yguan/scratch/websky_6144/logs
cp "${REPO}/run_am_v2.slurm" "${SUB}/"; cd "${SUB}"
for OCT in 000 001 010 011 100 101 110 111; do
    Z=${OCT:0:1}; Y=${OCT:1:1}; X=${OCT:2:1}
    [ "$X" = "1" ] && OBSX=2618.0 || OBSX=-2618.0
    [ "$Y" = "1" ] && OBSY=2618.0 || OBSY=-2618.0
    [ "$Z" = "1" ] && OBSZ=2618.0 || OBSZ=-2618.0
    jid=$(sbatch --parsable --job-name="amv2-${OCT}" \
        --export=ALL,OCT=${OCT},OBSX=${OBSX},OBSY=${OBSY},OBSZ=${OBSZ} \
        --output="${LOGS}/prod_amv2_oct${OCT}_%j.out" --error="${LOGS}/prod_amv2_oct${OCT}_%j.err" run_am_v2.slurm)
    echo "oct${OCT}: job ${jid} obs=(${OBSX},${OBSY},${OBSZ})"
done
