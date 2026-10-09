#!/bin/bash
# Production-v5 campaign ("v5fs") on Flatiron Rusty, end to end:
#   8 catalog jobs (gpu, 4×A100) -> one full-sky AM table from the 8 v5 raw catalogs (tail_N) ->
#   AM per octant -> halo paint + CIB per octant (genx) -> full-sky analysis + AM tails, and Tier-A.
#   8 field-map jobs (gpu, 1×A100) run alongside the catalogs; the analysis waits for them too.
# 43 jobs in total. The gpu QoS caps a user at 16 GPUs, so catalogs run 4 at a time (~2 rounds).
# Run by hand from a clean, committed tree:  bash validation/websky_6144/production_v5/submit_all.sh
# Every job loads the live working tree: do not edit src/ until the catalog jobs have started.
set -euo pipefail
TAIL_N="${TAIL_N:-10}"
export REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
source "${REPO}/validation/websky_6144/production_v5/env.sh"
P="${REPO}/validation/websky_6144/production_v5"
C="${WS_CATS}/catalogs_v5"; LOGS="${WS_ROOT}/logs"
OCTS=(000 001 010 011 100 101 110 111)

# run_gpu_octant.jl refuses a dirty tree (provenance); fail here rather than after the queue wait
git -C "${REPO}" diff --quiet HEAD -- || { echo "working tree is dirty: commit first" >&2; exit 1; }
for f in halos_10x10.pksc halo_mass_completion.txt kap_lt4.5.fits tsz_2048.fits; do
    [ -f "${WS_REF}/${f}" ] || { echo "missing ${WS_REF}/${f}" >&2; exit 1; }
done
mkdir -p "${C}" "${LOGS}" "${WS_ROOT}"/{fieldmaps,halomaps,cibmaps,fullsky}_v5fs
echo "commit $(git -C "${REPO}" rev-parse --short HEAD)  $(date)" | tee -a "${LOGS}/v5_submit.txt"

cats=""; CATS=""; flds=""
for o in "${OCTS[@]}"; do
    j=$(sbatch --parsable --job-name="v5-cat-${o}" --export=ALL,OCT=${o} \
        --output="${LOGS}/v5_cat_oct${o}_%j.out" --error="${LOGS}/v5_cat_oct${o}_%j.err" "${P}/run_catalog.slurm")
    f=$(sbatch --parsable --job-name="v5-fld-${o}" --export=ALL,OCT=${o} \
        --output="${LOGS}/v5_fld_oct${o}_%j.out" --error="${LOGS}/v5_fld_oct${o}_%j.err" "${P}/run_fieldmap.slurm")
    echo "oct${o}: catalog ${j}  fieldmap ${f}"
    cats="${cats}:${j}"; flds="${flds}:${f}"; CATS="${CATS} ${C}/catalog_websky_6144_v5_oct${o}.pksc"
done

# CPU steps go through run_jl.slurm. SCRIPT/SARGS are set in the environment rather than in
# --export because SARGS contains commas.
AMS=validation/websky_6144/apply_abundance_match_fullsky.jl; TABLE="${C}/am_table_fullsky_tail${TAIL_N}.txt"
tj=$(SCRIPT=${AMS} SARGS="table ${TABLE} ${TAIL_N} 4.5${CATS}" sbatch --parsable --export=ALL \
     --dependency=afterok${cats} --job-name=v5fs-amtable --cpus-per-task=8 --mem=64G --time=1:30:00 \
     --output="${LOGS}/v5fs_amtable_%j.out" "${P}/run_jl.slurm")
echo "AM table: ${tj}"
ams=""; maps=""
for o in "${OCTS[@]}"; do
    AMC="${C}/catalog_websky_6144_v5_oct${o}_AMfs.pksc"; T="v5fs_oct${o}_AM"
    aj=$(SCRIPT=${AMS} SARGS="apply ${TABLE} ${C}/catalog_websky_6144_v5_oct${o}.pksc ${AMC} ${o}" sbatch --parsable \
         --export=ALL --dependency=afterok:${tj} --job-name="v5fs-am-${o}" --cpus-per-task=8 --mem=64G --time=0:30:00 \
         --output="${LOGS}/v5fs_am_oct${o}_%j.out" "${P}/run_jl.slurm")
    pj=$(SCRIPT=validation/websky_6144/production/paint_octant.jl \
         SARGS="${o} ${AMC} ${WS_ROOT}/halomaps_v5fs 4096 kappa,tsz,ksz ${T}" sbatch --parsable --export=ALL \
         --dependency=afterok:${aj} --job-name="v5fs-paint-${o}" --cpus-per-task=32 --mem=48G --time=0:45:00 \
         --output="${LOGS}/v5fs_paint_oct${o}_%j.out" "${P}/run_jl.slurm")
    cj=$(JPROJECT="${XGPAINT}" SCRIPT=validation/websky_6144/production/paint_cib.jl \
         SARGS="${o} ${AMC} ${WS_ROOT}/cibmaps_v5fs 4096 100,143,217,353,545,857 wcut ${T}" sbatch --parsable --export=ALL \
         --dependency=afterok:${aj} --job-name="v5fs-cib-${o}" --cpus-per-task=32 --mem=128G --time=1:30:00 \
         --output="${LOGS}/v5fs_cib_oct${o}_%j.out" "${P}/run_jl.slurm")
    echo "oct${o}: am ${aj}  paint ${pj}  cib ${cj}"
    ams="${ams}:${aj}"; maps="${maps}:${pj}:${cj}"
done
an=$(sbatch --parsable --export=ALL --dependency=afterok${maps}${flds} \
     --output="${LOGS}/v5fs_analysis_%j.out" "${P}/run_analysis.slurm")
ta=$(SCRIPT=validation/paper/tierA_v3fs.jl TIERA_CAMP=v5 sbatch --parsable --export=ALL --dependency=afterok${ams} \
     --job-name=v5fs-tierA --cpus-per-task=16 --mem=96G --time=1:30:00 --output="${LOGS}/v5fs_tierA_%j.out" "${P}/run_jl.slurm")
echo "analysis: ${an}  tierA: ${ta}"
