#!/bin/bash
# Full-sky AM rerun of the v3 catalogs ("v3fs"): one table from all 8 raw octants (fsky=1, as in
# the reference Websky procedure) with tail_N (frozen correction above the top TAIL_N halos per
# z-bin), applied per octant, then repaint (halo + CIB) and the full-sky analysis.
# Field maps don't depend on AM and are symlinked from fieldmaps_v3. New outputs live on scratch;
# /project/.../{halomaps,cibmaps,fullsky}_v3fs are symlinks to scratch dirs (project quota).
set -euo pipefail
TAIL_N="${TAIL_N:-10}"
REPO=/home/yguan/work/PeakPatch.jl
S=/home/yguan/scratch/websky_6144; C=${S}/catalogs_v3; LOGS=${S}/logs; SUB=${S}/slurm_v3
D=/home/yguan/projects/aip-aspuru-ab/yguan/websky
OCTS=(000 001 010 011 100 101 110 111)
mkdir -p "${SUB}" "${S}/v3fs/halomaps" "${S}/v3fs/cibmaps" "${S}/v3fs/fullsky" "${D}/fieldmaps_v3fs"
for d in halomaps cibmaps fullsky; do [ -e "${D}/${d}_v3fs" ] || ln -s "${S}/v3fs/${d}" "${D}/${d}_v3fs"; done
for f in "${D}"/fieldmaps_v3/*_v3_oct*_nside4096.fits; do
    b=$(basename "$f"); ln -sf "$f" "${D}/fieldmaps_v3fs/${b/_v3_oct/_v3fs_oct}"; done
cp "${REPO}"/validation/websky_6144/production/run_paint.slurm "${REPO}"/validation/websky_6144/production/run_cib_prod.slurm \
   "${REPO}"/validation/paper/run_jl.slurm "${SUB}/"
cd "${SUB}"
TABLE=${C}/am_table_fullsky_tail${TAIL_N}.txt
CATS=""; for o in "${OCTS[@]}"; do CATS="${CATS} ${C}/catalog_websky_6144_v3_oct${o}.pksc"; done
AMS=validation/websky_6144/apply_abundance_match_fullsky.jl
# SARGS has spaces: pass SCRIPT/SARGS through the environment (--export=ALL), not --export=K=V
tj=$(SCRIPT=${AMS} SARGS="table ${TABLE} ${TAIL_N} 4.5${CATS}" sbatch --parsable --export=ALL \
     --job-name=v3fs-amtable --mem=120G --time=2:55:00 --output="${LOGS}/v3fs_amtable_%j.out" run_jl.slurm)
echo "table: ${tj}"
deps=""
for o in "${OCTS[@]}"; do
    AMC=${C}/catalog_websky_6144_v3_oct${o}_AMfs.pksc; T=v3fs_oct${o}_AM
    aj=$(SCRIPT=${AMS} SARGS="apply ${TABLE} ${C}/catalog_websky_6144_v3_oct${o}.pksc ${AMC} ${o}" sbatch --parsable \
         --export=ALL --dependency=afterok:${tj} --job-name="v3fs-am-${o}" --mem=160G --time=1:30:00 \
         --output="${LOGS}/v3fs_am_oct${o}_%j.out" run_jl.slurm)
    pj=$(sbatch --parsable --dependency=afterok:${aj} --job-name="v3fs-paint-${o}" \
         --export=ALL,OCT=${o},CAT=${AMC},OUTD=${D}/halomaps_v3fs,TAG=${T} \
         --output="${LOGS}/v3fs_paint_oct${o}_%j.out" --error="${LOGS}/v3fs_paint_oct${o}_%j.err" run_paint.slurm)
    cj=$(sbatch --parsable --dependency=afterok:${aj} --job-name="v3fs-cib-${o}" \
         --export=ALL,OCT=${o},CAT=${AMC},OUTD=${D}/cibmaps_v3fs,TAG=${T} \
         --output="${LOGS}/v3fs_cib_oct${o}_%j.out" --error="${LOGS}/v3fs_cib_oct${o}_%j.err" run_cib_prod.slurm)
    echo "oct${o}: am ${aj}  paint ${pj}  cib ${cj}"
    deps="${deps}:${pj}:${cj}"
done
sed 's/v3-analysis/v3fs-analysis/; s/CAMPAIGN=v3$/CAMPAIGN=v3fs/' "${REPO}"/validation/paper/v3_analysis.slurm > v3fs_analysis.slurm
an=$(sbatch --parsable --dependency=afterok${deps} --output="${LOGS}/v3fs_analysis_%j.out" v3fs_analysis.slurm)
echo "analysis: ${an}"
