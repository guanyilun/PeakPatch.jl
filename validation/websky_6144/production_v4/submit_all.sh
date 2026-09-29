#!/bin/bash
# Production-v4 campaign ("v4fs"), end to end:
#   8 catalog jobs (v4 configs: volume reduction + per-tile threshold) -> one full-sky AM table
#   (tail_N) from all 8 raw octants -> AM applied per octant -> halo paint + CIB per octant ->
#   full-sky analysis (CAMPAIGN=v4fs), tail check and Tier-A (TIERA_CAMP=v4).
# Field maps are halo-independent: symlinked from fieldmaps_v3. All new outputs live on scratch;
# /project/.../{halomaps,cibmaps,fullsky,fieldmaps}_v4fs are symlinks / symlink dirs (quota).
set -euo pipefail
TAIL_N="${TAIL_N:-10}"
REPO=/home/yguan/work/PeakPatch.jl
S=/home/yguan/scratch/websky_6144; C=${S}/catalogs_v4; LOGS=${S}/logs; SUB=${S}/slurm_v4
D=/home/yguan/projects/aip-aspuru-ab/yguan/websky
OCTS=(000 001 010 011 100 101 110 111)
mkdir -p "${SUB}" "${C}" "${S}/v4fs/halomaps" "${S}/v4fs/cibmaps" "${S}/v4fs/fullsky" "${D}/fieldmaps_v4fs"
for d in halomaps cibmaps fullsky; do [ -e "${D}/${d}_v4fs" ] || ln -s "${S}/v4fs/${d}" "${D}/${d}_v4fs"; done
for f in "${D}"/fieldmaps_v3/*_v3_oct*_nside4096.fits; do b=$(basename "$f"); ln -sf "$f" "${D}/fieldmaps_v4fs/${b/_v3_oct/_v4fs_oct}"; done
cp "${REPO}"/validation/websky_6144/production_v4/run_catalog.slurm "${REPO}"/validation/websky_6144/production/run_paint.slurm \
   "${REPO}"/validation/websky_6144/production/run_cib_prod.slurm "${REPO}"/validation/paper/run_jl.slurm "${SUB}/"
cd "${SUB}"
cats=""; CATS=""
for o in "${OCTS[@]}"; do
    j=$(sbatch --parsable --job-name="v4-cat-${o}" --export=ALL,OCT=${o} \
        --output="${LOGS}/v4_cat_oct${o}_%j.out" --error="${LOGS}/v4_cat_oct${o}_%j.err" run_catalog.slurm)
    echo "oct${o} catalog: ${j}"; cats="${cats}:${j}"; CATS="${CATS} ${C}/catalog_websky_6144_v4_oct${o}.pksc"
done
AMS=validation/websky_6144/apply_abundance_match_fullsky.jl; TABLE=${C}/am_table_fullsky_tail${TAIL_N}.txt
tj=$(SCRIPT=${AMS} SARGS="table ${TABLE} ${TAIL_N} 4.5${CATS}" sbatch --parsable --export=ALL --dependency=afterok${cats} \
     --job-name=v4fs-amtable --mem=120G --time=2:55:00 --output="${LOGS}/v4fs_amtable_%j.out" run_jl.slurm)
echo "AM table: ${tj}"
deps=""
for o in "${OCTS[@]}"; do
    AMC=${C}/catalog_websky_6144_v4_oct${o}_AMfs.pksc; T=v4fs_oct${o}_AM
    aj=$(SCRIPT=${AMS} SARGS="apply ${TABLE} ${C}/catalog_websky_6144_v4_oct${o}.pksc ${AMC} ${o}" sbatch --parsable \
         --export=ALL --dependency=afterok:${tj} --job-name="v4fs-am-${o}" --mem=160G --time=1:30:00 \
         --output="${LOGS}/v4fs_am_oct${o}_%j.out" run_jl.slurm)
    pj=$(sbatch --parsable --dependency=afterok:${aj} --job-name="v4fs-paint-${o}" \
         --export=ALL,OCT=${o},CAT=${AMC},OUTD=${D}/halomaps_v4fs,TAG=${T} \
         --output="${LOGS}/v4fs_paint_oct${o}_%j.out" --error="${LOGS}/v4fs_paint_oct${o}_%j.err" run_paint.slurm)
    cj=$(sbatch --parsable --dependency=afterok:${aj} --job-name="v4fs-cib-${o}" \
         --export=ALL,OCT=${o},CAT=${AMC},OUTD=${D}/cibmaps_v4fs,TAG=${T} \
         --output="${LOGS}/v4fs_cib_oct${o}_%j.out" --error="${LOGS}/v4fs_cib_oct${o}_%j.err" run_cib_prod.slurm)
    echo "oct${o}: am ${aj}  paint ${pj}  cib ${cj}"; deps="${deps}:${aj}:${pj}:${cj}"
done
sed 's/v3-analysis/v4fs-analysis/; s/CAMPAIGN=v3$/CAMPAIGN=v4fs/' "${REPO}"/validation/paper/v3_analysis.slurm > v4fs_analysis.slurm
an=$(sbatch --parsable --dependency=afterok${deps} --output="${LOGS}/v4fs_analysis_%j.out" v4fs_analysis.slurm)
cat > tails_v4fs.sh <<'EOT'
#!/bin/bash
module load julia/1.12.5
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"
cd /home/yguan/work/PeakPatch.jl
for o in 000 001 010 011 100 101 110 111; do TAG=v4fs OCT=$o julia --project=validation -t 16 validation/paper/check_am_tail.jl; done
EOT
tl=$(sbatch --parsable --dependency=afterok${deps} --account=aip-aspuru-ab --partition=gpubase_l40s_b1 --gpus-per-node=1 \
     --cpus-per-task=16 --mem=64G --time=1:30:00 --output="${LOGS}/v4fs_tails_%j.out" tails_v4fs.sh)
ta=$(SCRIPT=validation/paper/tierA_v3fs.jl TIERA_CAMP=v4 sbatch --parsable --export=ALL --dependency=afterok${deps} \
     --job-name=tierA-v4fs --cpus-per-task=16 --mem=96G --time=2:55:00 --output="${LOGS}/tierA_v4fs_%j.out" run_jl.slurm)
echo "analysis: ${an}  tails: ${tl}  tierA: ${ta}"
