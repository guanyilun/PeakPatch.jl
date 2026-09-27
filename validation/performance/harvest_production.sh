#!/bin/bash
# Reproduce the frozen-campaign timing tables in PERFORMANCE_2026-09-26.md from the SLURM
# logs, sacct and file mtimes. The catalog logs carry no per-stage timestamps, so the stage
# split uses: job start (log) → pipeline end (start + elapsed_min; JIT/startup lumped into
# the pipeline) → raw catalog mtime (merge + finalize + write) → AM catalog mtime (AM).
set -uo pipefail
L=/home/yguan/scratch/websky_6144/logs
D=/home/yguan/projects/aip-aspuru-ab/yguan/websky
printf "%-4s %-7s %8s %8s %9s %9s %9s %8s %10s\n" oct node job_h pipe_min postpipe_h am_min shell_pct gpu_s maxrss_GB
for f in "$L"/prod_cat_oct*_44381*.out; do
    o=$(basename "$f" | sed -E 's/prod_cat_oct([01]{3})_.*/\1/'); j=$(basename "$f" .out | sed 's/.*_//')
    e=${f%.out}.err
    node=$(grep -o 'Node: [a-z0-9]*' "$f" | cut -d' ' -f2)
    t0=$(date -d "$(grep -o 'Start: .*' "$f" | cut -d' ' -f2-)" +%s)
    pipe=$(grep -o 'elapsed_min = [0-9.]*' "$e" | head -1 | awk '{print $3}')
    shell=$(grep '09_shell_analysis' "$e" | head -1 | grep -o '( *[0-9.]*%)' | tr -d '( %)')
    gpus=$(grep 'TOTAL_stages' "$e" | head -1 | grep -o 'total-GPU-s= *[0-9.]*' | grep -o '[0-9.]*$')
    traw=$(stat -c %Y "$D/catalog_websky_6144_prod_oct${o}.pksc")
    tam=$(stat -c %Y "$D/catalog_websky_6144_prod_oct${o}_AM.pksc")
    el=$(sacct -j "$j" -X -n --format=Elapsed | head -1 | tr -d ' ')
    rss=$(sacct -j "$j.batch" -n --format=MaxRSS | head -1 | tr -d ' K')
    post=$(awk -v a="$traw" -v b="$t0" -v p="$pipe" 'BEGIN{printf "%.2f", (a-b-60*p)/3600}')
    am=$(awk -v a="$tam" -v b="$traw" 'BEGIN{d=(a-b)/60; if (d > 600) print "rerun"; else printf "%.1f", d}')
    jh=$(echo "$el" | awk -F: '{printf "%.2f", $1+$2/60+$3/3600}')
    printf "%-4s %-7s %8s %8s %9s %9s %9s %8s %10.1f\n" "$o" "$node" "$jh" "$pipe" "$post" "$am" "$shell" "$gpus" \
        "$(awk -v r="$rss" 'BEGIN{print r/1048576}')"
done
echo
printf "%-4s %-7s %8s %10s %10s\n" oct node job_h paint_min maxrss_GB
for f in "$L"/prod_fld_oct*_44381*.out; do
    o=$(basename "$f" | sed -E 's/prod_fld_oct([01]{3})_.*/\1/'); j=$(basename "$f" .out | sed 's/.*_//')
    node=$(grep -o 'Node: [a-z0-9]*' "$f" | cut -d' ' -f2)
    paint=$(grep -o 'elapsed_min = [0-9.]*' "${f%.out}.err" | head -1 | awk '{print $3}')
    el=$(sacct -j "$j" -X -n --format=Elapsed | head -1 | tr -d ' ')
    rss=$(sacct -j "$j.batch" -n --format=MaxRSS | head -1 | tr -d ' K')
    printf "%-4s %-7s %8s %10s %10.1f\n" "$o" "$node" "$(echo "$el" | awk -F: '{printf "%.2f", $1+$2/60+$3/3600}')" \
        "$paint" "$(awk -v r="$rss" 'BEGIN{print r/1048576}')"
done
