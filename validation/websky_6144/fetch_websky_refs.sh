#!/bin/bash
# Fetch the released Websky v0.0 products used by the paper comparisons
# (docs/paper_comparison_plan_2026-09.md, item A0) into persistent /project storage.
# Resumable (wget -c); writes md5sums.txt. The full halos.pksc is large — pass --halos to include it.
set -euo pipefail
DEST=/home/yguan/projects/aip-aspuru-ab/yguan/websky_ref
URL=https://mocks.cita.utoronto.ca/data/websky/v0.0
mkdir -p "$DEST"; cd "$DEST"
files=(kap.fits kap_lt4.5.fits kap_gt4.5.fits ksz.fits tsz.fits isw.fits
       cib_nu0100.fits cib_nu0143.fits cib_nu0217.fits cib_nu0353.fits cib_nu0545.fits cib_nu0857.fits)
[[ "${1:-}" == "--halos" ]] && files+=(halos.pksc)
for f in "${files[@]}"; do
    [[ -s "$f" ]] && { echo "have $f"; continue; }
    wget -c -q --show-progress "$URL/$f" || echo "FAILED $f"
done
md5sum *.fits *.pksc 2>/dev/null > md5sums.txt
ls -la
