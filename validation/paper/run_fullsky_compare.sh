#!/bin/bash
# Full-sky ours-vs-Websky auto spectra (A3 error model, octant jackknife).
#   CAMPAIGN=prod|v2 bash run_fullsky_compare.sh [product ...]   products: tsz ksz isw kappa cib
set -euo pipefail
C="${CAMPAIGN:-prod}"
F=/home/yguan/projects/aip-aspuru-ab/yguan/websky/fullsky_${C}
W=/home/yguan/projects/aip-aspuru-ab/yguan/websky_ref
S=validation/paper/compare_auto.jl
cd /home/yguan/work/PeakPatch.jl
J="julia --project=validation -t ${SLURM_CPUS_PER_TASK:-32}"
P=("$@"); [ ${#P[@]} -eq 0 ] && P=(tsz ksz isw kappa cib)
for p in "${P[@]}"; do case $p in
  tsz)   $J $S ${C}_tsz $F/tsz_y_${C}_fullsky_nside4096.fits $W/tsz_2048.fits 4096 1 1 2048 --jk ;;
  ksz)   $J $S ${C}_ksz $F/ksz_total_uK_${C}_fullsky_nside4096.fits $W/ksz.fits 8192 1 1 0 --jk ;;
  isw)   $J $S ${C}_isw $F/isw_uK_${C}_fullsky_nside4096.fits $W/isw.fits 1024 1 2.7255e6 0 --jk ;;
  kappa) $J $S ${C}_kappa_lt4.5 $F/kappa_lt4.5_${C}_fullsky_nside4096.fits $W/kap_lt4.5.fits 8192 1 1 0 --jk ;;
  cib)   for nu in 0100 0143 0217 0353 0545 0857; do
           $J $S ${C}_cib${nu} $F/cib_nu${nu}_wcut_${C}_fullsky_nside4096.fits $W/cib_nu${nu}.fits 8192 1 1 0 --jk
         done ;;
esac; done
