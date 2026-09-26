#!/bin/bash
#SBATCH --job-name=pp-bench
#SBATCH --account=aip-aspuru-ab
#SBATCH --nodes=1
#SBATCH --time=1:30:00
# Partition, --gpus-per-node, --cpus-per-task, --mem and log paths are set by submit_bench.sh.
#
# Strong-scaling point: CONFIG on NDEV GPUs (warm-up run first, then timed run).
# GPU memory: nvidia-smi polled every 500 ms in the background → gpumem_<jobid>.csv.
set -euo pipefail
module load julia/1.12.5 cuda/12.6
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"
unset CUDA_VISIBLE_DEVICES
: "${CONFIG:?}"; : "${NDEV:?}"
B=/home/yguan/scratch/websky_6144/bench
MEMLOG=${B}/gpumem_${SLURM_JOB_ID}.csv
cd /home/yguan/work/PeakPatch.jl
echo "Node: $(hostname)  CONFIG=${CONFIG} NDEV=${NDEV}  Start: $(date)"
nvidia-smi --query-gpu=timestamp,index,name,memory.used,memory.total,utilization.gpu \
    --format=csv,noheader,nounits -lms 500 > "${MEMLOG}" &
SMI=$!
T0=$(date +%s)
julia --project=validation --threads=${JULIA_NUM_THREADS} validation/performance/bench_octant.jl \
    "${CONFIG}" "${NDEV}" "${B}/out_${SLURM_JOB_ID}" validation/performance/configs/warmup.toml
T1=$(date +%s)
kill ${SMI} || true
echo "JOB_WALL_S=$((T1 - T0))"
awk -F', ' '{ if ($4 > m[$2]) m[$2] = $4; s[$2] += $6; c[$2]++ }
     END { for (d in m) printf("GPU_PEAK dev=%s mem_used_MiB=%d mean_util_pct=%.1f\n", d, m[d], s[d]/c[d]) }' "${MEMLOG}"
rm -rf "${B}/out_${SLURM_JOB_ID}"
echo "End: $(date)"
