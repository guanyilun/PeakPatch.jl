#!/bin/bash
#SBATCH --job-name=pp-memsweep
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b1
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=2:30:00
#SBATCH --output=/home/yguan/scratch/websky_6144/logs/bench_memsweep_%j.out
#SBATCH --error=/home/yguan/scratch/websky_6144/logs/bench_memsweep_%j.err
#
# GPU memory high-water mark vs tile size n (ntile=2, 8 tiles, production physics), one GPU,
# one fresh Julia process per n so the CUDA.jl pool starts empty. An OOM is a data point.
set -uo pipefail
module load julia/1.12.5 cuda/12.6
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"
unset CUDA_VISIBLE_DEVICES
B=/home/yguan/scratch/websky_6144/bench
cd /home/yguan/work/PeakPatch.jl
echo "Node: $(hostname)  Start: $(date)"
for n in 256 320 384 448 512 576; do
    CFG=validation/performance/configs/mem_n${n}.toml
    MEMLOG=${B}/gpumem_${SLURM_JOB_ID}_n${n}.csv
    nvidia-smi --query-gpu=timestamp,index,name,memory.used,memory.total,utilization.gpu \
        --format=csv,noheader,nounits -lms 250 > "${MEMLOG}" &
    SMI=$!
    T0=$(date +%s)
    julia --project=validation --threads=${JULIA_NUM_THREADS} validation/performance/bench_octant.jl \
        "${CFG}" 1 "${B}/out_${SLURM_JOB_ID}"
    RC=$?
    T1=$(date +%s)
    kill ${SMI} || true
    awk -F', ' -v n=${n} -v rc=${RC} -v w=$((T1 - T0)) '{ if ($4 > m) m = $4; tot = $5 }
         END { printf("MEMSWEEP n=%d rc=%d wall_s=%d peak_mem_used_MiB=%d total_MiB=%d\n", n, rc, w, m, tot) }' "${MEMLOG}"
done
rm -rf "${B}/out_${SLURM_JOB_ID}"
echo "End: $(date)"
