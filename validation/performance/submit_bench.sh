#!/bin/bash
# Submit the performance benchmark set (PERFORMANCE_2026-09-26.md). sbatch cannot run from
# /home, so the scripts are copied to scratch first. Run from anywhere:
#   bash validation/performance/submit_bench.sh
set -euo pipefail
R=/home/yguan/work/PeakPatch.jl/validation/performance
B=/home/yguan/scratch/websky_6144/bench
L=/home/yguan/scratch/websky_6144/logs
mkdir -p "$B" "$L"
cp "$R/bench.slurm" "$R/mem_sweep.slurm" "$B/"
cd "$B"
CFG=validation/performance/configs/scale.toml
# strong scaling on L40S: 8 CPU cores per GPU (production: 32 cores / 4 GPUs)
for nd in 1 2 4; do
    sbatch --partition=gpubase_l40s_b1 --gpus-per-node=$nd --cpus-per-task=$((8 * nd)) --mem=$((60 * nd))G \
        --output="$L/bench_scale_l40s_g${nd}_%j.out" --error="$L/bench_scale_l40s_g${nd}_%j.err" \
        --export=ALL,CONFIG=$CFG,NDEV=$nd bench.slurm
done
# one H100 point (same config, 1 GPU)
sbatch --partition=gpubase_h100_b1 --gpus-per-node=1 --cpus-per-task=6 --mem=120G \
    --output="$L/bench_scale_h100_g1_%j.out" --error="$L/bench_scale_h100_g1_%j.err" \
    --export=ALL,CONFIG=$CFG,NDEV=1 bench.slurm
# GPU memory vs tile size
sbatch mem_sweep.slurm
