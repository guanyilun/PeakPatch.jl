#!/bin/bash
#SBATCH --job-name=env-check
#SBATCH --account=aip-aspuru-ab
#SBATCH --partition=gpubase_l40s_b1
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/home/yguan/scratch/websky_6144/logs/env_check_%j.out
#SBATCH --error=/home/yguan/scratch/websky_6144/logs/env_check_%j.err
# Checks the environment split: (1) validation/ script env loads CUDA + CUDAExt on a GPU
# and passes the GPU multires tests; (2) the root package env (CUDA weak) passes Pkg.test().
set -uo pipefail
module load julia/1.12.5 cuda/12.6
export JULIA_DEPOT_PATH="/home/yguan/scratch/julia_depot_killarney:${HOME}/.julia"
cd /home/yguan/work/PeakPatch.jl
echo "=== [1] validation env: CUDA + extension ==="
julia --project=validation -e 'using PeakPatch, CUDA; println("functional=", CUDA.functional(), " ext=", Base.get_extension(PeakPatch, :CUDAExt) !== nothing, " dev=", CUDA.name(CUDA.device()))'
echo "=== [1b] GPU multires tests (validation env) ==="
julia --project=validation -t 16 test/test_multiresolution_gpu.jl; echo "exit gpu_tests=$?"
echo "=== [2] root package Pkg.test() ==="
julia --project=. -t 16 -e 'using Pkg; Pkg.test()'; echo "exit pkg_test=$?"
