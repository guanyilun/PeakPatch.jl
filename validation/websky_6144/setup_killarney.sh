#!/bin/bash
# One-time setup for Websky 6144^3 GPU runs on Killarney.
#
# Run on a login node (or compute node with Julia available):
#   module load julia/1.12.5 cuda/12.6
#   bash validation/websky_6144/setup_killarney.sh
#
# What this does:
#   1. Adds CUDA.jl to the project (needed for GPU pipeline)
#   2. Instantiates all dependencies
#   3. Generates data files (filters, collapse table, CAMB P(k))
#   4. Generates all 8 octant TOML configs
#   5. Creates output directories

set -euo pipefail

echo "=== Step 1: Instantiate the validation/ script env (PeakPatch + CUDA + Healpix) ==="
# CUDA stays a weak dep of the package; never Pkg.add it into the root env.
julia --project=validation -e 'using Pkg; Pkg.instantiate()'

echo ""
echo "=== Step 2: Generate data files (filters + collapse table) ==="
julia --project=validation validation/websky_6144/generate_data_files.jl

echo ""
echo "=== Step 3: Generate CAMB P(k) ==="
if [ ! -d "/home/yguan/scratch/camb_env" ]; then
    echo "Creating CAMB venv..."
    python3 -m venv /home/yguan/scratch/camb_env
    source /home/yguan/scratch/camb_env/bin/activate
    pip install camb
else
    source /home/yguan/scratch/camb_env/bin/activate
fi
python3 validation/websky_6144/generate_pk_camb.py
deactivate

echo ""
echo "=== Step 4: Generate octant configs ==="
python3 validation/websky_6144/generate_octant_configs.py

echo ""
echo "=== Step 5: Create output directories ==="
mkdir -p /home/yguan/scratch/websky_6144/logs
mkdir -p /home/yguan/projects/aip-aspuru-ab/yguan/websky

echo ""
echo "=== Setup complete ==="
echo "To submit all 8 octants:"
echo "  cp validation/websky_6144/run_websky_6144_killarney.slurm /home/yguan/scratch/websky_6144/slurm/"
echo "  cd /home/yguan/scratch/websky_6144/slurm"
echo "  for oct in 000 001 010 011 100 101 110 111; do"
echo "      sbatch --export=OCTANT=\$oct run_websky_6144_killarney.slurm"
echo "  done"
