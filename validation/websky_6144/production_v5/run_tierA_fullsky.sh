#!/bin/bash
# Submit the full-sky Tier-A (validation/paper/tierA_fullsky.jl) for v5: the same 305 caps in our 8 AM
# octants and in the FULL Websky catalog (websky_ref/halos.pksc, 862.9M halos). This is the test where v4
# was 5% too clustered (xi 3-15 Mpc/h W/ours 0.948 +- 0.002; TIERA_FULLSKY_2026-10-02.md).
# CPU only. v4 took 56 min; memory was not recorded (the patch Tier-A peaked at 96 GiB), so 192G.
# Works from any directory:  bash validation/websky_6144/production_v5/run_tierA_fullsky.sh
set -euo pipefail
export REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
source "${REPO}/validation/websky_6144/production_v5/env.sh"
SCRIPT=validation/paper/tierA_fullsky.jl TIERA_CAMP=v5 sbatch --export=ALL --chdir="${REPO}" \
    --job-name=v5fs-tierA-fullsky --cpus-per-task=32 --mem=192G --time=2:30:00 \
    --output="${WS_ROOT}/logs/v5fs_tierA_fullsky_%j.out" "${REPO}/validation/websky_6144/production_v5/run_jl.slurm"
