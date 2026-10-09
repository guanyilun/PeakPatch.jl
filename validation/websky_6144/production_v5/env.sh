# Shared environment for the production-v5 jobs on Flatiron Rusty; sourced by every job script
# and by submit_all.sh. Override any variable by exporting it before submitting.
#   WS_ROOT  maps: fieldmaps_v5fs/ halomaps_v5fs/ cibmaps_v5fs/ fullsky_v5fs/, and logs/
#   WS_CATS  catalogs: catalogs_v5/ (raw, _AMfs, AM table, .prov.toml sidecars)
#   WS_REF   released Websky files (halos_10x10.pksc, kap/tsz/ksz/isw/cib maps, halo_mass_completion.txt)
# Julia 1.12 lives in the newer module tree (the default tree only has 1.11.2, which cannot read
# the 1.12 Manifest). CUDA.jl brings its own runtime as an artifact, so no cuda module is needed.
module load modules/2.5-20261005 julia/1.12.7
export REPO="${REPO:-/mnt/home/yguan/work/uoft/peakpatch/PeakPatch.jl}"
export XGPAINT="${XGPAINT:-/mnt/home/yguan/work/uoft/peakpatch/XGPaint.jl}"
export WS_ROOT="${WS_ROOT:-/mnt/ceph/users/yguan/projects/uoft/peakpatch/websky_6144}"
export WS_CATS="${WS_CATS:-${WS_ROOT}}"
export WS_REF="${WS_REF:-${WS_ROOT}/websky_ref}"
export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
