module PeakPatch

# Re-export submodules for convenient access
include("Cosmology/Cosmology.jl")
include("InitialConditions/PowerSpectrum.jl")
include("InitialConditions/LCG.jl")
include("InitialConditions/RandomField.jl")
include("InitialConditions/LPT.jl")
include("InitialConditions/NonGaussian.jl")
include("HaloFinder/Filters.jl")
include("HaloFinder/PeakFind.jl")
include("IO/Parameters.jl")
include("IO/Catalog.jl")
include("EllipsoidalCollapse/EllipsoidalCollapse.jl")
include("EllipsoidalCollapse/CollapseTable.jl")
include("HaloFinder/RadialShell.jl")
include("Merger/Exclusion.jl")
include("Merger/Merger.jl")
include("Pipeline.jl")
include("MultiTile.jl")

# Re-export public API
using .Cosmology: CosmologyParams, E2, H, chi, growth_factor, growth_rate, delta_c,
    DlinearTables, Dlinear_tables, Dlinear_ab, Dfnofa,
    ChiToZTable, build_chi_to_z, chi_to_z, peak_redshift
using .PowerSpectrum: load_pk, load_pk_nongaussian
using .LCG: LCG
using .NonGaussian: apply_fnl_correlated!, apply_fnl_uncorrelated!
using .RandomField: generate_grf, generate_grf_lcg,
    fill_noise_threefry!, fill_noise_threefry_region!
using .LPT: displacements_1lpt, displacements_2lpt
using .Filters: gaussian_window, gaussian_window_fortran, tophat_window, read_filterbank
using .PeakFind: PeakCandidate, find_peaks
using .RadialShell: ShellCell, PeakGrid, PeakResult, no_collapse,
    hRinteg, atab4, precompute_shells, analyse_peak, normalize_strain!, normalize_strain,
    fsc_of_z, get_evals, reset_dump_counters!, get_dump_counts
using .Parameters: PipelineConfig, FortranParams, read_params_bin, write_params_bin
using .Catalog: HaloRecord, ExtHaloRecord, write_pksc, read_pksc
using .EllipsoidalCollapse: EllipsoidParams, evolve_ellipse_full,
    get_b_2, _elliptic_rd
using .CollapseTable: CollapseTableParams, CollapseTableInterp,
    make_table, make_table_threaded,
    write_homeltab, read_homeltab, interpolate
using .Exclusion: SpatialHash, build_hash, sphere_overlap,
    lagrangian_exclusion!, volume_reduction!
using .Merger: merge_catalog
using .Pipeline: run_tile
using .MultiTile: run_multitile, extract_tile, tile_center

# MPI extension stub — method defined in ext/MPIExt.jl when MPI+PencilFFTs are loaded
function run_multitile_mpi end

# HDF5 extension stub — method defined in ext/HDF5Ext.jl when HDF5 is loaded
function write_catalog_hdf5 end

export
    CosmologyParams, E2, H, chi, growth_factor, growth_rate, delta_c,
    DlinearTables, Dlinear_tables, Dlinear_ab, Dfnofa,
    ChiToZTable, build_chi_to_z, chi_to_z, peak_redshift,
    load_pk, load_pk_nongaussian,
    apply_fnl_correlated!, apply_fnl_uncorrelated!,
    generate_grf,
    fill_noise_threefry!, fill_noise_threefry_region!,
    displacements_1lpt, displacements_2lpt,
    gaussian_window, gaussian_window_fortran, tophat_window, read_filterbank,
    PeakCandidate, find_peaks,
    ShellCell, PeakGrid, PeakResult,
    hRinteg, atab4, precompute_shells, analyse_peak, normalize_strain!, normalize_strain,
    fsc_of_z, get_evals, reset_dump_counters!, get_dump_counts,
    PipelineConfig, FortranParams, read_params_bin, write_params_bin,
    HaloRecord, ExtHaloRecord, write_pksc, read_pksc,
    EllipsoidParams, evolve_ellipse_full, get_b_2,
    CollapseTableParams, CollapseTableInterp,
    make_table, make_table_threaded,
    write_homeltab, read_homeltab, interpolate,
    sphere_overlap, merge_catalog,
    run_tile,
    run_multitile, extract_tile, tile_center,
    run_multitile_mpi,
    write_catalog_hdf5

end # module PeakPatch
