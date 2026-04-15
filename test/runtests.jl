using PeakPatch
using Test

@testset "PeakPatch.jl" begin
    include("test_cosmology.jl")
    include("test_powerspectrum.jl")
    include("test_filters.jl")
    include("test_randomfield.jl")
    include("test_lpt.jl")
    include("test_parameters.jl")
    include("test_catalog.jl")
    include("test_ellipsoidalcollapse.jl")
    include("test_peakfind.jl")
    include("test_radialshell.jl")
    include("test_radialshell_fortran.jl")
    include("test_pipeline.jl")
    include("test_merger.jl")
    include("test_multitile.jl")
    include("test_shell_gpu.jl")
end
