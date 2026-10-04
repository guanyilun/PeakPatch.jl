using PeakPatch
using Test
using HDF5

@testset "Provenance and compact catalogs" begin
    root = pkgdir(PeakPatch)
    tmp = mktempdir()
    cfgp = joinpath(tmp, "cfg.toml")
    write(cfgp, """
    [cosmology]
    Om = 0.31
    OB = 0.049
    OL = 0.69
    h = 0.68
    [grid]
    n = 64
    boxsize = 54.5
    nbuff = 8
    [run]
    seed = 777
    [files]
    pk = "validation/websky_6144/data/pk_websky.dat"
    filterbank = "validation/websky_6144/data/filters_websky.dat"
    homeltab = "validation/websky_6144/data/HomelTab_websky.dat"
    """)
    p = capture_provenance(cfgp; allow_dirty=true, extra=Dict("env.gpu" => "test", "note" => "x"))
    @test p["recipe.seed"] == "777"
    @test p["recipe.config_toml"] == read(cfgp, String)
    @test p["extra.note"] == "x" && p["env.gpu"] == "test"
    @test all(p["input.$k.sha256"] != "missing" for k in ("pk", "filterbank", "homeltab"))
    @test length(p["code.commit"]) == 40 || p["code.commit"] == ""      # empty outside a git checkout
    @test p["code.dirty"] in ("true", "false")
    # verification against this checkout: inputs always match; code matches the captured state
    @test isempty(filter(m -> startswith(m, "input"), verify_provenance_inputs(p)))
    # sidecar round trip, with the product hash
    prod = joinpath(tmp, "cat.pksc"); write(prod, rand(UInt8, 1000))
    side = write_provenance_sidecar(prod, p)
    q = read_provenance_sidecar(side)
    @test q["recipe.config_toml"] == p["recipe.config_toml"] && haskey(q, "lineage.product_sha256")
    # post-processing stage
    add_stage!(q, "am"; inputs=Dict("table" => prod), allow_dirty=true, notes=Dict("octant" => "000"))
    @test q["lineage.am.octant"] == "000" && q["lineage.am.input.table.sha256"] == q["lineage.product_sha256"]
    # compact catalog: quantisation bounds and provenance round trip
    n = 2000
    hs = [ExtHaloRecord(Float32(100 * randn()), Float32(100 * randn()), Float32(100 * randn()),
                        Float32(300 * randn()), Float32(300 * randn()), Float32(300 * randn()),
                        Float32(1 + 5rand()), 0f0, 0f0, 0f0, 1.7f0,
                        Float32(0.3rand()), Float32(0.2rand() - 0.1), ntuple(_ -> 0f0, 7)..., Float32(3rand()),
                        ntuple(_ -> 0f0, 12)...) for _ in 1:n]
    out = joinpath(tmp, "c.h5")
    write_compact_catalog(out, hs, q; extra_fields=true)
    c = read_compact_catalog(out)
    @test maximum(abs.(c.x .- [h.x for h in hs])) <= 0.005 + 1e-4
    @test maximum(abs.(c.vz .- [h.vz for h in hs])) <= 0.5 + 1e-3
    Mt = [4π / 3 * 2.775e11 * 0.31 * Float64(h.RTHL)^3 for h in hs]
    @test maximum(abs.(c.M ./ Mt .- 1)) < 1.2e-4
    @test maximum(abs.(c.zform .- [h.zform for h in hs])) <= 1.0001e-4
    @test c.provenance == q
end
