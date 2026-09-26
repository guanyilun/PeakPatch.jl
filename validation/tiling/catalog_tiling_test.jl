# Catalog-level tiling-independence test, plan item A5(a) (docs/paper_comparison_plan_2026-09.md).
#
# Same seed, same global grid (N, Websky cellsize 0.85221 Mpc/h), same coarse grid M
# (block 12 = production cf32), same nbuff=16, finecell filter bank, 2LPT, z=0 snapshot;
# run_multitile_split with ntile = 1, 2, 4, 8 (tile sizes nmesh = NC/ntile + 2 nbuff).
# The white noise and the coarse field are bitwise identical across tilings by construction,
# so differences can only come from the per-tile isolated residual convolution and the
# per-tile peak search / exclusion at tile faces. Matches merged catalogs halo-by-halo.
#
#   julia --project=validation -t 16 validation/tiling/catalog_tiling_test.jl [outdir]
const USE_GPU = get(ENV, "CAT_GPU", "1") == "1"
USE_GPU && @eval using CUDA
using PeakPatch
using Printf, Statistics

const OUT = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "results")
mkpath(OUT)
const DATA = joinpath(@__DIR__, "..", "websky_6144", "data")
const N = parse(Int, get(ENV, "CAT_N", "480"))
const NBUFF = 16
const BLOCK = 12
const M = N ÷ BLOCK
const ALATT = 0.85221
const NC = N - 2NBUFF
const SEED = 12345
const NTILES = parse.(Int, split(get(ENV, "CAT_NTILES", "1,2,4,8"), ","))
const RHO_M = 2.775e11 * 0.31
mass(h) = 4π / 3 * Float64(h.RTHL)^3 * RHO_M

println("GPU: ", USE_GPU ? Base.invokelatest(() -> CUDA.name(CUDA.device())) : "none (CPU)", "   N=$N  NC=$NC  M=$M  block=$BLOCK  nbuff=$NBUFF  box=$(N*ALATT) Mpc/h")

function make_cfg(nmesh)
    PipelineConfig(Dict{String,Any}(
        "cosmology" => Dict{String,Any}("Om" => 0.31, "OB" => 0.049, "OL" => 0.69, "h" => 0.68),
        "grid" => Dict{String,Any}("n" => nmesh, "boxsize" => nmesh * ALATT, "nbuff" => NBUFF,
                                   "cenx" => 0.0, "ceny" => 0.0, "cenz" => 0.0),
        "run" => Dict{String,Any}("ievol" => 0, "z_max" => 0.0, "z_out" => 0.0, "ilpt" => 2,
                                  "ioutshear" => 0, "wsmooth" => 1, "rmax2rs" => 0.0),
        "files" => Dict{String,Any}("pk" => joinpath(DATA, "pk_websky.dat"),
                                    "filterbank" => joinpath(DATA, "filters_websky_finecell.dat"),
                                    "homeltab" => joinpath(DATA, "HomelTab_websky.dat"),
                                    "output" => joinpath(OUT, "unused.pksc"))))
end

using Random
# Canonical order: descending RTHL, ties broken by (x, y, z). merge_catalog's stable
# sortperm then visits tied radii in a fixed order, independent of tile completion order.
canon(hs) = hs[sortperm(hs; by=h -> (-h.RTHL, h.x, h.y, h.z))]
recset(hs) = Set((h.x, h.y, h.z, h.vx, h.vy, h.vz, h.RTHL, h.vx2, h.vy2, h.vz2) for h in hs)
function save(path, hs)
    open(path, "w") do io      # x y z vx vy vz RTHL vx2 vy2 vz2 per halo, Float32
        for h in hs; write(io, h.x, h.y, h.z, h.vx, h.vy, h.vz, h.RTHL, h.vx2, h.vy2, h.vz2); end
    end
end
runcat(cfg, nt) = Base.invokelatest(run_multitile_split, cfg; ntile=nt, seed=SEED, coarse_grid=M, use_gpu=USE_GPU)

cats = Dict{Int,Any}()
const REPEAT_NT = parse(Int, get(ENV, "CAT_REPEAT", "4"))
for nt in NTILES
    nsub = NC ÷ nt; nmesh = nsub + 2NBUFF
    cfg = make_cfg(nmesh)
    t0 = time()
    try
        raw = runcat(cfg, nt)
        t1 = time()
        merged = merge_catalog(canon(raw))
        merged_asis = merge_catalog(raw)
        cats[nt] = (; raw, merged, merged_asis, nsub, nmesh)
        save(joinpath(OUT, "raw_ntile$(nt).f32"), canon(raw)); save(joinpath(OUT, "merged_ntile$(nt).f32"), merged)
        nties = length(raw) - length(unique(h.RTHL for h in raw))
        @printf("ntile=%d nmesh=%d: %d raw → %d merged (canonical order), %d merged (as returned)  [RTHL ties in raw list: %d]  (field+find %.0f s)\n",
                nt, nmesh, length(raw), length(merged), length(merged_asis), nties, t1 - t0)
        # (ii) merge-order sensitivity on a fixed raw list: 3 random input permutations
        rs = recset(merged)
        for trial in 1:3
            mp = merge_catalog(shuffle(MersenneTwister(trial), raw))
            @printf("    merge of shuffled input #%d: %d halos, %d records differ from canonical merge\n",
                    trial, length(mp), length(symdiff(recset(mp), rs)))
        end
        @printf("    merge of as-returned order: %d records differ from canonical merge\n", length(symdiff(recset(merged_asis), rs)))
        # (iii) run-to-run determinism of the raw peak list (same tiling, same GPU)
        if nt == REPEAT_NT
            raw2 = runcat(cfg, nt)
            @printf("    repeat run ntile=%d: %d raw; raw peak lists bit-identical as sets: %s (%d records differ); same order: %s\n",
                    nt, length(raw2), recset(raw2) == recset(raw), length(symdiff(recset(raw2), recset(raw))),
                    all(raw2[i] === raw[i] for i in 1:min(length(raw), length(raw2))) && length(raw) == length(raw2))
            m2a = merge_catalog(raw2); m2c = merge_catalog(canon(raw2))
            @printf("    repeat run merged: as-returned order %d vs %d halos, %d records differ run-to-run; canonical order %d records differ run-to-run\n",
                    length(m2a), length(merged_asis), length(symdiff(recset(m2a), recset(merged_asis))), length(symdiff(recset(m2c), rs)))
        end
    catch e
        @printf("ntile=%d nmesh=%d FAILED: %s\n", nt, nmesh, sprint(showerror, e)[1:min(end, 300)])
    end
    flush(stdout); GC.gc(); USE_GPU && Base.invokelatest(CUDA.reclaim)
end

# ---- halo matching (Lagrangian peak position, nearest unclaimed within tolerance) ----
# tolerance max(1 cell, 0.2 R_TH): flat, heavily smoothed peaks of big halos can shift by a
# cell or two under a sub-percent field perturbation. Capped at 3 cells (search radius).
tol(h) = min(3.0, max(1.0, 0.2 * h.RTHL / ALATT))
key(h) = (floor(Int, h.x / ALATT), floor(Int, h.y / ALATT), floor(Int, h.z / ALATT))
function build_hash(hs)
    d = Dict{NTuple{3,Int},Vector{Int}}()
    for (i, h) in enumerate(hs); push!(get!(d, key(h), Int[]), i); end
    d
end
# distance (cells) to the nearest tile face of a tiling with core size nsub
function face_dist(h, nsub)
    minimum(begin
        u = c / ALATT + NC / 2 + 0.5; w = mod(u - 0.5, nsub); min(w, nsub - w)
    end for c in (h.x, h.y, h.z))
end
function match(ref, cand)
    hsh = build_hash(cand); used = falses(length(cand))
    mi = zeros(Int, length(ref)); md = fill(Inf, length(ref))
    # greedy by descending reference mass so big halos claim their partner first
    for i in sortperm(ref; by=h -> -h.RTHL)
        h = ref[i]; k = key(h); best = 0; bd = Inf
        for dz in -3:3, dy in -3:3, dx in -3:3
            for j in get(hsh, (k[1]+dx, k[2]+dy, k[3]+dz), Int[])
                used[j] && continue
                g = cand[j]
                d = sqrt((h.x-g.x)^2 + (h.y-g.y)^2 + (h.z-g.z)^2) / ALATT
                d < bd && (bd = d; best = j)
            end
        end
        if best > 0 && bd <= tol(h)
            mi[i] = best; md[i] = bd; used[best] = true
        end
    end
    mi, md
end

# ---- (i) pre-merge peak lists: exact comparison ----
# Peak positions are grid-cell centres, so a peak present in both lists has identical x,y,z.
println("\n(i) raw (pre-merge) peak lists vs ntile=$(minimum(keys(cats)))")
open(joinpath(OUT, "raw_compare.csv"), "w") do io
    println(io, "ntile,d_lo,d_hi,n_ref,pos_found_frac,bit_identical_frac,rthl_identical_frac,median_abs_dlnR")
    rnt = minimum(keys(cats)); rr = cats[rnt].raw
    for nt in sort(collect(keys(cats)))
        nt == rnt && continue
        cr = cats[nt].raw; nsub = cats[nt].nsub
        byp = Dict((h.x, h.y, h.z) => h for h in cr)
        @printf("  ntile=%d: %d vs %d raw peaks\n", nt, length(cr), length(rr))
        for (lo, hi) in ((0, 2), (2, 4), (4, 8), (8, 16), (16, 32), (32, 1000), (0, 1000))
            ir = [h for h in rr if lo <= face_dist(h, nsub) < hi]; isempty(ir) && continue
            f = [get(byp, (h.x, h.y, h.z), nothing) for h in ir]
            nf = count(!isnothing, f)
            nb = count(q -> f[q] !== nothing && f[q] === ir[q], eachindex(ir))
            nrr = count(q -> f[q] !== nothing && f[q].RTHL == ir[q].RTHL, eachindex(ir))
            dl = [abs(log(f[q].RTHL / ir[q].RTHL)) for q in eachindex(ir) if f[q] !== nothing]
            @printf("    face dist [%3d,%4d): n=%7d  position found %.4f  bit-identical record %.4f  identical RTHL %.4f  median|ΔlnR| %.2e\n",
                    lo, hi, length(ir), nf / length(ir), nb / length(ir), nrr / length(ir), isempty(dl) ? NaN : median(dl))
            @printf(io, "%d,%d,%d,%d,%.6f,%.6f,%.6f,%.4e\n", nt, lo, hi, length(ir), nf / length(ir), nb / length(ir), nrr / length(ir), isempty(dl) ? NaN : median(dl))
        end
    end
end

const MEDGES = [1e12, 3e12, 1e13, 3e13, 1e14, 1e16]
isempty(cats) && error("no successful runs")
refnt = minimum(keys(cats))
ref = cats[refnt].merged
println("\nreference: ntile=$refnt ($(length(ref)) merged halos)")
open(joinpath(OUT, "catalog_match_summary.csv"), "w") do io
    println(io, "ntile,ref_ntile,mass_lo,mass_hi,n_ref,n_cand,matched_frac,median_lnM_ratio,rms_lnM_ratio,rms_dpsi_over_sigma")
    for nt in sort(collect(keys(cats)))
        nt == refnt && continue
        cand = cats[nt].merged; nsub = cats[nt].nsub
        mi, md = match(ref, cand)
        sig_psi = std([h.vx for h in ref])
        @printf("\nntile=%d vs %d: %d vs %d halos (N ratio %.5f)\n", nt, refnt, length(cand), length(ref), length(cand) / length(ref))
        for b in 1:length(MEDGES)-1
            lo, hi = MEDGES[b], MEDGES[b+1]
            ir = [i for i in eachindex(ref) if lo <= mass(ref[i]) < hi]
            nc = count(h -> lo <= mass(h) < hi, cand)
            im = [i for i in ir if mi[i] > 0]
            lr = [log(mass(cand[mi[i]]) / mass(ref[i])) for i in im]
            dp = [(cand[mi[i]].vx - ref[i].vx) for i in im]
            fm = isempty(ir) ? NaN : length(im) / length(ir)
            medr = isempty(lr) ? NaN : median(lr); rmsr = isempty(lr) ? NaN : sqrt(mean(abs2, lr))
            rdp = isempty(dp) ? NaN : sqrt(mean(abs2, dp)) / sig_psi
            @printf("  M=[%.0e,%.0e): n_ref=%7d n_cand=%7d matched %.4f  ln(M ratio) median %+.2e rms %.2e  ψx diff rms/σ %.2e\n",
                    lo, hi, length(ir), nc, fm, medr, rmsr, rdp)
            @printf(io, "%d,%d,%.3e,%.3e,%d,%d,%.6f,%.4e,%.4e,%.4e\n", nt, refnt, lo, hi, length(ir), nc, fm, medr, rmsr, rdp)
        end
        # unmatched fraction vs distance to the nearest tile face (of the candidate tiling)
        dbins = [0, 2, 4, 8, 16, 32, 1000]
        println("  unmatched fraction of ref halos vs distance to nearest ntile=$nt tile face (cells):")
        open(joinpath(OUT, "catalog_unmatched_vs_face_ntile$(nt).csv"), "w") do io2
            println(io2, "d_lo,d_hi,n_ref,unmatched_frac,exact_mass_frac")
            for b in 1:length(dbins)-1
                ir = [i for i in eachindex(ref) if dbins[b] <= face_dist(ref[i], nsub) < dbins[b+1]]
                isempty(ir) && continue
                fu = count(i -> mi[i] == 0, ir) / length(ir)
                fe = count(i -> mi[i] > 0 && cand[mi[i]].RTHL == ref[i].RTHL, ir) / length(ir)
                @printf("    [%3d,%4d): n=%7d unmatched %.4f  identical RTHL %.4f\n", dbins[b], dbins[b+1], length(ir), fu, fe)
                @printf(io2, "%d,%d,%d,%.6f,%.6f\n", dbins[b], dbins[b+1], length(ir), fu, fe)
            end
        end
        pos = md[mi .> 0]
        @printf("  matched position offsets (cells): median %.3f, 99%% %.3f, exact(0) %.4f\n",
                median(pos), quantile(pos, 0.99), count(==(0.0), pos) / length(pos))
    end
end
# cumulative counts N(>M) per tiling
open(joinpath(OUT, "catalog_cumcounts.csv"), "w") do io
    println(io, "ntile,M,N_gt_M")
    for nt in sort(collect(keys(cats))), Mth in 10 .^ (12:0.25:15)
        @printf(io, "%d,%.4e,%d\n", nt, Mth, count(h -> mass(h) > Mth, cats[nt].merged))
    end
end
println("\nwrote results to $OUT")
