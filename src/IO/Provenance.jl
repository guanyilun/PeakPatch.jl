"""
    Provenance

Reproduction metadata that travels with every output product (catalogs, maps), so a file is self-contained:
the exact recipe (config text, seed), the exact code (git commit + dirty flag + diff), the exact inputs
(SHA-256 of P(k), filter bank, collapse table, …), and the environment that produced it.

A `Provenance` is a flat `Dict{String,String}` (portable: HDF5 attributes, FITS header cards, TOML).
Keys are namespaced: `recipe.*`, `code.*`, `input.*`, `env.*`, `run.*`, plus free-form `extra.*`.
"""
module Provenance

using SHA
using TOML
using Dates

export capture_provenance, provenance_toml, write_provenance_sidecar, read_provenance_sidecar,
       verify_provenance_inputs, add_stage!

const PP_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
# code paths whose state defines the build (data/validation edits do not change the catalog)
const CODE_PATHS = ["src", "ext", "Project.toml", "Manifest.toml"]

_sha256_file(path) = open(io -> bytes2hex(sha256(io)), path)
_sha256_str(s::AbstractString) = bytes2hex(sha256(s))

function _git(args::AbstractVector; repo=PP_ROOT)
    try
        return strip(read(pipeline(`git -C $repo $(String.(args))`; stderr=devnull), String))
    catch
        return ""
    end
end

"""
    capture_provenance(config_path; allow_dirty=ENV["PEAKPATCH_ALLOW_DIRTY"]=="1", extra=Dict())
        -> Dict{String,String}

Record everything needed to regenerate a product from `config_path` (a PipelineConfig TOML).

- `recipe.config_toml` holds the full config text (not a path). `recipe.seed` holds the seed.
- `code.commit` holds the git commit. `code.dirty` is true when `src/`, `ext/`, `Project.toml`, `Manifest.toml` or
  any of `scripts` (driver scripts, hashed as `code.script<i>.*`) differ from HEAD. In that case the run errors unless `allow_dirty`. When allowed, the full diff is embedded as
  `code.diff` together with its hash `code.diff_sha256`.
- `input.<name>.path` and `input.<name>.sha256` cover each file the config points to (pk, filterbank, homeltab).
- `env.*` holds the host, Julia and PeakPatch versions, the active Manifest hash, the thread count and the UTC
  date. Callers add GPU and driver details through `extra` (`env.gpu`, …).
"""
function capture_provenance(config_path::AbstractString;
                            allow_dirty::Bool=get(ENV, "PEAKPATCH_ALLOW_DIRTY", "0") == "1",
                            scripts::AbstractVector{<:AbstractString}=String[],
                            extra::AbstractDict=Dict{String,String}())
    p = Dict{String,String}()
    txt = read(config_path, String)
    cfg = TOML.parse(txt)
    p["recipe.config_path"] = abspath(config_path)
    p["recipe.config_toml"] = txt
    p["recipe.config_sha256"] = _sha256_str(txt)
    run = get(cfg, "run", Dict{String,Any}())
    haskey(run, "seed") && (p["recipe.seed"] = string(run["seed"]))
    # code state
    commit = _git(["rev-parse", "HEAD"])
    p["code.commit"] = commit
    p["code.describe"] = _git(["describe", "--always", "--dirty", "--tags"])
    # driver scripts (e.g. validation/websky_6144/run_gpu_octant.jl) also shape the product: hash them and
    # include them in the dirty check
    spaths = String[]
    for (i, sc) in enumerate(scripts)
        sp = abspath(sc)
        p["code.script$(i).path"] = relpath(sp, PP_ROOT)
        p["code.script$(i).sha256"] = isfile(sp) ? _sha256_file(sp) : "missing"
        push!(spaths, relpath(sp, PP_ROOT))
    end
    p["code.paths"] = join(vcat(CODE_PATHS, spaths), ",")
    diff = _git(["diff", "HEAD", "--", CODE_PATHS..., spaths...])
    dirty = !isempty(diff)
    p["code.dirty"] = string(dirty)
    if dirty
        allow_dirty || error("capture_provenance: the PeakPatch code tree ($(p["code.paths"])) differs from " *
                             "HEAD $(commit[1:min(end, 10)]); commit first, or set PEAKPATCH_ALLOW_DIRTY=1 to embed the diff")
        p["code.diff"] = diff
        p["code.diff_sha256"] = _sha256_str(diff)
    end
    # inputs referenced by the config (paths relative to the repo root, as the pipeline resolves them)
    files = get(cfg, "files", Dict{String,Any}())
    for key in ("pk", "filterbank", "homeltab")
        haskey(files, key) || continue
        f = files[key]; fp = isabspath(f) ? f : joinpath(PP_ROOT, f)
        p["input.$key.path"] = f
        p["input.$key.sha256"] = isfile(fp) ? _sha256_file(fp) : "missing"
    end
    # environment
    p["env.host"] = gethostname()
    p["env.julia_version"] = string(VERSION)
    p["env.peakpatch_version"] = string(pkgversion(parentmodule(@__MODULE__)))
    proj = Base.active_project()
    man = proj === nothing ? "" : joinpath(dirname(proj), "Manifest.toml")
    p["env.manifest_path"] = man
    p["env.manifest_sha256"] = isfile(man) ? _sha256_file(man) : "missing"
    p["env.julia_threads"] = string(Threads.nthreads())
    p["env.date_utc"] = string(Dates.now(Dates.UTC))
    haskey(ENV, "SLURM_JOB_ID") && (p["env.slurm_job_id"] = ENV["SLURM_JOB_ID"])
    for (k, v) in extra
        p[startswith(string(k), r"(recipe|code|input|env|run|extra|encoding|lineage)\.") ? string(k) : "extra.$(k)"] = string(v)
    end
    return p
end

"""TOML text of a provenance dict (sorted keys; multi-line values quoted by TOML)."""
provenance_toml(p::AbstractDict) = sprint(io -> TOML.print(io, Dict("provenance" => Dict(p)); sorted=true))

"""
    write_provenance_sidecar(product_path, p) -> sidecar path

For formats without metadata support (`.pksc`): `<product>.prov.toml` next to the product, including the
product's own SHA-256 so a separated sidecar can be matched back to its file.
"""
function write_provenance_sidecar(product_path::AbstractString, p::AbstractDict)
    q = Dict{String,String}(p)
    isfile(product_path) && (q["lineage.product_sha256"] = _sha256_file(product_path))
    q["lineage.product_file"] = basename(product_path)
    side = product_path * ".prov.toml"
    open(io -> write(io, provenance_toml(q)), side, "w")
    return side
end

read_provenance_sidecar(path::AbstractString) =
    Dict{String,String}(string(k) => string(v) for (k, v) in TOML.parsefile(endswith(path, ".prov.toml") ? path : path * ".prov.toml")["provenance"])

"""
    verify_provenance_inputs(p) -> Vector{String}

Re-hash the recorded input files and the code tree against `p`. Returns a list of mismatches (empty = this
checkout and these inputs reproduce the recorded recipe).
"""
function verify_provenance_inputs(p::AbstractDict)
    bad = String[]
    for (k, v) in p
        m = match(r"^input\.(.+)\.sha256$", k); m === nothing && continue
        f = p["input.$(m[1]).path"]; fp = isabspath(f) ? f : joinpath(PP_ROOT, f)
        (isfile(fp) && _sha256_file(fp) == v) || push!(bad, "input $(m[1]) ($f)")
    end
    _git(["rev-parse", "HEAD"]) == get(p, "code.commit", "") || push!(bad, "code.commit (checkout is $(_git(["rev-parse", "--short", "HEAD"])))")
    diff = _git(["diff", "HEAD", "--", split(get(p, "code.paths", join(CODE_PATHS, ",")), ",")...])
    if get(p, "code.dirty", "false") == "true"
        _sha256_str(diff) == get(p, "code.diff_sha256", "") || push!(bad, "code.diff (working tree differs from the recorded diff)")
    else
        isempty(diff) || push!(bad, "code tree is dirty")
    end
    return bad
end

"""
    add_stage!(p, stage; scripts=String[], inputs=Dict(), allow_dirty=…, notes=Dict()) -> p

Append a post-processing stage (e.g. `"am"` for abundance matching) to an existing provenance `p` (usually read
from the parent product's sidecar or file), under `lineage.<stage>.*`: the code state at that stage (commit,
dirty flag and diff over src/ext/Project/Manifest plus `scripts`), SHA-256 of each extra input
(`inputs = Dict("table" => path, …)`), host and date, and free-form `notes`. Errors on a dirty tree unless
`allow_dirty` (default: ENV PEAKPATCH_ALLOW_DIRTY == "1").
"""
function add_stage!(p::AbstractDict, stage::AbstractString;
                    scripts::AbstractVector{<:AbstractString}=String[], inputs::AbstractDict=Dict{String,String}(),
                    allow_dirty::Bool=get(ENV, "PEAKPATCH_ALLOW_DIRTY", "0") == "1",
                    notes::AbstractDict=Dict{String,String}())
    pre = "lineage.$stage."
    spaths = String[relpath(abspath(sc), PP_ROOT) for sc in scripts]
    for (i, sp) in enumerate(spaths)
        fp = joinpath(PP_ROOT, sp)
        p[pre * "script$(i).path"] = sp
        p[pre * "script$(i).sha256"] = isfile(fp) ? _sha256_file(fp) : "missing"
    end
    commit = _git(["rev-parse", "HEAD"])
    diff = _git(["diff", "HEAD", "--", CODE_PATHS..., spaths...])
    p[pre * "commit"] = commit
    p[pre * "dirty"] = string(!isempty(diff))
    if !isempty(diff)
        allow_dirty || error("add_stage!($stage): code tree differs from HEAD $(commit[1:min(end, 10)]); commit first, or set PEAKPATCH_ALLOW_DIRTY=1")
        p[pre * "diff"] = diff; p[pre * "diff_sha256"] = _sha256_str(diff)
    end
    for (k, f) in inputs
        p[pre * "input.$(k).path"] = string(f)
        p[pre * "input.$(k).sha256"] = isfile(f) ? _sha256_file(f) : "missing"
    end
    p[pre * "host"] = gethostname(); p[pre * "date_utc"] = string(Dates.now(Dates.UTC))
    haskey(ENV, "SLURM_JOB_ID") && (p[pre * "slurm_job_id"] = ENV["SLURM_JOB_ID"])
    for (k, v) in notes; p[pre * string(k)] = string(v); end
    return p
end

end # module Provenance
