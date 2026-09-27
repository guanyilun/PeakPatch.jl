#!/usr/bin/env julia
# Diagnostic: where does Julia lose peaks vs Fortran? Run the global (Fortran-
# equivalent) path on the A/B config with dump-counter instrumentation.
# Fortran reference (same inputs): 265858 peaks found -> 224181 kept, 41677 dumped.
using TOML, PeakPatch, Printf
import PeakPatch.RadialShell as RS

cfg = PipelineConfig(TOML.parsefile(joinpath(@__DIR__, "config_ab_test.toml")))
seed = 13579; ntile = 4

RS.reset_dump_counters!()
@info "Running GLOBAL path with dump instrumentation..."
halos = run_multitile(cfg; ntile=ntile, seed=seed, verbose=false)
d = RS.get_dump_counts()
kept = length(halos)
dumped = d.premask + d.no_fcrit + d.m0_one + d.no_collapse + d.rthl_neg
found = kept + dumped
@printf("\n==== JULIA global, %d^3 box=320.8 z=0 (vs Fortran) ====\n", cfg.n*0+256)
@printf("  peaks found (kept+dumped) = %8d   [Fortran 265858]\n", found)
@printf("  kept                      = %8d   [Fortran 224181]\n", kept)
@printf("  dumped total              = %8d   [Fortran  41677]\n", dumped)
@printf("  -- dump breakdown --\n")
@printf("     premask  (excluded)    = %8d\n", d.premask)
@printf("     no_collapse            = %8d\n", d.no_collapse)
@printf("     no_fcrit               = %8d\n", d.no_fcrit)
@printf("     m0_one                 = %8d\n", d.m0_one)
@printf("     rthl_neg               = %8d\n", d.rthl_neg)
