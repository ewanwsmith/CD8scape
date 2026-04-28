"""
runtests.jl — Top-level test runner for CD8scape Julia source files.

Run from the repo root:
    julia --project=src test/runtests.jl

Or run a single test file directly:
    julia --project=src test/test_path_utils.jl
"""

using Test

println("=" ^ 70)
println("CD8scape Julia Test Suite")
println("=" ^ 70)

test_files = [
    "test_path_utils.jl",
    "test_generate_peptides.jl",
    "test_clean_peptides.jl",
    "test_parse_aa_variants.jl",
    "test_parse_trajectories.jl",
    "test_parse_vcf.jl",
    "test_percentile.jl",
    "test_read_ncbi_frames.jl",
    "test_read_samfire_frames.jl",
    "test_simulate_variants.jl",
    "test_variant_fates.jl",
    "test_process_scores.jl",
]

@testset "CD8scape" begin
    for f in test_files
        path = joinpath(@__DIR__, f)
        println("\n--- Running: $f ---")
        include(path)
    end
end
