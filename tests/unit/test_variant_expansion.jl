using Test
using Random
include("../../src/variant_expansion.jl")
using .VariantExpansion

# deterministic RNG for sampling tests
rng = MersenneTwister(1234)

@testset "expand_sites - small enumeration" begin
    sites = [Site(100, "A", ["C"], nothing), Site(105, "G", ["T"], nothing)]
    combos = expand_sites(sites; max_alts_per_site=2, max_combinations=100)
    @test length(combos) == 4
    # combinations should be all pairs of [REF,ALT]
    expected = [["A","G"],["C","G"],["A","T"],["C","T"]]
    @test sort(combos) == sort(expected)
end

@testset "expand_sites - saturated sampling" begin
    # 5 sites with 3 options each => 3^5 = 243 combos; set max_combinations=10 to force sampling
    sites = [Site(i, "R", ["A","B"], [0.2,0.3]) for i in 1:5]
    combos = expand_sites(sites; max_alts_per_site=2, max_combinations=10, rng=rng, sample_on_saturate=true)
    @test length(combos) <= 10
    # ensure each combo length equals number of sites
    @test all(length(c) == length(sites) for c in combos)
end

@testset "expand_sites - saturated deterministic fallback" begin
    sites = [Site(i, "R", ["A","B"], nothing) for i in 1:6]
    combos = expand_sites(sites; max_alts_per_site=2, max_combinations=2, sample_on_saturate=false)
    @test length(combos) == 2
    @test combos[1] == [s.ref for s in sites]
end

@testset "expand_sites - saturation logging" begin
    tmplog_path, tmplog_io = Base.mktemp()
    close(tmplog_io)
    sites = [Site(i, "R", ["A","B","C"], nothing) for i in 1:8]
    # set max_combinations low to force saturation
    combos = expand_sites(sites; max_alts_per_site=3, max_combinations=5, sample_on_saturate=false, saturation_log=tmplog_path)
    @test isfile(tmplog_path)
    # read the file and ensure at least one line is present
    txt = read(tmplog_path, String)
    @test !isempty(strip(txt))
    Base.rm(tmplog_path)
end

@testset "expand_sites - stop on saturation" begin
    tmplog_path, tmplog_io = Base.mktemp()
    close(tmplog_io)
    sites = [Site(i, "R", ["A","B","C"], nothing) for i in 1:8]
    combos = expand_sites(sites; max_alts_per_site=3, max_combinations=5, sample_on_saturate=false, saturation_log=tmplog_path, stop_on_saturation=true)
    @test isempty(combos)
    @test isfile(tmplog_path)
    txt = read(tmplog_path, String)
    @test occursin("stop", txt)
    Base.rm(tmplog_path)
end
