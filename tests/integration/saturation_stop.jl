
include("../../src/variant_expansion.jl")
using .VariantExpansion

using Test, CSV, DataFrames

@testset "expand_sites - saturation logging" begin
    # tmp log path
    logpath = Base.mktemp()[1]
    Base.rm(logpath) # let expand_sites create it

    # build 6 sites with REF + 4 ALTs -> 5^6 combos (15625)
    sites = [Site(i, "R", ["A","B","C","D"], nothing) for i in 1:6]

    combos = expand_sites(sites; max_alts_per_site=4, max_combinations=1000, stop_on_saturation=true, saturation_log=logpath)

    @test isempty(combos)
    @test isfile(logpath)

    # Parse the saturation log CSV while ignoring comment lines starting with '#'
    txt = read(logpath, String)
    io = IOBuffer(txt)
    # CSV.File supports a `comment` keyword to skip comment lines
    df = DataFrame(CSV.File(io; comment="#"))
    @test nrow(df) >= 1
    row = df[end, :]

    # Strict assertions for fields
    @test parse(Int, string(row[:n_sites])) == 6
    @test parse(Int, string(row[:product])) == 15625
    @test parse(Int, string(row[:threshold])) == 1000
    @test string(row[:strategy]) == "stop"
    @test string(row[:site_positions]) == "1;2;3;4;5;6"

    # cleanup
    Base.rm(logpath)
end
