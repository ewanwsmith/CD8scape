"""
test_clean_peptides.jl — Tests for src/clean_peptides.jl

Covers:
  process_peptides – duplicate removal, stop-codon filtering,
                     file not found, case-insensitive filename matching
"""

using Test
using CSV
using DataFrames

include(joinpath(@__DIR__, "../src/clean_peptides.jl"))

# ---------------------------------------------------------------------------
function with_tmpdir(f::Function)
    dir = mktempdir()
    try; f(dir); finally; rm(dir; recursive=true, force=true); end
end

function write_pep(dir, lines; filename="Peptides.pep")
    path = joinpath(dir, filename)
    open(path, "w") do io
        println(io, "x")       # CSV header row (column name "x")
        for l in lines; println(io, l); end
    end
    return path
end

# ===========================================================================
@testset "process_peptides" begin

    @testset "removes duplicate peptides" begin
        with_tmpdir() do dir
            write_pep(dir, ["ACDEFGHIK", "ACDEFGHIK", "MNPQRSTVW"])
            process_peptides(dir)
            df = CSV.read(joinpath(dir, "Peptides.pep"), DataFrame)
            @test nrow(df) == 2
        end
    end

    @testset "removes peptides containing stop codon (*)" begin
        with_tmpdir() do dir
            write_pep(dir, ["ACDE*GHIK", "MNPQRSTVW"])
            process_peptides(dir)
            df = CSV.read(joinpath(dir, "Peptides.pep"), DataFrame)
            peps = string.(df[:, 1])
            @test !any(occursin("*", p) for p in peps)
        end
    end

    @testset "keeps valid peptides unchanged" begin
        with_tmpdir() do dir
            write_pep(dir, ["ACDEFGHIK", "MNPQRSTVWY"])
            process_peptides(dir)
            df = CSV.read(joinpath(dir, "Peptides.pep"), DataFrame)
            @test nrow(df) == 2
        end
    end

    @testset "all stop-codon peptides results in empty file" begin
        with_tmpdir() do dir
            write_pep(dir, ["A*CDEFG", "*MNPQRS"])
            process_peptides(dir)
            df = CSV.read(joinpath(dir, "Peptides.pep"), DataFrame)
            @test nrow(df) == 0
        end
    end

    @testset "all duplicates results in single row" begin
        with_tmpdir() do dir
            write_pep(dir, ["ACDEFGHIK", "ACDEFGHIK", "ACDEFGHIK"])
            process_peptides(dir)
            df = CSV.read(joinpath(dir, "Peptides.pep"), DataFrame)
            @test nrow(df) == 1
        end
    end

    @testset "missing Peptides.pep prints error and returns" begin
        with_tmpdir() do dir
            # No pep file — function should return without crashing
            @test_nowarn process_peptides(dir)
        end
    end

    @testset "case-insensitive filename match (PEPTIDES.PEP)" begin
        with_tmpdir() do dir
            write_pep(dir, ["ACDEFGHIK"]; filename="PEPTIDES.PEP")
            # Should not throw even if casing differs
            @test_nowarn process_peptides(dir)
        end
    end

    @testset "duplicate removal + stop codon removal combined" begin
        with_tmpdir() do dir
            write_pep(dir, ["ACDEFGHIK", "ACDEFGHIK", "A*CDEFGHI", "MNPQRSTVW"])
            process_peptides(dir)
            df = CSV.read(joinpath(dir, "Peptides.pep"), DataFrame)
            peps = string.(df[:, 1])
            @test length(peps) == 2
            @test !any(occursin("*", p) for p in peps)
        end
    end

    @testset "single valid peptide preserved" begin
        with_tmpdir() do dir
            write_pep(dir, ["ACDEFGHIK"])
            process_peptides(dir)
            df = CSV.read(joinpath(dir, "Peptides.pep"), DataFrame)
            @test nrow(df) == 1
        end
    end

end
