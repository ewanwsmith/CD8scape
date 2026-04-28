"""
test_process_scores.jl — Tests for src/process_scores.jl

Covers:
  reshape_hla_data     – renames HLA → MHC, drops original column
  process_and_join     – correct join, missing peptide reporting,
                         handles absent Locus column gracefully
"""

using Test
using DataFrames
using CSV

include(joinpath(@__DIR__, "../src/process_scores.jl"))

# ---------------------------------------------------------------------------
function with_tmpdir(f::Function)
    dir = mktempdir()
    try; f(dir); finally; rm(dir; recursive=true, force=true); end
end

# ===========================================================================
@testset "reshape_hla_data" begin

    @testset "renames HLA to MHC" begin
        df = DataFrame(Peptide=["AAA"], EL_Rank=[0.5], HLA=["HLA-A02:01"])
        result = reshape_hla_data(df)
        @test :MHC in propertynames(result)
        @test !(:HLA in propertynames(result))
    end

    @testset "MHC values preserved" begin
        df = DataFrame(Peptide=["AAA"], EL_Rank=[0.5], HLA=["HLA-A02:01"])
        result = reshape_hla_data(df)
        @test result.MHC[1] == "HLA-A02:01"
    end

    @testset "other columns untouched" begin
        df = DataFrame(Peptide=["PPP"], EL_Rank=[1.0], HLA=["HLA-B07:02"],
                       Extra=["x"])
        result = reshape_hla_data(df)
        @test :Peptide in propertynames(result)
        @test :EL_Rank in propertynames(result)
        @test :Extra in propertynames(result)
    end

    @testset "multiple rows all renamed" begin
        df = DataFrame(Peptide=["AAA","BBB"], EL_Rank=[0.5,1.0],
                       HLA=["HLA-A02:01","HLA-B07:02"])
        result = reshape_hla_data(df)
        @test nrow(result) == 2
        @test all(!ismissing, result.MHC)
    end

end

# ===========================================================================
@testset "process_and_join" begin

    function make_test_data(dir)
        # Minimal processed_output.csv (NetMHCpan-like output)
        mhc = DataFrame(
            Pos     = [1, 2, 3],
            Peptide = ["ACDEFGHIK", "MNPQRSTVW", "UNKNOWN___"],
            EL_Rank = [0.1, 0.5, 2.0],
            HLA     = ["HLA-A02:01", "HLA-A02:01", "HLA-A02:01"],
        )
        CSV.write(joinpath(dir, "processed_output.csv"), mhc)

        # Minimal peptides_labels.csv
        labels = DataFrame(
            Peptide      = ["ACDEFGHIK", "MNPQRSTVW"],
            Peptide_label = ["A1B_ORF1_1_A", "C2D_ORF1_1_D"],
            Locus        = [100, 200],
        )
        CSV.write(joinpath(dir, "peptides_labels.csv"), labels)
    end

    @testset "returns a DataFrame" begin
        with_tmpdir() do dir
            make_test_data(dir)
            df = process_and_join(dir)
            @test df isa DataFrame
        end
    end

    @testset "joined result has Locus column from labels" begin
        with_tmpdir() do dir
            make_test_data(dir)
            df = process_and_join(dir)
            @test :Locus in propertynames(df)
        end
    end

    @testset "matched peptides have non-missing Pos" begin
        with_tmpdir() do dir
            make_test_data(dir)
            df = process_and_join(dir)
            matched = filter(row -> row.Peptide in ["ACDEFGHIK", "MNPQRSTVW"], df)
            @test all(!ismissing, matched.Pos)
        end
    end

    @testset "unmatched peptide has missing Locus" begin
        with_tmpdir() do dir
            make_test_data(dir)
            df = process_and_join(dir)
            unmatched = filter(row -> row.Peptide == "UNKNOWN___", df)
            @test nrow(unmatched) == 1
            @test ismissing(unmatched.Locus[1])
        end
    end

    @testset "throws when required files missing" begin
        with_tmpdir() do dir
            @test_throws ErrorException process_and_join(dir)
        end
    end

    @testset "all rows of processed_output preserved (left join)" begin
        with_tmpdir() do dir
            make_test_data(dir)
            df = process_and_join(dir)
            @test nrow(df) == 3
        end
    end

end
