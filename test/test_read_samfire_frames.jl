"""
test_read_samfire_frames.jl — Tests for src/read_samfire_frames.jl

Covers:
  read_dat_file         – normal 3-line blocks, multiple blocks,
                          incomplete trailing block, bad start/end lines
  save_dataframe_as_csv – writes CSV to correct path
"""

using Test
using DataFrames
using CSV

include(joinpath(@__DIR__, "../src/read_samfire_frames.jl"))

# ---------------------------------------------------------------------------
function with_tmpdir(f::Function)
    dir = mktempdir()
    try; f(dir); finally; rm(dir; recursive=true, force=true); end
end

function write_dat(content::String)
    f, io = mktemp(); write(io, content); close(io); return f
end

# ===========================================================================
@testset "read_dat_file" begin

    @testset "single valid block produces 1-row DataFrame" begin
        f = write_dat("77 496\nATGCCCGAATTT\nMPEF\n")
        df = read_dat_file(f); rm(f)
        @test nrow(df) == 1
        @test df.Region[1] == "77,496"
        @test df.Consensus_sequence[1] == "ATGCCCGAATTT"
        @test df.Description[1] == "Frame_1"
    end

    @testset "two blocks produce 2-row DataFrame" begin
        content = "77 496\nATGCCC\nMPE\n" *
                  "606 980\nGAATTT\nEF\n"
        f = write_dat(content); df = read_dat_file(f); rm(f)
        @test nrow(df) == 2
        @test df.Region[2] == "606,980"
    end

    @testset "Description auto-labelled Frame_1, Frame_2, ..." begin
        content = "1 100\nAAA\nK\n" *
                  "200 300\nCCC\nP\n" *
                  "400 500\nGGG\nG\n"
        f = write_dat(content); df = read_dat_file(f); rm(f)
        @test df.Description[1] == "Frame_1"
        @test df.Description[2] == "Frame_2"
        @test df.Description[3] == "Frame_3"
    end

    @testset "output has Region, Consensus_sequence, Description columns" begin
        f = write_dat("1 100\nATGC\nM\n")
        df = read_dat_file(f); rm(f)
        @test :Region in propertynames(df)
        @test :Consensus_sequence in propertynames(df)
        @test :Description in propertynames(df)
    end

    @testset "incomplete trailing block is skipped (no error)" begin
        content = "1 100\nATGC\nM\n" *
                  "200 300\n"          # only 2 of 3 lines present
        f = write_dat(content); df = read_dat_file(f); rm(f)
        @test nrow(df) == 1
    end

    @testset "block with non-integer start/end is skipped" begin
        content = "one hundred\nATGC\nM\n" *
                  "1 100\nCCCC\nP\n"
        f = write_dat(content); df = read_dat_file(f); rm(f)
        @test nrow(df) == 1
        @test df.Region[1] == "1,100"
    end

    @testset "region formatted as 'start,end'" begin
        f = write_dat("10 50\nATGC\nM\n")
        df = read_dat_file(f); rm(f)
        @test df.Region[1] == "10,50"
    end

    @testset "empty file returns empty DataFrame" begin
        f = write_dat("")
        df = read_dat_file(f); rm(f)
        @test nrow(df) == 0
    end

    @testset "returns a DataFrame" begin
        f = write_dat("1 100\nATGC\nM\n")
        df = read_dat_file(f); rm(f)
        @test df isa DataFrame
    end

end

# ===========================================================================
@testset "save_dataframe_as_csv" begin

    @testset "writes CSV to directory of input file" begin
        with_tmpdir() do dir
            input_path = joinpath(dir, "Reading_Frames.dat")
            write(input_path, "dummy\n")
            df = DataFrame(Region=["1,100"], Consensus_sequence=["ATGC"], Description=["Frame_1"])
            save_dataframe_as_csv(df, input_path)
            out = joinpath(dir, "frames.csv")
            @test isfile(out)
        end
    end

    @testset "written CSV is readable and correct" begin
        with_tmpdir() do dir
            input_path = joinpath(dir, "Reading_Frames.dat")
            write(input_path, "dummy\n")
            df = DataFrame(Region=["1,100","200,300"],
                           Consensus_sequence=["ATGC","CCCC"],
                           Description=["Frame_1","Frame_2"])
            save_dataframe_as_csv(df, input_path)
            df2 = CSV.read(joinpath(dir, "frames.csv"), DataFrame)
            @test nrow(df2) == 2
            @test df2.Region[1] == "1,100"
        end
    end

    @testset "custom output filename honoured" begin
        with_tmpdir() do dir
            input_path = joinpath(dir, "Reading_Frames.dat")
            write(input_path, "dummy\n")
            df = DataFrame(Region=["1,100"], Consensus_sequence=["ATGC"], Description=["F1"])
            save_dataframe_as_csv(df, input_path, "custom_frames.csv")
            @test isfile(joinpath(dir, "custom_frames.csv"))
        end
    end

end
