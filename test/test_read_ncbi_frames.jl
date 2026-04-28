"""
test_read_ncbi_frames.jl — Tests for src/read_ncbi_frames.jl

Covers:
  parse_header           – typical header, join() syntax, description cleanup
  read_fasta_metadata    – single record, multiple records, blank-line handling
  read_fasta_sequence    – multi-line sequence, header skipped, strips whitespace
  extract_regions_sequence – single segment, multi-segment, boundary checks,
                             error on out-of-bounds
"""

using Test
using DataFrames

include(joinpath(@__DIR__, "../src/read_ncbi_frames.jl"))

# ===========================================================================
@testset "parse_header" begin

    @testset "typical NCBI header" begin
        region, desc = parse_header(">acc1:77..496 surface glycoprotein [HIV-1]")
        @test region == "77,496"
        @test occursin("surface glycoprotein", desc)
        @test !occursin("[HIV-1]", desc)   # bracket annotation stripped
    end

    @testset "header without description" begin
        region, desc = parse_header(">acc1:100..300")
        @test region == "100,300"
        @test strip(desc) == ""
    end

    @testset "join() multi-segment syntax" begin
        region, desc = parse_header(">join(acc1:77..496,acc2:606..980) polyprotein")
        parts = split(region, ";")
        @test length(parts) == 2
        @test parts[1] == "77,496"
        @test parts[2] == "606,980"
    end

    @testset "description leading pipe is stripped" begin
        _, desc = parse_header(">acc1:1..100 |some desc")
        @test !startswith(desc, "|")
    end

    @testset "bracket annotations removed from description" begin
        _, desc = parse_header(">acc1:1..100 gag [organism=HIV-1] [strain=HXB2]")
        @test !occursin("[organism=HIV-1]", desc)
        @test !occursin("[strain=HXB2]", desc)
    end

    @testset "returns two AbstractString values" begin
        # parse_header uses strip() which returns SubString{String}, a subtype
        # of AbstractString but not String itself.
        r, d = parse_header(">acc:1..10 desc")
        @test r isa AbstractString
        @test d isa AbstractString
    end

    @testset "single-position region (start==end)" begin
        region, _ = parse_header(">acc:50..50 single")
        @test region == "50,50"
    end

end

# ===========================================================================
@testset "read_fasta_metadata" begin

    function write_fasta(content)
        f, io = mktemp(); write(io, content); close(io); return f
    end

    @testset "single header returns 1-row DataFrame" begin
        f = write_fasta(">acc:1..30 gene\nATGCATGCATGCATGCATGCATGCATGCAT\n")
        df = read_fasta_metadata(f); rm(f)
        @test nrow(df) == 1
        @test df.Region[1] == "1,30"
    end

    @testset "multiple headers return multi-row DataFrame" begin
        content = ">acc1:1..30 geneA\nATGCATGCAT\n" *
                  ">acc2:50..80 geneB\nCCCGGGTTT\n"
        f = write_fasta(content); df = read_fasta_metadata(f); rm(f)
        @test nrow(df) == 2
        @test df.Region[1] == "1,30"
        @test df.Region[2] == "50,80"
    end

    @testset "returns DataFrame with Region and Description columns" begin
        f = write_fasta(">acc:1..10 x\nATGCATGCAT\n")
        df = read_fasta_metadata(f); rm(f)
        @test :Region in propertynames(df)
        @test :Description in propertynames(df)
    end

    @testset "blank lines between sequences handled" begin
        content = ">acc1:1..10 x\nATGCATGCAT\n\n>acc2:20..30 y\nCCCGGGTTTAA\n"
        f = write_fasta(content); df = read_fasta_metadata(f); rm(f)
        @test nrow(df) == 2
    end

    @testset "sequence lines (non-headers) are ignored" begin
        content = ">acc:1..12 gene\nATGCATGCATGC\nGGGGGGGGGGGG\n"
        f = write_fasta(content); df = read_fasta_metadata(f); rm(f)
        @test nrow(df) == 1
    end

end

# ===========================================================================
@testset "read_fasta_sequence" begin

    function write_fasta(content)
        f, io = mktemp(); write(io, content); close(io); return f
    end

    @testset "single-line sequence" begin
        f = write_fasta(">header\nATGCATGC\n")
        seq = read_fasta_sequence(f); rm(f)
        @test seq == "ATGCATGC"
    end

    @testset "multi-line sequence concatenated" begin
        f = write_fasta(">header\nATGC\nATGC\nATGC\n")
        seq = read_fasta_sequence(f); rm(f)
        @test seq == "ATGCATGCATGC"
    end

    @testset "header line excluded from sequence" begin
        f = write_fasta(">ATGCATGC\nGGGGGGGG\n")   # header looks like sequence
        seq = read_fasta_sequence(f); rm(f)
        @test seq == "GGGGGGGG"
    end

    @testset "strips whitespace from lines" begin
        f = write_fasta(">header\n  ATGC  \n  ATGC  \n")
        seq = read_fasta_sequence(f); rm(f)
        @test seq == "ATGCATGC"
    end

    @testset "empty file returns empty string" begin
        f, io = mktemp(); close(io)
        seq = read_fasta_sequence(f); rm(f)
        @test seq == ""
    end

end

# ===========================================================================
@testset "extract_regions_sequence" begin

    consensus = "ATGCCCGAATTTGGG"  # 15 nt

    @testset "single segment extracts correctly" begin
        @test extract_regions_sequence("1,6", consensus) == "ATGCCC"
    end

    @testset "single segment — end of sequence" begin
        @test extract_regions_sequence("13,15", consensus) == "GGG"
    end

    @testset "multi-segment concatenated in order" begin
        result = extract_regions_sequence("1,3;7,9", consensus)
        @test result == "ATGGAA"
    end

    @testset "1-based single nucleotide" begin
        @test extract_regions_sequence("4,4", consensus) == "C"
    end

    @testset "entire sequence" begin
        result = extract_regions_sequence("1,15", consensus)
        @test result == consensus
    end

    @testset "out-of-bounds throws error" begin
        @test_throws ErrorException extract_regions_sequence("1,100", consensus)
    end

    @testset "invalid region format throws error" begin
        @test_throws ErrorException extract_regions_sequence("1-6", consensus)
    end

    @testset "empty region string returns empty" begin
        @test extract_regions_sequence("", consensus) == ""
    end

end
