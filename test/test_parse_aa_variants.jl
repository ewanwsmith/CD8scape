"""
test_parse_aa_variants.jl — Tests for src/parse_aa_variants.jl

Covers:
  CODON_DICT           – completeness (shared with generate_peptides)
  build_aa_to_codon    – reverse mapping, uniqueness, no stop codons
  AA_TO_CODON          – spot-checks
  parse_aa_file        – normal records, blank-line separation, malformed
                         lines, invalid aa_pos, multi-word orf names
  frame_start          – single-segment, multi-segment region strings
"""

using Test

include(joinpath(@__DIR__, "../src/parse_aa_variants.jl"))

# ===========================================================================
@testset "build_aa_to_codon / AA_TO_CODON" begin

    @testset "no stop codon in reverse map" begin
        @test !haskey(AA_TO_CODON, "*")
    end

    @testset "all 20 standard amino acids present" begin
        standard = collect("ACDEFGHIKLMNPQRSTVWY")
        for aa in standard
            @test haskey(AA_TO_CODON, string(aa))
        end
    end

    @testset "values are valid 3-nt codons" begin
        for (aa, codon) in AA_TO_CODON
            @test length(codon) == 3
            @test all(c -> c in "ACGT", codon)
        end
    end

    @testset "each codon translates back to its amino acid" begin
        for (aa, codon) in AA_TO_CODON
            @test get(CODON_DICT, codon, "?") == aa
        end
    end

    @testset "M maps to ATG" begin
        @test AA_TO_CODON["M"] == "ATG"
    end

    @testset "W maps to TGG" begin
        @test AA_TO_CODON["W"] == "TGG"
    end

end

# ===========================================================================
@testset "parse_aa_file" begin

    function write_aa(content::String)
        f, io = mktemp()
        write(io, content)
        close(io)
        return f
    end

    @testset "parses single well-formed record" begin
        f = write_aa("ORF1 10\nK M\n")
        records = parse_aa_file(f)
        rm(f)
        @test length(records) == 1
        @test records[1].orf_name == "ORF1"
        @test records[1].aa_pos   == 10
        @test records[1].ancestral_aa == "K"
        @test records[1].derived_aa   == "M"
    end

    @testset "parses two records separated by blank line" begin
        f = write_aa("ORF1 5\nA T\n\nORF2 12\nG V\n")
        records = parse_aa_file(f)
        rm(f)
        @test length(records) == 2
        @test records[2].orf_name == "ORF2"
        @test records[2].aa_pos   == 12
    end

    @testset "multi-word orf name" begin
        f = write_aa("Capsid Protein 7\nI L\n")
        records = parse_aa_file(f)
        rm(f)
        @test records[1].orf_name == "Capsid Protein"
        @test records[1].aa_pos   == 7
    end

    @testset "lowercase aa codes are uppercased" begin
        f = write_aa("ORF1 3\nk m\n")
        records = parse_aa_file(f)
        rm(f)
        @test records[1].ancestral_aa == "K"
        @test records[1].derived_aa   == "M"
    end

    @testset "empty file returns empty vector" begin
        f = write_aa("")
        records = parse_aa_file(f)
        rm(f)
        @test isempty(records)
    end

    @testset "invalid aa_pos (zero) is skipped" begin
        f = write_aa("ORF1 0\nK M\n")
        records = parse_aa_file(f)
        rm(f)
        @test isempty(records)
    end

    @testset "invalid aa_pos (non-numeric) is skipped" begin
        f = write_aa("ORF1 abc\nK M\n")
        records = parse_aa_file(f)
        rm(f)
        @test isempty(records)
    end

    @testset "multi-character aa codes are skipped" begin
        f = write_aa("ORF1 5\nLys Met\n")
        records = parse_aa_file(f)
        rm(f)
        @test isempty(records)
    end

    @testset "multiple records all parsed" begin
        lines = join(["ORF$(i) $(i)\nA T\n" for i in 1:5], "\n")
        f = write_aa(lines)
        records = parse_aa_file(f)
        rm(f)
        @test length(records) == 5
        @test [r.aa_pos for r in records] == collect(1:5)
    end

    @testset "records returned as NamedTuples with correct fields" begin
        f = write_aa("ORF1 1\nA T\n")
        records = parse_aa_file(f)
        rm(f)
        r = records[1]
        @test hasproperty(r, :orf_name)
        @test hasproperty(r, :aa_pos)
        @test hasproperty(r, :ancestral_aa)
        @test hasproperty(r, :derived_aa)
    end

end

# ===========================================================================
@testset "frame_start" begin

    @testset "single segment '77,496'" begin
        @test frame_start("77,496") == 77
    end

    @testset "multi-segment '77,496;606,980'" begin
        @test frame_start("77,496;606,980") == 77
    end

    @testset "start of 1 (gene starts at position 1)" begin
        @test frame_start("1,300") == 1
    end

    @testset "large position" begin
        @test frame_start("10000,20000") == 10000
    end

    @testset "three segments, first start used" begin
        @test frame_start("50,100;200,300;400,500") == 50
    end

end
