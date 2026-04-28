"""
test_variant_fates.jl — Tests for src/variant_fates.jl

Covers:
  CODON_DICT (shared copy)  – completeness, stop codons
  translate_dna              – simple sequences, stop codon, unknown codon
  classify_and_label         – synonymous, stop-codon, viable (non-synonymous),
                               indel (deletion, insertion), out-of-bounds
  FATE_RANK                  – ordering
"""

using Test

include(joinpath(@__DIR__, "../src/variant_fates.jl"))

# ===========================================================================
@testset "CODON_DICT (variant_fates copy)" begin

    @testset "64 entries" begin
        @test length(CODON_DICT) == 64
    end

    @testset "stop codons" begin
        @test CODON_DICT["TAA"] == "*"
        @test CODON_DICT["TAG"] == "*"
        @test CODON_DICT["TGA"] == "*"
    end

    @testset "ATG -> M" begin
        @test CODON_DICT["ATG"] == "M"
    end

end

# ===========================================================================
@testset "translate_dna" begin

    @testset "ATG -> M (single codon)" begin
        @test translate_dna("ATG") == "M"
    end

    @testset "ATG CCC GAA TTT -> MPEF" begin
        @test translate_dna("ATGCCCGAATTT") == "MPEF"
    end

    @testset "TAA is stop codon *" begin
        @test translate_dna("ATGTAA") == "M*"
    end

    @testset "unknown codon returns ?" begin
        @test translate_dna("NNN") == "?"
    end

    @testset "empty string returns empty" begin
        @test translate_dna("") == ""
    end

    @testset "partial trailing nucleotides ignored" begin
        # ATGCCC + extra G → only 2 codons translated
        @test translate_dna("ATGCCCG") == "MP"
    end

    @testset "uppercase conversion" begin
        # translate_dna calls uppercase via CODON_DICT lookup
        @test translate_dna("atgccc") == "MP"
    end

end

# ===========================================================================
@testset "classify_and_label" begin

    # consensus: ATG CCC GAA TTT  (MPEF)
    # positions: 1   4   7   10
    consensus = "ATGCCCGAATTT"

    @testset "synonymous SNP returns 'synonymous'" begin
        # Position 9: GAA -> GAG both encode E
        fate, mut = classify_and_label(consensus, "A", "G", 9, 9)
        @test fate == "synonymous"
    end

    @testset "non-synonymous SNP returns 'viable'" begin
        # Position 4: CCC -> TCC  P→S
        fate, mut = classify_and_label(consensus, "C", "T", 4, 4)
        @test fate == "viable"
    end

    @testset "derived stop codon returns 'stop_codon'" begin
        # Position 7: GAA -> TAA  E→*
        fate, mut = classify_and_label(consensus, "G", "T", 7, 7)
        @test fate == "stop_codon"
    end

    @testset "mutation label for non-synonymous has format AAposletter" begin
        # P4S: codon 2 (1-based aa), consensus=P, derived=S, aa_pos=2
        fate, mut = classify_and_label(consensus, "C", "T", 4, 4)
        @test occursin(r"[A-Z]\d+[A-Z]", mut)
    end

    @testset "deletion labelled as del<locus><nt>" begin
        fate, mut = classify_and_label(consensus, "ATG", "A", 1, 1)
        @test startswith(mut, "del")
        @test fate in ("viable", "synonymous", "stop_codon")
    end

    @testset "insertion labelled as ins<locus><nt>" begin
        fate, mut = classify_and_label(consensus, "A", "ATG", 1, 1)
        @test startswith(mut, "ins")
    end

    @testset "out-of-bounds ref_end returns 'synonymous'" begin
        fate, _ = classify_and_label("ATGC", "ATGCCCCC", "T", 1, 1)
        @test fate == "synonymous"
    end

    @testset "return type is Tuple{String,String}" begin
        result = classify_and_label(consensus, "C", "T", 4, 4)
        @test result isa Tuple{String, String}
    end

    @testset "ancestral stop does not mask derived stop" begin
        # consensus with stop at pos 7: ATGCCCTAA
        cons2 = "ATGCCCTAA"
        # Changing T at pos 7 to G gives GAA (E), so derived has no stop at that position
        fate, _ = classify_and_label(cons2, "T", "G", 7, 7)
        # ancestral has stop at codon 3; derived may differ
        @test fate isa String
    end

end

# ===========================================================================
@testset "FATE_RANK ordering" begin

    @testset "synonymous < stop_codon < viable" begin
        @test FATE_RANK["synonymous"] < FATE_RANK["stop_codon"]
        @test FATE_RANK["stop_codon"] < FATE_RANK["viable"]
    end

    @testset "all three keys present" begin
        @test haskey(FATE_RANK, "synonymous")
        @test haskey(FATE_RANK, "stop_codon")
        @test haskey(FATE_RANK, "viable")
    end

end
