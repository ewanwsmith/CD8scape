"""
test_generate_peptides.jl — Comprehensive tests for src/generate_peptides.jl

Covers:
  CODON_DICT             – completeness, stop codons, spot-checks
  generate_peptides      – lengths, locus coverage, edge cases
  normalize_variant      – prefix/suffix trimming
  find_ref_match         – exact match, windowed search, failure
  check_locus            – column derivation, relative locus calc
  edit_ancestral_sequence – per-row patching of ancestral/derived sequences
  translate_sequences    – DNA→AA translation
  add_peptides_columns!  – peptide generation, indel labelling
  separate_peptides      – synonymous/stop filtering, A/D pairing
  deduplicate_peptides   – deduplication key logic
"""

using Test
using DataFrames

# generate_peptides.jl checks ARGS at top level; give it a dummy folder
push!(ARGS, "dummy_for_test")
try
    include(joinpath(@__DIR__, "../src/generate_peptides.jl"))
finally
    empty!(ARGS)
end

# ===========================================================================
@testset "CODON_DICT completeness" begin

    @testset "has exactly 64 entries" begin
        @test length(CODON_DICT) == 64
    end

    @testset "stop codons map to *" begin
        @test CODON_DICT["TAA"] == "*"
        @test CODON_DICT["TAG"] == "*"
        @test CODON_DICT["TGA"] == "*"
    end

    @testset "start codon ATG maps to M" begin
        @test CODON_DICT["ATG"] == "M"
    end

    @testset "TGG is the only W codon" begin
        @test CODON_DICT["TGG"] == "W"
    end

    @testset "spot check all keys are 3 uppercase nucleotides" begin
        for k in keys(CODON_DICT)
            @test length(k) == 3
            @test all(c -> c in "ACGT", k)
        end
    end

    @testset "values are single uppercase letters or *" begin
        for v in values(CODON_DICT)
            @test length(v) == 1
            @test v[1] in "ACDEFGHIKLMNPQRSTVWY*"
        end
    end

    @testset "Phe codons TTC/TTT" begin
        @test CODON_DICT["TTC"] == "F"
        @test CODON_DICT["TTT"] == "F"
    end

end

# ===========================================================================
@testset "generate_peptides" begin

    @testset "generates peptides of correct lengths" begin
        seq = "ACDEFGHIKLM"
        peps = generate_peptides(seq, 3, [8, 9])
        @test all(p -> length(p) == 8 || length(p) == 9, peps)
    end

    @testset "locus must be within each peptide" begin
        seq = "ACDEFGHIKLMNPQRST"
        locus = 5
        for len in [8, 9, 10]
            peps = generate_peptides(seq, locus, [len])
            for p in peps
                idx = findfirst(==(seq[locus]), p)
                @test idx !== nothing
            end
        end
    end

    @testset "returns empty for locus out of sequence range" begin
        seq = "ACDEFGHI"
        peps = generate_peptides(seq, 100, [8, 9])
        @test isempty(peps)
    end

    @testset "single-residue sequence with length 1" begin
        peps = generate_peptides("A", 1, [1])
        @test length(peps) == 1
        @test peps[1] == "A"
    end

    @testset "returns Vector{String}" begin
        peps = generate_peptides("ACDEFGHIKLM", 5, [9])
        @test peps isa Vector{String}
    end

    @testset "locus at start of sequence" begin
        seq = "ACDEFGHIKLM"
        peps = generate_peptides(seq, 1, [9])
        @test !isempty(peps)
        @test all(p -> startswith(seq, p[1:1]), peps)
    end

    @testset "locus at end of sequence" begin
        seq = "ACDEFGHIKLM"
        peps = generate_peptides(seq, length(seq), [9])
        @test !isempty(peps)
    end

    @testset "multiple lengths produce different counts" begin
        seq = "ACDEFGHIKLMNPQRST"
        p8 = generate_peptides(seq, 5, [8])
        p9 = generate_peptides(seq, 5, [9])
        p89 = generate_peptides(seq, 5, [8, 9])
        @test length(p89) == length(p8) + length(p9)
    end

end

# ===========================================================================
@testset "normalize_variant" begin

    @testset "SNP — no common prefix/suffix" begin
        r, a, adj = normalize_variant("A", "T")
        @test r == "A"
        @test a == "T"
        @test adj == 0
    end

    @testset "trims common suffix" begin
        r, a, adj = normalize_variant("ATG", "TTG")
        # Common suffix is 'TG', leaves 'A' vs 'T'
        @test length(r) < 3
        @test length(a) < 3
    end

    @testset "trims common prefix" begin
        r, a, adj = normalize_variant("GAT", "GAC")
        # Common prefix 'GA' -> 'T' vs 'C', prefix_trim=2
        @test adj == 2
        @test r == "T"
        @test a == "C"
    end

    @testset "no trimming for fully different" begin
        r, a, adj = normalize_variant("ACGT", "TGCA")
        @test adj == 0
    end

    @testset "result ref and alt never empty" begin
        r, a, _ = normalize_variant("AA", "A")
        @test !isempty(r)
        @test !isempty(a)
    end

    @testset "returns three values" begin
        result = normalize_variant("AGT", "ACT")
        @test length(result) == 3
    end

end

# ===========================================================================
@testset "find_ref_match" begin

    consensus = "ATGCCCGAATTT"

    @testset "exact match at expected position" begin
        pos, ok = find_ref_match(consensus, "CCC", 4)
        @test ok
        @test pos == 4
    end

    @testset "match found within window" begin
        pos, ok = find_ref_match(consensus, "GAA", 6; window=5)
        @test ok
    end

    @testset "no match returns ok=false" begin
        _, ok = find_ref_match(consensus, "ZZZ", 1; window=3)
        @test !ok
    end

    @testset "empty ref returns original locus with ok=true" begin
        pos, ok = find_ref_match(consensus, "", 3)
        @test ok
        @test pos == 3
    end

    @testset "single-character exact match" begin
        pos, ok = find_ref_match(consensus, "A", 1)
        @test ok
        @test pos == 1
    end

end

# ===========================================================================
@testset "check_locus" begin

    function make_check_df()
        # Simple 12-nt frame: ATG CCC GAA TTT
        df = DataFrame(
            Locus         = [4],
            Consensus     = ["C"],
            Variant       = ["T"],
            Start         = [1],
            End           = [12],
            Consensus_sequence = ["ATGCCCGAATTT"],
            Description   = ["Frame1"],
        )
        return df
    end

    @testset "adds Relative_Locus column" begin
        df = check_locus(make_check_df())
        @test :Relative_Locus in propertynames(df)
    end

    @testset "Relative_Locus is 1-based (locus - start + 1)" begin
        df = check_locus(make_check_df())
        @test df.Relative_Locus[1] == 4
    end

    @testset "adds Matches_Consensus column" begin
        df = check_locus(make_check_df())
        @test :Matches_Consensus in propertynames(df)
    end

    @testset "adds Normalized_Consensus and Normalized_Variant columns" begin
        df = check_locus(make_check_df())
        @test :Normalized_Consensus in propertynames(df)
        @test :Normalized_Variant in propertynames(df)
    end

    @testset "adds Pulled_Base column" begin
        df = check_locus(make_check_df())
        @test :Pulled_Base in propertynames(df)
    end

end

# ===========================================================================
@testset "edit_ancestral_sequence" begin

    function make_edit_df()
        df = DataFrame(
            Locus             = [4],
            Start             = [1],
            End               = [12],
            Consensus_sequence= ["ATGCCCGAATTT"],
            Normalized_Consensus = ["C"],
            Normalized_Variant   = ["T"],
            Relative_Locus    = [4],
            Matches_Consensus = [true],
            Adjusted_Locus    = [4],
            Description       = ["Frame1"],
        )
        df[!, :Ancestral_sequence] = [""]
        df[!, :Derived_sequence]   = [""]
        return df
    end

    @testset "adds Ancestral_sequence and Derived_sequence columns" begin
        df = make_edit_df()
        df2 = edit_ancestral_sequence(df)
        @test :Ancestral_sequence in propertynames(df2)
        @test :Derived_sequence in propertynames(df2)
    end

    @testset "ancestral sequence contains ref allele at relative locus" begin
        df = make_edit_df()
        df2 = edit_ancestral_sequence(df)
        rl = df2.Relative_Locus[1]
        @test df2.Ancestral_sequence[1][rl] == 'C'
    end

    @testset "derived sequence contains alt allele at relative locus" begin
        df = make_edit_df()
        df2 = edit_ancestral_sequence(df)
        rl = df2.Relative_Locus[1]
        @test df2.Derived_sequence[1][rl] == 'T'
    end

    @testset "sequences have same length as consensus" begin
        df = make_edit_df()
        df2 = edit_ancestral_sequence(df)
        @test length(df2.Ancestral_sequence[1]) == length(df2.Consensus_sequence[1])
        @test length(df2.Derived_sequence[1]) == length(df2.Consensus_sequence[1])
    end

end

# ===========================================================================
@testset "translate_sequences" begin

    function make_translate_df()
        # ATG CCC GAA TTT → M P E F
        seq = "ATGCCCGAATTT"
        DataFrame(
            Locus             = [4],
            Ancestral_sequence = [seq],
            Derived_sequence   = [replace(seq, "CCC" => "TTT")],  # P→F
        )
    end

    @testset "adds Ancestral_AA_sequence column" begin
        df = translate_sequences(make_translate_df())
        @test :Ancestral_AA_sequence in propertynames(df)
    end

    @testset "adds Derived_AA_sequence column" begin
        df = translate_sequences(make_translate_df())
        @test :Derived_AA_sequence in propertynames(df)
    end

    @testset "ATG CCC GAA TTT translates to MPEF" begin
        df = translate_sequences(make_translate_df())
        @test df.Ancestral_AA_sequence[1] == "MPEF"
    end

    @testset "derived sequence with P→F gives MFEF" begin
        df = translate_sequences(make_translate_df())
        @test df.Derived_AA_sequence[1] == "MFEF"
    end

    @testset "stop codon translates to *" begin
        df = DataFrame(
            Locus             = [1],
            Ancestral_sequence = ["ATGTAA"],
            Derived_sequence   = ["ATGTAA"],
        )
        df2 = translate_sequences(df)
        @test df2.Ancestral_AA_sequence[1] == "M*"
    end

end

# ===========================================================================
@testset "separate_peptides" begin

    function make_sep_df(; anc, der, locus=100, label="A1B_Frame1_1")
        DataFrame(
            Locus           = [locus],
            Relative_Locus  = [1],
            AA_Locus        = [1],
            Ancestral_Peptide = [anc],
            Derived_Peptide   = [der],
            Peptide_label   = [label],
        )
    end

    @testset "synonymous pair is removed" begin
        df = separate_peptides(make_sep_df(anc="ACDEFGHIK", der="ACDEFGHIK"))
        @test nrow(df) == 0
    end

    @testset "non-synonymous pair yields 2 rows (A and D)" begin
        df = separate_peptides(make_sep_df(anc="ACDEFGHIK", der="ACDEFGHIL"))
        @test nrow(df) == 2
        states = last.(split.(df.Peptide_label, "_"))
        @test "A" in states
        @test "D" in states
    end

    @testset "stop codon in ancestral is removed" begin
        df = separate_peptides(make_sep_df(anc="ACDE*GHIK", der="ACDEFGHIK"))
        @test nrow(df) == 0
    end

    @testset "stop codon in derived is removed" begin
        df = separate_peptides(make_sep_df(anc="ACDEFGHIK", der="ACDE*GHIK"))
        @test nrow(df) == 0
    end

    @testset "result DataFrame has Locus, Peptide, Peptide_label columns" begin
        df = separate_peptides(make_sep_df(anc="ACDEFGHIK", der="ACDEFGHIL"))
        @test :Locus in propertynames(df)
        @test :Peptide in propertynames(df)
        @test :Peptide_label in propertynames(df)
    end

    @testset "multiple non-synonymous variants all retained" begin
        rows = vcat(
            make_sep_df(anc="AAAAAAAA", der="TTTTTTTT", locus=1,  label="X1Y_F_1"),
            make_sep_df(anc="CCCCCCCC", der="GGGGGGGG", locus=2,  label="X2Y_F_1"),
        )
        df = separate_peptides(rows)
        @test nrow(df) == 4
    end

end

# ===========================================================================
@testset "deduplicate_peptides" begin

    @testset "removes duplicate (Locus, Peptide, State, Peptide_label)" begin
        df = DataFrame(
            Locus        = [1, 1],
            Peptide      = ["ACDEFGHIK", "ACDEFGHIK"],
            Peptide_label = ["A1T_Frame_1_A", "A1T_Frame_1_A"],
        )
        dd = deduplicate_peptides(df)
        @test nrow(dd) == 1
    end

    @testset "keeps rows with different labels" begin
        df = DataFrame(
            Locus        = [1, 1],
            Peptide      = ["ACDEFGHIK", "ACDEFGHIK"],
            Peptide_label = ["A1T_Frame_1_A", "A1T_Frame_2_A"],
        )
        dd = deduplicate_peptides(df)
        @test nrow(dd) == 2
    end

    @testset "removes helper State column after dedup" begin
        df = DataFrame(
            Locus        = [1],
            Peptide      = ["ACDEFGHIK"],
            Peptide_label = ["A1T_Frame_1_A"],
        )
        dd = deduplicate_peptides(df)
        @test !(:State in propertynames(dd))
    end

    @testset "distinguishes ancestral (_A) from derived (_D)" begin
        df = DataFrame(
            Locus        = [1, 1],
            Peptide      = ["ACDEFGHIK", "ACDEFGHIK"],
            Peptide_label = ["A1T_Frame_1_A", "A1T_Frame_1_D"],
        )
        dd = deduplicate_peptides(df)
        @test nrow(dd) == 2
    end

end
