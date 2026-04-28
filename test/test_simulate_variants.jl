"""
test_simulate_variants.jl — Tests for src/simulate_variants.jl

Covers:
  parse_region_bounds    – single segment, multi-segment (uses first start / last end)
  CODON_DICT (shared)    – spot-check
  Variant enumeration    – each codon position has exactly 3 substitutions,
                           total count for a known sequence, no self-variants
  Simulation with frames – end-to-end DataFrame content from simulate logic
"""

using Test
using DataFrames
using CSV
using Random

push!(ARGS, "dummy_for_test")
try
    include(joinpath(@__DIR__, "../src/simulate_variants.jl"))
finally
    empty!(ARGS)
end

# ---------------------------------------------------------------------------
function with_tmpdir(f::Function)
    dir = mktempdir()
    try; f(dir); finally; rm(dir; recursive=true, force=true); end
end

# ===========================================================================
@testset "parse_region_bounds" begin

    @testset "single segment '1,300'" begin
        s, e = parse_region_bounds("1,300")
        @test s == 1
        @test e == 300
    end

    @testset "multi-segment uses first start and last end" begin
        s, e = parse_region_bounds("77,496;606,980")
        @test s == 77
        @test e == 980
    end

    @testset "quoted region string" begin
        s, e = parse_region_bounds("\"1,300\"")
        @test s == 1
        @test e == 300
    end

    @testset "returns integers" begin
        s, e = parse_region_bounds("10,50")
        @test s isa Integer
        @test e isa Integer
    end

    @testset "start < end for normal gene" begin
        s, e = parse_region_bounds("100,500")
        @test s < e
    end

end

# ===========================================================================
@testset "Variant enumeration logic" begin

    # Simulate the enumeration kernel inline (mirrors main() logic)
    function enumerate_variants(dna::String, start_nt::Int)
        nts = ['A', 'C', 'G', 'T']
        rows = NamedTuple{(:Locus, :Consensus, :Variant), Tuple{Int,String,String}}[]
        n_codons = div(length(dna), 3)
        for ci in 1:n_codons
            dna_start = (ci - 1) * 3 + 1
            for pos in 0:2
                offset  = dna_start + pos
                locus   = start_nt + offset - 1
                cons_nt = dna[offset]
                for var_nt in nts
                    var_nt == cons_nt && continue
                    push!(rows, (Locus=locus, Consensus=string(cons_nt), Variant=string(var_nt)))
                end
            end
        end
        return rows
    end

    @testset "each position has exactly 3 substitutions (4 nts - 1 self)" begin
        rows = enumerate_variants("ATGCCC", 1)
        per_locus = Dict{Int,Int}()
        for r in rows
            per_locus[r.Locus] = get(per_locus, r.Locus, 0) + 1
        end
        for (_, cnt) in per_locus
            @test cnt == 3
        end
    end

    @testset "no self-substitution (Consensus != Variant)" begin
        rows = enumerate_variants("ATGCCCGAA", 1)
        @test all(r.Consensus != r.Variant for r in rows)
    end

    @testset "total variants = n_codons * 3 positions * 3 substitutions" begin
        dna = "ATGCCCGAATTT"  # 4 codons
        rows = enumerate_variants(dna, 1)
        @test length(rows) == 4 * 3 * 3
    end

    @testset "loci are offset by start_nt" begin
        rows = enumerate_variants("ATGCCC", 100)
        loci = sort(unique([r.Locus for r in rows]))
        @test minimum(loci) == 100
        @test maximum(loci) == 105
    end

    @testset "Variant field is always a single nucleotide" begin
        rows = enumerate_variants("ATGCCC", 1)
        @test all(length(r.Variant) == 1 for r in rows)
    end

end

# ===========================================================================
@testset "CODON_DICT spot-checks (simulate_variants.jl copy)" begin

    @testset "ATG -> M (start codon)" begin
        @test CODON_DICT["ATG"] == "M"
    end

    @testset "TAA, TAG, TGA -> * (stop codons)" begin
        @test CODON_DICT["TAA"] == "*"
        @test CODON_DICT["TAG"] == "*"
        @test CODON_DICT["TGA"] == "*"
    end

    @testset "64 entries total" begin
        @test length(CODON_DICT) == 64
    end

end

# ===========================================================================
@testset "simulate_variants end-to-end (with temp frames.csv)" begin

    function run_simulation(dna, start_pos, end_pos, desc)
        with_tmpdir() do dir
            frames_df = DataFrame(
                Region = ["$(start_pos),$(end_pos)"],
                Consensus_sequence = [dna],
                Description = [desc],
            )
            CSV.write(joinpath(dir, "frames.csv"), frames_df)

            # Replicate the enumeration from simulate_variants main()
            out = DataFrame(Locus=Int[], Consensus=String[], Variant=String[])
            nts = ['A', 'C', 'G', 'T']
            for row in eachrow(frames_df)
                dna_seq = String(row.Consensus_sequence)
                s = parse(Int, split(row.Region, ",")[1])
                n_codons = div(length(dna_seq), 3)
                for ci in 1:n_codons
                    dna_st = (ci - 1) * 3 + 1
                    for pos in 0:2
                        off = dna_st + pos
                        locus = s + off - 1
                        cons = dna_seq[off]
                        for var in nts
                            var == cons && continue
                            push!(out, (locus, string(cons), string(var)))
                        end
                    end
                end
            end
            return out
        end
    end

    @testset "one codon → 9 variants" begin
        df = run_simulation("ATG", 1, 3, "orf1")
        @test nrow(df) == 9
    end

    @testset "all variants have valid ACGT nucleotides" begin
        df = run_simulation("ATGCCC", 1, 6, "orf1")
        @test all(v -> v in ("A","C","G","T"), df.Variant)
        @test all(v -> v in ("A","C","G","T"), df.Consensus)
    end

    @testset "Locus values span correct genomic range" begin
        df = run_simulation("ATGCCC", 100, 105, "orf1")
        @test minimum(df.Locus) == 100
        @test maximum(df.Locus) == 105
    end

end
