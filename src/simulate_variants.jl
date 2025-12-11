#!/usr/bin/env julia
"""
simulate_variants.jl

Generates a synthetic variants.csv from frames.csv such that
for each codon in each frame, we enumerate all possible
single-nucleotide substitutions at each codon position (0,1,2).

Output variants.csv columns:
- Locus: genomic nucleotide position mutated (any of the 3 codon bases)
- Consensus: consensus nucleotide (A/C/G/T) at that genomic position
- Variant: variant nucleotide (A/C/G/T, excluding the consensus)

Usage:
    julia simulate_variants.jl <folder_path>

Requires:
    frames.csv in <folder_path> with columns: Region, Consensus_sequence, Description
"""

using CSV
using DataFrames
using Random

const AMINO_ACIDS = collect("ARNDCEQGHILKMFPSTWYV")  # 20 canonical AAs

# Global codon dictionary for translation
const CODON_DICT = Dict(
    "ATA"=>"I","ATC"=>"I","ATT"=>"I","ATG"=>"M",
    "ACA"=>"T","ACC"=>"T","ACG"=>"T","ACT"=>"T",
    "AAC"=>"N","AAT"=>"N","AAA"=>"K","AAG"=>"K",
    "AGC"=>"S","AGT"=>"S","AGA"=>"R","AGG"=>"R",
    "CTA"=>"L","CTC"=>"L","CTG"=>"L","CTT"=>"L",
    "CCA"=>"P","CCC"=>"P","CCG"=>"P","CCT"=>"P",
    "CAC"=>"H","CAT"=>"H","CAA"=>"Q","CAG"=>"Q",
    "CGA"=>"R","CGC"=>"R","CGG"=>"R","CGT"=>"R",
    "GTA"=>"V","GTC"=>"V","GTG"=>"V","GTT"=>"V",
    "GCA"=>"A","GCC"=>"A","GCG"=>"A","GCT"=>"A",
    "GAC"=>"D","GAT"=>"D","GAA"=>"E","GAG"=>"E",
    "GGA"=>"G","GGC"=>"G","GGG"=>"G","GGT"=>"G",
    "TCA"=>"S","TCC"=>"S","TCG"=>"S","TCT"=>"S",
    "TTC"=>"F","TTT"=>"F","TTA"=>"L","TTG"=>"L",
    "TAC"=>"Y","TAT"=>"Y","TAA"=>"*","TAG"=>"*",
    "TGC"=>"C","TGT"=>"C","TGA"=>"*","TGG"=>"W"
)

function parse_region_bounds(region::AbstractString)
    # Remove quotes if present
    region = replace(region, '"' => "")
    parts = split(region, ";")
    first_seg = split(parts[1], ",")
    start = parse(Int, first_seg[1])
    # Use last segment for end
    last_seg = split(parts[end], ",")
    end_ = parse(Int, last_seg[2])
    return start, end_
end

function translate_dna_to_protein(dna_sequence::String)::String
    prot = IOBuffer()
    up = uppercase(dna_sequence)
    for i in 1:3:length(up)-2
        aa = get(CODON_DICT, up[i:i+2], "?")
        print(prot, aa)
    end
    return String(take!(prot))
end

function main()
    if length(ARGS) < 1
        println("Usage: julia simulate_variants.jl <folder_path>")
        exit(1)
    end
    folder = ARGS[1]
    frames_path = joinpath(folder, "frames.csv")
    if !isfile(frames_path)
        println("Error: frames.csv not found in $folder")
        exit(1)
    end

    frames = CSV.read(frames_path, DataFrame)
    # Prepare variants rows
    out = DataFrame(Locus=Int[], Consensus=String[], Variant=String[])

    nucleotides = ['A', 'C', 'G', 'T']

    for row in eachrow(frames)
        region = replace(row.Region, '"' => "")
        region_bounds = split(region, ";")
        start_nt = parse(Int, split(region_bounds[1], ",")[1])
        end_nt = parse(Int, split(region_bounds[end], ",")[2])
        dna = String(row.Consensus_sequence)
        # For each codon, mutate each position (0,1,2) to capture all possible AA changes via single-nucleotide variants
        n_codons = div(length(dna), 3)
        for codon_idx in 1:n_codons
            dna_start = (codon_idx - 1) * 3 + 1
            for pos in 0:2
                offset = dna_start + pos
                locus = start_nt + offset - 1
                consensus_nt = dna[offset]
                for var_nt in nucleotides
                    if var_nt == consensus_nt
                        continue
                    end
                    push!(out, (locus, string(consensus_nt), string(var_nt)))
                end
            end
        end
    end

    # Parse additional arguments
    sample_mode = "t"
    prop = 0.1
    n = 100
    seed = 1320
    i = 2
    while i <= length(ARGS)
        arg = ARGS[i]
        if arg == "--sample"
            i += 1
            sample_mode = ARGS[i]
        elseif arg == "--prop"
            i += 1
            prop = parse(Float64, ARGS[i])
        elseif arg == "--n"
            i += 1
            n = parse(Int, ARGS[i])
        elseif arg == "--seed"
            i += 1
            seed = parse(Int, ARGS[i])
        end
        i += 1
    end
    Random.seed!(seed)

    variants_path = joinpath(folder, "variants.csv")
    println("Simulated variant rows: ", nrow(out))
    # Sampling logic
    if sample_mode == "p"
        if prop <= 0 || prop > 1.0
            println("Warning: --prop value ", prop, " is outside [0,1]. Defaulting to 1.0 (100%).")
            prop = 1.0
        end
        out = out[shuffle(1:nrow(out))[1:ceil(Int, prop * nrow(out))], :]
        println("Sampled proportion: ", prop)
        println("(output: ", nrow(out), " variants)")
    elseif sample_mode == "n"
        n = min(n, nrow(out))
        out = out[shuffle(1:nrow(out))[1:n], :]
        println("Sampled n: ", n)
        println("(output: ", nrow(out), " variants)")
    elseif sample_mode == "t"
        println("Total variants: ", nrow(out))
    else
        println("Unknown sample mode: ", sample_mode, ". Defaulting to total.")
    end
    # Always write the current 'out' DataFrame (sampled or total)
    CSV.write(variants_path, out)
    println("Simulated variants written to: $variants_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
