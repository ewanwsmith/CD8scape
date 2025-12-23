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
    julia simulate_variants.jl <folder_path> [--n <count>] [--p <proportion>] [--seed <int>] [--suffix <name>] [--latest]

Sampling options (standalone, no --sample needed):
- --n: sample an absolute count of variants (e.g., --n 500)
- --p: sample a proportion of variants in [0,1] (e.g., --p 0.25)
    If both --n and --p are provided, --n takes precedence.

Requires:
    frames.csv in <folder_path> with columns: Region, Consensus_sequence, Description
"""

using CSV
using DataFrames
using Random
include("path_utils.jl")

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

 

function main()
    if length(ARGS) < 1
        println("Usage: julia simulate_variants.jl <folder_path> [--n <count>] [--p <proportion>] [--seed <int>] [--suffix <name>] [--latest]")
        exit(1)
    end
    folder = ARGS[1]
    # Defaults for simulation (suffix token without leading underscore)
    suffix = "simulated"
    latest = true

    # Parse additional arguments
    # New behavior: --n and --p are sufficient; no --sample mode required
    # Defaults when flags are present without values:
    #   --n       -> 1000
    #   --p/--prop-> 0.1
    have_n = false
    have_p = false
    prop = 1.0
    n = 0
    seed = 1320
    i = 2
    while i <= length(ARGS)
        arg = ARGS[i]
        if arg == "--prop" || arg == "--p"
            # Support optional value; default to 0.1 if omitted
            # Example: `--p` is equivalent to `--p 0.1`
            if i + 1 <= length(ARGS) && !startswith(ARGS[i+1], "--")
                i += 1
                prop = parse(Float64, ARGS[i])
            else
                prop = 0.1
            end
            have_p = true
        elseif arg == "--n"
            # Support optional value; default to 1000 if omitted
            # Example: `--n` is equivalent to `--n 1000`
            if i + 1 <= length(ARGS) && !startswith(ARGS[i+1], "--")
                i += 1
                n = parse(Int, ARGS[i])
            else
                n = 1000
            end
            have_n = true
        elseif arg == "--seed"
            i += 1
            seed = parse(Int, ARGS[i])
        elseif arg == "--suffix"
            # Optional: override default output suffix token (no leading underscore)
            if i + 1 <= length(ARGS) && !startswith(ARGS[i+1], "--")
                i += 1
                suffix = ARGS[i]
                if startswith(suffix, "_")
                    suffix = suffix[2:end]
                end
            else
                suffix = "simulated"
            end
        elseif arg == "--latest"
            latest = true
        end
        i += 1
    end
    Random.seed!(seed)
    # Resolve frames.csv for reading (ignore suffix for input; allow latest fallback)
    frames_path = resolve_read(joinpath(folder, "frames.csv"); suffix="", latest=latest)
    if !isfile(frames_path)
        println("Error: frames.csv not found in $folder")
        exit(1)
    end

    # Read frames and build full variant enumeration
    frames = CSV.read(frames_path, DataFrame)
    out = DataFrame(Locus=Int[], Consensus=String[], Variant=String[])
    nucleotides = ['A', 'C', 'G', 'T']
    for row in eachrow(frames)
        region = replace(row.Region, '"' => "")
        region_bounds = split(region, ";")
        start_nt = parse(Int, split(region_bounds[1], ",")[1])
        end_nt = parse(Int, split(region_bounds[end], ",")[2])
        dna = String(row.Consensus_sequence)
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

    variants_path = resolve_write(joinpath(folder, "variants.csv"); suffix=suffix)
    println("Simulated variant rows: ", nrow(out))
    # Sampling logic: --n and --p are sufficient
    if have_n && have_p
        println("Note: both --n and --p provided; using --n and ignoring --p.")
    end
    if have_n
        if n <= 0
            println("Warning: --n should be > 0. Outputting total variants.")
        else
            n = min(n, nrow(out))
            out = out[shuffle(1:nrow(out))[1:n], :]
            println("Sampled n: ", n)
            println("(output: ", nrow(out), " variants)")
        end
    elseif have_p
        if prop <= 0 || prop > 1.0
            println("Warning: --p/--prop value ", prop, " is outside (0,1]. Defaulting to 1.0 (100%).")
            prop = 1.0
        end
        k = ceil(Int, prop * nrow(out))
        out = out[shuffle(1:nrow(out))[1:k], :]
        println("Sampled proportion: ", prop)
        println("(output: ", nrow(out), " variants)")
    else
        println("No sampling options provided; outputting total variants: ", nrow(out))
    end
    # Always write the current 'out' DataFrame (sampled or total)
    CSV.write(variants_path, out)
    println("Simulated variants written to: $variants_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
