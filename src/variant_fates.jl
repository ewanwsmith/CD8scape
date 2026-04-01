#!/usr/bin/env julia
"""
variant_fates.jl

Traces each variant in variants.csv through the pipeline's filtering stages and
writes variant_fates.csv with a `fate` column for every (Locus, Consensus, Variant).

Fate values (in pipeline order):
  out_of_frame  - locus does not fall within any reading frame
  synonymous    - amino acid change is synonymous in all matching frames
  stop_codon    - derived sequence introduces an earlier stop codon
  non_binding   - variant is non-synonymous but both HMBR_A and HMBR_D > 2
  passed        - variant survives all filters and appears in harmonic_mean_best_ranks.csv

Called automatically at the end of `run` and `run_supertype`.

Usage:
    julia variant_fates.jl <folder_path> [--suffix <name>] [--latest|--no-latest]
"""

using DataFrames, CSV
include("path_utils.jl")

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

translate_dna(seq::String) =
    join([get(CODON_DICT, uppercase(seq[i:i+2]), "?") for i in 1:3:(length(seq)-2)], "")

"""
    classify_in_frame(consensus_seq, ref_nt, alt_nt, rel_locus)

Returns one of "synonymous", "stop_codon", or "viable" for a variant
that is confirmed to fall within a reading frame.
"""
function classify_in_frame(consensus_seq::String, ref_nt::String, alt_nt::String, rel_locus::Int)
    n = length(consensus_seq)
    ref_end = rel_locus + length(ref_nt) - 1
    ref_end > n && return "synonymous"  # out-of-bounds: treat conservatively

    prefix = rel_locus > 1 ? consensus_seq[1:(rel_locus-1)] : ""
    suffix = ref_end < n  ? consensus_seq[(ref_end+1):end]  : ""

    anc_aa = translate_dna(string(prefix, ref_nt, suffix))
    der_aa = translate_dna(string(prefix, alt_nt, suffix))

    anc_stop = findfirst(==('*'), anc_aa)
    der_stop = findfirst(==('*'), der_aa)

    # Stop codon introduced or shifted earlier in the derived sequence
    if der_stop !== nothing && (anc_stop === nothing || der_stop < anc_stop)
        return "stop_codon"
    end

    anc_aa == der_aa ? "synonymous" : "viable"
end

# Fate priority (higher index = more favourable / further through pipeline)
const FATE_RANK = Dict("out_of_frame"=>0, "stop_codon"=>1, "synonymous"=>2,
                       "non_binding"=>3, "passed"=>4)

function main()
    folder_path = ARGS[1]
    suffix = ""
    latest = true
    i = 2
    while i <= length(ARGS)
        a = ARGS[i]
        if a == "--suffix" && i + 1 <= length(ARGS)
            i += 1; suffix = ARGS[i]
        elseif a == "--latest";    latest = true
        elseif a == "--no-latest"; latest = false
        end
        i += 1
    end

    variants_path = resolve_read(joinpath(folder_path, "variants.csv"); suffix=suffix, latest=latest)
    frames_path   = resolve_read(joinpath(folder_path, "frames.csv");   suffix="",     latest=latest)

    variants = CSV.read(variants_path, DataFrame)
    frames   = CSV.read(frames_path,   DataFrame)

    frames[!, :Start] = [parse(Int, split(split(r, ";")[1],   ",")[1]) for r in frames.Region]
    frames[!, :End]   = [parse(Int, split(split(r, ";")[end], ",")[2]) for r in frames.Region]

    # Loci that survived the binding filter
    passed_loci = Set{Int}()
    hmbr_path = resolve_read(joinpath(folder_path, "harmonic_mean_best_ranks.csv"); suffix=suffix, latest=latest)
    if isfile(hmbr_path)
        hmbr = CSV.read(hmbr_path, DataFrame)
        union!(passed_loci, hmbr.Locus)
    else
        println("Warning: harmonic_mean_best_ranks.csv not found; variants that reached the binding filter will be marked non_binding.")
    end

    fates = String[]
    for var_row in eachrow(variants)
        locus = var_row.Locus
        ref   = String(var_row.Consensus)
        alt   = String(var_row.Variant)

        matching = filter(f -> f.Start <= locus <= f.End, frames)

        if isempty(matching)
            push!(fates, "out_of_frame")
            continue
        end

        # Classify across all matching frames; keep the most favourable fate
        best = "out_of_frame"
        for f in eachrow(matching)
            rel = locus - f.Start + 1
            frame_class = classify_in_frame(String(f.Consensus_sequence), ref, alt, rel)
            candidate = if frame_class == "viable"
                locus in passed_loci ? "passed" : "non_binding"
            else
                frame_class  # "synonymous" or "stop_codon"
            end
            if get(FATE_RANK, candidate, 0) > get(FATE_RANK, best, 0)
                best = candidate
            end
        end

        push!(fates, best)
    end

    result = copy(variants)
    result[!, :fate] = fates

    out_path = resolve_write(joinpath(folder_path, "variant_fates.csv"); suffix=suffix)
    CSV.write(out_path, result)
    println("variant_fates.csv written to $out_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) == 1 && ARGS[1] == "--help"
        println("Usage: julia variant_fates.jl <folder_path> [--suffix <name>] [--latest|--no-latest]")
        exit(0)
    elseif length(ARGS) < 1
        println("Error: No data folder provided."); exit(1)
    end
    main()
end
