#!/usr/bin/env julia
"""
variant_fates.jl

Traces each variant in variants.csv through the pipeline's filtering stages and
writes variant_fates.csv in long format with columns:

  Frame     - protein / region label  (matches harmonic_mean_best_ranks.csv)
  Locus     - genomic position
  Mutation  - amino-acid change label (e.g. A17T, K50*, del100atg)
  source    - upstream pipeline node
  target    - downstream pipeline node

Each variant produces one row per stage it reaches.  Aggregating by
(source, target) and counting gives the link values for a Sankey diagram.

Stage nodes:
  Input          → Out of frame | In frame
  In frame       → Synonymous   | Stop codon | Non-synonymous
  Non-synonymous → Non-binding  | Passed

For variants that passed the binding filter, Frame and Mutation are taken
directly from harmonic_mean_best_ranks.csv (one set of journey rows is
emitted per (Frame, Mutation) entry at that Locus).  For variants filtered
earlier, Frame is the first matching reading frame and Mutation is derived
from translation.

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
    classify_and_label(consensus_seq, ref_nt, alt_nt, rel_locus, genomic_locus)

Returns (fate_class, mutation_label) where fate_class is one of
"synonymous", "stop_codon", or "viable".
Indel labels use genomic_locus to match the convention in harmonic_mean_best_ranks.csv.
"""
function classify_and_label(consensus_seq::String, ref_nt::String, alt_nt::String,
                             rel_locus::Int, genomic_locus::Int)
    n = length(consensus_seq)
    ref_end = rel_locus + length(ref_nt) - 1
    ref_end > n && return ("synonymous", "?")

    prefix = rel_locus > 1 ? consensus_seq[1:(rel_locus-1)] : ""
    suffix = ref_end < n  ? consensus_seq[(ref_end+1):end]  : ""

    anc_aa = translate_dna(string(prefix, ref_nt, suffix))
    der_aa = translate_dna(string(prefix, alt_nt, suffix))

    is_indel = length(ref_nt) != length(alt_nt)
    aa_pos   = ceil(Int, rel_locus / 3)

    # Build mutation label
    mutation = if is_indel
        if length(ref_nt) > length(alt_nt)
            "del$(genomic_locus)$(lowercase(ref_nt[(length(alt_nt)+1):end]))"
        else
            "ins$(genomic_locus)$(lowercase(alt_nt[(length(ref_nt)+1):end]))"
        end
    else
        anc_c = aa_pos <= length(anc_aa) ? string(anc_aa[aa_pos]) : "?"
        der_c = aa_pos <= length(der_aa) ? string(der_aa[aa_pos]) : "?"
        "$(anc_c)$(aa_pos)$(der_c)"
    end

    # Classification: stop codon introduced / moved earlier?
    anc_stop = findfirst(==('*'), anc_aa)
    der_stop = findfirst(==('*'), der_aa)
    if der_stop !== nothing && (anc_stop === nothing || der_stop < anc_stop)
        return ("stop_codon", mutation)
    end

    anc_aa == der_aa ? ("synonymous", mutation) : ("viable", mutation)
end

# Fate priority used when a variant maps to multiple frames with different fates
const FATE_RANK = Dict("synonymous"=>0, "stop_codon"=>1, "viable"=>2)

function main()
    folder_path = ARGS[1]
    suffix = ""
    latest = true
    i = 2
    while i <= length(ARGS)
        a = ARGS[i]
        if a == "--suffix" && i + 1 <= length(ARGS); i += 1; suffix = ARGS[i]
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

    # Build locus → [(Frame, Mutation)] lookup from HMBR for passed variants
    hmbr_by_locus = Dict{Int, Vector{NamedTuple{(:Frame,:Mutation),Tuple{String,String}}}}()
    hmbr_path = resolve_read(joinpath(folder_path, "harmonic_mean_best_ranks.csv");
                             suffix=suffix, latest=latest)
    if isfile(hmbr_path)
        for row in eachrow(CSV.read(hmbr_path, DataFrame))
            entry = (Frame=String(row.Frame), Mutation=String(row.Mutation))
            push!(get!(hmbr_by_locus, row.Locus, []), entry)
        end
    else
        println("Warning: harmonic_mean_best_ranks.csv not found; " *
                "binding-filter survivors will be marked Non-binding.")
    end

    # Accumulate output columns
    col_frame    = Union{String,Missing}[]
    col_locus    = Int[]
    col_mutation = Union{String,Missing}[]
    col_source   = String[]
    col_target   = String[]

    function emit!(frame, locus, mutation, source, target)
        push!(col_frame, frame); push!(col_locus, locus)
        push!(col_mutation, mutation)
        push!(col_source, source); push!(col_target, target)
    end

    for var_row in eachrow(variants)
        locus = var_row.Locus
        ref   = String(var_row.Consensus)
        alt   = String(var_row.Variant)

        matching = filter(f -> f.Start <= locus <= f.End, frames)

        if isempty(matching)
            emit!(missing, locus, missing, "Input", "Out of frame")
            continue
        end

        # Variants in HMBR passed the binding filter; emit one journey per (Frame, Mutation)
        if haskey(hmbr_by_locus, locus)
            for e in hmbr_by_locus[locus]
                emit!(e.Frame, locus, e.Mutation, "Input",          "In frame")
                emit!(e.Frame, locus, e.Mutation, "In frame",       "Non-synonymous")
                emit!(e.Frame, locus, e.Mutation, "Non-synonymous", "Passed")
            end
            continue
        end

        # Not in HMBR — classify from translation, take most favourable frame fate
        best_class    = "synonymous"
        best_mutation = ""
        best_frame    = String(matching[1, :Description])

        for f in eachrow(matching)
            rel = locus - f.Start + 1
            fc, mut = classify_and_label(String(f.Consensus_sequence), ref, alt, rel, locus)
            rank = get(FATE_RANK, fc, 0)
            if rank > get(FATE_RANK, best_class, 0)
                best_class    = fc
                best_mutation = mut
                best_frame    = String(f.Description)
            elseif best_mutation == ""
                best_mutation = mut
                best_frame    = String(f.Description)
            end
        end

        if best_class == "viable"
            # Reached binding filter but didn't pass
            emit!(best_frame, locus, best_mutation, "Input",          "In frame")
            emit!(best_frame, locus, best_mutation, "In frame",       "Non-synonymous")
            emit!(best_frame, locus, best_mutation, "Non-synonymous", "Non-binding")
        elseif best_class == "stop_codon"
            emit!(best_frame, locus, best_mutation, "Input",    "In frame")
            emit!(best_frame, locus, best_mutation, "In frame", "Stop codon")
        else
            emit!(best_frame, locus, best_mutation, "Input",    "In frame")
            emit!(best_frame, locus, best_mutation, "In frame", "Synonymous")
        end
    end

    result = DataFrame(Frame=col_frame, Locus=col_locus, Mutation=col_mutation,
                       source=col_source, target=col_target)
    out_path = resolve_write(joinpath(folder_path, "variant_fates.csv"); suffix=suffix)
    CSV.write(out_path, result)
    println("variant_fates.csv written to $out_path ($(nrow(result)) rows)")
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
