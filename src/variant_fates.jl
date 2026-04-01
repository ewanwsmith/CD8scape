#!/usr/bin/env julia
"""
variant_fates.jl

Traces each variant in variants.csv through the pipeline's filtering stages and
writes variant_fates.csv in wide format with columns:

  Frame          - protein / region label  (matches harmonic_mean_best_ranks.csv)
  Locus          - genomic position
  Mutation       - amino-acid change label (e.g. A17T, K50*, del100atg)
  frame_filter   - "In frame" | "Out of frame"
  peptide_filter - "Non-synonymous" | "Synonymous" | "Stop codon" | missing
  binding_filter - "Passed" | "Non-binding" | missing

One row is emitted per variant (or per (Frame, Mutation) entry in
harmonic_mean_best_ranks.csv when a locus appears in multiple ORFs).
Adjacent stage columns give the source→target pairs needed for a Sankey diagram.

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

    mutation = if is_indel
        length(ref_nt) > length(alt_nt) ?
            "del$(genomic_locus)$(lowercase(ref_nt[(length(alt_nt)+1):end]))" :
            "ins$(genomic_locus)$(lowercase(alt_nt[(length(ref_nt)+1):end]))"
    else
        anc_c = aa_pos <= length(anc_aa) ? string(anc_aa[aa_pos]) : "?"
        der_c = aa_pos <= length(der_aa) ? string(der_aa[aa_pos]) : "?"
        "$(anc_c)$(aa_pos)$(der_c)"
    end

    anc_stop = findfirst(==('*'), anc_aa)
    der_stop = findfirst(==('*'), der_aa)
    if der_stop !== nothing && (anc_stop === nothing || der_stop < anc_stop)
        return ("stop_codon", mutation)
    end

    anc_aa == der_aa ? ("synonymous", mutation) : ("viable", mutation)
end

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

    # Build locus → [(Frame, Mutation)] from HMBR (one entry per passing ORF)
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

    col_frame    = Union{String,Missing}[]
    col_locus    = Int[]
    col_mutation = Union{String,Missing}[]
    col_stage1   = String[]                      # frame_filter
    col_stage2   = Union{String,Missing}[]       # peptide_filter
    col_stage3   = Union{String,Missing}[]       # binding_filter

    function push_row!(frame, locus, mutation, s1, s2, s3)
        push!(col_frame, frame); push!(col_locus, locus); push!(col_mutation, mutation)
        push!(col_stage1, s1);   push!(col_stage2, s2);   push!(col_stage3, s3)
    end

    for var_row in eachrow(variants)
        locus = var_row.Locus
        ref   = String(var_row.Consensus)
        alt   = String(var_row.Variant)

        matching = filter(f -> f.Start <= locus <= f.End, frames)

        if isempty(matching)
            push_row!(missing, locus, missing, "Out of frame", missing, missing)
            continue
        end

        # Passed binding filter — one row per (Frame, Mutation) in HMBR
        if haskey(hmbr_by_locus, locus)
            for e in hmbr_by_locus[locus]
                push_row!(e.Frame, locus, e.Mutation, "In frame", "Non-synonymous", "Passed")
            end
            continue
        end

        # Classify from translation; take most favourable fate across matching frames
        best_class    = "synonymous"
        best_mutation = ""
        best_frame    = String(matching[1, :Description])

        for f in eachrow(matching)
            rel = locus - f.Start + 1
            fc, mut = classify_and_label(String(f.Consensus_sequence), ref, alt, rel, locus)
            if get(FATE_RANK, fc, 0) > get(FATE_RANK, best_class, 0)
                best_class = fc; best_mutation = mut; best_frame = String(f.Description)
            elseif best_mutation == ""
                best_mutation = mut; best_frame = String(f.Description)
            end
        end

        s2, s3 = if best_class == "viable"
            "Non-synonymous", "Non-binding"
        elseif best_class == "stop_codon"
            "Stop codon", missing
        else
            "Synonymous", missing
        end

        push_row!(best_frame, locus, best_mutation, "In frame", s2, s3)
    end

    result = DataFrame(
        Frame          = col_frame,
        Locus          = col_locus,
        Mutation       = col_mutation,
        frame_filter   = col_stage1,
        peptide_filter = col_stage2,
        binding_filter = col_stage3,
    )
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
