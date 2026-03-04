#!/usr/bin/env julia
"""
parse_aa_variants.jl

Parses an amino-acid-level variant file (.aa) and writes variants.csv.

Input file format (.aa):
    <orf_name> <aa_position>
    <ancestral_aa> <derived_aa>

    (blank lines between records are ignored)

Example:
    Orf3 23
    K M

    Orf1 45
    A T

The orf_name must match a Description value in frames.csv.
The aa_position is 1-based within the translated protein.
Ancestral and derived amino acids are single-letter codes.

Usage:
    julia parse_aa_variants.jl <folder_path> [--suffix <name>] [--latest|--no-latest]
"""

using DataFrames
using CSV
include("path_utils.jl")

const CODON_DICT = Dict(
    "ATA" => "I", "ATC" => "I", "ATT" => "I", "ATG" => "M",
    "ACA" => "T", "ACC" => "T", "ACG" => "T", "ACT" => "T",
    "AAC" => "N", "AAT" => "N", "AAA" => "K", "AAG" => "K",
    "AGC" => "S", "AGT" => "S", "AGA" => "R", "AGG" => "R",
    "CTA" => "L", "CTC" => "L", "CTG" => "L", "CTT" => "L",
    "CCA" => "P", "CCC" => "P", "CCG" => "P", "CCT" => "P",
    "CAC" => "H", "CAT" => "H", "CAA" => "Q", "CAG" => "Q",
    "CGA" => "R", "CGC" => "R", "CGG" => "R", "CGT" => "R",
    "GTA" => "V", "GTC" => "V", "GTG" => "V", "GTT" => "V",
    "GCA" => "A", "GCC" => "A", "GCG" => "A", "GCT" => "A",
    "GAC" => "D", "GAT" => "D", "GAA" => "E", "GAG" => "E",
    "GGA" => "G", "GGC" => "G", "GGG" => "G", "GGT" => "G",
    "TCA" => "S", "TCC" => "S", "TCG" => "S", "TCT" => "S",
    "TTC" => "F", "TTT" => "F", "TTA" => "L", "TTG" => "L",
    "TAC" => "Y", "TAT" => "Y", "TAA" => "*", "TAG" => "*",
    "TGC" => "C", "TGT" => "C", "TGA" => "*", "TGG" => "W"
)

"""Build a reverse codon dictionary: single-letter AA -> first alphabetical codon."""
function build_aa_to_codon()::Dict{String,String}
    d = Dict{String, String}()
    for (codon, aa) in sort(collect(CODON_DICT); by=first)
        if aa != "*" && !haskey(d, aa)
            d[aa] = codon
        end
    end
    return d
end

const AA_TO_CODON = build_aa_to_codon()

"""
    parse_aa_file(filepath) -> Vector{NamedTuple}

Parse the .aa file and return a list of records, each with fields:
  orf_name, aa_pos, ancestral_aa, derived_aa
"""
function parse_aa_file(filepath::String)
    records = NamedTuple{(:orf_name, :aa_pos, :ancestral_aa, :derived_aa),
                         Tuple{String,Int,String,String}}[]
    lines = filter(!isempty, [strip(l) for l in readlines(filepath)])

    i = 1
    while i <= length(lines)
        line1 = lines[i]
        if i + 1 > length(lines)
            @warn "Incomplete record at line $i of $filepath (missing aa line). Skipping."
            break
        end
        line2 = lines[i + 1]

        tokens1 = split(line1)
        tokens2 = split(line2)

        if length(tokens1) < 2
            @warn "Line $i ('$line1') does not have <orf_name> <aa_position>. Skipping."
            i += 1
            continue
        end
        if length(tokens2) < 2
            @warn "Line $(i+1) ('$line2') does not have <ancestral_aa> <derived_aa>. Skipping."
            i += 2
            continue
        end

        orf_name = join(tokens1[1:end-1], " ")
        aa_pos   = tryparse(Int, tokens1[end])
        anc_aa   = uppercase(String(tokens2[1]))
        der_aa   = uppercase(String(tokens2[2]))

        if isnothing(aa_pos) || aa_pos < 1
            @warn "Invalid aa_position '$(tokens1[end])' at line $i. Skipping."
            i += 2
            continue
        end

        if length(anc_aa) != 1 || length(der_aa) != 1
            @warn "Ancestral ('$anc_aa') and derived ('$der_aa') must be single-letter amino acid codes at line $(i+1). Skipping."
            i += 2
            continue
        end

        push!(records, (orf_name=orf_name, aa_pos=aa_pos, ancestral_aa=anc_aa, derived_aa=der_aa))
        i += 2
    end

    return records
end

"""
    frame_start(region_str) -> Int

Extract the start nucleotide position from a Region string (e.g. "77,496" or "77,496;606,980").
"""
function frame_start(region_str::String)::Int
    first_segment = split(region_str, ';')[1]
    return parse(Int, split(first_segment, ',')[1])
end

function main()
    if length(ARGS) < 1
        println("Usage: julia parse_aa_variants.jl <folder_path> [--suffix <name>] [--latest|--no-latest]")
        exit(1)
    end

    folder_path = ARGS[1]
    suffix = ""
    latest = true
    i = 2
    while i <= length(ARGS)
        a = ARGS[i]
        if a == "--suffix"
            if i + 1 <= length(ARGS) && !startswith(ARGS[i+1], "--")
                i += 1
                suffix = ARGS[i]
            end
        elseif a == "--latest"
            latest = true
        elseif a == "--no-latest"
            latest = false
        end
        i += 1
    end

    # Find the .aa file
    aa_files = filter(f -> endswith(f, ".aa"), readdir(folder_path; join=true))
    if isempty(aa_files)
        println("Error: No .aa file found in $folder_path.")
        exit(1)
    end
    if length(aa_files) > 1
        println("Warning: Multiple .aa files found. Using: $(aa_files[1])")
    end
    aa_filepath = aa_files[1]
    println("Using AA variant file: $aa_filepath")

    # Read frames.csv
    frames_path = resolve_read(joinpath(folder_path, "frames.csv"); suffix=suffix, latest=latest)
    frames = CSV.read(frames_path, DataFrame)
    println("Using frames file: $frames_path")

    # Parse the .aa file
    records = parse_aa_file(aa_filepath)
    if isempty(records)
        println("Error: No valid records found in $aa_filepath.")
        exit(1)
    end

    # Build variants, forcing canonical codons regardless of consensus
    out = DataFrame(Locus=Int[], Consensus=String[], Variant=String[])

    for rec in records
        # Find the frame whose Description matches orf_name
        matches = findall(==(rec.orf_name), frames.Description)
        if isempty(matches)
            @warn "No frame found with Description '$(rec.orf_name)'. Skipping variant $(rec.orf_name) $(rec.aa_pos) $(rec.ancestral_aa)→$(rec.derived_aa)."
            continue
        end
        if length(matches) > 1
            @warn "Multiple frames match Description '$(rec.orf_name)'. Using first match."
        end
        frame_idx = matches[1]
        frame_row = frames[frame_idx, :]

        region_str = String(frame_row.Region)
        if occursin(';', region_str)
            @warn "Frame '$(rec.orf_name)' has a multi-segment region ('$region_str'). Locus is computed from the first segment start; results may be inaccurate for AAs in later segments."
        end

        start_pos = frame_start(region_str)
        cons_seq  = String(frame_row.Consensus_sequence)

        # Relative nucleotide position within Consensus_sequence (1-based)
        rel_nt = (rec.aa_pos - 1) * 3 + 1

        if rel_nt + 2 > length(cons_seq)
            @warn "AA position $(rec.aa_pos) is out of bounds for frame '$(rec.orf_name)' (sequence length $(length(cons_seq)) nt = $(length(cons_seq) ÷ 3) aa). Skipping."
            continue
        end

        if !haskey(AA_TO_CODON, rec.ancestral_aa)
            @warn "No codon found for ancestral amino acid '$(rec.ancestral_aa)'. Skipping variant $(rec.orf_name) $(rec.aa_pos) $(rec.ancestral_aa)→$(rec.derived_aa)."
            continue
        end
        if !haskey(AA_TO_CODON, rec.derived_aa)
            @warn "No codon found for derived amino acid '$(rec.derived_aa)'. Skipping variant $(rec.orf_name) $(rec.aa_pos) $(rec.ancestral_aa)→$(rec.derived_aa)."
            continue
        end

        # Always use canonical codons (forced), regardless of what is in the consensus.
        # This allows substitutions to be specified even when sequence data is missing or
        # the consensus does not match the specified ancestral amino acid.
        ancestral_codon = AA_TO_CODON[rec.ancestral_aa]
        derived_codon   = AA_TO_CODON[rec.derived_aa]

        # Report and warn if the consensus codon differs from the canonical ancestral
        actual_codon      = uppercase(cons_seq[rel_nt:(rel_nt + 2)])
        translated_actual = get(CODON_DICT, actual_codon, "?")
        if translated_actual != rec.ancestral_aa
            @warn "Consensus mismatch at $(rec.orf_name) AA $(rec.aa_pos): consensus is '$translated_actual' (codon '$actual_codon'), expected '$(rec.ancestral_aa)'. Forcing ancestral codon to '$ancestral_codon' and overwriting frames."
        end

        if ancestral_codon == derived_codon
            @warn "Canonical codons for '$(rec.ancestral_aa)' and '$(rec.derived_aa)' are identical. Skipping $(rec.orf_name) $(rec.aa_pos) $(rec.ancestral_aa)→$(rec.derived_aa)."
            continue
        end

        # Force the canonical ancestral codon into the frame's Consensus_sequence so that
        # generate_peptides.jl produces the correct ancestral peptide even on mismatch.
        frames[frame_idx, :Consensus_sequence] = cons_seq[1:rel_nt-1] * ancestral_codon * cons_seq[rel_nt+3:end]

        # Absolute nucleotide locus (start of the codon in genome coordinates)
        locus = start_pos + rel_nt - 1

        push!(out, (locus, ancestral_codon, derived_codon))
        println("  $(rec.orf_name) AA $(rec.aa_pos): $(rec.ancestral_aa)→$(rec.derived_aa) | codon $ancestral_codon→$derived_codon | locus $locus")
    end

    if nrow(out) == 0
        println("Error: No variants could be parsed from $aa_filepath.")
        exit(1)
    end

    # Write back the (possibly modified) frames so generate_peptides.jl sees forced codons
    CSV.write(frames_path, frames)
    println("Frames written to: $frames_path")

    output_path = resolve_write(joinpath(folder_path, "variants.csv"); suffix=suffix)
    CSV.write(output_path, out)
    println("Variants written to: $output_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
