#!/usr/bin/env julia
"""
generate_peptides.jl

Generates ancestral and derived peptides for each locus from frames and variants data.

Usage:
    julia generate_peptides.jl <folder_path> [--suffix <name>] [--latest|--no-latest]

Arguments:
    <folder_path>   Path to the folder containing frames.csv and variants.csv.
"""

using DataFrames
using CSV
include("path_utils.jl")

# ------------------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------------------
"""
A dictionary mapping codons to their corresponding amino acids.
"""
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

# ------------------------------------------------------------------------------
# Function Definitions
# ------------------------------------------------------------------------------
"""
print_help()

Prints out usage and help information for this script.
"""
function print_help()
    println("Usage: julia generate_peptides.jl <data_folder>")
    println("")
    println("This script processes peptide data from frames.csv and variants.csv, ")
    println("found in the specified data_folder directory, and generates two output files:")
    println("  1. peptides_labels.csv: A CSV file containing labeled peptide sequences.")
    println("  2. Peptides.pep: A .pep file containing peptide sequences without headers.")
    println("")
    println("Arguments:")
    println("  <data_folder>: A directory containing frames.csv and variants.csv.")
    println("")
    println("For help, run: generate_peptides.jl --help")
end

"""
    generate_peptides(sequence::String, aa_locus::Int, lengths::Vector{Int}) :: Vector{String}

Generates peptide sequences of specified lengths from a given amino acid sequence.

# Arguments
- sequence::String: The amino acid sequence to generate peptides from.
- aa_locus::Int: The starting amino acid locus for generating peptides.
- lengths::Vector{Int}: A vector of peptide lengths to generate.

# Returns
- A vector of generated peptide sequences.
"""
function generate_peptides(sequence::String, aa_locus::Int, lengths::Vector{Int})::Vector{String}
    peptides = String[]
    for len in lengths
        for i in 1:(length(sequence) - len + 1)
            start_pos = i
            end_pos = i + len - 1
            if start_pos ≤ aa_locus ≤ end_pos
                push!(peptides, sequence[start_pos:end_pos])
            end
        end
    end
    return peptides
end

"""
    join_data(folder_path::String) :: DataFrame

Reads `frames.csv` and `variants.csv` from the specified folder, extracts start and
end positions from the 'Region' column, and performs a cross join followed by filtering
to keep only rows where the variant locus is within the frame.

# Arguments
- folder_path::String: The directory containing frames.csv and variants.csv

# Returns
- A DataFrame of joined data.
"""
function join_data(folder_path::String; suffix::AbstractString="", latest::Bool=true)::DataFrame
    # Prefer suffixed inputs if they exist; otherwise, fall back to latest discovery
    base_frames   = joinpath(folder_path, "frames.csv")
    base_variants = joinpath(folder_path, "variants.csv")

    candidate_frames = with_suffix(base_frames, suffix)
    frames_path = isfile(candidate_frames) ? candidate_frames : discover_path(base_frames; latest=latest)

    candidate_variants = with_suffix(base_variants, suffix)
    variants_path = isfile(candidate_variants) ? candidate_variants : discover_path(base_variants; latest=latest)
    
    frames = CSV.read(frames_path, DataFrame)
    variants = CSV.read(variants_path, DataFrame)
    
    frames[!, :Start] = [
        parse(Int, split(split(region, ";")[1], ",")[1]) for region in frames.Region
    ]
    frames[!, :End] = [
        parse(Int, split(split(region, ";")[end], ",")[2]) for region in frames.Region
    ]
    
    # Cross join variants and frames, then keep rows where the locus falls within the frame.
    result = crossjoin(variants, frames)
    filter!(row -> row.Start ≤ row.Locus ≤ row.End, result)
    
    return result
end

"""
    check_locus(df::DataFrame) :: DataFrame

Calculates the relative locus and extracts the pulled base from the ancestral (consensus) sequence.
Removed filtering on ancestral match so that all rows are retained here.

# Arguments
- df::DataFrame: The input DataFrame after join_data.

# Returns
- The same DataFrame with added Relative_Locus, Pulled_Base, and Matches_Consensus columns.
"""
function normalize_variant(ref::AbstractString, alt::AbstractString)
    r = String(ref)
    a = String(alt)
    # Trim common suffix while keeping at least 1 base in each (VCF-style)
    while length(r) > 1 && length(a) > 1 && last(r) == last(a)
        r = r[1:end-1]
        a = a[1:end-1]
    end
    prefix_trim = 0
    while length(r) > 1 && length(a) > 1 && first(r) == first(a)
        r = r[2:end]
        a = a[2:end]
        prefix_trim += 1
    end
    return r, a, prefix_trim
end

function find_ref_match(consensus::AbstractString, ref::AbstractString, rl::Int; window::Int=10, max_mismatch::Int=0)
    seq = String(consensus)
    r = String(ref)
    rlen = length(r)
    if rlen == 0
        return rl, true
    end
    min_pos = max(1, rl - window)
    max_pos = min(length(seq) - rlen + 1, rl + window)
    best = nothing
    for pos in min_pos:max_pos
        window_ref = seq[pos:(pos + rlen - 1)]
        mismatches = 0
        if window_ref == r
            mismatches = 0
        else
            for i in 1:rlen
                if window_ref[i] != r[i]
                    mismatches += 1
                    if mismatches > max_mismatch
                        break
                    end
                end
            end
        end
        if mismatches <= max_mismatch
            dist = abs(pos - rl)
            best = best === nothing ? (dist, mismatches, pos) :
                (dist < best[1] || (dist == best[1] && mismatches < best[2]) ? (dist, mismatches, pos) : best)
        end
    end
    return best === nothing ? (rl, false) : (best[3], true)
end

function check_locus(df::DataFrame)::DataFrame
    # Normalize variant representation (trim common prefix/suffix) and optionally realign
    norm_ref = String[]
    norm_alt = String[]
    adjusted_locus = Int[]
    matches = Bool[]

    for row in eachrow(df)
        ref = String(row.Consensus)
        alt = String(row.Variant)
        r, a, prefix_trim = normalize_variant(ref, alt)
        adj = row.Locus + prefix_trim
        # Compute Relative_Locus from adjusted position
        rl = adj - row.Start + 1
        # If ref does not match at the adjusted position, try to find nearby match
        ref_match = false
        if 1 ≤ rl && rl + length(r) - 1 ≤ length(row.Consensus_sequence)
            ref_match = row.Consensus_sequence[rl:(rl + length(r) - 1)] == r
        end
        if !ref_match
            # allow a small number of mismatches to tolerate minor ref discrepancies
            new_rl, ok = find_ref_match(row.Consensus_sequence, r, rl; window=10, max_mismatch=1)
            if ok
                adj = row.Start + new_rl - 1
                rl = new_rl
                if 1 ≤ rl && rl + length(r) - 1 ≤ length(row.Consensus_sequence)
                    r = row.Consensus_sequence[rl:(rl + length(r) - 1)]
                end
                ref_match = true
            else
                # normalize against consensus at the adjusted locus
                ref_len = length(r)
                if 1 ≤ rl && rl + ref_len - 1 ≤ length(row.Consensus_sequence)
                    r = row.Consensus_sequence[rl:(rl + ref_len - 1)]
                    ref_match = true
                    @warn "Reference allele mismatch at Locus $(row.Locus). Normalizing REF to consensus sequence."
                end
            end
        end
        push!(norm_ref, r)
        push!(norm_alt, a)
        push!(adjusted_locus, adj)
        push!(matches, ref_match)
    end

    df[!, :Normalized_Consensus] = norm_ref
    df[!, :Normalized_Variant] = norm_alt
    df[!, :Adjusted_Locus] = adjusted_locus
    df[!, :Matches_Consensus] = matches

    # Fix: make Relative_Locus 1-based (first nucleotide in frame is 1)
    df[!, :Relative_Locus] = df.Adjusted_Locus .- df.Start .+ 1

    df[!, :Pulled_Base] = [
        begin
            if ismissing(rl)
                missing
            else
                ref = String(cons)
                ref_len = length(ref)
                if 1 ≤ rl && rl + ref_len - 1 ≤ length(seq)
                    string(seq[rl:(rl + ref_len - 1)])
                else
                    missing
                end
            end
        end
        for (rl, seq, cons) in zip(df.Relative_Locus, df.Consensus_sequence, df.Normalized_Consensus)
    ]

    return df
end

"""
    edit_ancestral_sequence(df::DataFrame) :: DataFrame

Edits the ancestral (consensus) sequence at the relative locus to produce a derived sequence.

# Arguments
- df::DataFrame: Input DataFrame after check_locus.

# Returns
- A DataFrame with a new 'Derived_sequence' column.
"""
function edit_ancestral_sequence(df::DataFrame)::DataFrame
    df[!, :Derived_sequence] = similar(df.Consensus_sequence)
    for row in eachrow(df)
        rl = row.Relative_Locus
        consensus = row.Consensus_sequence
        variant_nt = String(row.Normalized_Variant)
        ref_nt = String(row.Normalized_Consensus)

        if !ismissing(rl) && 1 ≤ rl ≤ length(consensus)
            ref_len = length(ref_nt)
            seq_len = length(consensus)
            if rl + ref_len - 1 ≤ seq_len && consensus[rl:(rl + ref_len - 1)] == ref_nt
                prefix = rl > 1 ? consensus[1:(rl - 1)] : ""
                suffix_start = rl + ref_len
                suffix = suffix_start ≤ seq_len ? consensus[suffix_start:end] : ""
                row.Derived_sequence = string(prefix, variant_nt, suffix)
            else
                @warn "Reference allele mismatch or out of bounds at Locus $(row.Locus). Keeping consensus sequence."
                row.Derived_sequence = consensus
            end
        else
            row.Derived_sequence = consensus
        end
    end
    return df
end

"""
    translate_sequences(df::DataFrame) :: DataFrame

Translates both the ancestral and derived DNA sequences into amino acid sequences.

# Arguments
- df::DataFrame: Input DataFrame after edit_consensus_sequence.

# Returns
- A DataFrame with 'Ancestral_AA_sequence' and 'Derived_AA_sequence' columns added.
"""
function translate_sequences(df::DataFrame)::DataFrame
    translate_dna_to_protein(dna_sequence::String)::String = 
        join([get(CODON_DICT, uppercase(dna_sequence[i:i+2]), "?") 
              for i in 1:3:length(dna_sequence)-2], "")
    
    df[!, :Ancestral_AA_sequence] = [translate_dna_to_protein(seq) for seq in df.Consensus_sequence]
    df[!, :Derived_AA_sequence] = [translate_dna_to_protein(seq) for seq in df.Derived_sequence]
    return df
end

"""
    add_peptides_columns!(df::DataFrame, 
                          relative_locus_col::Symbol, 
                          consensus_col::Symbol, 
                          variant_col::Symbol, 
                          substr_lengths::Vector{Int}) :: DataFrame

Generates peptides from both ancestral and derived amino acid sequences and flattens
the data into a new DataFrame with corresponding annotations.

# Arguments
- df::DataFrame: Input DataFrame after translate_sequences.
- relative_locus_col::Symbol: Column symbol for Relative_Locus.
- consensus_col::Symbol: Column symbol for Ancestral amino acid sequences.
- variant_col::Symbol: Column symbol for Derived amino acid sequences.
- substr_lengths::Vector{Int}: Peptide lengths to generate.

# Returns
- A flattened DataFrame with columns for ancestral peptides, derived peptides, and labels.
"""
function add_peptides_columns!(df::DataFrame, rl::Symbol, cs::Symbol, vs::Symbol, lens::Vector{Int})::DataFrame
    # Fix: AA_Locus should be 1-based, so use ceil instead of floor
    df[!, :AA_Locus] = ceil.(Int, df[!, rl] ./ 3)

    out = DataFrame(
        Locus = Int[],
        Relative_Locus = Int[],
        AA_Locus = Int[],
        Ancestral_Peptide = Union{Missing,String}[],
        Derived_Peptide = Union{Missing,String}[],
        Peptide_label = String[]
    )

    for row in eachrow(df)
        cps = generate_peptides(row[cs], row.AA_Locus, lens)
        vps = generate_peptides(row[vs], row.AA_Locus, lens)

        ref_nt = String(row.Normalized_Consensus)
        alt_nt = String(row.Normalized_Variant)
        is_indel = length(ref_nt) != length(alt_nt)

        if length(cps) == length(vps) && !is_indel
            consensus_aa = (1 ≤ row.AA_Locus ≤ length(row.Ancestral_AA_sequence)) ? row.Ancestral_AA_sequence[row.AA_Locus] : "?"
            variant_aa = (1 ≤ row.AA_Locus ≤ length(row.Derived_AA_sequence)) ? row.Derived_AA_sequence[row.AA_Locus] : "?"

            # Include all variants (synonymous filtering happens in separate_peptides)
            change_label = "$(consensus_aa)$(row.AA_Locus)$(variant_aa)"
            base = "$(change_label)_$(replace(String(row.Description), " " => "_"))"

            for (i, (c, v)) in enumerate(zip(cps, vps))
                label = "$(base)_$(i)"
                push!(out, (row.Locus, row.Relative_Locus, row.AA_Locus, c, v, label))
            end
        elseif is_indel
            # Label indels as del/ins + Locus + only what was deleted or inserted
            if length(ref_nt) > length(alt_nt)
                # Deletion: show only the deleted nucleotides
                deleted = lowercase(ref_nt[length(alt_nt)+1:end])
                change_label = "del$(row.Locus)$(deleted)"
            else
                # Insertion: show only the inserted nucleotides
                inserted = lowercase(alt_nt[length(ref_nt)+1:end])
                change_label = "ins$(row.Locus)$(inserted)"
            end
            base = "$(change_label)_$(replace(String(row.Description), " " => "_"))"

            maxlen = max(length(cps), length(vps))
            for i in 1:maxlen
                label = "$(base)_$(i)"
                c = i <= length(cps) ? cps[i] : missing
                v = i <= length(vps) ? vps[i] : missing
                push!(out, (row.Locus, row.Relative_Locus, row.AA_Locus, c, v, label))
            end
        else
            consensus_aa = (1 ≤ row.AA_Locus ≤ length(row.Ancestral_AA_sequence)) ? row.Ancestral_AA_sequence[row.AA_Locus] : "?"
            variant_aa = (1 ≤ row.AA_Locus ≤ length(row.Derived_AA_sequence)) ? row.Derived_AA_sequence[row.AA_Locus] : "?"

            change_label = "$(consensus_aa)$(row.AA_Locus)$(variant_aa)"
            base = "$(change_label)_$(replace(String(row.Description), " " => "_"))"

            maxlen = max(length(cps), length(vps))
            for i in 1:maxlen
                label = "$(base)_$(i)"
                c = i <= length(cps) ? cps[i] : missing
                v = i <= length(vps) ? vps[i] : missing
                push!(out, (row.Locus, row.Relative_Locus, row.AA_Locus, c, v, label))
            end
            @warn "Mismatched peptide counts at Locus $(row.Locus). Emitting all available peptides (A=$(length(cps)), D=$(length(vps)))."
        end
    end

    return out
end

"""
    separate_peptides(df::DataFrame) :: DataFrame

Splits the flattened DataFrame into two rows per original row: one for the ancestral
peptide and one for the derived peptide.

Only non-synonymous variants (where ancestral and derived peptides differ) are retained.

# Arguments
- df::DataFrame: Input DataFrame after add_peptides_columns!.

# Returns
- A DataFrame with separate rows for ancestral and derived peptides.
"""
function separate_peptides(df::DataFrame)::DataFrame
    # Retain peptides where at least one non-missing peptide is valid and non-synonymous
    initial_rows = nrow(df)

    dropped_syn_df = filter(row -> !ismissing(row.Ancestral_Peptide) && !ismissing(row.Derived_Peptide) && row.Ancestral_Peptide == row.Derived_Peptide, df)
    dropped_stop_df = filter(row -> (!ismissing(row.Ancestral_Peptide) && occursin('*', row.Ancestral_Peptide)) ||
                                   (!ismissing(row.Derived_Peptide) && occursin('*', row.Derived_Peptide)), df)
    removed_syn_loci = length(unique(dropped_syn_df.Locus))
    removed_stop_loci = length(unique(dropped_stop_df.Locus))
    removed_syn_rows = nrow(dropped_syn_df)
    removed_stop_rows = nrow(dropped_stop_df)

    filtered_df = filter(row -> begin
            c = row.Ancestral_Peptide
            v = row.Derived_Peptide
            c_missing = ismissing(c)
            v_missing = ismissing(v)
            c_stop = !c_missing && occursin('*', c)
            v_stop = !v_missing && occursin('*', v)
            synonymous = !c_missing && !v_missing && c == v
            keep_c = !c_missing && !c_stop && !synonymous
            keep_v = !v_missing && !v_stop && !synonymous
            keep_c && keep_v
        end, df)

    println("Removed $removed_syn_rows peptides from $removed_syn_loci loci due to synonymity.")
    println("Removed $removed_stop_rows peptides from $removed_stop_loci loci due to stop codons.")
    
    transformed_df = DataFrame(
        Locus = Int[],
        Peptide = String[],
        Peptide_label = String[]
    )

    for row in eachrow(filtered_df)
        c = row.Ancestral_Peptide
        v = row.Derived_Peptide
        c_missing = ismissing(c)
        v_missing = ismissing(v)
        c_stop = !c_missing && occursin('*', c)
        v_stop = !v_missing && occursin('*', v)
        synonymous = !c_missing && !v_missing && c == v

        if !c_missing && !c_stop && !synonymous
            push!(transformed_df, (
                row.Locus,
                c,
                "$(row.Peptide_label)_A"
            ))
        end

        if !v_missing && !v_stop && !synonymous
            push!(transformed_df, (
                row.Locus,
                v,
                "$(row.Peptide_label)_D"
            ))
        end
    end

    return transformed_df
end

function deduplicate_peptides(df::DataFrame)::DataFrame
    # Preserve distinct AA changes per locus: include full Peptide_label in dedup key
    # so that ancestral peptides shared across multiple variants are not collapsed.
    df[!, :State] = last.(split.(df.Peptide_label, "_"))
    dedup_df = unique(df, [:Locus, :Peptide, :State, :Peptide_label])
    select!(dedup_df, Not(:State))  # remove helper column
    return dedup_df
end

"""
    write_peptides_file_no_headers(df::DataFrame, folder_path::String)

Writes the peptide sequences (only) to a .pep file in the specified folder.

# Arguments
- df::DataFrame: Input DataFrame from separate_peptides().
- folder_path::String: Output directory path.

# Output
- A file named 'Peptides.pep' containing peptide sequences.
"""
function write_peptides_file_no_headers(df::DataFrame, folder_path::String)
    file_path = joinpath(folder_path, "Peptides.pep")
    
    open(file_path, "w") do io
        for row in eachrow(df)
            println(io, row.Peptide)
        end
    end

    println("Peptides.pep file has been written to: $file_path")
end

# ------------------------------------------------------------------------------
# Main Script
# ------------------------------------------------------------------------------
if length(ARGS) == 1 && ARGS[1] == "--help"
    print_help()
    exit(0)
elseif length(ARGS) < 1
    println("Error: No data folder provided.")
    println("Run `julia generate_peptides.jl --help` for usage.")
    exit(1)
end

function main()
    # The first argument should be the data_folder
    data_folder = ARGS[1]
    # parse optional flags
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

    # Run the processing pipeline
    joined = join_data(data_folder; suffix=suffix, latest=latest)
    checked = check_locus(joined)
    edited = edit_ancestral_sequence(checked)
    translated = translate_sequences(edited)

    # Define locus-based peptide lengths to generate
    substr_lengths = [8, 9, 10, 11]

    flattened_peptides_df = add_peptides_columns!(
        translated, 
        :Relative_Locus, 
        :Ancestral_AA_sequence, 
        :Derived_AA_sequence, 
        substr_lengths
    )

    transformed_df = separate_peptides(flattened_peptides_df)
    transformed_df = deduplicate_peptides(transformed_df)

    peptide_df = transformed_df

    # Save to CSV
    csv_file_path = resolve_write(joinpath(data_folder, "peptides_labels.csv"); suffix=suffix)
    CSV.write(csv_file_path, peptide_df)
    println("peptides_labels.csv file has been written to: $csv_file_path")

    # Save peptides to .pep file
    write_peptides_file_no_headers(peptide_df, data_folder)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end