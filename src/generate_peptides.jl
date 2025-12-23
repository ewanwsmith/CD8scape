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
function check_locus(df::DataFrame)::DataFrame
    # Fix: make Relative_Locus 1-based (first nucleotide in frame is 1)
    df[!, :Relative_Locus] = df.Locus .- df.Start .+ 1
    
    df[!, :Pulled_Base] = [
        (1 ≤ rl ≤ length(seq) ? string(seq[rl]) : missing)
        for (rl, seq) in zip(df.Relative_Locus, df.Consensus_sequence)
    ]
    
    df[!, :Matches_Consensus] = coalesce.(df.Pulled_Base .== df.Consensus, false)
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
        variant_nt = row.Variant

        if !ismissing(rl) && 1 ≤ rl ≤ length(consensus)
            seq = collect(consensus)
            seq[rl] = variant_nt[1]
            row.Derived_sequence = join(seq)
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
        Ancestral_Peptide = String[],
        Derived_Peptide = String[],
        Peptide_label = String[]
    )

    for row in eachrow(df)
        cps = generate_peptides(row[cs], row.AA_Locus, lens)
        vps = generate_peptides(row[vs], row.AA_Locus, lens)

        if length(cps) == length(vps)
            consensus_aa = (1 ≤ row.AA_Locus ≤ length(row.Ancestral_AA_sequence)) ? row.Ancestral_AA_sequence[row.AA_Locus] : "?"
            variant_aa = (1 ≤ row.AA_Locus ≤ length(row.Derived_AA_sequence)) ? row.Derived_AA_sequence[row.AA_Locus] : "?"

            change_label = "$(consensus_aa)$(row.AA_Locus)$(variant_aa)"
            base = "$(change_label)_$(replace(String(row.Description), " " => "_"))"

            for (i, (c, v)) in enumerate(zip(cps, vps))
                label = "$(base)_$(i)"
                push!(out, (row.Locus, row.Relative_Locus, row.AA_Locus, c, v, label))
            end
        else
            @warn "Mismatched peptide counts at Locus $(row.Locus). Skipping."
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
    # Retain only rows where the consensus peptide and variant peptide differ.
    initial_rows = nrow(df)
    
    # Identify loci dropped for synonymity and stop codons separately
    dropped_syn_df = filter(row -> row.Ancestral_Peptide == row.Derived_Peptide, df)
    dropped_stop_df = filter(row -> occursin('*', row.Ancestral_Peptide) || occursin('*', row.Derived_Peptide), df)
    removed_syn_loci = length(unique(dropped_syn_df.Locus))
    removed_stop_loci = length(unique(dropped_stop_df.Locus))
    removed_syn_rows = nrow(dropped_syn_df)
    removed_stop_rows = nrow(dropped_stop_df)

    filtered_df = filter(row -> row.Ancestral_Peptide != row.Derived_Peptide &&
                                !occursin('*', row.Ancestral_Peptide) &&
                                !occursin('*', row.Derived_Peptide), df)

    println("Removed $removed_syn_rows peptides from $removed_syn_loci loci due to synonymity.")
    println("Removed $removed_stop_rows peptides from $removed_stop_loci loci due to stop codons.")
    
    transformed_df = DataFrame(
        Locus = Int[],
        Peptide = String[],
        Peptide_label = String[]
    )

    for row in eachrow(filtered_df)
        # Add ancestral peptide row
        push!(transformed_df, (
            row.Locus,
            row.Ancestral_Peptide,
            "$(row.Peptide_label)_A"
        ))

        # Add derived peptide row
        push!(transformed_df, (
            row.Locus,
            row.Derived_Peptide,
            "$(row.Peptide_label)_D"
        ))
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