#!/usr/bin/env julia

using DataFrames
using CSV

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
    join_data(folder_path::String) :: DataFrame

Reads `frames.csv` and `variants.csv` from the specified folder, extracts start and
end positions from the 'Region' column, and performs a cross join followed by filtering
to keep only rows where the variant locus is within the frame.

# Arguments
- folder_path::String: The directory containing frames.csv and variants.csv

# Returns
- A DataFrame of joined data.
"""
function join_data(folder_path::String)::DataFrame
    frames_path = joinpath(folder_path, "frames.csv")
    variants_path = joinpath(folder_path, "variants.csv")
    
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

Calculates the relative locus and extracts the pulled base from the consensus sequence.
Removed filtering on consensus match so that all rows are retained here.

# Arguments
- df::DataFrame: The input DataFrame after join_data.

# Returns
- The same DataFrame with added Relative_Locus, Pulled_Base, and Matches_Consensus columns.
"""
function check_locus(df::DataFrame)::DataFrame
    df[!, :Relative_Locus] = df.Locus .- df.Start
    
    df[!, :Pulled_Base] = [
        (1 ≤ rl ≤ length(seq) ? string(seq[rl]) : missing)
        for (rl, seq) in zip(df.Relative_Locus, df.Consensus_sequence)
    ]
    
    df[!, :Matches_Consensus] = coalesce.(df.Pulled_Base .== df.Consensus, false)
    return df
end

"""
    edit_consensus_sequence(df::DataFrame) :: DataFrame

Edits the consensus sequence at the relative locus to produce a variant sequence.

# Arguments
- df::DataFrame: Input DataFrame after check_locus.

# Returns
- A DataFrame with a new 'Variant_sequence' column.
"""
function edit_consensus_sequence(df::DataFrame)::DataFrame
    df[!, :Variant_sequence] = Vector{String}(undef, nrow(df))
    
    for row in eachrow(df)
        if 1 ≤ row.Relative_Locus ≤ length(row.Consensus_sequence)
            sequence_chars = collect(row.Consensus_sequence)
            sequence_chars[row.Relative_Locus] = row.Variant[1]
            row.Variant_sequence = join(sequence_chars)
        else
            row.Variant_sequence = row.Consensus_sequence
        end
    end
    
    return df
end

"""
    translate_sequences(df::DataFrame) :: DataFrame

Translates both the consensus and variant DNA sequences into amino acid sequences.

# Arguments
- df::DataFrame: Input DataFrame after edit_consensus_sequence.

# Returns
- A DataFrame with 'Consensus_AA_sequence' and 'Variant_AA_sequence' columns added.
"""
function translate_sequences(df::DataFrame)::DataFrame
    translate_dna_to_protein(dna_sequence::String)::String = 
        join([get(CODON_DICT, uppercase(dna_sequence[i:i+2]), "?") 
              for i in 1:3:length(dna_sequence)-2], "")
    
    df[!, :Consensus_AA_sequence] = [translate_dna_to_protein(seq) for seq in df.Consensus_sequence]
    df[!, :Variant_AA_sequence] = [translate_dna_to_protein(seq) for seq in df.Variant_sequence]
    return df
end

"""
    generate_peptides(sequence::String, locus::Int, substr_lengths::Vector{Int}) :: Vector{String}

Generates all possible peptide substrings of specified lengths that contain the specified locus.

# Arguments
- sequence::String: The amino acid sequence.
- locus::Int: Position of interest (1-based index).
- substr_lengths::Vector{Int}: List of peptide lengths to generate.

# Returns
- A vector of peptide strings.
"""
function generate_peptides(sequence::String, locus::Int, substr_lengths::Vector{Int})::Vector{String}
    peptides = String[]
    seq_length = length(sequence)
    
    for substr_len in substr_lengths
        start_min = max(1, locus - substr_len + 1)
        start_max = min(seq_length - substr_len + 1, locus)
        
        for start_pos in start_min:start_max
            end_pos = start_pos + substr_len - 1
            if end_pos ≤ seq_length && start_pos ≤ locus ≤ end_pos
                push!(peptides, sequence[start_pos:end_pos])
            end
        end
    end
    
    return peptides
end

"""
    add_peptides_columns!(df::DataFrame, 
                          relative_locus_col::Symbol, 
                          consensus_col::Symbol, 
                          variant_col::Symbol, 
                          substr_lengths::Vector{Int}) :: DataFrame

Generates peptides from both consensus and variant amino acid sequences and flattens
the data into a new DataFrame with corresponding annotations.

# Arguments
- df::DataFrame: Input DataFrame after translate_sequences.
- relative_locus_col::Symbol: Column symbol for Relative_Locus.
- consensus_col::Symbol: Column symbol for Consensus amino acid sequences.
- variant_col::Symbol: Column symbol for Variant amino acid sequences.
- substr_lengths::Vector{Int}: Peptide lengths to generate.

# Returns
- A flattened DataFrame with columns for consensus peptides, variant peptides, and labels.
"""
function add_peptides_columns!(
    df::DataFrame, 
    relative_locus_col::Symbol, 
    consensus_col::Symbol, 
    variant_col::Symbol, 
    substr_lengths::Vector{Int}
)::DataFrame
    # Calculate AA_Locus from Relative_Locus
    df[!, :AA_Locus] = ceil.(Int, df[!, relative_locus_col] / 3) .+ 1
    
    # Initialize the output DataFrame
    flattened_peptides_df = DataFrame(
        Locus = Int[],
        Relative_Locus = Int[],
        AA_Locus = Int[],
        Consensus_Peptide = String[],
        Variant_Peptide = String[],
        Peptide_label = String[]
    )
    
    # Iterate through each row of the DataFrame
    for row in eachrow(df)
        # Generate the consensus and variant peptides
        consensus_peptides = generate_peptides(row[consensus_col], row.AA_Locus, substr_lengths)
        variant_peptides = generate_peptides(row[variant_col], row.AA_Locus, substr_lengths)
        
        if length(consensus_peptides) == length(variant_peptides)
            counter = 1
            for (cons_pep, var_pep) in zip(consensus_peptides, variant_peptides)
                # Replace spaces with underscores in Description
                description_with_underscores = replace(row.Description, " " => "_")
                # Include the ORF information (with underscores) in the label
                peptide_label = "$(row.Consensus)$(row.Locus)$(row.Variant)_$(description_with_underscores)_$counter"
                push!(flattened_peptides_df, (
                    row.Locus, 
                    row.Relative_Locus, 
                    row.AA_Locus,
                    cons_pep, 
                    var_pep, 
                    peptide_label
                ))
                counter += 1
            end
        else
            @warn "Mismatched peptide lists in row with Locus $(row.Locus). Skipping this row."
        end
    end
    
    return flattened_peptides_df
end

"""
    separate_peptides(df::DataFrame) :: DataFrame

Splits the flattened DataFrame into two rows per original row: one for the consensus
peptide and one for the variant peptide.

Only non-synonymous variants (where consensus and variant peptides differ) are retained.

# Arguments
- df::DataFrame: Input DataFrame after add_peptides_columns!.

# Returns
- A DataFrame with separate rows for consensus and variant peptides.
"""
function separate_peptides(df::DataFrame)::DataFrame
    # Retain only rows where the consensus peptide and variant peptide differ.
    filtered_df = filter(row -> row.Consensus_Peptide != row.Variant_Peptide, df)
    
    transformed_df = DataFrame(
        Locus = Int[],
        Peptide = String[],
        Peptide_label = String[]
    )

    for row in eachrow(filtered_df)
        # Add consensus peptide row
        push!(transformed_df, (
            row.Locus,
            row.Consensus_Peptide,
            "$(row.Peptide_label)_C"
        ))

        # Add variant peptide row
        push!(transformed_df, (
            row.Locus,
            row.Variant_Peptide,
            "$(row.Peptide_label)_V"
        ))
    end

    return transformed_df
end

# ------------------------------------------------------------------------------
# Background peptide generation functions
# ------------------------------------------------------------------------------
"""
    generate_background_peptides(aa_sequence::String, region::String, description::String)

Generates background peptide sequences for the consensus amino acid sequence based on region and description.

# Arguments
- aa_sequence::String: The amino acid sequence.
- region::String: The nucleotide region in "start,end" format or "start1,end1;start2,end2".
- description::String: The description to include in the label.

# Returns
- A tuple of peptide sequences and corresponding labels.
"""
function generate_background_peptides(aa_sequence::String, region::String, description::String)
    peptide_lengths = [8, 9, 10, 11]
    peptides = String[]
    labels = String[]

    # Handle multiple region entries (e.g., "13468,13502;13600,13650")
    region_parts = split(region, ";")  # ["13468,13502", "13600,13650"]

    start_nt, _ = parse.(Int, split(region_parts[1], ","))  # First pair
    _, end_nt = parse.(Int, split(region_parts[end], ","))  # Last pair

    # Replace spaces with underscores in Description
    description_clean = replace(description, " " => "_")

    for len in peptide_lengths
        for i in 1:(length(aa_sequence) - len + 1)
            peptide = aa_sequence[i:(i+len-1)]
            
            # Compute the nucleotide loci for this peptide
            peptide_start_nt = start_nt + (i - 1) * 3
            peptide_end_nt = start_nt + (i + len - 2) * 3

            push!(peptides, peptide)
            
            # Create label, incorporating ORF information (Description) with underscores
            label = "$(peptide_start_nt)-$(peptide_end_nt)_$(description_clean)_A"
            push!(labels, label)
        end
    end

    return peptides, labels
end

"""
    create_background_peptide_dataframe(frames::DataFrame)

Creates a DataFrame for background peptides using the provided frames DataFrame.

# Arguments
- frames::DataFrame: The input frames DataFrame containing amino acid sequences.

# Returns
- A DataFrame containing background peptides and their labels.
"""
function create_background_peptide_dataframe(frames::DataFrame)
    peptide_list = String[]
    label_list = String[]

    for row in eachrow(frames)
        if !ismissing(row.Consensus_AA_sequence) && row.Consensus_AA_sequence != "missing"
            # Pass Description along with other data to generate_background_peptides
            peptides, labels = generate_background_peptides(row.Consensus_AA_sequence, String(row.Region), String(row.Description))
            append!(peptide_list, peptides)
            append!(label_list, labels)
        end
    end

    return DataFrame(Peptide = peptide_list, Peptide_label = label_list)
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
    println("Run `julia script.jl --help` for usage.")
    exit(1)
end

# The first argument should be the data_folder
data_folder = ARGS[1]

# Run the processing pipeline
joined = join_data(data_folder)
checked = check_locus(joined)
edited = edit_consensus_sequence(checked)
translated = translate_sequences(edited)

# Define locus-based peptide lengths to generate
substr_lengths = [8, 9, 10, 11]

flattened_peptides_df = add_peptides_columns!(
    translated, 
    :Relative_Locus, 
    :Consensus_AA_sequence, 
    :Variant_AA_sequence, 
    substr_lengths
)

transformed_df = separate_peptides(flattened_peptides_df)

# generate background peptides
background_peptides_df = create_background_peptide_dataframe(translated)

# Create a new "Locus" column for background_peptides_df with all values set to "A"
background_peptides_df.Locus .= "A"

# Stack the dataframes vertically
peptide_df = vcat(transformed_df, background_peptides_df)

# Save to CSV
csv_file_path = joinpath(data_folder, "peptides_labels.csv")
CSV.write(csv_file_path, peptide_df)
println("peptides_labels.csv file has been written to: $csv_file_path")

# Save peptides to .pep file
write_peptides_file_no_headers(peptide_df, data_folder)