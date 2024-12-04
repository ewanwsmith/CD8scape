# Include environment configurations
include("./env.jl")

using DataFrames
using CSV

# ------------------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------------------

# Codon to amino acid mapping
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
    join_data(folder_path::String) -> DataFrame

Reads `frames.csv` and `trajectories.csv` from the specified folder, processes the
regions to extract start and end positions, performs a cross join, and filters
rows based on locus conditions.

# Arguments
- `folder_path::String`: Path to the folder containing the CSV files.

# Returns
- `DataFrame`: The joined and filtered DataFrame.
"""
function join_data(folder_path::String)::DataFrame
    # Construct file paths
    frames_path = joinpath(folder_path, "frames.csv")
    trajectories_path = joinpath(folder_path, "trajectories.csv")
    
    # Read CSV files
    frames = CSV.read(frames_path, DataFrame)
    trajectories = CSV.read(trajectories_path, DataFrame)
    
    # Extract Start and End from the Region column
    frames[!, :Start] = [
        parse(Int, split(split(region, ";")[1], ",")[1]) for region in frames.Region
    ]
    frames[!, :End] = [
        parse(Int, split(split(region, ";")[end], ",")[2]) for region in frames.Region
    ]
    
    # Perform cross join and filter based on locus conditions
    result = crossjoin(trajectories, frames)
    filter!(row -> row.Start ≤ row.Locus ≤ row.End, result)
    
    return result
end

"""
    check_locus(df::DataFrame) -> DataFrame

Generates `Relative_Locus`, extracts `Pulled_Base` from `Consensus_sequence`,
checks for consensus matches, filters matching rows, and removes unnecessary columns.

# Arguments
- `df::DataFrame`: The input DataFrame.

# Returns
- `DataFrame`: The filtered DataFrame with consensus matches.
"""
function check_locus(df::DataFrame)::DataFrame
    # Calculate Relative_Locus
    df[!, :Relative_Locus] = df.Locus .- df.Start
    
    # Extract Pulled_Base from Consensus_sequence
    df[!, :Pulled_Base] = [
        (1 ≤ rl ≤ length(seq) ? string(seq[rl]) : missing)
        for (rl, seq) in zip(df.Relative_Locus, df.Consensus_sequence)
    ]
    
    # Check for consensus matches
    df[!, :Matches_Consensus] = df.Pulled_Base .== df.Consensus
    
    # Replace `missing` with `false` before filtering
    df[!, :Matches_Consensus] = coalesce.(df[!, :Matches_Consensus], false)
    
    # Filter rows where Matches_Consensus is true
    df = filter(:Matches_Consensus => identity, df)
    
    # Remove Pulled_Base and Matches_Consensus columns
    select!(df, Not([:Pulled_Base, :Matches_Consensus]))
    
    return df
end

"""
    edit_consensus_sequence(df::DataFrame) -> DataFrame

Modifies the `Consensus_sequence` at the `Relative_Locus` to match the `Variant`,
creating a new `Variant_sequence` column.

# Arguments
- `df::DataFrame`: The input DataFrame.

# Returns
- `DataFrame`: The DataFrame with the `Variant_sequence` column added.
"""
function edit_consensus_sequence(df::DataFrame)::DataFrame
    # Initialize Variant_sequence column
    df[!, :Variant_sequence] = Vector{String}(undef, nrow(df))
    
    # Iterate over each row to modify the sequence
    for row in eachrow(df)
        if 1 ≤ row.Relative_Locus ≤ length(row.Consensus_sequence)
            # Convert to character array for modification
            sequence_chars = collect(row.Consensus_sequence)
            # Replace character at Relative_Locus with Variant
            sequence_chars[row.Relative_Locus] = row.Variant[1]
            # Update Variant_sequence
            row.Variant_sequence = join(sequence_chars)
        else
            # If out of bounds, retain the original Consensus_sequence
            row.Variant_sequence = row.Consensus_sequence
        end
    end
    
    return df
end

"""
    translate_sequences(df::DataFrame) -> DataFrame

Translates DNA sequences in `Consensus_sequence` and `Variant_sequence` to
amino acid sequences using the provided codon table.

# Arguments
- `df::DataFrame`: The input DataFrame.

# Returns
- `DataFrame`: The DataFrame with translated amino acid sequences.
"""
function translate_sequences(df::DataFrame)::DataFrame
    # Function to translate a DNA sequence to a protein sequence
    translate_dna_to_protein(dna_sequence::String)::String = 
        join([get(CODON_DICT, uppercase(dna_sequence[i:i+2]), "?") 
              for i in 1:3:length(dna_sequence)-2], "")
    
    # Translate Consensus and Variant sequences
    df[!, :Consensus_AA_sequence] = [translate_dna_to_protein(seq) for seq in df.Consensus_sequence]
    df[!, :Variant_AA_sequence] = [translate_dna_to_protein(seq) for seq in df.Variant_sequence]
    
    return df
end

"""
    generate_peptides(sequence::String, locus::Int, substr_lengths::Vector{Int}) -> Vector{String}

Generates peptides of specified lengths that include the locus position.

# Arguments
- `sequence::String`: The amino acid sequence.
- `locus::Int`: The position of interest within the sequence.
- `substr_lengths::Vector{Int}`: Desired lengths of peptides.

# Returns
- `Vector{String}`: A vector of generated peptide sequences.
"""
function generate_peptides(sequence::String, locus::Int, substr_lengths::Vector{Int})::Vector{String}
    peptides = String[]
    seq_length = length(sequence)
    
    for substr_len in substr_lengths
        # Determine valid start positions
        start_min = max(1, locus - substr_len + 1)
        start_max = min(seq_length - substr_len + 1, locus)
        
        for start_pos in start_min:start_max
            end_pos = start_pos + substr_len - 1
            # Ensure the locus is within the peptide
            if end_pos <= seq_length && start_pos ≤ locus ≤ end_pos
                push!(peptides, sequence[start_pos:end_pos])
            end
        end
    end
    
    return peptides
end

"""
    add_peptides_columns!(df::DataFrame, relative_locus_col::Symbol, 
                         consensus_col::Symbol, variant_col::Symbol, 
                         substr_lengths::Vector{Int}) -> DataFrame

Generates peptides from consensus and variant amino acid sequences and adds them
to a new flattened DataFrame with relevant annotations.

# Arguments
- `df::DataFrame`: The input DataFrame.
- `relative_locus_col::Symbol`: Column name for Relative_Locus.
- `consensus_col::Symbol`: Column name for Consensus amino acid sequences.
- `variant_col::Symbol`: Column name for Variant amino acid sequences.
- `substr_lengths::Vector{Int}`: Desired lengths of peptides.

# Returns
- `DataFrame`: A flattened DataFrame containing peptide pairs and annotations.
"""
function add_peptides_columns!(
    df::DataFrame, 
    relative_locus_col::Symbol, 
    consensus_col::Symbol, 
    variant_col::Symbol, 
    substr_lengths::Vector{Int}
)::DataFrame
    # Calculate AA_Locus (1-based indexing)
    df[!, :AA_Locus] = ceil.(Int, df[!, relative_locus_col] / 3) .+ 1
    
    # Initialize the flattened peptides DataFrame
    flattened_peptides_df = DataFrame(
        Locus = Int[],
        Relative_Locus = Int[],
        AA_Locus = Int[],
        Consensus_Peptide = String[],
        Variant_Peptide = String[],
        Peptide_label = String[]
    )
    
    # Iterate through each row to generate peptides
    for row in eachrow(df)
        consensus_peptides = generate_peptides(row[consensus_col], row.AA_Locus, substr_lengths)
        variant_peptides = generate_peptides(row[variant_col], row.AA_Locus, substr_lengths)
        
        # Ensure both peptide lists are of the same length
        if length(consensus_peptides) == length(variant_peptides)
            for (cons_pep, var_pep) in zip(consensus_peptides, variant_peptides)
                # Calculate site position within the peptide
                site_position = row.AA_Locus - (length(row[consensus_col]) - length(cons_pep))
                
                # Create peptide label
                peptide_label = "$(row.Consensus)$(row.Locus)$(row.Variant) "
                
                # Append to the flattened DataFrame
                push!(flattened_peptides_df, (
                    row.Locus, 
                    row.Relative_Locus, 
                    row.AA_Locus,
                    cons_pep, 
                    var_pep, 
                    peptide_label
                ))
            end
        else
            @warn "Mismatched peptide lists in row with Locus $(row.Locus). Skipping this row."
        end
    end
    
    return flattened_peptides_df
end

# ------------------------------------------------------------------------------
# Data Processing Pipeline
# ------------------------------------------------------------------------------

# Define the folder path containing the CSV files
data_folder = "/Users/e.smith.5/Documents/PhD/CD8scape/data/RSV_example"

# Join Data
joined = join_data(data_folder)

# Check Locus
checked = check_locus(joined)

# Edit Consensus Sequence
edited = edit_consensus_sequence(checked)

# Translate Sequences
translated = translate_sequences(edited)

# Define desired peptide lengths
substr_lengths = [8, 9, 10, 11]

# Generate Flattened Peptides DataFrame
flattened_peptides_df = add_peptides_columns!(
    translated, 
    :Relative_Locus, 
    :Consensus_AA_sequence, 
    :Variant_AA_sequence, 
    substr_lengths
)

# Display the resulting DataFrame
display(flattened_peptides_df)