#!/usr/bin/env julia
"""
generate_peptides.jl

Generates consensus and variant peptides for each locus from frames and variants data.

Usage:
    julia generate_peptides.jl <folder_path>

Arguments:
    <folder_path>   Path to the folder containing frames.csv and variants.csv.
"""

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
    edit_consensus_sequence(df::DataFrame) :: DataFrame

Edits the consensus sequence at the relative locus to produce a variant sequence.

# Arguments
- df::DataFrame: Input DataFrame after check_locus.

# Returns
- A DataFrame with a new 'Variant_sequence' column.
"""
function edit_consensus_sequence(df::DataFrame)::DataFrame
    df[!, :Variant_sequence] = similar(df.Consensus_sequence)
    for row in eachrow(df)
        rl = row.Relative_Locus
        consensus = row.Consensus_sequence
        consensus_nt = row.Consensus
        variant_nt = row.Variant

        if !ismissing(rl) && 1 ≤ rl ≤ length(consensus)
            # Handle different types of variants
            if length(consensus_nt) == length(variant_nt)
                # SNV: Simple substitution
                seq = collect(consensus)
                seq[rl] = variant_nt[1]
                row.Variant_sequence = join(seq)
            elseif length(consensus_nt) > length(variant_nt)
                # Deletion: Remove nucleotides
                deletion_length = length(consensus_nt) - length(variant_nt)
                seq = collect(consensus)
                # Replace the REF sequence with ALT sequence at the specified position
                seq_before = seq[1:(rl-1)]
                seq_after = seq[(rl + length(consensus_nt)):end]
                new_seq = vcat(seq_before, collect(variant_nt), seq_after)
                row.Variant_sequence = join(new_seq)
            else
                # Insertion: Add nucleotides
                seq = collect(consensus)
                # Replace the REF sequence with ALT sequence
                seq_before = seq[1:(rl-1)]
                seq_after = seq[(rl + length(consensus_nt)):end]
                new_seq = vcat(seq_before, collect(variant_nt), seq_after)
                row.Variant_sequence = join(new_seq)
            end
        else
            row.Variant_sequence = consensus
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
    
    df[!, :Consensus_AA_sequence] = [translate_dna_to_protein(String(seq)) for seq in df.Consensus_sequence]
    df[!, :Variant_AA_sequence] = [translate_dna_to_protein(String(seq)) for seq in df.Variant_sequence]
    return df
end

"""
    generate_peptides_covering_variant(consensus_seq::String, variant_seq::String, aa_locus::Int, lengths::Vector{Int}) :: Tuple{Vector{String}, Vector{String}}

Generates all possible 8-11mer peptides that cover the variant amino acid position from both
consensus and variant sequences.

# Arguments
- consensus_seq::String: The consensus amino acid sequence
- variant_seq::String: The variant amino acid sequence  
- aa_locus::Int: The amino acid position of the variant
- lengths::Vector{Int}: Peptide lengths to generate (should be [8,9,10,11])

# Returns
- A tuple of (consensus_peptides, variant_peptides) vectors
"""
function generate_peptides_covering_variant(consensus_seq::String, variant_seq::String, aa_locus::Int, lengths::Vector{Int})::Tuple{Vector{String}, Vector{String}}
    consensus_peptides = String[]
    variant_peptides = String[]
    
    for len in lengths
        # Generate all possible peptides of this length that cover the variant position
        # from the consensus sequence
        for i in 1:(length(consensus_seq) - len + 1)
            start_pos = i
            end_pos = i + len - 1
            if start_pos ≤ aa_locus ≤ end_pos
                consensus_peptide = consensus_seq[start_pos:end_pos]
                push!(consensus_peptides, consensus_peptide)
            end
        end
        
        # Generate all possible peptides of this length that cover the variant position
        # from the variant sequence
        for i in 1:(length(variant_seq) - len + 1)
            start_pos = i
            end_pos = i + len - 1
            if start_pos ≤ aa_locus ≤ end_pos
                variant_peptide = variant_seq[start_pos:end_pos]
                push!(variant_peptides, variant_peptide)
            end
        end
    end
    
    return (consensus_peptides, variant_peptides)
end

"""
    add_peptides_columns!(df::DataFrame, 
                          relative_locus_col::Symbol, 
                          consensus_col::Symbol, 
                          variant_col::Symbol, 
                          substr_lengths::Vector{Int}) :: DataFrame

Generates peptides from both consensus and variant amino acid sequences and flattens
the data into a new DataFrame with corresponding annotations, adjusting for indels.

# Arguments
- df::DataFrame: Input DataFrame after translate_sequences.
- relative_locus_col::Symbol: Column symbol for Relative_Locus.
- consensus_col::Symbol: Column symbol for Consensus amino acid sequences.
- variant_col::Symbol: Column symbol for Variant amino acid sequences.
- substr_lengths::Vector{Int}: Peptide lengths to generate.

# Returns
- A flattened DataFrame with columns for consensus peptides, variant peptides, and labels.
"""
function add_peptides_columns!(df::DataFrame, rl::Symbol, cs::Symbol, vs::Symbol, lens::Vector{Int})::DataFrame
    # Fix: AA_Locus should be 1-based, so use ceil instead of floor
    df[!, :AA_Locus] = ceil.(Int, df[!, rl] ./ 3)

    out = DataFrame(
        Locus = Int[],
        Relative_Locus = Int[],
        AA_Locus = Int[],
        Consensus_Peptide = String[],
        Variant_Peptide = String[],
        Peptide_label = String[]
    )

    for row in eachrow(df)
        # Check if this is an indel (consensus and variant nucleotides have different lengths)
        is_indel = length(row.Consensus) != length(row.Variant)
        
        if is_indel
            # For indels, we need to handle the affected region differently
            # Calculate the size difference to determine the affected amino acid range
            nt_diff = abs(length(row.Variant) - length(row.Consensus))
            aa_diff = ceil(Int, nt_diff / 3)  # Number of amino acids potentially affected
            
            # Generate consensus peptides from original sequence covering the variant locus
            cps, _ = generate_peptides_covering_variant(row[cs], row[cs], row.AA_Locus, lens)
            
            # For variant peptides, we need to cover a broader region that includes the entire affected area
            # Start from the variant locus and extend by the amino acid difference
            extended_start = max(1, row.AA_Locus - aa_diff)
            extended_end = min(length(row[vs]), row.AA_Locus + aa_diff)
            
            # Generate all 8-11mers that cover any part of the extended affected region
            vps = String[]
            for aa_pos in extended_start:extended_end
                for len in lens
                    for start_pos in max(1, aa_pos - len + 1):min(length(row[vs]) - len + 1, aa_pos)
                        if start_pos + len - 1 <= length(row[vs])
                            peptide = row[vs][start_pos:start_pos + len - 1]
                            if !(peptide in vps)  # Avoid duplicates
                                push!(vps, peptide)
                            end
                        end
                    end
                end
            end
        else
            # For SNVs, use the existing logic
            cps, vps = generate_peptides_covering_variant(row[cs], row[vs], row.AA_Locus, lens)
        end

        # Create consensus/variant peptide pairs
        if is_indel
            # For indels, create nucleotide position-based notation
            if length(row.Consensus) > length(row.Variant)
                # Deletion
                deleted_length = length(row.Consensus) - length(row.Variant)
                if deleted_length == 1
                    # Single nucleotide deletion
                    change_label = "nt$(row.Relative_Locus)del"
                else
                    # Multi-nucleotide deletion
                    start_nt = row.Relative_Locus
                    end_nt = start_nt + deleted_length - 1
                    change_label = "nt$(start_nt)_$(end_nt)del"
                end
            else
                # Insertion
                inserted_length = length(row.Variant) - length(row.Consensus)
                # For insertions, use flanking positions and show inserted sequence
                pos1 = row.Relative_Locus - 1
                pos2 = row.Relative_Locus
                # Get the inserted sequence (what's in variant but not in consensus)
                inserted_seq = row.Variant[length(row.Consensus)+1:end]
                change_label = "nt$(pos1)_$(pos2)ins$(inserted_seq)"
            end
        else
            # For SNVs, use the traditional amino acid substitution format
            consensus_aa = (1 ≤ row.AA_Locus ≤ length(row.Consensus_AA_sequence)) ? row.Consensus_AA_sequence[row.AA_Locus] : "?"
            variant_aa = (1 ≤ row.AA_Locus ≤ length(row.Variant_AA_sequence)) ? row.Variant_AA_sequence[row.AA_Locus] : "?"
            change_label = "$(consensus_aa)$(row.AA_Locus)$(variant_aa)"
        end
        base = "$(change_label)_$(replace(String(row.Description), " " => "_"))"

        # Add consensus peptides
        for (i, c) in enumerate(cps)
            label = "$(base)_C_$(i)"
            push!(out, (row.Locus, row.Relative_Locus, row.AA_Locus, c, "", label))
        end
        
        # Add variant peptides
        for (i, v) in enumerate(vps)
            label = "$(base)_V_$(i)"
            push!(out, (row.Locus, row.Relative_Locus, row.AA_Locus, "", v, label))
        end
    end

    return out
end

"""
    separate_peptides(df::DataFrame) :: DataFrame

Processes the DataFrame containing consensus and variant peptides separately,
filtering out peptides with stop codons and creating the final peptide list.

# Arguments
- df::DataFrame: Input DataFrame after add_peptides_columns!.

# Returns
- A DataFrame with separate rows for consensus and variant peptides.
"""
function separate_peptides(df::DataFrame)::DataFrame
    # Filter out rows with stop codons
    initial_rows = nrow(df)
    
    dropped_df = filter(row -> 
        (row.Consensus_Peptide != "" && occursin('*', row.Consensus_Peptide)) ||
        (row.Variant_Peptide != "" && occursin('*', row.Variant_Peptide)), df)
    removed_loci = length(unique(dropped_df.Locus))
    
    filtered_df = filter(row -> 
        (row.Consensus_Peptide == "" || !occursin('*', row.Consensus_Peptide)) &&
        (row.Variant_Peptide == "" || !occursin('*', row.Variant_Peptide)), df)
    
    final_rows = nrow(filtered_df)
    removed_rows = initial_rows - final_rows
    println("Removed $removed_rows peptides from $removed_loci loci due to stop codons.")
    
    transformed_df = DataFrame(
        Locus = Int[],
        Peptide = String[],
        Peptide_label = String[]
    )

    for row in eachrow(filtered_df)
        if row.Consensus_Peptide != ""
            # This is a consensus peptide
            push!(transformed_df, (
                row.Locus,
                row.Consensus_Peptide,
                row.Peptide_label
            ))
        elseif row.Variant_Peptide != ""
            # This is a variant peptide
            push!(transformed_df, (
                row.Locus,
                row.Variant_Peptide,
                row.Peptide_label
            ))
        end
    end

    return transformed_df
end

function deduplicate_peptides(df::DataFrame)::DataFrame
    df[!, :State] = last.(split.(df.Peptide_label, "_"))
    dedup_df = unique(df, [:Locus, :Peptide, :State])
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

# Calculate AA_Locus (1-based, from Relative_Locus)
translated[!, :AA_Locus] = ceil.(Int, translated[!, :Relative_Locus] ./ 3)


# Remove loci with no amino acid change anywhere in the sequence
syn_loci = filter(row -> begin
    # Check if this is an indel
    is_indel = length(row.Consensus) != length(row.Variant)
    
    if is_indel
        # For indels, check if any amino acid changes in the entire sequence
        # This handles frameshifts that affect downstream amino acids
        consensus_aa = row.Consensus_AA_sequence
        variant_aa = row.Variant_AA_sequence
        
        # If sequences have different lengths, they definitely create changes
        if length(consensus_aa) != length(variant_aa)
            return false  # Not synonymous, has changes
        end
        
        # Check if any amino acid position differs
        for i in 1:min(length(consensus_aa), length(variant_aa))
            if consensus_aa[i] != variant_aa[i]
                return false  # Not synonymous, has changes
            end
        end
        return true  # All amino acids same, is synonymous
    else
        # For SNVs, use the original logic - check only the variant position
        locus = row.AA_Locus isa Int ? row.AA_Locus : tryparse(Int, row.AA_Locus)
        if locus === nothing || locus < 1 || locus > min(length(row.Consensus_AA_sequence), length(row.Variant_AA_sequence))
            false
        else
            row.Consensus_AA_sequence[locus] == row.Variant_AA_sequence[locus]
        end
    end
end, translated)
filtered_loci = filter(row -> begin
    # Check if this is an indel
    is_indel = length(row.Consensus) != length(row.Variant)
    
    if is_indel
        # For indels, check if any amino acid changes in the entire sequence
        consensus_aa = row.Consensus_AA_sequence
        variant_aa = row.Variant_AA_sequence
        
        # If sequences have different lengths, they definitely create changes
        if length(consensus_aa) != length(variant_aa)
            return true  # Has changes, not synonymous
        end
        
        # Check if any amino acid position differs
        for i in 1:min(length(consensus_aa), length(variant_aa))
            if consensus_aa[i] != variant_aa[i]
                return true  # Has changes, not synonymous
            end
        end
        return false  # All amino acids same, is synonymous
    else
        # For SNVs, use the original logic - check only the variant position
        locus = row.AA_Locus isa Int ? row.AA_Locus : tryparse(Int, row.AA_Locus)
        if locus === nothing || locus < 1 || locus > min(length(row.Consensus_AA_sequence), length(row.Variant_AA_sequence))
            false
        else
            row.Consensus_AA_sequence[locus] != row.Variant_AA_sequence[locus]
        end
    end
end, translated)
println("Dropped $(nrow(syn_loci)) synonymous loci.")

# Define locus-based peptide lengths to generate
substr_lengths = [8, 9, 10, 11]

flattened_peptides_df = add_peptides_columns!(
    filtered_loci,
    :Relative_Locus,
    :Consensus_AA_sequence,
    :Variant_AA_sequence,
    substr_lengths
)

transformed_df = separate_peptides(flattened_peptides_df)
transformed_df = deduplicate_peptides(transformed_df)

peptide_df = transformed_df

# Save to CSV
csv_file_path = joinpath(data_folder, "peptides_labels.csv")
CSV.write(csv_file_path, peptide_df)
println("peptides_labels.csv file has been written to: $csv_file_path")

# Save peptides to .pep file
write_peptides_file_no_headers(peptide_df, data_folder)