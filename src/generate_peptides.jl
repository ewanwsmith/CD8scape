include("./env.jl")

using DataFrames
using CSV

# Define a function to perform the operation with the folder path as an argument
function join_data(folder_path::String)
    # Construct file paths for frames and trajectories CSV files
    frames_path = joinpath(folder_path, "frames.csv")
    trajectories_path = joinpath(folder_path, "trajectories.csv")
    
    # Read the frames and trajectories data
    frames = CSV.read(frames_path, DataFrame)
    trajectories = CSV.read(trajectories_path, DataFrame)
    
    # Create Start and End columns by extracting the first Start and last End values
    frames[!, :Start] = [parse(Int, split(split(region, ";")[1], ",")[1]) for region in frames.Region]
    frames[!, :End] = [parse(Int, split(split(region, ";")[end], ",")[2]) for region in frames.Region]
    
    # Perform a cross join and then filter based on the condition
    result = crossjoin(trajectories, frames)
    result = filter(row -> row.Start <= row.Locus <= row.End, result)

    return result
end

function check_locus(df)
    # Generate the `Relative_Locus` column
    df[!, :Relative_Locus] = df.Locus .- df.Start

    # Create `Pulled_Base` by converting the extracted character to a String and then to String1
    df[!, :Pulled_Base] = [
        1 <= row.Relative_Locus <= length(row.Consensus_sequence) ? 
        String1(string(row.Consensus_sequence[row.Relative_Locus])) : missing
        for row in eachrow(df)
    ]

    # Check if Pulled_Base matches Consensus
    df[!, :Matches_Consensus] = df.Pulled_Base .== df.Consensus

    # Drop rows where Matches_Consensus is false (in-place)
    df = filter(row -> row.Matches_Consensus == true, df)

    # Drop Pulled_Base and Matches_Consensus columns from df
    select!(df, Not([:Pulled_Base, :Matches_Consensus]))

    return df
end

function edit_consensus_sequence(df)
    # Add an empty Variant_sequence column to hold the modified sequences
    df[!, :Variant_sequence] = fill("", nrow(df))

    # Iterate over each row in the DataFrame
    for row in eachrow(df)
        # Ensure Relative_Locus is within bounds
        if 1 <= row.Relative_Locus <= length(row.Consensus_sequence)
            # Convert Consensus_sequence to a character array for easy modification
            sequence_chars = collect(row.Consensus_sequence)
            
            # Modify the character at Relative_Locus to match Variant
            sequence_chars[row.Relative_Locus] = row.Variant[1]  # Assuming Variant is a single-character string or String1
            
            # Update Variant_sequence with the modified sequence
            row.Variant_sequence = join(sequence_chars)
        else
            # If out of bounds, set Variant_sequence to Consensus_sequence unmodified
            row.Variant_sequence = row.Consensus_sequence
        end
    end
    return df
end

function translate_sequences(df::DataFrame)
    # Codon to amino acid dictionary
    codon_table = Dict(
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

    # Function to translate a DNA sequence to an amino acid sequence
    function translate_dna_to_protein(dna_sequence::String)
        protein_sequence = ""
        for i in 1:3:length(dna_sequence)-2
            codon = dna_sequence[i:i+2]
            protein_sequence *= get(codon_table, codon, "?")  # Use "?" for unknown codons
        end
        return protein_sequence
    end

    # Add translated sequences as new columns
    df.Consensus_AA_sequence = [translate_dna_to_protein(seq) for seq in df.Consensus_sequence]
    df.Variant_AA_sequence = [translate_dna_to_protein(seq) for seq in df.Variant_sequence]

    return df
end

joined = join_data("/Users/e.smith.5/Documents/PhD/CD8scape/data/RSV_example")
checked = check_locus(joined)
edited = edit_consensus_sequence(checked)
translated = translate_sequences(edited)

function generate_peptides(sequence::String, locus::Int, substr_lengths::Vector{Int})
    peptides = String[]
    for substr_len in substr_lengths
        # Calculate safe indexing bounds using `eachindex`
        start_pos_min = max(first(eachindex(sequence)), locus - substr_len + 1)
        start_pos_max = min(last(eachindex(sequence)) - substr_len + 1, locus)

        for start_pos in start_pos_min:start_pos_max
            end_pos = start_pos + substr_len - 1
            
            # Ensure the substring includes the locus character
            if end_pos <= last(eachindex(sequence)) && (locus >= start_pos && locus <= end_pos)
                push!(peptides, sequence[start_pos:end_pos])
            end
        end
    end
    return peptides
end

function add_peptides_columns!(df::DataFrame, relative_locus_col::Symbol, consensus_col::Symbol, variant_col::Symbol, substr_lengths::Vector{Int})
    # Calculate AA_Locus from Relative_Locus with +1 adjustment
    df.AA_Locus = [Int(fld(locus, 3)) + 1 for locus in df[!, relative_locus_col]]

    # Initialize an empty DataFrame for storing the flattened peptide pairs
    flattened_peptides_df = DataFrame(
        Locus = Int[], Relative_Locus = Int[], AA_Locus = Int[],
        Consensus_Peptide = String[], Variant_Peptide = String[], Peptide_label = String[]
    )

    # Iterate through each row of the original DataFrame
    for row in eachrow(df)
        consensus_peptides = generate_peptides(row[consensus_col], row[:AA_Locus], substr_lengths)
        variant_peptides = generate_peptides(row[variant_col], row[:AA_Locus], substr_lengths)

        # Ensure the lists are of equal length
        if length(consensus_peptides) == length(variant_peptides)
            for j in eachindex(consensus_peptides)
                consensus_peptide = consensus_peptides[j]
                variant_peptide = variant_peptides[j]

                # Calculate the peptide length safely using `eachindex`
                peptide_length = last(eachindex(consensus_peptide))
                sequence_length = last(eachindex(row[consensus_col]))

                # Calculate the position of the locus within the peptide (1-based index)
                site_position = row[:AA_Locus] - (sequence_length - peptide_length)

                # Construct the Peptide_label
                peptide_label = "$(row.Consensus)$(row.Locus)$(row.Variant)_$site_position"

                # Add the new row to the flattened DataFrame
                push!(flattened_peptides_df, (
                    row.Locus, row.Relative_Locus, row.AA_Locus,
                    consensus_peptide, variant_peptide, peptide_label
                ))
            end
        else
            @warn "Mismatched peptide lists in row with Locus $(row.Locus). Skipping this row."
        end
    end

    return flattened_peptides_df
end

# Define the desired lengths for peptides
substr_lengths = [8, 9, 10, 11]

# Generate the flattened DataFrame with peptide pairs
flattened_peptides_df = add_peptides_columns!(translated, :Relative_Locus, :Consensus_AA_sequence, :Variant_AA_sequence, substr_lengths)

# Display the resulting flattened DataFrame
display(flattened_peptides_df)