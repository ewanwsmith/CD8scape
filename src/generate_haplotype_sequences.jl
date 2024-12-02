#!/usr/bin/env julia

# Include environment settings
include("./env.jl")

using DataFrames
using CSV
using ArgParse
# using FilePathsBase  # Removed to prevent path handling conflicts

# Function to read the locus file and return a DataFrame
function read_locus_file(filepath::String)::DataFrame
    """
    Reads a locus file and returns a DataFrame with columns:
    - `locus` (Int): The position of the locus.
    - `original_base` (String): The original base at the locus.
    - `variant_base` (String): The variant base at the locus.
    """
    # Initialize arrays to store the data
    loci = Int[]
    original_bases = String[]
    variant_bases = String[]

    # Open and read the file line by line
    open(filepath, "r") do f
        for line in eachline(f)
            tokens = split(strip(line))
            if length(tokens) >= 3
                # Parse the first three columns
                locus = tryparse(Int, tokens[1])
                if locus === nothing
                    @warn "Invalid locus value in line: $line"
                    continue
                end
                original_base = tokens[2]
                variant_base = tokens[3]

                # Append the parsed data to the arrays
                push!(loci, locus)
                push!(original_bases, original_base)
                push!(variant_bases, variant_base)
            else
                @warn "Line with insufficient data: $line"
            end
        end
    end

    # Create a DataFrame from the collected data
    df = DataFrame(
        locus = loci,
        original_base = original_bases,
        variant_base = variant_bases
    )
    return df
end

# Function to read the haplotypes file and return a DataFrame
function read_haps_file(filepath::String)::DataFrame
    """
    Reads a haplotypes file and returns a DataFrame with one column:
    - `haps` (String): The haplotype sequences.
    """
    # Initialize an array to store the haplotypes
    haps = String[]
    
    # Define a regular expression pattern to match lines starting with haplotype sequences
    pattern = r"^[ACGT]+"

    # Open and read the file line by line
    open(filepath, "r") do f
        for line in eachline(f)
            line = strip(line)
            if isempty(line)
                continue  # Skip empty lines
            end
            if occursin(pattern, line)
                # Split the line into tokens
                tokens = split(line)
                # The first token is the haplotype
                hap = tokens[1]
                # Verify that the rest of the tokens are numbers
                if all(x -> tryparse(Float64, x) !== nothing, tokens[2:end])
                    push!(haps, hap)
                else
                    @warn "Skipping line with non-numeric data: $line"
                end
            else
                # Line doesn't start with a haplotype sequence; skip it
                continue
            end
        end
    end

    # Create a DataFrame with the collected haplotypes
    df = DataFrame(haps = haps)
    return df
end

# Function to construct the bases DataFrame from loci and haplotypes
function construct_bases_dataframe(loci_df::DataFrame, haps_df::DataFrame)::DataFrame
    """
    Constructs a DataFrame combining loci information with haplotypes.

    Parameters:
    - `loci_df` (DataFrame): DataFrame containing loci information.
    - `haps_df` (DataFrame): DataFrame containing haplotype sequences.

    Returns:
    - `bases_df` (DataFrame): DataFrame where each row corresponds to a haplotype, and columns represent the base at each locus.
    """
    # Extract the unique loci and sort them
    loci = unique(loci_df.locus)
    sort!(loci)
    num_loci = length(loci)

    # Initialize a DataFrame with columns for each locus
    # Include a column for haplotype identifiers
    bases_df = DataFrame()

    # Add the haplotype strings as an identifier column
    bases_df[!, :Haplotype] = haps_df.haps  # Renamed 'haps' to 'Haplotype'

    # Check that each haplotype string has the same length as the number of loci
    for (i, hap_string) in enumerate(haps_df.haps)
        hap_length = length(hap_string)
        if hap_length != num_loci
            error("Length of haplotype string (length $hap_length) does not match number of loci ($num_loci) for haplotype at row $(i).")
        end
    end

    # Add columns for each locus, filling in the bases from each haplotype
    # The position in the haplotype string corresponds to the sorted loci positions
    for (idx, locus) in enumerate(loci)
        # Column name based on the locus
        column_name = Symbol("locus_$(locus)")
        # Extract the base at this position for each haplotype
        bases_df[!, column_name] = [hap[idx] for hap in haps_df.haps]
    end

    return bases_df
end

# Function to read frames.csv into a DataFrame
function read_frames_csv(file_path::String)::DataFrame
    """
    Reads a CSV file specified by `file_path` into a DataFrame.

    Parameters:
    - `file_path` (String): The path to the frames.csv file.

    Returns:
    - `df` (DataFrame): The DataFrame containing the data from frames.csv.
    """
    # Check if the file exists
    if !isfile(file_path)
        error("The file '$file_path' does not exist.")
    end

    # Read the CSV file into a DataFrame
    df = CSV.read(file_path, DataFrame)

    return df
end

# Function to join frames and bases DataFrames with multiple rows per haplotype-frame combination
function join_frames_bases_multiple_rows(frames_df::DataFrame, bases_df::DataFrame)::DataFrame
    """
    Joins the frames DataFrame with the bases DataFrame based on locus positions,
    producing multiple rows per haplotype and frame combination.
    In each row, only the loci that match the frame are present; others are missing.
    The 'Region' column from frames_df is included in the final DataFrame.

    Parameters:
    - `frames_df` (DataFrame): DataFrame containing frames information with columns:
        - `Description` (String)
        - `Region` (String): "start,end" format
        - `Consensus_sequence` (String)
    - `bases_df` (DataFrame): DataFrame containing haplotypes and locus data with columns:
        - `Haplotype` (String): Haplotype identifiers
        - `locus_XXX` (String): Base at locus XXX

    Returns:
    - `result_df` (DataFrame): DataFrame combining information from both DataFrames,
      with multiple rows per haplotype and frame, and loci outside the frame set to missing.
    """

    # Extract locus columns and their positions
    locus_cols = filter(col -> startswith(col, "locus_"), names(bases_df))
    # Extract locus numbers
    locus_numbers = parse.(Int, replace.(locus_cols, "locus_" => ""))

    # Parse 'Region' in frames_df into 'start' and 'end'
    frames_df = deepcopy(frames_df)
    frames_df[!, :start] = parse.(Int, first.(split.(frames_df.Region, ",")))
    frames_df[!, :end] = parse.(Int, last.(split.(frames_df.Region, ",")))

    # Initialize an empty DataFrame to collect results
    result_rows = DataFrame()

    # Iterate over each haplotype
    for hap_row in eachrow(bases_df)
        haplotype = hap_row.Haplotype  # Use 'Haplotype' column
        # Iterate over each frame
        for frame_row in eachrow(frames_df)
            frame_start = frame_row.start
            frame_end = frame_row.end
            # Determine which loci fall within the frame's region
            locus_in_frame = (locus_numbers .>= frame_start) .& (locus_numbers .<= frame_end)
            if any(locus_in_frame)
                # Create a new row
                new_row = Dict{Symbol, Any}()
                # Add haplotype identifier
                new_row[:Haplotype] = haplotype
                # Add frame information
                new_row[:Description] = frame_row.Description
                new_row[:Region] = frame_row.Region  # Include the 'Region' column
                new_row[:Consensus_sequence] = frame_row.Consensus_sequence
                # For each locus, assign base or missing based on whether it is in the frame
                for (i, locus_col) in enumerate(locus_cols)
                    if locus_in_frame[i]
                        new_row[Symbol(locus_col)] = hap_row[Symbol(locus_col)]
                    else
                        new_row[Symbol(locus_col)] = missing
                    end
                end
                # Append new_row to result_rows
                append!(result_rows, DataFrame(new_row))
            end
        end
    end

    return result_rows
end

# Function to adjust locus positions based on the 'start' value in the 'Region' column and reshape back to wide format
function adjust_locus_positions(df::DataFrame)::DataFrame
    """
    Adjusts the locus positions in the DataFrame by subtracting the 'start' value from the 'Region' column.
    The DataFrame is reshaped into long format to accommodate per-row adjustments, and then reshaped back to wide format.

    Parameters:
    - `df` (DataFrame): The DataFrame to process, which includes 'Region' column and 'locus_' columns.

    Returns:
    - `adjusted_df` (DataFrame): The adjusted DataFrame in wide format with adjusted locus positions as column names.
    """
    # Extract locus columns
    locus_cols = filter(col -> startswith(col, "locus_"), names(df))

    # Reshape the DataFrame into long format
    df_long = stack(df, locus_cols, variable_name="locus_col", value_name="base")

    # Extract locus number from 'locus_col'
    df_long[!, :locus] = parse.(Int, replace.(df_long.locus_col, "locus_" => ""))

    # Extract 'start' from 'Region' column
    df_long[!, :start] = parse.(Int, first.(split.(df_long.Region, ",")))

    # Compute relative position
    df_long[!, :relative_position] = df_long.locus .- df_long.start

    # Create a new adjusted locus column name
    df_long[!, :adjusted_locus_col] = "locus_" .* string.(df_long.relative_position)

    # Now, pivot the DataFrame back to wide format
    adjusted_df = unstack(df_long, [:Haplotype, :Description, :Region, :Consensus_sequence], :adjusted_locus_col, :base)

    return adjusted_df
end

# Function to create Haplotype_Sequence by applying locus substitutions to Consensus_sequence
function create_haplotype_sequences(df::DataFrame)::DataFrame
    """
    For each row in the DataFrame, creates a new sequence by applying the base substitutions
    from the locus columns to the Consensus_sequence, resulting in Haplotype_Sequence.

    Parameters:
    - `df` (DataFrame): The DataFrame with adjusted locus positions and locus columns.

    Returns:
    - `df` (DataFrame): The DataFrame with an additional column 'Haplotype_Sequence'.
    """
    # Copy the DataFrame to avoid modifying the original
    df = deepcopy(df)

    # Get the locus columns
    locus_cols = filter(col -> startswith(col, "locus_"), names(df))

    # Extract the adjusted positions from the locus column names
    # Map from column names to positions
    adjusted_positions = Dict{String, Int}()
    for col in locus_cols
        # Extract the position from the column name
        # If the column name is like "locus_0", extract 0
        pos_str = replace(col, "locus_" => "")
        pos = parse(Int, pos_str)
        adjusted_positions[col] = pos
    end

    # Prepare a vector to store the Haplotype_Sequence
    haplotype_sequences = Vector{String}(undef, nrow(df))

    # For each row, apply the substitutions
    for (i, row) in enumerate(eachrow(df))
        # Get the consensus sequence
        consensus_seq = row.Consensus_sequence
        # Convert to a mutable vector of characters
        haplotype_seq = collect(consensus_seq)

        # For each locus column
        for col in locus_cols
            base = row[col]
            if !ismissing(base)
                # Get the adjusted position
                pos = adjusted_positions[col]
                # Julia strings are 1-based, so adjust position accordingly
                index = pos + 1
                # Check if index is within bounds
                if 1 <= index <= length(haplotype_seq)
                    # Replace the base at this position
                    haplotype_seq[index] = base
                else
                    @warn "Index $index is out of bounds for sequence of length $(length(haplotype_seq)) in row $i"
                end
            end
        end
        # Convert back to string and store in the vector
        haplotype_sequences[i] = join(haplotype_seq)
    end

    # Add the new column to the DataFrame
    df.Haplotype_Sequence = haplotype_sequences

    return df
end

# Function to create 'Label' column by joining 'Haplotype' and modified 'Description'
function create_label_column(df::DataFrame)::DataFrame
    """
    Creates a new column 'Label' by joining 'Haplotype' and 'Description' columns.
    Spaces in 'Description' are replaced with underscores before joining.

    Parameters:
    - `df` (DataFrame): The DataFrame with 'Haplotype' and 'Description' columns.

    Returns:
    - `df` (DataFrame): The DataFrame with an additional column 'Label'.
    """
    # Copy the DataFrame to avoid modifying the original
    df = deepcopy(df)

    # Prepare the 'Label' column
    labels = Vector{String}(undef, nrow(df))

    for (i, row) in enumerate(eachrow(df))
        # Get the Haplotype and Description
        haplotype = row.Haplotype
        description = row.Description

        # Replace spaces with underscores in Description
        description_no_spaces = replace(description, " " => "_")

        # Join Haplotype and modified Description with '-'
        label = string(haplotype, "-", description_no_spaces)

        # Store in the labels vector
        labels[i] = label
    end

    # Add the 'Label' column to the DataFrame
    df.Label = labels

    return df
end

# Function to create output DataFrame 'output_df' with 'Sequence' and 'Label' columns
function create_output_df(df::DataFrame)::DataFrame
    """
    Creates an output DataFrame 'output_df' with columns 'Sequence' and 'Label'.
    For each unique 'Description', it adds the 'Consensus_sequence' with a Label 'Consensus-<Description>'.
    Then, it adds all 'Haplotype_Sequence's with their corresponding 'Label's.

    Parameters:
    - `df` (DataFrame): The final DataFrame with 'Haplotype_Sequence', 'Consensus_sequence', 'Description', and 'Label' columns.

    Returns:
    - `output_df` (DataFrame): The output DataFrame with 'Sequence' and 'Label' columns.
    """
    # Initialize an empty DataFrame
    output_df = DataFrame(Sequence=String[], Label=String[])

    # Get unique Descriptions
    unique_descriptions = unique(df.Description)

    # For each unique Description, add the Consensus_sequence
    for desc in unique_descriptions
        # Filter rows with this Description
        desc_rows = df[df.Description .== desc, :]
        if nrow(desc_rows) > 0
            # Get the first Consensus_sequence (they should be the same for the same Description)
            consensus_sequence = desc_rows.Consensus_sequence[1]

            # Create the Label "Consensus-<Description>" with spaces replaced by underscores
            desc_label = replace(desc, " " => "_")
            consensus_label = "Consensus-" * desc_label

            # Append to output_df
            push!(output_df, (Sequence=consensus_sequence, Label=consensus_label))
        end
    end

    # Now, add all Haplotype_Sequences with their Labels
    haplotype_sequences = df.Haplotype_Sequence
    labels = df.Label  # Use 'Label' with uppercase 'L'

    # Append these to output_df
    for i in 1:length(haplotype_sequences)
        push!(output_df, (Sequence=haplotype_sequences[i], Label=labels[i]))
    end

    return output_df
end

# Function to write sequences and labels to a FASTA file
function write_fasta_file(output_df::DataFrame, filepath::String)
    """
    Writes the sequences and labels from output_df to a FASTA file.

    Parameters:
    - `output_df` (DataFrame): The DataFrame with 'Sequence' and 'Label' columns.
    - `filepath` (String): The path to the output FASTA file.
    """
    # Open the file for writing
    open(filepath, "w") do io
        for row in eachrow(output_df)
            # Write the header line with the label
            println(io, ">", row.Label)
            # Write the sequence
            println(io, row.Sequence)
        end
    end
end

# Function to parse command-line arguments
function get_args()::Dict{Symbol, Any}
    """
    Parses command-line arguments using ArgParse.

    Returns:
    - args (Dict{Symbol, Any}): A dictionary of parsed arguments with Symbol keys.
    """
    s = ArgParseSettings(description="Process locus, haplotypes, and frames files to generate a FASTA file.")

    @add_arg_table s begin
        "--locus", "-l"
        help = "Path to the locus file."
        arg_type = String
        required = true

        "--haps", "-p"  # Changed short option from '-h' to '-p'
        help = "Path to the haplotypes file."
        arg_type = String
        required = true

        "--frames", "-f"
        help = "Path to the frames CSV file."
        arg_type = String
        required = true

        "--output", "-o"
        help = "Path to the output FASTA file. If not provided, defaults to the frames file directory with name 'haplotype_sequences.fa'."
        arg_type = String
        required = false
    end

    # Parse the arguments using ArgParse.parse_args
    args = ArgParse.parse_args(ARGS, s)

    # Convert String keys to Symbol keys
    args_sym = Dict{Symbol, Any}(Symbol.(keys(args)) .=> values(args))

    return args_sym
end

# Main function to run the script
function main()
    # Parse command-line arguments
    args = get_args()

    # Assign arguments to variables using Symbol keys
    locus_file_path = args[:locus]
    haps_file_path = args[:haps]
    frames_csv_path = args[:frames]

    # Determine the output FASTA file path
    if args[:output] !== nothing
        output_fasta_path = args[:output]
    else
        # Extract the directory path from frames_csv_path
        frames_directory = dirname(frames_csv_path)
        println("Determined frames_directory: $frames_directory")
        # Construct the output FASTA file path in the same directory with the desired name
        output_fasta_path = joinpath(frames_directory, "haplotype_sequences.fa")  # Default file name
    end

    # Debugging Statements
    println("Output FASTA File: $output_fasta_path")
    println("Type of output_fasta_path: ", typeof(output_fasta_path))

    # Validate input files
    for filepath in [locus_file_path, haps_file_path, frames_csv_path]
        if !isfile(filepath)
            error("Input file '$filepath' does not exist.")
        end
    end

    println("\nInput Files:")
    println("  Locus File: $locus_file_path")
    println("  Haplotypes File: $haps_file_path")
    println("  Frames CSV File: $frames_csv_path")
    println("Output FASTA File: $output_fasta_path\n")

    # Read data from files
    loci_df = read_locus_file(locus_file_path)
    haps_df = read_haps_file(haps_file_path)
    bases_df = construct_bases_dataframe(loci_df, haps_df)
    frames_df = read_frames_csv(frames_csv_path)

    # Display the Bases DataFrame
    println("Constructed Bases DataFrame:")
    display(first(bases_df, 5))  # Display first 5 rows for brevity
    println("... (Total Rows: $(nrow(bases_df)))\n")

    # Join frames and bases
    joined_df = join_frames_bases_multiple_rows(frames_df, bases_df)

    # Display the Joined DataFrame
    println("Joined DataFrame:")
    display(first(joined_df, 5))  # Display first 5 rows for brevity
    println("... (Total Rows: $(nrow(joined_df)))\n")

    # Adjust locus positions based on 'Region' start value and reshape back to wide format
    adjusted_df = adjust_locus_positions(joined_df)

    # Display the Adjusted DataFrame
    println("Adjusted DataFrame with Adjusted Locus Positions:")
    display(first(adjusted_df, 5))  # Display first 5 rows for brevity
    println("... (Total Rows: $(nrow(adjusted_df)))\n")

    # Create haplotype sequences
    final_df = create_haplotype_sequences(adjusted_df)

    # Create the Label column
    final_df = create_label_column(final_df)

    # Display the final DataFrame
    println("Final DataFrame with Haplotype Sequences and Labels:")
    display(first(final_df, 5))  # Display first 5 rows for brevity
    println("... (Total Rows: $(nrow(final_df)))\n")

    # Create the output DataFrame
    output_df = create_output_df(final_df)

    # Display the output DataFrame
    println("Output DataFrame:")
    display(first(output_df, 5))  # Display first 5 rows for brevity
    println("... (Total Rows: $(nrow(output_df)))\n")

    # Write the sequences to a FASTA file in the specified output path
    write_fasta_file(output_df, output_fasta_path)

    println("FASTA file has been written to: $output_fasta_path")
end

# Execute the main function
main()