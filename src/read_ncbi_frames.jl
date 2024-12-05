"""
A Julia script to process an NCBI FASTA file and extract regions.

Usage:
    julia script_name.jl <path_to_ncbi_fasta> [path_to_consensus_fasta]

    - If the path to the consensus FASTA file is not provided, the script will
      look for 'Consensus0.fa' in the same directory as the NCBI FASTA file.
"""

# Include environment settings
include("./env.jl")

using DataFrames
using CSV
using FilePathsBase  # For path manipulations

# Function to parse the FASTA header
function parse_header(header::AbstractString)
    # Remove the initial '>' if present
    header = strip(header, ['>'])

    # Initialize variables
    accession = ""
    region_str = ""  # String representation of regions
    description = ""

    # Check if header contains 'join('
    if occursin("join(", header)
        # Extract the 'join(...)' part
        join_part, rest_of_header = split(header, ')', limit=2)
        region_str = join_part[6:end]  # Remove 'join('
        rest_of_header = strip(rest_of_header)

        # Extract accession from the first region
        first_region = first(split(region_str, ','))
        acc_and_coords = split(first_region, ':')
        if !isempty(acc_and_coords)
            accession = acc_and_coords[1]
        end

        # Process region_str to remove accession numbers and reformat
        region_parts = split(region_str, ',')
        regions_formatted = []
        for region in region_parts
            # Remove accession number if present
            acc_and_coords = split(region, ':')
            coords = length(acc_and_coords) == 2 ? acc_and_coords[2] : acc_and_coords[1]
            # Split start and end positions
            if occursin("..", coords)
                start_pos, end_pos = split(coords, "..")
                formatted_region = "$start_pos,$end_pos"
            else
                formatted_region = "$coords,$coords"
            end
            push!(regions_formatted, formatted_region)
        end
        # Join regions with ';' between them
        region_str = join(regions_formatted, ';')

        description = rest_of_header
    else
        # Header does not contain 'join('
        parts = split(header, ' ', limit=2)
        acc_and_coords = split(parts[1], ':')
        accession = acc_and_coords[1]
        if length(acc_and_coords) == 2
            coords = acc_and_coords[2]
            # Split start and end positions
            if occursin("..", coords)
                start_pos, end_pos = split(coords, "..")
                region_str = "$start_pos,$end_pos"
            else
                region_str = "$coords,$coords"
            end
        else
            region_str = ""
        end

        description = length(parts) > 1 ? parts[2] : ""
    end

    # **Modification starts here**

    # Remove leading '|' from the description
    description = lstrip(description, '|')

    # Remove any content inside square brackets, including the brackets
    description = replace(description, r"\[.*?\]" => "")

    # Remove any extra whitespace
    description = strip(description)

    return accession, region_str, description
end

# Function to read FASTA metadata
function read_fasta_metadata(file_path::String)
    # Initialize an empty DataFrame with only the desired columns
    df = DataFrame(
        Description = String[],
        Region = String[]
    )

    # Open the file for reading
    open(file_path, "r") do io
        header = ""
        for line in eachline(io)
            line = strip(line)
            if isempty(line)
                continue
            elseif startswith(line, '>')
                # If there is a previous header, process it
                if !isempty(header)
                    accession, region_str, description = parse_header(header)
                    # Add the row to the DataFrame
                    push!(df, (
                        description,
                        region_str
                    ))
                end
                # Update the header
                header = line
            else
                # Skip sequence lines, as we don't need them
                continue
            end
        end
        # Process the last header
        if !isempty(header)
            accession, region_str, description = parse_header(header)
            push!(df, (
                description,
                region_str
            ))
        end
    end

    return df
end

# Function to read the consensus sequence from a FASTA file
function read_fasta_sequence(file_path::String)::String
    sequence = ""
    open(file_path, "r") do io
        for line in eachline(io)
            line = strip(line)
            if isempty(line) || startswith(line, '>')
                continue  # Skip headers and empty lines
            else
                sequence *= line  # Accumulate sequence lines
            end
        end
    end
    return sequence
end

# Function to extract and concatenate sequences based on regions
function extract_regions_sequence(region_str::String, sequence::String)::String
    # Split the region string into individual regions
    regions = split(region_str, ';')
    subsequences = []
    for region in regions
        start_end = split(region, ',')
        if length(start_end) != 2
            error("Invalid region format: $region")
        end
        start_pos = parse(Int, start_end[1])
        end_pos = parse(Int, start_end[2])
        # Check that positions are within bounds
        if start_pos < 1 || end_pos > length(sequence)
            error("Region positions out of bounds: $start_pos to $end_pos")
        end
        if start_pos > end_pos
            error("Start position greater than end position in region: $region")
        end
        # Extract the substring from the consensus sequence
        subseq = sequence[start_pos:end_pos]
        push!(subsequences, subseq)
    end
    # Concatenate subsequences if multiple regions
    concatenated_seq = join(subsequences, "")
    return concatenated_seq
end

# Main function to execute the script logic
function main()
    # Check if the script is called with the correct number of arguments
    if length(ARGS) < 1
        println("Usage: julia script_name.jl <path_to_ncbi_fasta> [path_to_consensus_fasta]")
        return
    end

    # Get the NCBI FASTA file path from the command-line arguments
    ncbi_file_path = ARGS[1]

    # Get the Consensus FASTA file path if provided
    consensus_fasta = ""
    if length(ARGS) >= 2
        consensus_fasta = ARGS[2]
    else
        # Use default value in the same directory as ncbi_file_path
        input_dir = dirname(abspath(ncbi_file_path))
        consensus_fasta = joinpath(input_dir, "Consensus0.fa")
        println("No consensus FASTA provided. Using default path: $consensus_fasta")
    end

    # Verify that the NCBI FASTA file exists
    if !isfile(ncbi_file_path)
        @error "NCBI FASTA file does not exist: $ncbi_file_path"
        return
    end

    # Verify that the Consensus FASTA file exists
    if !isfile(consensus_fasta)
        @error "Consensus FASTA file does not exist: $consensus_fasta"
        return
    end

    # Read metadata from the NCBI FASTA file
    regions_df = read_fasta_metadata(ncbi_file_path)

    # Read the consensus sequence
    consensus_seq = read_fasta_sequence(consensus_fasta)

    # Add a new column "Consensus_sequence" with the extracted sequences
    regions_df[!, :Consensus_sequence] = [extract_regions_sequence(regions_df.Region[i], consensus_seq) for i in 1:nrow(regions_df)]

    # Display the updated DataFrame
    println("Generated DataFrame:")
    println(regions_df)

    # Define the output CSV path in the same directory as the NCBI file
    input_dir = dirname(abspath(ncbi_file_path))
    output_csv_path = joinpath(input_dir, "frames.csv")

    # Save the DataFrame to frames.csv in the input directory
    CSV.write(output_csv_path, regions_df)

    println("DataFrame has been saved to $output_csv_path")
end

# Execute the main function
main()