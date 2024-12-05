# Include environment settings
include("./env.jl")

using DataFrames
using CSV
using FilePathsBase

# Function to read the .dat file and create the DataFrame with "Region", "Consensus_sequence", and "Description" columns
function read_dat_file(filepath::String)::DataFrame
    # Read all lines from the file
    lines = readlines(filepath)
    
    # Initialize vectors to store the data
    regions = String[]
    consensus = String[]
    descriptions = String[]
    
    # Total number of lines
    n = length(lines)
    
    # Initialize frame counter
    frame_counter = 1
    
    # Iterate over the lines in steps of 2 (assuming 2 lines per entry)
    i = 1
    while i <= n
        # Ensure that we have at least two lines left for a complete block
        if i + 1 > n
            @warn "Incomplete block starting at line $i. Skipping."
            break
        end
        
        # Parse the Start and End from the first line
        start_end = split(strip(lines[i]))
        if length(start_end) < 2
            @warn "Line $i does not contain two elements. Skipping this block."
            i += 2
            continue
        end
        start_val = tryparse(Int, start_end[1])
        end_val = tryparse(Int, start_end[2])
        
        if isnothing(start_val) || isnothing(end_val)
            @warn "Invalid Start or End values at line $i. Skipping this block."
            i += 2
            continue
        end
        
        # Combine Start and End into "Region" as "Start,End"
        region = "$(start_val),$(end_val)"
        push!(regions, region)
        
        # Read the Consensus_sequence from the next line
        consensus_seq = strip(lines[i + 1])
        push!(consensus, consensus_seq)
        
        # Generate the Description as "Frame_n"
        description = "Frame_$(frame_counter)"
        push!(descriptions, description)
        
        # Increment the frame counter
        frame_counter += 1
        
        # Move to the next block
        i += 2
    end
    
    # Create the DataFrame with "Region", "Consensus_sequence", and "Description"
    frames_df = DataFrame(
        Region = regions,
        Consensus_sequence = consensus,
        Description = descriptions
    )
    
    return frames_df
end

# Function to save DataFrame as CSV in the same directory as the input file
function save_dataframe_as_csv(df::DataFrame, input_filepath::String, output_filename::String="frames.csv")
    # Extract the directory from the input file path
    directory = dirname(input_filepath)
    
    # Construct the output file path
    output_filepath = joinpath(directory, output_filename)
    
    # Write the DataFrame to CSV
    CSV.write(output_filepath, df)
    
    println("DataFrame successfully saved to $output_filepath")
end

# Main function to orchestrate the processing
function main()
    # Check if the script is called with the correct number of arguments
    if length(ARGS) < 1
        println("Usage: julia script_name.jl <path_to_frames.dat>")
        return
    end

    # Get the file path from the command-line arguments
    frames_filepath = ARGS[1]

    # Check if the provided file exists
    if !isfile(frames_filepath)
        @error "The file '$frames_filepath' does not exist."
        return
    end
    
    # Read the .dat file and create the DataFrame
    frames_df = read_dat_file(frames_filepath)
    
    # Display the DataFrame (optional)
    println("Generated DataFrame:")
    println(frames_df)
    
    # Save the DataFrame as frames.csv in the same directory as frames_filepath
    save_dataframe_as_csv(frames_df, frames_filepath)
end

# Execute the main function
main()